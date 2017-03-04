//============================================================================//
// PRad Cluster Reconstruction Method                                         //
// Reconstruct the cluster using island algorithm from PrimEx                 //
//                                                                            //
// Ilya Larin, adapted GAMS island algorithm for HyCal in PrimEx.             //
// Maxime Levillain & Weizhi Xiong, developed a new method based on PrimEx    //
//                                  method, has a much better performance on  //
//                                  splitting, but has a slightly worse       //
//                                  resolution 08/01/2016                     //
// Weizhi Xiong, adapted island algorithm for PRad, reduced the discretized   //
//               energy value from 10 MeV to 0.1 MeV, wrote a C++ wrapper for //
//               the fortran code. 09/28/2016                                 //
// Chao Peng, rewrote the whole fortran code into C++. 11/21/2016             //
//============================================================================//

#include <cmath>
#include <algorithm>
#include <iostream>
#include "PRadIslandCluster.h"
#include "PRadClusterProfile.h"
#include "PRadHyCalDetector.h"
#include "PRadADCChannel.h"
#include "PRadTDCChannel.h"

const PRadClusterProfile &__ic_prof = PRadClusterProfile::Instance();

PRadIslandCluster::PRadIslandCluster(const std::string &path)
{
    // configuration
    Configure(path);
}

PRadIslandCluster::~PRadIslandCluster()
{
    // place holder
}

PRadHyCalCluster *PRadIslandCluster::Clone()
const
{
    return new PRadIslandCluster(*this);
}

void PRadIslandCluster::Configure(const std::string &path)
{
    PRadHyCalCluster::Configure(path);

    bool verbose = !path.empty();

    split_iter = getDefConfig<unsigned int>("Split Iteration", 6, verbose);
    least_share = getDefConfig<float>("Least Split Fraction", 0.01, verbose);
    bool corner = getDefConfig<bool>("Corner Connection", false, verbose);
    if(corner)
        adj_dist = CORNER_ADJACENT;
    else
        adj_dist = SIDE_ADJACENT;

    // set the min module energy for all the module type
    float univ_min_energy = getDefConfig<float>("Min Module Energy", 0., false);
    min_module_energy.resize(PRadHyCalModule::Max_Type, univ_min_energy);

    // update the min module energy if some type is specified
    // the key is "Min Module Energy [typename]"
    for(unsigned int i = 0; i < min_module_energy.size(); ++i)
    {
        // determine key name
        std::string type = PRadHyCalModule::get_module_type_name(i);
        std::string key = "Min Module Energy [" + type + "]";
        auto value = GetConfigValue(key);
        if(!value.IsEmpty())
            min_module_energy[i] = value.Float();
    }
}



//============================================================================//
// Method based on the code from I. Larin for PrimEx                          //
//============================================================================//
#ifdef ISLAND_FINE_SPLIT

void PRadIslandCluster::FormCluster(std::vector<ModuleHit> &hits,
                                    std::vector<ModuleCluster> &clusters)
const
{
    // clear container first
    clusters.clear();

    std::vector<std::vector<ModuleHit*>> groups;

    // group adjacent hits
    groupHits(hits, groups);

    // try to split the group
    for(auto &group : groups)
    {
        splitCluster(group, clusters);
    }
}

// group adjacent hits into raw clusters
// list is suitable for the merge process, but it has poor performance while loop
// over all the elements, the test with real data prefers vector as container
void PRadIslandCluster::groupHits(std::vector<ModuleHit> &hits,
                                  std::vector<std::vector<ModuleHit*>> &groups)
const
{
    // roughly combine all adjacent hits
    for(auto &hit : hits)
    {
        if(hit.energy < min_module_energy.at(hit.geo.type))
            continue;

        // not belong to any existing cluster
        if(!fillClusters(hit, groups)) {
            std::vector<ModuleHit*> new_group;
            new_group.reserve(50);
            new_group.push_back(&hit);
            groups.emplace_back(std::move(new_group));
        }
    }

    // merge adjacent groups

    for(auto it = groups.begin(); it != groups.end(); ++it)
    {
        auto it_next = it;
        while(++it_next != groups.end())
        {
            if(checkAdjacent(*it, *it_next)) {
                it_next->insert(it_next->end(), it->begin(), it->end());
                groups.erase(it--);
                break;
            }
        }
    }
}

bool PRadIslandCluster::fillClusters(ModuleHit &hit,
                                     std::vector<std::vector<ModuleHit*>> &groups)
const
{
    for(auto &group : groups)
    {
        for(auto &prev_hit : group)
        {
            // it belongs to a existing cluster
            if(PRadHyCalDetector::hit_distance(hit, *prev_hit) < adj_dist) {
                group.push_back(&hit);
                return true;
            }
        }
    }

    return false;
}

inline bool PRadIslandCluster::checkAdjacent(const std::vector<ModuleHit*> &g1,
                                             const std::vector<ModuleHit*> &g2)
const
{
    for(auto &m1 : g1)
    {
        for(auto &m2 : g2)
        {
            if(PRadHyCalDetector::hit_distance(*m1, *m2) < adj_dist) {
                return true;
            }
        }
    }

    return false;
}

// some global container and function to help splitting and improve performance
#define SPLIT_MAX_HITS 100
#define SPLIT_MAX_CLUSTERS 10
float __ic_frac[SPLIT_MAX_HITS][SPLIT_MAX_CLUSTERS];
float __ic_tot_frac[SPLIT_MAX_HITS];

inline void __ic_sum_frac(size_t hits, size_t maximums)
{
    for(size_t i = 0; i < hits; ++i)
    {
        __ic_tot_frac[i] = 0;
        for(size_t j = 0; j < maximums; ++j)
            __ic_tot_frac[i] += __ic_frac[i][j];
    }
}

// split one group into several clusters
void PRadIslandCluster::splitCluster(const std::vector<ModuleHit*> &group,
                                     std::vector<ModuleCluster> &clusters)
const
{
    // find local maximum
    auto maximums = findMaximums(group);

    // no cluster center found
    if(maximums.empty())
        return;

    if((maximums.size() == 1) ||                    // only 1 cluster
       (group.size() >= SPLIT_MAX_HITS) ||          // too many hits
       (maximums.size() >= SPLIT_MAX_CLUSTERS)) {   // too many local maximums
        // create cluster based on the center
        clusters.emplace_back(*maximums.front());
        auto &cluster = clusters.back();

        for(auto &hit : group)
            cluster.AddHit(*hit);
    // split hits between several maximums
    } else {
        splitHits(maximums, group, clusters);
    }
}

// find local maximums in a group of adjacent hits
std::vector<ModuleHit*> PRadIslandCluster::findMaximums(const std::vector<ModuleHit*> &hits)
const
{
    std::vector<ModuleHit*> local_max;
    local_max.reserve(20);
    for(auto it = hits.begin(); it != hits.end(); ++it)
    {
        auto &hit1 = *it;
        if(hit1->energy < min_center_energy)
            continue;

        bool maximum = true;
        for(auto &hit2 : hits)
        {
            if(hit1 == hit2)
                continue;

            // we count corner in
            if((PRadHyCalDetector::hit_distance(*hit1, *hit2) < CORNER_ADJACENT) &&
               (hit2->energy > hit1->energy)) {
                maximum = false;
                break;
            }
        }

        if(maximum) {
            local_max.push_back(hit1);
        }
    }

    return local_max;
}

// split hits between several local maximums inside a cluster group
void PRadIslandCluster::splitHits(const std::vector<ModuleHit*> &maximums,
                                  const std::vector<ModuleHit*> &hits,
                                  std::vector<ModuleCluster> &clusters)
const
{
    // initialize fractions
    for(size_t i = 0; i < maximums.size(); ++i)
    {
        auto &center = *maximums.at(i);
        for(size_t j = 0; j < hits.size(); ++j)
        {
            auto &hit = *hits.at(j);
            __ic_frac[j][i] = __ic_prof.GetProfile(center, hit).frac*center.energy;
        }
    }

    // do iteration to evaluate the share of hits between several maximums
    evalFraction(hits, maximums, split_iter);

    // done iteration, add cluster according to the final share of energy
    for(size_t i = 0; i < maximums.size(); ++i)
    {
        clusters.emplace_back(*maximums.at(i));
        auto &cluster = clusters.back();

        for(size_t j = 0; j < hits.size(); ++j)
        {
            if(__ic_frac[j][i] == 0.)
                continue;

            // too small share, treat as zero
            if(__ic_frac[j][i]/__ic_tot_frac[j] < least_share) {
                __ic_tot_frac[j] -= __ic_frac[j][i];
                continue;
            }

            ModuleHit new_hit(*hits.at(j));
            new_hit.energy *= __ic_frac[j][i]/__ic_tot_frac[j];
            cluster.AddHit(new_hit);

            // update the center energy
            if(new_hit == cluster.center)
                cluster.center.energy = new_hit.energy;
        }
    }
}

inline void PRadIslandCluster::evalFraction(const std::vector<ModuleHit*> &hits,
                                            const std::vector<ModuleHit*> &maximums,
                                            size_t iters)
const
{
    // temp containers for reconstruction
    BaseHit temp[POS_RECON_HITS];

    // iterations to refine the split energies
    while(iters-- > 0)
    {
        __ic_sum_frac(hits.size(), maximums.size());
        for(size_t i = 0; i < maximums.size(); ++i)
        {
            // cluster center reconstruction
            auto &center = *maximums.at(i);
            float tot_E = 0.;
            int count = 0;
            for(size_t j = 0; j < hits.size(); ++j)
            {
                auto &hit = *hits.at(j);
                if(__ic_frac[j][i] == 0.)
                    continue;

                // using 3x3 to reconstruct hit position
                if(PRadHyCalDetector::hit_distance(center, hit) < CORNER_ADJACENT) {
                    temp[count].x = hit.geo.x;
                    temp[count].y = hit.geo.y;
                    temp[count].E = hit.energy*__ic_frac[j][i]/__ic_tot_frac[j];
                    tot_E += temp[count].E;
                    count++;
                }
            }

            BaseHit recon;
            PRadHyCalCluster::reconstructPos(temp, count, &recon);

            // update profile with the reconstructed center
            for(size_t j = 0; j < hits.size(); ++j)
            {
                auto &hit = *hits.at(j);
                __ic_frac[j][i] = __ic_prof.GetProfile(recon.x, recon.y, hit).frac*tot_E;
            }
        }
    }
    __ic_sum_frac(hits.size(), maximums.size());
}



//============================================================================//
// Method based on code from M. Levillain and W. Xiong                        //
//============================================================================//
#else

void PRadIslandCluster::FormCluster(std::vector<ModuleHit> &hits,
                                    std::vector<ModuleCluster> &clusters)
const
{
    // clear container first
    clusters.clear();

    // form clusters with high energy hit seed
    groupHits(hits, clusters);
}

void PRadIslandCluster::groupHits(std::vector<ModuleHit> &hits,
                                  std::vector<ModuleCluster> &clusters)
const
{
    // sort hits by energy
    std::sort(hits.begin(), hits.end(),
              [] (const ModuleHit &m1, const ModuleHit &m2)
              {
                  return m1.energy > m2.energy;
              });

    // loop over all hits
    for(auto &hit : hits)
    {
        // less than min module energy, ignore this hit
        if(hit.energy < min_module_energy.at(hit.geo.type))
            continue;

        // not belongs to any cluster, and the energy is larger than center threshold
        if(!fillClusters(hit, clusters) && (hit.energy > min_center_energy))
        {
            clusters.emplace_back(hit);
            clusters.back().AddHit(hit);
        }
    }
}

bool PRadIslandCluster::fillClusters(ModuleHit &hit, std::vector<ModuleCluster> &c)
const
{
    std::vector<unsigned int> indices;
    indices.reserve(5);

    for(unsigned int i = 0; i < c.size(); ++i)
    {
        for(auto &prev_hit : c.at(i).hits)
        {
            if(PRadHyCalDetector::hit_distance(hit, prev_hit) < CORNER_ADJACENT) {
                indices.push_back(i);
                break;
            }
        }
    }

    // it belongs to no cluster
    if(indices.empty())
        return false;

    // it belongs to single cluster
    if(indices.size() == 1) {
        c.at(indices.front()).AddHit(hit);
        return true;
    }

    // it belongs to several clusters
    return splitHit(hit, c, indices);
}

// split hit that belongs to several clusters
bool PRadIslandCluster::splitHit(ModuleHit &hit,
                                 std::vector<ModuleCluster> &clusters,
                                 std::vector<unsigned int> &indices)
const
{
    // energy fraction
    float frac[indices.size()];
    float total_frac = 0.;

    // rough splitting, only use the center position
    // to refine it, use a reconstruction position or position from other detector
    for(unsigned int i = 0; i < indices.size(); ++i)
    {
        auto &center = clusters.at(indices.at(i)).center;
        // we are comparing the relative amount of energy to be shared, so use of
        // center energy should be equivalent to total cluster energy
        frac[i] = __ic_prof.GetProfile(center, hit).frac * center.energy;
        total_frac += frac[i];
    }

    // this hit is too far away from all clusters, discard it
    if(total_frac == 0.) {
        return false;
    }

    // add hit to all clusters
    for(unsigned int i = 0; i < indices.size(); ++i)
    {
        // this cluster has no share of the hit
        if(frac[i] == 0.)
            continue;

        auto &cluster = clusters.at(indices.at(i));
        ModuleHit shared_hit(hit);
        shared_hit.energy *= frac[i]/total_frac;
        cluster.AddHit(shared_hit);
    }

    return true;
}

#endif


//============================================================================//
// PRad Cluster Reconstruction Method                                         //
// Reconstruct the cluster within a square area                               //
//                                                                            //
// Weizhi Xiong, Chao Peng                                                    //
// 06/10/2016                                                                 //
//============================================================================//

#include "PRadSquareCluster.h"
#include "PRadHyCalReconstructor.h"
#include "PRadHyCalModule.h"
#include <algorithm>
#include <cmath>



PRadSquareCluster::PRadSquareCluster(PRadHyCalReconstructor *r)
: rec(r)
{
    // place holder
}

PRadSquareCluster::~PRadSquareCluster()
{
    // place holder
}

PRadHyCalCluster *PRadSquareCluster::Clone()
const
{
    return new PRadSquareCluster(*this);
}

inline bool PRadSquareCluster::checkBelongs(const ModuleHit &center,
                                            const ModuleHit &hit,
                                            float factor)
const
{
    float dist_x = factor*center->GetSizeX();
    float dist_y = factor*center->GetSizeY();

    if((fabs(center->GetX() - hit->GetX()) > dist_x) ||
       (fabs(center->GetY() - hit->GetY()) > dist_y))
        return false;

    return true;
}

void PRadSquareCluster::FormCluster(std::vector<ModuleHit> &hs, std::vector<ModuleCluster> &cls)
const
{
    // clear container first
    cls.clear();

    // form clusters with high energy hit seed
    groupHits(hs, cls);
}

void PRadSquareCluster::groupHits(std::vector<ModuleHit> &hits,
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
        // not belongs to any cluster, and the energy is larger than center threshold
        if(!fillClusters(hit, clusters) && (hit.energy > rec->config.min_center_energy))
        {
            clusters.emplace_back(hit, hit->GetLayoutFlag());
            clusters.back().AddHit(hit);
        }
    }
}

bool PRadSquareCluster::fillClusters(ModuleHit &hit, std::vector<ModuleCluster> &c)
const
{
    std::vector<unsigned int> indices;
    indices.reserve(50);

    // check how many clusters the hit belongs to
    for(unsigned int i = 0; i < c.size(); ++i)
    {
        const auto &center = c.at(i).center;
        // within the square range
        if(checkBelongs(center, hit, float(rec->config.square_size)/2.)) {
            indices.push_back(i);
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
bool PRadSquareCluster::splitHit(ModuleHit &hit,
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
        frac[i] = rec->getProf(center, hit).frac * center.energy;
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


//============================================================================//
// PRad Cluster Reconstruction Method                                         //
// Reconstruct the cluster using island algorithm from PrimEx                 //
//                                                                            //
// Ilya Larin, adapted GAMS island algorithm for HyCal in PrimEx.             //
// Maxime Levillain & Weizhi Xiong, developed a new method based on PrimEx    //
//                                  method, has a better performance but a    //
//                                  slightly worse result, 08/01/2016         //
// Weizhi Xiong, adapted island algorithm for PRad, reduced the discretized   //
//               energy value from 10 MeV to 0.1 MeV, wrote a C++ wrapper for //
//               the fortran code. 09/28/2016                                 //
// Chao Peng, rewrote the whole fortran code into C++. 11/21/2016             //
//            improve the island clustering with DFS method, remove LX's      //
//            non-splitting method since it is no longer used, 03/17/2018     //
//============================================================================//

#include <cmath>
#include <algorithm>
#include <iostream>
#include "PRadIslandCluster.h"
#include "PRadHyCalReconstructor.h"
#include "PRadHyCalDetector.h"
#include "PRadADCChannel.h"
#include "PRadTDCChannel.h"



PRadIslandCluster::PRadIslandCluster(PRadHyCalReconstructor *r)
: rec(r)
{
    // place holder
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



//============================================================================//
// Method based on the code from I. Larin for PrimEx                          //
//============================================================================//

void PRadIslandCluster::FormCluster(std::vector<ModuleHit> &hs, std::vector<ModuleCluster> &cls)
const
{
    // clear container first
    cls.clear();

    std::vector<std::vector<ModuleHit*>> groups;

    // group adjacent hits
    groupHits(hs, groups);

    // try to split the group
    for(auto &group : groups)
    {
        splitCluster(group, cls);
    }
}


// recursive function for the DFS grouping
void dfs_group(std::vector<ModuleHit*> &container, int idx,
               std::vector<ModuleHit> &hits, std::vector<bool> &visits, bool corner)
{
    auto &hit = hits[idx];
    container.push_back(&hit);
    visits[idx] = true;
    for(size_t i = 0; i < hits.size(); ++i)
    {
        if(visits[i] || ! hit->IsNeighbor(hits[i].id, corner))
                continue;
        dfs_group(container, i, hits, visits, corner);
    }
}

// group adjacent hits into raw clusters
// list is suitable for the merge process, but it has poor performance while loop
// over all the elements, the test with real data prefers vector as container
void PRadIslandCluster::groupHits(std::vector<ModuleHit> &hits,
                                  std::vector<std::vector<ModuleHit*>> &groups)
const
{
    std::vector<bool> visits(hits.size(), false);
    for(size_t i = 0; i < hits.size(); ++i)
    {
        // already in a group
        if(visits[i])
            continue;

        // create a new group and reserve some space for the possible hits
        groups.emplace_back();
        groups.back().reserve(ISLAND_GROUP_RESERVE);

        // group all the possible hits
        dfs_group(groups.back(), i, hits, visits, rec->config.corner_conn);
    }

}

// split one group into several clusters
void PRadIslandCluster::splitCluster(const std::vector<ModuleHit*> &group,
                                     std::vector<ModuleCluster> &clusters)
const
{
    // find local maximum
    auto maxima = findMaximums(group);

    // no cluster center found
    if(maxima.empty())
        return;

    if((maxima.size() == 1) ||                    // only 1 cluster
       (group.size() >= SPLIT_MAX_HITS) ||          // too many hits
       (maxima.size() >= SPLIT_MAX_MAXIMA)) {     // too many local maxima
        // create cluster based on the center
        clusters.emplace_back(*maxima.front(), (*maxima.front())->GetLayoutFlag());
        auto &cluster = clusters.back();

        for(auto &hit : group)
            cluster.AddHit(*hit);
    // split hits between several maxima
    } else {
        splitHits(maxima, group, clusters);
    }
}

// find local maxima in a group of adjacent hits
std::vector<ModuleHit*> PRadIslandCluster::findMaximums(const std::vector<ModuleHit*> &hits)
const
{
    std::vector<ModuleHit*> local_max;
    local_max.reserve(20);
    for(auto it = hits.begin(); it != hits.end(); ++it)
    {
        auto &hit1 = *it;
        // not a qualified center
        if(hit1->energy < rec->config.min_center_energy)
            continue;

        bool maximum = true;
        for(auto &hit2 : hits)
        {
            if(hit1 == hit2)
                continue;

            // we count corner in
            if((*hit1)->IsNeighbor(hit2->id, true) &&
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

// split hits between several local maxima inside a cluster group
void PRadIslandCluster::splitHits(const std::vector<ModuleHit*> &maxima,
                                  const std::vector<ModuleHit*> &hits,
                                  std::vector<ModuleCluster> &clusters)
const
{
    static SplitContainer split;

    // initialize fractions
    for(size_t i = 0; i < maxima.size(); ++i)
    {
        auto &center = *maxima.at(i);
        for(size_t j = 0; j < hits.size(); ++j)
        {
            auto &hit = *hits.at(j);
            split.frac[j][i] = rec->getProf(center, hit).frac*center.energy;
        }
    }

    // do iteration to evaluate the share of hits between several maxima
    evalFraction(hits, maxima, split);

    // done iteration, add cluster according to the final share of energy
    for(size_t i = 0; i < maxima.size(); ++i)
    {
        clusters.emplace_back(*maxima[i], (*maxima[i])->GetLayoutFlag());
        auto &cluster = clusters.back();

        for(size_t j = 0; j < hits.size(); ++j)
        {
            if(split.frac[j][i] == 0.)
                continue;

            // too small share, treat as zero
            if(split.norm_frac(i, j) < rec->config.least_split) {
                split.total[j] -= split.frac[j][i];
                continue;
            }

            ModuleHit new_hit(*hits.at(j));
            new_hit.energy *= split.norm_frac(i, j);
            cluster.AddHit(new_hit);

            // update the center energy
            if(new_hit == cluster.center)
                cluster.center.energy = new_hit.energy;

            // set flag to mark the splitted clusters
            SET_BIT(cluster.flag, kSplit);
        }
    }
}

inline void PRadIslandCluster::evalFraction(const std::vector<ModuleHit*> &hits,
                                            const std::vector<ModuleHit*> &maxima,
                                            SplitContainer &split)
const
{
    // temp containers for reconstruction
    BaseHit temp[POS_RECON_HITS];

    // iterations to refine the split energies
    size_t iters = rec->config.split_iter;
    while(iters-- > 0)
    {
        split.sum_frac(hits.size(), maxima.size());
        for(size_t i = 0; i < maxima.size(); ++i)
        {
            // cluster center reconstruction
            auto &center = *maxima.at(i);
            float tot_E = center.energy;
            int count = 0;
            for(size_t j = 0; j < hits.size(); ++j)
            {
                auto &hit = *hits.at(j);

                if(hit.id == center.id || split.frac[j][i] == 0.)
                    continue;

                // using 3x3 to reconstruct hit position
                double dx, dy;
                center->QuantizedDist(hit.ptr, dx, dy);
                if(std::abs(dx) < 1.01 && std::abs(dy) < 1.01) {
                    temp[count].x = dx;
                    temp[count].y = dy;
                    temp[count].E = hit.energy*split.norm_frac(i, j);
                    tot_E += temp[count].E;
                    count++;
                }
            }

            BaseHit recon;
            rec->reconstructPos(center, temp, count, &recon);

            // update profile with the reconstructed center
            for(size_t j = 0; j < hits.size(); ++j)
            {
                auto &hit = *hits.at(j);
                split.frac[j][i] = rec->getProf(recon, hit).frac*tot_E;
            }
        }
    }
    split.sum_frac(hits.size(), maxima.size());
}


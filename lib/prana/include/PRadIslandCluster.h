#ifndef PRAD_ISLAND_CLUSTER_H
#define PRAD_ISLAND_CLUSTER_H

#include <vector>
#include "PRadHyCalCluster.h"



// reserve some space for grouping hits
#define ISLAND_GROUP_RESERVE 50
// do not split the group if the number of hits in it exceed this number
#define SPLIT_MAX_HITS 100
// do not split the group if the number of maxima in it exceed this number
#define SPLIT_MAX_MAXIMA 10

struct SplitContainer
{
    float frac[SPLIT_MAX_HITS][SPLIT_MAX_MAXIMA], total[SPLIT_MAX_HITS];

    // a helper function to sum 2d array
    void sum_frac(size_t hits, size_t maximums)
    {
        for(size_t i = 0; i < hits; ++i)
        {
            total[i] = 0;
            for(size_t j = 0; j < maximums; ++j)
                total[i] += frac[i][j];
        }
    }

    // get the normalized fraction
    inline float norm_frac(size_t i, size_t j)
    {
        return frac[j][i]/total[j];
    }
};

class PRadIslandCluster : public PRadHyCalCluster
{
public:
    PRadIslandCluster(class PRadHyCalReconstructor *r);
    virtual ~PRadIslandCluster();
    PRadHyCalCluster *Clone() const;

    void FormCluster(std::vector<ModuleHit> &hs, std::vector<ModuleCluster> &cls) const;

protected:
    void groupHits(std::vector<ModuleHit> &hits,
                   std::vector<std::vector<ModuleHit*>> &groups) const;
    bool fillClusters(ModuleHit &hit, std::vector<std::vector<ModuleHit*>> &groups) const;
    bool checkAdjacent(const std::vector<ModuleHit*> &g1, const std::vector<ModuleHit*> &g2) const;
    void splitCluster(const std::vector<ModuleHit*> &grp, std::vector<ModuleCluster> &c) const;
    std::vector<ModuleHit*> findMaximums(const std::vector<ModuleHit*> &g) const;
    void splitHits(const std::vector<ModuleHit*> &maximums,
                   const std::vector<ModuleHit*> &hits,
                   std::vector<ModuleCluster> &clusters) const;
    void evalFraction(const std::vector<ModuleHit*> &maximums,
                      const std::vector<ModuleHit*> &hits,
                      SplitContainer &split) const;

protected:
    class PRadHyCalReconstructor *rec;
};

#endif

#ifndef PRAD_ISLAND_CLUSTER_H
#define PRAD_ISLAND_CLUSTER_H

#include <vector>
#include "PRadHyCalCluster.h"

// use the original island method, do iterations on splitting clusters
// undefine it will switch the method to do a coarse splitting with an
// improved performance
// SUGGEST to be defined since performance should not be an issue for
// this adapted C++ island code
#define ISLAND_FINE_SPLIT
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

    inline float norm_frac(size_t i, size_t j)
    {
        return frac[j][i]/total[j];
    }
};

class PRadIslandCluster : public PRadHyCalCluster
{
public:
    PRadIslandCluster(const std::string &path = "");
    virtual ~PRadIslandCluster();
    PRadHyCalCluster *Clone() const;

    void Configure(const std::string &path);
    void FormCluster(std::vector<ModuleHit> &hits,
                     std::vector<ModuleCluster> &clusters) const;

protected:
// primex method, do iterations for splitting
#ifdef ISLAND_FINE_SPLIT
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

// M. Levillain and W. Xiong method, a quick but slightly rough splitting
#else
    void groupHits(std::vector<ModuleHit> &hits,
                   std::vector<ModuleCluster> &clusters) const;
    bool fillClusters(ModuleHit &hit, std::vector<ModuleCluster> &clusters) const;
    bool splitHit(ModuleHit &hit,
                  std::vector<ModuleCluster> &clusters,
                  std::vector<unsigned int> &indices) const;
#endif

protected:
    // parameters for reconstruction
    bool corner_conn;
    unsigned int split_iter;
    float least_share;
    std::vector<float> min_module_energy;
};

#endif

#ifndef PRAD_SQUARE_CLUSTER_H
#define PRAD_SQUARE_CLUSTER_H

#include <vector>
#include "PRadHyCalCluster.h"

class PRadSquareCluster : public PRadHyCalCluster
{
public:
    PRadSquareCluster();
    virtual ~PRadSquareCluster();
    PRadHyCalCluster *Clone() const;

    void FormCluster(PRadHyCalReconstructor *r);

protected:
    void groupHits(std::vector<ModuleHit> &hits,
                   std::vector<ModuleCluster> &clusters) const;
    bool fillClusters(ModuleHit &hit, std::vector<ModuleCluster> &clusters) const;
    bool splitHit(ModuleHit &hit,
                  std::vector<ModuleCluster> &clusters,
                  std::vector<unsigned int> &indices) const;
    bool checkBelongs(const ModuleHit &center, const ModuleHit &hit, float factor) const;

protected:
    PRadHyCalReconstructor *rec;
};

#endif

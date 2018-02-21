#ifndef PRAD_HYCAL_CLUSTER_H
#define PRAD_HYCAL_CLUSTER_H

#include <vector>
#include "PRadEventStruct.h"



class PRadHyCalCluster
{
public:
    virtual ~PRadHyCalCluster();
    virtual PRadHyCalCluster *Clone() const;

    virtual void FormCluster(std::vector<ModuleHit> &hs, std::vector<ModuleCluster> &cls) const;

protected:
    PRadHyCalCluster();
};

#endif

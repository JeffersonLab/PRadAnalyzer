#ifndef PRAD_HYCAL_CLUSTER_H
#define PRAD_HYCAL_CLUSTER_H

#include "PRadHyCalReconstructor.h"



class PRadHyCalCluster
{
public:
    virtual ~PRadHyCalCluster();
    virtual PRadHyCalCluster *Clone() const;

    virtual void FormCluster(PRadHyCalReconstructor *recon);

protected:
    PRadHyCalCluster();
};

#endif

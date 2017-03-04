#ifndef PRAD_HYCAL_CLUSTER_H
#define PRAD_HYCAL_CLUSTER_H

#include <vector>
#include <string>
#include <iostream>
#include "ConfigObject.h"
#include "PRadEventStruct.h"
#include "PRadHyCalDetector.h"

// we use 3x3 adjacent hits to reconstruct position
// here gives a larger volume to save the information
#define POS_RECON_HITS 15

class PRadHyCalCluster : public ConfigObject
{
public:
    virtual ~PRadHyCalCluster();
    virtual PRadHyCalCluster *Clone();
    virtual void Configure(const std::string &path);
    virtual void FormCluster(std::vector<ModuleHit> &hits,
                             std::vector<ModuleCluster> &clusters) const;
    virtual bool CheckCluster(const ModuleCluster &hit) const;
    virtual void LeakCorr(ModuleCluster &cluster, const std::vector<ModuleHit> &dead) const;

    void ReadVModuleList(const std::string &path);
    float GetWeight(const float &E, const float &E0) const;
    float GetShowerDepth(int module_type, const float &E) const;
    void AddVirtHits(ModuleCluster &cluster, const std::vector<ModuleHit> &dead) const;
    void CorrectVirtHits(ModuleCluster &cluster, float x, float y) const;

    HyCalHit Reconstruct(const ModuleCluster &cluster, const float &alpE = 1.) const;

protected:
    PRadHyCalCluster();
    int fillHits(BaseHit *temp,
                 int max_hits,
                 const ModuleHit &center,
                 const std::vector<ModuleHit> &hits) const;
    void reconstructPos(BaseHit *temp, int count, BaseHit *recon) const;

protected:
    bool depth_corr;
    bool leak_corr;
    bool linear_corr;
    float log_weight_thres;
    float min_cluster_energy;
    float min_center_energy;
    float least_leak;
    float linear_corr_limit;
    unsigned int min_cluster_size;
    unsigned int leak_iters;
    std::vector<ModuleHit> inner_virtual;
    std::vector<ModuleHit> outer_virtual;
};

#endif

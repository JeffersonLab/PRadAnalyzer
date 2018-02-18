#ifndef PRAD_HYCAL_CLUSTER_H
#define PRAD_HYCAL_CLUSTER_H

#include <vector>
#include <string>
#include <iostream>
#include "ConfigObject.h"
#include "PRadEventStruct.h"
#include "PRadHyCalDetector.h"
#include "PRadClusterProfile.h"



// we use 3x3 adjacent hits to reconstruct position
// here gives a larger volume to save the information
#define POS_RECON_HITS 15

class PRadHyCalCluster : public ConfigObject
{
public:
    virtual ~PRadHyCalCluster();
    virtual PRadHyCalCluster *Clone() const;
    virtual void Configure(const std::string &path);
    virtual void CollectHits(PRadHyCalDetector *det);
    virtual void Reconstruct(PRadHyCalDetector *det, PRadClusterProfile *prof);
    virtual void FormCluster(std::vector<ModuleHit> &hits,
                             std::vector<ModuleCluster> &clusters) const;
    virtual bool CheckCluster(const ModuleCluster &hit) const;
    virtual void LeakCorr(ModuleCluster &cluster) const;

    inline void AddHit(const ModuleHit &hit) {module_hits.emplace_back(hit);}
    inline void AddHit(ModuleHit &&hit) {module_hits.emplace_back(hit);}
    inline void ClearHits() {module_hits.clear();}
    const std::vector<ModuleHit> &GetHits() const {return module_hits;}
    const std::vector<ModuleCluster> &GetClusters() const {return module_clusters;}

    void ReadVModuleList(const std::string &path);
    float GetWeight(const float &E, const float &E0) const;
    float GetShowerDepth(int module_type, const float &E) const;
    void CorrectVirtHits(BaseHit &hit, std::vector<ModuleHit> &vhits,
                         const ModuleCluster &cluster) const;

    HyCalHit ReconstructHit(const ModuleCluster &cluster, const float &alpE = 1.) const;

protected:
    PRadHyCalCluster();
    int fillHits(BaseHit *temp,
                 int max_hits,
                 const ModuleHit &center,
                 const std::vector<ModuleHit> &hits) const;
    int reconstructPos(const ModuleHit &center, BaseHit *temp, int count, BaseHit *hit) const;
    int reconstructPos(const ModuleCluster &cluster, BaseHit *hit) const;
    PRadClusterProfile::Value getProf(double cx, double cy, double cE, const ModuleHit &hit) const;
    PRadClusterProfile::Value getProf(const BaseHit &center, const ModuleHit &hit) const;
    PRadClusterProfile::Value getProf(const ModuleHit &center, const ModuleHit &hit) const;
    double evalCluster(const BaseHit &center, const ModuleCluster &cluster) const;
    inline double hitDistance(const ModuleHit &m1, const ModuleHit &m2)
    const
    {
        return detector->QuantizedDist(m1.ptr, m2.ptr);
    }

protected:
    PRadHyCalDetector *detector;
    PRadClusterProfile *profile;
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
    std::vector<ModuleHit> module_hits;
    std::vector<ModuleCluster> module_clusters;
};

#endif

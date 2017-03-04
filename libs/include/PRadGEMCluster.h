#ifndef PRAD_GEM_CLUSTER_H
#define PRAD_GEM_CLUSTER_H

#include <string>
#include <unordered_map>
#include "PRadGEMPlane.h"
#include "PRadEventStruct.h"
#include "ConfigObject.h"

class PRadGEMDetector;

class PRadGEMCluster : public ConfigObject
{
public:
    PRadGEMCluster(const std::string &c_path = "");
    virtual ~PRadGEMCluster();

    // functions that to be overloaded
    void Configure(const std::string &path = "");

    void FormClusters(std::vector<StripHit> &hits,
                      std::vector<StripCluster> &clusters) const;
    void CartesianReconstruct(const std::vector<StripCluster> &x_cluster,
                              const std::vector<StripCluster> &y_cluster,
                              std::vector<GEMHit> &container,
                              int det_id) const;

protected:
    void groupHits(std::vector<StripHit> &h, std::vector<StripCluster> &c) const;
    void splitCluster(std::vector<StripCluster> &c) const;
    bool splitCluster_sub(StripCluster &c, StripCluster &c1) const;
    void filterCluster(std::vector<StripCluster> &c) const;
    bool filterCrossTalk(const StripCluster &cluster,
                         const std::vector<StripCluster> &clusters) const;
    void reconstructCluster(std::vector<StripCluster> &c) const;

protected:
    // parameters
    unsigned int min_cluster_hits;
    unsigned int max_cluster_hits;
    float split_cluster_diff;
    float cross_talk_width;

    // cross talk characteristic distances
    std::vector<double> charac_distance;
};

#endif

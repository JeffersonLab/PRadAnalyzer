#ifndef PRAD_HYCAL_RECONSTRUCTOR_H
#define PRAD_HYCAL_RECONSTRUCTOR_H

#include <vector>
#include <string>
#include "PRadEventStruct.h"
#include "PRadClusterProfile.h"
#include "ConfigParser.h"
#include "ConfigObject.h"


// we use 3x3 adjacent hits to reconstruct position
// here gives a larger volume to save the information
#define POS_RECON_HITS 15

class PRadHyCalDetector;

// class to manage the various hycal clustering methods
class PRadHyCalReconstructor : public ConfigObject
{
public:
    friend class PRadIslandCluster;
    friend class PRadSquareCluster;
    // clustering method enum
    enum MethodEnum
    {
        Undefined = -1,
        Island = 0,
        Square,
        Max_Methods,
    };
    // macro in ConfigParser.h
    ENUM_MAP(MethodEnum, "Island|Square");

public:
    PRadHyCalReconstructor(const std::string &conf_path = "");

    // copy/move constructors
    PRadHyCalReconstructor(const PRadHyCalReconstructor &that);
    PRadHyCalReconstructor(PRadHyCalReconstructor &&that);

    // destructor
    virtual ~PRadHyCalReconstructor();

    // copy/move assignment operators
    PRadHyCalReconstructor &operator =(const PRadHyCalReconstructor &rhs);
    PRadHyCalReconstructor &operator =(PRadHyCalReconstructor &&rhs);

    // configuration
    void Configure(const std::string &path);

    // core functions
    void Reconstruct(PRadHyCalDetector *det);
    void Reconstruct(PRadHyCalDetector *det, const EventData &event);
    void CollectHits(PRadHyCalDetector *det);
    void CollectHits(PRadHyCalDetector *det, const EventData &event);
    void ReconstructHits(PRadHyCalDetector *det);
    void AddTiming(PRadHyCalDetector *det);
    void AddTiming(PRadHyCalDetector *det, const EventData &event);

    // help functions
    HyCalHit Cluster2Hit(const ModuleCluster &cl) const;
    void LeakCorr(ModuleCluster &cluster) const;
    void CorrectVirtHits(BaseHit &hit, std::vector<ModuleHit> &vhits,
                         const ModuleCluster &cluster) const;
    bool CheckCluster(const ModuleCluster &cluster) const;
    double EvalCluster(const BaseHit &c, const ModuleCluster &cl) const;

    // profile related
    PRadClusterProfile *GetProfile() {return &profile;}
    void LoadProfile(int t, const std::string &path) {profile.Load(t, path);}

    // methods information
    bool SetMethod(const std::string &name);
    bool SetMethod(MethodEnum newtype);
    class PRadHyCalCluster *GetMethod() const {return method;}
    MethodEnum GetMethodType() const {return type;}
    std::string GetMethodName() const {return MethodEnum2str(type);}
    std::vector<std::string> GetMethodNames() const;

    // containers
    const std::vector<ModuleHit> &GetHits() const {return module_hits;}
    const std::vector<ModuleCluster> &GetClusters() const {return module_clusters;}

protected:
    float getWeight(const float &E, const float &E0) const;
    float getShowerDepth(int module_type, const float &E) const;
    int reconstructPos(const ModuleHit &center, BaseHit *temp, int count, BaseHit *hit) const;
    int reconstructPos(const ModuleCluster &cl, BaseHit *hit) const;
    int fillHits(BaseHit *temp, int max_hits, const ModuleHit &center,
                 const std::vector<ModuleHit> &hits) const;
    PRadClusterProfile::Value getProf(const ModuleHit &c, const ModuleHit &hit) const;
    PRadClusterProfile::Value getProf(double cx, double cy, double cE, const ModuleHit &hit) const;
    PRadClusterProfile::Value getProf(const BaseHit &c, const ModuleHit &hit) const;


private:
    MethodEnum type;
    class PRadHyCalCluster *method;
    PRadClusterProfile profile;

    std::vector<ModuleHit> module_hits;
    std::vector<ModuleCluster> module_clusters;

    bool depth_corr;
    bool leak_corr;
    bool linear_corr;
    bool corner_conn;
    float log_weight_thres;
    float min_cluster_energy;
    float min_center_energy;
    float least_leak;
    float least_share;
    float linear_corr_limit;
    unsigned int min_cluster_size;
    unsigned int leak_iters;
    unsigned int split_iter;
    unsigned int square_size;
    std::vector<float> min_module_energy;
};

#endif // PRAD_HYCAL_RECONSTRUCTOR

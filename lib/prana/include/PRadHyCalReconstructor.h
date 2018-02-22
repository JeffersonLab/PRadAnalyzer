#ifndef PRAD_HYCAL_RECONSTRUCTOR_H
#define PRAD_HYCAL_RECONSTRUCTOR_H

#include <vector>
#include <string>
#include "PRadEventStruct.h"
#include "PRadClusterProfile.h"
#include "PRadClusterDensity.h"
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
    enum ClMethod
    {
        Undefined_ClMethod = -1,
        Island = 0,
        Square,
        Max_ClMethods,
    };
    // macro in ConfigParser.h
    ENUM_MAP(ClMethod, 0, "Island|Square");

    // position reconstruction method
    enum PosMethod
    {
        Undefined_PosMethod = -1,
        Logarithmic = 0,
        Linear,
        Max_PosMethods,
    };
    ENUM_MAP(PosMethod, 0, "Logarithmic|Linear");

    // configuraiton data package
    struct Config
    {
        // general
        bool depth_corr, leak_corr, linear_corr, den_corr, sene_corr;
        float log_weight_thres, min_cluster_energy, min_center_energy;
        float least_leak, linear_corr_limit;
        unsigned int min_cluster_size, leak_iters;
        std::vector<float> min_module_energy;

        // for island
        bool corner_conn;
        unsigned int split_iter;
        float least_split;
        // for square
        unsigned int square_size;
    };

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
    PRadClusterDensity *GetDensityParams() {return &density;}
    void LoadDensityParams(int t, const std::string &p_path, const std::string &e_path)
    {density.Load(t, p_path, e_path);}
    void ChooseDensitySet(PRadClusterDensity::SetEnum i) {density.ChooseSet(i);}

    // methods information
    bool SetClusterMethod(const std::string &name);
    bool SetClusterMethod(ClMethod newtype);
    class PRadHyCalCluster *GetClusterMethod() const {return cluster;}
    ClMethod GetClusterMethodType() const {return cltype;}
    std::string GetClusterMethodName() const {return ClMethod2str(cltype);}
    std::vector<std::string> GetClusterMethodNames() const;

    bool SetPositionMethod(const std::string &name);
    bool SetPositionMethod(PosMethod newtype);
    PosMethod GetPositionMethodType() const {return postype;}
    std::string GetPositionMethodName() const {return PosMethod2str(postype);}
    std::vector<std::string> GetPositionMethodNames() const;


    // containers
    const std::vector<ModuleHit> &GetHits() const {return module_hits;}
    const std::vector<ModuleCluster> &GetClusters() const {return module_clusters;}

protected:
    float getWeight(const float &E, const float &E0) const;
    float getPosBias(const std::vector<float> &pars, const float &dx) const;
    float getShowerDepth(int module_type, const float &E) const;
    int reconstructPos(const ModuleHit &center, BaseHit *temp, int count, BaseHit *hit) const;
    int reconstructPos(const ModuleCluster &cl, BaseHit *hit) const;
    int fillHits(BaseHit *temp, int max_hits, const ModuleHit &center,
                 const std::vector<ModuleHit> &hits) const;
    PRadClusterProfile::Value getProf(const ModuleHit &c, const ModuleHit &hit) const;
    PRadClusterProfile::Value getProf(double cx, double cy, double cE, const ModuleHit &hit) const;
    PRadClusterProfile::Value getProf(const BaseHit &c, const ModuleHit &hit) const;


private:
    PRadClusterProfile profile;
    PRadClusterDensity density;
    ClMethod cltype;
    class PRadHyCalCluster *cluster;
    PosMethod postype;
    Config config;

    std::vector<ModuleHit> module_hits;
    std::vector<ModuleCluster> module_clusters;
};

#endif // PRAD_HYCAL_RECONSTRUCTOR

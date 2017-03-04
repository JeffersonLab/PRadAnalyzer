#ifndef PRAD_CLUSTER_PROFILE_H
#define PRAD_CLUSTER_PROFILE_H

#include <string>
#include "PRadHyCalDetector.h"
#include "PRadEventStruct.h"

class PRadClusterProfile
{
public:
    struct Profile
    {
        float frac;
        float err;

        Profile() : frac(0), err(0) {};
        Profile(float f, float e) : frac(f), err(e) {};
    };

public:
    static PRadClusterProfile &Instance()
    {
        static PRadClusterProfile instance;

        return instance;
    }

    // copy/move constructors
    PRadClusterProfile(const PRadClusterProfile &that) = delete;
    PRadClusterProfile(PRadClusterProfile &&that) = delete;

    virtual ~PRadClusterProfile();

    // copy/move assignment operators
    PRadClusterProfile &operator =(const PRadClusterProfile &rhs) = delete;
    PRadClusterProfile &operator =(PRadClusterProfile &&rhs) = delete;

    void Resize(int type, int xsize, int ysize);
    void Clear();
    void LoadProfile(int type, const std::string &path);
    const Profile &GetProfile(int type, int x, int y) const;
    const Profile &GetProfile(const ModuleHit &m1, const ModuleHit &m2) const;
    const Profile &GetProfile(const float &x, const float &y, const ModuleHit &hit) const;
    float EvalEstimator(const BaseHit &hit, const ModuleCluster &cluster) const;

private:
    PRadClusterProfile(int type = 2, int xsize = 501, int ysize = 501);
    void reserve();
    void release();

private:
    int types;
    int x_steps;
    int y_steps;
    Profile ***profiles;
    Profile empty_prof;
    Profile trans_prof;
};

#endif

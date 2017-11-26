#ifndef PRAD_CLUSTER_PROFILE_H
#define PRAD_CLUSTER_PROFILE_H

#include <vector>
#include <string>
#include "PRadEventStruct.h"

#define CLPROF_MIN_ENE 200      // max energy in the profile, MeV
#define CLPROF_MAX_ENE 2100     // min energy in the profile, MeV
#define CLPROF_MAX_DIST 5       // max distance in the profile, module size
#define CLPROF_STEP_ENE 100     // each step in energy
#define CLPROF_STEP_DIST 0.001  // each step in distance

class PRadHyCalDetector;

class PRadClusterProfile
{
public:
    struct Value
    {
        double frac, err;

        Value(double f = 0., double e = 0.) : frac(f), err(e) {}
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

    void LoadProfile(int type, const std::string &path);
    Value GetProfile(int type, double dist, double energy) const;

private:
    PRadClusterProfile();

private:
    std::vector<std::vector<std::vector<Value>>> profiles;
};

#endif

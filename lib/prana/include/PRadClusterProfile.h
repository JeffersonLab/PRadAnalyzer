#ifndef PRAD_CLUSTER_PROFILE_H
#define PRAD_CLUSTER_PROFILE_H

#include <vector>
#include <string>
#include "PRadEventStruct.h"


class PRadHyCalDetector;

class PRadClusterProfile
{
public:
    struct Value
    {
        double frac, err;

        Value(double f = 0., double e = 0.) : frac(f), err(e) {}
    };

    struct Profile
    {
        double min_ene, max_ene, step_ene;
        double max_dist, step_dist;
        std::vector<std::vector<Value>> values;

        void Resize(int Ne, int Nd)
        {
            // energy grids
            values.resize(Ne);
            // distance grids
            for(auto &e_prof : values)
                e_prof.resize(Nd);
        }

        typedef std::vector<std::vector<Value>>::size_type size_type;
        inline std::vector<Value> &operator [] (size_type i) {return values[i];}
        inline const std::vector<Value> &operator [] (size_type i) const {return values[i];}
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
    std::vector<Profile> profiles;
};

#endif

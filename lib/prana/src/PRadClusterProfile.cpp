//============================================================================//
// A class to store the information about the HyCal cluster profile           //
// It is a singleton class to be shared among different clustering methods    //
//                                                                            //
// Chao Peng                                                                  //
// 11/22/2016                                                                 //
//                                                                            //
// UPDATE - Chao Peng, Xinzhan Bai, 11/25/2017                                //
// Cluster profile is replaced with the one from GEANT4 simulation PRadSim    //
// Raw profile from simulation is obtained by Xinzhan Bai                     //
// The raw profile was then smoothed by a general fit based on neural network //
//============================================================================//

#include "PRadClusterProfile.h"
#include "PRadHyCalModule.h"
#include "ConfigParser.h"
#include <cmath>



PRadClusterProfile::PRadClusterProfile()
{
    profiles.resize(static_cast<size_t>(PRadHyCalModule::Max_Types));
}

PRadClusterProfile::~PRadClusterProfile()
{
    // place holder
}

void PRadClusterProfile::LoadProfile(int type, const std::string &path)
{
    if((size_t)type >=  profiles.size()) {
        std::cerr << "PRad Cluster Profile Error: Exceed current capacity, "
                  << "only has " << profiles.size() << " types."
                  << std::endl;
        return;
    }

    ConfigParser parser;
    if(!parser.OpenFile(path)) {
        std::cerr << "PRad Cluster Profile Error: File"
                  << " \"" << path << "\" "
                  << "cannot be opened."
                  << std::endl;
        return;
    }

    auto &profile = profiles[type];
    // read configurations
    if(parser.ParseLine() && parser.NbofElements() == 5) {
        parser >> profile.min_ene >> profile.max_ene >> profile.step_ene
               >> profile.max_dist >> profile.step_dist;
    } else {
        std::cerr << "PRad Cluster Profile Error: Cannot find correct configuration "
                  << "from profile " << " \"" << path << "\", "
                  << std::endl;
    }

    // prepare enough space
    int Ne = (profile.max_ene - profile.min_ene)/profile.step_ene + 1;
    int Nd = profile.max_dist/profile.step_dist + 1;
    profile.Resize(Ne, Nd);

    // read profiles
    int ie, id;
    double val, err;
    while(parser.ParseLine())
    {
        if(!parser.CheckElements(4))
            continue;

        parser >> ie >> id >> val >> err;

        if(id >= 0 && id < Nd && ie >= 0 && ie < Ne) {
            profile[ie][id] = Value(val, err);
        } else {
            std::cout << "PRad Cluster Profile Warning: "
                      << "Step (" << ie << ", " << id << ") is out of range"
                      << std::endl;
        }
    }
}

// get the profile value by distance and energy
PRadClusterProfile::Value PRadClusterProfile::GetProfile(int type, double dist, double energy)
const
{
    if((size_t)type >=  profiles.size())
        return Value();

    auto &profile = profiles[type];

    // out of range, should be 0
    if(dist >= profile.max_dist)
        return Value();

    // just round the step value for distance since the step size is small enough
    int id = int(dist/profile.step_dist + 0.5);

    double norm_e = (energy - profile.min_ene)/profile.step_ene;
    int ie = int(norm_e);

    // lower than the limit
    if(ie < 0) {
        return profile.values.front().at(id);
    // higher than the limit
    } else if(ie + 1 >= (int)profile.values.size()) {
        return profile.values.back().at(id);
    }

    // no need to interpolate
    double res = norm_e - ie, res2 = 1. - res;
    if(res < 0.05)
        return profile[ie][id];
    if(res2 < 0.05)
        return profile[ie + 1][id];

    // interpolation for energy
    auto &val1 = profile[ie][id];
    auto &val2 = profile[ie + 1][id];

    return Value(res*val2.frac + res2*val1.frac, res*val2.err + res2*val1.err);
}

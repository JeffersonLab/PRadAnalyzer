#ifndef PRAD_CLUSTER_DENSITY_H
#define PRAD_CLUSTER_DENSITY_H

#include <vector>
#include <string>
#include <unordered_map>
#include "PRadEventStruct.h"
#include "ConfigParser.h"



#define NB_ECORR_PARS 8


class PRadClusterDensity
{
public:
    enum SetEnum
    {
        Set_1GeV = 0,
        Set_2GeV,
        Max_SetEnums,
    };
    // macro in ConfigParser.h
    ENUM_MAP(SetEnum, 0, "Set_1GeV|Set_2GeV");

    struct Params
    {
        std::vector<float> x, y;
    };

    struct ParamsSet
    {
        float beam_energy;
        std::vector<float> energy_range;
        std::vector<Params> ppars;
        std::unordered_map<int, Params> epars_ee, epars_ep;

        void Resize(int ng, int ne, int np)
        {
            ppars.resize(ng*ne);
            for(auto &par : ppars)
            {
                par.x.resize(np, 0.);
                par.y.resize(np, 0.);
            }
        }
    };

public:
    PRadClusterDensity();
    virtual ~PRadClusterDensity();

    bool Load(int is, const std::string &p_path, const std::string &e_path);
    void CorrectBias(const ModuleHit &ctr, HyCalHit &hit, bool pos, bool ene) const;
    float GetPosBias(const std::vector<float> &pars, const float &dx) const;
    float GetEneBias(const std::vector<float> &pars, const float &dx,
                     const float &dy, const float &E0) const;

    void ChooseSet(SetEnum iset) {cur_set = iset;}

private:
    int getEnergyIndex(float energy, float theta, float res) const;
    int getGeometryIndex(const ModuleHit &center) const;
    bool processPosPars(ConfigParser &c_parser, ParamsSet &pset);
    bool processEnePars(ConfigParser &c_parser, ParamsSet &pset);


private:
    SetEnum cur_set;
    ParamsSet psets[Max_SetEnums];
};

#endif // PRAD_CLUSTER_DENSITY_H

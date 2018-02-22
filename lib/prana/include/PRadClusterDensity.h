#ifndef PRAD_CLUSTER_DENSITY_H
#define PRAD_CLUSTER_DENSITY_H

#include <vector>
#include <string>
#include <iostream>
#include <map>
#include "PRadEventStruct.h"
#include "ConfigParser.h"



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
        std::vector<Params> pars;

        void Resize(int ng, int ne, int np)
        {
            pars.resize(ng*ne);
            for(auto &par : pars)
            {
                par.x.resize(np, 0.);
                par.y.resize(np, 0.);
            }
        }
    };

public:
    PRadClusterDensity();
    virtual ~PRadClusterDensity();

    bool Load(int is, const std::string &path, float wthres = 0.05);
    bool CorrectBias(const ModuleHit &ctr, BaseHit &hit, float energy) const;
    float GetPosBias(const std::vector<float> &pars, const float &dx) const;

    void ChooseSet(SetEnum iset) {cur_set = iset;}

private:
    int getEnergyIndex(float energy, float theta, float res) const;
    int getGeometryIndex(const ModuleHit &center) const;


private:
    SetEnum cur_set;
    ParamsSet psets[Max_SetEnums];
};

#endif // PRAD_CLUSTER_DENSITY_H

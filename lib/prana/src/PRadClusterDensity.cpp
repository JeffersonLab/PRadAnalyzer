//============================================================================//
// A container class for PRad HyCal cluster density correction                //
// This correction eliminates the S shape bias in the reconstructed  position //
// Method based on                                                            //
// A. Lednev, NIM Physics Research A 366 (1995) 292-297                       //
//                                                                            //
// Maxime Lavillain, developer                                                //
// Chao Peng, integrated it into PRadHyCalReconstructor                       //
// 02/21/2018                                                                 //
//============================================================================//

#include "PRadClusterDensity.h"
#include "PRadHyCalModule.h"
#include "PRadHyCalDetector.h"
#include "PRadCoordSystem.h"
#include "ConfigObject.h"
#include "canalib.h"


// constructor
PRadClusterDensity::PRadClusterDensity()
: cur_set(Set_1GeV)
{
    // place holder
}

// destructor
PRadClusterDensity::~PRadClusterDensity()
{
    // place holder
}

bool PRadClusterDensity::Load(int is, const std::string &path, float wthres)
{
    ConfigParser c_parser;
    if(!c_parser.ReadFile(path)) {
        std::cerr << "PRad Cluster Density Error: Cannot open data file "
                  << "\"" << path << "\"." << std::endl;
        return false;
    }

    if(is < 0 || is >= static_cast<int>(Max_SetEnums)) {
        std::cerr << "PRad Cluster Density Error: Invalid parameter set enum = " << is
                  << ", skip reading parameters." << std::endl;
        return false;
    }

    auto &pset = psets[is];

    // configuration object
    ConfigObject conf_obj;
    std::string label;

    // read configuration
    while(c_parser.ParseLine())
    {
        // name
        c_parser >> label;
        // param start label
        if(label == "params")
            break;

        // configuration strings before the param label
        conf_obj.ReadConfigString(c_parser.CurrentLine());
    }

    // configurate pset
    pset.beam_energy = conf_obj.GetConfig<float>("Beam Energy");
    auto rangestr = conf_obj.GetConfig<std::string>("Energy Range");
    pset.energy_range = ConfigParser::stofs(rangestr, ",", " \t");
    int n_pars = conf_obj.GetConfig<int>("Number Of Params");
    int n_geo = conf_obj.GetConfig<int>("Geometry Groups");
    int n_ene = conf_obj.GetConfig<int>("Energy Groups");

    // sanity check, (range size - 1) + ep + Moller
    if(n_ene != static_cast<int>(pset.energy_range.size()) + 1) {
        std::cerr << "PRad HyCal Density Error: Unmatched energy range setting and "
                  << "number of energy groups in file "
                  << "\"" << path << "\"." << std::endl;
        return false;
    }

    pset.Resize(n_geo, n_ene, n_pars);

    float *pars;
    size_t index = 0, size = pset.pars.size();
    while(c_parser.ParseLine() && index < 2*size)
    {
        if(c_parser.NbofElements() < n_pars + 1)
            continue;

        c_parser >> label;

        if(index < size) {
            pars = &pset.pars[index].x[0];
        } else {
            pars = &pset.pars[index - size].y[0];
        }

        for(int i = 0; i < n_pars; ++i)
        {
            c_parser >> pars[i];
        }

        index ++;
    }

    return true;
}

// correct S shape position reconstruction
bool PRadClusterDensity::CorrectBias(const ModuleHit &ctr, BaseHit &hit, float energy)
const
{
    // geometrical index
    int ig = getGeometryIndex(ctr);
    float angle = PRadCoordSystem::GetPolarAngle(Point(hit.x, hit.y, PRadCoordSystem::hycal_z()));
    float res = ctr->GetDetector()->GetEneRes(ctr.ptr, energy);

    // energy index
    int ie = getEnergyIndex(energy, angle*cana::deg2rad, 6.*res);

    auto &pars = psets[static_cast<int>(cur_set)].pars;
    int idx = ie + ig*5;
    if(ie < 0 || ig < 0 || idx >= (int)pars.size()) {
        std::cerr << "PRad Cluster Density Error: Cannot find corresponding parameters." << std::endl;
        return false;
    }

    float dx = (hit.x - ctr->GetX())/ctr->GetSizeX();
    float dy = (hit.y - ctr->GetY())/ctr->GetSizeY();

    hit.x += GetPosBias(pars[idx].x, dx)*ctr->GetSizeX();
    hit.y += GetPosBias(pars[idx].y, dy)*ctr->GetSizeY();

    return true;
}

// get S shape bias in position reconstruction
float PRadClusterDensity::GetPosBias(const std::vector<float> &pars, const float &dx)
const
{
    // formula is
    // dx*(dx^2 - 0.25)*c0*(dx^4 + c1*dx^2 + c2)*(dx^2 - c3)
    float dx2 = dx*dx, dx4 = dx2*dx2;
    return dx*(dx2 - 0.25)*pars[0]*(dx4 + pars[1]*dx2 + pars[2])*(dx2 - pars[3]);
}

// helper function, determine Moller event energy
inline float moller_energy(float E, float theta)
{
    float a = (E - cana::ele_mass)/(E + cana::ele_mass);
    float cth = cos(theta);
    return cana::ele_mass*(1. + a*cth*cth)/(1. - a*cth*cth);
}

// helper function, determine ep elastic event energy
inline float mott_energy(float E, float theta)
{
    float sin2 = cana::pow2(sin(theta/2.));
    return E/(1. + 2.*E*sin2/cana::proton_mass);
}

// return energy index for searching parameters
// priority for ep or moller energy range
int PRadClusterDensity::getEnergyIndex(float energy, float theta, float res)
const
{
    auto &pset = psets[static_cast<int>(cur_set)];
    // ep energy
    float Eep = mott_energy(pset.beam_energy, theta);
    if(std::abs(energy/Eep - 1.) < res) return pset.energy_range.size() - 1;

    // moller energy
    float Eml = moller_energy(pset.beam_energy, theta);
    if(std::abs(energy/Eml - 1.) < res) return pset.energy_range.size();

    // energy range defined in the parameter data file
    for(int i = 0; i < (int)pset.energy_range.size() - 1; ++i)
    {
        if(energy > pset.energy_range[i] && energy <= pset.energy_range[i + 1])
            return i;
    }

    // cannot find corresponding range
    return -1;
}

// return geometrical index for searching parameters
// the density correction study was done by grouping several adjacent modules
// so this function finds the group number
int PRadClusterDensity::getGeometryIndex(const ModuleHit &center)
const
{
    if(center->GetType() == PRadHyCalModule::PbWO4) {
        int row = center->GetRow();
        int col = center->GetColumn();
        if(row < 17 && col < 17) return 32 + row/2*9 + col/2;
        if(row < 17) return 32 + row/2*9 + (33 - col)/2;
        if(col < 17) return 32 + (33 - row)/2*9 + col/2;
        return 32 + (33 - row)/2*9 + (33 - col)/2;
    } else {
        int row = (center->GetID() - 1)/30;
        int col = (center->GetID() - 1)%30;
        if(row < 6 && col < 24) return row/3*8 + col/3;
        if(row < 24 && col >= 24) return 16 + (col - 24)/3*8 + row/3;
        if(row >= 24 && col >= 6) return (29 - row)/3*8 + (29 - col)/3;
        if(row >= 6 && col < 6) return 16 + (5 - col)/3*8 + (29 - row)/3;
    }

    // error
    return -1;
}

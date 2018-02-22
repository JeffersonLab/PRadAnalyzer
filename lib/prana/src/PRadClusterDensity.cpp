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

// load parameters
bool PRadClusterDensity::Load(int is, const std::string &p_path, const std::string &e_path)
{
    if(is < 0 || is >= static_cast<int>(Max_SetEnums)) {
        std::cerr << "PRad Cluster Density Error: Invalid parameter set enum = " << is
                  << ", skip reading parameters." << std::endl;
        return false;
    }

    auto &pset = psets[is];

    // position density part
    ConfigParser c_parser;
    if(!c_parser.ReadFile(p_path)) {
        std::cerr << "PRad Cluster Density Error: Cannot open data file "
                  << "\"" << p_path << "\"." << std::endl;
        return false;
    }

    bool pos_success = processPosPars(c_parser, pset);

    // energy density part
    if(!c_parser.ReadFile(e_path)) {
        std::cerr << "PRad Cluster Density Error: Cannot open data file "
                  << "\"" << e_path << "\"." << std::endl;
        return false;
    }

    bool ene_success = processEnePars(c_parser, pset);

    return pos_success&ene_success;
}

// correct S shape position reconstruction
void PRadClusterDensity::CorrectBias(const ModuleHit &ctr, HyCalHit &hit, bool pos_corr, bool ene_corr)
const
{
    if(!pos_corr && !ene_corr) return;

    // required infromation
    auto &pset = psets[static_cast<int>(cur_set)];

    // energy index
    float angle = PRadCoordSystem::GetPolarAngle(Point(hit.x, hit.y, PRadCoordSystem::hycal_z() + hit.z));
    int ie = getEnergyIndex(hit.E, angle*cana::deg2rad, 6.*hit.sig_ene);

    // position correction
    if(pos_corr && !TEST_BIT(hit.flag, kDenCorr)) {
        // geometrical index
        int ig = getGeometryIndex(ctr);
        int idx = ie + ig*5;
        if(ie >= 0 && ig >= 0 && idx < (int) pset.ppars.size()) {
            float dx = (hit.x - ctr->GetX())/ctr->GetSizeX();
            float dy = (hit.y - ctr->GetY())/ctr->GetSizeY();
            auto &pars = pset.ppars[idx];
            hit.x += GetPosBias(pars.x, dx)*ctr->GetSizeX();
            hit.y += GetPosBias(pars.y, dy)*ctr->GetSizeY();
            SET_BIT(hit.flag, kDenCorr);
        }
    }

    // energy correction
    if(ene_corr && !TEST_BIT(hit.flag, kSEneCorr)) {
        int ep_set = pset.energy_range.size() - 1;
        int ee_set = ep_set + 1;
        // this will be affected by pos correction, it is intended
        float dx = (hit.x - ctr->GetX())/ctr->GetSizeX();
        float dy = (hit.y - ctr->GetY())/ctr->GetSizeY();

        // ep index
        if(ie == ep_set) {
            auto pit = pset.epars_ep.find(ctr.id);
            if(pit != pset.epars_ep.end()) {
                hit.E_Scorr = GetEneBias(pit->second.x, dx, dy, hit.E);
                SET_BIT(hit.flag, kSEneCorr);
                hit.E += hit.E_Scorr;
            }
        // ee index
        } else if(ie == ee_set) {
            auto pit = pset.epars_ee.find(ctr.id);
            if(pit != pset.epars_ee.end()) {
                hit.E_Scorr = GetEneBias(pit->second.x, dx, dy, hit.E);
                SET_BIT(hit.flag, kSEneCorr);
                hit.E += hit.E_Scorr;
            }
        }
    }

    return;
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

float PRadClusterDensity::GetEneBias(const std::vector<float> &pars, const float &dx,
                                     const float &dy, const float &E0)
const
{
    // formula is
    // E_0/c0 /(1 + c1*dx^2 + c2*dy^2 + c3*dx^2*dy^2 + c4*dx^4 + c5*dy^4 + c6*dx + c7 * dy)
    float dx2 = dx*dx, dx4 = dx2*dx2;
    float dy2 = dy*dy, dy4 = dy2*dy2;
    return E0/pars[0]/(1. + pars[1]*dx2 + pars[2]*dy2 + pars[3]*dx2*dy2
                       + pars[4]*dx4 + pars[5]*dy4 + pars[6]*dx + pars[7]*dy);
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

// helper function to process position density file
bool PRadClusterDensity::processPosPars(ConfigParser &c_parser, ParamsSet &pset)
{
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
                  << "number of energy groups in position parameter file."
                  << std::endl;
        return false;
    }

    pset.Resize(n_geo, n_ene, n_pars);

    // read parameters for position density
    float *pars;
    size_t index = 0, size = pset.ppars.size();
    while(c_parser.ParseLine() && index < 2*size)
    {
        c_parser >> label;

        if(index < size) {
            pars = &pset.ppars[index].x[0];
        } else {
            pars = &pset.ppars[index - size].y[0];
        }

        for(int i = 0; i < n_pars; ++i)
        {
            c_parser >> pars[i];
        }

        index ++;
    }

    return true;
}

// helper function to process energy density file
bool PRadClusterDensity::processEnePars(ConfigParser &c_parser, ParamsSet &pset)
{
    bool ep_or_ee = true;

    std::string name;
    Params par;
    // number of parameters required for energy correction
    par.x.resize(NB_ECORR_PARS, 0.);
    while(c_parser.ParseLine())
    {
        // check labels
        c_parser >> name;
        if(name == "EE_PARAMS") {
            ep_or_ee = false;
            continue;
        } if(name == "EP_PARAMS") {
            ep_or_ee = true;
            continue;
        }

        // read parameters
        for(int i = 0; i < NB_ECORR_PARS; ++i)
        {
            c_parser >> par.x[i];
        }

        // set map
        int id = PRadHyCalModule::name_to_id(name);

        if(ep_or_ee) {
            pset.epars_ep[id] = par;
        } else {
            pset.epars_ee[id] = par;
        }
    }

    return true;
}


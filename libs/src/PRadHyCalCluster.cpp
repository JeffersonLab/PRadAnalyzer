//============================================================================//
// Basic PRad Cluster Reconstruction Class For HyCal                          //
// Different reconstruction methods can be implemented accordingly            //
//                                                                            //
// Chao Peng, Weizhi Xiong                                                    //
// 09/28/2016                                                                 //
//============================================================================//

#include <cmath>
#include <iostream>
#include <iomanip>
#include "PRadHyCalCluster.h"
#include "PRadClusterProfile.h"


const PRadClusterProfile &__hc_prof = PRadClusterProfile::Instance();

PRadHyCalCluster::PRadHyCalCluster()
: depth_corr(true), leak_corr(true), linear_corr(true),
  log_weight_thres(3.6), min_cluster_energy(30.), min_center_energy(10.),
  least_leak(0.05), linear_corr_limit(0.6), min_cluster_size(1), leak_iters(3)
{
    // place holder
}

PRadHyCalCluster::~PRadHyCalCluster()
{
    // place holder
}

PRadHyCalCluster* PRadHyCalCluster::Clone()
{
    return new PRadHyCalCluster(*this);
}

void PRadHyCalCluster::Configure(const std::string &path)
{
    bool verbose = false;

    if(!path.empty()) {
        ConfigObject::Configure(path);
        verbose = true;
    }

    depth_corr = getDefConfig<bool>("Shower Depth Correction", true, verbose);
    leak_corr = getDefConfig<bool>("Leakage Correction", true, verbose);
    linear_corr = getDefConfig<bool>("Non Linearity Correction", true, verbose);
    log_weight_thres = getDefConfig<float>("Log Weight Threshold", 3.6, verbose);
    min_cluster_energy = getDefConfig<float>("Minimum Cluster Energy", 50., verbose);
    min_center_energy = getDefConfig<float>("Minimum Center Energy", 10., verbose);
    min_cluster_size = getDefConfig<unsigned int>("Minimum Cluster Size", 1, verbose);
    least_leak = getDefConfig<float>("Least Leakage Fraction", 0.05, verbose);
    leak_iters = getDefConfig<unsigned int>("Leakage Iterations", 3, verbose);
    linear_corr_limit = getDefConfig<float>("Non Linearity Limit", 0.6, verbose);

    ReadVModuleList(GetConfig<std::string>("Virtual Module List"));
}

void PRadHyCalCluster::ReadVModuleList(const std::string &path)
{
    if(path.empty())
        return;

    ConfigParser c_parser;
    c_parser.ReadFile(path);

    inner_virtual.clear();
    outer_virtual.clear();

    std::string name;
    std::string type, sector;
    PRadHyCalModule::Geometry geo;

    // some info that is not read from list
    while (c_parser.ParseLine())
    {
        if(!c_parser.CheckElements(8))
            continue;

        c_parser >> name >> type
                 >> geo.size_x >> geo.size_y >> geo.size_z
                 >> geo.x >> geo.y >> geo.z;

        geo.type = PRadHyCalModule::get_module_type(type.c_str());

        if(ConfigParser::str_upper(name) == "INNER") {
            ModuleHit inner_hit(false);
            inner_hit.id = -1;
            inner_hit.geo = geo;
            inner_hit.sector = 0;
            inner_virtual.push_back(inner_hit);
        } else if(ConfigParser::str_upper(name) == "OUTER") {
            ModuleHit outer_hit(false);
            outer_hit.id = -1;
            outer_hit.geo = geo;
            outer_hit.sector = 1;
            outer_virtual.push_back(outer_hit);
        }
    }
}

void PRadHyCalCluster::FormCluster(std::vector<ModuleHit> &,
                                   std::vector<ModuleCluster> &)
const
{
    // to be implemented by methods
}

float PRadHyCalCluster::GetWeight(const float &E, const float &E0)
const
{
    float w = log_weight_thres + log(E/E0);
    if(w < 0.)
        return 0.;
    return w;
}

// get shower depth, unit is in MeV
float PRadHyCalCluster::GetShowerDepth(int module_type, const float &E)
const
{
    if(depth_corr && E > 0.) {
        // here all the values are hard coded, because these are all physical
        // values corresponding to the material, so no need to change
        // it returns the maximum shower depth that
        // t = X0*(ln(E0/Ec) - Cf),
        // where X0 is radiation length, Ec is critical energy, Cf = -0.5 for
        // electron induced shower and 0.5 for photon
        // units are in mm and MeV
        if(module_type == PRadHyCalModule::PbWO4)
            return 8.6*(log(E/1.1) - 0.5);

        // -101.2 is the surface difference between Lead Glass and Lead Tungstate modules
        if(module_type == PRadHyCalModule::PbGlass)
            return 26.7*(log(E/2.84) - 0.5);
    }

    return 0.;
}

// check if the cluster is good enough
bool PRadHyCalCluster::CheckCluster(const ModuleCluster &cluster)
const
{
    if((cluster.energy < min_cluster_energy) ||
       (cluster.hits.size() < min_cluster_size))
            return false;

    return true;
}

// reconstruct cluster
HyCalHit PRadHyCalCluster::Reconstruct(const ModuleCluster &cluster, const float &alpE)
const
{
    // initialize the hit
    HyCalHit hycal_hit(cluster.center.id,       // center id
                       cluster.center.flag,     // module flag
                       cluster.energy,          // total energy
                       cluster.leakage);        // energy from leakage corr

    // do non-linearity energy correction
    if(linear_corr && fabs(alpE) < linear_corr_limit) {
        float corr = 1./(1 + alpE);
        // save the correction factor, not alpha(E)
        hycal_hit.lin_corr = corr;
        hycal_hit.E *= corr;
    }

    // count modules
    hycal_hit.nblocks = cluster.hits.size();

    // fill 3x3 hits around center into temp container for position reconstruction
    BaseHit cl[POS_RECON_HITS];
    int count = fillHits(cl, POS_RECON_HITS, cluster.center, cluster.hits);

    // record how many hits participated in position reconstruction
    hycal_hit.npos = count;

    // reconstruct position
    reconstructPos(cl, count, (BaseHit*)&hycal_hit);
    hycal_hit.z = cluster.center.geo.z;

    // z position will need a depth correction
    hycal_hit.z += GetShowerDepth(cluster.center.geo.type, cluster.energy);

    return hycal_hit;
}

// leakage correction, dead module hits will be provided by hycal detector
void PRadHyCalCluster::LeakCorr(ModuleCluster &cluster, const std::vector<ModuleHit> &dead)
const
{
    if(!leak_corr)
        return;

    if(TEST_BIT(cluster.center.flag, kDeadNeighbor))
        AddVirtHits(cluster, dead);

    if(TEST_BIT(cluster.center.flag, kInnerBound))
        AddVirtHits(cluster, inner_virtual);

    if(TEST_BIT(cluster.center.flag, kOuterBound))
        AddVirtHits(cluster, outer_virtual);
}

// add virtual hits to correct energy leakage
void PRadHyCalCluster::AddVirtHits(ModuleCluster &cluster, const std::vector<ModuleHit> &dead)
const
{
    if(dead.empty())
        return;

    const auto &center = cluster.center;

    // temp container to reconstruct position
    BaseHit cl[POS_RECON_HITS], temp_hit(center.geo.x, center.geo.y, 0., cluster.energy);

    // estimator to check if virtual hits will improve the profile
    float estimator = __hc_prof.EvalEstimator(temp_hit, cluster);

    // this cluster is too bad
    if(estimator > 5.)
        return;

    // temporty container for dead hits energies
    float dead_energy[dead.size()], temp_energy[dead.size()];
    // initialize
    for(unsigned int i = 0; i < dead.size(); ++i)
    {
        dead_energy[i] = 0.;
    }

    // iteration to correct leakage
    for(unsigned int iter = 0; iter < leak_iters; ++iter)
    {
        // check profile to update dead hits' energies
        for(unsigned int i = 0; i < dead.size(); ++i)
        {
            float frac = __hc_prof.GetProfile(temp_hit.x, temp_hit.y, dead.at(i)).frac;
            // full correction would be frac/(1 - frac), but it may result in divergence
            temp_energy[i] = cluster.energy*frac;
        }

        // reconstruct position using cluster hits and dead modules
        // fill existing cluster hits
        int count = fillHits(cl, POS_RECON_HITS, center, cluster.hits);

        // fill virtual hits for dead modules
        temp_hit.E = cluster.energy;
        for(unsigned int i = 0; i < dead.size(); ++i)
        {
            if(temp_energy[i] == 0.)
                continue;

            const auto &hit = dead.at(i);
            if(PRadHyCalDetector::hit_distance(center, hit) < CORNER_ADJACENT) {
                cl[count].x = hit.geo.x;
                cl[count].y = hit.geo.y;
                cl[count].E = temp_energy[i];
                temp_hit.E += temp_energy[i];
                count++;
            }
        }

        // reconstruct position
        reconstructPos(cl, count, &temp_hit);

        // check if the correction helps improve the cluster profile
        float new_est = __hc_prof.EvalEstimator(temp_hit, cluster);
        // not improving, stop!
        if(new_est > estimator)
            break;

        // improved! apply changes
        estimator = new_est;
        for(unsigned int i = 0; i < dead.size(); ++i)
        {
            dead_energy[i] = temp_energy[i];
        }
    }

    // leakage correction to cluster
    for(unsigned int i = 0; i < dead.size(); ++i)
    {
        // leakage is large enough
        if(dead_energy[i] >= least_leak*cluster.energy) {

            // add virtual hit
            ModuleHit vhit(dead.at(i));
            vhit.energy = dead_energy[i];
            cluster.AddHit(vhit);

            // record leakage correction
            cluster.leakage += vhit.energy;
        }
    }
}

// only use the center 3x3 to fill the temp container
inline int PRadHyCalCluster::fillHits(BaseHit *temp,
                                      int max_hits,
                                      const ModuleHit &center,
                                      const std::vector<ModuleHit> &hits)
const
{
    int count = 0;
    for(auto &hit : hits)
    {
        if(count >= max_hits) {
            std::cout << "PRad HyCal Cluster Warning: Exceeds the "
                      << "hits limit (" << max_hits << ") "
                      << "for  position reconstruction."
                      << std::endl;
            break;
        }

        if(PRadHyCalDetector::hit_distance(center, hit) < CORNER_ADJACENT) {
            temp[count].x = hit.geo.x;
            temp[count].y = hit.geo.y;
            temp[count].E = hit.energy;
            count++;
        }
    }
    return count;
}

// reconstruct position from the temp container
void PRadHyCalCluster::reconstructPos(BaseHit *temp, int count, BaseHit *recon)
const
{
    // get total energy
    float energy = 0;
    for(int i= 0; i < count; ++i)
    {
        energy += temp[i].E;
    }

    // reconstruct position
    float wx = 0, wy = 0, wtot = 0;
    for(int i = 0; i < count; ++i)
    {
        float weight = GetWeight(temp[i].E, energy);
        wx += temp[i].x*weight;
        wy += temp[i].y*weight;
        wtot += weight;
    }

    recon->x = wx/wtot;
    recon->y = wy/wtot;
}

// correct virtual hits energy if we know the real positon (from other detector)
void PRadHyCalCluster::CorrectVirtHits(ModuleCluster &cluster, float x, float y)
const
{
    // no need to correct
    if(cluster.leakage == 0.)
        return;

    // re-calculate cluster energy, zero leakage energy
    cluster.energy -= cluster.leakage;
    cluster.leakage = 0.;

    // change virtual hits
    for(auto it = cluster.hits.begin(); it != cluster.hits.end(); ++it)
    {
        // real hit, no need to correct
        if(it->real)
            continue;

        // check profile
        float frac = __hc_prof.GetProfile(x, y, *it).frac;

        // full energy correction because we trust the position
        if(frac > 0. && frac < 1.) {
            it->energy = cluster.energy*frac/(1 - frac);
            cluster.leakage += it->energy;
        // remove virtual hit if its energy should be zero
        } else {
            cluster.hits.erase(it--);
        }
    }

    // update energy
    cluster.energy += cluster.leakage;
}

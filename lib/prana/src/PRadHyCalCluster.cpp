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



PRadHyCalCluster::PRadHyCalCluster()
: detector(nullptr), depth_corr(true), leak_corr(true), linear_corr(true),
  log_weight_thres(3.6), min_cluster_energy(30.), min_center_energy(10.),
  least_leak(0.05), linear_corr_limit(0.6), min_cluster_size(1), leak_iters(3)
{
    // place holder
}

PRadHyCalCluster::~PRadHyCalCluster()
{
    // place holder
}

PRadHyCalCluster* PRadHyCalCluster::Clone() const
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

}

void PRadHyCalCluster::CollectHits(PRadHyCalDetector *det)
{
    // collect hits
    module_hits.clear();

    for(auto &module : det->GetModuleList())
    {
        float energy = module->GetEnergy();
        if(energy > 0) {
            module_hits.emplace_back(module, module->GetID(), energy);
        }
    }
}

void PRadHyCalCluster::Reconstruct(PRadHyCalDetector *det, PRadClusterProfile *p)
{
    detector = det;
    profile = p;

    // form clusters
    FormCluster(module_hits, module_clusters);

    // reconstruct hit from cluster
    det->ClearHits();
    for(auto &cluster : module_clusters)
    {
        // discard cluster that does not satisfy certain conditions
        if(!CheckCluster(cluster))
            continue;

        // leakage correction for dead modules
        LeakCorr(cluster);

        // reconstruct hit the position based on the cluster
        HyCalHit hit = reconstructHit(cluster);

        // final hit reconstructed
        det->AddHit(std::move(hit));
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

// leakage correction, dead module hits will be provided by hycal detector
void PRadHyCalCluster::LeakCorr(ModuleCluster &cluster)
const
{
    if(!leak_corr ||                            // correction disabled
       TEST_BIT(cluster.flag, kLeakCorr) ||     // already corrected
       cluster.hits.size() < 4)                 // insufficient hits to constrain
        return;

    const auto &vnbrs = cluster.center->GetVirtNeighbors();

    // no need to correct
    if(vnbrs.empty())
        return;

    // add virtual hits for each virtual neighbor module
    std::vector<ModuleHit> vhits;
    vhits.reserve(vnbrs.size());
    for(auto &vnbr : vnbrs)
    {
        vhits.emplace_back(vnbr.ptr, vnbr->GetID(), 0., false);
    }

    // reconstruct hit position
    BaseHit pos(0., 0., 0., cluster.energy);
    reconstructPos(cluster, &pos);
    // get estimator for the current cluster
    double est = evalCluster(pos, cluster);

    // iteration to correct virtual hits
    int iter = leak_iters;
    std::vector<double> vhits_e(vhits.size());
    while(iter-- > 0) {
        // save current status of vhits
        for(size_t i = 0; i < vhits.size(); ++i)
        {
            vhits_e[i] = vhits[i].energy;
        }

        // correct virtual hits, reconstruct hit position
        CorrectVirtHits(pos, vhits, cluster);
        // re-evalate the cluster profile
        double new_est = evalCluster(pos, cluster);

        // worse
        if(new_est >= est) {
            // restore energies
            for(size_t i = 0; i < vhits.size(); ++i)
            {
                vhits[i].energy = vhits_e[i];
            }
            // done
            break;
        } else {
            est = new_est;
        }
    }

    // apply the virtual hits
    for(auto &vhit : vhits)
    {
        if(vhit.energy <= 0.) continue;
        cluster.hits.push_back(vhit);
        cluster.energy += vhit.energy;
        cluster.leakage += vhit.energy;
    }

    // including leakage correction
    SET_BIT(cluster.flag, kLeakCorr);
}

// correct virtual hits energy with given hit, and update the hit after correction
void PRadHyCalCluster::CorrectVirtHits(BaseHit &hit, std::vector<ModuleHit> &vhits,
                                       const ModuleCluster &cluster)
const
{
    // update virtual hit energy
    double tote = cluster.energy;
    for(auto &vhit : vhits)
    {
        // check profile
        float frac = getProf(hit.x, hit.y, cluster.energy, vhit).frac;

        double ene;
        if(frac > least_leak && frac < 1.) {
            ene = hit.E*frac;
        } else {
            ene = 0.;
        }
        vhit.energy = ene;
        tote += ene;
    }

    // reconstruct the position
    BaseHit temp[POS_RECON_HITS];
    int count = fillHits(temp, POS_RECON_HITS, cluster.center, cluster.hits);
    count += fillHits(&temp[count], POS_RECON_HITS - count, cluster.center, vhits);
    reconstructPos(cluster.center, temp, count, &hit);
    hit.E = tote;
}

HyCalHit PRadHyCalCluster::ReconstructHit(const ModuleCluster &cluster, PRadHyCalDetector *det)
{
    detector = det;
    return reconstructHit(cluster);
}

double PRadHyCalCluster::EvalCluster(const BaseHit &hit, const ModuleCluster &cluster, PRadClusterProfile *prof)
{
    profile = prof;
    return evalCluster(hit, cluster);
}

// reconstruct hit from cluster
HyCalHit PRadHyCalCluster::reconstructHit(const ModuleCluster &cluster)
const
{

    // initialize the hit
    HyCalHit hycal_hit(cluster.center.id,               // center id
                       cluster.flag,                    // cluster flag
                       cluster.energy,                  // total energy
                       cluster.leakage);                // energy from leakage corr

    // count modules
    hycal_hit.nblocks = cluster.hits.size();

    // reconstruct position
    hycal_hit.npos = reconstructPos(cluster, (BaseHit*)&hycal_hit);

    // get non-linear correction factor
    float alpE = cluster.center->GetCalibConst().NonLinearCorr(cluster.energy);

    // do non-linearity energy correction
    if(linear_corr && fabs(alpE) < linear_corr_limit) {
        float corr = 1./(1. + alpE);
        // save the correction factor, not alpha(E)
        hycal_hit.lin_corr = corr;
        hycal_hit.E *= corr;
    }

    // z position will need a depth correction
    hycal_hit.z += GetShowerDepth(cluster.center->GetType(), cluster.energy);

    return hycal_hit;
}

// only use the center 3x3 to fill the temp container
inline int PRadHyCalCluster::fillHits(BaseHit *temp, int max_hits,
                                      const ModuleHit &center,
                                      const std::vector<ModuleHit> &hits)
const
{
    int count = 0;
    for(auto &hit : hits)
    {
        if(center.id == hit.id) continue;
        if(count >= max_hits) {
            std::cout << "PRad HyCal Cluster Warning: Exceeds the "
                      << "hits limit (" << max_hits << ") "
                      << "for  position reconstruction."
                      << std::endl;
            break;
        }

        double dx, dy;
        detector->QuantizedDist(center.ptr, hit.ptr, dx, dy);
        if(std::abs(dx) < 1.01 && std::abs(dy) < 1.01) {
            temp[count].x = dx;
            temp[count].y = dy;
            temp[count].E = hit.energy;
            count++;
        }
    }
    return count;
}

// reconstruct position from the temp container
int PRadHyCalCluster::reconstructPos(const ModuleHit &center,
                                     BaseHit *temp, int count, BaseHit *hit)
const
{
    // get total energy
    float energy = center.energy;
    for(int i= 0; i < count; ++i)
    {
        energy += temp[i].E;
    }

    // reconstruct position
    float wx = 0, wy = 0, wtot = GetWeight(center.energy, energy);
    // center only contains a little portion of the total energy
    // possibly cosmic
    if(wtot == 0.) {
        hit->x = center->GetX();
        hit->y = center->GetY();
        hit->z = center->GetZ();
        return 1;
    }

    int phits = 0;
    for(int i = 0; i < count; ++i)
    {
        float weight = GetWeight(temp[i].E, energy);
        if(weight > 0.) {
            wx += temp[i].x*weight;
            wy += temp[i].y*weight;
            wtot += weight;
            phits++;
        }
    }

    hit->x = center->GetX() + wx/wtot*center->GetSizeX();
    hit->y = center->GetY() + wy/wtot*center->GetSizeY();
    hit->z = center->GetZ();

    return phits;
}

// reconstruct position from cluster
int PRadHyCalCluster::reconstructPos(const ModuleCluster &cl, BaseHit *hit)
const
{
    BaseHit temp[POS_RECON_HITS];
    int count = fillHits(temp, POS_RECON_HITS, cl.center, cl.hits);
    return reconstructPos(cl.center, temp, count, hit);
}

// get profile values from PRadClusterProfile
typedef PRadClusterProfile::Value ProfVal;

ProfVal PRadHyCalCluster::getProf(const ModuleHit &c, const ModuleHit &hit)
const
{
    // magic number 0.78, the center module contains about 78% of the total energy
    return profile->GetProfile(c->GetType(), hitDistance(c, hit), c.energy/0.78);
}

ProfVal PRadHyCalCluster::getProf(double cx, double cy, double cE, const ModuleHit &hit)
const
{
    int sid = detector->GetSectorID(cx, cy);
    int type = detector->GetSectorInfo().at(sid).mtype;
    double dist = detector->QuantizedDist(cx, cy, sid,
                                          hit->GetX(), hit->GetY(), hit->GetSectorID());
    return profile->GetProfile(type, dist, cE);
}

ProfVal PRadHyCalCluster::getProf(const BaseHit &c, const ModuleHit &hit)
const
{
    int sid = detector->GetSectorID(c.x, c.y);
    int type = detector->GetSectorInfo().at(sid).mtype;
    double dist = detector->QuantizedDist(c.x, c.y, sid,
                                          hit->GetX(), hit->GetY(), hit->GetSectorID());
    return profile->GetProfile(type, dist, c.E);
}

// evaluate how well this cluster can be described by the profile
double PRadHyCalCluster::evalCluster(const BaseHit &c, const ModuleCluster &cl)
const
{
    double est = 0.;

    // determine energy resolution
    double res = 0.026;  // 2.6% for PbWO4
    if(TEST_BIT(cl.flag, kPbGlass))
        res = 0.065;    // 6.5% for PbGlass
    if(TEST_BIT(cl.flag, kTransition))
        res = 0.050;    // 5.0% for transition
    res /= sqrt(c.E/1000.);

    int count = 0;
    for(auto &hit : cl.hits)
    {
        auto prof = getProf(c, hit);
        if(prof.frac < 0.01)
          continue;

        ++count;

        double diff = hit.energy - c.E*prof.frac;
        double sigma2 = 0.816*hit.energy + res*c.E*prof.err;

        // log likelyhood for double exponential distribution
        est += fabs(diff)/sqrt(sigma2);
    }

    return est/count;
}


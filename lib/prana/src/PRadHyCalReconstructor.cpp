//============================================================================//
// Class to manage HyCal clustering method and reconstruct hits for HyCal     //
// Different clustering methods can be switched                               //
//                                                                            //
// Chao Peng                                                                  //
// 02/18/2018                                                                 //
//============================================================================//

#include "PRadHyCalReconstructor.h"
#include "PRadSquareCluster.h"
#include "PRadIslandCluster.h"
#include "PRadHyCalDetector.h"
#include "PRadHyCalSystem.h"
#include "PRadClusterProfile.h"
#include <unordered_map>
#include <iostream>



// constructor
PRadHyCalReconstructor::PRadHyCalReconstructor(const std::string &conf_path)
: type(Undefined), method(nullptr)
{
    // default island method
    SetMethod(Island);
    Configure(conf_path);
}

// copy/move constructor
PRadHyCalReconstructor::PRadHyCalReconstructor(const PRadHyCalReconstructor &that)
: type(that.type), profile(that.profile)
{
    method = that.method->Clone();
}

PRadHyCalReconstructor::PRadHyCalReconstructor(PRadHyCalReconstructor &&that)
: type(that.type), profile(std::move(that.profile))
{
    method = that.method;
    that.method = nullptr;
}

// destructor
PRadHyCalReconstructor::~PRadHyCalReconstructor()
{
    delete method, method = nullptr;
}

// copy/move assignment operator
PRadHyCalReconstructor &PRadHyCalReconstructor::operator =(const PRadHyCalReconstructor &rhs)
{
    type = rhs.type;
    method = rhs.method->Clone();
    profile = rhs.profile;
    return *this;
}

PRadHyCalReconstructor &PRadHyCalReconstructor::operator =(PRadHyCalReconstructor &&rhs)
{
    type = rhs.type;
    method = rhs.method;
    rhs.method = nullptr;
    profile = std::move(rhs.profile);
    return *this;
}


void PRadHyCalReconstructor::Configure(const std::string &path)
{
    bool verbose = false;

    if(!path.empty()) {
        ConfigObject::Configure(path);
        verbose = true;
    }

    depth_corr = getDefConfig<bool>("Shower Depth Correction", true, verbose);
    leak_corr = getDefConfig<bool>("Leakage Correction", true, verbose);
    linear_corr = getDefConfig<bool>("Non Linearity Correction", true, verbose);
    corner_conn = getDefConfig<bool>("Corner Connection", false, verbose);

    log_weight_thres = getDefConfig<float>("Log Weight Threshold", 3.6, verbose);
    min_cluster_energy = getDefConfig<float>("Minimum Cluster Energy", 50., verbose);
    min_center_energy = getDefConfig<float>("Minimum Center Energy", 10., verbose);
    min_cluster_size = getDefConfig<unsigned int>("Minimum Cluster Size", 1, verbose);
    least_leak = getDefConfig<float>("Least Leakage Fraction", 0.05, verbose);
    leak_iters = getDefConfig<unsigned int>("Leakage Iterations", 3, verbose);
    linear_corr_limit = getDefConfig<float>("Non Linearity Limit", 0.6, verbose);
    split_iter = getDefConfig<unsigned int>("Split Iteration", 6, verbose);
    least_share = getDefConfig<float>("Least Split Fraction", 0.01, verbose);

    square_size = getDefConfig<unsigned int>("Square Size", 5, verbose);

    // set the min module energy for all the module type
    float univ_min_energy = getDefConfig<float>("Min Module Energy", 0., false);
    min_module_energy.resize(PRadHyCalModule::Max_Types, univ_min_energy);

    // update the min module energy if some type is specified
    // the key is "Min Module Energy [typename]"
    for(unsigned int i = 0; i < min_module_energy.size(); ++i)
    {
        // determine key name
        std::string type = PRadHyCalModule::Type2str(i);
        std::string key = "Min Module Energy [" + type + "]";
        auto value = GetConfigValue(key);
        if(!value.IsEmpty())
            min_module_energy[i] = value.Float();
    }
}

// reconstruct the event to clusters
void PRadHyCalReconstructor::Reconstruct(PRadHyCalDetector *hycal, const EventData &event)
{
    // cannot reconstruct without necessary objects
    if(!hycal || !method) {
        std::cerr << "PRad HyCal Reconstructor Error: undefined method or null "
                  << "detector pointer. Abort event reconstruction."
                  << std::endl;
        return;
    }

    // no need to reconstruct non-physics event
    if(!event.is_physics_event())
        return;

    // collect hits
    CollectHits(hycal, event);

    ReconstructHits(hycal);

    // add timing information
    AddTiming(hycal, event);
}

void PRadHyCalReconstructor::Reconstruct(PRadHyCalDetector *hycal)
{
    // cannot reconstruct without necessary objects
    if(!hycal || !method) {
        std::cerr << "PRad HyCal Reconstructor Error: undefined method or null detector pointer. "
                  << "Abort event reconstruction."
                  << std::endl;
        return;
    }

    // collect hits
    CollectHits(hycal);

    // reconstruct
    ReconstructHits(hycal);

    // add timing information
    AddTiming(hycal);
}

// steps to do cluster reconstruction
void PRadHyCalReconstructor::ReconstructHits(PRadHyCalDetector *det)
{
    // form clusters
    method->FormCluster(this);

    // hits container from the detector
    auto &container = det->GetHits();
    container.clear();

    for(auto &cluster : module_clusters)
    {
        // discard cluster that does not satisfy certain conditions
        if(!CheckCluster(cluster))
            continue;

        // leakage correction for dead modules
        LeakCorr(cluster);

        // reconstruct hit the position based on the cluster
        container.emplace_back(Cluster2Hit(cluster));
    }
}

// collect hits from detector
void PRadHyCalReconstructor::CollectHits(PRadHyCalDetector *det)
{
    module_hits.clear();

    for(auto &module : det->GetModuleList())
    {
        float energy = module->GetEnergy();
        if(energy > min_module_energy[module->GetType()]) {
            module_hits.emplace_back(module, module->GetID(), energy);
        }
    }
}

// collect hits from event
void PRadHyCalReconstructor::CollectHits(PRadHyCalDetector *det, const EventData &event)
{
    module_hits.clear();

    auto sys = det->GetSystem();

    if(!sys) {
        std::cerr << "PRad HyCal Reconstructor Error: cannot reconstruct a given "
                  << "event without a DAQ system connected to the detector. "
                  << "No module hits collected."
                  << std::endl;
        return;
    }

    for(auto adc : event.get_adc_data())
    {
        auto channel = sys->GetADCChannel(adc.channel_id);
        if(!channel) continue;

        auto module = channel->GetModule();
        if(!module) continue;

        double val = (double)adc.value - channel->GetPedestal().mean;
        double energy = module->GetEnergy(val);

        if(energy > min_module_energy[module->GetType()]) {
            module_hits.emplace_back(module, module->GetID(), energy);
        }
    }
}

// add timing information from detector
void PRadHyCalReconstructor::AddTiming(PRadHyCalDetector *det)
{
    // add timing information
    for(auto &hit : det->GetHits())
    {
        auto center = det->GetModule(hit.cid);
        if(!center) continue;

        auto tdc = center->GetTDC();
        if(tdc) hit.set_time(tdc->GetTimeMeasure());
    }
}

// add timing information from event
void PRadHyCalReconstructor::AddTiming(PRadHyCalDetector *det, const EventData &event)
{
    // build map for tdc information
    std::unordered_map<uint16_t, std::vector<uint16_t>> tdc_info;
    for(auto &tdc : event.tdc_data)
    {
        auto it = tdc_info.find(tdc.channel_id);
        if(it == tdc_info.end()) {
            std::vector<uint16_t> info;
            info.push_back(tdc.value);
            tdc_info[tdc.channel_id] = info;
        } else {
            it->second.push_back(tdc.value);
        }
    }

    // add timing information
    for(auto &hit : det->GetHits())
    {
        auto center = det->GetModule(hit.cid);
        if(!center) continue;

        auto tdc = center->GetTDC();
        if(!tdc) continue;

        auto it = tdc_info.find(tdc->GetID());

        if(it != tdc_info.end()) {
            hit.set_time(it->second);
        }
    }
}

// set method
bool PRadHyCalReconstructor::SetMethod(MethodEnum newtype)
{
    if(newtype == type && method) {
        return true;
    }

    // create corresponding method
    PRadHyCalCluster *newone = nullptr;

    switch(newtype)
    {
    case Island: newone = new PRadIslandCluster(); break;
    case Square: newone = new PRadSquareCluster(); break;
    default: break;
    }

    // warn the failure of set method
    if(newone == nullptr) {
        auto method_names = GetMethodNames();
        std::cout << "PRad HyCal System Warning: Failed to set clustering method. \n"
                  << "Available methods are: \n";

        for(auto &n : method_names)
        {
            std::cout << "\t" << n << "\n";
        }
        std::cout << std::endl;

        return false;
    }

    // free previous method
    delete method;
    method = newone;
    type = newtype;
    return true;
}

// set method
bool PRadHyCalReconstructor::SetMethod(const std::string &name)
{
    return SetMethod(str2MethodEnum(name.c_str()));
}

// show available methods
std::vector<std::string> PRadHyCalReconstructor::GetMethodNames()
const
{
    std::vector<std::string> res;
    for(int i = 0; i < static_cast<int>(Max_Methods); ++i)
    {
        res.emplace_back(MethodEnum2str(i));
    }

    return res;
}

// check if the cluster is good enough
bool PRadHyCalReconstructor::CheckCluster(const ModuleCluster &cluster)
const
{
    if((cluster.energy < min_cluster_energy) ||
       (cluster.hits.size() < min_cluster_size))
            return false;

    return true;
}

// reconstruct hit from cluster
HyCalHit PRadHyCalReconstructor::Cluster2Hit(const ModuleCluster &cluster)
const
{

    // initialize the hit
    HyCalHit hycal_hit(cluster.center.id,               // center id
                       cluster.flag,                    // cluster flag
                       cluster.energy,                  // total energy
                       cluster.leakage);                // energy from leakage corr

    // count modules
    hycal_hit.nblocks = cluster.hits.size();

    // position from module hits, record the number of participants
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
    hycal_hit.z += getShowerDepth(cluster.center->GetType(), cluster.energy);

    return hycal_hit;
}

// leakage correction, dead module hits will be provided by hycal detector
void PRadHyCalReconstructor::LeakCorr(ModuleCluster &cluster)
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
    double est = EvalCluster(pos, cluster);

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
        double new_est = EvalCluster(pos, cluster);

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
void PRadHyCalReconstructor::CorrectVirtHits(BaseHit &hit, std::vector<ModuleHit> &vhits,
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

// evaluate how well this cluster can be described by the profile
double PRadHyCalReconstructor::EvalCluster(const BaseHit &c, const ModuleCluster &cl)
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

// get position weight
float PRadHyCalReconstructor::getWeight(const float &E, const float &E0)
const
{
    float w = log_weight_thres + log(E/E0);
    if(w < 0.)
        return 0.;
    return w;
}

// get shower depth, unit is in MeV
float PRadHyCalReconstructor::getShowerDepth(int module_type, const float &E)
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

// reconstruct position from the temp container
int PRadHyCalReconstructor::reconstructPos(const ModuleHit &center,
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
    float wx = 0, wy = 0, wtot = getWeight(center.energy, energy);
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
        float weight = getWeight(temp[i].E, energy);
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
int PRadHyCalReconstructor::reconstructPos(const ModuleCluster &cl, BaseHit *hit)
const
{
    BaseHit temp[POS_RECON_HITS];
    int count = fillHits(temp, POS_RECON_HITS, cl.center, cl.hits);
    return reconstructPos(cl.center, temp, count, hit);
}

// only use the center 3x3 to fill the temp container
int PRadHyCalReconstructor::fillHits(BaseHit *temp, int max_hits,
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
        center->GetDetector()->QuantizedDist(center.ptr, hit.ptr, dx, dy);
        if(std::abs(dx) < 1.01 && std::abs(dy) < 1.01) {
            temp[count].x = dx;
            temp[count].y = dy;
            temp[count].E = hit.energy;
            count++;
        }
    }
    return count;
}


// get profiles
typedef PRadClusterProfile::Value ProfVal;

ProfVal PRadHyCalReconstructor::getProf(const ModuleHit &c, const ModuleHit &hit)
const
{
    double dist = c->GetDetector()->QuantizedDist(c.ptr, hit.ptr);
    // magic number 0.78, the center module contains about 78% of the total energy
    return profile.Get(c->GetType(), dist, c.energy/0.78);
}

ProfVal PRadHyCalReconstructor::getProf(double cx, double cy, double cE, const ModuleHit &hit)
const
{
    auto detector = hit->GetDetector();
    int sid = detector->GetSectorID(cx, cy);
    int type = detector->GetSectorInfo().at(sid).mtype;
    double dist = detector->QuantizedDist(cx, cy, sid,
                                          hit->GetX(), hit->GetY(), hit->GetSectorID());
    return profile.Get(type, dist, cE);
}

ProfVal PRadHyCalReconstructor::getProf(const BaseHit &c, const ModuleHit &hit)
const
{
    auto detector = hit->GetDetector();
    int sid = detector->GetSectorID(c.x, c.y);
    int type = detector->GetSectorInfo().at(sid).mtype;
    double dist = detector->QuantizedDist(c.x, c.y, sid,
                                          hit->GetX(), hit->GetY(), hit->GetSectorID());
    return profile.Get(type, dist, c.E);
}


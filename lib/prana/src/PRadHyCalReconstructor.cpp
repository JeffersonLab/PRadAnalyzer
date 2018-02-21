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



//============================================================================//
// Constructor, Destructor, Assignment Operators                              //
//============================================================================//

// constructor
PRadHyCalReconstructor::PRadHyCalReconstructor(const std::string &conf_path)
: cltype(Undefined_ClMethod), cluster(nullptr), postype(Logarithmic)
{
    // default island method
    SetClusterMethod(Island);
    Configure(conf_path);
}

// copy/move constructor
PRadHyCalReconstructor::PRadHyCalReconstructor(const PRadHyCalReconstructor &that)
: ConfigObject(that), profile(that.profile),
  cltype(that.cltype), postype(that.postype), config(that.config)
{
    cluster = that.cluster->Clone();
}

PRadHyCalReconstructor::PRadHyCalReconstructor(PRadHyCalReconstructor &&that)
: ConfigObject(that), profile(std::move(that.profile)),
  cltype(that.cltype), postype(that.postype), config(std::move(that.config))
{
    cluster = that.cluster;
    that.cluster = nullptr;
}

// destructor
PRadHyCalReconstructor::~PRadHyCalReconstructor()
{
    delete cluster, cluster = nullptr;
}

// copy/move assignment operator
PRadHyCalReconstructor &PRadHyCalReconstructor::operator =(const PRadHyCalReconstructor &rhs)
{
    if(this == &rhs)
        return *this;

    PRadHyCalReconstructor that(rhs); // copy constructor
    *this = std::move(that); // move assignment operator
    return *this;
}

PRadHyCalReconstructor &PRadHyCalReconstructor::operator =(PRadHyCalReconstructor &&rhs)
{
    if(this == &rhs)
        return *this;

    ConfigObject::operator =(rhs);
    profile = std::move(rhs.profile);
    cltype = rhs.cltype;
    cluster = rhs.cluster;
    rhs.cluster = nullptr;
    postype = rhs.postype;
    config = std::move(rhs.config);
    return *this;
}


//============================================================================//
// Public Member Functions                                                    //
//============================================================================//

// configuration
void PRadHyCalReconstructor::Configure(const std::string &path)
{
    bool verbose = false;

    if(!path.empty()) {
        ConfigObject::Configure(path);
        verbose = true;
    }

    // general
    CONF_CONN(config.depth_corr, "Shower Depth Correction", true, verbose);
    CONF_CONN(config.leak_corr, "Leakage Correction", true, verbose);
    CONF_CONN(config.linear_corr, "Non Linearity Correction", true, verbose);
    CONF_CONN(config.pos_s_corr, "S Shape Correction", false, verbose);

    CONF_CONN(config.log_weight_thres, "Log Weight Threshold", 3.6, verbose);
    CONF_CONN(config.min_cluster_energy, "Minimum Cluster Energy", 50., verbose);
    CONF_CONN(config.min_center_energy, "Minimum Center Energy", 10., verbose);
    CONF_CONN(config.min_cluster_size, "Minimum Cluster Size", 1, verbose);
    CONF_CONN(config.least_leak, "Least Leakage Fraction", 0.05, verbose);
    CONF_CONN(config.leak_iters, "Leakage Iterations", 3, verbose);
    CONF_CONN(config.linear_corr_limit, "Non Linearity Limit", 0.6, verbose);

    // square
    CONF_CONN(config.square_size, "Square Size", 5, verbose);

    // island
    CONF_CONN(config.corner_conn, "Corner Connection", false, verbose);
    CONF_CONN(config.split_iter, "Split Iteration", 6, verbose);
    CONF_CONN(config.least_split, "Least Split Fraction", 0.01, verbose);

    // default min module energy
    config.min_module_energy.resize(static_cast<int>(PRadHyCalModule::Max_Types), 0.);
    // update the min module energy if some type is specified
    // the key is "Min Module Energy [typename]"
    for(int i = 0; i < (int)config.min_module_energy.size(); ++i)
    {
        // determine key name
        std::string type = PRadHyCalModule::Type2str(i);
        std::string key = "Min Module Energy [" + type + "]";
        auto value = GetConfigValue(key);
        if(!value.IsEmpty())
            config.min_module_energy[i] = value.Float();
    }

    // s shape correction parameters
    config.pos_s_pars.resize(static_cast<int>(Max_PosMethods));
    int max_pars = 4;
    for(int i = 0; i < (int)config.pos_s_pars.size(); ++i)
    {
        std::string type = PosMethod2str(i);
        std::string key = "S Shape Parameters [" + type + "]";
        auto valstr = GetConfigValue(key);
        auto vals = ConfigParser::split(valstr, ",");  // split input string

        auto &par_pack = config.pos_s_pars.at(i);
        par_pack.resize(max_pars, 0.);
        for(int j = 0; j < max_pars && !vals.empty(); ++j)
        {
            // trim off white spaces
    	    ConfigValue val = ConfigParser::trim(vals.front(), " \t");
    	    vals.pop_front();
            par_pack[j] = val.Float();
        }
    }
}

// reconstruct the event to clusters
void PRadHyCalReconstructor::Reconstruct(PRadHyCalDetector *hycal, const EventData &event)
{
    // cannot reconstruct without necessary objects
    if(!hycal || !cluster) {
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
    if(!hycal || !cluster) {
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
    cluster->FormCluster(module_hits, module_clusters);

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
        if(energy > config.min_module_energy[module->GetType()]) {
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

        if(energy > config.min_module_energy[module->GetType()]) {
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

// set cluster method by enum
bool PRadHyCalReconstructor::SetClusterMethod(ClMethod newtype)
{
    if(newtype == cltype && cluster) {
        return true;
    }

    // create corresponding method
    PRadHyCalCluster *newone = nullptr;

    switch(newtype)
    {
    case Island: newone = new PRadIslandCluster(this); break;
    case Square: newone = new PRadSquareCluster(this); break;
    default: break;
    }

    // warn the failure of set method
    if(newone == nullptr) {
        std::cout << "PRad HyCal Reconstructor Warning: Failed to set clustering "
                  << "method. \nAvailable methods are: \n";

        auto method_names = GetClusterMethodNames();
        for(auto &n : method_names)
        {
            std::cout << "\t" << n << "\n";
        }
        std::cout << std::endl;

        return false;
    }

    // free previous method
    delete cluster;
    cluster = newone;
    cltype = newtype;
    return true;
}

// set cluster method by name
bool PRadHyCalReconstructor::SetClusterMethod(const std::string &name)
{
    return SetClusterMethod(str2ClMethod(name.c_str()));
}

// show available cluster methods
std::vector<std::string> PRadHyCalReconstructor::GetClusterMethodNames()
const
{
    std::vector<std::string> res;
    for(int i = 0; i < static_cast<int>(Max_ClMethods); ++i)
    {
        res.emplace_back(ClMethod2str(i));
    }

    return res;
}

// set position reconstruction method by enum
bool PRadHyCalReconstructor::SetPositionMethod(PosMethod newtype)
{
    int id = static_cast<int>(newtype);
    if(id < 0 || id >= static_cast<int>(Max_PosMethods)) {
        std::cout << "PRad HyCal Reconstructor Warning: Failed to set clustering "
                  << "method. \nAvailable methods are: \n";

        auto method_names = GetPositionMethodNames();
        for(auto &n : method_names)
        {
            std::cout << "\t" << n << "\n";
        }
        std::cout << std::endl;

        return false;
    }

    postype = newtype;
    return true;
}

// set position reconstruction method by name
bool PRadHyCalReconstructor::SetPositionMethod(const std::string &name)
{
    return SetPositionMethod(str2PosMethod(name.c_str()));
}

// show available position methods
std::vector<std::string> PRadHyCalReconstructor::GetPositionMethodNames()
const
{
    std::vector<std::string> res;
    for(int i = 0; i < static_cast<int>(Max_PosMethods); ++i)
    {
        res.emplace_back(PosMethod2str(i));
    }

    return res;
}

// check if the cluster is good enough
bool PRadHyCalReconstructor::CheckCluster(const ModuleCluster &cluster)
const
{
    if((cluster.energy < config.min_cluster_energy) ||
       (cluster.hits.size() < config.min_cluster_size))
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
    if(config.linear_corr && fabs(alpE) < config.linear_corr_limit) {
        float corr = 1./(1. + alpE);
        // save the correction factor, not alpha(E)
        hycal_hit.lin_corr = corr;
        hycal_hit.E *= corr;
    }

    // z position will need a depth correction
    hycal_hit.z += getShowerDepth(cluster.center->GetType(), cluster.energy);

    // add resolution information
    auto detector = cluster.center->GetDetector();
    hycal_hit.sig_ene = detector->GetEneRes(cluster.center.ptr, hycal_hit.E);
    hycal_hit.sig_pos = detector->GetPosRes(cluster.center.ptr, hycal_hit.E);

    return hycal_hit;
}

// leakage correction, dead module hits will be provided by hycal detector
void PRadHyCalReconstructor::LeakCorr(ModuleCluster &cluster)
const
{
    if(!config.leak_corr ||                     // correction disabled
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
    int iter = config.leak_iters;
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
        if(frac > config.least_leak && frac < 1.) {
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
        double dE = res*c.E;
        double sigma2 = c.E*c.E*prof.err*prof.err + dE*dE*prof.frac*prof.frac;

        // log likelyhood for double exponential distribution
        est += fabs(diff)/sqrt(sigma2);
    }

    return est/count;
}

//============================================================================//
// Protected/Private Functions                                                //
//============================================================================//

// get position weight
float PRadHyCalReconstructor::getWeight(const float &E, const float &E0)
const
{
    float w;
    switch(postype)
    {
    default:
    case Logarithmic:
        w = config.log_weight_thres + std::log(E/E0);
        break;
    case Linear:
        w = E/E0;
        break;
    }

    if(w < 0.)
        return 0.;
    return w;
}

// get shower depth, unit is in MeV
float PRadHyCalReconstructor::getShowerDepth(int module_type, const float &E)
const
{
    if(config.depth_corr && E > 0.) {
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

// correct S shape bias in position reconstruction
float PRadHyCalReconstructor::getPosBias(const float &dx)
const
{
    // formula is
    // dx*(dx^2 - 0.25)*(c0 + c1*dx^2 + c2*dx^4 + c3*dx^6)
    float dx2 = dx*dx;
    float res = 0.;
    float fac = 1.;
    for(auto &par : config.pos_s_pars[static_cast<size_t>(postype)])
    {
        res += fac*par;
        fac *= dx2;
    }

    return dx*res*(dx2 - 0.25);
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

    float dx = wx/wtot;
    float dy = wy/wtot;

    hit->x = center->GetX() + dx*center->GetSizeX();
    hit->y = center->GetY() + dy*center->GetSizeY();
    hit->z = center->GetZ();

    // currently supports crystal modules
    if(config.pos_s_corr && center->GetType() == PRadHyCalModule::PbWO4) {
        hit->x += getPosBias(dx);
        hit->y += getPosBias(dy);
    }

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


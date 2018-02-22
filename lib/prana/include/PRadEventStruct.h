#ifndef PRAD_EVENT_STRUCT_H
#define PRAD_EVENT_STRUCT_H

// class, structure and enum related the reconstructed event

#include <string>
#include <vector>
#include <deque>
#include <utility>
#include <cstdint>
#include "generalstruct.h"
#include "datastruct.h"

// some discriminator related settings
#define REF_CHANNEL 7
#define REF_PULSER_FREQ 500000

#define FCUP_CHANNEL 6
#define FCUP_OFFSET 100.0
#define FCUP_SLOPE 906.2

//============================================================================//
// *BEGIN* RUN INFORMATION STRUCTURE                                          //
//============================================================================//
struct RunInfo
{
    int32_t run_number;
    double beam_charge;
    double dead_count;
    double ungated_count;

    RunInfo()
    : run_number(0), beam_charge(0.), dead_count(0.), ungated_count(0.)
    {}

    RunInfo(const int &run, const double &c, const double &d, const double &ug)
    : run_number(run), beam_charge(c), dead_count(d), ungated_count(ug)
    {}

    void reset()
    {
        run_number = 0;
        beam_charge = 0.;
        dead_count = 0.;
        ungated_count = 0.;
    }
};
//============================================================================//
// *END* RUN INFORMATION STRUCTURE                                            //
//============================================================================//



//============================================================================//
// *BEGIN* ONLINE INFORMATION STRUCTURE                                       //
//============================================================================//
struct TriggerChannel
{
    std::string name;
    uint32_t id;
    double freq;

    TriggerChannel()
    : name("undefined"), id(0), freq(0.)
    {}
    TriggerChannel(const std::string &n, const uint32_t &i)
    : name(n), id(i), freq(0.)
    {}
};

// online information
struct OnlineInfo
{
    double live_time;
    double beam_current;
    std::vector<TriggerChannel> trigger_info;

    OnlineInfo()
    : live_time(0.), beam_current(0.)
    {}

    void add_trigger(const std::string &n, const uint32_t &i)
    {
        trigger_info.emplace_back(n, i);
    }

    void reset()
    {
        live_time = 0.;
        beam_current = 0.;
    }

    void clear()
    {
        reset();
        trigger_info.clear();
    }
};
//============================================================================//
// *END* ONLINE INFORMATION STRUCTURE                                         //
//============================================================================//



//============================================================================//
// *BEGIN* RAW EPICS DATA STRUCTURE                                           //
//============================================================================//
struct EpicsData
{
    int32_t event_number;
    std::vector<float> values;

    EpicsData()
    {}
    EpicsData(const int &ev, const std::vector<float> &val)
    : event_number(ev), values(val)
    {}

    void clear()
    {
        event_number = 0;
        values.clear();
    }

    bool operator <(const int &evt) const {return event_number < evt;}
    bool operator >(const int &evt) const {return event_number > evt;}
    bool operator <=(const int &evt) const {return event_number <= evt;}
    bool operator >=(const int &evt) const {return event_number >= evt;}
    bool operator ==(const int &evt) const {return event_number == evt;}
    bool operator !=(const int &evt) const {return event_number != evt;}
};
//============================================================================//
// *END* RAW EPICS DATA STRUCTURE                                             //
//============================================================================//



//============================================================================//
// *BEGIN* RAW EVENT DATA COMPONENTS                                          //
//============================================================================//
typedef struct ChannelData
{
    uint16_t channel_id;
    uint16_t value;

    ChannelData()
    : channel_id(0), value(0)
    {}
    ChannelData(const uint16_t &i, const uint16_t &v)
    : channel_id(i), value(v)
    {}

} TDC_Data, ADC_Data;

struct DSC_Data
{
    uint32_t gated_count;
    uint32_t ungated_count;

    DSC_Data()
    : gated_count(0), ungated_count(0)
    {}
    DSC_Data(const uint32_t &g1, const uint32_t &g2)
    : gated_count(g1), ungated_count(g2)
    {}
};

struct GEMChannelAddress
{
    uint8_t fec;
    uint8_t adc;
    uint8_t strip;

    GEMChannelAddress() {}
    GEMChannelAddress(const uint8_t &f,
                      const uint8_t &a,
                      const uint8_t &s)
    : fec(f), adc(a), strip(s)
    {}
};

struct GEM_Data
{
    GEMChannelAddress addr;
    std::vector<float> values;

    GEM_Data() {}
    GEM_Data(const uint8_t &f,
             const uint8_t &a,
             const uint8_t &s)
    : addr(f, a, s)
    {}

    void set_address (const uint8_t &f,
                      const uint8_t &a,
                      const uint8_t &s)
    {
        addr.fec = f;
        addr.adc = a;
        addr.strip = s;
    }

    void add_value(const float &v)
    {
        values.push_back(v);
    }
};

//============================================================================//
// *END* RAW EVENT DATA COMPONENTS                                            //
//============================================================================//



//============================================================================//
// *BEGIN* RAW EVENT DATA STRUCTURE                                           //
//============================================================================//
struct EventData
{
    // event info
    int32_t event_number;
    uint8_t type;
    uint8_t trigger;
    uint64_t timestamp;

    // data banks
    std::vector< ADC_Data > adc_data;
    std::vector< TDC_Data > tdc_data;
    std::vector< GEM_Data > gem_data;
    std::vector< DSC_Data > dsc_data;

    // constructors
    EventData()
    : event_number(0), type(0), trigger(0), timestamp(0)
    {}
    EventData(const uint8_t &t)
    : event_number(0), type(t), trigger(0), timestamp(0)
    {}
    EventData(const uint8_t &t,
              const PRadTriggerType &trg,
              std::vector<ADC_Data> &adc,
              std::vector<TDC_Data> &tdc,
              std::vector<GEM_Data> &gem,
              std::vector<DSC_Data> &dsc)
    : event_number(0), type(t), trigger((uint8_t)trg), timestamp(0),
      adc_data(adc), tdc_data(tdc), gem_data(gem), dsc_data(dsc)
    {}

    void clear()
    {
        event_number = 0;
        type = 0;
        trigger = 0;
        timestamp = 0;
        adc_data.clear();
        tdc_data.clear();
        gem_data.clear();
        dsc_data.clear();
    }

    void update_type(const uint8_t &t) {type = t;}
    void update_trigger(const uint8_t &t) {trigger = t;}
    void update_time(const uint64_t &t) {timestamp = t;}

    uint32_t get_type() const {return type;}
    uint32_t get_trigger() const {return trigger;}
    uint64_t get_time() const {return timestamp;}

    void add_adc(const ADC_Data &a) {adc_data.emplace_back(a);}
    void add_tdc(const TDC_Data &t) {tdc_data.emplace_back(t);}
    void add_gemhit(const GEM_Data &g) {gem_data.emplace_back(g);}
    void add_dsc(const DSC_Data &d) {dsc_data.emplace_back(d);}

    void add_adc(ADC_Data &&a) {adc_data.emplace_back(a);}
    void add_tdc(TDC_Data &&t) {tdc_data.emplace_back(t);}
    void add_gemhit(GEM_Data &&g) {gem_data.emplace_back(g);}
    void add_dsc(DSC_Data &&d) {dsc_data.emplace_back(d);}

    std::vector<ADC_Data> &get_adc_data() {return adc_data;}
    std::vector<TDC_Data> &get_tdc_data() {return tdc_data;}
    std::vector<GEM_Data> &get_gem_data() {return gem_data;}
    std::vector<DSC_Data> &get_dsc_data() {return dsc_data;}

    const std::vector<ADC_Data> &get_adc_data() const {return adc_data;}
    const std::vector<TDC_Data> &get_tdc_data() const {return tdc_data;}
    const std::vector<GEM_Data> &get_gem_data() const {return gem_data;}
    const std::vector<DSC_Data> &get_dsc_data() const {return dsc_data;}

    bool is_physics_event()
    const
    {
        return ( (trigger == PHYS_LeadGlassSum) ||
                 (trigger == PHYS_TotalSum)     ||
                 (trigger == PHYS_TaggerE)      ||
                 (trigger == PHYS_Scintillator) );
    }

    bool is_monitor_event()
    const
    {
        return ( (trigger == LMS_Led) ||
                 (trigger == LMS_Alpha) );
    }

    bool is_sync_event()
    const
    {
        return type == CODA_Sync;
    }

    double get_beam_time()
    const
    {
        double elapsed_time = 0.;
        if(dsc_data.size() > REF_CHANNEL)
        {
            elapsed_time = (double)dsc_data.at(REF_CHANNEL).ungated_count
                          /(double)REF_PULSER_FREQ;
        }

        return elapsed_time;
    }

    double get_live_time()
    const
    {
        double live_time = 1.;
        if(dsc_data.size() > REF_CHANNEL)
        {
            live_time -= (double)dsc_data.at(REF_CHANNEL).gated_count
                        /(double)dsc_data.at(REF_CHANNEL).ungated_count;
        }

        return live_time;
    }

    double get_beam_charge()
    const
    {
        double beam_charge = 0.;
        if(dsc_data.size() > FCUP_CHANNEL)
        {
            beam_charge = ((double)dsc_data.at(FCUP_CHANNEL).ungated_count - FCUP_OFFSET)/FCUP_SLOPE;
        }

        return beam_charge;
    }

    double get_beam_current()
    const
    {
        double beam_time = get_beam_time();
        if(beam_time > 0.)
            return get_beam_charge()/beam_time;
        else
            return 0.;
    }

    DSC_Data get_dsc_channel(const uint32_t &idx)
    const
    {
        if(dsc_data.size() <= idx)
            return DSC_Data();
        else
            return dsc_data.at(idx);
    }

    DSC_Data get_ref_channel()
    const
    {
        return get_dsc_channel(REF_CHANNEL);
    }

    DSC_Data get_dsc_scaled_by_ref(const uint32_t &idx)
    const
    {
        if(idx >= dsc_data.size())
            return DSC_Data();

        uint64_t ref_pulser = get_ref_channel().ungated_count;
        uint64_t ungated_scaled = (dsc_data.at(idx).ungated_count*ref_pulser)/REF_PULSER_FREQ;
        uint64_t gated_scaled = (dsc_data.at(idx).gated_count*ref_pulser)/REF_PULSER_FREQ;

        return DSC_Data((uint32_t)gated_scaled, (uint32_t)ungated_scaled);
    }

    bool operator == (const int &ev) const {return ev == event_number;}
    bool operator != (const int &ev) const {return ev != event_number;}
    bool operator > (const int &ev) const {return ev > event_number;}
    bool operator < (const int &ev) const {return ev < event_number;}
    bool operator > (const EventData &other) const {return other.event_number > event_number;}
    bool operator < (const EventData &other) const {return other.event_number < event_number;}
};
//============================================================================//
// *END* RAW EVENT DATA STRUCTURE                                             //
//============================================================================//



//============================================================================//
// *BEGIN* HYCAL MODULE HIT AND CLUSTER STRUCTURE                             //
//============================================================================//

struct ModuleHit
{
    class PRadHyCalModule *ptr;     // pointer to the hycal module
    int32_t id;                     // module id
    float energy;                   // participated energy, may be splitted
    bool real;                      // false for virtual hit to correct leakage

    ModuleHit(bool r = true)
    : ptr(nullptr), id(0), energy(0), real(r)
    {}

    ModuleHit(PRadHyCalModule *p, int i, float e, bool r = true)
    : ptr(p), id(i), energy(e), real(r)
    {}

    bool operator ==(const ModuleHit &rhs) const {return id == rhs.id;}
    PRadHyCalModule *operator->() {return ptr;}
    const PRadHyCalModule *operator->() const {return ptr;}
};

struct ModuleCluster
{
    ModuleHit center;               // center hit
    std::vector<ModuleHit> hits;    // hits group
    float energy;                   // cluster energy
    float leakage;                  // energy leakage
    uint32_t flag;                  // cluster flag

    ModuleCluster()
    : energy(0), leakage(0), flag(0)
    {
        hits.reserve(100);
    }

    ModuleCluster(const ModuleHit &hit, uint32_t f)
    : center(hit), energy(0), leakage(0), flag(f)
    {
        hits.reserve(100);
    }

    void AddHit(const ModuleHit &hit)
    {
        hits.emplace_back(hit);
        energy += hit.energy;
    }

    void Merge(const ModuleCluster &that)
    {
        hits.reserve(hits.size() + that.hits.size());
        hits.insert(hits.end(), that.hits.begin(), that.hits.end());

        energy += that.energy;
        leakage += that.leakage;

        if(center.energy < that.center.energy) {
            center = that.center;
            flag = that.flag;
        }
    }

    void FindCenter()
    {
        if(hits.empty())
            return;

        float max_e = center.energy;
        ModuleHit *cptr = nullptr;
        for(auto &hit : hits)
        {
            if(hit.energy > max_e) {
                max_e = hit.energy;
                cptr = &hit;
            }
        }

        if(cptr)
            center = *cptr;
    }
};

//============================================================================//
// *END* HYCAL MODULE HIT AND CLUSTER STRUCTURE                               //
//============================================================================//



//============================================================================//
// *BEGIN* GEM STRIP HIT AND CLUSTER STRUCTURE                                //
//============================================================================//
struct StripHit
{
    int32_t strip;
    float charge;
    float position;
    bool cross_talk;
    APVAddress apv_addr;

    StripHit()
    : strip(0), charge(0.), position(0.), cross_talk(false), apv_addr(-1, -1)
    {}
    StripHit(int s, float c, float p, bool f = false, int fec = -1, int adc = -1)
    : strip(s), charge(c), position(p), cross_talk(f), apv_addr(fec, adc)
    {}
};

struct StripCluster
{
    float position;
    float peak_charge;
    float total_charge;
    std::vector<StripHit> hits;

    StripCluster()
    : position(0.), peak_charge(0.), total_charge(0.)
    {}

    StripCluster(const std::vector<StripHit> &p)
    : position(0.), peak_charge(0.), total_charge(0.), hits(p)
    {}

    StripCluster(std::vector<StripHit> &&p)
    : position(0.), peak_charge(0.), total_charge(0.), hits(std::move(p))
    {}

    bool IsCrossTalk()
    const
    {
        for(auto &hit : hits)
        {
            if(!hit.cross_talk)
                return false;
        }
        return true;
    }
};

//============================================================================//
// *END* GEM STRIP HIT AND CLUSTER STRUCTURE                                  //
//============================================================================//



//============================================================================//
// *BEGIN* DETECTOR HIT STRUCTURE                                             //
//============================================================================//

// status enums require bitwise manipulation
// the following defintions are copies from Rtypes.h in root (cern)
#define SET_BIT(n,i)  ( (n) |= (1ULL << i) )
#define CLEAR_BIT(n,i)  ( (n) &= ~(1ULL << i) )
#define TEST_BIT(n,i)  ( (bool)( n & (1ULL << i) ) )

// common part between gem and hycal hits
class BaseHit
{
public:
    float x;            // Cluster's x-position (mm)
    float y;            // Cluster's y-position (mm)
    float z;            // Cluster's z-position (mm)
    float E;            // Cluster's energy (MeV)

    BaseHit()
    : x(0.), y(0.), z(0.), E(0.)
    {}

    BaseHit(float xi, float yi, float zi, float Ei)
    : x(xi), y(yi), z(zi), E(Ei)
    {}
};

// hycal hit status
enum HyCalHitStatus
{
    kPbGlass = 0,       // cluster center at lead glass region
    kPbWO4,             // cluster center at lead tungstate region
    kTransition,        // cluster center at transition region
    kSplit,             // cluster after splitting
    kDeadModule,        // cluster center is a dead module (only possible from leakage correction)
    kDeadNeighbor,      // cluster center is near a dead module
    kInnerBound,        // cluster near the inner hole of HyCal
    kOuterBound,        // cluster near the outer boundary of HyCal
    kLeakCorr,          // cluster with leakage correction
    kDenCorr,           // cluster with density correction (for position bias)
    kSEneCorr,          // cluster with S shape energy correction (for energy bias)
};

// hycal reconstructed hit
class HyCalHit : public BaseHit
{
public:
#define TIME_MEASURE_SIZE 3
    uint32_t flag;          // overall status of the cluster
    int16_t type;           // Cluster types: 0,1,2,3,4;-1
    int16_t status;         // Spliting status
    int16_t nblocks;        // Number of blocks in a cluster
    int16_t npos;           // Number of blocks participated in position reconstruction
    int16_t cid;            // Cluster's central cell ID
    float E_leak;           // Leakage correction on energy (MeV)
    float lin_corr;         // Non Linearity factor for energy correction E_f = E_i*lin_corr
    float E_Scorr;          // S-shape correction on energy (MeV)
    float sig_ene;          // energy resolution
    float sig_pos;          // position resolution
    uint16_t time[TIME_MEASURE_SIZE];      // time information from central TDC group

    HyCalHit()
    : flag(0), type(0), status(0), nblocks(0), npos(0), cid(0), E_leak(0.),
      lin_corr(1.), E_Scorr(0.), sig_ene(0.), sig_pos(0.)
    {
        clear_time();
    }

    HyCalHit(int16_t id, uint32_t f, float ene, float leak)
    : BaseHit(0., 0., 0., ene), flag(f), type(0), status(0), nblocks(0), npos(0),
      cid(id), E_leak(leak), lin_corr(1.), E_Scorr(0.), sig_ene(0.), sig_pos(0.)
    {
        clear_time();
    }

    void clear_time()
    {
        for(int i = 0; i < TIME_MEASURE_SIZE; ++i)
            time[i] = 0;
    }

    void set_time(const std::vector<uint16_t> &t)
    {
        for(int i = 0; i < TIME_MEASURE_SIZE; ++i)
        {
            if(i < (int)t.size())
                time[i] = t[i];
            else
                time[i] = 0;
        }
    }
};

// gem reconstructed hit
class GEMHit : public BaseHit
{
public:
    int32_t det_id;         // which GEM detector it belongs to
    float x_charge;         // x charge
    float y_charge;         // y charge
    float x_peak;           // x peak charge
    float y_peak;           // y peak charge
    int32_t x_size;         // x hits size
    int32_t y_size;         // y hits size
    float sig_pos;          // position resolution

    GEMHit()
    : det_id(-1), x_charge(0.), y_charge(0.), x_peak(0.), y_peak(0.),
      x_size(0), y_size(0), sig_pos(0.)
    {}

    GEMHit(float xx, float yy, float zz, int d, float xc, float yc,
           float xp, float yp, int xs, int ys, float sig)
    : BaseHit(xx, yy, zz, 0.), det_id(d), x_charge(xc), y_charge(yc),
      x_peak(xp), y_peak(yp), x_size(xs), y_size(ys), sig_pos(sig)
    {}

};


// hit matching status
enum MatchHitStatus
{
    // since this status flag is originally from HyCalHitStatus
    // start with its original value to avoid mismatch with old code

    kGEM1Match = 8,     // found a matched GEM hit on PRadGEM1
    kGEM2Match,         // found a matched GEM hit on PRadGEM2
};

class MatchHit : public BaseHit
{
public:
    HyCalHit hycal;
    GEMHit gem;
    std::vector<GEMHit> gem1;
    std::vector<GEMHit> gem2;
    uint32_t mflag;
    // this index is kept because of example/physCalib is using it
    // TODO revamp physCalib and remove this member
    uint32_t hycal_idx;

    MatchHit(const HyCalHit &hit)
    : BaseHit(hit.x, hit.y, hit.z, hit.E), hycal(hit), mflag(0)
    {}

    MatchHit(const HyCalHit &hit, std::vector<GEMHit> &&v1, std::vector<GEMHit> &&v2)
    : BaseHit(hit.x, hit.y, hit.z, hit.E), hycal(hit), gem1(v1), gem2(v2), mflag(0)
    {}

    MatchHit(const HyCalHit &hit, const std::vector<GEMHit> &v1, const std::vector<GEMHit> &v2)
    : BaseHit(hit.x, hit.y, hit.z, hit.E), hycal(hit), gem1(v1), gem2(v2), mflag(0)
    {}

    void SubstituteCoord(const BaseHit &h)
    {
        x = h.x;
        y = h.y;
        z = h.z;
    }
};

// a function to get the hit with the greatest energy
template<typename Iter>
inline Iter MostEnergeticHit(Iter beg, Iter end)
{
    if(beg == end) return end;

    Iter res = beg;
    for(++beg; beg != end; ++beg)
    {
        if(beg->E > res->E)
            res = beg;
    }

    return res;
}

//============================================================================//
// *END* DETECTOR HIT STRUCTURE                                               //
//============================================================================//

#endif

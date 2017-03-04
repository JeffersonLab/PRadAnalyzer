#ifndef PRAD_EVENT_STRUCT_H
#define PRAD_EVENT_STRUCT_H

#include <string>
#include <vector>
#include <deque>
#include <utility>
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
    int run_number;
    double beam_charge;
    double dead_count;
    double ungated_count;

    RunInfo()
    : run_number(0), beam_charge(0.), dead_count(0.), ungated_count(0.)
    {};

    RunInfo(const int &run, const double &c, const double &d, const double &ug)
    : run_number(run), beam_charge(c), dead_count(d), ungated_count(ug)
    {};

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
    {};
    TriggerChannel(const std::string &n, const uint32_t &i)
    : name(n), id(i), freq(0.)
    {};
};

// online information
struct OnlineInfo
{
    double live_time;
    double beam_current;
    std::vector<TriggerChannel> trigger_info;

    OnlineInfo()
    : live_time(0.), beam_current(0.)
    {};

    void add_trigger(const std::string &n, const uint32_t &i)
    {
        trigger_info.emplace_back(n, i);
    };

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
    int event_number;
    std::vector<float> values;

    EpicsData()
    {};
    EpicsData(const int &ev, const std::vector<float> &val)
    : event_number(ev), values(val)
    {};

    void clear()
    {
        event_number = 0;
        values.clear();
    };

    bool operator <(const int &evt) const {return event_number < evt;};
    bool operator >(const int &evt) const {return event_number > evt;};
    bool operator <=(const int &evt) const {return event_number <= evt;};
    bool operator >=(const int &evt) const {return event_number >= evt;};
    bool operator ==(const int &evt) const {return event_number == evt;};
    bool operator !=(const int &evt) const {return event_number != evt;};
};
//============================================================================//
// *END* RAW EPICS DATA STRUCTURE                                             //
//============================================================================//

//============================================================================//
// *BEGIN* RAW EVENT DATA COMPONENTS                                          //
//============================================================================//
typedef struct ChannelData
{
    unsigned short channel_id;
    unsigned short value;

    ChannelData()
    : channel_id(0), value(0)
    {};
    ChannelData(const unsigned short &i, const unsigned short &v)
    : channel_id(i), value(v)
    {};

} TDC_Data, ADC_Data;

struct DSC_Data
{
    unsigned int gated_count;
    unsigned int ungated_count;

    DSC_Data()
    : gated_count(0), ungated_count(0)
    {};
    DSC_Data(const unsigned int &g1, const unsigned int &g2)
    : gated_count(g1), ungated_count(g2)
    {};
};

struct GEMChannelAddress
{
    unsigned char fec;
    unsigned char adc;
    unsigned char strip;

    GEMChannelAddress() {};
    GEMChannelAddress(const unsigned char &f,
                      const unsigned char &a,
                      const unsigned char &s)
    : fec(f), adc(a), strip(s)
    {};
};

struct GEM_Data
{
    GEMChannelAddress addr;
    std::vector<float> values;

    GEM_Data() {};
    GEM_Data(const unsigned char &f,
             const unsigned char &a,
             const unsigned char &s)
    : addr(f, a, s)
    {};

    void set_address (const unsigned char &f,
                      const unsigned char &a,
                      const unsigned char &s)
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
    int event_number;
    unsigned char type;
    unsigned char trigger;
    uint64_t timestamp;

    // data banks
    std::vector< ADC_Data > adc_data;
    std::vector< TDC_Data > tdc_data;
    std::vector< GEM_Data > gem_data;
    std::vector< DSC_Data > dsc_data;

    // constructors
    EventData()
    : event_number(0), type(0), trigger(0), timestamp(0)
    {};
    EventData(const unsigned char &t)
    : event_number(0), type(t), trigger(0), timestamp(0)
    {};
    EventData(const unsigned char &t,
              const PRadTriggerType &trg,
              std::vector<ADC_Data> &adc,
              std::vector<TDC_Data> &tdc,
              std::vector<GEM_Data> &gem,
              std::vector<DSC_Data> &dsc)
    : event_number(0), type(t), trigger((unsigned char)trg), timestamp(0),
      adc_data(adc), tdc_data(tdc), gem_data(gem), dsc_data(dsc)
    {};

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
    };

    void update_type(const unsigned char &t) {type = t;};
    void update_trigger(const unsigned char &t) {trigger = t;};
    void update_time(const uint64_t &t) {timestamp = t;};

    unsigned int get_type() const {return type;};
    unsigned int get_trigger() const {return trigger;};
    uint64_t get_time() const {return timestamp;};

    void add_adc(const ADC_Data &a) {adc_data.emplace_back(a);};
    void add_tdc(const TDC_Data &t) {tdc_data.emplace_back(t);};
    void add_gemhit(const GEM_Data &g) {gem_data.emplace_back(g);};
    void add_dsc(const DSC_Data &d) {dsc_data.emplace_back(d);};

    void add_adc(ADC_Data &&a) {adc_data.emplace_back(a);};
    void add_tdc(TDC_Data &&t) {tdc_data.emplace_back(t);};
    void add_gemhit(GEM_Data &&g) {gem_data.emplace_back(g);};
    void add_dsc(DSC_Data &&d) {dsc_data.emplace_back(d);};

    std::vector<ADC_Data> &get_adc_data() {return adc_data;};
    std::vector<TDC_Data> &get_tdc_data() {return tdc_data;};
    std::vector<GEM_Data> &get_gem_data() {return gem_data;};
    std::vector<DSC_Data> &get_dsc_data() {return dsc_data;};

    const std::vector<ADC_Data> &get_adc_data() const {return adc_data;};
    const std::vector<TDC_Data> &get_tdc_data() const {return tdc_data;};
    const std::vector<GEM_Data> &get_gem_data() const {return gem_data;};
    const std::vector<DSC_Data> &get_dsc_data() const {return dsc_data;};

    bool is_physics_event()
    const
    {
        return ( (trigger == PHYS_LeadGlassSum) ||
                 (trigger == PHYS_TotalSum)     ||
                 (trigger == PHYS_TaggerE)      ||
                 (trigger == PHYS_Scintillator) );
    };

    bool is_monitor_event()
    const
    {
        return ( (trigger == LMS_Led) ||
                 (trigger == LMS_Alpha) );
    };

    bool is_sync_event()
    const
    {
        return type == CODA_Sync;
    };

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
    };

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
    };

    double get_beam_charge()
    const
    {
        double beam_charge = 0.;
        if(dsc_data.size() > FCUP_CHANNEL)
        {
            beam_charge = ((double)dsc_data.at(FCUP_CHANNEL).ungated_count - FCUP_OFFSET)/FCUP_SLOPE;
        }

        return beam_charge;
    };

    double get_beam_current()
    const
    {
        double beam_time = get_beam_time();
        if(beam_time > 0.)
            return get_beam_charge()/beam_time;
        else
            return 0.;
    };

    DSC_Data get_dsc_channel(const uint32_t &idx)
    const
    {
        if(dsc_data.size() <= idx)
            return DSC_Data();
        else
            return dsc_data.at(idx);
    };

    DSC_Data get_ref_channel()
    const
    {
        return get_dsc_channel(REF_CHANNEL);
    };

    DSC_Data get_dsc_scaled_by_ref(const uint32_t &idx)
    const
    {
        if(idx >= dsc_data.size())
            return DSC_Data();

        uint64_t ref_pulser = get_ref_channel().ungated_count;
        uint64_t ungated_scaled = (dsc_data.at(idx).ungated_count*ref_pulser)/REF_PULSER_FREQ;
        uint64_t gated_scaled = (dsc_data.at(idx).gated_count*ref_pulser)/REF_PULSER_FREQ;

        return DSC_Data((unsigned int)gated_scaled, (unsigned int)ungated_scaled);
    };

    bool operator == (const int &ev) const {return ev == event_number;};
    bool operator != (const int &ev) const {return ev != event_number;};
    bool operator > (const int &ev) const {return ev > event_number;};
    bool operator < (const int &ev) const {return ev < event_number;};
    bool operator > (const EventData &other) const {return other.event_number > event_number;};
    bool operator < (const EventData &other) const {return other.event_number < event_number;};
};
//============================================================================//
// *END* RAW EVENT DATA STRUCTURE                                             //
//============================================================================//

//============================================================================//
// *BEGIN* CLUSTER STRUCTURE                                                  //
//============================================================================//

// common part between gem and hycal cluster
class BaseHit
{
public:
    float x;            // Cluster's x-position (mm)
    float y;            // Cluster's y-position (mm)
    float z;            // Cluster's z-position (mm)
    float E;            // Cluster's energy (MeV)

    BaseHit()
    : x(0.), y(0.), z(0.), E(0.)
    {};

    BaseHit(float xi, float yi, float zi, float Ei)
    : x(xi), y(yi), z(zi), E(Ei)
    {};
};

// hycal cluster status
enum HyCalHitStatus
{
// this enum is for bitwise manipulation
// the following defintions are from Rtypes.h in root (cern)
// copied here for a clearer reference
#define SET_BIT(n,i)  ( (n) |= (1ULL << i) )
#define CLEAR_BIT(n,i)  ( (n) &= ~(1ULL << i) )
#define TEST_BIT(n,i)  ( (bool)( n & (1ULL << i) ) )

    kPbGlass = 0,       //cluster center at lead glass region
    kPbWO4,             //cluster center at lead tungstate region
    kTransition,        //cluster center at transition region
    kSplit,             //cluster after splitting
    kDeadModule,        //cluster center is a dead module (only possible from leakage correction)
    kDeadNeighbor,      //cluster center is near a dead module
    kInnerBound,        //cluster near the inner hole of HyCal
    kOuterBound,        //cluster near the outer boundary of HyCal
    kGEM1Match,         //cluster found a match GEM hit on GEM 1
    kGEM2Match,         //cluster found a match GEM hit on GEM 2
};

// hycal cluster
class HyCalHit : public BaseHit
{
public:
#define TIME_MEASURE_SIZE 3
    unsigned int flag;  // overall status of the cluster
    short type;         // Cluster types: 0,1,2,3,4;-1
    short status;       // Spliting status
    short nblocks;      // Number of blocks in a cluster
    short npos;         // Number of blocks participated in position reconstruction
    short cid;          // Cluster's central cell ID
    float E_leak;       // Leakage correction on energy (MeV)
    float lin_corr;     // Non Linearity factor for energy correction E_f = E_i*lin_corr
    unsigned short time[TIME_MEASURE_SIZE];      // time information from central TDC group

    HyCalHit()
    : flag(0), type(0), status(0), nblocks(0), npos(0), cid(0), E_leak(0.), lin_corr(1.)
    {
        clear_time();
    }

    HyCalHit(short id, unsigned int f, float ene, float leak)
    : BaseHit(0., 0., 0., ene), flag(f), type(0), status(0), nblocks(0), npos(0),
      cid(id), E_leak(leak), lin_corr(1.)
    {
        clear_time();
    }

    void clear_time()
    {
        for(int i = 0; i < TIME_MEASURE_SIZE; ++i)
            time[i] = 0;
    }

    void set_time(const std::vector<unsigned short> &t)
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

// gem cluster
class GEMHit : public BaseHit
{
public:
    int det_id;       // which GEM detector it belongs to
    float x_charge;   // x charge
    float y_charge;   // y charge
    float x_peak;     // x peak charge
    float y_peak;     // y peak charge
    int x_size;       // x hits size
    int y_size;       // y hits size

    GEMHit()
    : det_id(-1), x_charge(0.), y_charge(0.), x_peak(0.), y_peak(0.),
      x_size(0), y_size(0)
    {};

    GEMHit(float x, float y, float z,
           int d, float xc, float yc, float xp, float yp, int xs, int ys)
    : BaseHit(x, y, z, 0.), det_id(d), x_charge(xc), y_charge(yc),
      x_peak(xp), y_peak(yp), x_size(xs), y_size(ys)
    {};

};

class MatchHit : public BaseHit
{
public:
    HyCalHit hycal;
    std::vector<GEMHit> gem1;
    std::vector<GEMHit> gem2;
    // this index is kept because of decoder/physCalib is using it
    // TODO revamp physCalib and remove this member
    uint32_t hycal_idx;

    MatchHit(const HyCalHit &hit)
    : BaseHit(hit.x, hit.y, hit.z, hit.E), hycal(hit)
    {};

    MatchHit(const HyCalHit &hit, std::vector<GEMHit> &&v1, std::vector<GEMHit> &&v2)
    : BaseHit(hit.x, hit.y, hit.z, hit.E), hycal(hit), gem1(v1), gem2(v2)
    {};

    MatchHit(const HyCalHit &hit, const std::vector<GEMHit> &v1, const std::vector<GEMHit> &v2)
    : BaseHit(hit.x, hit.y, hit.z, hit.E), hycal(hit), gem1(v1), gem2(v2)
    {};

    void SubstituteCoord(const BaseHit &h)
    {
        x = h.x;
        y = h.y;
        z = h.z;
    }
};

//============================================================================//
// *END* CLUSTER STRUCTURE                                                    //
//============================================================================//

#endif

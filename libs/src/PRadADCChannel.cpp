//============================================================================//
// Basic DAQ channel unit                                                     //
//                                                                            //
// Chao Peng                                                                  //
// 02/17/2016                                                                 //
//============================================================================//

#include "PRadADCChannel.h"
#include "PRadTDCChannel.h"
#include "PRadHyCalModule.h"
#include <iomanip>
#include <utility>



//============================================================================//
// Constructors, Destructors, Assignment Operators                            //
//============================================================================//

// constructor
PRadADCChannel::PRadADCChannel(const std::string &name, const ChannelAddress &daqAddr)
: PRadDAQChannel(name, daqAddr),
  module(nullptr), tdc_group(nullptr), occupancy(0), sparsify(0), adc_value(0)
{
    // initialize histograms
    trg_hist.resize(MAX_Trigger, nullptr);
    std::vector<std::string> hn= {"Physics", "Pedestal", "LMS"};

    AddHist(hn[0], new TH1I((ch_name+"_"+hn[0]).c_str(), hn[0].c_str(), 2048, 0, 8191));
    AddHist(hn[1], new TH1I((ch_name+"_"+hn[1]).c_str(), hn[1].c_str(), 1024, 0, 1023));
    AddHist(hn[2], new TH1I((ch_name+"_"+hn[2]).c_str(), hn[2].c_str(), 2048, 0, 8191));

    // default hist-trigger mapping
    //    MapHist("PED", TI_Error);
    if(ch_name.find("LMS") != std::string::npos) {
        MapHist("Physics", LMS_Alpha);
        MapHist("Pedestal", PHYS_LeadGlassSum);
        MapHist("Pedestal", PHYS_TotalSum);
        MapHist("Pedestal", PHYS_TaggerE);
        MapHist("Pedestal", PHYS_Scintillator);
    } else {
        MapHist("Pedestal", LMS_Alpha);
        MapHist("Physics", PHYS_LeadGlassSum);
        MapHist("Physics", PHYS_TotalSum);
        MapHist("Physics", PHYS_TaggerE);
        MapHist("Physics", PHYS_Scintillator);
    }
    MapHist("LMS", LMS_Led);
}

// copy constructor
PRadADCChannel::PRadADCChannel(const PRadADCChannel &that)
: PRadDAQChannel(that),
  module(nullptr), tdc_group(nullptr), pedestal(that.pedestal),
  occupancy(that.occupancy), sparsify(that.sparsify), adc_value(that.adc_value)
{
    for(auto &it : that.hist_map)
    {
        hist_map[it.first] = (TH1*) it.second->Clone();
    }

    trg_hist.resize(that.trg_hist.size(), nullptr);
    for(size_t i = 0; i < that.trg_hist.size(); ++i)
    {
        if(that.trg_hist[i] != nullptr)
            MapHist(that.trg_hist[i]->GetTitle(), i);
    }
}

// move constructor
PRadADCChannel::PRadADCChannel(PRadADCChannel &&that)
: PRadDAQChannel(that),
  module(nullptr), tdc_group(nullptr), pedestal(that.pedestal),
  occupancy(that.occupancy), sparsify(that.sparsify), adc_value(that.adc_value),
  trg_hist(std::move(that.trg_hist)), hist_map(std::move(that.hist_map))
{
    // place holder
}

// destructor
PRadADCChannel::~PRadADCChannel()
{
    for(auto &ele : hist_map)
    {
        delete ele.second, ele.second = nullptr;
    }

    UnsetTDC();
    UnsetModule();
}

// copy assignment operator
PRadADCChannel &PRadADCChannel::operator =(const PRadADCChannel &rhs)
{
    if(this == &rhs)
        return *this;

    PRadADCChannel that(rhs);
    *this = std::move(that);
    return *this;
}

// move assignment operator
PRadADCChannel &PRadADCChannel::operator =(PRadADCChannel &&rhs)
{
    if(this == &rhs)
        return *this;

    // release memories
    for(auto &it : hist_map)
        delete it.second;
    hist_map.clear();
    trg_hist.clear();

    PRadDAQChannel::operator =(rhs);
    pedestal = rhs.pedestal;
    occupancy = rhs.occupancy;
    sparsify = rhs.sparsify;
    adc_value = rhs.adc_value;
    hist_map = std::move(hist_map);
    trg_hist = std::move(trg_hist);

    return *this;
}



//============================================================================//
// Public Member Functions                                                    //
//============================================================================//

// set a tdc channel
void PRadADCChannel::SetTDC(PRadTDCChannel *t, bool force_set)
{
    if(t == tdc_group)
        return;

    if(!force_set)
        UnsetTDC();

    tdc_group = t;
}

// unset current tdc channel
void PRadADCChannel::UnsetTDC(bool force_unset)
{
    if(!tdc_group)
        return;

    if(!force_unset)
        tdc_group->DisconnectChannel(ch_id, true);

    tdc_group = nullptr;
}

// set a module to this channel
void PRadADCChannel::SetModule(PRadHyCalModule *m, bool force_set)
{
    if(m == module)
        return;

    if(!force_set)
        UnsetModule();

    module = m;
}

// unset current module
void PRadADCChannel::UnsetModule(bool force_unset)
{
    if(!module)
        return;

    if(!force_unset)
        module->UnsetChannel(true);

    module = nullptr;
}

// add a histogram
bool PRadADCChannel::AddHist(const std::string &n, TH1 *hist)
{
    if(hist == nullptr)
        return false;

    std::string key = ConfigParser::str_lower(n);
    auto it = hist_map.find(key);
    if(it != hist_map.end()) {
        std::cerr << "PRad DAQ Channel Error: Failed to add hist " << n
                  << ", a hist with the same name exists."
                  << std::endl;
        return false;
    }
    hist_map[key] = hist;
    return true;

}

// map a histogram to trigger type
bool PRadADCChannel::MapHist(const std::string &name, int trg)
{
    if((size_t)trg >= trg_hist.size())
        return false;

    TH1 *hist = GetHist(name);
    if(hist == nullptr)
        return false;

    size_t index = (size_t) trg;
    trg_hist[index] = hist;

    return true;
}

// get histogram
TH1 *PRadADCChannel::GetHist(const std::string &n)
const
{
    auto it = hist_map.find(ConfigParser::str_lower(n));
    if(it == hist_map.end()) {
        return nullptr;
    }
    return it->second;
}

// get the histogram list
std::vector<TH1*> PRadADCChannel::GetHistList()
const
{
    std::vector<TH1*> hlist;

    for(auto &it : hist_map)
    {
        hlist.push_back(it.second);
    }

    return hlist;
}

// set pedestal
void PRadADCChannel::SetPedestal(const Pedestal &p)
{
    pedestal = p;

    sparsify = (unsigned short)(pedestal.mean + 5.*pedestal.sigma + 0.5); // round
}

// set pedestal
void PRadADCChannel::SetPedestal(const double &m, const double &s)
{
    SetPedestal(Pedestal(m, s));
}

// reset current data
void PRadADCChannel::Reset()
{
    occupancy = 0;
    adc_value = 0;
    ResetHists();
}

// reset histograms
void PRadADCChannel::ResetHists()
{
    for(auto &it : hist_map)
    {
        if(it.second)
            it.second->Reset();
    }
}

// erase histograms
void PRadADCChannel::ClearHists()
{
    for(auto &it : hist_map)
        delete it.second;

    hist_map.clear();

    for(auto &hist : trg_hist)
        hist = nullptr;
}

// zero suppression, triggered when adc value is statistically
// above pedestal (5 sigma)
bool PRadADCChannel::Sparsify(const unsigned short &adcVal)
{
    if(adcVal < sparsify)
        return false;

    ++occupancy;
    return true;
}

double PRadADCChannel::GetEnergy()
const
{
    if(!module)
        return 0;

    return module->GetEnergy((double)adc_value - pedestal.mean);
}

double PRadADCChannel::GetEnergy(const unsigned short &adcVal)
const
{
    if(!module)
        return 0;

    return module->GetEnergy((double)adcVal - pedestal.mean);
}

std::ostream &operator <<(std::ostream &os, const PRadADCChannel::Pedestal &ped)
{
    return os << std::setw(12) << ped.mean
              << std::setw(12) << ped.sigma;
}

std::ostream &operator <<(std::ostream &os, const PRadADCChannel &ch)
{
    return os << (PRadDAQChannel)ch << ch.GetPedestal();
}

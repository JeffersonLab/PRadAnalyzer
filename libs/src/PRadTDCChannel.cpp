//============================================================================//
// a base TDC Channel class, it can connect to multiple ADC Channels          //
//                                                                            //
// Chao Peng                                                                  //
// 12/11/2016                                                                 //
//============================================================================//

#include "PRadTDCChannel.h"
#include "PRadADCChannel.h"
#include "TH1I.h"



//============================================================================//
// Constructors, Destructor, Assignment Operators                             //
//============================================================================//

// constructor
PRadTDCChannel::PRadTDCChannel(const std::string &name, const ChannelAddress &addr)
: PRadDAQChannel(name, addr)
{
    std::string tdc_name = "TDC_" + name;
    tdc_hist = new TH1I(tdc_name.c_str(), "Time Measure", 20000, 0, 19999);
}

// copy/move constructors
// connections to ADC channels won't be copied
PRadTDCChannel::PRadTDCChannel(const PRadTDCChannel &that)
: PRadDAQChannel(that), time_measure(that.time_measure)
{
    if(that.tdc_hist)
        tdc_hist = (TH1*) that.tdc_hist->Clone();
    else
        tdc_hist = nullptr;
}

PRadTDCChannel::PRadTDCChannel(PRadTDCChannel &&that)
: PRadDAQChannel(that), time_measure(std::move(that.time_measure))
{
    tdc_hist = that.tdc_hist;
    that.tdc_hist = nullptr;
}

// destructor
PRadTDCChannel::~PRadTDCChannel()
{
    DisconnectChannels();
    delete tdc_hist;
}

// copy/move assignment operators
PRadTDCChannel &PRadTDCChannel::operator =(const PRadTDCChannel &rhs)
{
    if(this == &rhs)
        return *this;

    delete tdc_hist, tdc_hist = nullptr;

    PRadDAQChannel::operator =(rhs);
    if(rhs.tdc_hist)
        tdc_hist = (TH1*) rhs.tdc_hist->Clone();
    time_measure = rhs.time_measure;

    return *this;
}

PRadTDCChannel &PRadTDCChannel::operator =(PRadTDCChannel &&rhs)
{
    if(this == &rhs)
        return *this;

    delete tdc_hist;

    PRadDAQChannel::operator =(rhs);
    tdc_hist = rhs.tdc_hist;
    rhs.tdc_hist = nullptr;
    time_measure = std::move(rhs.time_measure);

    return *this;
}



//============================================================================//
// Public Member Functions                                                    //
//============================================================================//

// connect adc channel
void PRadTDCChannel::ConnectChannel(PRadADCChannel *ch)
{
    if(!ch)
        return;

    auto it = group_map.find(ch->GetID());
    if(it != group_map.end()) {
        std::cerr << "PRad TDC Channel Error: Failed to connect ADC channel "
                  << ch->GetName() << " (" << ch->GetID() << ") , a channel "
                  << "with same id exists."
                  << std::endl;
        return;
    }

    group_map[ch->GetID()] = ch;
    ch->SetTDC(this);
}

// disconnect adc channel
void PRadTDCChannel::DisconnectChannel(int id, bool force_disconn)
{
    auto it = group_map.find(id);

    if(it == group_map.end())
        return;

    if(!force_disconn)
        it->second->UnsetTDC(true);

    group_map.erase(it);
}

// disconnect all channels
void PRadTDCChannel::DisconnectChannels()
{
    for(auto &it : group_map)
        it.second->UnsetTDC(true);

    group_map.clear();
}

// add a time measure value
void PRadTDCChannel::AddTimeMeasure(const unsigned short &count)
{
    time_measure.push_back(count);
}

// add multiple time measure values
void PRadTDCChannel::AddTimeMeasure(const std::vector<unsigned short> &counts)
{
    time_measure.insert(time_measure.end(), counts.begin(), counts.end());
}

// set time measure values
void PRadTDCChannel::SetTimeMeasure(const std::vector<unsigned short> &counts)
{
    time_measure = counts;
}

// fill histogram
void PRadTDCChannel::FillHist(const unsigned short &time)
{
    tdc_hist->Fill(time);
}

// reset all data
void PRadTDCChannel::Reset()
{
    ResetHists();
    ClearTimeMeasure();
}

// reset histogram
void PRadTDCChannel::ResetHists()
{
    tdc_hist->Reset();
}

// clear time measure
void PRadTDCChannel::ClearTimeMeasure()
{
    time_measure.clear();
}

// get adc channel
PRadADCChannel* PRadTDCChannel::GetADCChannel(int id)
const
{
    auto it = group_map.find(id);
    if(it == group_map.end())
        return nullptr;
    return it->second;
}

std::vector<PRadADCChannel*> PRadTDCChannel::GetChannelList()
const
{
    std::vector<PRadADCChannel*> res;

    for(auto it : group_map)
    {
        res.push_back(it.second);
    }

    return res;
}

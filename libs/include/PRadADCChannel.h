#ifndef PRAD_ADC_CHANNEL_H
#define PRAD_ADC_CHANNEL_H

#include <string>
#include <vector>
#include <iostream>
#include <unordered_map>
#include "PRadDAQChannel.h"
#include "TH1.h"


class PRadHyCalModule;
class PRadTDCChannel;

class PRadADCChannel : public PRadDAQChannel
{
public:
    struct Pedestal
    {
        double mean;
        double sigma;

        Pedestal()
        : mean(0), sigma(0)
        {};
        Pedestal(const double &m, const double &s)
        : mean(m), sigma(s)
        {};
    };

public:
    // constructor
    PRadADCChannel(const std::string &name, const ChannelAddress &daqAddr);

    // copy/move constructor
    PRadADCChannel(const PRadADCChannel &that);
    PRadADCChannel(PRadADCChannel &&that);

    // destructor
    virtual ~PRadADCChannel();

    // copy/move assignment operator
    PRadADCChannel &operator =(const PRadADCChannel &rhs);
    PRadADCChannel &operator =(PRadADCChannel &&rhs);

    // public member functions
    // set members
    void SetTDC(PRadTDCChannel *t, bool force_set = false);
    void UnsetTDC(bool force_unset = false);
    void SetModule(PRadHyCalModule *m, bool force_set = false);
    void UnsetModule(bool force_unset = false);
    void SetPedestal(const Pedestal &ped);
    void SetPedestal(const double &m, const double &s);
    void SetValue(const unsigned short &adc) {adc_value = adc;};
    // reset data
    void Reset();

    // check if adc passed threshold
    void Sparsify() {occupancy++;};
    bool Sparsify (const unsigned short &adcVal);
    int GetOccupancy() const {return occupancy;};
    unsigned short GetValue() const {return adc_value;};
    double GetReducedValue() const {return (double)adc_value - pedestal.mean;};
    double GetEnergy() const;
    double GetEnergy(const unsigned short &adcVal) const;
    const Pedestal &GetPedestal() const {return pedestal;};
    PRadHyCalModule *GetModule() const {return module;};
    PRadTDCChannel *GetTDC() const {return tdc_group;};

    // histograms manipulations
    void ResetHists();
    void ClearHists();
    bool AddHist(const std::string &name, TH1 *hist);
    bool MapHist(const std::string &name, int trg);
    void RemoveHist(const std::string &name);
    template<typename T>
    void FillHist(const T& t, int trg)
    {
        if(trg_hist[trg]) {
            trg_hist[trg]->Fill(t);
        }
    };
    TH1 *GetHist(const std::string &name = "Physics") const;
    TH1 *GetHist(PRadTriggerType type) const {return trg_hist[(int)type];};
    std::vector<TH1*> GetHistList() const;

protected:
    PRadHyCalModule *module;
    PRadTDCChannel *tdc_group;
    Pedestal pedestal;
    int occupancy;
    unsigned short sparsify;
    unsigned short adc_value;
    // histograms
    std::vector<TH1*> trg_hist;
    std::unordered_map<std::string, TH1*> hist_map;
};

std::ostream &operator <<(std::ostream &os, const PRadADCChannel::Pedestal &ped);
std::ostream &operator <<(std::ostream &os, const PRadADCChannel &ch);
#endif

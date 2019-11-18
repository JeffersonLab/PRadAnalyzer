#ifndef PRAD_HYCAL_SYSTEM_H
#define PRAD_HYCAL_SYSTEM_H

#include <string>
#include <ostream>
#include <unordered_map>
#include "PRadEventStruct.h"
#include "PRadHyCalDetector.h"
#include "PRadHyCalReconstructor.h"
#include "PRadTDCChannel.h"
#include "PRadADCChannel.h"
#include "PRadCalibConst.h"
#include "ConfigObject.h"


// adc searching speed is important, thus reserve buckets to have unordered_map
// better formed
#define ADC_BUCKETS 2000

class TH1D;

class PRadHyCalSystem : public ConfigObject
{
public:
    // data structure for finding the calibration period
    struct CalPeriod
    {
        int begin, end;
        int main, sub;

        CalPeriod() : begin(0), end(0), main(0), sub(0) {}
        CalPeriod(int b, int e, int m, int s) : begin(b), end(e), main(m), sub(s) {}

        bool operator ==(int run) const {return (run >= begin) && (run <= end);}
        bool operator !=(int run) const {return (run < begin) || (run > end);}
        bool operator <(int run) const {return end < run;}
        bool operator >(int run) const {return begin > run;}
    };

public:
    // constructor
    PRadHyCalSystem(const std::string &path = "");

    // copy/move constructors
    PRadHyCalSystem(const PRadHyCalSystem &that);
    PRadHyCalSystem(PRadHyCalSystem &&that);

    // destructor
    virtual ~PRadHyCalSystem();

    // copy/move assignment operators
    PRadHyCalSystem &operator =(const PRadHyCalSystem &rhs);
    PRadHyCalSystem &operator =(PRadHyCalSystem &&rhs);

    // configuration
    void Configure(const std::string &path);
    bool ReadChannelList(const std::string &path);
    bool ReadRunInfoFile(const std::string &path);
    bool ReadTriggerEffFile(const std::string &path);
    bool ReadCalPeriodFile(const std::string &path);
    void ChooseRun(const std::string &path, bool verbose = true);
    void ChooseRun(int run, bool verbose = true);
    void UpdateRunFiles(bool verbose = true);

    // connections
    void BuildConnections();

    // events related
    void ChooseEvent(const EventData &data);
    void Reset();
    inline void Reconstruct() {return recon.Reconstruct(hycal);}
    inline void Reconstruct(const EventData &data) {return recon.Reconstruct(hycal, data);}
    PRadHyCalReconstructor *GetReconstructor() {return &recon;}

    // detector related
    void SetDetector(PRadHyCalDetector *h);
    void RemoveDetector();
    void DisconnectDetector(bool force_disconn = false);
    double GetEnergy(const EventData &event) const;
    PRadHyCalModule *GetModule(const int &id) const;
    PRadHyCalModule *GetModule(const std::string &name) const;
    std::vector<PRadHyCalModule*> GetModuleList() const;
    PRadHyCalDetector *GetDetector() const {return hycal;}

    // daq related
    bool AddADCChannel(PRadADCChannel *adc);
    bool AddTDCChannel(PRadTDCChannel *tdc);
    void ClearADCChannel();
    void ClearTDCChannel();
    PRadADCChannel *GetADCChannel(const int &id) const;
    PRadADCChannel *GetADCChannel(const std::string &name) const;
    PRadADCChannel *GetADCChannel(const ChannelAddress &addr) const;
    PRadTDCChannel *GetTDCChannel(const int &id) const;
    PRadTDCChannel *GetTDCChannel(const std::string &name) const;
    PRadTDCChannel *GetTDCChannel(const ChannelAddress &addr) const;
    const std::vector<PRadADCChannel*> &GetADCList() const {return adc_list;}
    const std::vector<PRadTDCChannel*> &GetTDCList() const {return tdc_list;}
    void Sparsify(const EventData &event);

    // histogram related
    void FillHists(const EventData &event);
    void FillEnergyHist();
    void FillEnergyHist(const double &e);
    void FillEnergyHist(const EventData &event);
    void ResetEnergyHist();
    TH1 *GetEnergyHist() const {return energy_hist;}
    void SaveHists(const std::string &path) const;
    std::vector<double> FitHist(const std::string &channel,
                                const std::string &hist_name,
                                const std::string &fit_function,
                                const double &range_min,
                                const double &range_max,
                                const bool &verbose) const;
    void FitPedestal();
    void CorrectGainFactor(int ref);

private:
    PRadHyCalDetector *hycal;
    PRadHyCalReconstructor recon;
    TH1D *energy_hist;

    std::vector<CalPeriod> cal_period;

    // channel lists
    std::vector<PRadADCChannel*> adc_list;
    std::vector<PRadTDCChannel*> tdc_list;

    // channel maps
    std::unordered_map<ChannelAddress, PRadADCChannel*> adc_addr_map;
    std::unordered_map<std::string, PRadADCChannel*> adc_name_map;
    std::unordered_map<ChannelAddress, PRadTDCChannel*> tdc_addr_map;
    std::unordered_map<std::string, PRadTDCChannel*> tdc_name_map;
};

#endif

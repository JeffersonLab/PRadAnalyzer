#ifndef PRAD_HYCAL_SYSTEM_H
#define PRAD_HYCAL_SYSTEM_H

#include <string>
#include <ostream>
#include <unordered_map>
#include "PRadEventStruct.h"
#include "PRadHyCalDetector.h"
#include "PRadHyCalCluster.h"
#include "PRadSquareCluster.h"
#include "PRadIslandCluster.h"
#include "PRadClusterProfile.h"
#include "PRadTDCChannel.h"
#include "PRadADCChannel.h"
#include "PRadCalibConst.h"
#include "ConfigObject.h"

#ifdef USE_PRIMEX_METHOD
#include "PRadPrimexCluster.h"
#endif

// adc searching speed is important, thus reserve buckets to have unordered_map
// better formed
#define ADC_BUCKETS 2000

// data structure for finding the calibration period
struct CalPeriod
{
    int begin;
    int end;
    int main;
    int sub;

    CalPeriod()
    : begin(0), end(0), main(0), sub(0)
    {};
    CalPeriod(int b, int e, int m, int s)
    : begin(b), end(e), main(m), sub(s)
    {};

    bool operator ==(const int &run)
    const
    {
        return (run >= begin) && (run <= end);
    }

    bool operator !=(const int &run)
    const
    {
        return (run < begin) || (run > end);
    }

    bool operator <(const int &run)
    const
    {
        return end < run;
    }

    bool operator >(const int &run)
    const
    {
        return begin > run;
    }
};

class TH1D;

class PRadHyCalSystem : public ConfigObject
{
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
    void ReadChannelList(const std::string &path);
    void ReadRunInfoFile(const std::string &path);
    void ReadTriggerEffFile(const std::string &path);
    void ReadCalPeriodFile(const std::string &path);
    void ChooseRun(const std::string &path, bool verbose = true);
    void ChooseRun(int run, bool verbose = true);
    void UpdateRunFiles(bool verbose = true);

    // connections
    void BuildConnections();

    // events related
    void ChooseEvent(const EventData &data);
    void Reconstruct();
    void Reconstruct(const EventData &data);
    void Reset();

    // detector related
    void SetDetector(PRadHyCalDetector *h);
    void RemoveDetector();
    void DisconnectDetector(bool force_disconn = false);
    double GetEnergy(const EventData &event) const;
    PRadHyCalModule *GetModule(const int &id) const;
    PRadHyCalModule *GetModule(const std::string &name) const;
    std::vector<PRadHyCalModule*> GetModuleList() const;
    PRadHyCalDetector *GetDetector() const {return hycal;};

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
    const std::vector<PRadADCChannel*> &GetADCList() const {return adc_list;};
    const std::vector<PRadTDCChannel*> &GetTDCList() const {return tdc_list;};
    void Sparsify(const EventData &event);

    // clustering method related
    bool AddClusterMethod(const std::string &name, PRadHyCalCluster *c);
    void RemoveClusterMethod(const std::string &name);
    void ClearClusterMethods();
    void SetClusterMethod(const std::string &name);
    PRadHyCalCluster *GetClusterMethod() const {return recon;};
    PRadHyCalCluster *GetClusterMethod(const std::string &name) const;
    std::string GetClusterMethodName() const;
    std::vector<std::string> GetClusterMethodNames() const;

    // histogram related
    void FillHists(const EventData &event);
    void FillEnergyHist();
    void FillEnergyHist(const double &e);
    void FillEnergyHist(const EventData &event);
    void ResetEnergyHist();
    TH1 *GetEnergyHist() const {return energy_hist;};
    void SaveHists(const std::string &path) const;
    std::vector<double> FitHist(const std::string &channel,
                                const std::string &hist_name,
                                const std::string &fit_function,
                                const double &range_min,
                                const double &range_max,
                                const bool &verbose) const throw(PRadException);
    void FitPedestal();
    void CorrectGainFactor(int ref);

private:
    PRadHyCalDetector *hycal;
    PRadHyCalCluster *recon;
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

    // clustering method map
    std::unordered_map<std::string, PRadHyCalCluster*> recon_map;
};

#endif

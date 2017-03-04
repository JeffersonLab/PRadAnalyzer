#ifndef PRAD_TDC_CHANNEL_H
#define PRAD_TDC_CHANNEL_H

#include <vector>
#include <string>
#include <unordered_map>
#include "PRadDAQChannel.h"

class PRadADCChannel;
class TH1;

class PRadTDCChannel : public PRadDAQChannel
{
public:
    // constructor
    PRadTDCChannel(const std::string &n, const ChannelAddress &addr);

    // copy/move constructors
    PRadTDCChannel(const PRadTDCChannel &that);
    PRadTDCChannel(PRadTDCChannel &&that);

    // destructor
    virtual ~PRadTDCChannel();

    // copy/move assignment operators
    PRadTDCChannel &operator =(const PRadTDCChannel &rhs);
    PRadTDCChannel &operator =(PRadTDCChannel &&rhs);

    void ConnectChannel(PRadADCChannel *ch);
    void DisconnectChannel(int id, bool force_disconn = false);
    void DisconnectChannels();
    void AddTimeMeasure(const unsigned short &count);
    void AddTimeMeasure(const std::vector<unsigned short> &counts);
    void SetTimeMeasure(const std::vector<unsigned short> &counts);
    void FillHist(const unsigned short &time);
    void Reset();
    void ResetHists();
    void ClearTimeMeasure();

    PRadADCChannel* GetADCChannel(int id) const;
    TH1 *GetHist() const {return tdc_hist;};
    std::vector<PRadADCChannel*> GetChannelList() const;
    const std::vector<unsigned short> &GetTimeMeasure() const {return time_measure;};

private:
    std::unordered_map<int, PRadADCChannel*> group_map;
    std::vector<unsigned short> time_measure;
    TH1 *tdc_hist;
};

#endif

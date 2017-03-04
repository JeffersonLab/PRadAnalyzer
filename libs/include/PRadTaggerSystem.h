#ifndef PRAD_TAGGER_SYSTEM_H
#define PRAD_TAGGER_SYSTEM_H

#include "PRadEventStruct.h"
#include "datastruct.h"

#define TAGGER_CHANID 30000 // Tagger tdc id will start from this number
#define TAGGER_T_CHANID 1000 // Start from TAGGER_CHANID, more than 1000 will be t channel

class TH2I;

class PRadTaggerSystem
{
public:
    // constructor
    PRadTaggerSystem();

    // copy/move constructors
    PRadTaggerSystem(const PRadTaggerSystem &that);
    PRadTaggerSystem(PRadTaggerSystem &&that);

    // destructor
    virtual ~PRadTaggerSystem();

    // copy/move assignment operators
    PRadTaggerSystem &operator =(const PRadTaggerSystem &rhs);
    PRadTaggerSystem &operator =(PRadTaggerSystem &&rhs);

    void Reset();

    // fill hists
    void FeedTaggerHits(const TDCV1190Data &data, EventData &event);
    void FillHists(const EventData &event);

    // get hists
    TH2I *GetECounterHist() const {return hist_E;};
    TH2I *GetTCounterHist() const {return hist_T;};

private:
    TH2I *hist_E;
    TH2I *hist_T;
};

#endif

#ifndef PRAD_DATA_HANDLER_H
#define PRAD_DATA_HANDLER_H

#include <deque>
#include <thread>
#include <unordered_map>
#include "PRadEvioParser.h"
#include "PRadDSTParser.h"
#include "PRadEventStruct.h"
#include "PRadException.h"

// PMT 0 - 2
#define DEFAULT_REF_PMT 2

class PRadHyCalSystem;
class PRadGEMSystem;
class PRadEPICSystem;
class PRadTaggerSystem;
class TH2I;

class PRadDataHandler
{
public:
    // constructor
    PRadDataHandler();

    // copy/move constructors
    PRadDataHandler(const PRadDataHandler &that);
    PRadDataHandler(PRadDataHandler &&that);

    // destructor
    virtual ~PRadDataHandler();

    // copy/move assignment operators
    PRadDataHandler &operator =(const PRadDataHandler &rhs);
    PRadDataHandler &operator =(PRadDataHandler &&rhs);

    // mode change
    void SetOnlineMode(const bool &mode);

    // set systems
    void SetHyCalSystem(PRadHyCalSystem *hycal) {hycal_sys = hycal;};
    void SetGEMSystem(PRadGEMSystem *gem) {gem_sys = gem;};
    void SetEPICSystem(PRadEPICSystem *epics) {epic_sys = epics;};
    void SetTaggerSystem(PRadTaggerSystem *tagger) {tagger_sys = tagger;};
    PRadHyCalSystem *GetHyCalSystem() const {return hycal_sys;};
    PRadGEMSystem *GetGEMSystem() const {return gem_sys;};
    PRadEPICSystem *GetEPICSystem() const {return epic_sys;};

    // file reading and writing
    void Decode(const void *buffer);
    void ReadFromDST(const std::string &path, unsigned int mode = 0);
    void ReadFromEvio(const std::string &path, int evt = -1, bool verbose = false);
    void ReadFromSplitEvio(const std::string &path, int split = -1, bool verbose = true);
    void WriteToDST(const std::string &path);
    void Replay(const std::string &r_path, int split = -1, const std::string &w_path = "");

    // data handler
    void Clear();
    void StartofNewEvent(const unsigned char &tag);
    void EndofThisEvent(const unsigned int &ev);
    void EndProcess(EventData *data);
    void FillHistograms(const EventData &data);
    void UpdateTrgType(const unsigned char &trg);


    // feeding data
    void FeedData(const JLabTIData &tiData);
    void FeedData(const JLabDSCData &dscData);
    void FeedData(const ADC1881MData &adcData);
    void FeedData(const TDCV767Data &tdcData);
    void FeedData(const TDCV1190Data &tdcData);
    void FeedData(const GEMRawData &gemData);
    void FeedData(const std::vector<GEMZeroSupData> &gemData);
    void FeedData(const EPICSRawData &epicsData);


    // show data
    void ChooseEvent(const int &idx = -1);
    void ChooseEvent(const EventData &event);
    int GetCurrentEventNb() const {return current_event;};
    unsigned int GetEventCount() const {return event_data.size();};
    const EventData &GetEvent(const unsigned int &index) const throw (PRadException);
    const std::deque<EventData> &GetEventData() const {return event_data;};

    // analysis tools
    void InitializeByData(const std::string &path = "", int ref = DEFAULT_REF_PMT);
    void RefillEnergyHist();
    int FindEvent(int event_number) const;

private:
    void waitEventProcess();

private:
    PRadEvioParser parser;
    PRadDSTParser dst_parser;
    PRadEPICSystem *epic_sys;
    PRadTaggerSystem *tagger_sys;
    PRadHyCalSystem *hycal_sys;
    PRadGEMSystem *gem_sys;
    bool onlineMode;
    bool replayMode;
    int current_event;
    std::thread end_thread;

    // data related
    std::deque<EventData> event_data;
    EventData *new_event;
    EventData *proc_event;
};

#endif

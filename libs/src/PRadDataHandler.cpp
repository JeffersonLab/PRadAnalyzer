//============================================================================//
// The data handler and container class                                       //
// Dealing with the data from all the channels                                //
//                                                                            //
// Chao Peng                                                                  //
// 02/07/2016                                                                 //
//============================================================================//

#include <iostream>
#include <iomanip>
#include <algorithm>
#include "PRadDataHandler.h"
#include "PRadInfoCenter.h"
#include "PRadEPICSystem.h"
#include "PRadTaggerSystem.h"
#include "PRadHyCalSystem.h"
#include "PRadGEMSystem.h"
#include "PRadBenchMark.h"
#include "ConfigParser.h"
#include "canalib.h"
#include "TH2.h"



//============================================================================//
// Constructors, Destructor, Assignment Operators                             //
//============================================================================//

// constructor
PRadDataHandler::PRadDataHandler()
: parser(this), dst_parser(this),
  epic_sys(nullptr), tagger_sys(nullptr), hycal_sys(nullptr), gem_sys(nullptr),
  onlineMode(false), replayMode(false), current_event(0),
  new_event(new EventData), proc_event(new EventData)
{
    // place holder
}

// copy/move constructors
PRadDataHandler::PRadDataHandler(const PRadDataHandler &that)
: parser(this), dst_parser(this),
  epic_sys(nullptr), tagger_sys(nullptr), hycal_sys(nullptr), gem_sys(nullptr),
  onlineMode(that.onlineMode), replayMode(that.replayMode),
  current_event(that.current_event), event_data(that.event_data),
  new_event(new EventData(*that.new_event)), proc_event(new EventData(*that.proc_event))
{
    // place holder
}

PRadDataHandler::PRadDataHandler(PRadDataHandler &&that)
: parser(this), dst_parser(this),
  epic_sys(nullptr), tagger_sys(nullptr), hycal_sys(nullptr), gem_sys(nullptr),
  onlineMode(that.onlineMode), replayMode(that.replayMode),
  current_event(that.current_event), event_data(std::move(that.event_data)),
  new_event(new EventData(std::move(*that.new_event))),
  proc_event(new EventData(std::move(*that.proc_event)))
{
    // place holder
}

// destructor
PRadDataHandler::~PRadDataHandler()
{
    delete new_event;
    delete proc_event;
}

// copy/move assignment operators
PRadDataHandler &PRadDataHandler::operator =(const PRadDataHandler &rhs)
{
    if(this == &rhs)
        return *this;

    PRadDataHandler that(rhs);
    *this = std::move(that);
    return *this;
}

PRadDataHandler &PRadDataHandler::operator =(PRadDataHandler &&rhs)
{
    if(this == &rhs)
        return *this;

    delete new_event;
    delete proc_event;

    new_event = rhs.new_event;
    rhs.new_event = nullptr;
    proc_event = rhs.proc_event;
    rhs.proc_event = nullptr;
    onlineMode = rhs.onlineMode;
    replayMode = rhs.replayMode;
    current_event = rhs.current_event;
    event_data = std::move(rhs.event_data);

    return *this;
}

//============================================================================//
// Public Member Functions                                                    //
//============================================================================//

// change online mode
void PRadDataHandler::SetOnlineMode(const bool &mode)
{
    onlineMode = mode;
    Clear();
}

// decode an event buffer
void PRadDataHandler::Decode(const void *buffer)
{
    parser.ReadEventBuffer(buffer);

    waitEventProcess();
}

// read from DST format file
void PRadDataHandler::ReadFromDST(const std::string &path, unsigned int mode)
{
    try {
        dst_parser.OpenInput(path);

        dst_parser.SetMode(mode);

        std::cout << "Data Handler: Reading events from DST file "
                  << "\"" << path << "\""
                  << std::endl;

        while(dst_parser.Read())
        {
            switch(dst_parser.EventType())
            {
            case PRadDSTParser::Type::event:
                // fill histogram
                FillHistograms(dst_parser.GetEvent());
                // count occupancy
                if(hycal_sys)
                    hycal_sys->Sparsify(dst_parser.GetEvent());
                // save data
                event_data.emplace_back(std::move(dst_parser.event));
                break;
            case PRadDSTParser::Type::epics:
                if(epic_sys)
                    epic_sys->AddEvent(std::move(dst_parser.epics_event));
                break;
            default:
                break;
            }
        }

    } catch(PRadException &e) {
        std::cerr << e.FailureType() << ": "
                  << e.FailureDesc() << std::endl
                  << "Read from DST Aborted!" << std::endl;
    } catch(std::exception &e) {
        std::cerr << e.what() << std::endl
                  << "Read from DST Aborted!" << std::endl;
    }
    dst_parser.CloseInput();
 }


// read fro evio file
int PRadDataHandler::ReadFromEvio(const std::string &path, int evt, bool verbose)
{
    int count = parser.ReadEvioFile(path.c_str(), evt, verbose);
    waitEventProcess();
    return count;
}

// read from splitted evio file
int PRadDataHandler::ReadFromSplitEvio(const std::string &path, int split, bool verbose)
{
    if(split < 0) {// default input, no split
        return ReadFromEvio(path.c_str(), -1, verbose);
    } else {
        int count = 0;
        for(int i = 0; i <= split; ++i)
        {
            std::string split_path = path + "." + std::to_string(i);
            count += ReadFromEvio(split_path.c_str(), -1, verbose);
        }
        return count;
    }
}

// erase the data container and all the connected systems
void PRadDataHandler::Clear()
{
    // used memory won't be released, but it can be used again for new data file
    event_data = std::deque<EventData>();
    parser.SetEventNumber(0);

    PRadInfoCenter::Instance().Reset();

    if(epic_sys)
        epic_sys->Reset();

    if(tagger_sys)
        tagger_sys->Reset();

    if(hycal_sys)
        hycal_sys->Reset();

    if(gem_sys)
        gem_sys->Reset();
}

// signal of new event
void PRadDataHandler::StartofNewEvent(const unsigned char &tag)
{
    new_event->update_type(tag);
}

// update trigger type
void PRadDataHandler::UpdateTrgType(const unsigned char &trg)
{
    if(new_event->trigger && (new_event->trigger != trg)) {
        std::cerr << "ERROR: Trigger type mismatch at event "
                  << parser.GetEventNumber()
                  << ", was " << (int) new_event->trigger
                  << " now " << (int) trg
                  << std::endl;
    }
    new_event->trigger = trg;
}

// feed JLab TI Data
void PRadDataHandler::FeedData(const JLabTIData &tiData)
{
    new_event->timestamp = tiData.time_high;
    new_event->timestamp <<= 32;
    new_event->timestamp |= tiData.time_low;
}

// feed JLab discriminator data
void PRadDataHandler::FeedData(const JLabDSCData &dscData)
{
    for(uint32_t i = 0; i < dscData.size; ++i)
    {
        new_event->dsc_data.emplace_back(dscData.gated_buf[i], dscData.ungated_buf[i]);
    }
}

// feed ADC1881M data
void PRadDataHandler::FeedData(const ADC1881MData &adcData)
{
    if(!hycal_sys)
        return;

    // get the channel
    PRadADCChannel *channel = hycal_sys->GetADCChannel(adcData.addr);

    if(!channel)
        return;

    if(new_event->is_physics_event()) {
        if(channel->Sparsify(adcData.val)) {
            new_event->add_adc(ADC_Data(channel->GetID(), adcData.val)); // store this data word
        }
    } else if (new_event->is_monitor_event()) {
        new_event->add_adc(ADC_Data(channel->GetID(), adcData.val));
    }

}

// feed TDC CAEN v767 data
void PRadDataHandler::FeedData(const TDCV767Data &tdcData)
{
    if(!hycal_sys)
        return;

    PRadTDCChannel *tdc = hycal_sys->GetTDCChannel(tdcData.addr);

    if(!tdc)
        return;

    new_event->tdc_data.push_back(TDC_Data(tdc->GetID(), tdcData.val));
}

// feed TDC CAEN v1190 data
void PRadDataHandler::FeedData(const TDCV1190Data &tdcData)
{
    if(!hycal_sys)
        return;

    // tagger hits
    if(tdcData.addr.crate == PRadTagE) {
        if(tagger_sys)
            tagger_sys->FeedTaggerHits(tdcData, *new_event);
        return;
    }

    PRadTDCChannel *tdc = hycal_sys->GetTDCChannel(tdcData.addr);

    if(!tdc)
        return;

    new_event->add_tdc(TDC_Data(tdc->GetID(), tdcData.val));
}

// feed GEM data
void PRadDataHandler::FeedData(const GEMRawData &gemData)
{
    if(gem_sys)
        gem_sys->FillRawData(gemData, *new_event);
}

// feed GEM data which has been zero-suppressed
void PRadDataHandler::FeedData(const std::vector<GEMZeroSupData> &gemData)
{
    if(gem_sys)
        gem_sys->FillZeroSupData(gemData, *new_event);
}

// feed EPICS data
void PRadDataHandler::FeedData(const EPICSRawData &epicsData)
{
    if(epic_sys)
        epic_sys->FillRawData(epicsData.buf);
}

// Fill the event information into histograms
void PRadDataHandler::FillHistograms(const EventData &data)
{
    if(hycal_sys)
        hycal_sys->FillHists(data);

    if(tagger_sys)
        tagger_sys->FillHists(data);
}

// signal of event end, save event or discard event in online mode
void PRadDataHandler::EndofThisEvent(const unsigned int &ev)
{
    new_event->event_number = ev;
    // wait for the process thread
    waitEventProcess();

    // swap the pointers
    EventData *tmp = new_event;
    new_event = proc_event;
    proc_event = tmp;

    end_thread = std::thread(&PRadDataHandler::EndProcess, this, proc_event);
}

// wait for the end process to finish
inline void PRadDataHandler::waitEventProcess()
{
    if(end_thread.joinable())
        end_thread.join();
}

void PRadDataHandler::EndProcess(EventData *ev)
{
    if(ev->get_type() == EPICS_Info) {

        if(epic_sys) {
            if(replayMode)
                dst_parser.WriteEPICS(EpicsData(ev->event_number, epic_sys->GetCurrentValues()));
            else
                epic_sys->SaveData(ev->event_number, onlineMode);
        }

    } else { // event or sync event

        FillHistograms(*ev);
        PRadInfoCenter::Instance().UpdateInfo(*ev);

        // online mode only saves the last event, to reduce usage of memory
        if(onlineMode && event_data.size())
            event_data.pop_front();

        if(replayMode)
            dst_parser.WriteEvent(*ev);
        else
            event_data.emplace_back(std::move(*ev)); // save event

    }

    // clear the event for the future usage
    ev->clear();
}

// show the event to event viewer
void PRadDataHandler::ChooseEvent(const int &idx)
{
    if (event_data.size()) { // offline mode, pick the event given by console
        if((unsigned int) idx >= event_data.size())
            ChooseEvent(event_data.back());
        else
            ChooseEvent(event_data.at(idx));
    }
}

void PRadDataHandler::ChooseEvent(const EventData &event)
{
    if(hycal_sys)
        hycal_sys->ChooseEvent(event);
    if(gem_sys)
        gem_sys->ChooseEvent(event);

    current_event = event.event_number;
}

// get the event by index
const EventData &PRadDataHandler::GetEvent(const unsigned int &index)
const
throw (PRadException)
{
    if(!event_data.size())
        throw PRadException("PRad Data Handler Error", "Empty data bank!");

    if(index >= event_data.size()) {
        return event_data.back();
    } else {
        return event_data.at(index);
    }
}

// Refill energy hist after correct gain factos
void PRadDataHandler::RefillEnergyHist()
{
    if(!hycal_sys)
        return;

    hycal_sys->ResetEnergyHist();

    for(auto &event : event_data)
    {
        if(!event.is_physics_event())
            continue;

        hycal_sys->FillEnergyHist(event);
    }
}

// try to get all the needed information from monitor events
void PRadDataHandler::InitializeByData(const std::string &path, int ref)
{
    PRadBenchMark timer;

    if(!hycal_sys || !gem_sys) {
        std::cerr << "Data Handler: HyCal System or GEM System missing, "
                  << "abort initialization."
                  << std::endl;
        return;
    }

    std::cout << "Data Handler: Initializing from Data File "
              << "\"" << path << "\"."
              << std::endl;

    if(!path.empty()) {
        PRadInfoCenter::SetRunNumber(path);
        gem_sys->SetPedestalMode(true);
        parser.ReadEvioFile(path.c_str(), 20000);
    }

    std::cout << "Data Handler: Fitting Pedestal for HyCal." << std::endl;
    hycal_sys->FitPedestal();

    std::cout << "Data Handler: Correct HyCal Gain Factor. " << std::endl;
    hycal_sys->CorrectGainFactor(ref);

    std::cout << "Data Handler: Fitting Pedestal for GEM." << std::endl;
    gem_sys->FitPedestal();

    std::cout << "Data Handler: Releasing Memeory." << std::endl;
    gem_sys->SetPedestalMode(false);

    // save run number
    int run_number = PRadInfoCenter::GetRunNumber();
    Clear();
    PRadInfoCenter::SetRunNumber(run_number);

    std::cout << "Data Handler: Done initialization, took "
              << timer.GetElapsedTime()/1000. << " s"
              << std::endl;
}

// find event by its event number
// it is assumed the files decoded are all from 1 single run and they are loaded in order
// otherwise this function will not work properly
int PRadDataHandler::FindEvent(int evt)
const
{
    auto it = cana::binary_search(event_data.begin(), event_data.end(), evt);

    if(it == event_data.end())
        return -1;

    return it - event_data.begin();
}

// replay the raw data file, do zero suppression and save it in DST format
void PRadDataHandler::Replay(const std::string &r_path, int split, const std::string &w_path)
{
    if(w_path.empty()) {
        std::string file = "prad_" + std::to_string(PRadInfoCenter::GetRunNumber()) + ".dst";
        dst_parser.OpenOutput(file);
    } else {
        dst_parser.OpenOutput(w_path);
    }

    std::cout << "Replay started!" << std::endl;
    PRadBenchMark timer;

    dst_parser.WriteHyCalInfo(hycal_sys);
    dst_parser.WriteGEMInfo(gem_sys);
    dst_parser.WriteEPICSMap(epic_sys);

    replayMode = true;

    int count = ReadFromSplitEvio(r_path, split);

    dst_parser.WriteRunInfo();

    replayMode = false;

    std::cout << "Replay done, took "
              << timer.GetElapsedTime()/1000. << " s! "
              << "Replayed " << count << " events."
              << std::endl;
    dst_parser.CloseOutput();
}

// write the current data bank to DST file
void PRadDataHandler::WriteToDST(const std::string &path)
{
    try {
        dst_parser.OpenOutput(path);

        std::cout << "Data Handler: Saving DST file "
                  << "\"" << path << "\""
                  << std::endl;

        if(hycal_sys)
            dst_parser.WriteHyCalInfo(hycal_sys);
        if(gem_sys)
            dst_parser.WriteGEMInfo(gem_sys);
        if(epic_sys)
            dst_parser.WriteEPICSMap(epic_sys);

        if(epic_sys) {
            for(auto &epics : epic_sys->GetEventData())
            {
                dst_parser.WriteEPICS(epics);
            }
        }

        for(auto &event : event_data)
        {
            dst_parser.WriteEvent(event);
        }

        dst_parser.WriteRunInfo();

    } catch(PRadException &e) {
        std::cerr << e.FailureType() << ": "
                  << e.FailureDesc() << std::endl
                  << "Write to DST Aborted!" << std::endl;
    } catch(std::exception &e) {
        std::cerr << e.what() << std::endl
                  << "Write to DST Aborted!" << std::endl;
    }

    dst_parser.CloseOutput();
}

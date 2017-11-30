//============================================================================//
// A singleton class that stores and manages the run information or online    //
// information, it can be shared through all the classes                      //
// TODO merge live_scaled_charge into RunInfo                                 //
// It may change existing DST file parsing, so need a careful treatment       //
//                                                                            //
// Chao Peng                                                                  //
// 11/19/2016                                                                 //
//============================================================================//

#include "PRadInfoCenter.h"
#include "ConfigParser.h"
#include <iostream>

// add the trigger channels
PRadInfoCenter::PRadInfoCenter()
: live_scaled_charge(0)
{
    online_info.add_trigger("Lead Glass Sum", 0);
    online_info.add_trigger("Total Sum", 1);
    online_info.add_trigger("LMS Led", 2);
    online_info.add_trigger("LMS Alpha Source", 3);
    online_info.add_trigger("Tagger Master OR", 4);
    online_info.add_trigger("Scintillator", 5);
}

// clear all the information
void PRadInfoCenter::Reset()
{
    run_info.reset();
    online_info.reset();
    live_scaled_charge = 0.;
}

// update information from event data
void PRadInfoCenter::UpdateInfo(const EventData &event)
{
    // only synchronization event contains the following information
    // this is by the design of our DAQ system
    if(event.get_type() != CODA_Sync)
        return;

    // online information update, update to the latest values
    // update triggers
    for(auto trg_ch : online_info.trigger_info)
    {
        if(trg_ch.id < event.dsc_data.size())
        {
            // get ungated trigger counts
            unsigned int counts = event.get_dsc_channel(trg_ch.id).ungated_count;

            // calculate the frequency
            trg_ch.freq = (double)counts / event.get_beam_time();
        } else {
            std::cerr << "PRad Info Center Erro: Unmatched discriminator data"
                      << " from event " << event.event_number
                      << ", expect trigger " << trg_ch.name
                      << " at channel " << trg_ch.id
                      << ", but the event only has " << event.dsc_data.size()
                      << " dsc channels."
                      << std::endl;
        }
    }

    // update live time
    online_info.live_time = event.get_live_time();

    //update beam current
    online_info.beam_current = event.get_beam_current();

    // run information, accu
    // only collect run information for physics events
    if(!event.is_physics_event())
        return;

    double beam_charge = event.get_beam_charge();
    unsigned int dead_count = event.get_ref_channel().gated_count;
    unsigned int total_count = event.get_ref_channel().ungated_count;
    run_info.beam_charge += beam_charge;
    run_info.ungated_count += total_count;
    run_info.dead_count += dead_count;
    live_scaled_charge += beam_charge * (1. - (double)dead_count/(double)total_count);
}

// set run number
bool PRadInfoCenter::SetRunNumber(int run)
{
    if(Instance().run_info.run_number != run) {
        Instance().run_info.run_number = run;
        return true;
    }

    return false;
}

// get run number
int PRadInfoCenter::GetRunNumber()
{
    return Instance().run_info.run_number;
}

// get beam charge
double PRadInfoCenter::GetBeamCharge()
{
    return Instance().run_info.beam_charge;
}

// get beam charge scaled by live time
double PRadInfoCenter::GetLiveBeamCharge()
{
    return Instance().live_scaled_charge;
}

// get live time
double PRadInfoCenter::GetLiveTime()
{
    const RunInfo &run = Instance().run_info;
    if(!run.ungated_count)
        return 0.;
    return (1. - (double)run.dead_count/(double)run.ungated_count);
};

// set run number from file path
bool PRadInfoCenter::SetRunNumber(const std::string &path)
{
    std::string file_name = ConfigParser::decompose_path(path).name;
    int run = ConfigParser::find_integer(file_name);

    if(run > 0 && Instance().run_info.run_number != run) {
        Instance().run_info.run_number = run;
        return true;
    }

    return false;
}

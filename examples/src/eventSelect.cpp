//============================================================================//
// An example showing how to cut off events from bad events list              //
//                                                                            //
// Chao Peng                                                                  //
// 11/12/2016                                                                 //
//============================================================================//

#include "PRadDSTParser.h"
#include "PRadEventFilter.h"
#include "PRadBenchMark.h"
#include "ConfigParser.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#define PROGRESS_COUNT 10000

using namespace std;

void EventSelect(const string &file, const string &bad_file);
void WriteEvents(PRadDSTParser &dst, const PRadEventFilter &filter, const vector<EventData> &data);

// expecting two input string, event file path and bad events list path
int main(int argc, char *argv[])
{
    if(argc != 3) {
        cout << "usage: eventSelect <file> <bad_event_list>" << endl;
        return 0;
    }

    string file = argv[1];
    string bad_file = argv[2];

    EventSelect(file, bad_file);
}

void EventSelect(const string &file, const string &bad_file)
{
    cout << "Event selection of " << "\"" << file << "\"."
         << endl
         << "Based on the list of " << "\"" << bad_file << "\"."
         << endl;

    PRadEventFilter filter(bad_file);
    PRadDSTParser dst_parser, dst_parser2;

    dst_parser.OpenInput(file);
    auto path_info = ConfigParser::decompose_path(file);
    string fname = path_info.name;
    string dir = path_info.dir;

    dst_parser.OpenOutput(fname + "_sel.dst");
    dst_parser2.OpenOutput(fname + "_mon.dst");

    PRadBenchMark timer;
    double time = 0.;

    // a package for a period (from sync event to sync event
    vector<EventData> event_pack;
    event_pack.reserve(150000);

    unsigned int m_count = 0, p_count = 0, t_count = 0;
    bool first_sync = true;
    while(dst_parser.Read())
    {
        if(dst_parser.EventType() == PRadDSTParser::Type::event) {
            auto event = dst_parser.GetEvent();
            t_count++;

            // write the monitoring event
            if(event.is_monitor_event()) {
                m_count++;
                dst_parser2.WriteEvent(event);
                continue;
            }

            // throw away everything before the first sync
            if(!first_sync)
                event_pack.push_back(event);

            // sync event is the beginning of this package
            if(event.is_sync_event()) {
                // don't save the events before the first sync
                if(first_sync) {
                    first_sync = false;
                    timer.Reset();
                    cout << "Found the first sync event at "
                         << event.event_number
                         << ", start recording."
                         << endl;
                } else {
                    WriteEvents(dst_parser, filter, event_pack);

                    // show progress
                    double t = timer.GetElapsedTime();
                    time += t;
                    timer.Reset();
                    p_count += event_pack.size();

                    cout << "---[ ev " << p_count << " ]---"
                         << "---[ " << t << " ms ]---"
                         << "\r" << flush;
                }

                // clear event package after writing
                event_pack.clear();
            }
        } else if(dst_parser.EventType() == PRadDSTParser::Type::epics) {
            dst_parser.WriteEPICS();
        }
    }

    time += timer.GetElapsedTime();
    dst_parser.CloseInput();
    dst_parser.CloseOutput();
    dst_parser2.CloseOutput();

    cout << endl;
    cout << "TIMER: Finished event selection, select "
         << m_count << " monitor events and "
         << p_count << " physics events from "
         << t_count << " events, took "
         << time/1000. << " s."
         << endl;
}

void WriteEvents(PRadDSTParser &dst, const PRadEventFilter &filter, const vector<EventData> &data)
{
    if(data.empty())
        return;

    if(filter.IsBadPeriod(data.front(), data.back()))
        return;

    for(auto &event : data)
    {
        dst.WriteEvent(event);
    }
}

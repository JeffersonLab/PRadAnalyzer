//============================================================================//
// An example showing how to get beam charge information from one run         //
//                                                                            //
// Chao Peng                                                                  //
// 11/12/2016                                                                 //
//============================================================================//

#include "PRadDSTParser.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

using namespace std;

void showLiveTime(const string &file);
inline double liveTime(double gated, double ungated)
{
    return 1. - gated/ungated;
}

int main(int argc, char *argv[])
{
    if(argc != 2)
    {
        cout << "usage: " << argv[0] << " <file>" << endl;
        return -1;
    }

    showLiveTime(argv[1]);

    return 0;
}

void showLiveTime(const string &file)
{
    PRadDSTParser dst_parser;
    dst_parser.OpenInput(file);

    cout << "Live time for " << "\"" << file << "\"." << endl;

    int total_sum_cnt = 0, lg_sum_cnt = 0;
    int previous_period = 0;
    vector<int> counts_total(6, 0), counts_period(6, 0);

    while(dst_parser.Read())
    {
        // only interested in events
        if(dst_parser.EventType() != PRadDSTParser::Type::event)
            continue;

        auto event = dst_parser.GetEvent();

        // only interested in physics events
        if(!event.is_physics_event())
            continue;

        // count trigger types
        switch(event.trigger)
        {
        case PHYS_LeadGlassSum: lg_sum_cnt++; break;
        case PHYS_TotalSum: total_sum_cnt++; break;
        default: break;
        }


        // show gated counter and ungated counter for triggers
        if(event.is_sync_event()) {
            counts_period[0] = event.get_trg_channel(PHYS_LeadGlassSum).gated_count;
            counts_period[1] = event.get_trg_channel(PHYS_LeadGlassSum).ungated_count;
            counts_period[2] = event.get_trg_channel(PHYS_TotalSum).gated_count;
            counts_period[3] = event.get_trg_channel(PHYS_TotalSum).ungated_count;
            counts_period[4] = event.get_ref_channel().gated_count;
            counts_period[5] = event.get_ref_channel().ungated_count;

            for(size_t i = 0; i < counts_total.size(); ++i)
            {
                counts_total[i] += counts_period[i];
            }

            std::cout << "======================================================\n"
                      << "Current period: " << previous_period << " - " << event.event_number << "\n"
                      << "TRG(LG_SUM) = " << lg_sum_cnt << "\n"
                      << "TRG(TOTAL_SUM) = " << total_sum_cnt << "\n"
                      << "LT(LG_SUM) = " << liveTime(counts_period[0], counts_period[1]) << "\n"
                      << "LT(TOTAL_SUM) = " << liveTime(counts_period[2], counts_period[3]) << "\n"
                      << "LT(REF_PULSER) = " << liveTime(counts_period[4], counts_period[5])
                      << std::endl;

            lg_sum_cnt = 0;
            total_sum_cnt = 0;
            previous_period = event.event_number;
        }

    }

    std::cout << "======================================================\n"
              << "Overal livetime: \n"
              << "LT(LG_SUM) = " << liveTime(counts_total[0], counts_total[1]) << "\n"
              << "LT(TOTAL_SUM) = " << liveTime(counts_total[2], counts_total[3]) << "\n"
              << "LT(REF_PULSER) = " << liveTime(counts_total[4], counts_total[5]) << std::endl;

}

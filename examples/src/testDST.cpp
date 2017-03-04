//============================================================================//
// An example showing how to use DST Parser to read and save selected events  //
//                                                                            //
// Chao Peng                                                                  //
// 11/20/2016                                                                 //
//============================================================================//

#include "PRadDSTParser.h"
#include "PRadHyCalSystem.h"
#include "PRadGEMSystem.h"
#include "PRadEPICSystem.h"
#include <iostream>
#include <string>
#include <vector>

using namespace std;

PRadEPICSystem *epics;
PRadHyCalSystem *hycal;
void SelectEvent(const string &file);

int main(int argc, char *argv[])
{
    if(argc < 2) {
        cout << "usage: testDST <file1> <file2> ..." << endl;
        return 0;
    }

    epics = new PRadEPICSystem("config/epics_channels.conf");
    hycal = new PRadHyCalSystem("config/hycal.conf");

    for(int i = 1; i < argc; ++i)
    {
        string file = argv[i];
        SelectEvent(file);
    }

    return 0;
}

void SelectEvent(const string &file)
{
    PRadDSTParser *dst_parser = new PRadDSTParser();
    dst_parser->OpenInput(file);
    dst_parser->OpenOutput(ConfigParser::decompose_path(file).name + "_select.dst");

    float beam_energy = 0.;
    int beam_energy_ch = epics->GetChannel("MBSY2C_energy");

    int count = 0;
    while(dst_parser->Read())
    {
        if(dst_parser->EventType() == PRadDSTParser::Type::event) {
            // you can push this event into data handler
            // handler->GetEventData().push_back(dst_parser->GetEvent()
            // or you can just do something with this event and discard it
            if(++count%10000 == 0)
                cout << count << "\r" << flush;

            auto event = dst_parser->GetEvent();
            if(!event.is_physics_event())
                continue;
            hycal->Reconstruct(event);
            auto &hits = hycal->GetDetector()->GetHits();

            if((hits.size() != 1) ||
               (!TEST_BIT(hits.front().flag, kTransition)))
                continue;

            float energy = 0;
            for(auto &hit : hits)
            {
                energy += hit.E;
            }

            // only save the event with the energy close to beam energy
            if(fabs((energy - beam_energy)/beam_energy) <= 0.3)
                dst_parser->WriteEvent();
        } else if (dst_parser->EventType() == PRadDSTParser::Type::epics) {
            auto epics_ev = dst_parser->GetEPICSEvent();
            dst_parser->WriteEPICS(epics_ev);
            beam_energy = epics_ev.values.at(beam_energy_ch);
        }
    }
    cout << endl;
    dst_parser->CloseInput();
    dst_parser->CloseOutput();
    delete dst_parser;
}


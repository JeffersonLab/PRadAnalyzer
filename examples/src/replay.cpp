//============================================================================//
// An application of replay raw data file and save the replayed data into DST //
// file. This is the 1st-level replay, it only discards the pedestal data     //
//                                                                            //
// Chao Peng                                                                  //
// 10/04/2016                                                                 //
//============================================================================//

#include "PRadDataHandler.h"
#include "PRadEPICSystem.h"
#include "PRadTaggerSystem.h"
#include "PRadHyCalSystem.h"
#include "PRadGEMSystem.h"
#include "PRadInfoCenter.h"
#include "PRadEvioParser.h"
#include "PRadBenchMark.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

using namespace std;

void print_instruction()
{
    cout << "usage: " << endl
         << setw(10) << "-i : " << "input file path" << endl
         << setw(10) << "-o : " << "output file path" << endl
         << setw(10) << "-h : " << "show options" << endl
         << endl;
}

int main(int argc, char * argv[])
{
    if(argc < 2) {
        print_instruction();
        return 0;
    }

    char *ptr;
    string output, input;

    // -i input_file -o output_file
    for(int i = 1; i < argc; ++i)
    {
        ptr = argv[i];
        if(*(ptr++) == '-') {
            switch(*(ptr++))
            {
            case 'o':
                output = argv[++i];
                break;
            case 'i':
                input = argv[++i];
                break;
            case 'h':
                print_instruction();
                break;
            default:
                cout << "Unkown option! check with -h" << endl;
                exit(1);
            }
        }
    }

    PRadDataHandler *handler = new PRadDataHandler();
    PRadEPICSystem *epics = new PRadEPICSystem("config/epics_channels.conf");
    PRadHyCalSystem *hycal = new PRadHyCalSystem("config/hycal.conf");
    PRadGEMSystem *gem = new PRadGEMSystem("config/gem.conf");
    PRadTaggerSystem *tagger = new PRadTaggerSystem;

    handler->SetEPICSystem(epics);
    handler->SetTaggerSystem(tagger);
    handler->SetHyCalSystem(hycal);
    handler->SetGEMSystem(gem);

    PRadBenchMark timer;
//    handler->ReadFromDST("/work/hallb/prad/replay/prad_001292.dst");
//    handler->ReadFromDST("prad_1292.dst");
//    handler->ReadFromEvio("/work/prad/xbai/1323/prad_001323.evio.1");
//    handler->ReadFromSplitEvio("/work/prad/xbai/1323/prad_001323.evio", 10);
    handler->InitializeByData(input+".0");
    handler->Replay(input, 1500, output);
//    handler->GetSRS()->SavePedestal("gem_ped.txt");


    cout << "TIMER: Finished, took " << timer.GetElapsedTime() << " ms" << endl;
    cout << "Read " << handler->GetEventCount() << " events and "
         << epics->GetEventCount() << " EPICS events from file."
         << endl;
    cout << PRadInfoCenter::GetBeamCharge() << endl;
    cout << PRadInfoCenter::GetLiveTime() << endl;

//    handler->WriteToDST("prad_001323_0-10.dst");
    //handler->PrintOutEPICS();
    return 0;
}


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
#include "ConfigOption.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

using namespace std;

void print_instruction()
{
    cout << "usage: replay <in_file> <out_file>" << endl
         << "\t" << "-s <value>: spliting file number, default -1 (no splitting)\n"
         << "\t" << "--init-evio: initialize from evio.0 file\n"
         << "\t" << "--init-database: initialize from database\n"
         << "\t" << "-r <value>: only has effect if --init-database is set, default -1 (determined from file name)\n"
         << "\t" << "-h " << "show options\n"
         << endl;
}

int main(int argc, char * argv[])
{
    ConfigOption conf_opt;
    conf_opt.AddOpt('s', ConfigOption::arg_require);
    conf_opt.AddOpt('h', ConfigOption::arg_none);
    conf_opt.AddOpt("init-evio", ConfigOption::arg_none, 'e');
    conf_opt.AddOpt("init-database", ConfigOption::arg_none, 'd');
    conf_opt.AddOpt('r', ConfigOption::arg_require);

    if(!conf_opt.ParseArgs(argc, argv)) {
        print_instruction();
        return -1;
    }

    bool evio_database = true;

    int split = -1, run = -1;
    for(auto &opt : conf_opt.GetOptions())
    {
        switch(opt.mark)
        {
        case 's':
            split = opt.var.Int();
            break;
        case 'e':
            evio_database = true;
            break;
        case 'd':
            evio_database = false;
            break;
        case 'r':
            run = opt.var.Int();
            break;
        case 'h':
        default:
            print_instruction();
            return -1;
        }
    }

    if(conf_opt.NbofArgs() != 2) {
        std::cerr << "Wrong number of arguments, require 2." << std::endl;
        return -1;
    }

    string input = conf_opt.GetArgument(0).String();
    string output = conf_opt.GetArgument(1).String();

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
    if(evio_database) {
        handler->InitializeByData(input + ".0");
    } else {
        if(run < 0)
            hycal->ChooseRun(input);
        else
            hycal->ChooseRun(run);
    }
    handler->Replay(input, split, output);


    cout << "TIMER: Finished, took " << timer.GetElapsedTime() << " ms" << endl;
    cout << PRadInfoCenter::GetBeamCharge() << endl;
    cout << PRadInfoCenter::GetLiveTime() << endl;

    return 0;
}


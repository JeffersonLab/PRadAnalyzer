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

int main(int argc, char * argv[])
{
    ConfigOption conf_opt;
    conf_opt.AddOpt(ConfigOption::arg_require, 's');
    conf_opt.AddLongOpt(ConfigOption::arg_none, "init-evio", 'e');
    conf_opt.AddLongOpt(ConfigOption::arg_none, "init-database", 'd');
    conf_opt.AddOpt(ConfigOption::arg_require, 'r');
    conf_opt.AddOpt(ConfigOption::help_message, 'h');

    conf_opt.SetDesc("usage: replay <in_file> <out_file>");
    conf_opt.SetDesc('s', "spliting file number, default -1 (no splitting).");
    conf_opt.SetDesc('r', "set run number, only valid for --init-database, default -1 (determined from file name).");
    conf_opt.SetDesc('e', "initialize from evio.0 file");
    conf_opt.SetDesc('d', "initialize from database");
    conf_opt.SetDesc('h', "show instruction.");

    if(!conf_opt.ParseArgs(argc, argv) || conf_opt.NbofArgs() != 2) {
        std::cout << conf_opt.GetInstruction() << std::endl;
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
        default:
            std::cout << conf_opt.GetInstruction() << std::endl;
            return -1;
        }
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


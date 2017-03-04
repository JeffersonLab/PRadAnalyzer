//============================================================================//
// An example showing how to get the cluster from HyCal and check how well it //
// could be described by cluster profile                                      //
//                                                                            //
// Chao Peng                                                                  //
// 11/12/2016                                                                 //
//============================================================================//

#include "cosmicEval.h"
#include "PRadHyCalSystem.h"
#include "PRadDSTParser.h"
#include "PRadBenchMark.h"
#include "TFile.h"
#include "TH1.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

using namespace std;

#define PROGRESS_COUNT 10000

double beam_energy = 1100.;

// declaration of functions
void Evaluate(const string &file);
PRadHyCalSystem *sys;

int main(int argc, char *argv[])
{
    if(argc < 2)
    {
        cout << "usage: cosmicCheck <file1> <file2> ..." << endl;
        return 0;
    }

    sys = new PRadHyCalSystem("config/hycal.conf");

    for(int i = 1; i < argc; ++i)
    {
        string file = argv[i];
        cout << "Analyzing File " << file << "..." << endl;
        Evaluate(file);          // a normal production run
    }

    return 0;
}

void Evaluate(const string &file)
{
    // remove directory and affix
    const string &fname = ConfigParser::decompose_path(file).name;

    PRadDSTParser dst_parser;
    dst_parser.OpenInput(file);
    dst_parser.OpenOutput(fname + "_sav.dst");
    PRadDSTParser dst_parser2;
    dst_parser2.OpenOutput(fname + "_rej.dst");

    TFile f((fname + "_prof.root").c_str(), "RECREATE");
    TH1I hist_size("Number of Clusters", "Number of Clusters", 50, 0, 50);
    TH1I hist_msize("Max Cluster Size", "Max Cluster Size", 50, 0, 50);
    TH1F hist_mene("Max Hit Energy", "Max Hit Energy", 1000, 0, 3000);
    TH1F hist_unif("Uniformity", "Energy Uniformity", 100, 0., 10.);
    TH1F hist_rsq("R Square", "R Square", 100, 0., 1.);
    TH1F hist_chisq("Chi Square", "Chi Square", 100, 0., 10.);

    int count = 0;
    PRadBenchMark timer;
    while(dst_parser.Read())
    {
        if(dst_parser.EventType() == PRadDSTParser::Type::event) {

            auto &event = dst_parser.GetEvent();

            // save sync event no matter what it is
            if(event.is_sync_event())
            {
                dst_parser.WriteEvent(event);
                continue;
            }

            // discard non-physics events
            if(!event.is_physics_event())
                continue;

            count++;
            if(count%PROGRESS_COUNT == 0) {
                cout <<"----------event " << count
                     << "-------[ " << timer.GetElapsedTimeStr() << " ]------"
                     << "\r" << flush;
            }

            auto param = AnalyzeEvent(sys, event, 3.0);

            hist_size.Fill(param.group_size);
            hist_msize.Fill(param.max_group_size);
            hist_mene.Fill(param.group_energy.maximum);
            hist_unif.Fill(param.max_group_energy.uniform);
            hist_rsq.Fill(param.max_group_line.rsq);
            hist_chisq.Fill(param.max_group_line.chisq);

            if((param.group_energy.maximum > 1.3 * beam_energy) ||
               ((param.max_group_size == 1) && (param.group_energy.maximum > 0.2 * beam_energy)) ||
               (param.max_group_size == 0))
                dst_parser2.WriteEvent(event);
            else
                dst_parser.WriteEvent(event);

        } else if (dst_parser.EventType() == PRadDSTParser::Type::epics) {
            dst_parser.WriteEPICS();
        }
    }

    cout <<"----------event " << count
         << "-------[ " << timer.GetElapsedTimeStr() << " ]------"
         << endl;

    dst_parser.CloseInput();
    dst_parser.CloseOutput();
    dst_parser2.CloseOutput();
    hist_size.Write();
    hist_msize.Write();
    hist_mene.Write();
    hist_unif.Write();
    hist_rsq.Write();
    hist_chisq.Write();
    f.Close();
}


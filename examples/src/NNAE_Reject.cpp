//============================================================================//
// An example showing how to read-in NNID root file and reject the event that //
// is categorized as cosmic ray (or background) by the Neural Network         //
//                                                                            //
// Chao Peng                                                                  //
// 07/05/2017                                                                 //
//============================================================================//


#include "PRadHyCalSystem.h"
#include "PRadDSTParser.h"
#include "PRadBenchMark.h"
#include "PRadCoordSystem.h"
#include "PRadInfoCenter.h"
#include "ConfigOption.h"
#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#define PROGRESS_COUNT 10000

using namespace std;
void RejectCosmic(const char *path, const char *nnid_file, const char *out_dir,
                  double ang_min, double ang_max, double ene_min, double ene_max);

int main(int argc, char *argv[])
{
    ConfigOption conf_opt;
    conf_opt.AddOpt(ConfigOption::help_message, 'h');
    conf_opt.AddLongOpt(ConfigOption::arg_require, "cut-angle-min", 'a');
    conf_opt.AddLongOpt(ConfigOption::arg_require, "cut-angle-max", 'b');
    conf_opt.AddLongOpt(ConfigOption::arg_require, "cut-energy-min", 'c');
    conf_opt.AddLongOpt(ConfigOption::arg_require, "cut-energy-max", 'd');

    conf_opt.SetDesc("usage: NNAE_Prepare <data_file> <root_file> <out_dir>");
    conf_opt.SetDesc('a', "cut minimum angle of the most energetic cluster.");
    conf_opt.SetDesc('b', "cut maximum angle of the most energetic cluster.");
    conf_opt.SetDesc('c', "cut minimum energy of the most energetic cluster.");
    conf_opt.SetDesc('d', "cut minimum energy of the most energetic cluster.");

    if(!conf_opt.ParseArgs(argc, argv) || conf_opt.NbofArgs() != 3) {
        std::cout << conf_opt.GetInstruction() << std::endl;
        return -1;
    }

    double ang_min = -1., ang_max = -1., ene_min = -1., ene_max = -1.;
    for(auto &opt : conf_opt.GetOptions())
    {
        switch(opt.mark)
        {
        case 'a': ang_min = opt.var.Double(); break;
        case 'b': ang_max = opt.var.Double(); break;
        case 'c': ene_min = opt.var.Double(); break;
        case 'd': ene_max = opt.var.Double(); break;
        default :
            std::cout << conf_opt.GetInstruction() << std::endl;
            return -1;
        }
    }

    RejectCosmic(conf_opt.GetArgument(0).c_str(),
                 conf_opt.GetArgument(1).c_str(),
                 conf_opt.GetArgument(2).c_str(),
                 ang_min, ang_max, ene_min, ene_max);

    return 0;
}

void RejectCosmic(const char *path, const char *nnid_file, const char *out_dir,
                  double ang_min, double ang_max, double ene_min, double ene_max)
{
    string prad_path = getenv("PRAD_PATH");
    if(prad_path.size()) {
        prad_path += "/";
    }

    // read-in root file and save the rejected event numbers
    std::unordered_set<int> rej_list;
    TFile *myfile= new TFile(nnid_file,"READ");
    TTree *tree= (TTree*) myfile->Get("Tid");
    double nn_id;
    int EvNb;
    tree->SetBranchAddress("nn_id",&nn_id);
    tree->SetBranchAddress("EvNb",&EvNb);
    int entries= tree->GetEntries();
    for(int i = 0; i < entries; i++)
    {
        tree->GetEntry(i);
        // rejection criterium nn_id < 0
        if(nn_id < 0.) {
            rej_list.insert(EvNb);
        }
    }

    // reject events
    PRadDSTParser dst_parser, dst_parser2;
    PRadHyCalSystem sys(prad_path + "config/hycal.conf");
    sys.ChooseRun(path);
    PRadCoordSystem coord(prad_path + "database/coordinates.dat");
    coord.ChooseCoord(PRadInfoCenter::GetRunNumber());

    string name = ConfigParser::decompose_path(path).name;
    string save_path = ConfigParser::form_path(out_dir, name + "_sav.dst");
    string rej_path = ConfigParser::form_path(out_dir, name + "_rej.dst");

    dst_parser.OpenInput(path);
    dst_parser.OpenOutput(save_path);
    dst_parser2.OpenOutput(rej_path);
    int count = 0;
    PRadBenchMark timer;

    while(dst_parser.Read())
    {
        if(dst_parser.EventType() == PRadDSTParser::Type::event) {

            if(++count % PROGRESS_COUNT == 0) {
                cout << "----------event " << count
                     << "-------[ " << timer.GetElapsedTimeStr() << " ]------"
                     << "\r" << flush;
            }

            auto event = dst_parser.GetEvent();
            sys.Reconstruct(event);
            auto &hits = sys.GetDetector()->GetHits();
            coord.TransformHits(sys.GetDetector());

            if(hits.empty()) continue;

            auto &emax_hit = *MostEnergeticHit(hits.begin(), hits.end());
            if(ene_min > 0. && emax_hit.E < ene_min) continue;
            if(ene_max > 0. && emax_hit.E > ene_max) continue;
            if(ang_min > 0. && PRadCoordSystem::GetPolarAngle(emax_hit) < ang_min) continue;
            if(ang_max > 0. && PRadCoordSystem::GetPolarAngle(emax_hit) > ang_max) continue;

            // in the rejected list
            if(rej_list.find(event.event_number) != rej_list.end()) {
                dst_parser2.Write(event);
            } else {
                dst_parser.Write(event);
            }

        } else if(dst_parser.EventType() == PRadDSTParser::Type::epics) {
            dst_parser.WriteEPICS();
        }
    }

    cout << "----------event " << count
         << "-------[ " << timer.GetElapsedTimeStr() << " ]------"
         << endl;

    dst_parser.CloseInput();
    dst_parser.CloseOutput();
    dst_parser2.CloseOutput();
}

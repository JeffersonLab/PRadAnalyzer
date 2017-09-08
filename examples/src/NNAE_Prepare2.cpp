//============================================================================//
// An example showing how to extract HyCal information of each event in root  //
// file. The event is pre-analyzed and some features are extracted.           //
// Output can be used for external Neural Network application                 //
//                                                                            //
// Chao Peng                                                                  //
// 07/05/2017                                                                 //
//============================================================================//


#include "PRadHyCalSystem.h"
#include "PRadGEMSystem.h"
#include "PRadDetMatch.h"
#include "PRadDSTParser.h"
#include "PRadBenchMark.h"
#include "PRadCoordSystem.h"
#include "PRadInfoCenter.h"
#include "cosmicEval.h"
#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#define PROGRESS_COUNT 10000

using namespace std;
void SavePixels(const char *path, const char *save_path);

int main(int argc, char *argv[])
{

    if(argc != 3) {
        cerr << "usage: NNAE_Prepare <data_file> <root_file>" << endl;
        return -1;
    }

    SavePixels(argv[1], argv[2]);
    return 0;
}


void SavePixels(const char *path, const char *save_path)
{
    string prad_path = getenv("PRAD_PATH");
    if(prad_path.size()) {
        prad_path += "/";
    }

    PRadDSTParser dst_parser;
    PRadHyCalSystem sys(prad_path + "config/hycal.conf");
    sys.ChooseRun(path);
    PRadCoordSystem coord(prad_path + "database/coordinates.dat");
    coord.ChooseCoord(PRadInfoCenter::GetRunNumber());
    PRadGEMSystem gem(prad_path + "config/gem.conf");
    PRadDetMatch det_match(prad_path + "config/det_match.conf");

    string name = ConfigParser::decompose_path(path).name + "_params";

    int32_t EvNb;
    bool discharge;
    EventParam param;

    TFile file(save_path,"RECREATE");
    TTree T("T","T");
    T.Branch("EvNb", &EvNb, "EvNb/I");
    T.Branch("discharge", &discharge, "discharge/O");
    T.Branch("Nhits", &param.group_size, "Nhits/i");
    T.Branch("MCl_Nhits", &param.max_group_size, "MCl_Nhits/i");
    T.Branch("energy_param", &param.group_energy, "total/D:maximum:uniform");
    T.Branch("MCl_energy_param", &param.max_group_energy, "total/D:maximum:uniform");
    T.Branch("line_param", &param.group_line, "success/D:k:b:rsq:chisq");
    T.Branch("MCl_line_param", &param.max_group_line, "success/D:k:b:rsq:chisq");
    T.Branch("Ngroups", &param.n_groups, "Ngroups/i");
    T.Branch("group_hits", &param.group_hits, "group_hits[Ngroups]/i");

    dst_parser.OpenInput(path);
    int count = 0;
    PRadBenchMark timer;

    while(dst_parser.Read())
    {
        if(dst_parser.EventType() == PRadDSTParser::Type::event) {

            auto event = dst_parser.GetEvent();
            if(!event.is_physics_event())
                continue;

            if((++count)%PROGRESS_COUNT == 0) {
                cout << "----------event " << count
                     << "-------[ " << timer.GetElapsedTimeStr() << " ]------"
                     << "\r" << flush;
            }
            EvNb = event.event_number;
            param = AnalyzeEvent(&sys, event);

            double max_e = 0., next_e = 0.;
            for(auto &module : sys.GetDetector()->GetModuleList())
            {
                double e = module->GetEnergy();
                if(e > max_e) {
                    next_e = max_e;
                    max_e = e;
                }
            }

            if(max_e > 300. && next_e < 5.)
                discharge = true;
            else
                discharge = false;

            T.Fill();
            /*
            std::cout << param.group_size << ", "
                      << param.max_group_size << ", "
                      << param.group_energy.total << ", "
                      << param.group_energy.maximum << ", "
                      << param.group_energy.uniform << ", "
                      << param.group_line.k << ", "
                      << param.group_line.b << ", "
                      << param.group_line.rsq << ", "
                      << param.group_line.chisq << ", "
                      << std::endl;
            T.Show(T.GetEntries()-1);
            */
        }
    }

    cout << "----------event " << count
         << "-------[ " << timer.GetElapsedTimeStr() << " ]------"
         << endl;

    dst_parser.CloseInput();
    T.Write();
    file.Save();
}

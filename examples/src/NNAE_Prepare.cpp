//============================================================================//
// An example showing how to extract HyCal information of each event in root  //
// file. The output can be used for external Neural Network application       //
//                                                                            //
// Chao Peng                                                                  //
// 07/05/2017                                                                 //
//============================================================================//


#include "PRadHyCalSystem.h"
#include "PRadDSTParser.h"
#include "PRadBenchMark.h"
#include "PRadCoordSystem.h"
#include "PRadInfoCenter.h"
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

    string name = ConfigParser::decompose_path(path).name;

    unsigned int Nhits;
    int EvNb;
    double theta;
    double id[1728], E[1728];

    TFile file(save_path,"RECREATE");
    TTree T("T","T");
    T.Branch("Nhits", &Nhits, "Nhits/i");
    T.Branch("EvNb", &EvNb, "EvNb/I");
    T.Branch("Theta", &theta, "Theta/D");
    T.Branch("ID", &id[0], "ID[Nhits]/D");
    T.Branch("Energy", &E[0], "Energy[Nhits]/D");

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

            // update event information to all HyCal modules
            sys.ChooseEvent(event);
            sys.GetDetector()->CollectHits();

            Nhits = sys.GetDetector()->GetModuleHits().size();
            for(unsigned int i = 0; i < Nhits; ++i)
            {
                id[i] = sys.GetDetector()->GetModuleHits().at(i).id;
                E[i] = sys.GetDetector()->GetModuleHits().at(i).energy;
            }

            // reconstruct event and get the most energized cluster's theta angle
            sys.Reconstruct(event);
            auto &hits = sys.GetDetector()->GetHits();
            coord.Transform(PRadDetector::HyCal, hits.begin(), hits.end());
            if(hits.empty()) {
                theta = 0.;
            } else {
                double max_energy = hits.front().E;
                theta = PRadCoordSystem::GetPolarAngle(hits.front());
                for(auto &hit : hits)
                {
                    if(hit.E > max_energy) {
                        max_energy = hit.E;
                        theta = PRadCoordSystem::GetPolarAngle(hit);
                    }
                }
            }

            T.Fill();
        }
    }

    cout << "----------event " << count
         << "-------[ " << timer.GetElapsedTimeStr() << " ]------"
         << endl;

    dst_parser.CloseInput();
    T.Write();
    file.Save();
}

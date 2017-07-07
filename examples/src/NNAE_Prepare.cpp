//============================================================================//
// An example showing how to extract HyCal information of each event in root  //
// file. The output can be used for external Neural Network application       //
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

    string name = ConfigParser::decompose_path(path).name;

    unsigned int Nclusters, Nhits;
    int EvNb;
    double theta, good;
    double id[1728], E[1728];
    double ClAngle[1000], ClEnergy[1000];

    TFile file(save_path,"RECREATE");
    TTree T("T","T");
    T.Branch("Nclusters", &Nclusters, "Nclusters/i");
    T.Branch("ClAngle", &ClAngle, "ClAngle[Nclusters]/D");
    T.Branch("ClEnergy", &ClEnergy, "ClEnergy[Nclusters]/D");
    T.Branch("Nhits", &Nhits, "Nhits/i");
    T.Branch("EvNb", &EvNb, "EvNb/I");
    T.Branch("Theta", &theta, "Theta/D");
    T.Branch("Good", &good, "Good/D");
    T.Branch("ID", &id[0], "ID[Nhits]/D");
    T.Branch("Energy", &E[0], "Energy[Nhits]/D");

    dst_parser.OpenInput(path);
    int count = 0;
    PRadBenchMark timer;

    PRadGEMDetector *gem1 = gem.GetDetector("PRadGEM1");
    PRadGEMDetector *gem2 = gem.GetDetector("PRadGEM2");

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
            gem.Reconstruct(event);
            auto &hits = sys.GetDetector()->GetHits();

            coord.TransformHits(sys.GetDetector());
            coord.TransformHits(gem1);
            coord.TransformHits(gem2);

            Nclusters = hits.size();
            for(unsigned int i = 0; i < Nclusters; ++i)
            {
                ClEnergy[i] = hits[i].E;
                ClAngle[i] = PRadCoordSystem::GetPolarAngle(hits[i]);
            }

            // hits matching, return matched index
            auto matched = det_match.Match(hits, gem1->GetHits(), gem2->GetHits());

            if(!matched.empty())
                good = 1.;
            else
                good = -1.;

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

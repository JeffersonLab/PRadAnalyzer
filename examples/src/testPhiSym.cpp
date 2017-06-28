//============================================================================//
// An example showing how to use HyCAL clustering method to reconstruct data  //
//                                                                            //
// Chao Peng                                                                  //
// 11/12/2016                                                                 //
//============================================================================//

#include "PRadHyCalSystem.h"
#include "PRadDSTParser.h"
#include "PRadBenchMark.h"
#include "PRadInfoCenter.h"
#include "PRadCoordSystem.h"
#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#define PROGRESS_COUNT 10000
using namespace std;

int main(int argc, char *argv[])
{
    if(argc != 2)
        return 0;

    string file = argv[1];
    double beam_energy = 2142;

    PRadHyCalSystem *sys = new PRadHyCalSystem("config/hycal.conf");
    sys->ChooseRun(file);
    PRadDSTParser *dst_parser = new PRadDSTParser();
    PRadCoordSystem *coord = new PRadCoordSystem("database/coordinates.dat");
    coord->ChooseCoord(PRadInfoCenter::GetRunNumber());

    dst_parser->OpenInput(file);
    TFile *f = new TFile((ConfigParser::decompose_path(file).name + "_phisym.root").c_str(), "RECREATE");
    TH1F *hist1 = new TH1F("phi_pwo", "phi_sym", 18, -180, 180);
    TH1F *hist2 = new TH1F("phi_lg", "phi_sym", 18, -180, 180);
    TH1F *hist3 = new TH1F("phi_trans", "phi_sym", 18, -180, 180);

    PRadBenchMark timer;
    int count = 0;

    while(dst_parser->Read())
    {
        if(dst_parser->EventType() == PRadDSTParser::Type::event) {

            auto event = dst_parser->GetEvent();
            if(!event.is_physics_event())
                continue;

            if(++count%PROGRESS_COUNT == 0) {
                cout <<"------[ ev " << count << " ]---"
                     << "---[ " << timer.GetElapsedTimeStr() << " ]---"
                     << "---[ " << timer.GetElapsedTime()/(double)count << " ms/ev ]------"
                     << "\r" << flush;
            }

            sys->Reconstruct(event);
            auto& hits = sys->GetDetector()->GetHits();
            coord->Transform(PRadDetector::HyCal, hits.begin(), hits.end());
            for(auto &hit : hits)
            {
                if(std::abs(1. - hit.E/beam_energy) >= 0.2 ||
                   PRadCoordSystem::GetPolarAngle(hit) < 0.7)
                    continue;
                if(TEST_BIT(hit.flag, kTransition))
                    hist3->Fill(PRadCoordSystem::GetAzimuthalAngle(hit));
                else if(TEST_BIT(hit.flag, kPbWO4))
                    hist1->Fill(PRadCoordSystem::GetAzimuthalAngle(hit));
                else if(TEST_BIT(hit.flag, kPbGlass))
                    hist2->Fill(PRadCoordSystem::GetAzimuthalAngle(hit));
            }
        }
    }

    hist1->Write();
    hist2->Write();
    hist3->Write();
    f->Close();
    dst_parser->CloseInput();
    cout <<"------[ ev " << count << " ]---"
         << "---[ " << timer.GetElapsedTimeStr() << " ]---"
         << "---[ " << timer.GetElapsedTime()/(double)count << " ms/ev ]------"
         << endl;
    cout << "Finished." << endl;

    return 0;
}

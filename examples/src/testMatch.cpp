//============================================================================//
// An example showing how to reconstruct, transform and match the detectors   //
//                                                                            //
// Chao Peng                                                                  //
// 10/24/2016                                                                 //
//============================================================================//

#include "PRadDataHandler.h"
#include "PRadDSTParser.h"
#include "PRadInfoCenter.h"
#include "PRadBenchMark.h"
#include "PRadEPICSystem.h"
#include "PRadTaggerSystem.h"
#include "PRadHyCalSystem.h"
#include "PRadGEMSystem.h"
#include "PRadCoordSystem.h"
#include "PRadDetMatch.h"
#include "canalib.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include <iostream>
#include <string>
#include <vector>

#define PROGRESS_COUNT 10000

using namespace std;

PRadEPICSystem *epics;
PRadHyCalSystem *hycal;
PRadGEMSystem *gem;
PRadCoordSystem *coord_sys;
PRadDetMatch *det_match;
void testMatch(const string &path);
ostream &operator <<(ostream &os, const PRadBenchMark &timer);

int main(int argc, char *argv[])
{
    if(argc < 2) {
        cout << "usage: testMatch <file1> <file2> ..." << endl;
        return 0;
    }

    // initialize objects
    epics = new PRadEPICSystem("config/epics_channels.conf");
    gem = new PRadGEMSystem("config/gem.conf");
    hycal = new PRadHyCalSystem("config/hycal.conf");
    coord_sys = new PRadCoordSystem("database/coordinates.dat");
    det_match = new PRadDetMatch("config/det_match.conf");

    for(int i = 1; i < argc; ++i)
    {
        string file = argv[i];
        testMatch(file);
    }

    return 0;
}

void testMatch(const string &path)
{
    hycal->ChooseRun(path);
    coord_sys->ChooseCoord(PRadInfoCenter::GetRunNumber());

    PRadDSTParser *dst_parser = new PRadDSTParser();
    dst_parser->OpenInput(path);
    dst_parser->OpenOutput(ConfigParser::decompose_path(path).name + "_match.dst");

    string outfile = ConfigParser::decompose_path(path).name + "_match.root";
    TFile f(outfile.c_str(), "RECREATE");
    TH1F *hist[3];
    hist[0] = new TH1F("PbGlass R Diff", "Diff in R", 1000, -100, 100);
    hist[1] = new TH1F("PbWO4 R Diff", "Diff in R", 1000, -100, 100);
    hist[2] = new TH1F("Trans R Diff", "Diff in R", 1000, -100, 100);
    TH2F *hist2d = new TH2F("R Diff", "HyCal - GEM", 800, 0, 800, 200, -100, 100);
    TH2F *histev = new TH2F("E vs theta", "Event Distribution", 200, 0, 8, 200, 0, 1400);

    PRadHyCalDetector *hycal_det = hycal->GetDetector();
    PRadGEMDetector *gem_det1 = gem->GetDetector("PRadGEM1");
    PRadGEMDetector *gem_det2 = gem->GetDetector("PRadGEM2");

    PRadBenchMark timer;
    int count = 0;

    while(dst_parser->Read())
    {
        if(dst_parser->EventType() == PRadDSTParser::Type::event) {
            auto &event = dst_parser->GetEvent();

            // only interested in physics event
            if(!event.is_physics_event())
                continue;

            if((++count)%PROGRESS_COUNT == 0) {
                cout <<"------[ ev " << count << " ]---"
                     << "---[ " << timer.GetElapsedTimeStr() << " ]---"
                     << "---[ " << timer.GetElapsedTime()/(double)count << " ms/ev ]------"
                     << "\r" << flush;
            }

            // update run information
            PRadInfoCenter::Instance().UpdateInfo(event);

            // reconstruct
            hycal->ChooseEvent(event);
            hycal->Reconstruct();
            gem->Reconstruct(event);

            // get reconstructed clusters
            auto &hycal_hit = hycal_det->GetHits();
            auto &gem1_hit = gem_det1->GetHits();
            auto &gem2_hit = gem_det2->GetHits();

            // coordinates transform, projection
            coord_sys->Transform(PRadDetector::HyCal, hycal_hit.begin(), hycal_hit.end());
            coord_sys->Transform(PRadDetector::PRadGEM1, gem1_hit.begin(), gem1_hit.end());
            coord_sys->Transform(PRadDetector::PRadGEM2, gem2_hit.begin(), gem2_hit.end());

            // hits matching, return matched index
            auto matched = det_match->Match(hycal_hit, gem1_hit, gem2_hit);

            // project to HyCal surface
            coord_sys->Projection(matched.begin(), matched.end());

            for(auto &hit : matched)
            {
                int hidx = 0;
                if(TEST_BIT(hit.hycal.flag, kPbWO4))
                    hidx = 1;
                if(TEST_BIT(hit.hycal.flag, kTransition))
                    hidx = 2;

                float r = sqrt(hit.hycal.x*hit.hycal.x + hit.hycal.y*hit.hycal.y);
                float r2 = sqrt(hit.x*hit.x + hit.y*hit.y);
                hist[hidx]->Fill(r - r2);
                hist2d->Fill(r, r - r2);

                float angle = coord_sys->GetPolarAngle(hit);
                histev->Fill(angle, hit.E);
            }

            if(matched.size() >= 1 && matched.size() <= 2)
                dst_parser->WriteEvent(event);

        } else if(dst_parser->EventType() == PRadDSTParser::Type::epics) {
            // save epics into handler, otherwise get epicsvalue won't work
            epics->AddEvent(dst_parser->GetEPICSEvent());
            dst_parser->WriteEPICS();
        }
    }

    dst_parser->CloseInput();
    dst_parser->CloseOutput();

    cout <<"------[ ev " << count << " ]---"
         << "---[ " << timer.GetElapsedTimeStr() << " ]---"
         << "---[ " << timer.GetElapsedTime()/(double)count << " ms/ev ]------"
         << endl;
    cout << "TIMER: Finished, took " << timer.GetElapsedTime() << " ms" << endl;
    cout << PRadInfoCenter::GetBeamCharge() << endl;
    cout << PRadInfoCenter::GetLiveTime() << endl;

    for(int i = 0; i < 3; ++i)
    {
        hist[i]->Write();
    }
    hist2d->Write();
    histev->Write();
    f.Close();
}

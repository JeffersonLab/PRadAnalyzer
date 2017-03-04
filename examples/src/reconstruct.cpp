//============================================================================//
// An example showing how to reconstruct, transform and match the detectors   //
//                                                                            //
// Chao Peng                                                                  //
// 10/24/2016                                                                 //
//============================================================================//

#include "PRadDSTParser.h"
#include "PRadInfoCenter.h"
#include "PRadBenchMark.h"
#include "PRadEPICSystem.h"
#include "PRadHyCalSystem.h"
#include "PRadGEMSystem.h"
#include "PRadCoordSystem.h"
#include "PRadDetMatch.h"
#include "canalib.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include <iostream>
#include <string>
#include <vector>

#define PROGRESS_COUNT 10000

using namespace std;

void reconstruct(const char *path);

int main(int argc, char *argv[])
{
    if(argc != 2) {
        cout << "usage: reconstruct <file>" << endl;
        return 0;
    }

    reconstruct(argv[1]);

    return 0;
}

void reconstruct(const char *path)
{
    PRadEPICSystem epics("config/epics_channels.conf");
    PRadGEMSystem gem("config/gem.conf");
    PRadHyCalSystem hycal("config/hycal.conf");
    PRadCoordSystem coord_sys("database/coordinates.dat");
    PRadDetMatch det_match("config/det_match.conf");

    PRadHyCalDetector *hycal_det = hycal.GetDetector();
    PRadGEMDetector *gem_det1 = gem.GetDetector("PRadGEM1");
    PRadGEMDetector *gem_det2 = gem.GetDetector("PRadGEM2");

    // get run number from the file name and update the calibration
    // constants and run info accordingly
    PRadInfoCenter::SetRunNumber(path);
    hycal.UpdateRunFiles();
    coord_sys.ChooseCoord(PRadInfoCenter::GetRunNumber());

    PRadDSTParser dst_parser;
    dst_parser.OpenInput(path);

    // root file
    string rootfile = ConfigParser::decompose_path(path).name + "_recon.root";
    TFile f(rootfile.c_str(), "RECREATE");
    TTree t("T", "T");

    const size_t max_hits = 100;
    int N, Nm;
    double E[max_hits], x[max_hits], y[max_hits], z[max_hits];
    bool match[max_hits];
    t.Branch("Nhits", &N, "Nhits/I");
    t.Branch("Nmatch", &N, "Nmatch/I");
    t.Branch("E", &E, "E/D");
    t.Branch("X", &x[0], "X[Nhits]/D");
    t.Branch("Y", &y[0], "Y[Nhits]/D");
    t.Branch("Z", &z[0], "Z[Nhits]/D");
    t.Branch("Match", &match[0], "Match[Nhits]/O");

    PRadBenchMark timer;
    int count = 0;
    while(dst_parser.Read())
    {
        if(dst_parser.EventType() == PRadDSTParser::Type::event) {

            auto &event = dst_parser.GetEvent();
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
            hycal.Reconstruct(event);
            gem.Reconstruct(event);

            // get reconstructed clusters
            auto &hycal_hit = hycal_det->GetHits();
            auto &gem1_hit = gem_det1->GetHits();
            auto &gem2_hit = gem_det2->GetHits();

            // fill in arrays for writing to root file
            N = (hycal_hit.size() > max_hits)? max_hits : hycal_hit.size();
            for(int i = 0; i < N; ++i)
            {
                E[i] = hycal_hit[i].E;
                x[i] = hycal_hit[i].x;
                y[i] = hycal_hit[i].y;
                z[i] = hycal_hit[i].z;
                match[i] = false;
            }

            // coordinates transform, projection
            coord_sys.Transform(PRadDetector::HyCal, hycal_hit.begin(), hycal_hit.end());
            coord_sys.Transform(PRadDetector::PRadGEM1, gem1_hit.begin(), gem1_hit.end());
            coord_sys.Transform(PRadDetector::PRadGEM2, gem2_hit.begin(), gem2_hit.end());

            coord_sys.Projection(hycal_hit.begin(), hycal_hit.end());
            coord_sys.Projection(gem1_hit.begin(), gem1_hit.end());
            coord_sys.Projection(gem2_hit.begin(), gem2_hit.end());

            // hits matching, return matched index
            auto matched = det_match.Match(hycal_hit, gem1_hit, gem2_hit);

            // replace the hycal hits with the matched information
            Nm = 0;
            for(auto &m : matched)
            {
                if((int)m.hycal_idx >= N)
                    continue;

                x[m.hycal_idx] = m.x;
                y[m.hycal_idx] = m.y;
                z[m.hycal_idx] = m.z;
                match[m.hycal_idx] = true;
                Nm ++;
            }

            t.Fill();

        } else if(dst_parser.EventType() == PRadDSTParser::Type::epics) {
            // save epics into handler, otherwise get epicsvalue won't work
            epics.AddEvent(dst_parser.GetEPICSEvent());
        }
    }

    dst_parser.CloseInput();

    cout <<"------[ ev " << count << " ]---"
         << "---[ " << timer.GetElapsedTimeStr() << " ]---"
         << "---[ " << timer.GetElapsedTime()/(double)count << " ms/ev ]------"
         << endl;
    cout << "TIMER: Finished, took " << timer.GetElapsedTime() << " ms" << endl;
    cout << PRadInfoCenter::GetLiveBeamCharge() << endl;
    cout << PRadInfoCenter::GetLiveTime() << endl;

    t.Write();
    f.Close();
}

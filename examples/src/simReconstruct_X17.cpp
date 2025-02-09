//============================================================================//
// An example showing how to reconstruct, transform and match the detectors   //
// For simulation                                                             //
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
#include "TRandom2.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include <iostream>
#include <string>
#include <vector>

#define PROGRESS_COUNT 10000
#define NMODULES 1728

static TRandom2 *RandGen = new TRandom2();
static double fResoPar[3000][3];

// we would get the exact position from GUN.Z, but to include the target profile effects, here it is a fixed position
static const double target_center = -3000. + 89.0;
static const double hycal_pwo_surf = 2735.15;
static const double hycal_lg_diff = - 97.3;

void simReconstruct(TChain *t, const char *outf);
void setModuleEnergy(PRadHyCalModule &m, double edep);

int main(int argc, char *argv[])
{
    if(argc < 3) {
        std::cout << "usage: reconstruct <file1> [<file2> <file3> ...] <out_file>" << std::endl;
        return -1;
    }

    // read simulation files
    auto t = new TChain("T");
    for (int i = 1; i < argc - 1; ++i) {
        t->Add(argv[i]);
    }

    // read resolution parameters
    ConfigParser parser;
    if (!parser.ReadFile("./database/hycal_resolution_curve_2terms.dat")) {
        std::cout << "cannot find hycal_resolution_curve.dat" << std::endl;
    //    exit(0);
    }

    double input[3];
    while (parser.ParseLine()) {
        std::string ch_name = parser.TakeFirst();
        for (int i = 0; i < 3; i++) { input[i] = parser.TakeFirst().Double(); }
        int id = PRadHyCalModule::name_to_id(ch_name);
        if (id <= 0) { continue; }
        for (int i = 0; i < 3; i++) { fResoPar[id - 1][i] = input[i]; }
    }

    simReconstruct(t, argv[argc - 1]);
    return 0;
}

void simReconstruct(TChain *tch, const char *outf)
{
    TFile *f = new TFile(outf, "RECREATE");

    // just copy these branches
    std::vector<std::string> selections = {
        "GUN.N", "GUN.PID", "GUN.X", "GUN.Y", "GUN.Z", "GUN.E", "GUN.Theta", "GUN.Phi",
//        "GEM.N", "GEM.X", "GEM.Y", "GEM.Z", "GEM.DID",
//        "ST.N", "ST.X", "ST.Y", "ST.Z",
    };

    tch->SetBranchStatus("*", 0);
    for (auto &br : selections) {
        tch->SetBranchStatus(br.c_str(), 1);
    }
    auto t = tch->CloneTree(0);

    // hits reconstruction
    int Nhits, CID[100], match[100];
    double hE[100], hX[100], hY[100], hZ[100];
    // retrieve part of the cluster information
    t->Branch("Hit.N", &Nhits, "Hit.N/I");
    t->Branch("Hit.X", hX, "Hit.X[Hit.N]/D");
    t->Branch("Hit.Y", hY, "Hit.Y[Hit.N]/D");
    t->Branch("Hit.Z", hZ, "Hit.Z[Hit.N]/D");
    t->Branch("Hit.E", hE, "Hit.E[Hit.N]/D");
    t->Branch("Hit.CID", CID, "Hit.CID[Hit.N]/I");
    t->Branch("Hit.Match", match, "Hit.Match[Hit.N]/I");

    // setup hycal
    PRadHyCalSystem hycal("config/hycal.conf");
    PRadHyCalDetector *hycal_det = hycal.GetDetector();
    hycal_det->SortModuleList();
    auto mlist = hycal_det->GetModuleList();

    // set branch addresses
    tch->SetBranchStatus("HC.ModuleEdep", 1);
    tch->SetBranchStatus("GEM.*", 1);
    int Ngem, gdid[100];
    double medep[NMODULES], gx[100], gy[100], gz[100];
    tch->SetBranchAddress("HC.ModuleEdep", medep);
    tch->SetBranchAddress("GEM.N", &Ngem);
    tch->SetBranchAddress("GEM.X", gx);
    tch->SetBranchAddress("GEM.Y", gy);
    tch->SetBranchAddress("GEM.Z", gz);
    tch->SetBranchAddress("GEM.DID", gdid);

    PRadBenchMark timer;
    Point target_point = Point(0., 0., target_center);

    std::vector<std::string> first_layer = {"W526", "W527", "W528", "W529",
                                            "W560", "W563",
                                            "W594", "W597",
                                            "W628", "W629", "W630", "W631"};
    std::vector<std::string> second_layer = {"W491", "W492", "W493", "W494", "W495", "W496",
                                            "W525", "W530",
                                            "W559", "W564",
                                            "W593", "W598",
                                            "W627", "W632",
                                            "W661", "W662", "W663", "W664", "W665", "W666"};
    std::vector<PRadHyCalModule*> turnedoffs;
    for (auto &name : first_layer) { turnedoffs.push_back(hycal.GetModule(name)); }
    for (auto &name : second_layer) { turnedoffs.push_back(hycal.GetModule(name)); }
    // loop over all events
    for (int i = 0; i < tch->GetEntries(); ++i) {
        if ((i + 1) % PROGRESS_COUNT == 0) {
            std::cout <<"------[ ev " << i + 1 << " ]---"
                      << "---[ " << timer.GetElapsedTimeStr() << " ]---"
                      << "---[ " << timer.GetElapsedTime()/(double)(i + 1) << " ms/ev ]------"
                      << "\r" << std::flush;
        }
        tch->GetEntry(i);
        // Update HyCal energy
        for (size_t k = 0; k < mlist.size(); ++k) {
            setModuleEnergy(*mlist[k], medep[k]);
        }
        for (auto turnedoff : turnedoffs) {
            setModuleEnergy(*turnedoff, 0.);
        }

        hycal.Reconstruct();
        // match hits
        auto &hyhits = hycal_det->GetHits();
        Nhits = hyhits.size();
        for (int j = 0; j < Nhits; ++j) {
            auto &hit = hyhits[j];
            CID[j] = hit.cid;
            hE[j] = hit.E;
            hX[j] = hit.x;
            hY[j] = hit.y;
            hZ[j] = hit.z;
            match[j] = 0;
/*
            match[j] = 0;
            Point hpos(hit.x, hit.y, hit.z + hycal_pwo_surf);
            double best_dist = 1000.;
            Point best_match = hpos;

            // simulation output has an inverted x-axis as compared to the analyzer
            for (int k = 0; k < Ngem; ++k) {
                Point gpos(-gx[k], gy[k], gz[k]);
                auto dist = PRadCoordSystem::ProjectionDistance(gpos, hpos, target_point, hycal_pwo_surf);
                if (dist < 5.0*2.6/std::sqrt(hit.E/1000.)) {
                    SET_BIT(match[j], gdid[k]);
                    if (dist < best_dist) {
                        best_dist = dist;
                        best_match = gpos;
                    }
                }
            }
            PRadCoordSystem::Projection(best_match, target_point, hycal_pwo_surf);
            hX[j] = best_match.x;
            hY[j] = best_match.y;
            hZ[j] = best_match.z;
*/
        }
        t->Fill();
    }

    std::cout <<"------[ ev " << tch->GetEntries() << " ]---"
              << "---[ " << timer.GetElapsedTimeStr() << " ]---"
              << "---[ " << timer.GetElapsedTime()/(double)tch->GetEntries() << " ms/ev ]------"
              << std::endl;
    t->Write();
    f->Close();
}

void setModuleEnergy(PRadHyCalModule &m, double edep)
{
    double ped = RandGen->Gaus(m.GetChannel()->GetPedestal().mean, m.GetChannel()->GetPedestal().sigma);
    auto ch = m.GetChannel();
    if (ch->IsDead() || edep > 2e4 || edep < 1e-3) {
        ch->SetValue(ped);
        return;
    }
    double reso = std::sqrt(0.76) * std::sqrt(pow(fResoPar[m.GetID() - 1][0] / std::sqrt(edep / 1000.), 2) +
                  pow(fResoPar[m.GetID() - 1][1] / (edep / 1000.), 2) +
                  pow(fResoPar[m.GetID() - 1][2], 2));
    double val = ped + (RandGen->Gaus(edep,  edep * reso)) / m.GetCalibrationFactor();
    ch->SetValue(val);
}


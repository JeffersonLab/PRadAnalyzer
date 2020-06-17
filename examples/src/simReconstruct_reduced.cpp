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
#include "TRandom3.h"
#include "TH1.h"
#include "TH2.h"
#include <iostream>
#include <string>
#include <vector>

#define PROGRESS_COUNT 10000
#define NMODULES 1728

static TRandom2 *RandGen = new TRandom2();
static double fResoPar[3000][3];

static const double target_center = -3000. + 88.9;
static const double hycal_pwo_surf = 2950.;
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
        "ST.N", "ST.X", "ST.Y", "ST.Z",
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
    t->Branch("Hits.N", &Nhits, "Hits.N/I");
    t->Branch("Hits.X", hX, "Hits.X[Hits.N]/D");
    t->Branch("Hits.Y", hY, "Hits.Y[Hits.N]/D");
    t->Branch("Hits.Z", hZ, "Hits.Z[Hits.N]/D");
    t->Branch("Hits.E", hE, "Hits.E[Hits.N]/D");
    t->Branch("Hits.CID", CID, "Hits.CID[Hits.N]/I");
    t->Branch("Hits.match", match, "Hits.match[Hits.N]/I");

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
    // tch->SetBranchAddress("HC.ModuleEdep", medep);
    tch->SetBranchAddress("GEM.N", &Ngem);
    tch->SetBranchAddress("GEM.X", gx);
    tch->SetBranchAddress("GEM.Y", gy);
    tch->SetBranchAddress("GEM.Z", gz);
    tch->SetBranchAddress("GEM.DID", gdid);

    int Nhyc;
    double cx[300], cy[300], cz[300], cE[300];
    tch->SetBranchAddress("HC.N", &Nhyc);
    tch->SetBranchAddress("HC.P", cE);
    tch->SetBranchAddress("HC.X", cx);
    tch->SetBranchAddress("HC.Y", cy);
    tch->SetBranchAddress("HC.Z", cz);


    PRadBenchMark timer;
    TRandom2 rng;
    Point target_point = Point(0., 0., target_center);
    // loop over all events
    for (int i = 0; i < tch->GetEntries(); ++i) {
        if ((i + 1) % PROGRESS_COUNT == 0) {
            std::cout <<"------[ ev " << i + 1 << " ]---"
                      << "---[ " << timer.GetElapsedTimeStr() << " ]---"
                      << "---[ " << timer.GetElapsedTime()/(double)(i + 1) << " ms/ev ]------"
                      << "\r" << std::flush;
        }
        tch->GetEntry(i);

        // match hits
        Nhits = 0;
        for (int j = 0; j < Nhyc; ++j, ++Nhits) {
            // probably a secondary hit
            // if ( cz[j] > hycal_pwo_surf + 1.0) { continue; }

            auto module = hycal_det->GetModule(cx[j], cy[j]);
            // use id to check resolution
            int cid = 1000;
            if (module) { cid = module->GetID(); }

            // smear resolution
            double res = (cid < 1000.) ? 0.065 : 0.026/std::sqrt(cE[j]/1000.);
            BaseHit hit(cx[j] + rng.Gaus(0., res*100.),
                        cy[j] + rng.Gaus(0., res*100.),
                        cz[j],
                        rng.Gaus(1.0, res)*cE[j]);
            CID[j] = cid;

            hE[j] = hit.E;
            Point hpos(hit.x, hit.y, hit.z + hycal_pwo_surf);
            match[j] = 0;
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
        }
        t->Fill();
    }

    std::cout <<"------[ ev " << tch->GetEntries() << " ]---"
              << "---[ " << timer.GetElapsedTimeStr() << " ]---"
              << "---[ " << timer.GetElapsedTime()/(double)tch->GetEntries() << " ms/ev ]------"
              << std::endl;
    t->Write();
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


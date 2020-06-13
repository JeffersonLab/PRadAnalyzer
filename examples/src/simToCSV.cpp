//============================================================================//
// An example showing how to convert root files to a csv file                 //
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

void toCSV(TChain *t, const char *outf);

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

    toCSV(t, argv[argc - 1]);
    return 0;
}

void toCSV(TChain *tch, const char *outf)
{
    auto path = ConfigParser::decompose_path(outf);
    path.name += "_gun";

    std::ofstream output(outf);
    std::ofstream outgun(ConfigParser::compose_path(path));

    // header
    char sep = ',';
    output << "Event Number" << sep << "x" << sep << "y" << sep << "z" << sep << "E" << sep
           << "match" << sep << "scint" << sep << "cid" << std::endl;
    outgun << "Event Number" << sep << "x" << sep << "y" << sep << "z" << sep << "E" << sep
           << "pid" << sep << "theta" << sep << "phi" << std::endl;

    // set branch addresses
    int Nhits, Nsc, match[100], cid[100];
    double hx[100], hy[100], hz[100], hE[100], sx[100], sy[100];
    tch->SetBranchAddress("Hits.N", &Nhits);
    tch->SetBranchAddress("Hits.X", hx);
    tch->SetBranchAddress("Hits.Y", hy);
    tch->SetBranchAddress("Hits.Z", hz);
    tch->SetBranchAddress("Hits.E", hE);
    tch->SetBranchAddress("Hits.CID", cid);
    tch->SetBranchAddress("Hits.match", match);
    tch->SetBranchAddress("ST.N", &Nsc);
    tch->SetBranchAddress("ST.X", sx);
    tch->SetBranchAddress("ST.Y", sy);

    int Ngun, gpid[100];
    double gx[100], gy[100], gz[100], gE[100], gth[100], gph[100];
    tch->SetBranchAddress("GUN.N", &Ngun);
    tch->SetBranchAddress("GUN.X", gx);
    tch->SetBranchAddress("GUN.Y", gy);
    tch->SetBranchAddress("GUN.Z", gz);
    tch->SetBranchAddress("GUN.E", gE);
    tch->SetBranchAddress("GUN.PID", gpid);
    tch->SetBranchAddress("GUN.Theta", gth);
    tch->SetBranchAddress("GUN.Phi", gph);



    PRadBenchMark timer;
    // loop over all events
    for (int iev = 0; iev < tch->GetEntries(); ++iev) {
        if ((iev + 1) % PROGRESS_COUNT == 0) {
            std::cout <<"------[ ev " << iev + 1 << " ]---"
                      << "---[ " << timer.GetElapsedTimeStr() << " ]---"
                      << "---[ " << timer.GetElapsedTime()/(double)(iev + 1) << " ms/ev ]------"
                      << "\r" << std::flush;
        }
        tch->GetEntry(iev);

        // get scintillator status
        int scint = 0;
        for (int j = 0; j < Nsc; ++j) {
            if (sx[j] <= 20 && sx[j] >= -20) {
                SET_BIT(scint, (sy[j] > 0 ? 1 : 3));
            } else if (sy[j] <= 20 && sy[j] >= -20) {
                SET_BIT(scint, (sy[j] > 0 ? 2 : 4));
            }
        }

        for (int j = 0; j < Nhits; ++j) {
            output << iev << sep << hx[j] << sep << hy[j] << sep << hz[j] << sep << hE[j] << sep
                   << match[j] << sep << scint << sep << cid << std::endl;
        }

        for (int j = 0; j < Ngun; ++j) {
            outgun << iev << sep << gx[j] << sep << gy[j] << sep << gz[j] << sep << gE[j] << sep
                   << gpid[j] << sep << gth[j] << sep << gph[j] << std::endl;
        }
    }

    std::cout <<"------[ ev " << tch->GetEntries() << " ]---"
              << "---[ " << timer.GetElapsedTimeStr() << " ]---"
              << "---[ " << timer.GetElapsedTime()/(double)tch->GetEntries() << " ms/ev ]------"
              << std::endl;
    output.close();
    outgun.close();
}


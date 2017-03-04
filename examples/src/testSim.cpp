//============================================================================//
// An example showing how to read and analyze simulation output from PRadSim  //
//                                                                            //
// Chao Peng                                                                  //
// 10/12/2016                                                                 //
//============================================================================//

#include "PRadDataHandler.h"
#include "PRadHyCalSystem.h"
#include "PRadDSTParser.h"
#include "PRadEvioParser.h"
#include "PRadBenchMark.h"
#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <string>
#include <vector>

using namespace std;

int main(int /*argc*/, char * /*argv*/ [])
{
    // simulation data is more like raw evio data with HyCal information only,
    // so we only need hycal system to connected to the handler
    PRadDataHandler *handler = new PRadDataHandler();
    PRadHyCalSystem *hycal = new PRadHyCalSystem("config/hycal.conf");

    handler->SetHyCalSystem(hycal);

    PRadBenchMark timer;

    //    handler->ReadFromEvio("/home/chao/Desktop/prad_001287.evio.0");
    // read simulation output
    handler->ReadFromEvio("/home/chao/geant4/PRadSim/output/simrun_47.evio");

    TFile *f = new TFile("testSim.root", "RECREATE");
    TTree *t = new TTree("T", "T");

    int N;
    double E[100], x[100], y[100]; // maximum number of clusters, 100 is enough
    // retrieve part of the cluster information
    t->Branch("NClusters", &N, "NClusters/I");
    t->Branch("ClusterE", &E, "ClusterE/D");
    t->Branch("ClusterX", &x[0], "ClusterX[NClusters]/D");
    t->Branch("ClusterY", &y[0], "ClusterY[NClusters]/D");


    for(auto &event : handler->GetEventData())
    {
        hycal->Reconstruct(event);
        auto &hits = hycal->GetDetector()->GetHits();
        N = (int)hits.size();
        for(size_t i = 0; i < hits.size(); ++i)
        {
            E[i] = hits[i].E;
            x[i] = hits[i].x;
            y[i] = hits[i].y;
        }
        t->Fill();
    }

    cout << "TIMER: Finished, took " << timer.GetElapsedTime() << " ms." << endl;
    cout << "Read " << handler->GetEventCount() << " events." << endl;

    f->cd();
    t->Write();
    f->Close();

    handler->WriteToDST("simrun_47.dst");
    return 0;
}

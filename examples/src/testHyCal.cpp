//============================================================================//
// An example showing how to use HyCAL clustering method to reconstruct data  //
//                                                                            //
// Chao Peng                                                                  //
// 11/12/2016                                                                 //
//============================================================================//

#include "PRadHyCalSystem.h"
#include "PRadDataHandler.h"
#include "PRadDSTParser.h"
#include "PRadBenchMark.h"
#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

using namespace std;

int main(int /*argc*/, char * /*argv*/ [])
{
    PRadHyCalSystem *sys = new PRadHyCalSystem("config/hycal.conf");
    PRadDSTParser *dst_parser = new PRadDSTParser();

    // show adc and tdc channels
    /*
    for(auto adc : sys->GetADCList())
    {
        cout << "ADC: "
             << *adc;
        if(adc->GetTDC())
            cout << *adc->GetTDC();
        cout << endl;
    }
    for(auto tdc : sys->GetTDCList())
    {
        cout << "TDC: "
             << setw(6) << tdc->GetName() << "  "
             << tdc->GetAddress()
             << endl;
        for(auto adc : tdc->GetChannelList())
        {
            cout << setw(6) << adc->GetName() << ",";
        }
        cout << endl;
    }
    */
    // show the modules
    /*
    PRadHyCalDetector *hycal = sys->GetDetector();
    for(auto module : hycal->GetModuleList())
    {
        cout << setw(4) << module->GetID() << ": " << *module;
        cout << setw(12) << module->GetNonLinearConst();
        if(module->GetChannel())
            cout << module->GetChannel()->GetAddress();
        cout << endl;
    }
    */

    // test reconstruction performance
//    dst_parser->OpenInput("/work/hallb/prad/replay/prad_001288.dst");
    dst_parser->OpenInput("prad_1310_select.dst");
    dst_parser->OpenOutput("prad_1310_leak.dst");
    TFile *f = new TFile("prad_1310_cluster.root", "RECREATE");
    TTree *t = new TTree("T", "T");

    int N;
    double E[100], x[100], y[100], corr[100], leak[100];
    // retrieve part of the cluster information
    t->Branch("NClusters", &N, "NClusters/I");
    t->Branch("ClusterE", &E, "ClusterE/D");
    t->Branch("ClusterX", &x[0], "ClusterX[NClusters]/D");
    t->Branch("ClusterY", &y[0], "ClusterY[NClusters]/D");
    t->Branch("NonLinearCorr", &corr[0], "NonLinearCorr[NClusters]/D");
    t->Branch("LeakageCorr", &leak[0], "LeakageCorr[NClusters]/D");

    PRadBenchMark timer;

    while(dst_parser->Read())
    {
        if(dst_parser->EventType() == PRadDSTParser::Type::event) {

            auto event = dst_parser->GetEvent();
            if(!event.is_physics_event())
                continue;

            sys->Reconstruct(event);
            auto hits = sys->GetDetector()->GetHits();
            N = hits.size();
            bool save = false;
            for(size_t i = 0; i < hits.size(); ++i)
            {
                const auto &hit = hits.at(i);
                x[i] = hit.x;
                y[i] = hit.y;
                E[i] = hit.E;
                corr[i] = hit.lin_corr;
                leak[i] = hit.E_leak;
                if(hit.E_leak > 100)
                    save = true;
            }
            t->Fill();

            if(save)
                dst_parser->WriteEvent(event);
        }
    }

    t->Write();
    f->Close();
    dst_parser->CloseInput();
    dst_parser->CloseOutput();
    cout << "TIMER: Finished, took " << timer.GetElapsedTime() << " ms" << endl;

    return 0;
}

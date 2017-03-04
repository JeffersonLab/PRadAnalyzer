//============================================================================//
// An example showing how to use GEM clustering method to reconstruct gem data//
//                                                                            //
// Chao Peng                                                                  //
// 10/07/2016                                                                 //
//============================================================================//

#include "PRadDataHandler.h"
#include "PRadDSTParser.h"
#include "PRadEvioParser.h"
#include "PRadBenchMark.h"
#include "PRadGEMSystem.h"
#include "PRadInfoCenter.h"
#include <iostream>
#include <string>
#include <vector>

#define PROGRESS_COUNT 10000

using namespace std;

int main(int argc, char *argv [])
{

    if(argc != 2) {
        cout << "usage: testGEM <dst_file>" << endl;
        return 0;
    }

    string dst_file = argv[1];

    PRadGEMSystem *gem_srs = new PRadGEMSystem("config/gem.conf");
    PRadDSTParser *dst_parser = new PRadDSTParser();

    // show the FECs and APVs
    for(auto &fec : gem_srs->GetFECList())
    {
        cout << "FEC " << fec->GetID() << ": " << endl;
        for(auto &apv : fec->GetAPVList())
        {
            cout << "     APV: " << apv->GetFECID() << ", " << apv->GetADCChannel() << endl;
        }
    }

    auto det_list = gem_srs->GetDetectorList();

    // show the planes and the strips that connected to APVs
    for(auto &detector: det_list)
    {
        cout << "Detector: " << detector->GetName() << endl;
        for(auto &plane : detector->GetPlaneList())
        {
            cout << "    " << "Plane: " << plane->GetName() << endl;
            for(auto &apv : plane->GetAPVList())
            {
                cout << "    " << "    "
                     << "APV: " << apv->GetPlaneIndex()
                     << ", " << apv->GetAddress();

                int min = apv->GetPlaneStripNb(0);
                int max = apv->GetPlaneStripNb(0);
                for(size_t i = 1; i < apv->GetTimeSampleSize(); ++i)
                {
                    int strip = apv->GetPlaneStripNb(i);
                    if(min > strip) min = strip;
                    if(max < strip) max = strip;
                }

                cout << ", " << min
                     << " ~ " << max
                     << endl;
            }
        }
    }

    //dst_parser->OpenInput("/work/hallb/prad/replay/prad_001288.dst");
    dst_parser->OpenInput(dst_file);

    // test reconstruction performance
    PRadBenchMark timer;

    int count = 0;
    double time = 0;

    while(dst_parser->Read())
    {
        if(dst_parser->EventType() == PRadDSTParser::Type::event) {
            auto event = dst_parser->GetEvent();
            if(!event.is_physics_event())
                continue;

            ++count;
            gem_srs->Reconstruct(event);

            if(count%PROGRESS_COUNT == 0) {
                double t = timer.GetElapsedTime();
                time += t;
                timer.Reset();

                cout <<"------[ ev " << count << " ]---"
                     << "---[ " << t << " ms ]---"
                     << "---[ " << time/(double)count << " ms/ev ]------"
                     << "\r" << flush;
            }
            // show strip clusters
            // detectors from GEM system

	    //-------------------------------------------------
	    /*
            for(auto &detector: gem_srs->GetDetectorList())
            {
                cout << "Detector: " << detector->GetName() << endl;
                // planes from a detector
                for(auto &plane : detector->GetPlaneList())
                {
                    cout << "    " << "Plane: " << plane->GetName() << endl;
                    // clusters from a plane
                    for(auto &cluster : plane->GetStripClusters())
                    {
                        cout << "    " << "    "
                             << "Cluster: "
                             << cluster.position << ", "
                             << cluster.peak_charge
                             << endl;
                        // hits from a cluster
                        for(auto &hit : cluster.hits)
                        {
                            cout << "    " << "    " << "     "
                                 << hit.strip << ", " << hit.charge << endl;
                        }
                    }
                }
            }
	    */
	    //-------------------------------------------------
        }
    }

    dst_parser->CloseInput();
    time += timer.GetElapsedTime();
    cout << endl;
    cout << "TIMER: Finished, read and reconstructed " << count << " events, "
         << "took " << time/1000. << " s."
         << endl;

    return 0;
}

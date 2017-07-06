//============================================================================//
// An example showing how to read-in NNID root file and reject the event that //
// is categorized as cosmic ray (or background) by the Neural Network         //
//                                                                            //
// Chao Peng                                                                  //
// 07/05/2017                                                                 //
//============================================================================//


#include "PRadHyCalSystem.h"
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
void RejectCosmic(const char *path, const char *nnid_file, const char *out_dir);

int main(int argc, char *argv[])
{

    if(argc != 4) {
        cerr << "usage: NNAE_Prepare <data_file> <root_file> <out_directory>" << endl;
        return -1;
    }

    RejectCosmic(argv[1], argv[2], argv[3]);
    return 0;
}

void RejectCosmic(const char *path, const char *nnid_file, const char *out_dir)
{
    string prad_path = getenv("PRAD_PATH");
    if(prad_path.size()) {
        prad_path += "/";
    }

    // read-in root file and save the rejected event numbers
    std::unordered_set<int> rej_list;
    TFile *myfile= new TFile(nnid_file,"READ");
    TTree *tree= (TTree*) myfile->Get("Tid");
    double nn_id;
    unsigned int EvNb;
    tree->SetBranchAddress("nn_id",&nn_id);
    tree->SetBranchAddress("EvNb",&EvNb);
    int entries= tree->GetEntries();
    for(int i = 0; i < entries; i++)
    {
        tree->GetEntry(i);
        // rejection criterium nn_id < 0
        if(nn_id < 0.) {
            rej_list.insert(EvNb);
        }
    }

    // reject events
    PRadDSTParser dst_parser, dst_parser2;
    PRadHyCalSystem sys(prad_path + "config/hycal.conf");
    sys.ChooseRun(path);
    PRadCoordSystem coord(prad_path + "database/coordinates.dat");
    coord.ChooseCoord(PRadInfoCenter::GetRunNumber());

    string name = ConfigParser::decompose_path(path).name;
    string save_path = ConfigParser::form_path(out_dir, name + "_sav.dst");
    string rej_path = ConfigParser::form_path(out_dir, name + "_rej.dst");

    dst_parser.OpenInput(path);
    dst_parser.OpenOutput(save_path);
    dst_parser2.OpenOutput(rej_path);
    int count = 0;
    PRadBenchMark timer;

    while(dst_parser.Read())
    {
        if(dst_parser.EventType() == PRadDSTParser::Type::event) {

            if(++count % PROGRESS_COUNT == 0) {
                cout << "----------event " << count
                     << "-------[ " << timer.GetElapsedTimeStr() << " ]------"
                     << "\r" << flush;
            }

            auto event = dst_parser.GetEvent();

            // in the rejected list
            if(rej_list.find((unsigned int)event.event_number) != rej_list.end()) {
                dst_parser2.WriteEvent(event);
            } else {
                dst_parser.WriteEvent(event);
            }

        } else if(dst_parser.EventType() == PRadDSTParser::Type::epics) {
            dst_parser.WriteEPICS();
        }
    }

    cout << "----------event " << count
         << "-------[ " << timer.GetElapsedTimeStr() << " ]------"
         << endl;

    dst_parser.CloseInput();
    dst_parser.CloseOutput();
    dst_parser2.CloseOutput();
}

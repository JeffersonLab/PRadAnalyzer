//============================================================================//
// List all the neighbor modules of the input module                          //
//                                                                            //
// Chao Peng                                                                  //
// 11/27/2017                                                                 //
//============================================================================//

#include "PRadHyCalDetector.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

using namespace std;

int main(int argc, char *argv [])
{
    if(argc < 2) {
        cout << "usage " << argv[0] << " <module_id>" << endl;
        return -1;
    }

    int id = atoi(argv[1]);

    string prad_root = getenv("PRAD_PATH");
    if(prad_root.empty()) prad_root = ".";

    PRadHyCalDetector hycal;
    hycal.ReadModuleList(prad_root + "/database/hycal_module.txt");

    auto module = hycal.GetModule(id);
    if(module == nullptr) {
        cout << "cannot find module " << id << endl;
        return -1;
    }

    cout << "list all the neighbors of module " << id << " and the quantized distance between them: \n";
    for(auto nbr : module->GetNeighbors())
    {
        cout << nbr->GetID()
             << ", dist = " << nbr.dist
             << ", dx = " << nbr.dx
             << ", dy = " << nbr.dy
             << endl;
    }

    return 0;
}

//============================================================================//
// Show the quantized distance between any two modules                        //
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
    if(argc < 3) {
        cout << "usage " << argv[0] << " <module_id1> <module_id2>" << endl;
        return -1;
    }

    int id1 = atoi(argv[1]);
    int id2 = atoi(argv[2]);

    string prad_root = getenv("PRAD_PATH");
    if(prad_root.empty()) prad_root = ".";

    PRadHyCalDetector hycal;
    hycal.ReadModuleList(prad_root + "/database/hycal_module.txt");

    auto m1 = hycal.GetModule(id1);
    auto m2 = hycal.GetModule(id2);

    if(m1 == nullptr || m2 == nullptr) {
        cout << "cannot find modules of " << id1 << ", " << id2 << endl;
        return -1;
    }

    cout << "the quantized distance between " << id1 << " and " << id2 << " is "
         << hycal.QuantizedDist(m1, m2) << endl;

    return 0;
}

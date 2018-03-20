//============================================================================//
// An example to test the performance of reconstruction methods               //
//                                                                            //
// Chao Peng                                                                  //
// 11/12/2016                                                                 //
//============================================================================//

#include "PRadBenchMark.h"
#include "ConfigParser.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

using namespace std;

int main(int argc, char *argv[])
{
    if(argc < 2) {
        cout << "usage: " << argv[0] << " <file>" << endl;
        return -1;
    }

    PRadBenchMark timer;

    ConfigParser parser;
    int count = 100;

    for(int i = 1; i <= count; ++i)
    {
        parser.OpenFile(argv[1]);
        while(parser.ParseLine()) {;}

        cout <<"------[ " << i << " ]---"
             << "---[ " << timer.GetElapsedTimeStr() << " ]---"
             << "---[ " << timer.GetElapsedTime()/(double)i << " ms/time ]------"
             << "\r" << flush;
    }

    cout <<"------[ ev " << count << " ]---"
         << "---[ " << timer.GetElapsedTimeStr() << " ]---"
         << "---[ " << timer.GetElapsedTime()/(double)count << " ms/time ]------"
         << endl;
    cout << "Finished." << endl;
}

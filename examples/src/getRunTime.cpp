//============================================================================//
// An example to get the run time from a dst file                             //
//                                                                            //
// Chao Peng                                                                  //
// 12/06/2017                                                                 //
//============================================================================//

#include "PRadDSTParser.h"
#include <iostream>

using namespace std;

int main(int argc, char *argv[])
{
    if(argc != 2) {
        cout << "usage: " << argv[0] << " <dst_file> " << endl;
        return -1;
    }

    PRadDSTParser dst;
    dst.OpenInput(argv[1]);
    if(!dst.ReadMap()) return -1;

    auto evmap = dst.GetInputMap().GetType(PRadDSTParser::Type::event);

    dst.Read(evmap.front());
    if(dst.EventType() != PRadDSTParser::Type::event) return -1;
    auto begin_time = dst.GetEvent().timestamp;

    dst.Read(evmap.back());
    if(dst.EventType() != PRadDSTParser::Type::event) return -1;
    auto end_time = dst.GetEvent().timestamp;

    // 4 ns to 1 sec
    double sec = double(end_time - begin_time)/2.5e8;

    std::cout << "run time: " << sec << " s." << std::endl;
}


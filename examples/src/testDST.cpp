//============================================================================//
// An example showing the features of PRadDSTParser                           //
//                                                                            //
// Chao Peng                                                                  //
// 07/21/2017                                                                 //
//============================================================================//

#include "PRadDSTParser.h"
#include "PRadBenchMark.h"
#include "ConfigParser.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#define PROGRESS_COUNT 10000

using namespace std;

void TestDST(const string &inf, const string &outf);
void ListEPICS(const string &inf);

int main(int argc, char *argv[])
{
    if(argc != 3) {
        cout << "usage: testDST <in_file> <out_file>" << endl;
        return 0;
    }

    TestDST(argv[1], argv[2]);
    ListEPICS(argv[2]);
}

void TestDST(const string &inf, const string &outf)
{
    PRadDSTParser dst_parser;

    dst_parser.OpenInput(inf);
    dst_parser.OpenOutput(outf);

    while(dst_parser.Read())
    {
        if(dst_parser.EventType() == PRadDSTParser::Type::event) {
            dst_parser.WriteEvent();
        } else if(dst_parser.EventType() == PRadDSTParser::Type::epics) {
            dst_parser.WriteEPICS();
        }
    }

    dst_parser.CloseInput();
    dst_parser.CloseOutput();
}

void ListEPICS(const string &inf)
{
    PRadDSTParser dst_parser;
    dst_parser.OpenInput(inf);

    if(dst_parser.ReadEventMap()) {
        auto dst_map = dst_parser.GetInputMap();

        for(auto &pos : dst_map.epics_pos)
        {
            if(dst_parser.Read(pos)) {
                if(dst_parser.EventType() == PRadDSTParser::Type::epics) {
                    cout << dst_parser.GetEPICSEvent().event_number << endl;
                } else {
                    cout << "Wrong type!" << endl;
                }
            } else {
                cout << "Read failure at " << pos << endl;
            }
        }
    } else {
        cout << "No event map information." << endl;
    }

    dst_parser.CloseInput();
}

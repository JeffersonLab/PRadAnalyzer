//============================================================================//
// An example showing the features of PRadDSTParser                           //
//                                                                            //
// Chao Peng                                                                  //
// 07/21/2017                                                                 //
//============================================================================//

#include "PRadDSTParser.h"
#include "DSTReaderV1.h"
#include "PRadBenchMark.h"
#include "ConfigParser.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include "DSTReaderV1.h"

#define PROGRESS_COUNT 10000

using namespace std;

void UpdateDST(const string &inf, const string &outf);
void ListEvents(const string &inf);

int main(int argc, char *argv[])
{
    if(argc != 3) {
        cout << "usage: updateDST <in_file> <out_file>" << endl;
        return 0;
    }

    UpdateDST(argv[1], argv[2]);
//    ListEvents(argv[2]);
}

void UpdateDST(const string &inf, const string &outf)
{
    DSTReaderV1 dst_old;
    PRadDSTParser dst_new;

    dst_old.OpenInput(inf);
    dst_new.OpenOutput(outf);

    while(dst_old.Read())
    {
        if(dst_old.EventType() == DSTReaderV1::Type::event) {
            dst_new.WriteEvent(dst_old.GetEvent());
        } else if(dst_old.EventType() == DSTReaderV1::Type::epics) {
            dst_new.WriteEPICS(dst_old.GetEPICSEvent());
        }
    }

    dst_old.CloseInput();
    dst_new.CloseOutput();
}

void ListEvents(const string &inf)
{
    PRadDSTParser dst_parser;
    dst_parser.OpenInput(inf);

    auto type = PRadDSTParser::Type::event;

    if(dst_parser.ReadMap()) {
        auto dst_map = dst_parser.GetInputMap();

        for(auto &pos : dst_map.GetType(type))
        {
            if(dst_parser.Read(pos)) {
                if(dst_parser.EventType() != type) {
                    cout << "Wrong type!" << endl;
                } else {
                    cout << dst_parser.GetEvent().event_number << endl;
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

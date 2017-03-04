//============================================================================//
// An example showing how to get beam charge information from one run         //
//                                                                            //
// Chao Peng                                                                  //
// 11/12/2016                                                                 //
//============================================================================//

#include "PRadDSTParser.h"
#include "PRadEventFilter.h"
#include "PRadInfoCenter.h"
#include "ConfigParser.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

#define PROGRESS_COUNT 10000

using namespace std;

void ChargeCount(const string &file, ofstream &out);

int main(int argc, char *argv[])
{
    if(argc < 2)
    {
        cout << "usage: beamChargeCount <file1> <file2> ..." << endl;
        return 0;
    }

    ofstream output("beam_charge.dat", ios::app);
    for(int i = 1; i < argc; ++i)
    {
        string file = argv[i];
        ChargeCount(file, output);
    }
    output.close();
}

void ChargeCount(const string &file, ofstream &out)
{
    cout << "Beam Charge Counting for " << "\"" << file << "\"."
         << endl;

    PRadInfoCenter::Instance().Reset();
    PRadInfoCenter::SetRunNumber(file);

    PRadDSTParser dst_parser;
    dst_parser.OpenInput(file);

    while(dst_parser.Read())
    {
        if(dst_parser.EventType() == PRadDSTParser::Type::event) {
            auto event = dst_parser.GetEvent();
            PRadInfoCenter::Instance().UpdateInfo(event);
        }
    }
    dst_parser.CloseInput();

    cout << "TIMER: Finished beam charge counting for run "
         << PRadInfoCenter::GetRunNumber() << "."
         << endl
         << "Total beam charge: "
         << PRadInfoCenter::GetBeamCharge() << " nC."
         << endl
         << "Average live time: "
         << PRadInfoCenter::GetLiveTime()*100. << "%."
         << endl
         << "Accepted beam charge: "
         << PRadInfoCenter::GetLiveBeamCharge() << " nC."
         << endl;

    out << setw(6) << PRadInfoCenter::GetRunNumber()
        << setw(12) << PRadInfoCenter::GetBeamCharge()
        << setw(12) << PRadInfoCenter::GetLiveTime()
        << setw(12) << PRadInfoCenter::GetLiveBeamCharge()
        << endl;
}

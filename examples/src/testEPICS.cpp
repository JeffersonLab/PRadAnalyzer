//============================================================================//
// An example showing how to get the EPICS information from data files        //
//                                                                            //
// Chao Peng                                                                  //
// 10/04/2016                                                                 //
//============================================================================//

#include "PRadDataHandler.h"
#include "PRadEPICSystem.h"
#include "PRadDSTParser.h"
#include "PRadEvioParser.h"
#include "PRadBenchMark.h"
#include <iostream>
#include <string>
#include <vector>

using namespace std;

void testEPICS(const char *path, const char *channel);

int main(int argc, char *argv[])
{
    if(argc < 3) {
        cout << "usage: testEPICS <file> <channel_name>" << endl;
        return 0;
    }

    testEPICS(argv[1], argv[2]);

}

void testEPICS(const char *path, const char *channel)
{
    PRadEPICSystem *epics = new PRadEPICSystem("config/epics_channels.conf");
    PRadDSTParser *dst_parser = new PRadDSTParser();

    dst_parser->OpenInput(path);

    while(dst_parser->Read())
    {
        if(dst_parser->EventType() == PRadDSTParser::Type::epics) {
            // save event into epics system, otherwise find epicsvalue won't work
            epics->AddEvent(dst_parser->GetEPICSEvent());
            epics->GetValue(channel);
            cout << channel << ", event = "
                 << epics->GetEventNumber() << ", value = "
                 << epics->GetValue(channel)
                 << endl;
        }
    }

    dst_parser->CloseInput();

    delete epics;
    delete dst_parser;

}

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

int main(int /*argc*/, char * /*argv*/ [])
{
    // either read the epics channel map from text file
    PRadEPICSystem *epics = new PRadEPICSystem("config/epics_channels.conf");
    PRadDSTParser *dst_parser = new PRadDSTParser();

    // or read the saved epics channel information from dst file, 
    // to do so, you need set a handler with the epics system to the dst parser,
    // as
    //PRadDataHandler *handler = new PRadDataHandler;
    //handler->SetEPICSystem(epics);
    //dst_parser->SetHandler(handler);

    PRadBenchMark timer;
    // here shows an example how to read DST file while not saving all the events
    // in memory
    dst_parser->OpenInput("/work/hallb/prad/replay/prad_001288.dst");

    int count = 0;
    while(dst_parser->Read() && count < 300)
    {
        if(dst_parser->EventType() == PRadDSTParser::Type::event) {
            ++count;
            auto event = dst_parser->GetEvent();
            cout << event.event_number << ", energy is "
                 << epics->FindValue(event.event_number, "MBSY2C_energy") << endl;
        } else if(dst_parser->EventType() == PRadDSTParser::Type::epics) {
            // save event into epics system, otherwise find epicsvalue won't work
            epics->AddEvent(dst_parser->GetEPICSEvent());
            cout << dst_parser->GetEPICSEvent().event_number << ", "
                 << dst_parser->GetEPICSEvent().values.at(0) << endl;
        }
    }

    dst_parser->CloseInput();

    cout << "TIMER: Finished, took " << timer.GetElapsedTime() << " ms" << endl;

    return 0;
}

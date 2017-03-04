//============================================================================//
// An example showing how to get the LMS data from the replayed DST file      //
//                                                                            //
// Chao Peng                                                                  //
// 02/27/2016                                                                 //
//============================================================================//

#include "PRadHyCalSystem.h"
#include "PRadDSTParser.h"
#include "PRadBenchMark.h"
#include <iostream>
#include <string>
#include <vector>

using namespace std;

int main(int /*argc*/, char * /*argv*/ [])
{
    PRadHyCalSystem *hycal = new PRadHyCalSystem;
    PRadDSTParser *dst_parser = new PRadDSTParser();

    // read configuration files
    hycal->Configure("config/hycal.conf");

    PRadBenchMark timer;

    // here shows an example how to read DST file while not saving all the events
    // in memory
    dst_parser->OpenInput("/work/hallb/prad/replay/prad_001288.dst");

    int count = 0;
    while(dst_parser->Read() && count < 20000)
    {
        if(dst_parser->EventType() != PRadDSTParser::Type::event)
            continue;

        ++count;
        hycal->FillHists(dst_parser->GetEvent());
    }

    dst_parser->CloseInput();
    cout << "TIMER: Finished reading, took "
         << timer.GetElapsedTime() << " ms"
         << endl;

    timer.Reset();

    // channel name, histogram name, fit function, min range, max range, verbose
    auto pars1 = hycal->FitHist("W1115", "Pedestal", "gaus", 0, 8200, true);
    auto pars2 = hycal->FitHist("W1115", "LMS", "gaus", 0, 8200, true);

    // for Gaussian, par0 is amplitude, par1 is mean, par2 is sigma
    cout << "Pedestal for W1115: mean " << pars1.at(1)
         << ", sigma " << pars1.at(2) << endl;
    cout << "LMS for W1115: mean " << pars2.at(1)
         << ", sigma " << pars2.at(2) << endl;

    cout << "TIMER: Finished fitting, took "
         << timer.GetElapsedTime() << " ms"
         << endl;

//    handler->WriteToDST("prad_001323_0-10.dst");
    //handler->PrintOutEPICS();
    return 0;
}

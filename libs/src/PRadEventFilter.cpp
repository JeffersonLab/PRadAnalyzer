//============================================================================//
// An class used to cut off bad events for PRad                               //
// Now it reads a list of bad events, expecting more criterion implemented    //
//                                                                            //
// Maxime Levillain, Chao Peng                                                //
// 10/17/2016                                                                 //
//============================================================================//

#include "PRadEventFilter.h"
#include "ConfigParser.h"

PRadEventFilter::PRadEventFilter(const std::string &path)
{
    LoadBadEventsList(path);
}

PRadEventFilter::~PRadEventFilter()
{
    // place holder
}

void PRadEventFilter::LoadBadEventsList(const std::string &path, bool clear_exist)
{
    if(clear_exist)
        bad_events_list.clear();

    ConfigParser c_parser;
    // remove *, space and tab at both ends of each element
    c_parser.SetWhiteSpace(" \t*");

    if (!c_parser.ReadFile(path)) {
        std::cerr << "PRad Event Filter: Cannot open bad event list "
                  << "\"" << path << "\"."
                  << "No bad events list loaded."
                  << std::endl;
        return;
    }

    int val1, val2;

    while(c_parser.ParseLine())
    {
        if(c_parser.NbofElements() != 2)
            continue;

        c_parser >> val1 >> val2;
        bad_events_list.emplace_back(val1, val2);
    }

}

void PRadEventFilter::ClearBadEventsList()
{
    bad_events_list.clear();
}

bool PRadEventFilter::IsBadEvent(const EventData &event)
const
{
    // loop the whole list
    for(auto &interval : bad_events_list)
    {
        if((event.event_number >= interval.begin) &&
           (event.event_number <= interval.end))
        {
            // found in the listed bad events intervals
            return true;
        }
    }

    // not in the list
    return false;
}

bool PRadEventFilter::IsBadPeriod(const EventData &begin, const EventData &end)
const
{
    // loop the whole list
    for(auto &interval : bad_events_list)
    {
        // firstly, either begin or end itself is a bad event
        if(IsBadEvent(begin) || IsBadEvent(end))
            return true;

        // secondly, if the whole period contains a bad period
        if((begin.event_number < interval.begin) &&
           (end.event_number > interval.end))
            return true;
    }

    return false;
}


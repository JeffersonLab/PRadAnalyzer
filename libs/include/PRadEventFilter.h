#ifndef PRAD_EVENT_FILTER
#define PRAD_EVENT_FILTER

#include <vector>
#include <string>
#include "PRadEventStruct.h"

class PRadEventFilter
{
private:
    struct ev_interval
    {
        int begin;
        int end;

        ev_interval() {};
        ev_interval(int b, int e) : begin(b), end(e) {};
    };

public:
    PRadEventFilter(const std::string &path);
    virtual ~PRadEventFilter();

    void LoadBadEventsList(const std::string &path, bool clear_exist = true);
    void ClearBadEventsList();
    bool IsBadEvent(const EventData &event) const;
    bool IsBadPeriod(const EventData &begin, const EventData &end) const;

private:
    std::vector<ev_interval> bad_events_list;
};

#endif

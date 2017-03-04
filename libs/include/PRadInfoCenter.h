#ifndef PRAD_INFO_CENTER_H
#define PRAD_INFO_CENTER_H

#include <string>
#include "datastruct.h"
#include "PRadEventStruct.h"

class PRadInfoCenter
{
public:
    static PRadInfoCenter &Instance()
    {
        static PRadInfoCenter instance;

        return instance;
    }

    PRadInfoCenter(const PRadInfoCenter &)  = delete;
    void operator=(const PRadInfoCenter &)  = delete;

    static bool SetRunNumber(int run);
    static bool SetRunNumber(const std::string &path);
    static int GetRunNumber();
    static double GetBeamCharge();
    static double GetLiveBeamCharge();
    static double GetLiveTime();

    void Reset();
    void UpdateInfo(const EventData &event);
    void SetRunInfo(const RunInfo &info) {run_info = info;};
    void SetOnlineInfo(const OnlineInfo &info) {online_info = info;};

    const RunInfo &GetRunInfo() const {return run_info;};
    const OnlineInfo &GetOnlineInfo() const {return online_info;};

private:
    RunInfo run_info;
    OnlineInfo online_info;
    double live_scaled_charge;

    PRadInfoCenter();
};

#endif

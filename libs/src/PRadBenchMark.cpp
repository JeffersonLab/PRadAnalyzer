//============================================================================//
// A timer class to test the performance                                      //
//                                                                            //
// Chao Peng                                                                  //
// 02/12/2016                                                                 //
//============================================================================//

#include "PRadBenchMark.h"



PRadBenchMark::PRadBenchMark()
{
    Reset();
}

PRadBenchMark::~PRadBenchMark()
{
    // place holder
}

void PRadBenchMark::Reset()
{
    time_point = std::chrono::high_resolution_clock::now();
}

unsigned int PRadBenchMark::GetElapsedTime()
const
{
    auto time_end = std::chrono::high_resolution_clock::now();
    auto int_ms = std::chrono::duration_cast<std::chrono::milliseconds>(time_end - time_point);
    return int_ms.count();
}

std::string PRadBenchMark::GetElapsedTimeStr(bool show_msec)
const
{
    int time = GetElapsedTime();
    int t_sec = time/1000;
    int hour = t_sec/3600;
    int min = (t_sec%3600)/60;
    int sec = (t_sec%3600)%60;


    std::string time_str;
    time_str.reserve(30);

    if(hour < 10)
        time_str += "0";
    time_str += std::to_string(hour) + ":";
    if(min < 10)
        time_str += "0";
    time_str += std::to_string(min) + ":";
    if(sec < 10)
        time_str += "0";
    time_str += std::to_string(sec);

    if(show_msec) {
        time_str += ".";
        int msec = time%1000;
        if(msec < 10)
            time_str += "0";
        if(msec < 100)
            time_str += "0";
        time_str += std::to_string(msec);
    }

    return time_str;
}


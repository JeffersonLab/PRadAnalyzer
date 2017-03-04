//============================================================================//
// A simple benchmark class                                                   //
//                                                                            //
// Chao Peng                                                                  //
// 02/17/2016                                                                 //
//============================================================================//

#ifndef PRAD_BENCH_MARK_H
#define PRAD_BENCH_MARK_H

#include <chrono>
#include <string>

class PRadBenchMark
{
public:
    PRadBenchMark();
    virtual ~PRadBenchMark();

    void Reset();
    unsigned int GetElapsedTime() const;
    std::string GetElapsedTimeStr(bool show_msec = true) const;

private:
    std::chrono::time_point<std::chrono::high_resolution_clock> time_point;
};

#endif

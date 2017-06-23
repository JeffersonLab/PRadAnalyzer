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
#include <vector>

class PRadBenchMark
{
public:
    PRadBenchMark();
    virtual ~PRadBenchMark();

    void Reset();
    void Save();
    void ClearSaved() {saved.clear();}
    unsigned int GetElapsedTime() const;
    std::string GetElapsedTimeStr(bool show_msec = true) const;
    unsigned int GetLastSaved() const {if(saved.empty()) return 0; else return saved.back();}
    const std::vector<unsigned int> &GetAllSaved() const {return saved;}

private:
    std::chrono::time_point<std::chrono::high_resolution_clock> time_point;
    std::vector<unsigned int> saved;
};

#endif

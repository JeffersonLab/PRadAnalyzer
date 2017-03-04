#ifndef PRAD_EPIC_SYSTEM_H
#define PRAD_EPIC_SYSTEM_H

#include <string>
#include <vector>
#include <deque>
#include <unordered_map>
#include "PRadEventStruct.h"
#include "PRadException.h"

// epics channel
struct EPICSChannel
{
    std::string name;
    uint32_t id;
    float value;

    EPICSChannel(const std::string &n, const uint32_t &i, const float &v)
    : name(n), id(i), value(v)
    {};
};

class PRadEPICSystem
{
public:
    PRadEPICSystem(const std::string &path = "");
    virtual ~PRadEPICSystem();

    void Reset();
    void ReadMap(const std::string &path);
    void SaveMap(const std::string &path) const;
    void AddChannel(const std::string &name);
    void AddChannel(const std::string &name, uint32_t id, float value);
    void UpdateChannel(const std::string &name, const float &value);
    void AddEvent(EpicsData &&data);
    void AddEvent(const EpicsData &data);
    void FillRawData(const char *buf);
    void SaveData(const int &event_number, bool online = false);

    std::vector<EPICSChannel> GetSortedList() const;
    const std::vector<float> &GetValues() const {return epics_values;};
    float GetValue(const std::string &name) const;
    int GetChannel(const std::string &name) const;
    const EpicsData &GetEvent(const unsigned int &index) const throw(PRadException);
    const std::deque<EpicsData> &GetEventData() const {return epics_data;};
    unsigned int GetEventCount() const {return epics_data.size();};
    float FindValue(int event_number, const std::string &name) const;
    int FindEvent(int event_number) const;


private:
    // data related
    std::unordered_map<std::string, uint32_t> epics_map;
    std::vector<float> epics_values;
    std::deque<EpicsData> epics_data;
};

#endif

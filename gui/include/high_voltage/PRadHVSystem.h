#ifndef PRAD_HV_SYSTEM_H
#define PRAD_HV_SYSTEM_H

#include "CAENHVSystem.h"
#include <vector>
#include <string>
#include <thread>
#include <mutex>
#include <unordered_map>
#include "datastruct.h"

class PRadEventViewer;
class PRadHVSystem
{
public:
    struct Voltage
    {
       float Vmon;
       float Vset;
       bool ON;
       Voltage() : Vmon(0), Vset(0), ON(false) {};
       Voltage(float vm, float vs, bool o = false) : Vmon(vm), Vset(vs), ON(o) {};
    };


public:
    PRadHVSystem(PRadEventViewer *p);
    virtual ~PRadHVSystem();
    void AddCrate(const std::string &name,
                  const std::string &ip,
                  const unsigned char &id,
                  const CAENHV::CAENHV_SYSTEM_TYPE_t &type = CAENHV::SY1527,
                  const int &linkType = LINKTYPE_TCPIP,
                  const std::string &username = "admin",
                  const std::string &password = "admin");
    void Connect();
    void Disconnect();
    void StartMonitor();
    void StopMonitor();
    void ReadVoltage();
    void CheckStatus();
    void SaveCurrentSetting(const std::string &path);
    void RestoreSetting(const std::string &path);
    void SetVoltage(const ChannelAddress &addr, const float &Vset);
    void SetPower(const bool &on_off);
    void SetPower(const ChannelAddress &addr, const bool &on_off);

    CAEN_Crate *GetCrate(const std::string &name) const;
    CAEN_Crate *GetCrate(int id) const;
    CAEN_Board *GetBoard(const std::string &name, int slot) const;
    CAEN_Board *GetBoard(int id, int slot) const;
    CAEN_Channel *GetChannel(const std::string &name, int slot, int channel) const;
    CAEN_Channel *GetChannel(int id, int slot, int channel) const;
    CAEN_Channel *GetChannel(const ChannelAddress &addr) const;
    Voltage GetVoltage(const std::string &name, int slot, int channel) const;
    Voltage GetVoltage(int id, int slot, int channel) const;
    Voltage GetVoltage(const ChannelAddress &addr) const;

private:
    PRadEventViewer *console;
    std::vector<CAEN_Crate*> crateList;
    volatile bool alive;
    std::thread queryThread;
    std::mutex locker;
    std::unordered_map<int, CAEN_Crate*> crate_id_map;
    std::unordered_map<std::string, CAEN_Crate*> crate_name_map;
    void queryLoop();
};

#endif

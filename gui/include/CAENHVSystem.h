#ifndef CAEN_HV_SYSTEM_H
#define CAEN_HV_SYSTEM_H

#include "CAENHVWrapper.h"
#include <string>
#include <vector>
#include <iostream>

class CAEN_Channel;
class CAEN_Board;
class CAEN_Crate;

std::ostream &operator<<(std::ostream &os, CAEN_Board const &b);
void CAEN_ShowError(const std::string &prefix, const int &err, const bool &ShowSuccess = false);
float CAEN_VoltageLimit(const std::string &name);
void CAEN_ShowChError(const std::string &n, const unsigned int &bitmap);

class CAEN_Channel
{
private:
    CAEN_Board *mother;
    unsigned short channel;
    std::string name;
    bool on_off;
    float Vmon;
    float Vset;
    float limit;

public:
    // constructor
    CAEN_Channel(CAEN_Board *m)
    : mother(m), channel(-1), name(""), on_off(false), Vmon(0), Vset(0),
      limit(CAEN_VoltageLimit(name))
    {};
    CAEN_Channel(CAEN_Board *m, const unsigned short &c, const std::string &n)
    : mother(m), channel(c), name(n), on_off(false), Vmon(0), Vset(0),
      limit(CAEN_VoltageLimit(name))
    {};
    CAEN_Channel(CAEN_Board *m, const unsigned short &c, const std::string &n,
                 const bool &o, const float &vm, const float vs)
    : mother(m), channel(c), name(n), on_off(o), Vmon(vm), Vset(vs),
      limit(CAEN_VoltageLimit(name))
    {};

    virtual ~CAEN_Channel();

    void SetPower(const bool &on);
    void SetVoltage(const float &v);
    void SetName(const std::string &n);
    void SetLimit(const float &l);
    void ReadVoltage();
    void CheckStatus();
    void UpdateVoltage(const bool &pw, const float &vm, const float &vs);
    const std::string &GetName() {return name;};
    const float &GetVSet() {return Vset;};
    const float &GetVMon() {return Vmon;};
    const bool &IsTurnedOn() {return on_off;};
    const unsigned short &GetChannel() {return channel;};
    CAEN_Board *GetMOther() {return mother;};
};

class CAEN_Board
{
private:
    CAEN_Crate *mother;
    std::string model;
    std::string desc;
    unsigned short slot;
    unsigned short nChan;
    unsigned short serNum;
    unsigned char fmwLSB;
    unsigned char fmwMSB;
    int primary;
    std::vector<CAEN_Channel*> channelList;

public:
    // constructor
    CAEN_Board(CAEN_Crate *mo)
    : mother(mo), slot(-1), nChan(0), serNum(0), fmwLSB(0), fmwMSB(0), primary(-1)
    {};
    CAEN_Board(CAEN_Crate *mo, std::string m, std::string d, unsigned short s, unsigned short n,
               unsigned short ser, unsigned char lsb, unsigned char msb)
    : mother(mo), model(m), desc(d), slot(s), nChan(n), serNum(ser), fmwLSB(lsb), fmwMSB(msb), primary(-1)
    {};
    CAEN_Board(CAEN_Crate *mo, char* m, char* d, unsigned short s, unsigned short n,
               unsigned short ser, unsigned char lsb, unsigned char msb)
    : mother(mo), model(m), desc(d), slot(s), nChan(n), serNum(ser), fmwLSB(lsb), fmwMSB(msb), primary(-1)
    {};

    virtual ~CAEN_Board();

    void ReadBoardMap();
    void ReadVoltage();
    void CheckStatus();
    void SetPower(const bool &on_off);
    void SetPower(const std::vector<unsigned int> &on_off);
    void SetVoltage(const std::vector<float> &Vset);
    int GetHandle();
    const unsigned short &GetSlot() {return slot;};
    CAEN_Crate *GetMother() {return mother;};
    CAEN_Channel *GetPrimaryChannel();
    CAEN_Channel *GetChannel(int i);
    std::vector<CAEN_Channel*> &GetChannelList() {return channelList;};
    const std::string &GetModel() {return model;};
    const std::string &GetDescription() {return desc;};
    const unsigned short &GetSize() {return nChan;};
    const unsigned short &GetSerialNum() {return serNum;};
    unsigned short GetFirmware();
    const int &GetPrimaryChannelNumber() {return primary;};
 };

class CAEN_Crate
{
private:
    unsigned char id;
    std::string name;
    std::string ip;
    CAENHV::CAENHV_SYSTEM_TYPE_t sys_type;
    int link_type;
    std::string username;
    std::string password;
    int handle;
    bool mapped;
    unsigned short slot_map[50];
    std::vector<CAEN_Board*> boardList;

public:
    // constructor
    CAEN_Crate(const unsigned char &i,
               const std::string &n,
               const std::string &p,
               const CAENHV::CAENHV_SYSTEM_TYPE_t &type,
               const int &link,
               const std::string &user,
               const std::string &pwd)
    : id(i), name(n), ip(p), sys_type(type), link_type(link),
      username(user), password(pwd), handle(-1), mapped(false)
    {};

    virtual ~CAEN_Crate();

    bool Initialize();
    bool DeInitialize();
    void ReadCrateMap();
    void HeartBeat();
    void Clear();
    void ReadVoltage();
    void CheckStatus();
    void SetPower(const bool &on_off);
    const int &GetHandle() {return handle;};
    const std::string &GetName() {return name;};
    const std::string &GetIP() {return ip;};
    std::vector<CAEN_Board*> &GetBoardList() {return boardList;};
    CAEN_Board *GetBoard(const unsigned short &slot);
};

#endif

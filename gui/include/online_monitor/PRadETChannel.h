#ifndef PRAD_ET_CHANNEL_H
#define PRAD_ET_CHANNEL_H

#include <unordered_map>
#include <string>
#include <stdint.h>
#include "et.h"
#include "PRadException.h"
#include "PRadETStation.h"

#define ET_CHUNK_SIZE 500

class PRadETChannel
{

public:
    class Configuration
    {
    public:
        Configuration();
        virtual ~Configuration();

        et_openconfig &Get() {return config;}

        // wrapper functions
        void Initialize();
        void SetWait(int val);
        void SetTimeOut(struct timespec val);
        void SetHost(const char *val);
        void SetCast(int val);
        void SetTTL(int val);
        void SetPort(unsigned short val);
        void SetServerPort(unsigned short val);
        void AddBroadCast(const char *val);
        void RemoveBroadCast(const char *val);
        void AddMultiCast(const char *val);
        void RemoveMultiCast(const char *val);
        void SetPolicy(int val);
        void SetMode(int val);
        void SetDebugDefault(int val);
        void SetInterface(const char *val);
        void SetTCP(int rBufSize, int sBufSize, int noDelay);

    private:
        et_openconfig config;
    };


public:
    PRadETChannel(size_t size = 1048576);
    virtual ~PRadETChannel();
    void Open(const char *ipAddr, int tcpPort, const char *etFile) throw(PRadException);
    void NewStation(const std::string &name);
    void SwitchStation(const std::string &name);
    void RemoveStation(const std::string &name) throw(PRadException);
    void AttachStation() throw(PRadException);
    void DetachStation();
    void ForceClose();
    bool Read() throw(PRadException);
    void *GetBuffer() {return (void*) buffer;}
    size_t GetBufferLength() {return bufferSize;}
    Configuration &GetConfig() {return config;}
    et_sys_id &GetID() {return et_id;}
    PRadETStation *GetCurrentStation() {return curr_stat;}
    PRadETStation *GetStation(const std::string &name);

private:
    Configuration config;
    PRadETStation *curr_stat;
    std::unordered_map<std::string, PRadETStation*> stations;
    et_sys_id et_id;
    et_event *etEvent;
    uint32_t *buffer;
    size_t bufferSize;
    void copyEvent();
};

#endif

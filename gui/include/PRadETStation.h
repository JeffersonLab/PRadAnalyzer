#ifndef PRAD_ET_STATION_H
#define PRAD_ET_STATION_H

#include "et.h"
#include "PRadException.h"
#include <string>

class PRadETChannel;

class PRadETStation
{
public:
    class Configuration
    {
    public:
        Configuration();
        virtual ~Configuration();

        et_statconfig &Get() {return config;};

        //wrapper functions
        void Initialize();
        void SetBlock(int val);
        void SetFlow(int val);
        void SetSelect(int val);
        void SetUser(int val);
        void SetRestore(int val);
        void SetCUE(int val);
        void SetPrescale(int val);
        void SetSelectWords(int val[]);
        void SetFunction(const char *val);
        void SetLib(const char *val);
        void SetClass(const char *val);

    private:
        et_statconfig config;
    };

public:
    PRadETStation(PRadETChannel *p, std::string n, int mode = 2);
    virtual ~PRadETStation();
    Configuration &GetConfig() {return config;};
    et_stat_id &GetID() {return station_id;};
    et_att_id &GetAttachID() {return attach_id;};
    std::string GetName() {return name;};
    void PreSetting(int mode) throw(PRadException);
    void Create() throw(PRadException);
    void Attach() throw(PRadException);
    void Detach() throw(PRadException);
    void Remove() throw(PRadException);

private:
    PRadETChannel *et_system;
    std::string name;
    et_att_id attach_id;
    et_stat_id station_id;
    Configuration config;
};

#endif

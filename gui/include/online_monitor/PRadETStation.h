#ifndef PRAD_ET_STATION_H
#define PRAD_ET_STATION_H

#include "et.h"
#include "PRadException.h"
#include <memory>
#include <string>
#include <vector>

class PRadETChannel;

class PRadETStation
{
public:
    class Configuration
    {
    public:
        Configuration() {
            void *ptr;
            et_station_config_init(&ptr);
            conf = std::shared_ptr<void>(ptr, [] (void *p) { et_station_config_destroy(p); });
        }

        et_statconfig Get() {return conf.get();}

        //wrapper functions
        void SetBlock(int val) { et_station_config_setblock(conf.get(), val); }
        int GetBlock() { int val; et_station_config_getblock(conf.get(), &val); return val; }

        void SetFlow(int val) { et_station_config_setflow(conf.get(), val); }
        int GetFlow() { int val; et_station_config_getflow(conf.get(), &val); return val; }

        void SetSelect(int val) { et_station_config_setselect(conf.get(), val); }
        int GetSelect() {int val; et_station_config_getselect(conf.get(), &val); return val; }

        void SetUser(int val) { et_station_config_setuser(conf.get(), val); }
        int GetUser() { int val; et_station_config_getuser(conf.get(), &val); return val; }

        void SetRestore(int val) { et_station_config_setrestore(conf.get(), val); }
        int GetRestore() { int val; et_station_config_getrestore(conf.get(), &val); return val; }

        void SetCUE(int val) { et_station_config_setcue(conf.get(), val); }
        int GetCUE() { int val; et_station_config_getcue(conf.get(), &val); return val; }

        void SetPrescale(int val) { et_station_config_setprescale(conf.get(), val); }
        int GetPrescale() { int val; et_station_config_getprescale(conf.get(), &val); return val; }

        void SetSelectWords(std::vector<int> val) { et_station_config_setselectwords(conf.get(), &val[0]); }
        std::vector<int> GetSelectWords() { int vals[4]; et_station_config_getselectwords(conf.get(), vals); return std::vector<int>(vals, vals + 4); }

        void SetFunction(const std::string &val) { et_station_config_setfunction(conf.get(), val.c_str()); }
        std::string GetFunction() { char val[1024]; et_station_config_getfunction(conf.get(), val); return std::string(val); }

        void SetLib(const std::string &val) { et_station_config_setlib(conf.get(), val.c_str()); }
        std::string GetLib() { char val[1024]; et_station_config_getlib(conf.get(), val); return std::string(val); }

        void SetClass(const std::string &val) { et_station_config_setclass(conf.get(), val.c_str()); }
        std::string GetClass() { char val[1024]; et_station_config_getclass(conf.get(), val); return std::string(val); }

    private:
        std::shared_ptr<void> conf;
    };

public:
    PRadETStation(PRadETChannel *p, std::string n, int mode = 2);
    virtual ~PRadETStation();
    Configuration &GetConfig() {return config;}
    et_stat_id &GetID() {return station_id;}
    et_att_id &GetAttachID() {return attach_id;}
    std::string GetName() {return name;}
    void PreSetting(int mode);
    void Create();
    void Attach();
    void Detach();
    void Remove();

private:
    PRadETChannel *et_system;
    std::string name;
    et_att_id attach_id;
    et_stat_id station_id;
    Configuration config;
};

#endif

#ifndef CONFIG_OBJECT_H
#define CONFIG_OBJECT_H

#include <string>
#include <utility>
#include <vector>
#include <unordered_map>
#include "ConfigParser.h"


#define CONF_CONN(val, str, def, warn) val=GetConfigValue<decltype(val)>(str, def, warn)
#define CONF_CONN2(val, def, warn) val=GetConfiValue<decltype(val)>(str, def, warn)
#define GET_CONF(obj, val, str, def, warn) val=obj.GetConfiValue<decltype(val)>(str, def, warn)
#define GET_CONF2(obj, val, def, warn) val=GetConfiValue<decltype(val)>(#val, def, warn)


class ConfigObject
{
public:
    // constructor, desctructor
    ConfigObject(const std::string &spliiter = ":=", const std::string &ignore = " _\t", bool case_ins = true);

    virtual ~ConfigObject();

    // public member functions
    void ClearConfig();
    void SaveConfig(const std::string &path = "") const;
    void ListKeys() const;
    bool HasKey(const std::string &name) const;

    bool ReadConfigFile(const std::string &path);
    void ReadConfigString(const std::string &content);
    void SetConfigValue(const std::string &var_name, const ConfigValue &c_value);
    void SetIgnoreChars(const std::string &ignore) {ignore_chars = ignore;}
    void SetSplitChars(const std::string &splitter) {split_chars = splitter;}
    void SetReplacePair(const std::string &open, const std::string &close)
    {
        replace_pair = std::make_pair(open, close);
    }

    // get members
    ConfigValue GetConfigValue(const std::string &var_name) const;
    template<typename T>
    T GetConfigValue(const std::string &var_name)
    const
    {
        return GetConfigValue(var_name).Convert<T>();
    }

    ConfigValue GetConfigValue(const std::string &var_name, const ConfigValue &def_value, bool verbose = true);
    template<typename T>
    T GetConfigValue(const std::string &var_name, const T &val, bool verbose = true)
    {
        return GetConfigValue(var_name, ConfigValue(val), verbose).Convert<T>();
    }

    const std::string &GetConfigPath() const {return config_path;}
    const std::string &GetSplitChars() const {return split_chars;}
    const std::string &GetSpaceChars() const {return ignore_chars;}
    const std::pair<std::string, std::string> &GetReplacePair() const {return replace_pair;}
    std::vector<std::string> GetKeyList() const;


    // functions that to be overloaded
    virtual void Configure(const std::string &path = "");

protected:
    // protected member functions
    void reform(std::string &input, const std::string &open, const std::string &close) const;
    std::string formKey(const std::string &raw_key) const;

private:
    void parserProcess(ConfigParser &p, const std::string &source);
    void parseControl(const std::string &control_word);
    void parseTerm(std::string &&var_name, std::string &&var_value);

protected:
    std::string split_chars;
    std::string ignore_chars;
    bool case_insensitive;
    std::pair<std::string, std::string> replace_pair;
    std::string config_path;
    std::unordered_map<std::string, std::string> config_map;

    // return this reference when there is no value found in the map
    ConfigValue __empty_value;
};

#endif

//============================================================================//
// A class based on the support from ConfigParser and ConfigValue             //
// It provides a simple way to read text file as configuration file, and read //
// or modify a parameter in the inherited class                               //
// The Configure() function should be overloaded according to specialized     //
// requirements, and be called after the parameters being configured          //
//                                                                            //
// 10/31/2016                                                                 //
//============================================================================//

#include <iostream>
#include <iomanip>
#include <algorithm>
#include "ConfigObject.h"



//============================================================================//
// Constructor, Destructor                                                    //
//============================================================================//

// constructor
ConfigObject::ConfigObject(const std::string &splitter,
                           const std::string &ignore)
: split_chars(splitter), ignore_chars(ignore), __empty_value("")
{
    // set default replace bracket
    replace_pair = std::make_pair("{", "}");
}

// destructor
ConfigObject::~ConfigObject()
{
    // place holder
}


//============================================================================//
// Public Member Function                                                     //
//============================================================================//

// configure the cluster method
void ConfigObject::Configure(const std::string &path)
{
    // save the path
    config_path = path;

    // clear the map
    config_map.clear();

    // read configuration file in
    ReadConfigFile(path);
}

// clear all the loaded configuration values
void ConfigObject::ClearConfig()
{
    config_path = "";
    config_map.clear();
}

// read configuration file and build the configuration map
bool ConfigObject::ReadConfigFile(const std::string &path)
{
    ConfigParser c_parser(split_chars); // self-defined splitters

    if(!c_parser.ReadFile(path)) {
        std::cerr << "Cannot open configuration file "
                  << "\"" << path << "\""
                  << std::endl;
        return false;
    }

    // current directory
    std::string this_dir = ConfigParser::decompose_path(path).dir;

    while(c_parser.ParseLine())
    {
        // possible control words
        if(c_parser.NbofElements() == 1) {
            std::string control = c_parser.TakeFirst();
            size_t pos = control.find("{THIS_DIR}");
            if(pos != std::string::npos)
                control.replace(pos, 10, this_dir);
            parseControl(control);
        }
        // var_name and var_value
        else if (c_parser.NbofElements() == 2) {
            std::string var_name, key, var_value;
            c_parser >> var_name >> var_value;
            size_t pos = var_value.find("{THIS_DIR}");
            if(pos != std::string::npos)
                var_value.replace(pos, 10, this_dir);
            parseTerm(std::move(var_name), std::move(var_value));
        }
        // unsupported format
        else {
            std::cout << "Warning: Unsupported format in file "
                      << "\"" << path << "\" "
                      << "at line " << c_parser.LineNumber()
                      << std::endl
                      << "\"" << c_parser.CurrentLine() << "\""
                      << std::endl;
        }
    }
    return true;
}

// check if a certain term is configured
bool ConfigObject::HasKey(const std::string &var_name)
const
{
    std::string key = ConfigParser::str_lower(ConfigParser::str_remove(var_name, ignore_chars));
    if(config_map.find(key) != config_map.end())
        return true;
    return false;
}

// list all the existing configuration keys
void ConfigObject::ListKeys()
const
{
    for(auto &it : config_map)
    {
        std::cout << it.first << std::endl;
    }
}

// get all the existing configuration keys
std::vector<std::string> ConfigObject::GetKeyList()
const
{
    std::vector<std::string> res;

    for(auto &it : config_map)
        res.push_back(it.first);

    return res;
}

// save current configuration into a file
void ConfigObject::SaveConfig(const std::string &path)
const
{
    std::string save_path;

    if(path.empty())
        save_path = config_path;
    if(save_path.empty())
        save_path = "latest.conf";

    std::ofstream save(save_path);
    for(auto &it : config_map)
    {
        save << it.first
             << " = "
             << it.second
             << std::endl;
    }
}

// get configuration value by its name/key
ConfigValue ConfigObject::GetConfigValue(const std::string &var_name)
const
{
    // convert to lower case and remove uninterested characters
    std::string key = ConfigParser::str_lower(ConfigParser::str_remove(var_name, ignore_chars));

    auto it = config_map.find(key);
    if(it == config_map.end())
        return __empty_value;
    else
        return form(it->second, replace_pair.first, replace_pair.second);
}

// set configuration value by its name/key
void ConfigObject::SetConfigValue(const std::string &var_name, const ConfigValue &c_value)
{
    // convert to lower case and remove uninterested characters
    std::string key = ConfigParser::str_lower(ConfigParser::str_remove(var_name, ignore_chars));

    config_map[key] = c_value;
}



//============================================================================//
// Protected Member Function                                                  //
//============================================================================//

// get configuration value from the map
// if no such config value exists, it will fill the default value in
ConfigValue ConfigObject::getDefConfig(const std::string &name,
                                       const ConfigValue &def_value,
                                       bool verbose)
{
    std::string key = ConfigParser::str_lower(ConfigParser::str_remove(name, ignore_chars));

    auto it = config_map.find(key);
    if(it == config_map.end())
    {
        if(def_value.IsEmpty())
            return __empty_value;

        if(verbose) {
            std::cout << name << " (key: " << key << ")"
                      << " not defined in configuration file, set to default value "
                      << def_value
                      << std::endl;
        }
        config_map[key] = def_value;
        return def_value;
    }

    return form(it->second, replace_pair.first, replace_pair.second);
}

// replace the contents inside replace_pair with the configuration value
ConfigValue ConfigObject::form(const std::string &input,
                               const std::string &op,
                               const std::string &cl)
const
{
    std::string result = input;
    int pos1, pos2;

    while(ConfigParser::find_pair(result, op, cl, pos1, pos2))
    {
        size_t beg = pos1 + op.size();
        size_t end = pos2 - cl.size();
        size_t size = end - beg + 1;

        std::string var = result.substr(beg, size);
        std::string val;

        if(pos1 > 0 && result.at(pos1 - 1) == '$') {
            pos1 --;
            val = std::getenv(var.c_str());
        } else {
            val = GetConfig<std::string>(var);
        }

        result.replace(pos1, pos2 - pos1 + 1, val);
    }

    return ConfigValue(std::move(result));
}



//============================================================================//
// Private Member Function                                                    //
//============================================================================//

// parse the control word and respond
void ConfigObject::parseControl(const std::string &word)
{
    if(ConfigParser::str_upper(word.substr(0, 7)) == "INCLUDE") {
        // need the most outer pair
        auto pairs = ConfigParser::find_pairs(word, "(", ")");
        // not find pair
        if(pairs.empty()) {
            std::cout << "Unsupported control word format: " << word << "."
                      << "Expected: INCLUDE(path)"
                      << std::endl;
            return;
        }
        int begin = pairs.back().first + 1;
        int length = pairs.back().second - begin;
        ReadConfigFile(form(word.substr(begin, length), replace_pair.first, replace_pair.second)._value);
    }
    else {
        std::cout << "Unsupported control word: " << word << std::endl;
    }
}

// parse the configuration term
void ConfigObject::parseTerm(std::string &&var_name, std::string &&var_value)
{
    // convert to lower case and remove uninterested characters
    std::string key = ConfigParser::str_lower(ConfigParser::str_remove(var_name, ignore_chars));

    if(key.back() == '+') {
        key.pop_back();
        auto it = config_map.find(key);
        if(it != config_map.end())
            it->second += var_value;
        else
            config_map[key] = var_value;
    } else {
        config_map[key] = var_value;
    }
}

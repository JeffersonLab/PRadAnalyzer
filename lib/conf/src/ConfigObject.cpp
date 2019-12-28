//============================================================================//
// A class based on the support from ConfigParser and ConfigValue             //
// It provides a simple way to read text file as configuration file, and read //
// or modify a parameter in the inherited class                               //
// The Configure() function should be overloaded according to specialized     //
// requirements, and be called after the parameters being configured          //
//                                                                            //
// Chao Peng                                                                  //
// 10/31/2016                                                                 //
//============================================================================//

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include "ConfigObject.h"


//============================================================================//
// Constructor, Destructor                                                    //
//============================================================================//

// constructor
ConfigObject::ConfigObject(const std::string &splitter, const std::string &ignore, const std::string &var_open,
        const std::string &var_close, bool case_ins)
: split_chars(splitter), ignore_chars(ignore), variable_pair({var_open, var_close}), case_insensitive(case_ins),
  __empty_value("")
{
    // place holder
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
    ConfigParser c_parser;
    c_parser.SetSplitters(split_chars); // self-defined splitters

    if (c_parser.ReadFile(path)) {
        // current directory
        parserProcess(c_parser, path);
        return true;
    } else {
        std::cerr << "Cannot open configuration file "
                  << "\"" << path << "\""
                  << std::endl;
        return false;
    }
}

// read the configuration string directly
void ConfigObject::ReadConfigString(const std::string &content, const std::string &path)
{
    ConfigParser c_parser;
    c_parser.SetSplitters(split_chars);

    c_parser.ReadBuffer(content.c_str());

    parserProcess(c_parser, path);
}

// continue parse the terms
void ConfigObject::parserProcess(ConfigParser &c_parser, const std::string &source)
{
    char abs_path[2048];
    realpath(source.c_str(), abs_path);
    std::string cur_dir = ConfigParser::decompose_path(abs_path).dir;

    // THIS_DIR needs special treatment as many files in different dirs may be loaded into one instance
    std::string dir_key = variable_pair.first + "THIS_DIR" + variable_pair.second;

    while (c_parser.ParseLine()) {
        // possible control words
        if (c_parser.NbofElements() == 1) {
            std::string control = c_parser.TakeFirst();
            size_t pos = control.find(dir_key);
            if(pos != std::string::npos) { control.replace(pos, dir_key.size(), cur_dir); }
            parseControl(control);
        // var_name and var_value
        } else if (c_parser.NbofElements() == 2) {
            std::string var_name, key, var_value;
            c_parser >> var_name >> var_value;
            size_t pos = var_value.find(dir_key);
            if(pos != std::string::npos) { var_value.replace(pos, dir_key.size(), cur_dir); }
            parseTerm(std::move(var_name), std::move(var_value));
        // unsupported format
        } else {
            std::cout << "Warning: Unsupported format at line "
                      << c_parser.LineNumber() << " from " << source << "\n"
                      << "\"" << c_parser.CurrentLine() << "\""
                      << std::endl;
        }
    }
}

// check if a certain term is configured
bool ConfigObject::HasKey(const std::string &var_name)
const
{
    auto key = formKey(var_name);

    if (config_map.find(key) != config_map.end())
        return true;
    return false;
}

// list all the existing configuration keys
void ConfigObject::ListKeys()
const
{
    for (auto &it : config_map) {
        std::cout << it.first << std::endl;
    }
}

// get all the existing configuration keys
std::vector<std::string> ConfigObject::GetKeyList()
const
{
    std::vector<std::string> res;

    for (auto &it : config_map) {
        res.push_back(it.first);
    }

    return res;
}

// save current configuration into a file
void ConfigObject::SaveConfig(const std::string &path)
const
{
    std::string save_path;

    if (path.empty()) 
        save_path = config_path;
    if(save_path.empty())
        save_path = "latest.conf";

    std::ofstream save(save_path);
    for (auto &it : config_map) {
        save << it.first
             << " " << split_chars.front() << " "
             << it.second
             << std::endl;
    }
}

// get configuration value by its name/key
ConfigValue ConfigObject::Value(const std::string &var_name)
const
{
    // convert to lower case and remove uninterested characters
    auto key = formKey(var_name);

    auto it = config_map.find(key);
    if (it == config_map.end()) {
        return __empty_value;
    } else {
        ConfigValue result(it->second);
        reform(result._value, variable_pair.first, variable_pair.second);
        return result;
    }
}

// set configuration value by its name/key
void ConfigObject::SetConfigValue(const std::string &var_name, const ConfigValue &c_value)
{
    // convert to lower case and remove uninterested characters
    auto key = formKey(var_name);

    config_map[key] = c_value;
}


// get configuration value from the map
// if no such config value exists, it will fill the default value in
ConfigValue ConfigObject::Value(const std::string &var_name, const ConfigValue &def_value, bool verbose)
{
    auto key = formKey(var_name);

    auto it = config_map.find(key);
    if (it == config_map.end()) {
        if (def_value.IsEmpty())
            return __empty_value;

        if (verbose) {
            std::cout << var_name << " (key: " << key << ")"
                      << " not defined in configuration file, set to default value "
                      << def_value
                      << std::endl;
        }
        config_map[key] = def_value;
        return def_value;
    }

    ConfigValue result(it->second);
    reform(result._value, variable_pair.first, variable_pair.second);
    return result;
}



//============================================================================//
// Protected Member Function                                                  //
//============================================================================//

// build the key
std::string ConfigObject::formKey(const std::string &rawKey)
const
{
    std::string key = ConfigParser::str_remove(rawKey, ignore_chars);
    if (case_insensitive) {
        key = ConfigParser::str_lower(key);
    }
    return key;
}

// replace the contents inside variable_pair with the configuration value
void ConfigObject::reform(std::string &input, const std::string &op, const std::string &cl)
const
{
    // loop until no pair found
    while (true) {
        auto rpair = ConfigParser::find_pair(input, op, cl);
        if (rpair.first != std::string::npos && rpair.second != std::string::npos) {
            // get content inside the bracket
            std::string var = input.substr(rpair.first + op.size(),
                                           rpair.second - op.size() - rpair.first);
            // deal with nested structure
            reform(var, op, cl);

            // replace content
            std::string val = HasKey(var) ? Value(var)._value : std::getenv(var.c_str());

            // replace variable with configuration value
            input.replace(rpair.first, rpair.second - rpair.first + cl.size(), val);
        } else {
            // no pair found any more
            return;
        }
    }
}




//============================================================================//
// Private Member Function                                                    //
//============================================================================//

// parse the control word and respond
void ConfigObject::parseControl(const std::string &word)
{
    if (ConfigParser::str_upper(word.substr(0, 7)) == "INCLUDE") {
        // need the most outer pair
        auto p = ConfigParser::find_pair(word, "(", ")");
        // not find pair
        if (p.first == std::string::npos || p.second == std::string::npos) {
            std::cout << "Unsupported control word format: " << word << "."
                      << "Expected: INCLUDE(path)"
                      << std::endl;
            return;
        }

        int begin = p.first + 1;
        int length = p.second - begin;

        std::string new_path = word.substr(begin, length);
        reform(new_path, variable_pair.first, variable_pair.second);

        ReadConfigFile(new_path);
    }
    else {
        std::cout << "Unsupported control word: " << word << std::endl;
    }
}

// parse the configuration term
void ConfigObject::parseTerm(std::string &&var_name, std::string &&var_value)
{
    // convert to lower case and remove uninterested characters
    auto key = formKey(var_name);

    if (key.back() == '+') {
        key.pop_back();
        auto it = config_map.find(key);
        if (it != config_map.end())
            it->second += var_value;
        else
            config_map[key] = var_value;
    } else {
        config_map[key] = var_value;
    }
}

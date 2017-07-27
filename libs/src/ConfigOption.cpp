//============================================================================//
// A class to simplify the argument parsing procedure                         //
//                                                                            //
// Chao Peng                                                                  //
// 07/26/2017                                                                 //
//============================================================================//

#include "ConfigOption.h"
#include "ConfigParser.h"
#include <iostream>
#include <cstring>



ConfigOption::ConfigOption()
{
    // place holder
}

ConfigOption::~ConfigOption()
{
    // place holder
}

void ConfigOption::AddOpt(OptType type, char s_term)
{
    AddOpt(type, s_term, s_term);
}

void ConfigOption::AddLongOpt(OptType type, const char *l_term)
{
    AddLongOpt(type, l_term, *l_term);
}

void ConfigOption::AddOpts(OptType type, char s_term, const char *l_term)
{
    AddOpts(type, s_term, l_term, s_term);
}

void ConfigOption::AddOpt(OptType type, char s_term, char mark)
{
    if(s_opt_map.find(s_term) != s_opt_map.end()) {
        std::cerr << "ConfigOption: Short option " << s_term
                  << " has already been added, abort AddOpt."
                  << std::endl;
        return;
    }

    s_opt_map[s_term] = Opt(mark, type);
}

void ConfigOption::AddLongOpt(OptType type, const char *l_term, char mark)
{
    if(l_opt_map.find(l_term) != l_opt_map.end()) {
        std::cerr << "ConfigOption: Long option " << l_term
                  << " has already been added, abort AddOpt."
                  << std::endl;
        return;
    }

    l_opt_map[l_term] = Opt(mark, type);
}

void ConfigOption::AddOpts(OptType type, char s_term, const char *l_term, char mark)
{
    AddOpt(type, s_term, mark);
    AddLongOpt(type, l_term, mark);
}

bool ConfigOption::parseLongOpt(const char *arg)
{
    auto vars = ConfigParser::split(arg, strlen(arg), "=");

    std::string key = std::move(vars.front());
    vars.pop_front();
    auto it = l_opt_map.find(key);

    if(it == l_opt_map.end()) {
        std::cerr << "Unknown option --" << key << std::endl;
        return false;
    }

    switch(it->second.type)
    {
    case arg_require:
        if(vars.empty()) {
            std::cerr << "Lack of argument for option --" << key << std::endl;
            return false;
        }
        opt_pack.emplace_back(it->second.mark, it->second.type, std::move(vars.front()));
        break;
    case arg_none:
        opt_pack.emplace_back(it->second.mark, it->second.type);
        break;
    case help_message:
        return false;
    }

    return true;
}

bool ConfigOption::parseShortOpt(char key, int argc, char *argv[], int &idx)
{
    auto it = s_opt_map.find(key);
    if(it == s_opt_map.end()) {
        std::cerr << "Unknown option -" << key << std::endl;
        return false;
    }

    switch(it->second.type)
    {
    case arg_require:
        if((idx + 1) >= argc || *argv[idx + 1] == '-') {
            std::cerr << "Lack of argument for option -" << key << std::endl;
            return false;
        }
        opt_pack.emplace_back(it->second.mark, it->second.type, ConfigValue(argv[++idx]));
        break;
    case arg_none:
        opt_pack.emplace_back(it->second.mark, it->second.type);
        break;
    case help_message:
        return false;
    }

    return true;
}

bool ConfigOption::ParseArgs(int argc, char *argv[])
{
    bool success = true;
    arg_pack.clear();
    opt_pack.clear();

    for(int i = 1; i < argc; ++i)
    {
        char* ptr = argv[i];
        // option
        if(*ptr == '-') {
            // long option
            if(*(ptr + 1) == '-') {
                success &= parseLongOpt(ptr + 2);
            // short option
            } else {
                success &= parseShortOpt(*(ptr + 1), argc, argv, i);
            }
        // arguments
        } else {
            arg_pack.emplace_back(ptr);
        }
    }

    return success;
}

void ConfigOption::SetDesc(const char *desc)
{
    base_desc = desc;
}

void ConfigOption::SetDesc(char mark, const char *desc)
{
    std::string option;
    bool found_mark = false;
    for(auto &it : s_opt_map)
    {
        if(it.second.mark != mark) continue;
        option = " - ";
        option.back() = it.first;
        if(it.second.type == arg_require)
            option += " <arg>";
        option += ",";
        found_mark = true;
    }

    for(auto &it : l_opt_map)
    {
        if(it.second.mark != mark) continue;
        option += " --" + it.first;
        if(it.second.type == arg_require)
            option += "=<arg>";
        option += ",";
        found_mark = true;
    }

    if(!found_mark) {
        std::cerr << "Config Option: Cannot find any option with mark <"
                  << mark << ">, abort adding description." << std::endl;
        return;
    }

    option.back() = ':';
    option += " ";
    option += desc;
    option_desc.emplace_back(option);
}

std::string ConfigOption::GetInstruction()
{
    std::string res = base_desc + "\noptions:\n";
    for(auto &option : option_desc)
    {
        res += "\t" + option + "\n";
    }
    return res;
}


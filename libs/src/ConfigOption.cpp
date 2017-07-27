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

void ConfigOption::AddOpt(char s_term, OptType type)
{
    AddOpt(s_term, type, s_term);
}

void ConfigOption::AddOpt(char s_term, OptType type, char mark)
{
    if(s_opt_map.find(s_term) != s_opt_map.end()) {
        std::cerr << "ConfigOption: Short option " << s_term
                  << " has already been added, abort AddOpt."
                  << std::endl;
        return;
    }

    s_opt_map[s_term] = Opt(mark, type);
}

void ConfigOption::AddOpt(const char *l_term, OptType type, char mark)
{
    if(l_opt_map.find(l_term) != l_opt_map.end()) {
        std::cerr << "ConfigOption: Long option " << l_term
                  << " has already been added, abort AddOpt."
                  << std::endl;
        return;
    }

    l_opt_map[l_term] = Opt(mark, type);
}

void ConfigOption::AddOpt(char s_term, const char *l_term, OptType type)
{
    AddOpt(s_term, type, s_term);
    AddOpt(l_term, type, s_term);
}

void ConfigOption::AddOpt(char s_term, const char *l_term, OptType type, char mark)
{
    AddOpt(s_term, type, mark);
    AddOpt(l_term, type, mark);
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

    if(it->second.type == arg_require) {
        if(vars.empty()) {
            std::cerr << "Lack of argument for option --" << key << std::endl;
            return false;
        }
        opt_pack.emplace_back(it->second.mark, it->second.type, std::move(vars.front()));
    } else {
        opt_pack.emplace_back(it->second.mark, it->second.type);
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

    if(it->second.type == arg_require) {
        if((idx + 1) >= argc || *argv[idx + 1] == '-') {
            std::cerr << "Lack of argument for option -" << key << std::endl;
            return false;
        }

        opt_pack.emplace_back(it->second.mark, it->second.type, ConfigValue(argv[++idx]));
    } else {
        opt_pack.emplace_back(it->second.mark, it->second.type);
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


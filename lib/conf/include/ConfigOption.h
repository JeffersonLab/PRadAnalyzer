#ifndef CONFIG_OPTION_H
#define CONFIG_OPTION_H

#include "ConfigValue.h"
#include <vector>
#include <string>
#include <unordered_map>

class ConfigOption
{
public:
    enum OptType : int
    {
        arg_none = 0,
        arg_require,
        help_message,
    };

    struct Opt
    {
        char mark;
        OptType type;
        ConfigValue var;

        Opt() : mark(-1), type(arg_none) {}
        Opt(char m, OptType t) : mark(m), type(t) {}
        Opt(char m, OptType t, ConfigValue v) : mark(m), type(t), var(v) {}
    };

public:
    ConfigOption();
    virtual ~ConfigOption();

    void AddOpt(OptType type, char s_term);
    void AddOpt(OptType type, char s_term, char mark);
    void AddLongOpt(OptType type, const char *l_term);
    void AddLongOpt(OptType type, const char *l_term, char mark);
    void AddOpts(OptType type, char s_term, const char *l_term);
    void AddOpts(OptType type, char s_term, const char *l_term, char mark);
    void SetDesc(const char *desc);
    void SetDesc(char mark, const char *desc);
    std::string GetInstruction();
    bool ParseArgs(int argc, char *argv[]);
    size_t NbofArgs() const {return arg_pack.size();}
    size_t NbofOpts() const {return opt_pack.size();}
    const ConfigValue &GetArgument(size_t i) {return arg_pack.at(i);}
    const std::vector<ConfigValue> &GetArguments() {return arg_pack;}
    const std::vector<Opt> &GetOptions() {return opt_pack;}

private:
    bool parseLongOpt(const char *arg);
    bool parseShortOpt(char key, int argc, char *argv[], int &idx);

private:
    std::unordered_map<char, Opt> s_opt_map;
    std::unordered_map<std::string, Opt> l_opt_map;
    std::vector<Opt> opt_pack;
    std::vector<ConfigValue> arg_pack;
    std::string base_desc;
    std::vector<std::string> option_desc;
};


#endif // CONFIG_OPTION_H

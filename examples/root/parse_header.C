#include "ConfigParser.h"
#include <iostream>
#include <string>
#include <set>


void parse_header(const char *cpath = "tstyle_header.conf")
{
    ConfigParser parser;
    parser.SetSplitters(" \t");
    parser.SetQuotePair("(", ")");
    parser.ReadFile(cpath);

    std::set<std::string> enums = { "EPaperSize" };
    auto char_type = [](const std::string &type) { return (type == "const char*") || (type == "Option_t"); };

    std::string ret, func, input;
    std::vector<std::string> normals, specials, defaults;
    struct ArgAttr { std::string type, def; };
    while (parser.ParseLine()) {
        parser >> ret >> func >> input;
        auto args = ConfigParser::split(input, ",");
        std::vector<ArgAttr> attr;
        for (auto &arg : args) {
            std::string type, def;
            auto vals = ConfigParser::split(arg, "=");
            if (vals.size() > 1) { def = vals[1]; }

            vals = ConfigParser::split(arg, " ");
            type = vals[0];
            if (type == "const") {
                type += " " + vals[1];
                if (vals[2][0] == '*') {
                    type += '*';
                }
            }
            if (enums.find(type) != enums.end()) {
                type = "Int_t";
            }
            attr.push_back({type, def});
        }

        // default values
        std::string dline = func.substr(3) + " = ";
        for (size_t i = 0; i < attr.size(); ++i) {
            if (attr[i].def.empty()) { continue; }
            if (i > 0) { dline += ", "; }
            dline += attr[i].def;
        }
        defaults.push_back(dline);


        // special cases
        if (attr.size() == 1) {
            if ( char_type(attr[0].type) ) {
                normals.push_back("SET_ATTR_CSTR(style, " + func.substr(3) + ", obj);");
            } else {
                normals.push_back("SET_ATTR(style, " + func.substr(3) + ", " + attr[0].type + ", obj);");
            }
        } else if ((attr.size() == 2) && char_type(attr[1].type) && !char_type(attr[0].type)) {
            normals.push_back("SET_ATTR_CSTR2(style, " + func.substr(3) + ", " + attr[0].type + ", obj);");
        } else {
            std::cout << parser.CurrentLine() << std::endl;
            for (auto &a : attr) {
                std::cout << a.type << " = " << a.def << std::endl;
            }
        }
    }

    for (auto &line : normals) {
        std::cout << line << std::endl;
    }
    return 0;
}

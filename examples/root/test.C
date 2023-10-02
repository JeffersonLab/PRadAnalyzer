#include "ConfigParser.h"
#include "ConfigObject.h"
#include "canalib.h"

using namespace std;

// helper function
inline bool compare_str(const char *buf1, const char *buf2, size_t n)
{
    for (size_t i = 0; i < n; ++i) {
        if (buf1[i] != buf2[i]) {
            return false;
        }
    }
    return true;
}

// helper function
inline std::string break_line(const std::string &buf, const std::string &delim, size_t &i, const std::string &glue)
{
    size_t beg = i;
    if (delim.empty()) {
        i = buf.size();
        return buf.substr(beg);
    }

    for (; i <= buf.size() - delim.size(); ++i) {
        if (compare_str(&buf[i], delim.c_str(), delim.size())) {
            std::string res = buf.substr(beg, i - beg);
            i += delim.size();
            if (glue.size() && (res.size() > glue.size())) {
                size_t pos = res.size() - glue.size();
                if (compare_str(&res[pos], glue.c_str(), glue.size())) {
                    return res.substr(0, pos) + break_line(buf, delim, i, glue);
                }
            }
            return res;
        }
    }
    i = buf.size();
    return buf.substr(beg);
}

void test_find_pair()
{
    string tstr = "{|}asvs{asdasdss}asdasd}";
    string op = "{";
    string cl = "}";
    auto p = ConfigParser::find_pair(tstr, op, cl);

    if(p.first == string::npos || p.second == string::npos) {
        cout << "not found" << endl;
    } else {
        cout << tstr.substr(p.first + op.size(), p.second - op.size() - p.first) << endl;
    }
}

void test_reform()
{
    ConfigObject conf_obj;

    conf_obj.ReadConfigString("Var1 = 2 \n"
             "Var2 = 1{Var1}3456 \n"
             "Var3 = {Var{Var1}}789 \n"
             "Var4 = {Var{Var1}}{Var1}{Var1} \n"
             "Path = variable_test.conf \n"
             "INCLUDE({Path})");

    string var1 = conf_obj.Value<string>("Var1");
    string var2 = conf_obj.Value<string>("Var2");
    string var3 = conf_obj.Value<string>("Var3");
    string var4 = conf_obj.Value<string>("Var4");
    string var5 = conf_obj.Value<string>("Var5");

    cout << var1 << endl
         << var2 << endl
         << var3 << endl
         << var4 << endl
         << var5 << endl;
}

void test_block_read(const string &path = "block_test.conf")
{
    // read in file
    string buffer = ConfigParser::file_to_string(path);

    // remove comments
    ConfigParser::comment_between(buffer, "/*", "*/");
    ConfigParser::comment_line(buffer, "//", "\n");
    ConfigParser::comment_line(buffer, "#", "\n");

    // get content blocks
    auto text = ConfigParser::break_into_blocks(buffer, "{", "}");

    cout << "residual: " << text.residual << endl;

    for(auto &b : text.blocks)
    {
        cout << "label: " << b.label << endl;
        cout << "content: " << b.content << endl;
    }
}


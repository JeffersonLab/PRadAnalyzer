//============================================================================//
// A simple parser class to read text file and deal with strings              //
//                                                                            //
// Chao Peng                                                                  //
// 06/07/2016                                                                 //
//============================================================================//

#include "ConfigParser.h"
#include <cstring>
#include <climits>
#include <algorithm>

using namespace std;



//============================================================================//
// Constructors, Destructor, Assignment Operators                             //
//============================================================================//

// constructor, with format input
ConfigParser::ConfigParser(Format f)
: form(f), line_number(0)
{
    // place holder
}

// copy constructor, only copy format
ConfigParser::ConfigParser(const ConfigParser &that)
: form(that.form), line_number(0)
{
    // place holder
}

// move constructor, only move format
ConfigParser::ConfigParser(ConfigParser &&that)
: form(that.form), line_number(0)
{
    // place holder
}

// desctructor
ConfigParser::~ConfigParser()
{
    CloseFile();
}

// copy assignment operator
ConfigParser &ConfigParser::operator = (const ConfigParser &rhs)
{
    form = rhs.form;
    return *this;
}

// move assignment operator
ConfigParser &ConfigParser::operator = (ConfigParser &&rhs)
{
    form = rhs.form;
    return *this;
}

//============================================================================//
// Public Member Function                                                     //
//============================================================================//

// open a file for future parsing
bool ConfigParser::OpenFile(const string &path, size_t cap)
{
    Clear();

    infile.open(path);

    if(!infile.is_open())
        return false;

    buf.data.resize(cap);

    // success!
    return true;
}

// read the whole file into a buffer and break it into lines
bool ConfigParser::ReadFile(const string &path)
{
    Clear();

    infile.open(path);

    if(!infile.is_open())
        return false;

    infile.seekg(0, ios::end);
    buf.end = infile.tellg();
    buf.data.resize(buf.end);
    infile.seekg(0, ios::beg);

    infile.read(&buf.data[0], buf.end);
    infile.close();

    return true;
}

// close file
void ConfigParser::CloseFile()
{
    return infile.close();
}

// clear stored lines
void ConfigParser::Clear()
{
    buf.Reset();

    // reset line
    line_number = 0;
    cur_line.Reset();

    // close file
    CloseFile();
}

// read a buffer in
void ConfigParser::ReadBuffer(const char *buf_in)
{
    Clear();

    buf.end = strlen(buf_in);
    buf.data.resize(buf.end + 2);

    strncpy(&buf.data[0], buf_in, buf.end);
}

// parse a line from the file or buffer
// if the line is empty (all white spaces or comments), it will be skipped
// return false if reached the end
bool ConfigParser::ParseLine()
{
    elements.clear();

    while(elements.empty())
    {
        if(!getLine(cur_line))
            return false;

        // count the line number
        ++line_number;

        parseBuffer(cur_line);
    }

    return true;
}

// parse the whole file or buffer
// return false if nothing was found
bool ConfigParser::ParseAll()
{
    elements.clear();

    while(true)
    {
        if(!getLine(cur_line))
            return !elements.empty();

        // count the line number
        ++line_number;

        parseBuffer(cur_line);
    }
}

// parse an input string, split the string into elements
// the trail white spaces in the elements will be trimmed
int ConfigParser::ParseString(const string &line)
{
    deque<string> eles = split(line.c_str(), line.size(), form.split);

    int count = 0;
    for(auto &ele : eles)
    {
        string trim_ele = trim(ele, form.white);
        if(!trim_ele.empty()) {
            elements.emplace_back(move(trim_ele));
            count++;
        }
    }

    return count;
}

// check if the current elemtns number is in the range [num, num + optional]
// output a warning message if not
bool ConfigParser::CheckElements(int num, int optional)
{
    string num_str;

    if(optional > 0) {
        if((elements.size() >= (size_t)num) &&
           (elements.size() <= (size_t)(num + optional))) {
            return true;
        }

        num_str = to_string(num) + " - " + to_string(num + optional);

    } else if(optional == 0) {

        if(elements.size() == (size_t)num) {
            return true;
        }

        num_str = to_string(num);

    } else { // optional < 0
        if(elements.size() >= (size_t)num) {
            return true;
        }

        num_str = " >= " + to_string(num);
    }


    cout << "Config Parser Warning: Wrong format at line "
         << line_number
         << ", expecting " << num_str << " elements. "
         << endl
         << "\"" << cur_line.String() << "\""
         << endl;
    return false;
}


// take the first element
ConfigValue ConfigParser::TakeFirst()
{
    if(elements.empty()) {
        cout << "Config Parser Warning: Trying to take elements while there is "
             << "nothing, 0 value returned." << endl;
        return ConfigValue("0");
    }

    ConfigValue output(move(elements.front()));
    elements.pop_front();

    return output;
}



//============================================================================//
// Private Member Function                                                    //
//============================================================================//

// get buffer from the file or the input buffer
// return false if reached input end
bool ConfigParser::getBuffer()
{
    if(buf.begin < buf.end)
        return true;

    if(!infile.is_open() || infile.bad() || infile.eof())
        return false;

    infile.read(&buf.data[0], buf.data.size());

    buf.begin = 0;
    buf.end = infile.gcount();

    return true;
}

// trim white spaces
inline void trimbuf(const vector<char> &buf, size_t &begin, size_t &end, const string &w)
{
    while(begin < end)
    {
        if(w.find(buf[begin]) == string::npos)
            break;

        begin++;
    }

    while(end > begin)
    {
        if(w.find(buf[end - 1]) == string::npos)
            break;

        end--;
    }
}

inline bool compare(const char ch, const string &str, size_t &c1)
{
    if(str.empty())
        return false;

    c1 = (ch == str[c1]) ? (c1 + 1) : 0;

    return (c1 >= str.size());
}

inline bool rcompare(const ConfigParser::CharBuffer &buf, const string &str)
{
    if(str.empty() || buf.end <= buf.begin || buf.end - buf.begin < str.size())
        return false;

    for(size_t i = 1; i <= str.size(); ++i)
    {
        if(str[str.size() - i] != buf[buf.end - i])
            return false;
    }

    return true;
}

// a helper structure to check context status
struct TextStatus
{
    int val;
    size_t cmt1, cmt2, delim;

    TextStatus() : val(0), cmt1(0), cmt2(0), delim(0) {}
    inline void Set(int i) {val = i; cmt1 = 0; cmt2 = 0;}
};

// get a line from the file or buffer
// it deals with comments, white spaces
// return false if reached the end
bool ConfigParser::getLine(CharBuffer &line_buf, bool recursive)
{
    if(!recursive)
        line_buf.Reset();
    bool success = false;
    TextStatus stat;

    while(getBuffer())
    {
        success = true;

        while(buf.begin < buf.end)
        {
            auto &ch = buf[buf.begin++];
            switch(stat.val)
            {
            default:
            case 0:
                line_buf.Add(ch);
                // check if it is the end
                if(compare(ch, form.delim, stat.delim)) {
                    line_buf.end -= form.delim.size();
                    trimbuf(line_buf.data, line_buf.begin, line_buf.end, form.white);
                    // glue lines
                    if(rcompare(line_buf, form.glue)) {
                        line_buf.end -= form.glue.size();
                        return getLine(line_buf, true);
                    } else {
                        return success;
                    }
                } else if(compare(ch, form.cmtopen, stat.cmt1)) {
                    stat.Set(1);
                    line_buf.end -= form.cmtopen.size();
                } else if(compare(ch, form.cmtmark, stat.cmt2)) {
                    stat.Set(2);
                    line_buf.end -= form.cmtmark.size();
                }
                break;
            case 1:
                if(compare(ch, form.cmtclose, stat.cmt1)) {
                    stat.Set(0);
                }
                break;
            case 2:
                if(ch == '\n') {
                    stat.Set(0);
                    buf.begin -= 1;
                }
                break;
            }
        }
    }

    trimbuf(line_buf.data, line_buf.begin, line_buf.end, form.white);
    return success;
}

// parse an input char buffer, split the string into elements
// the trail white spaces in the elements will be trimmed
int ConfigParser::parseBuffer(const CharBuffer &line)
{
    if(line.begin >= line.end)
        return 0;

    size_t ele_begin = line.begin;
    int count = 0;

    // intended to visit i == line.end, so the rest of the string get parsed
    for(size_t i = 0; i <= line.end; ++i)
    {
        if(i == line.end || form.split.find(line[i]) != string::npos) {
            size_t ele_end = i;
            trimbuf(line.data, ele_begin, ele_end, form.white);
            if(ele_begin < ele_end) {
                elements.emplace_back(&line[ele_begin], ele_end - ele_begin);
                count++;
            }
            ele_begin = i + 1;
        }
    }

    return count;
}

//============================================================================//
// Public Static Function                                                     //
//============================================================================//

// comment out a string, remove chars from the comment mark to the line break
void ConfigParser::comment_line(string &str, const string &c, const string &b)
{
    // no need to continue
    if(str.empty() || c.empty() || b.empty())
        return;

    // loop until no marks found
    while(true)
    {
        size_t c_begin = str.find(c);
        if(c_begin != string::npos) {
            size_t c_end = str.find(b, c_begin + c.size());
            // found, comment out until the line break
            if(c_end != string::npos) {
                // do not remove line break
                str.erase(c_begin, c_end - c_begin);
            // not found, comment out until the end
            } else {
                str.erase(c_begin);
                // can stop now, since everything afterwards is removed
                return;
            }
        } else {
            // comment marks not found
            return;
        }
    }
}

// comment out between a pair of comment marks
// NOTICE: does not support nested structure of comment marks
void ConfigParser::comment_between(string &str, const string &open, const string &close)
{
    // no need to continue
    if(str.empty() || open.empty() || close.empty())
        return;

    while(true)
    {
        // find the openning comment mark
        size_t pos1 = str.find(open);
        if(pos1 != string::npos) {
            size_t pos2 = str.find(close, pos1 + open.size());
            // found pair
            if(pos2 != string::npos) {
                // remove everything between, including this pair
                str.erase(pos1, pos2 + close.size() - pos1);
            // comment pair not found
            } else {
                return;
            }
        } else {
            // comment pair not found
            return;
        }
    }
}

// trim all the characters defined as white space at both ends
string ConfigParser::trim(const string &str, const string &w)
{

    const auto strBegin = str.find_first_not_of(w);
    if (strBegin == string::npos)
        return ""; // no content

    const auto strEnd = str.find_last_not_of(w);

    const auto strRange = strEnd - strBegin + 1;
    return str.substr(strBegin, strRange);
}

// split a string into several pieces by all the characters defined as splitter
deque<string> ConfigParser::split(const string &str, const string &s)
{
    deque<string> eles;

    char *cstr = new char[str.length() + 1];

    strcpy(cstr, str.c_str());

    char *pch = strtok(cstr, s.c_str());

    while(pch != nullptr)
    {
        eles.emplace_back(pch);
        pch = strtok(nullptr, s.c_str());
    }

    delete[] cstr;

    return eles;
}

// split a char array into several pieces
deque<string> ConfigParser::split(const char* str, const size_t &len, const string &s)
{
    deque<string> eles;

    char *str_cpy = new char[len + 1];

    strncpy(str_cpy, str, len);
    // end of C string
    str_cpy[len] = '\0';

    char *pch = strtok(str_cpy, s.c_str());

    while(pch != nullptr)
    {
        eles.emplace_back(pch);
        pch = strtok(nullptr, s.c_str());
    }

    delete[] str_cpy;

    return eles;
}

// split a string and convert all parts to float numbers
vector<int> ConfigParser::stois(const string &str, const string &s, const string &w)
{
    vector<int> res;
    for(auto &val : split(str, s))
    {
        res.push_back(stoi(trim(val, w)));
    }
    return res;
}

// split a string and convert all parts to float numbers
vector<float> ConfigParser::stofs(const string &str, const string &s, const string &w)
{
    vector<float> res;
    for(auto &val : split(str, s))
    {
        res.push_back(stof(trim(val, w)));
    }
    return res;
}

// split a string and convert all parts to double numbers
vector<double> ConfigParser::stods(const string &str, const string &s, const string &w)
{
    vector<double> res;
    for(auto &val : split(str, s))
    {
        res.push_back(stod(trim(val, w)));
    }
    return res;
}

// get the split part at num
string ConfigParser::get_split_part(int num, const char *str, const char &s)
{
    // unavailable
    if(num < 0) return "";

    int beg = 0, cur = 0;
    while(str[cur] != '\0')
    {
        if(str[cur] == s) {
            // number macthed
            if(num-- == 0) {
                return string(&str[beg], cur - beg);
            // update segment
            } else {
                beg = cur + 1;
            }
        }
        ++cur;
    }

    // last element
    if(num == 0)
        return string(&str[beg], cur - beg);

    return "";
}

// check if the short string is the same with the first part of a long string
bool cstr_cmp_helper(const char *cmp, const char *str, int size)
{
    for(int i = 0; i < size; ++i)
    {
        if(cmp[i] != str[i] || cmp[i] == '\0' || str[i] == '\0')
            return false;
    }

    if(cmp[size] != '\0')
        return false;

    return true;
}

// split a long string and find if a short string is belong to its elements
int ConfigParser::get_part_count(const char *cmp, const char *str, const char &s)
{
    int cnt = 0, beg = 0, cur = 0;
    while(str[cur] != '\0')
    {
        if(str[cur] == s) {
            if(cstr_cmp_helper(cmp, &str[beg], cur - beg)) {
                return cnt;
            }

            ++cnt;
            beg = cur + 1;
        }
        ++cur;
    }

    if(cstr_cmp_helper(cmp, &str[beg], cur-beg))
        return cnt;

    return -1;
}

// find the integer in a string
int ConfigParser::find_integer(const string &str, const size_t &pos)
{
    vector<int> integers = find_integers(str);
    if(pos >= integers.size())
    {
        cerr << "Config Parser: Cannot find " << pos + 1 << " integers from "
             << "\"" << str << "\"."
             << endl;
        return 0;
    }

    return integers.at(pos);
}

// find all the integers in a string
vector<int> ConfigParser::find_integers(const string &str)
{
    vector<int> result;

    find_integer_helper(str, result);

    return result;
}

// helper function for finding a integer
void ConfigParser::find_integer_helper(const string &str, vector<int> &result)
{
   if(str.empty())
       return;

   int negative = 1;
   auto numBeg = str.find_first_of("-0123456789");
   if(numBeg == string::npos)
       return;

   // check negative sign
   string str2 = str.substr(numBeg);

   if(str2.at(0) == '-')
   {
       negative = -1;
       int num_check;

       do {
           str2.erase(0, 1);

           if(str2.empty())
               return;

           num_check = str2.at(0) - '0';
       } while (num_check > 9 || num_check < 0);
   }

   auto numEnd = str2.find_first_not_of("0123456789");
   if(numEnd == string::npos)
       numEnd = str2.size();

   int num = 0;
   size_t i = 0;

   for(; i < numEnd; ++i)
   {
       if( (num > INT_MAX/10) ||
           (num == INT_MAX/10 && ((str2.at(i) - '0') > (INT_MAX - num*10))) )
       {
           ++i;
           break;
       }

       num = num*10 + str2.at(i) - '0';
   }

   result.push_back(negative*num);
   find_integer_helper(str2.substr(i), result);
}

// return the lower case of this string
string ConfigParser::str_lower(const string &str)
{
    string res = str;
    for(auto &c : res)
    {
        c = tolower(c);
    }
    return res;
}

// return the upper case of this string
string ConfigParser::str_upper(const string &str)
{
    string res = str;
    for(auto &c : res)
    {
        c = toupper(c);
    }
    return res;
}

// remove characters in ignore list
string ConfigParser::str_remove(const string &str, const string &iignore)
{
    string res = str;

    for(auto &c : iignore)
    {
        res.erase(remove(res.begin(), res.end(), c), res.end());
    }
    return res;
}

// replace characters in the list with certain char
string ConfigParser::str_replace(const string &str, const string &list, const char &rc)
{
    if(list.empty())
        return str;

    string res = str;

    for(auto &c : res)
    {
        if(list.find(c) != string::npos)
            c = rc;
    }

    return res;
}

// compare two strings, can be case insensitive
bool ConfigParser::case_ins_equal(const string &str1, const string &str2)
{
    if(str1.size() != str2.size()) {
        return false;
    }

    for(auto c1 = str1.begin(), c2 = str2.begin(); c1 != str1.end(); ++c1, ++c2)
    {
        if(tolower(*c1) != tolower(*c2)) {
            return false;
        }
    }

    return true;
}

// find the first pair position in a string
// it will return the most outer pair if the first pair was in a nested structure
pair<size_t, size_t> ConfigParser::find_pair(const string &str,
                                             const string &open,
                                             const string &close,
                                             size_t pos)
{
    pair<size_t, size_t> res(string::npos, string::npos);

    if(open.empty() || close.empty() || str.size() <= pos)
        return res;

    res.first = str.find(open, pos);

    // pair not found
    if(res.first == string::npos) {
        return res;
    }

    int open_bracket = 1;
    size_t search_beg = res.first + open.size();

    // loop for nested structure
    while(open_bracket > 0)
    {
        size_t next_close = str.find(close, search_beg);

        // pair not found
        if(next_close == string::npos) {
            // change back to npos for the not-found indication
            res.first = string::npos;
            return res;
        }

        // check for nested structure
        size_t next_open = str.find(open, search_beg);

        // the comparison is based on the definition of string::npos
        // npos for not found is size_t = -1, which is the biggest size_t value
        // find another open before close
        if(next_open < next_close) {
            open_bracket++;
            search_beg = next_open + open.size();
        // else cases
        // 1. close mark found before open mark
        // 2. close mark found, open mark not
        // 3. close mark is the same as open mark, so the position is the same
        } else {
            open_bracket--;
            search_beg = next_close + close.size();
            res.second = next_close;
        }
    }

    return res;
}

// get file name and directory from a path
ConfigParser::PathInfo ConfigParser::decompose_path(const string &path)
{
    PathInfo res;
    if(path.empty()) return res;

    // find directory
    auto dir_pos = path.find_last_of("/");

    if(dir_pos != string::npos) {
        res.dir = path.substr(0, dir_pos);
        res.name = path.substr(dir_pos + 1);
    } else {
        res.name = path;
    }

    // find extension
    auto ext_pos = res.name.find_last_of(".");
    if(ext_pos != string::npos) {
        res.ext = res.name.substr(ext_pos + 1);
        res.name = res.name.substr(0, ext_pos);
    }

    return res;
}

// form the path
string ConfigParser::compose_path(const ConfigParser::PathInfo &path)
{
    string res(path.dir);
    res.reserve(path.dir.size() + path.name.size() + path.ext.size() + 2);

    if(!res.empty() && res.back() != '/')
        res += '/';

    res += path.name;

    if(!path.ext.empty())
        res += "." + path.ext;

    return res;
}

// form a path from given directory and file name, automatically add / if it is missing
string ConfigParser::form_path(const string &dir, const string &file)
{
    string file_path;
    file_path.reserve(dir.size() + file.size() + 1);

    file_path = dir;
    if(file_path.size() && file_path.back() != '/') file_path += "/";
    file_path += file;

    return file_path;
}

// read a file and return its content in a char string
string ConfigParser::file_to_string(const string &path)
{
    ifstream inf(path);

    if(!inf.is_open())
        return "";

    // read the whole file in
    string str;

    inf.seekg(0, ios::end);
    str.reserve(inf.tellg());
    inf.seekg(0, ios::beg);

    str.assign((istreambuf_iterator<char>(inf)), istreambuf_iterator<char>());
    inf.close();

    return str;
}

// break text file into several blocks in the format
// <label> <open_mark> <content> <close_mark>
// return extracted <residual> {<label> <content>}
ConfigParser::TextBlocks ConfigParser::break_into_blocks(const string &buf,
                                                         const string &open,
                                                         const string &close,
                                                         const string &seps)
{
    TextBlocks result;

    if(buf.empty() || open.empty() || close.empty())
        return result;

    size_t last_end = 0;
    // loop until no blocks found
    while(true)
    {
        // find the contents in block brackets
        auto p = find_pair(buf, open, close, last_end);

        // no pair found anymore
        if(p.first == string::npos || p.second == string::npos)
            break;

        // add content
        TextBlock block;
        block.content = trim(buf.substr(p.first + open.size(), p.second - p.first - open.size()), seps);

        // find label
        string head = buf.substr(last_end, p.first - last_end);
        if(head.empty()) {
            block.label = "";
        } else {
            // find end of label
            auto end = head.find_last_not_of(seps);
            if(end == string::npos) end = head.size() - 1;
            // find begin of label
            auto beg = head.find_last_of(seps, end);
            if(beg == string::npos) beg = 0;
            // add label
            block.label = trim(head.substr(beg, end - beg + 1), seps);
            // other content goes to residual
            result.residual += head.substr(0, beg);
        }
        // combine blocks
        result.blocks.emplace_back(move(block));
        last_end = p.second + close.size();
    }

    // trim
    result.residual = trim(result.residual, seps);

    return result;
}

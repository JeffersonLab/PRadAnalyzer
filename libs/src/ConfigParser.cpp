//============================================================================//
// A simple parser class to read text file                                    //
//                                                                            //
// Chao Peng                                                                  //
// 06/07/2016                                                                 //
//============================================================================//

#include "ConfigParser.h"
#include <streambuf>
#include <cstring>
#include <climits>
#include <algorithm>
#include <list>

using namespace std;



//============================================================================//
// Constructor, Destructor                                                    //
//============================================================================//
ConfigParser::ConfigParser(const string &s,
                           const string &w,
                           const vector<string> &c,
                           const string_pair &p,
                           const string &g)
: splitters(s), white_spaces(w), comment_marks(c), comment_pair(p), line_glues(g),
  line_number(0), in_comment_pair(false)
{
    // place holder
}

ConfigParser::~ConfigParser()
{
    // place holder
}

//============================================================================//
// Public Member Function                                                     //
//============================================================================//

// add comment mark
void ConfigParser::AddCommentMark(const string &c)
{
    auto it = find(comment_marks.begin(), comment_marks.end(), c);
    if(it == comment_marks.end())
        comment_marks.push_back(c);
}

// remove certain comment mark
void ConfigParser::RemoveCommentMark(const string &c)
{
    auto it = find(comment_marks.begin(), comment_marks.end(), c);
    if(it != comment_marks.end())
        comment_marks.erase(it);
}

// clear all comment marks
void ConfigParser::EraseCommentMarks()
{
    comment_marks.clear();
}

// open a file for future parsing
bool ConfigParser::OpenFile(const string &path)
{
    Clear();

    infile.open(path);

    return infile.is_open();
}

// read the whole file into a buffer and break it into lines
bool ConfigParser::ReadFile(const string &path)
{
    Clear();

    ifstream f(path);

    if(!f.is_open())
        return false;

    string buffer;

    f.seekg(0, ios::end);
    buffer.reserve(f.tellg());
    f.seekg(0, ios::beg);

    buffer.assign((istreambuf_iterator<char>(f)), istreambuf_iterator<char>());
    f.close();

    // remove comments and break the buffers into lines
    bufferProcess(buffer);
    return true;
}

// read a buffer and break it into lines
void ConfigParser::ReadBuffer(const char *buf)
{
    Clear();

    string buffer = buf;

    // remove comments and break the buffers into lines
    bufferProcess(buffer);
}

// clear stored lines
void ConfigParser::Clear()
{
    line_number = 0;
    infile.close();
    lines.clear();
}

// close file
void ConfigParser::CloseFile()
{
    infile.close();
}

// take a line
string ConfigParser::TakeLine()
{
    if(lines.size()) {
        string out = lines.front();
        lines.pop_front();
        return out;
    }

    return "";
}

// parse a line, take the line either from the opening file (if it is opened)
// or from the stored lines
// return false if empty
bool ConfigParser::ParseLine()
{
    elements.clear();

    if(infile.is_open())
        return parseFile();
    else
        return parseBuffer();
}

// parse the whole file all buffer
bool ConfigParser::ParseAll()
{
    elements.clear();
    if(infile.is_open()) {
        while(parseFile()) {;}
    } else {
        while(parseBuffer()) {;}
    }

    return !(elements.empty());
}

// parse a input string, trim the white space and split the string into elements
// by defined splitters, the elements will also be trimmed
// NOTICE: comments between a pair of marks are erased before input
int ConfigParser::ParseString(const string &line)
{
    deque<string> eles = split(line.c_str(), getCommentPoint(line), splitters);

    int count = 0;
    for(auto &ele : eles)
    {
        string trim_ele = trim(ele, white_spaces);
        if(trim_ele.size()) {
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
         << "\"" << current_line << "\""
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


// break text file into several blocks in the format
// <label> <open_mark> <content> <close_mark>
// return extracted <content>
vector<string> ConfigParser::ReadFileInBlocks(const string &path,
                                              const string &open,
                                              const string &close,
                                              const string &label,
                                              bool case_sensitive)
const
{
    vector<string> result;

    if(open.empty() || close.empty())
        return result;

    ifstream inf(path);

    if(!inf.is_open()) {
        return result;
    }

    // read the whole file in
    string buf;

    inf.seekg(0, ios::end);
    buf.reserve(inf.tellg());
    inf.seekg(0, ios::beg);

    buf.assign((istreambuf_iterator<char>(inf)), istreambuf_iterator<char>());
    inf.close();

    // remove comments
    comment_between(buf, comment_pair.first, comment_pair.second);
    for(auto &mark : comment_marks)
    {
        comment_out(buf, mark, "\n");
    }

    // find the contents in block brackets
    auto pairs = find_pairs(buf, open, close);

    int last_end = 0;
    // we are not parse it line by line, so add the line break into white spaces
    string whites = white_spaces + "\n";

    // lambda to check character
    auto is_white = [](const char &c, const string &ws)
                    {
                        for(auto &w : ws)
                        {
                            if(c == w)
                                return true;
                        }
                        return false;
                    };

    for(auto &p : pairs)
    {
        bool find_content = false;

        // no label specified
        if(label.empty()) {
            if(p.first > (int)last_end)
                find_content = true;
        // find the label
        } else {
            int beg = last_end, end = p.first - 1;
            bool find_end = false;
            for(int i = end; i >= beg; --i)
            {
                // find last not of white spaces
                if(!find_end) {
                    if(is_white(buf.at(i), whites)) {
                        end--;
                    } else {
                        find_end = true;
                    }
                } else {
                    if(is_white(buf.at(i), whites))
                        beg = i + 1;
                }
            }

            string lb = (end > beg) ? buf.substr(beg, end - beg + 1) : "";
            find_content = (case_sensitive) ? (lb == label) : strcmp_case_insensitive(lb, label);
        }

        // add contents into result container
        if(find_content) {
            size_t beg = p.first + open.size();
            size_t size = p.second - p.first - open.size();
            result.push_back(buf.substr(beg, size));
            last_end = p.second + close.size();
        }
    }

    return result;
}


//============================================================================//
// Private Member Function                                                    //
//============================================================================//

// process the read-in buffer
void ConfigParser::bufferProcess(string &buffer)
{
    if(comment_pair.first.size() && comment_pair.second.size()) {
        // comment by pair marks
        comment_between(buffer, comment_pair.first, comment_pair.second);
    }

    string line;
    for(auto &c : buffer)
    {
        if(c != '\n') {
            line += c;
        } else {
            lines.push_back(line);
            line.clear();
        }
    }
}

// take a line from opened file and parse it
// comment with comment pair will be removed here
bool ConfigParser::parseFile()
{
    int count = 0;
    // comment pair needs to be taken care here
    while(!count)
    {
        if(infile.bad() || infile.eof())
            return false;

        string parse_string;
        getLineFromFile(parse_string);

        if(parse_string.empty())
            continue;

        // no need to take care comment pair
        if(comment_pair.first.empty() || comment_pair.second.empty()) {
            count = ParseString(parse_string);
        } else {
            // if we had comment pair opened
            if(in_comment_pair) {
                // see if there is an end to the comment pair
                auto c_end = parse_string.find(comment_pair.second);
                // end of a comment pair
                if(c_end != string::npos) {
                    count = ParseString(parse_string.substr(c_end + comment_pair.second.size()));
                    in_comment_pair = false;
                }
            // if no previous comment pair opened
            } else {
                // remove complete comment pair in one line
                comment_between(parse_string, comment_pair.first, comment_pair.second);
                // see if there is any comment pair opening
                auto c_beg = parse_string.find(comment_pair.first);
                // find comment pair openning
                if(c_beg != string::npos) {
                    count = ParseString(parse_string.substr(0, c_beg));
                    in_comment_pair = true;
                } else {
                    // no special treatment
                    count = ParseString(parse_string);
                }
            }
        }
    }

    return true;
}

// get a line from ifstream, it takes care of empty line and concatenated lines
inline void ConfigParser::getLineFromFile(string &to_be_parsed)
{
    bool continuation;

    do {
        continuation = false;
        // end or error reached
        if(!getline(infile, current_line))
            break;

        // count the line number
        ++line_number;

        // trim white spaces at both ends
        current_line = trim(current_line, white_spaces);

        // should continue to read next line
        if(current_line.empty()) {
            continuation = true;
        } else {
            for(auto &c : line_glues)
            {
                if(current_line.back() == c) {
                    current_line.pop_back();
                    continuation = true;
                    break;
                }
            }
        }

        // add current line to be parsed
        to_be_parsed.append(current_line);

    } while(continuation);
}

// take a line from the stored lines buffer and parse it
// comment by comment pair has been already removed before going into the lines
bool ConfigParser::parseBuffer()
{
    // comment pair has been taken care before breaking into lines,
    // so it's simple here
    int count = 0;
    while(!count)
    {
        if(lines.empty())
            return false;

        string parse_string;
        getLineFromBuffer(parse_string);
        if(parse_string.empty())
            continue;

        count = ParseString(parse_string);
    }

    return true; // parsed a line
}

inline void ConfigParser::getLineFromBuffer(string &to_be_parsed)
{
    bool continuation;

    do {
        continuation = false;

        // end reached
        if(lines.empty())
            break;

        // take line
        current_line = move(lines.front());
        lines.pop_front();

        // count the line number
        ++line_number;

        // trim white spaces at both ends
        current_line = trim(current_line, white_spaces);

        // should continue to read next line
        if(current_line.empty()) {
            continuation = true;
        } else {
            for(auto &c : line_glues)
            {
                if(current_line.back() == c) {
                    current_line.pop_back();
                    continuation = true;
                    break;
                }
            }
        }

        // add current line to be parsed
        to_be_parsed.append(current_line);

    } while(continuation);
}

// comment out the characters with certain mark
size_t ConfigParser::getCommentPoint(const string &str)
{
    size_t res = str.size();

    for(auto &mark : comment_marks)
    {
        const auto c_begin = str.find(mark);
        if(c_begin != string::npos && c_begin < res) {
            res = c_begin;
        }
    }
    return res;
}

//============================================================================//
// Public Static Function                                                     //
//============================================================================//

// comment out a string, remove chars from the comment mark to the line break
void ConfigParser::comment_out(string &str, const string &c, const string &b)
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
deque<string> ConfigParser::split(const char* str, const size_t &size, const string &s)
{
    deque<string> eles;

    char *str_cpy = new char[size + 1];

    // end of C string
    str_cpy[size] = '\0';

    strncpy(str_cpy, str, size);

    char *pch = strtok(str_cpy, s.c_str());

    while(pch != nullptr)
    {
        eles.emplace_back(pch);
        pch = strtok(nullptr, s.c_str());
    }

    delete[] str_cpy;

    return eles;
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
bool ConfigParser::strcmp_case_insensitive(const string &str1, const string &str2)
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
// it will return the most inner pair if the first pair was in a nested structure
bool ConfigParser::find_pair(const string &str,
                             const string &open, const string &close,
                             int &open_pos, int &close_pos)
{
    auto pairs = find_pairs(str, open ,close);
    if(pairs.empty())
        return false;

    open_pos = pairs.front().first;
    close_pos = pairs.front().second;
    return true;
}

// find all pair positions in a string
vector<pair<int, int>> ConfigParser::find_pairs(const string &str,
                                                const string &op,
                                                const string &cl)
{
    vector<pair<int, int>> result;

    if(op == cl) {
        list<int> opens;
        int end_pos = str.size() - op.size() + 1;
        for(int i = 0; i < end_pos; ++i)
        {
            if(str.substr(i, op.size()) == op) {
                opens.push_back(i);
                i += op.size() - 1;
            }
        }
        while(opens.size() > 1)
        {
            int beg = opens.front();
            opens.pop_front();
            int end = opens.front();
            opens.pop_front();
            result.emplace_back(beg, end);
        }
        return result;
    }

    if(op.empty() || cl.empty())
        return result;

    int mark_size = op.size() + cl.size();

    // there won't be any
    if(mark_size > (int)str.size())
        return result;

    list<int> opens;
    int end_pos = str.size() - op.size() + 1;
    for(int i = 0; i < end_pos; ++i)
    {
        if(str.substr(i, op.size()) == op) {
            opens.push_back(i);
            i += op.size() - 1;
        }
    }

    list<int> closes;
    end_pos = str.size() - cl.size() + 1;
    for(int i = 0; i < end_pos; ++i)
    {
        if(str.substr(i, cl.size()) == cl)
        {
            bool share_char = false;
            for(auto &idx : opens)
            {
                if(abs(i - idx) < (int)op.size())
                    share_char = true;
            }

            if(!share_char) {
                closes.push_back(i);
                i += cl.size() - 1;
            }
        }
    }

    for(auto cl_it = closes.begin(); cl_it != closes.end(); ++cl_it)
    {
        for(auto op_it = opens.rbegin(); op_it != opens.rend(); ++op_it)
        {
            if((*cl_it - *op_it) >= (int)op.size())
            {
                result.emplace_back(*op_it, *cl_it);
                opens.erase(next(op_it).base());
                break;
            }
        }
    }

    return result;
}

// get file name and directory from a path
ConfigParser::PathInfo ConfigParser::decompose_path(const string &path)
{
    // find directory and suffix
    auto dir_pos = path.find_last_of("/");
    auto suf_pos = path.find_first_of(".");

    PathInfo res;
    if(dir_pos == string::npos) {
        res.dir = ".";
        res.name = path.substr(0, suf_pos);
        if(suf_pos != string::npos)
            res.suffix = path.substr(suf_pos + 1);
    } else {
        res.dir = path.substr(0, dir_pos);
        res.name = path.substr(dir_pos + 1, suf_pos - dir_pos - 1);
        if(suf_pos != string::npos && suf_pos > dir_pos)
            res.suffix = path.substr(suf_pos + 1);
    }

    return res;
}

// form the path
string ConfigParser::compose_path(const ConfigParser::PathInfo &path)
{
    string res(path.dir);
    res.reserve(path.dir.size() + path.name.size() + path.suffix.size() + 2);

    if(!res.empty() && res.back() != '/')
        res += '/';

    res += path.name;

    if(!path.suffix.empty())
        res += "." + path.suffix;

    return res;
}

string ConfigParser::form_path(const string &dir, const string &file)
{
    string file_path;
    file_path.reserve(dir.size() + file.size() + 1);

    file_path = dir;
    if(file_path.size() && file_path.back() != '/') file_path += "/";
    file_path += file;

    return file_path;
}

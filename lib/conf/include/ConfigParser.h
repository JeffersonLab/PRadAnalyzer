#ifndef CONFIG_PARSER_H
#define CONFIG_PARSER_H

#include <string>
#include <vector>
#include <queue>
#include <deque>
#include <fstream>
#include "ConfigValue.h"


// a macro to auto generate enum2str and str2enum
// name mapping begins at bias and continuously increase, split by '|'
// an example:
// enum ABC {a = 3, b, c};
// ENUM_MAP(ABC, 3, "a|b|c")
// ABC2str(3) = "a"
// str2ABC("b") = 4
#define ENUM_MAP(type, bias, strings) \
    static std::string type ## 2str(int T) \
    { \
        return ConfigParser::get_split_part(T - (bias), strings, '|'); \
    }; \
    static type str2 ## type(const char *str) \
    { \
        return static_cast<type>(bias + ConfigParser::get_part_count(str, strings, '|')); \
    }

// config parser class
class ConfigParser
{
    typedef std::pair<std::string, std::string> string_pair;

public:
    struct Format
    {
        struct StrPair { std::string open, close; };
        std::string split;              // element splitters (chars)
        std::string white;              // white spaces (chars)
        std::string delim;              // line delimiter (string)
        std::string glue;               // line glue (string)
        std::string linecmt;            // comment marks (string), until the line breaker '\n'
        StrPair blockcmt;               // block commenting marks (open string, close string)
        StrPair quote;                  // quote marks (open string, close string)

        static Format Basic() { return {" \t,", " \t", "\n", "", "", {"", ""}, {"\"", "\""}}; }
        static Format BashLike() { return {" \t,", " \t", "\n", "\\", "#", {"\'", "\'"}, {"\"", "\""}}; }
        static Format CLike() { return {" \t,\n", " \t\n", ";", "", "//", {"/*", "*/"}, {"\"", "\""}}; }
    };

public:
    ConfigParser(Format f = Format::BashLike());

    // format related
    inline void SetFormat(Format &&f) { fmt = f; }
    inline void SetFormat(const Format &f) { fmt = f; }
    inline void SetSplitters(std::string s) { fmt.split = s; }
    inline void SetWhiteSpaces(std::string w) { fmt.white = w; }
    inline void SetCommentMark(std::string c) { fmt.linecmt = c; }
    inline void SetCommentPair(std::string o, std::string c) { fmt.blockcmt = {o, c}; }
    inline void SetLineGlues(std::string g) { fmt.glue = g; }
    inline void SetLineBreaks(std::string b) { fmt.delim = b; }
    inline void SetQuotePair(std::string o, std::string c) { fmt.quote = {o, c}; }

    const Format &GetFormat() const {return fmt;}

    // dealing with file/buffer
    bool ReadFile(const std::string &path);
    void ReadBuffer(const char*);
    void Clear();

    // parse line, return false if no more line to parse
    bool ParseLine();
    // parse the whole file or buffer, return false if no elements found
    bool ParseAll();
    // parse a string, trim and split it into elements
    int ParseString(const std::string &line);

    // get current parsing status
    bool CheckElements(int num, int optional = 0);
    int NbofElements() const { return elements.size(); }
    int LineNumber() const { return line_number; }
    std::string CurrentLine() const;

    // take the elements
    ConfigValue TakeFirst();

    template<typename T>
    T TakeFirst()
    {
        return TakeFirst().Convert<T>();
    }

    template<typename T>
    ConfigParser &operator >>(T &t)
    {
        t = (*this).TakeFirst().Convert<T>();
        return *this;
    }

    template<class BidirIt>
    int Take(BidirIt first, BidirIt last)
    {
        int count = 0;
        for(auto it = first; it != last; ++it, ++count)
        {
            if(elements.empty())
                break;

            *it = elements.front();
            elements.pop_front();
        }
        return count;
    }

    template<template<class, class> class Container>
    Container<ConfigValue, std::allocator<ConfigValue>> TakeAll()
    {
        Container<ConfigValue, std::allocator<ConfigValue>> res;
        while(elements.size())
        {
            res.emplace_back(std::move(elements.front()));
            elements.pop_front();
        }
        return res;
    }

    template<template<class, class> class Container, class T>
    Container<T, std::allocator<T>> TakeAll()
    {
        Container<T, std::allocator<T>> res;
        while(elements.size())
        {
            ConfigValue tmp(std::move(elements.front()));
            elements.pop_front();
            res.emplace_back(tmp.Convert<T>());
        }
        return res;
    }


private:
    // private functions
    void toLines(std::string buf);
    void parseBuffer();
    void retrieveLine();

private:
    // private members
    Format fmt;
    int line_number;
    std::string curr_line;
    std::vector<std::string> quotes;
    std::deque<std::string> lines, elements;


public:
    // static functions
    static void comment_line(std::string &str, const std::string &cmt, const std::string &brk);
    static void comment_line(std::string &str, const std::string &cmt, const std::string &brk,
                             const std::string &qmark);
    static void comment_between(std::string &str, const std::string &open, const std::string &close);
    static void comment_between(std::string &str, const std::string &open, const std::string &close,
                                const std::string &qmark);
    static void tokenize(std::string &str, std::vector<std::string> &contents, const std::string &token,
                         const std::string &open, const std::string &close);
    static inline void tokenize(std::string &str, std::vector<std::string> &contents, const std::string &token,
                         const std::string &qmark) { tokenize(str, contents, token, qmark, qmark); }
    static void untokenize(std::string &str, const std::vector<std::string> &contents, const std::string &token,
                           const std::string &open, const std::string &close);
    static inline void untokenize(std::string &str, const std::vector<std::string> &contents, const std::string &token,
                                  const std::string &qmark = "") { untokenize(str, contents, token, qmark, qmark); }
    static std::string trim(const std::string &str, const std::string &w);
    static std::deque<std::string> split(const std::string &str, const std::string &s);
    static std::deque<std::string> split(const char* str, const size_t &len, const std::string &s);
    static std::string get_split_part(int num, const char *str, const char &s);
    static int get_part_count(const char *cmp, const char *str, const char &s);
    static std::vector<int> stois(const std::string &str, const std::string &s, const std::string &w);
    static std::vector<float> stofs(const std::string &str, const std::string &s, const std::string &w);
    static std::vector<double> stods(const std::string &str, const std::string &s, const std::string &w);
    static std::string str_remove(const std::string &str, const std::string &ignore);
    static std::string str_replace(const std::string &str, const std::string &ignore, const char &rc = ' ');
    static std::string str_lower(const std::string &str);
    static std::string str_upper(const std::string &str);
    static std::pair<size_t, size_t> find_pair(const std::string &str,
                                               const std::string &open,
                                               const std::string &close,
                                               size_t pos = 0);
    static bool case_ins_equal(const std::string &str1, const std::string &str2);
    static int find_integer(const std::string &str, const size_t &pos = 0);
    static std::vector<int> find_integers(const std::string &str);
    static void find_integer_helper(const std::string &str, std::vector<int> &result);
    struct PathInfo { std::string dir, name, ext; };
    static PathInfo decompose_path(const std::string &path);
    static std::string compose_path(const PathInfo &path);
    static std::string form_path(const std::string &dir, const std::string &file);
    static std::string file_to_string(const std::string &path);
    // break text file into several blocks in the format
    // <label> <open_mark> <content> <close_mark>, this structure can be separated by sep characters
    // return extracted <residual> {<label> <content>} with white characters trimmed
    struct TextBlock {std::string label, content;};
    struct TextBlocks {std::string residual; std::vector<TextBlock> blocks;};
    static TextBlocks break_into_blocks(const std::string &buf,
                                        const std::string &open = "{",
                                        const std::string &close = "}",
                                        const std::string &seps = " \t\n");

};

#endif

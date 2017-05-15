#ifndef CONFIG_PARSER_H
#define CONFIG_PARSER_H

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <deque>
#include "ConfigValue.h"

// config parser class
class ConfigParser
{
    typedef std::pair<std::string, std::string> string_pair;

public:
    ConfigParser(const std::string &s = " ,\t",                     // splitters
                 const std::string &w = " \t",                      // white_space
                 const std::vector<std::string> &c = {"#", "//"},   // comment mark
                 const string_pair &p = std::make_pair("/*", "*/"), // comment pair
                 const std::string &g = "\\");                      // line glue chars
    virtual ~ConfigParser();

    // set members
    void SetSplitters(const std::string &s) {splitters = s;}
    void SetWhiteSpaces(const std::string &w) {white_spaces = w;}
    void SetCommentMarks(const std::vector<std::string> &c) {comment_marks = c;}
    void SetCommentPair(const std::string &o, const std::string &c)
    {comment_pair = std::make_pair(o, c);}
    void SetLineGlues(const std::string &g) {line_glues = g;}
    void AddCommentMark(const std::string &c);
    void RemoveCommentMark(const std::string &c);
    void EraseCommentMarks();

    // dealing with file/buffer
    bool OpenFile(const std::string &path);
    bool ReadFile(const std::string &path);
    void ReadBuffer(const char*);
    void CloseFile();
    void Clear();

    // parse line, return false if no more line to parse
    bool ParseLine();
    // parse the whole file or buffer, return false if no elements found
    bool ParseAll();
    // parse a string, trim and split it into elements
    int ParseString(const std::string &line);

    // get current parsing status
    bool CheckElements(int num, int optional = 0);
    int NbofElements() const {return elements.size();}
    int NbofLines() const {return lines.size();}
    int LineNumber() const {return line_number;}
    const std::string &CurrentLine() const {return current_line;}

    // take the lines/elements
    std::string TakeLine();
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

    // get members
    const std::string &GetSplitters() const {return splitters;}
    const std::string &GetWhiteSpaces() const {return white_spaces;}
    const std::vector<std::string> &GetCommentMarks() const {return comment_marks;}
    const string_pair &GetCommentPair() const {return comment_pair;}
    const std::string &GetLineGlues() const {return line_glues;}


private:
    // private functions
    void bufferProcess(std::string &buffer);
    bool parseFile();
    bool parseBuffer();
    size_t getCommentPoint(const std::string &str);
    void getLineFromFile(std::string &to_be_parsed);
    void getLineFromBuffer(std::string &to_be_parsed);

private:
    // private members
    std::string splitters;
    std::string white_spaces;
    std::vector<std::string> comment_marks;
    string_pair comment_pair;
    std::string line_glues;
    std::deque<std::string> lines;
    std::deque<std::string> elements;
    std::string current_line;
    int line_number;
    bool in_comment_pair;
    std::ifstream infile;

public:
    // static functions
    static void comment_line(std::string &str, const std::string &cmt, const std::string &brk);
    static void comment_between(std::string &str, const std::string &open, const std::string &close);
    static std::string trim(const std::string &str, const std::string &w);
    static std::deque<std::string> split(const std::string &str, const std::string &s);
    static std::deque<std::string> split(const char* str, const size_t &size, const std::string &s);
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
    struct PathInfo { std::string dir, name, suffix; };
    static PathInfo decompose_path(const std::string &path);
    static std::string compose_path(const PathInfo &path);
    static std::string form_path(const std::string &dir, const std::string &file);
    static std::string file_to_string(const std::string &path);
    // break text file into several blocks in the format
    // <label> <open_mark> <content> <close_mark>
    // return extracted <label> <content>
    struct TextBlock {std::string label, content;};
    static std::vector<TextBlock> break_into_blocks(const std::string &buf,
                                                    const std::string &open = "{",
                                                    const std::string &close = "}",
                                                    const std::string &sep = " \t\n");

};

#endif

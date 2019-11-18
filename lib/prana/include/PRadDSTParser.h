#ifndef PRAD_DST_PARSER_H
#define PRAD_DST_PARSER_H

#include <vector>
#include <fstream>
#include <string>
#include "PRadException.h"
#include "PRadEventStruct.h"



class PRadDSTParser
{
public:
    enum HeaderType : uint16_t
    {
        // headers
        FileHeader = 0xa1b1,
        EventHeader = 0xa2b2,
        MapHeader = 0xa3b3,
    };

    enum class Type : uint16_t
    {
        // event types
        event = 0,
        epics,
        max_type,
    };

    struct Header
    {
        uint16_t htype, etype;
        uint32_t length;

        Header() : htype(0), etype(0), length(0) {}
        Header(uint16_t h, uint16_t e, uint32_t l = 0)
        : htype(h), etype(e), length(l) {}
        Header(HeaderType h, Type e, uint32_t l = 0)
        : htype(static_cast<uint16_t>(h)), etype(static_cast<uint16_t>(e)), length(l) {}

        inline bool Check(uint16_t h, uint16_t e)
        const
        {
            return (htype == h) && (etype == e);
        }

        inline bool Check(HeaderType h, Type e)
        const
        {
            return Check(static_cast<uint16_t>(h), static_cast<uint16_t>(e));
        }

        inline Type GetType(uint16_t h)
        const
        {
            if(htype == h)
                return static_cast<Type>(etype);
            return Type::max_type;
        }
    };

    struct Map
    {
        std::vector<std::vector<int64_t>> maps;

        Map()
        {
            maps.resize(static_cast<size_t>(Type::max_type));
        }

        void Clear()
        {
            for(auto &m : maps) m.clear();
        }

        bool Empty()
        const
        {
            bool empty = true;
            for(auto &m : maps) empty = empty&&m.empty();
            return empty;
        }

        void Add(Type t, int64_t pos)
        {
            size_t i = static_cast<size_t>(t);
            if(i < static_cast<size_t>(Type::max_type))
                maps[i].emplace_back(pos);
        }

        const std::vector<int64_t> &GetType(Type t)
        const
        {
            return maps[static_cast<size_t>(t)];
        }

        std::vector<int64_t> &GetType(Type t)
        {
            return maps[static_cast<size_t>(t)];
        }
    };

public:
    // constructor
    PRadDSTParser(uint32_t size = 1000000);

    // copy/move constructors
    PRadDSTParser(const PRadDSTParser &that) = delete;
    PRadDSTParser(PRadDSTParser &&that) = delete;

    // detructor
    virtual ~PRadDSTParser();

    // copy/move assignment operators
    PRadDSTParser &operator =(const PRadDSTParser &rhs) = delete;
    PRadDSTParser &operator =(PRadDSTParser &&rhs) = delete;

    // public member functions
    void OpenOutput(const std::string &path,
                    std::ios::openmode mode = std::ios::out | std::ios::binary);
    void OpenInput(const std::string &path,
                   std::ios::openmode mode = std::ios::in | std::ios::binary);
    void CloseOutput();
    void CloseInput();
    void ResizeBuffer(uint32_t size);

    // write information
    void WriteEvent();
    void WriteEPICS();
    void Write(const EventData &data);
    void Write(const EpicsData &data);

    // get information from current buffer
    bool Read(int64_t pos = -1);
    bool ReadMap();
    Type EventType() const {return cur_evh.GetType(EventHeader);}
    uint32_t GetBufSize() const {return buf_size;}
    EventData GetEvent() const;
    EpicsData GetEPICS() const;
    const Map &GetInputMap() const {return in_map;}
    const Map &GetOutputMap() const {return out_map;}


private:
    void writeMap();
    void saveBuffer(std::ofstream &ofs, Header evh);
    Header getBuffer(std::ifstream &ifs);

private:
    Map in_map, out_map;
    Header cur_evh;
    int64_t content_length;
    std::ofstream dst_out;
    std::ifstream dst_in;
    char *in_buf, *out_buf;
    uint32_t buf_size;
};

#endif

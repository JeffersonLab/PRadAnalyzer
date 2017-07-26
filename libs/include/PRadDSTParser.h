#ifndef PRAD_DST_PARSER_H
#define PRAD_DST_PARSER_H

#include <vector>
#include <fstream>
#include <string>
#include "PRadException.h"
#include "PRadEventStruct.h"

#define DST_BUF_SIZE 1000000

class PRadDSTParser
{
public:
    friend class PRadDataHandler;

    enum Header : unsigned int
    {
        // headers
        FileHeader = 0xc0c0c0,
        EventHeader = 0xe0e0e0,
        MapHeader = 0xa0a0a0,
    };

    enum class Type : unsigned int
    {
        // event types
        event = 0,
        epics,
        max_type,
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
    PRadDSTParser();

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
    bool Read(int64_t pos = -1);
    bool ReadMap();
    Type EventType() const {return ev_type;}
    const EventData &GetEvent() const {return event;}
    const EpicsData &GetEPICSEvent() const {return epics_event;}
    const Map &GetInputMap() const {return in_map;}
    const Map &GetOutputMap() const {return out_map;}

    // write information
    void WriteEvent() throw(PRadException);
    void WriteEPICS() throw(PRadException);
    void WriteEvent(const EventData &data) throw(PRadException);
    void WriteEPICS(const EpicsData &data) throw(PRadException);

private:
    void readEvent(EventData &data) throw(PRadException);
    void readEPICS(EpicsData &data) throw(PRadException);
    void saveBuffer(std::ofstream &ofs, uint32_t htype, uint32_t etype) throw(PRadException);
    Type getBuffer(std::ifstream &ifs) throw (PRadException);
    void writeMap() throw(PRadException);

private:
    Map in_map, out_map;
    std::ofstream dst_out;
    std::ifstream dst_in;
    int64_t content_length;
    EventData event;
    EpicsData epics_event;
    Type ev_type;
    char in_buf[DST_BUF_SIZE], out_buf[DST_BUF_SIZE];
    uint32_t in_idx, out_idx, in_bufl;
};

#endif

#ifndef PRAD_DST_PARSER_H
#define PRAD_DST_PARSER_H

#include <fstream>
#include <string>
#include "PRadException.h"
#include "PRadEventStruct.h"

#define DST_BUF_SIZE 1000000

class PRadDataHandler;
class PRadEPICSystem;
class PRadHyCalSystem;
class PRadGEMSystem;

class PRadDSTParser
{
public:
    friend class PRadDataHandler;

    enum PRadDSTHeader : unsigned int
    {
        // headers
        FileHeader = 0xc0c0c0,
        EventHeader = 0xe0e0e0,
    };

    enum class Type : unsigned int
    {
        // event types
        event = 0,
        epics,
        epics_map,
        run_info,
        hycal_info,
        gem_info,
        undefined,
    };

    enum class Mode : unsigned int
    {
        // by default it updates all info
        update_gem_info,
        update_hycal_ped,
        update_hycal_cal,
        update_run_info,
        update_epics_map,
    };

public:
    // constructor
    PRadDSTParser(PRadDataHandler *h = nullptr);

    // copy/move constructors
    PRadDSTParser(const PRadDSTParser &that) = delete;
    PRadDSTParser(PRadDSTParser &&that) = delete;

    // detructor
    virtual ~PRadDSTParser();

    // copy/move assignment operators
    PRadDSTParser &operator =(const PRadDSTParser &rhs) = delete;
    PRadDSTParser &operator =(PRadDSTParser &&rhs) = delete;

    // public member functions
    void SetHandler(PRadDataHandler *h) {handler = h;};
    PRadDataHandler *GetHandler() const {return handler;};
    void OpenOutput(const std::string &path,
                    std::ios::openmode mode = std::ios::out | std::ios::binary);
    void OpenInput(const std::string &path,
                   std::ios::openmode mode = std::ios::in | std::ios::binary);
    void CloseOutput();
    void CloseInput();
    void SetMode(uint32_t bit_word) {mode = bit_word;};
    void EnableMode(Mode m) {SET_BIT(mode, static_cast<uint32_t>(m));};
    void DisableMode(Mode m) {CLEAR_BIT(mode, static_cast<uint32_t>(m));};
    bool Read();
    Type EventType() const {return ev_type;};
    const EventData &GetEvent() const {return event;};
    const EpicsData &GetEPICSEvent() const {return epics_event;};

    // write information
    void WriteRunInfo() throw(PRadException);
    void WriteEvent() throw(PRadException);
    void WriteEPICS() throw(PRadException);
    void WriteEvent(const EventData &data) throw(PRadException);
    void WriteEPICS(const EpicsData &data) throw(PRadException);
    void WriteEPICSMap(const PRadEPICSystem *epics) throw(PRadException);
    void WriteHyCalInfo(const PRadHyCalSystem *hycal) throw(PRadException);
    void WriteGEMInfo(const PRadGEMSystem *gem) throw(PRadException);

private:
    void readRunInfo() throw(PRadException);
    void readEvent(EventData &data) throw(PRadException);
    void readEPICS(EpicsData &data) throw(PRadException);
    void readEPICSMap(PRadEPICSystem *epics) throw(PRadException);
    void readHyCalInfo(PRadHyCalSystem *hycal) throw(PRadException);
    void readGEMInfo(PRadGEMSystem *gem) throw(PRadException);
    void writeBuffer(char *ptr, uint32_t size);
    void readBuffer(char *ptr, uint32_t size);
    void saveBuffer(std::ofstream &ofs, uint32_t htype, uint32_t info) throw(PRadException);
    Type getBuffer(std::ifstream &ifs) throw (PRadException);

private:
    PRadDataHandler *handler;
    std::ofstream dst_out;
    std::ifstream dst_in;
    int64_t input_length;
    EventData event;
    EpicsData epics_event;
    Type ev_type;
    char in_buf[DST_BUF_SIZE];
    char out_buf[DST_BUF_SIZE];
    uint32_t in_idx;
    uint32_t out_idx;
    uint32_t in_bufl;
    uint32_t mode;
    bool old_ver;
};

#endif

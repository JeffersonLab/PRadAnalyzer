//===========================================================================//
//  Binary data file reading and writing related functions                   //
//  OBSOLETED! THIS CLASS IS FOR UPDATING OLD DST FILES TO NEW ONES          //
//                                                                           //
//  Chao Peng                                                                //
//  07/04/2016                                                               //
//===========================================================================//
#ifndef DST_READER_V1_H
#define DST_READER_V1_H

#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "PRadException.h"
#include "PRadEventStruct.h"
#include "PRadHyCalSystem.h"
#include "PRadGEMSystem.h"
#include "PRadEPICSystem.h"

#define DST_READER_BUF_SIZE 1000000
#define DST_FILE_VERSION 0x20  // current version
#define DST_FILE_VERSION_OLD 0x13 // supported old version


class DSTReaderV1
{
public:
    friend class PRadDataHandler;

    enum Header : unsigned int
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

public:
    // constructor
    DSTReaderV1();

    // copy/move constructors
    DSTReaderV1(const DSTReaderV1 &that) = delete;
    DSTReaderV1(DSTReaderV1 &&that) = delete;

    // detructor
    virtual ~DSTReaderV1();

    // copy/move assignment operators
    DSTReaderV1 &operator =(const DSTReaderV1 &rhs) = delete;
    DSTReaderV1 &operator =(DSTReaderV1 &&rhs) = delete;

    // public member functions
    void OpenInput(const std::string &path,
                   std::ios::openmode mode = std::ios::in | std::ios::binary);
    void CloseInput();
    bool Read(int64_t pos = -1);
    Type EventType() const {return ev_type;}
    const EventData &GetEvent() const {return event;}
    const EpicsData &GetEPICSEvent() const {return epics_event;}

private:
    void readRunInfo() throw(PRadException);
    void readEvent(EventData &data) throw(PRadException);
    void readEPICS(EpicsData &data) throw(PRadException);
    void readEPICSMap() throw(PRadException);
    void readHyCalInfo() throw(PRadException);
    void readGEMInfo() throw(PRadException);
    void readBuffer(char *ptr, uint32_t size);
    Type getBuffer(std::ifstream &ifs) throw (PRadException);

private:
    std::ifstream dst_in;
    int64_t input_length;
    EventData event;
    EpicsData epics_event;
    Type ev_type;
    char in_buf[DST_READER_BUF_SIZE];
    uint32_t in_idx, in_bufl;
    bool old_ver;
};

// helper functions
inline std::string __dst_ver_str(uint32_t ver)
{
    return std::to_string(ver >> 4) + "." + std::to_string(ver & 0xf);
}

inline bool __dst_check_htype(uint32_t header_word, uint32_t type)
{
    if((header_word >> 8) == type)
        return true;

    return false;
}

inline uint32_t __dst_form_header(uint32_t htype, uint32_t info)
{
    return (htype << 8) | info;
}

inline uint32_t __dst_get_ver(uint32_t header_word)
{
    if(__dst_check_htype(header_word, DSTReaderV1::FileHeader))
        return (header_word & 0xff);

    return 0;
}

inline DSTReaderV1::Type __dst_get_type(uint32_t header_word)
{
    if(__dst_check_htype(header_word, DSTReaderV1::EventHeader))
        return static_cast<DSTReaderV1::Type>(header_word & 0xff);

    return DSTReaderV1::Type::undefined;
}

inline uint32_t __dst_buf_to_header(char *buf, uint32_t index)
{
    if(index < 3)
        return 0;

    uint32_t header = buf[index]&0xff;
    header <<= 8;
    header |= buf[index-1]&0xff;
    header <<= 8;
    header |= buf[index-2]&0xff;
    header <<= 8;
    header |= buf[index-3]&0xff;

    return header;
}

template<typename T>
inline T __read(std::istream &is)
{
    T temp;
    is.read((char *) &temp, sizeof(T));
    return temp;
}

// constructor
DSTReaderV1::DSTReaderV1()
: input_length(0), ev_type(Type::undefined), in_idx(0), in_bufl(0), old_ver(false)
{
    // place holder
}

DSTReaderV1::~DSTReaderV1()
{
    // place holder
}

void DSTReaderV1::OpenInput(const std::string &path, std::ios::openmode mode)
{
    CloseInput();
    dst_in.open(path, mode);

    if(!dst_in.is_open()) {
        std::cerr << "DST Parser: Cannot open input file "
                  << "\"" << path << "\"!"
                  << std::endl;
        return;
    }

    dst_in.seekg(0, dst_in.end);
    input_length = dst_in.tellg();
    dst_in.seekg(0, dst_in.beg);

    uint32_t header;
    dst_in.read((char*) &header, sizeof(header));
    uint32_t ver = __dst_get_ver(header);

    if(ver == DST_FILE_VERSION_OLD) {
        old_ver = true;
    } else if(ver == DST_FILE_VERSION) {
        old_ver = false;
    } else {
        std::cerr << "DST Parser: Version mismatch between the file and library. "
                  << std::endl
                  << "Expected version " << __dst_ver_str(DST_FILE_VERSION)
                  << ", the file version is " << __dst_ver_str(ver) << "."
                  << std::endl;
        dst_in.close();
    }
}

void DSTReaderV1::CloseInput()
{
    dst_in.close();
}

void DSTReaderV1::readEvent(EventData &data)
throw(PRadException)
{
    data.clear();

    // event information
    readBuffer((char*) &data.event_number, sizeof(data.event_number));
    readBuffer((char*) &data.type        , sizeof(data.type));
    readBuffer((char*) &data.trigger     , sizeof(data.trigger));
    readBuffer((char*) &data.timestamp   , sizeof(data.timestamp));

    uint32_t adc_size, tdc_size, gem_size, value_size, dsc_size;
    ADC_Data adc;
    TDC_Data tdc;
    DSC_Data dsc;

    readBuffer((char*) &adc_size, sizeof(adc_size));
    for(uint32_t i = 0; i < adc_size; ++i)
    {
        readBuffer((char*) &adc, sizeof(adc));
        data.add_adc(adc);
    }

    readBuffer((char*) &tdc_size, sizeof(tdc_size));
    for(uint32_t i = 0; i < tdc_size; ++i)
    {
        readBuffer((char*) &tdc, sizeof(tdc));
        data.add_tdc(tdc);
    }

    float value;
    readBuffer((char*) &gem_size, sizeof(gem_size));
    for(uint32_t i = 0; i < gem_size; ++i)
    {
        GEM_Data gemhit;
        readBuffer((char*) &gemhit.addr, sizeof(gemhit.addr));
        readBuffer((char*) &value_size, sizeof(value_size));
        for(uint32_t j = 0; j < value_size; ++j)
        {
            readBuffer((char*) &value, sizeof(value));
            gemhit.add_value(value);
        }
        data.add_gemhit(gemhit);
    }

    readBuffer((char*) &dsc_size, sizeof(dsc_size));
    for(uint32_t i = 0; i < dsc_size; ++i)
    {
        readBuffer((char*) &dsc, sizeof(dsc));
        data.add_dsc(dsc);
    }
}

void DSTReaderV1::readEPICS(EpicsData &data)
throw(PRadException)
{
    data.clear();

    readBuffer((char*) &data.event_number, sizeof(data.event_number));

    uint32_t value_size;
    readBuffer((char*) &value_size, sizeof(value_size));

    for(uint32_t i = 0; i < value_size; ++i)
    {
        float value;
        readBuffer((char*) &value, sizeof(value));
        data.values.push_back(value);
    }
}

void DSTReaderV1::readRunInfo()
throw(PRadException)
{
    RunInfo runInfo;

    readBuffer((char*) &runInfo, sizeof(runInfo));
}

void DSTReaderV1::readEPICSMap()
throw(PRadException)
{
    uint32_t ch_size, str_size, id;
    std::string str;
    float value;

    readBuffer((char*) &ch_size, sizeof(ch_size));
    for(uint32_t i = 0; i < ch_size; ++i)
    {
        str = "";
        readBuffer((char*) &str_size, sizeof(str_size));
        for(uint32_t j = 0; j < str_size; ++j)
        {
            char c;
            readBuffer(&c, sizeof(c));
            str.push_back(c);
        }
        readBuffer((char*) &id, sizeof(id));
        readBuffer((char*) &value, sizeof(value));
    }
}

void DSTReaderV1::readHyCalInfo()
throw(PRadException)
{
    uint32_t ch_size;
    readBuffer((char*) &ch_size, sizeof(ch_size));

    for(uint32_t i = 0; i < ch_size; ++i)
    {
        PRadADCChannel::Pedestal ped;
        readBuffer((char*) &ped, sizeof(ped));

        double factor, base_factor;
        readBuffer((char*) &factor, sizeof(factor));
        readBuffer((char*) &base_factor, sizeof(base_factor));

        double gain;
        uint32_t gain_size;
        readBuffer((char*) &gain_size, sizeof(gain_size));

        for(uint32_t j = 0; j < gain_size; ++j)
        {
            readBuffer((char*) &gain, sizeof(gain));
        }
    }
}

void DSTReaderV1::readGEMInfo()
throw(PRadException)
{
    uint32_t apv_size, ped_size;
    readBuffer((char*) &apv_size, sizeof(apv_size));

    for(uint32_t i = 0; i < apv_size; ++i)
    {
        APVAddress addr;
        readBuffer((char*) &addr, sizeof(addr));

        readBuffer((char*) &ped_size, sizeof(ped_size));

        for(uint32_t j = 0; j < ped_size; ++j)
        {
            PRadGEMAPV::Pedestal ped;
            readBuffer((char*) &ped, sizeof(ped));
        }
    }
}

//============================================================================//
// Return type:  false. file end or error                                     //
//               true. successfully read                                      //
//============================================================================//
bool DSTReaderV1::Read(int64_t pos)
{
    if(pos > 0) dst_in.seekg(pos);
    if(dst_in.eof() || dst_in.tellg() >= input_length) return false;

    try {
        ev_type = getBuffer(dst_in);

        // reset in_buf index
        in_idx = 0;
        switch(ev_type)
        {
        case Type::event:
            readEvent(event);
            break;
        case Type::epics:
            readEPICS(epics_event);
            break;
        case Type::epics_map:
            readEPICSMap();
            break;
        case Type::run_info:
            readRunInfo();
            break;
        case Type::hycal_info:
            readHyCalInfo();
            break;
        case Type::gem_info:
            readGEMInfo();
            break;
        default:
            std::cerr << "READ DST ERROR: Undefined buffer type = 0x"
                      << std::hex << std::setw(8) << std::setfill('0') << static_cast<int>(ev_type)
                      << ", incorrect format or corrupted file."
                      << std::endl;
            return false;
        }
        return true;
    } catch(PRadException &e) {
        std::cerr << e.FailureType() << ": " << e.FailureDesc()
                  << std::endl
                  << "Read from DST Aborted!"
                  << std::endl;
        return false;
    } catch(std::exception &e) {
        std::cerr << e.what()
                  << std::endl
                  << "Read from DST Aborted!"
                  << std::endl;
        return false;
    }
}

inline void DSTReaderV1::readBuffer(char *ptr, uint32_t size)
{
    // read directly from dst_in if it is old version
    if(old_ver) {
        dst_in.read(ptr, size);
        return;
    }

    if(in_idx + size >= in_bufl) {
        std::cerr << "exceeds read-in buffer range! "
                  << in_idx + size << ", "
                  << in_bufl << ", "
                  << static_cast<uint32_t>(ev_type)
                  << std::endl;
        in_idx += size;
        exit(-1);
    }

    for(uint32_t i = 0; i < size; ++i, ++in_idx, ++ptr)
    {
        *ptr = in_buf[in_idx];
    }
}

inline DSTReaderV1::Type DSTReaderV1::getBuffer(std::ifstream &ifs)
throw (PRadException)
{
    if(!ifs.is_open())
        throw PRadException("READ DST", "input file is not opened!");

    // read header first
    uint32_t header;
    ifs.read((char*) &header, sizeof(header));
    Type buf_type = __dst_get_type(header);

    // old version does not have buffer length information
    // so we need read byte by byte from file
    if(!old_ver) {
        // read buffer length
        ifs.read((char*) &in_bufl, sizeof(in_bufl));
        if(in_bufl > DST_READER_BUF_SIZE) {
            throw PRadException("READ DST", "read-in buffer size (" + std::to_string(in_bufl) + ") exceeds limit!");
        }
        ifs.read(in_buf, in_bufl);
    }

    // return buffer type
    return buf_type;
}

#endif // DST_READER_V1


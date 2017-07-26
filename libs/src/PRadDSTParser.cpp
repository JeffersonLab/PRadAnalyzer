//===========================================================================//
//  Binary data file reading and writing related functions                   //
//                                                                           //
//  Chao Peng                                                                //
//  07/04/2016                                                               //
//===========================================================================//

#include <iostream>
#include <iomanip>
#include "PRadDSTParser.h"

#define DST_FILE_VERSION 0x21  // current version



//============================================================================//
// Inline helper functions                                                    //
//============================================================================//

// convert version type to human readable string
inline std::string dst_ver_str(uint32_t ver)
{
    return std::to_string(ver >> 4) + "." + std::to_string(ver & 0xf);
}

// check header type, return false if it is different with the input type
inline bool dst_check_htype(uint32_t header_word, uint32_t type)
{
    if((header_word >> 8) == type)
        return true;

    return false;
}

// get version from header
inline uint32_t dst_get_ver(uint32_t header_word)
{
    if(dst_check_htype(header_word, PRadDSTParser::FileHeader))
        return (header_word & 0xff);

    return 0;
}

// form header, combine header type and event type
inline uint32_t dst_form_header(uint32_t htype, uint32_t type)
{
    return (htype << 8) | type;
}

// get event type from designated header type
inline PRadDSTParser::Type dst_get_type(uint32_t header_word, uint32_t htype = PRadDSTParser::EventHeader)
{
    if(dst_check_htype(header_word, htype))
        return static_cast<PRadDSTParser::Type>(header_word & 0xff);

    return PRadDSTParser::Type::max_type;
}

// read from istream
template<typename T>
inline T ist_read(std::istream &is)
{
    T temp;
    is.read((char *) &temp, sizeof(T));
    return temp;
}

// read from istream, auto determine type
template<typename T>
inline void ist_read(std::istream &is, T &val)
{
    val = ist_read<T>(is);
}

// write to ostream
template<typename T>
inline void ost_write(std::ostream &os, T val)
{
    os.write((char *) &val, sizeof(T));
}

// read vector from istream
template<typename T>
inline std::vector<T> read_vector(std::istream &is)
{
    std::vector<T> res;
    uint32_t size = ist_read<uint32_t>(is);
    res.reserve(size);

    for(uint32_t i = 0; i < size; ++i)
        res.push_back(ist_read<T>(is));

    return res;
}

// read vector from istream, auto determine type
template<typename T>
inline void read_vector(std::istream &is, std::vector<T> &vec)
{
    vec = read_vector<T>(is);
}

// write vector to ostream
template<typename T>
inline void write_vector(std::ostream &os, const std::vector<T> &vec)
{
    ost_write(os, (uint32_t)vec.size());
    for(auto &ele : vec)
        ost_write(os, ele);
}

// read from char array
template<typename T>
inline T buf_read(const char *buf, uint32_t max_size, uint32_t &idx)
{
    uint32_t size = sizeof(T);
    if(idx + size >= max_size) {
        std::cerr << "exceeds read-in buffer range! "
                  << idx + size << " >= " << max_size
                  << std::endl;
        return T();
    }

    T result;
    char *ptr = (char*) &result;
    for(uint32_t i = 0; i < size; ++i, ++idx, ++ptr)
    {
        *ptr = buf[idx];
    }

    return result;
}

// read from char array, auto determine type
template<typename T>
inline void buf_read(const char *buf, uint32_t max_size, uint32_t &idx, T &val)
{
    val = buf_read<T>(buf, max_size, idx);
}

// write to char array
template<typename T>
inline void buf_write(char *buf, uint32_t max_size, uint32_t &idx, T val)
{
    uint32_t size = sizeof(T);
    if(idx + size >= max_size) {
        std::cerr << "exceeds buffer size limit! "
                  << idx + size << " > " << max_size
                  << std::endl;
        return;
    }

    char *ptr = (char*) &val;
    for(uint32_t i = 0; i < size; ++i, ++idx, ++ptr)
    {
        buf[idx] = *ptr;
    }
}

// read vector from char array
template<typename T>
inline std::vector<T> read_vector(const char *buf, uint32_t max_size, uint32_t &idx)
{
    std::vector<T> res;
    uint32_t size = buf_read<uint32_t>(buf, max_size, idx);
    res.reserve(size);

    for(uint32_t i = 0; i < size; ++i)
        res.push_back(buf_read<T>(buf, max_size, idx));

    return res;
}

// read vector from char array, auto determine type
template<typename T>
inline void read_vector(const char *buf, uint32_t max_size, uint32_t &idx, std::vector<T> &vec)
{
    vec = read_vector<T>(buf, max_size, idx);
}

// write vector to char array
template<typename T>
inline void write_vector(char *buf, uint32_t max_size, uint32_t &idx, const std::vector<T> &vec)
{
    buf_write(buf, max_size, idx, (uint32_t)vec.size());
    for(auto &ele : vec)
        buf_write(buf, max_size, idx, ele);
}



//============================================================================//
// Constructor, desctructor                                                   //
//============================================================================//

// constructor
PRadDSTParser::PRadDSTParser()
: content_length(0), ev_type(Type::max_type), in_idx(0), out_idx(0),
  in_bufl(0)
{
    // place holder
}

PRadDSTParser::~PRadDSTParser()
{
    // place holder
    CloseOutput();
    CloseInput();
}



//============================================================================//
// Public member functions                                                    //
//============================================================================//

// open output file
void PRadDSTParser::OpenOutput(const std::string &path, std::ios::openmode mode)
{
    CloseOutput();

    dst_out.open(path, mode);

    if(!dst_out.is_open()) {
        std::cerr << "DST Parser: Cannot open output file "
                  << "\"" << path << "\"!"
                  << std::endl;
        return;
    }

    // save header information
    ost_write(dst_out, dst_form_header(FileHeader, DST_FILE_VERSION));
    // save content length
    ost_write(dst_out, (int64_t)dst_out.tellp());
}

// close output file, save file related information (length, map)
void PRadDSTParser::CloseOutput()
{
    if(!dst_out.is_open())
        return;

    // write content length
    int64_t content_end = dst_out.tellp();
    // skip header
    dst_out.seekp(sizeof(uint32_t));
    // write length
    ost_write(dst_out, content_end);
    // back to the position
    dst_out.seekp(content_end);

    // write event map
    writeMap();
    dst_out.close();
}

// open input file
void PRadDSTParser::OpenInput(const std::string &path, std::ios::openmode mode)
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
    int64_t file_length = dst_in.tellg();
    dst_in.seekg(0, dst_in.beg);

    // check version
    uint32_t ver = dst_get_ver(ist_read<uint32_t>(dst_in));
    if(ver != DST_FILE_VERSION) {
        std::cerr << "DST Parser: Version mismatch between the file and library. "
                  << std::endl
                  << "Expected version " << dst_ver_str(DST_FILE_VERSION)
                  << ", the file version is " << dst_ver_str(ver) << "."
                  << std::endl;
        dst_in.close();
    }

    // check length
    content_length = ist_read<int64_t>(dst_in);
    if(content_length > file_length) {
        std::cerr << "DST Parser: Incorrect content length <" << content_length
                  << ">, it exceeds file length <" << file_length << ">."
                  << std::endl;
        dst_in.close();
    }
}

// close input file
void PRadDSTParser::CloseInput()
{
    in_map.Clear();
    dst_in.close();
}

// Read from the opened input
// Return type:  false. file end or error
//               true. successfully read
bool PRadDSTParser::Read(int64_t pos)
{
    if(pos > 0) dst_in.seekg(pos);
    if(dst_in.eof() || dst_in.tellg() >= content_length) return false;

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


//============================================================================//
// Read and write data structures                                             //
//============================================================================//

// write current event
void PRadDSTParser::WriteEvent()
throw(PRadException)
{
    try {
        WriteEvent(event);
    } catch(...) {
        throw;
    }
}

// write given event
void PRadDSTParser::WriteEvent(const EventData &ev)
throw(PRadException)
{
    // event information
    buf_write(out_buf, DST_BUF_SIZE, out_idx, ev.event_number);
    buf_write(out_buf, DST_BUF_SIZE, out_idx, ev.type);
    buf_write(out_buf, DST_BUF_SIZE, out_idx, ev.trigger);
    buf_write(out_buf, DST_BUF_SIZE, out_idx, ev.timestamp);

    // data
    write_vector(out_buf, DST_BUF_SIZE, out_idx, ev.adc_data);
    write_vector(out_buf, DST_BUF_SIZE, out_idx, ev.tdc_data);

    buf_write(out_buf, DST_BUF_SIZE, out_idx, (uint32_t)ev.gem_data.size());
    for(auto &gem : ev.gem_data)
    {
        buf_write(out_buf, DST_BUF_SIZE, out_idx, gem.addr);
        write_vector(out_buf, DST_BUF_SIZE, out_idx, gem.values);
    }

    write_vector(out_buf, DST_BUF_SIZE, out_idx, ev.dsc_data);

    // save buffer to file
    try {
        saveBuffer(dst_out, EventHeader, static_cast<uint32_t>(Type::event));
    } catch(...) {
        throw;
    }
}

// read event
void PRadDSTParser::readEvent(EventData &ev)
throw(PRadException)
{
    ev.clear();

    // event information
    buf_read(in_buf, in_bufl, in_idx, ev.event_number);
    buf_read(in_buf, in_bufl, in_idx, ev.type);
    buf_read(in_buf, in_bufl, in_idx, ev.trigger);
    buf_read(in_buf, in_bufl, in_idx, ev.timestamp);

    // read all data
    read_vector(in_buf, in_bufl, in_idx, ev.adc_data);
    read_vector(in_buf, in_bufl, in_idx, ev.tdc_data);

    // gem data structure is more complicated, it has a vector inside
    uint32_t gem_size = buf_read<uint32_t>(in_buf, in_bufl, in_idx);
    for(uint32_t i = 0; i < gem_size; ++i)
    {
        GEM_Data gemhit;
        buf_read(in_buf, in_bufl, in_idx, gemhit.addr);
        read_vector(in_buf, in_bufl, in_idx, gemhit.values);
        ev.add_gemhit(gemhit);
    }

    read_vector(in_buf, in_bufl, in_idx, ev.dsc_data);
}

// write current epics event
void PRadDSTParser::WriteEPICS()
throw(PRadException)
{
    try {
        WriteEPICS(epics_event);
    } catch(...) {
        throw;
    }
}

// write given epics event
void PRadDSTParser::WriteEPICS(const EpicsData &ep)
throw(PRadException)
{
    buf_write(out_buf, DST_BUF_SIZE, out_idx, ep.event_number);
    write_vector(out_buf, DST_BUF_SIZE, out_idx, ep.values);

    // save buffer to file
    try {
        saveBuffer(dst_out, EventHeader, static_cast<uint32_t>(Type::epics));
    } catch(...) {
        throw;
    }

}

// read epics event
void PRadDSTParser::readEPICS(EpicsData &ep)
throw(PRadException)
{
    ep.clear();

    buf_read(in_buf, in_bufl, in_idx, ep.event_number);
    read_vector(in_buf, in_bufl, in_idx, ep.values);
}

// write file map
void PRadDSTParser::writeMap()
throw(PRadException)
{
    // no need to write map
    if(out_map.Empty()) return;

    // write event map
    ost_write(dst_out, dst_form_header(MapHeader, static_cast<uint32_t>(Type::event)));
    write_vector(dst_out, out_map.GetType(Type::event));

    // write epics map
    ost_write(dst_out, dst_form_header(MapHeader, static_cast<uint32_t>(Type::epics)));
    write_vector(dst_out, out_map.GetType(Type::epics));

    out_map.Clear();
}

// read file map
bool PRadDSTParser::ReadMap()
{
    if(!dst_in.is_open())
        return false;

    in_map.Clear();

    // remember current position
    int64_t cur_pos = dst_in.tellg();

    // to content end, where the map is supposed to be
    dst_in.seekg(content_length);

    // check header
    Type map_type = dst_get_type(ist_read<uint32_t>(dst_in), MapHeader);

    switch(map_type)
    {
    case Type::event:
    case Type::epics:
        read_vector(dst_in, in_map.GetType(map_type));
        break;
    default:
        dst_in.seekg(cur_pos);
        return false;
    }

    dst_in.seekg(cur_pos);
    return true;
}



//============================================================================//
// Private member functions                                                   //
//============================================================================//

inline void PRadDSTParser::saveBuffer(std::ofstream &ofs, uint32_t htype, uint32_t etype)
throw (PRadException)
{
    if(!ofs.is_open())
        throw PRadException("WRITE DST", "output file is not opened!");

    // save current event position
    out_map.Add(static_cast<Type>(etype), ofs.tellp());

    // write header
    uint32_t header = dst_form_header(htype, etype);
    ofs.write((char*) &header, sizeof(header));

    // write buffer length
    ++out_idx;
    ofs.write((char*) &out_idx, sizeof(out_idx));

    // write buffer
    ofs.write(out_buf, out_idx);
    out_idx = 0;
}

inline PRadDSTParser::Type PRadDSTParser::getBuffer(std::ifstream &ifs)
throw (PRadException)
{
    if(!ifs.is_open())
        throw PRadException("READ DST", "input file is not opened!");

    // read header first
    uint32_t header = ist_read<uint32_t>(ifs);
    Type buf_type = dst_get_type(header);

    // read buffer length
    in_bufl = ist_read<uint32_t>(ifs);
    if(in_bufl > DST_BUF_SIZE) {
        throw PRadException("READ DST", "read-in buffer size (" + std::to_string(in_bufl) + ") exceeds limit!");
    }

    // read the whole buffer
    ifs.read(in_buf, in_bufl);

    // return buffer type
    return buf_type;
}

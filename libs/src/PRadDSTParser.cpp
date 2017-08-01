//===========================================================================//
//  Binary data file reading and writing related functions                   //
//                                                                           //
//  Chao Peng                                                                //
//  07/04/2016                                                               //
//===========================================================================//

#include <iostream>
#include <iomanip>
#include "PRadDSTParser.h"

#define DST_FILE_VERSION 0x22  // current version



//============================================================================//
// Inline helper functions                                                    //
//============================================================================//

// convert version type to human readable string
inline std::string dst_ver_str(uint16_t ver)
{
    return std::to_string(ver >> 4) + "." + std::to_string(ver & 0xf);
}

template<typename T>
inline std::streamsize vec_buf_size(const std::vector<T> &vec)
{
    return static_cast<std::streamsize>(vec.size()*sizeof(T));
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

// read from char array
template<typename T>
inline T buf_read(const char *buf, uint32_t max_size, uint32_t &idx)
{
    uint32_t size = sizeof(T);
    if(idx + size > max_size) {
        std::cerr << "exceeds read-in buffer range! "
                  << idx + size << " > " << max_size
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
PRadDSTParser::PRadDSTParser(uint32_t size)
: content_length(0), buf_size(size)
{
    in_buf = new char[size];
    out_buf = new char[size];
}

PRadDSTParser::~PRadDSTParser()
{
    CloseOutput();
    CloseInput();

    delete [] in_buf;
    delete [] out_buf;
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
    ost_write(dst_out, Header(FileHeader, DST_FILE_VERSION, sizeof(int64_t)));

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
    dst_out.seekp(sizeof(Header));
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

    // read file header
    Header header = ist_read<Header>(dst_in);

    if(!header.Check(FileHeader, DST_FILE_VERSION)) {
        std::cerr << "DST Parser: Unsupported format from \"" << path << "\".\n"
                  << "Expected version " << dst_ver_str(DST_FILE_VERSION)
                  << ", the file version is " << dst_ver_str(header.etype) << "."
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
        cur_evh = getBuffer(dst_in);
        Type ev_type = cur_evh.GetType(EventHeader);

        // reset in_buf index
        switch(ev_type)
        {
        case Type::event:
        case Type::epics:
            return true;
        default:
            std::cerr << "READ DST ERROR: Undefined buffer type = 0x"
                      << std::hex << std::setw(8) << std::setfill('0') << static_cast<int>(ev_type)
                      << ", incorrect format or corrupted file."
                      << std::endl;
            return false;
        }
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

// resize buffer
void PRadDSTParser::ResizeBuffer(uint32_t size)
{
    delete [] in_buf;
    delete [] out_buf;
    in_buf = new char[size];
    out_buf = new char[size];
    buf_size = size;
}


//============================================================================//
// Read and write data structures                                             //
//============================================================================//

// write current event
void PRadDSTParser::WriteEvent()
throw(PRadException)
{
    try {
        if(!cur_evh.Check(EventHeader, Type::event)) {
            Write(GetEvent());
        } else {
            throw PRadException("WRITE DST", "Current buffer has no event.");
        }
    } catch(...) {
        throw;
    }
}

// write given event
void PRadDSTParser::Write(const EventData &ev)
throw(PRadException)
{
    uint32_t out_idx = 0;
    // event information
    buf_write(out_buf, buf_size, out_idx, ev.event_number);
    buf_write(out_buf, buf_size, out_idx, ev.type);
    buf_write(out_buf, buf_size, out_idx, ev.trigger);
    buf_write(out_buf, buf_size, out_idx, ev.timestamp);

    // data
    write_vector(out_buf, buf_size, out_idx, ev.adc_data);
    write_vector(out_buf, buf_size, out_idx, ev.tdc_data);

    buf_write(out_buf, buf_size, out_idx, (uint32_t)ev.gem_data.size());
    for(auto &gem : ev.gem_data)
    {
        buf_write(out_buf, buf_size, out_idx, gem.addr);
        write_vector(out_buf, buf_size, out_idx, gem.values);
    }

    write_vector(out_buf, buf_size, out_idx, ev.dsc_data);

    // save buffer to file
    try {
        saveBuffer(dst_out, Header(EventHeader, Type::event, out_idx));
    } catch(...) {
        throw;
    }
}

// read event
EventData PRadDSTParser::GetEvent()
const
{
    EventData ev;
    if(!cur_evh.Check(EventHeader, Type::event)) {
        std::cerr << "DST Parser: Current buffer has no event." << std::endl;
        return ev;
    }

    uint32_t in_idx = 0;
    // event information
    buf_read(in_buf, cur_evh.length, in_idx, ev.event_number);
    buf_read(in_buf, cur_evh.length, in_idx, ev.type);
    buf_read(in_buf, cur_evh.length, in_idx, ev.trigger);
    buf_read(in_buf, cur_evh.length, in_idx, ev.timestamp);

    // read all data
    read_vector(in_buf, cur_evh.length, in_idx, ev.adc_data);
    read_vector(in_buf, cur_evh.length, in_idx, ev.tdc_data);

    // gem data structure is more complicated, it has a vector inside
    uint32_t gem_size = buf_read<uint32_t>(in_buf, cur_evh.length, in_idx);
    for(uint32_t i = 0; i < gem_size; ++i)
    {
        GEM_Data gemhit;
        buf_read(in_buf, cur_evh.length, in_idx, gemhit.addr);
        read_vector(in_buf, cur_evh.length, in_idx, gemhit.values);
        ev.add_gemhit(gemhit);
    }

    read_vector(in_buf, cur_evh.length, in_idx, ev.dsc_data);

    return ev;
}

// write current epics event
void PRadDSTParser::WriteEPICS()
throw(PRadException)
{
    try {
        if(!cur_evh.Check(EventHeader, Type::epics)) {
            Write(GetEPICS());
        } else {
            throw PRadException("WRITE DST", "Current buffer has no EPICS event.");
        }
    } catch(...) {
        throw;
    }
}

// write given epics event
void PRadDSTParser::Write(const EpicsData &ep)
throw(PRadException)
{
    uint32_t out_idx = 0;

    buf_write(out_buf, buf_size, out_idx, ep.event_number);
    write_vector(out_buf, buf_size, out_idx, ep.values);

    // save buffer to file
    try {
        saveBuffer(dst_out, Header(EventHeader, Type::epics, out_idx));
    } catch(...) {
        throw;
    }

}

// read epics event
EpicsData PRadDSTParser::GetEPICS()
const
{
    EpicsData ep;
    if(!cur_evh.Check(EventHeader, Type::epics)) {
        std::cerr << "DST Parser: Current buffer has no EPICS event." << std::endl;
        return ep;
    }

    uint32_t in_idx = 0;
    buf_read(in_buf, cur_evh.length, in_idx, ep.event_number);
    read_vector(in_buf, cur_evh.length, in_idx, ep.values);

    return ep;
}

// write file map
void PRadDSTParser::writeMap()
throw(PRadException)
{
    // directly write to ofstream since map can be huge
    try {
        for(uint16_t i = 0; i < out_map.maps.size(); ++i)
        {
            const auto &cur_map = out_map.GetType(static_cast<Type>(i));
            uint32_t size = cur_map.size();
            ost_write(dst_out, Header(MapHeader, i, vec_buf_size(cur_map) + sizeof(size)));
            ost_write(dst_out, size);
            dst_out.write((char*) &cur_map[0], vec_buf_size(cur_map));
        }

        out_map.Clear();
    } catch(...) {
        out_map.Clear();
        throw;
    }
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

    Header header = ist_read<Header>(dst_in);
    Type map_type = header.GetType(MapHeader);

    switch(map_type)
    {
    case Type::event:
    case Type::epics:
        {
            auto &this_map = in_map.GetType(map_type);
            uint32_t size = ist_read<uint32_t>(dst_in);
            this_map.resize(size);
            dst_in.read((char*) &this_map[0], vec_buf_size(this_map));
        }
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

inline void PRadDSTParser::saveBuffer(std::ofstream &ofs, Header evh)
throw (PRadException)
{
    if(!ofs.is_open())
        throw PRadException("WRITE DST", "output file is not opened!");

    // save current event position
    out_map.Add(static_cast<Type>(evh.etype), ofs.tellp());

    // write header
    ost_write(ofs, evh);

    // write buffer
    ofs.write(out_buf, evh.length);
}

inline PRadDSTParser::Header PRadDSTParser::getBuffer(std::ifstream &ifs)
throw (PRadException)
{
    if(!ifs.is_open())
        throw PRadException("READ DST", "input file is not opened!");

    // read header first
    Header eh = ist_read<Header>(ifs);

    if(eh.length > buf_size) {
        throw PRadException("READ DST", "read-in buffer size (" + std::to_string(cur_evh.length) + ") exceeds limit!");
    }

    // read the whole buffer
    ifs.read(in_buf, eh.length);

    // return buffer type
    return eh;
}

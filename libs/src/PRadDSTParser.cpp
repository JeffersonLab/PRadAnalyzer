//===========================================================================//
//  Binary data file reading and writing related functions                   //
//                                                                           //
//  Chao Peng                                                                //
//  07/04/2016                                                               //
//===========================================================================//

#include <iostream>
#include <iomanip>
#include "PRadDSTParser.h"
#include "PRadDataHandler.h"
#include "PRadEPICSystem.h"
#include "PRadHyCalSystem.h"
#include "PRadGEMSystem.h"
#include "PRadInfoCenter.h"


#define DST_FILE_VERSION 0x20  // current version
#define DST_FILE_VERSION_OLD 0x13 // supported old version

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
    if(__dst_check_htype(header_word, PRadDSTParser::FileHeader))
        return (header_word & 0xff);

    return 0;
}

inline PRadDSTParser::Type __dst_get_type(uint32_t header_word)
{
    if(__dst_check_htype(header_word, PRadDSTParser::EventHeader))
        return static_cast<PRadDSTParser::Type>(header_word & 0xff);

    return PRadDSTParser::Type::undefined;
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

// constructor
PRadDSTParser::PRadDSTParser(PRadDataHandler *h)
: handler(h), input_length(0), ev_type(Type::undefined), in_idx(0), out_idx(0),
  in_bufl(0), mode(0), old_ver(false)
{
    // place holder
}

PRadDSTParser::~PRadDSTParser()
{
    // place holder
}

void PRadDSTParser::OpenOutput(const std::string &path, std::ios::openmode mode)
{
    dst_out.open(path, mode);

    if(!dst_out.is_open()) {
        std::cerr << "DST Parser: Cannot open output file "
                  << "\"" << path << "\"!"
                  << std::endl;
        return;
    }

    // save header information
    uint32_t header = __dst_form_header(FileHeader, DST_FILE_VERSION);
    dst_out.write((char*) &header, sizeof(header));
}

void PRadDSTParser::CloseOutput()
{
    dst_out.close();
}

void PRadDSTParser::OpenInput(const std::string &path, std::ios::openmode mode)
{
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

void PRadDSTParser::CloseInput()
{
    dst_in.close();
}

void PRadDSTParser::WriteEvent()
throw(PRadException)
{
    try {
        WriteEvent(event);
    } catch(...) {
        throw;
    }
}

void PRadDSTParser::WriteEvent(const EventData &data)
throw(PRadException)
{
    // event information
    writeBuffer((char*) &data.event_number, sizeof(data.event_number));
    writeBuffer((char*) &data.type        , sizeof(data.type));
    writeBuffer((char*) &data.trigger     , sizeof(data.trigger));
    writeBuffer((char*) &data.timestamp   , sizeof(data.timestamp));

    // all data banks
    uint32_t adc_size = data.adc_data.size();
    uint32_t tdc_size = data.tdc_data.size();
    uint32_t gem_size = data.gem_data.size();
    uint32_t dsc_size = data.dsc_data.size();

    writeBuffer((char*) &adc_size, sizeof(adc_size));
    for(auto &adc : data.adc_data)
        writeBuffer((char*) &adc, sizeof(adc));

    writeBuffer((char*) &tdc_size, sizeof(tdc_size));
    for(auto &tdc : data.tdc_data)
        writeBuffer((char*) &tdc, sizeof(tdc));

    writeBuffer((char*) &gem_size, sizeof(gem_size));
    for(auto &gem : data.gem_data)
    {
        writeBuffer((char*) &gem.addr, sizeof(gem.addr));
        uint32_t hit_size = gem.values.size();
        writeBuffer((char*) &hit_size, sizeof(hit_size));
        for(auto &value : gem.values)
            writeBuffer((char*) &value, sizeof(value));
    }

    writeBuffer((char*) &dsc_size, sizeof(dsc_size));
    for(auto &dsc : data.dsc_data)
        writeBuffer((char*) &dsc, sizeof(dsc));

    // save buffer to file
    try {
        saveBuffer(dst_out, EventHeader, static_cast<uint32_t>(Type::event));
    } catch(...) {
        throw;
    }

}

void PRadDSTParser::readEvent(EventData &data)
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

void PRadDSTParser::WriteEPICS()
throw(PRadException)
{
    try {
        WriteEPICS(epics_event);
    } catch(...) {
        throw;
    }
}

void PRadDSTParser::WriteEPICS(const EpicsData &data)
throw(PRadException)
{
    writeBuffer((char*) &data.event_number, sizeof(data.event_number));

    uint32_t value_size = data.values.size();
    writeBuffer((char*) &value_size, sizeof(value_size));

    for(auto value : data.values)
        writeBuffer((char*) &value, sizeof(value));

    // save buffer to file
    try {
        saveBuffer(dst_out, EventHeader, static_cast<uint32_t>(Type::epics));
    } catch(...) {
        throw;
    }

}

void PRadDSTParser::readEPICS(EpicsData &data)
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

void PRadDSTParser::WriteRunInfo()
throw(PRadException)
{
    auto runInfo = PRadInfoCenter::Instance().GetRunInfo();
    writeBuffer((char*) &runInfo, sizeof(runInfo));

    // save buffer to file
    try {
        saveBuffer(dst_out, EventHeader, static_cast<uint32_t>(Type::run_info));
    } catch(...) {
        throw;
    }
}

void PRadDSTParser::readRunInfo()
throw(PRadException)
{
    RunInfo runInfo;

    readBuffer((char*) &runInfo, sizeof(runInfo));

    if(TEST_BIT(mode, static_cast<uint32_t>(Mode::update_run_info)))
        PRadInfoCenter::Instance().SetRunInfo(runInfo);
}

void PRadDSTParser::WriteEPICSMap(const PRadEPICSystem *epics)
throw(PRadException)
{
    if(!epics)
        throw PRadException("WRITE DST", "EPICS System does not exist!");

    std::vector<EPICSChannel> channels = epics->GetSortedList();

    uint32_t ch_size = channels.size();
    writeBuffer((char*) &ch_size, sizeof(ch_size));

    for(auto &ch : channels)
    {
        uint32_t str_size = ch.name.size();
        writeBuffer((char*) &str_size, sizeof(str_size));
        for(auto &c : ch.name)
            writeBuffer((char*) &c, sizeof(c));
        writeBuffer((char*) &ch.id, sizeof(ch.id));
        writeBuffer((char*) &ch.value, sizeof(ch.value));
    }

    // save buffer to file
    try {
        saveBuffer(dst_out, EventHeader, static_cast<uint32_t>(Type::epics_map));
    } catch(...) {
        throw;
    }
}

void PRadDSTParser::readEPICSMap(PRadEPICSystem *epics)
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

        if(epics && TEST_BIT(mode, static_cast<uint32_t>(Mode::update_epics_map)))
            epics->AddChannel(str, id, value);
    }
}

void PRadDSTParser::WriteHyCalInfo(const PRadHyCalSystem *hycal)
throw(PRadException)
{
    if(!hycal)
        throw PRadException("WRITE DST", "HyCal system does not exist!");

    uint32_t ch_size = hycal->GetADCList().size();
    writeBuffer((char*) &ch_size, sizeof(ch_size));

    for(auto channel : hycal->GetADCList())
    {
        PRadADCChannel::Pedestal ped = channel->GetPedestal();
        writeBuffer((char*) &ped, sizeof(ped));

        PRadHyCalModule *module = channel->GetModule();
        PRadCalibConst cal;
        if(module)
            cal = module->GetCalibConst();
        uint32_t gain_size = cal.base_gains.size();
        writeBuffer((char*) &cal.factor, sizeof(cal.factor));
        writeBuffer((char*) &cal.base_factor, sizeof(cal.base_factor));
        writeBuffer((char*) &gain_size, sizeof(gain_size));
        for(auto &gain : cal.base_gains)
            writeBuffer((char*) &gain, sizeof(gain));
    }

    // save buffer to file
    try {
        saveBuffer(dst_out, EventHeader, static_cast<uint32_t>(Type::hycal_info));
    } catch(...) {
        throw;
    }
}

void PRadDSTParser::readHyCalInfo(PRadHyCalSystem *hycal)
throw(PRadException)
{
    uint32_t ch_size;
    readBuffer((char*) &ch_size, sizeof(ch_size));

    for(uint32_t i = 0; i < ch_size; ++i)
    {
        PRadADCChannel::Pedestal ped;
        readBuffer((char*) &ped, sizeof(ped));

        PRadCalibConst cal;
        readBuffer((char*) &cal.factor, sizeof(cal.factor));
        readBuffer((char*) &cal.base_factor, sizeof(cal.base_factor));

        double gain;
        uint32_t gain_size;
        readBuffer((char*) &gain_size, sizeof(gain_size));

        for(uint32_t j = 0; j < gain_size; ++j)
        {
            readBuffer((char*) &gain, sizeof(gain));
            cal.base_gains.push_back(gain);
        }

        if(!hycal || hycal->GetADCChannel(i))
            continue;

        PRadADCChannel *adc = hycal->GetADCChannel(i);

        if(TEST_BIT(mode, static_cast<uint32_t>(Mode::update_hycal_ped)))
            adc->SetPedestal(ped);
        if(adc->GetModule() && TEST_BIT(mode, static_cast<uint32_t>(Mode::update_hycal_cal)))
            adc->GetModule()->SetCalibConst(cal);
    }
}

void PRadDSTParser::WriteGEMInfo(const PRadGEMSystem *gem)
throw(PRadException)
{
    if(!gem)
        throw PRadException("WRITE DST", "GEM system does not exist!");

    std::vector<PRadGEMAPV *> apv_list = gem->GetAPVList();

    uint32_t apv_size = apv_list.size();
    writeBuffer((char*) &apv_size, sizeof(apv_size));

    for(auto apv : apv_list)
    {
        APVAddress addr = apv->GetAddress();
        writeBuffer((char*) &addr, sizeof(addr));

        std::vector<PRadGEMAPV::Pedestal> ped_list = apv->GetPedestalList();
        uint32_t ped_size = ped_list.size();
        writeBuffer((char*) &ped_size, sizeof(ped_size));

        for(auto &ped : ped_list)
            writeBuffer((char*) &ped, sizeof(ped));
    }

    // save buffer to file
    try {
        saveBuffer(dst_out, EventHeader, static_cast<uint32_t>(Type::gem_info));
    } catch (...) {
        throw;
    }
}

void PRadDSTParser::readGEMInfo(PRadGEMSystem *gem)
throw(PRadException)
{
    uint32_t apv_size, ped_size;
    readBuffer((char*) &apv_size, sizeof(apv_size));

    for(uint32_t i = 0; i < apv_size; ++i)
    {
        APVAddress addr;
        readBuffer((char*) &addr, sizeof(addr));

        PRadGEMAPV *apv = nullptr;
        if(gem)
            apv = gem->GetAPV(addr);

        readBuffer((char*) &ped_size, sizeof(ped_size));

        for(uint32_t j = 0; j < ped_size; ++j)
        {
            PRadGEMAPV::Pedestal ped;
            readBuffer((char*) &ped, sizeof(ped));
            if(apv && TEST_BIT(mode, static_cast<uint32_t>(Mode::update_gem_info)))
                apv->UpdatePedestal(ped, j);
        }
    }
}

//============================================================================//
// Return type:  false. file end or error                                     //
//               true. successfully read                                      //
//============================================================================//
bool PRadDSTParser::Read()
{
    try {
        if(dst_in.tellg() < input_length && dst_in.tellg() != -1)
        {
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
                if(handler)
                    readEPICSMap(handler->GetEPICSystem());
                else
                    readEPICSMap(nullptr);
                break;
            case Type::run_info:
                readRunInfo();
                break;
            case Type::hycal_info:
                if(handler)
                    readHyCalInfo(handler->GetHyCalSystem());
                else
                    readHyCalInfo(nullptr);
                break;
            case Type::gem_info:
                if(handler)
                    readGEMInfo(handler->GetGEMSystem());
                else
                    readGEMInfo(nullptr);
                break;
            default:
                std::cerr << "READ DST ERROR: Undefined buffer type, incorrect "
                          << "format or corrupted file."
                          << std::endl;
                return false;
            }

            return true;
        } else {
            // file end
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

inline void PRadDSTParser::writeBuffer(char *ptr, uint32_t size)
{
    if(out_idx + size >= DST_BUF_SIZE) {
        std::cerr << "exceeds buffer size limit!" << std::endl;
        return;
    }

    for(uint32_t i = 0; i < size; ++i, ++out_idx, ++ptr)
    {
        out_buf[out_idx] = *ptr;
    }
}

inline void PRadDSTParser::readBuffer(char *ptr, uint32_t size)
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

inline void PRadDSTParser::saveBuffer(std::ofstream &ofs, uint32_t htype, uint32_t info)
throw (PRadException)
{
    if(!ofs.is_open())
        throw PRadException("WRITE DST", "output file is not opened!");

    // write header
    uint32_t header = __dst_form_header(htype, info);
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
    uint32_t header;
    ifs.read((char*) &header, sizeof(header));
    Type buf_type = __dst_get_type(header);

    // old version does not have buffer length information
    // so we need read byte by byte from file
    if(!old_ver) {
        // read buffer length
        ifs.read((char*) &in_bufl, sizeof(in_bufl));
        if(in_bufl > DST_BUF_SIZE) {
            std::cout << std::hex << in_bufl << std::endl;
            throw PRadException("READ DST", "read-in buffer exceeds size limit!");
        }
        ifs.read(in_buf, in_bufl);
    }

    // return buffer type
    return buf_type;
}

//============================================================================//
// A class to parse data buffer from file or ET                               //
// Make sure the endianness is correct or change the code befire using it     //
//                                                                            //
// Chao Peng                                                                  //
// 02/27/2016                                                                 //
//============================================================================//

#include "PRadEvioParser.h"
#include "PRadDataHandler.h"
#include <sstream>
#include <iostream>
#include <iomanip>

#ifdef MULTI_THREAD
#include <thread>
#define ROC_THREAD_THRES 5000     // open a new thread for large roc buffer size
#endif

#define MAX_BUFFER_SIZE 100000    // buffer to store a evio block
#define BLOCK_HEADER_SIZE 8       // evio block header size


using namespace std;



//============================================================================//
// Constuctor, Destructor                                                     //
//============================================================================//

// constructor
PRadEvioParser::PRadEvioParser(PRadDataHandler *handler)
: myHandler(handler), event_number(0)
{
    // place holder
}

// destructor
PRadEvioParser::~PRadEvioParser()
{
    // place holder
}



//============================================================================//
// Public Member Functions                                                    //
//============================================================================//

// simple binary reading for evio format files
void PRadEvioParser::ReadEvioFile(const char *filepath, int evt, bool verbose)
{
    // evio file is written in binary
    ifstream evio_in(filepath, ios::binary | ios::in);

    if(!evio_in.is_open()) {
        cerr << "Cannot open evio file "
             << "\"" << filepath << "\""
             << endl;
        return;
    }

    // get the total length of file
    evio_in.seekg(0, evio_in.end);
    int64_t length = evio_in.tellg();
    evio_in.seekg(0, evio_in.beg);

    // buffer is to store current event block, make sure it is large enough
    uint32_t *buffer = new uint32_t[MAX_BUFFER_SIZE];

    if(verbose) {
        cout << "Reading evio file " << filepath << endl;
    }

    // parse block, stop when read enough event
    // if evt <= 0, it reads all events
    int count = 0;
    while(evio_in.tellg() < length && evio_in.tellg() != -1)
    {
        try {
            count += parseEvioBlock(evio_in, buffer, evt-count);
        } catch (PRadException &e) {
            cerr << e.FailureType() << ": "
                 << e.FailureDesc() << endl;
            cerr << "Abort reading from file " << filepath << endl;
            break;
        }

        if(evt > 0 && count >= evt)
            break;
    }

    delete [] buffer;

    evio_in.close();
}

// read a event buffer, return its type
int PRadEvioParser::ReadEventBuffer(const void *buf)
{
    return parseEvent((const PRadEventHeader *)buf);
}


//============================================================================//
// Private Member Functions                                                   //
//============================================================================//

// parse a evio block data
int PRadEvioParser::parseEvioBlock(ifstream &in, uint32_t *buf, int max_evt)
throw(PRadException)
{
    streamsize buf_size = sizeof(uint32_t);

    // read the block size
    in.read((char*) &buf[0], buf_size);

    if(buf[0] > MAX_BUFFER_SIZE)
    {
        throw PRadException("Read Evio Block", "buffer size is not enough for the block (size " + to_string(buf[0]) + ")");
    }

    // read the whole block in
    in.read((char*) &buf[1], buf_size*(buf[0] - 1));

    // skip the block header
    uint32_t index = BLOCK_HEADER_SIZE;

    int buffer_cnt = 0;
    // inside a block
    while(index < buf[0])
    {
        int type = parseEvent((const PRadEventHeader *) &buf[index]);

        // only count physics event
        if((type == CODA_Event) ||
           (type == CODA_Sync))
            buffer_cnt ++;

        // to next event buffer
        index += buf[index] + 1;

        if(max_evt > 0 && buffer_cnt >= max_evt)
            break;
    }

    return  buffer_cnt;
}

// parse an evio event
int PRadEvioParser::parseEvent(const PRadEventHeader *header)
{
    // first check event type
    switch(header->tag)
    {
    case CODA_Event:
    case CODA_Sync:
    case EPICS_Info:
        // go on to parse event data
        break;

    case CODA_Prestart:
    case CODA_Go:
    case CODA_End:
    default:
        // no need to parse these events
        return header->tag;
    }

    // inform handler the start of a new event
    myHandler->StartofNewEvent(header->tag);

    // skip current header
    const uint32_t buf_size = header->length - 1;
    const uint32_t *buf = (const uint32_t*) &header[1];
    uint32_t index = 0;

#ifdef MULTI_THREAD
    vector<thread> roc_threads;
#endif

    // parse ROC data
    while(index < buf_size)
    {
        const PRadEventHeader *roc_header = (PRadEventHeader*)&buf[index];
        // skip header size and data size 2 + (length - 1)
        index += roc_header->length + 1;
#ifdef MULTI_THREAD
        // open a new thread for large roc data bank
        if(buf[index] > ROC_THREAD_THRES) {
            roc_threads.emplace_back(&PRadEvioParser::parseROCBank, this, roc_header);
        } else {
            parseROCBank(roc_header);
        }
#else
        parseROCBank(roc_header);
#endif
    }

#ifdef MULTI_THREAD
    for(auto &roc : roc_threads)
    {
        if(roc.joinable()) roc.join();
    }
#endif
    // inform handler the end of this event
    myHandler->EndofThisEvent(event_number);

    // return parsed event type
    return header->tag;
}

// parse ROC data
void PRadEvioParser::parseROCBank(const PRadEventHeader *roc_header)
{
    const uint32_t *buf = (const uint32_t*) &roc_header[1]; // skip current header

    switch(roc_header->tag)
    {
    case PRadTagE:  // Tagger E, ROC id 2
    case PRadSRS_2: // SRS, ROC id 8
    case PRadSRS_1: // SRS, ROC id 7
    case PRadROC_3: // Fastbus, ROC id 6
    case PRadROC_2: // Fastbus, ROC id 5
    case PRadROC_1: // Fastbus, ROC id 4
    case PRadTS:    // VME, ROC id 1
    case EPICS_IOC:
        break; // Interested in ROCs, to next header
    case EVINFO_BANK: // special bank
        event_number = buf[0]; // then skip
    default: // unrecognized ROC, skip
        return;
    }

    uint32_t roc_size = roc_header->length - 1;
    uint32_t index = 0;

    while(index < roc_size)
    {
        const PRadEventHeader *bank_header = (PRadEventHeader*)&buf[index];
        // skip bank header size and data size 2 + (length - 1)
        index += bank_header->length + 1;

        // parse bank
        parseDataBank(bank_header);
    }
}

// parse data banks
void PRadEvioParser::parseDataBank(const PRadEventHeader *data_header)
{
    const uint32_t *buffer = (const uint32_t*) &data_header[1]; // skip current header
    uint32_t dataSize = data_header->length - 1;

    // check the header, skip uninterested ones
    switch(data_header->tag)
    {
    default:
    case LIVE_BANK: // bank contains the live time
    case CONF_BANK: // configuration information
        break;
   case TI_BANK: // Bank 0x4, TI data, contains live time and event type information
        parseTIData(buffer, dataSize, data_header->num);
        break;
    case TDC_BANK:
    case TAG_BANK:
        parseTDCV1190(buffer, dataSize, data_header->num);
        break;
    case DSC_BANK:
        parseDSCData(buffer, dataSize);
        break;
    case FASTBUS_BANK: // Bank 0x7, Fastbus data
        parseADC1881M(buffer);
        break;
    case GEM_BANK: // Bank 0x8, gem data, single FEC right now
        parseGEMData(buffer, dataSize, data_header->num);
        break;
    case EPICS_BANK: // epics information
        parseEPICS(buffer);
        break;
    }
}

// Fastbus ADC1881M data
void PRadEvioParser::parseADC1881M(const uint32_t *data)
{
    // Self defined crate data header
    if((data[0]&0xff0fff00) != ADC1881M_DATABEG) {
        cerr << "Incorrect Fastbus bank header!"
             << "0x" << hex << setw(8) << setfill('0')
             <<  data[0] << endl;
        return;
    }

    // number of boards given by the self defined info word in CODA readout list
    const unsigned char boardNum = data[0]&0xFF;
    unsigned int index = 1, wordCount;
    ADC1881MData adcData;

    adcData.addr.crate = (data[0]>>20)&0xF;

    // parse the data for all boards
    for(unsigned char i = 0; i < boardNum; ++i)
    {
        if(data[index] == ADC1881M_ALIGNMENT) // 64 bit alignment, skip
            index++;
        else if(data[index] == ADC1881M_DATAEND) // self defined, end of crate word
            break;

        adcData.addr.slot = (data[index]>>27)&0x1F;
        wordCount = (data[index]&0x7F) + index;
        while(++index < wordCount)
        {
            if(((data[index]>>27)&0x1F) == (unsigned int)adcData.addr.slot) {
                adcData.addr.channel = (data[index]>>17)&0x3F;
                adcData.val = data[index]&0x3FFF;
                myHandler->FeedData(adcData); // feed data to handler
            } else { // show the error message
                cerr << "*** MISMATCHED CRATE ADDRESS ***" << endl;
                cerr << "GEOGRAPHICAL ADDRESS = "
                     << "0x" << hex << setw(8) << setfill('0') // formating
                     << adcData.addr.slot
                     << endl;
                cerr << "BOARD ADDRESS = "
                     << "0x" << hex << setw(8) << setfill('0')
                     << ((data[index]&0xf8000000)>>27)
                     << endl;
                cerr << "DATA WORD = "
                     << "0x" << hex << setw(8) << setfill('0')
                     << data[index]
                     << endl;
            }
        }
    }

}

// GEM data
void PRadEvioParser::parseGEMData(const uint32_t *data, const uint32_t &size,  const int &fec_id)
{
    // pre-zero-suppressed GEM data are in bank 99
    if(fec_id == 99) {
        parseGEMZeroSupData(data, size);
        return;
    }

    // parse raw GEM data
    GEMRawData gemData;
    uint32_t i = 0;

    while(i < size)
    {
        if((data[i]&0xffffff00) == GEMDATA_APVBEG) {
            gemData.addr.adc_ch = data[i]&0xff;
            gemData.addr.fec_id = (data[i+1] >> 16)&0xff;
            gemData.buf = &data[i+2];
            gemData.size = getAPVDataSize(gemData.buf);

            myHandler->FeedData(gemData);

            i += gemData.size;
        } else {
            ++i;
        }
    }
}

// parse zero-suppressed GEM data
void PRadEvioParser::parseGEMZeroSupData(const uint32_t *data, const uint32_t &size)
{
    // data word structure (32 bit word)
    // detector: 1 bit
    // plane: 1 bit
    // fec id: 4 bit
    // adc channel id: 4 bit
    // strip number: 7 bit
    // time sample number: 3 bit
    // polarity: 1 bit
    // adc value: 11 bit

    if((data[0]&0xffffff00) != GEMDATA_ZEROSUP) {
        cerr << "Unrecognized GEM zero suppressed data header word: "
             << "0x" << hex << setw(8) << setfill('0') << data[0]
             << endl;
    }

    vector<GEMZeroSupData> gemDataPack;
    for(uint32_t i = 1; i < size; ++i)
    {
        GEMZeroSupData gemData;
        gemData.addr.fec_id = (data[i] >> 26)&0xf;
        gemData.addr.adc_ch = (data[i] >> 22)&0xf;
        gemData.channel = (data[i] >> 15)&0x7f;
        gemData.time_sample = (data[i] >> 12)&0x7;
        gemData.adc_value = data[i]&0x7ff;

        gemDataPack.push_back(gemData);
    }

    myHandler->FeedData(gemDataPack);
}

// a helper function to determine the APV data size
uint32_t PRadEvioParser::getAPVDataSize(const uint32_t *data)
{
    uint32_t idx = 0;

    while((data[idx]&0xffffff00) != GEMDATA_APVBEG)
    {
        if(data[idx] == GEMDATA_FECEND) {
            return idx;
        }
        ++idx;
    }

    return idx - 1;
}

// parse CAEN V767 Data
void PRadEvioParser::parseTDCV767(const uint32_t *data, const uint32_t &size, const int &roc_id)
{
    if(!(data[0]&V767_HEADER_BIT)) {
        cerr << "Unrecognized V767 header word: "
             << "0x" << hex << setw(8) << setfill('0') << data[0]
             << endl;
        return;
    }
    if(!(data[size-1]&V767_END_BIT)) {
        cerr << "Unrecognized V767 EOB word: "
             << "0x" << hex << setw(8) << setfill('0') << data[size-1]
             << endl;
        return;
    }

    TDCV767Data tdcData;
    tdcData.addr.crate = roc_id;
    tdcData.addr.slot = data[0]>>27;
    for(uint32_t i = 1; i < size - 1; ++i)
    {
        if(data[i]&V767_INVALID_BIT) {
            cerr << "Event: "<< dec << event_number
                 << ", invalid data word: "
                 << "0x" << hex << setw(8) << setfill('0') << data[i]
                 << endl;
            continue;
        }
        tdcData.addr.channel = (data[i]>>24)&0x7f;
        tdcData.val = data[i]&0xfffff;
        myHandler->FeedData(tdcData);
    }
}

// parse CAEN V1190 Data
void PRadEvioParser::parseTDCV1190(const uint32_t *data, const uint32_t &size, const int &roc_id)
{
    TDCV1190Data tdcData;
    tdcData.addr.crate = roc_id;

    for(uint32_t i = 0; i < size; ++i)
    {
        switch(data[i]>>27)
        {
        case V1190_GLOBAL_HEADER:
            if(roc_id == PRadTS) {
                tdcData.addr.slot = 0; // geo address not supported in this crate
            } else {
                tdcData.addr.slot = data[i]&0x1f;
            }
            break;
        case V1190_TDC_MEASURE:
            tdcData.addr.channel = (data[i]>>19)&0x7f;
            tdcData.val = (data[i]&0x7ffff);
            myHandler->FeedData(tdcData);
            break;
        case V1190_TDC_ERROR:
/*
            cerr << "V1190 Error Word: "
		 << "0x" << hex << setw(8) << setfill('0') << data[i]
                 << endl;
            break;
*/
        case V1190_GLOBAL_TRAILER:
        case V1190_FILLER_WORD:
        default:
            break;
        }
    }
}

// parse JLab distriminator data
void PRadEvioParser::parseDSCData(const uint32_t *data, const uint32_t &size)
{
#define GATED_TDC_GROUP 3
#define GATED_TRG_GROUP 19
#define UNGATED_TDC_GROUP 35
#define UNGATED_TRG_GROUP 51

    if(size < 72) {
        cerr << "Unexpected scalar data bank size: " << size << endl;
        return;
    }

    JLabDSCData dscData;
    dscData.size = 8;
    dscData.gated_buf = &data[GATED_TRG_GROUP];
    dscData.ungated_buf = &data[UNGATED_TRG_GROUP];

    myHandler->FeedData(dscData);
}

// parse JLab TI data
void PRadEvioParser::parseTIData(const uint32_t *data, const uint32_t &size, const int &roc_id)
{
    // update trigger type
    myHandler->UpdateTrgType(bit_to_trigger(data[2]>>24));

    if(roc_id == PRadTS) {// we will be more interested in the TI-master
        // check block header first
        if((data[0] >> 27) != 0x10) {
            cerr << "Unexpected TI block header: "
                 << "0x" << hex << setw(8) << setfill('0') << data[0]
                 << endl;
            return;
        }
        // then check the bank size
        if(size != 9) {
            cerr << "Unexpected TI-master bank size: "
                 << size << endl;
            return;
        }
        JLabTIData tiData;
        tiData.time_low = data[4];
        tiData.time_high = data[5] & 0xffff;
        tiData.latch_word = data[6] & 0xff;
        tiData.lms_phase = (data[8] >> 16) & 0xff;
        myHandler->FeedData(tiData);
    }
}

// parse EPICS data (in textual format
void PRadEvioParser::parseEPICS(const uint32_t *data)
{
    EPICSRawData epics_data;
    epics_data.buf = (const char*) data;

    myHandler->FeedData(epics_data);
}



//============================================================================//
// Public Static Member Functions                                             //
//============================================================================//

PRadTriggerType PRadEvioParser::bit_to_trigger(const unsigned int &bit)
{
    int trg = 0;
    for(; (bit >> trg) > 0; ++trg)
    {
        if(trg >= MAX_Trigger) {
            return Undefined;
        }
    }

    return (PRadTriggerType) trg;
}

unsigned int PRadEvioParser::trigger_to_bit(const PRadTriggerType &trg)
{
    if(trg == NotFromTI)
        return 0;
    else
        return 1 << (int) trg;
}


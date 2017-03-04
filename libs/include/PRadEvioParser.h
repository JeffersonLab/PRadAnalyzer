#ifndef PRAD_EVIO_PARSER_H
#define PRAD_EVIO_PARSER_H

#include <fstream>
#include <cstdint>
#include "datastruct.h"
#include "PRadException.h"

class PRadDataHandler;

class PRadEvioParser
{
public:
    // constructor, destructor
    PRadEvioParser(PRadDataHandler* handler);
    virtual ~PRadEvioParser();

    // public member functions
    void ReadEvioFile(const char *filepath, int evt = -1, bool verbose = false);
    int ReadEventBuffer(const void *buf);

    void SetHandler(PRadDataHandler *h) {myHandler = h;};
    void SetEventNumber(const unsigned int &ev) {event_number = ev;};
    unsigned int GetEventNumber() const {return event_number;};

public:
    // static functions
    static PRadTriggerType bit_to_trigger(const unsigned int &bit);
    static unsigned int trigger_to_bit(const PRadTriggerType &trg);

private:
    // private member functions
    int parseEvioBlock(std::ifstream &s, uint32_t *buf, int max_evt) throw(PRadException);
    int parseEvent(const PRadEventHeader *evt_header);
    void parseROCBank(const PRadEventHeader *roc_header);
    void parseDataBank(const PRadEventHeader *data_header);
    void parseADC1881M(const uint32_t *data);
    void parseGEMData(const uint32_t *data, const uint32_t &size, const int &fec_id);
    void parseGEMZeroSupData(const uint32_t *data, const uint32_t &size);
    void parseTDCV767(const uint32_t *data, const uint32_t &size, const int &roc_id);
    void parseTDCV1190(const uint32_t *data, const uint32_t &size, const int &roc_id);
    void parseDSCData(const uint32_t *data, const uint32_t &size);
    void parseTIData(const uint32_t *data, const uint32_t &size, const int &roc_id);
    void parseEPICS(const uint32_t *data);
    uint32_t getAPVDataSize(const uint32_t *data);

private:
    PRadDataHandler *myHandler;
    unsigned int event_number;
};

#endif

#ifndef PRAD_DATA_STRUCT_H
#define PRAD_DATA_STRUCT_H

#include <cstddef>
#include <cstdint>

enum PRadEventType
{
    CODA_Unknown = 0x0,
    EPICS_Info = 0x1f,
    CODA_Event = 0x81,
    CODA_Prestart = 0x11,
    CODA_Go = 0x12,
    CODA_Sync = 0xc1,
    CODA_End = 0x20,
};

enum PRadTriggerType
{
    NotFromTI = 0,
    PHYS_LeadGlassSum,
    PHYS_TotalSum,
    LMS_Led,
    LMS_Alpha,
    PHYS_TaggerE,
    PHYS_Scintillator,
    MAX_Trigger,
    Undefined,
};

enum EvioBankType
{
    Unknown_32bit = 0x0,
    UnsignedInt_32bit = 0x01,
    Float_32bit = 0x02,
    CharString_8bit = 0x03,
    SignedShort_16bit = 0x04,
    UnsignedShort_16bit = 0x05,
    SignedChar_8bit = 0x06,
    UnsignedChar_8bit = 0x07,
    Double_64bit = 0x08,
    SignedInt_64bit = 0x09,
    UnsignedInt_64bit = 0x0a,
    SignedInt_32bit = 0x0b,
    EvioTagSegment = 0x0c,
    EvioSegment_B = 0x0d,
    EvioBank_B = 0x0e,
    EvioComposite = 0x0f,
    EvioBank = 0x10,
    EvioSegment = 0x20,
    EvioHollerit = 0x21,
    EvioNValue = 0x22,
};

enum PRadROCID
{
    PRadTS = 1,
    PRadTagE = 2,
    PRadROC_1 = 4,
    PRadROC_2 = 5,
    PRadROC_3 = 6,
    PRadSRS_1 = 7,
    PRadSRS_2 = 8,
    EPICS_IOC = 129,
};

enum PRadBankID
{
    TI_BANK = 0xe10a,
    TAG_BANK = 0xe10b,
    CONF_BANK = 0xe10e,
    LIVE_BANK = 0xe112,
    EPICS_BANK = 0xe114,
    DSC_BANK = 0xe115,
    GEM_BANK = 0xe11f,
    FASTBUS_BANK = 0xe120,
    TDC_BANK = 0xe121,
    EVINFO_BANK = 0xc000,
};

struct ChannelAddress
{
    unsigned int crate;
    unsigned int slot;
    unsigned int channel;

    ChannelAddress() {};
    ChannelAddress(const unsigned int &c,
                   const unsigned int &s,
                   const unsigned int &ch)
    : crate(c), slot(s), channel(ch)
    {};

    bool operator < (const ChannelAddress &rhs) const
    {
        if( crate != rhs.crate )
            return crate < rhs.crate ;
        else if( slot != rhs.slot )
            return slot < rhs.slot ;
        else if( channel != rhs.channel )
            return channel < rhs.channel ;
        else
        return false ;
    }

    bool operator == (const ChannelAddress &rhs) const
    {
        if( (crate != rhs.crate) ||
            (slot != rhs.slot)   ||
            (channel != rhs.channel) )
            return false;
        else
            return true;
    }

};

// a simple hash function for DAQ configuration
namespace std
{
    template<>
    struct hash<ChannelAddress>
    {
        size_t operator()(const ChannelAddress &addr)
        const
        {
            // crate id is 1-6, slot is 1-26, channel is 0-63
            // thus they can be filled in a 14 bit word
            // [ 0 0 0 | 0 0 0 0 0 | 0 0 0 0 0 0 ]
            // [ crate |    slot   |   channel   ]
            // this simple hash ensures no collision for current setup
            return ((addr.crate << 11) | (addr.slot << 6) | addr.channel);
        }
    };
}

struct APVAddress
{
    int fec_id;
    int adc_ch;

    APVAddress() {};
    APVAddress(const int &f, const int &a)
    : fec_id(f), adc_ch(a)
    {};

    bool operator < (const APVAddress &rhs) const
    {
        if( fec_id != rhs.fec_id)
            return fec_id < rhs.fec_id;
        else if( adc_ch != rhs.adc_ch)
            return adc_ch < adc_ch;
        else
            return false;
    }

    bool operator == (const APVAddress &rhs) const
    {
        if( (fec_id != rhs.fec_id) ||
            (adc_ch != rhs.adc_ch) )
            return false;
        else
            return true;
    }
};


// some words defined in readout list
#define ADC1881M_DATABEG 0xdc0adc00 //&0xff0fff00
#define ADC1881M_DATAEND 0xfabc0005
#define ADC1881M_ALIGNMENT 0x00000000 // 64 bit data alignment

#define GEMDATA_APVBEG 0x41444300 //&0xffffff00
#define GEMDATA_FECEND 0xfafafafa
#define GEMDATA_ZEROSUP 0xfecfec00 //&xffffff00

#define V767_HEADER_BIT  (1 << 22)
#define V767_END_BIT     (1 << 21)
#define V767_INVALID_BIT (V767_HEADER_BIT | V767_END_BIT)

// v1190 type check (data >> 27)
enum V1190WordType {
    V1190_GLOBAL_HEADER = 0x08,  // 01000
    V1190_GLOBAL_TRAILER = 0x10, // 10000
    V1190_GLOBAL_TIMETAG = 0x11, // 10001
    V1190_TDC_HEADER = 0x01,     // 00001
    V1190_TDC_TRAILER = 0x03,    // 00011
    V1190_TDC_ERROR = 0x04,      // 00100
    V1190_TDC_MEASURE = 0x00,    // 00000
    V1190_FILLER_WORD = 0x18,    // 11000
};

/* 32 bit event header structure
 * -------------------
 * |     length      |
 * -------------------
 * |  tag  |type| num|
 * -------------------
 */
struct PRadEventHeader
{
    unsigned int length;
    unsigned char num;
    unsigned char type;
    unsigned short tag;
};

struct JLabTIData
{
    unsigned char latch_word;
    unsigned char lms_phase;
    unsigned short time_high;
    unsigned int time_low;
};

struct JLabDSCData
{
    unsigned int size;
    const uint32_t *gated_buf;
    const uint32_t *ungated_buf;
};

struct ADC1881MData
{
    ChannelAddress addr;
    unsigned short val;
};

struct TDCV767Data
{
    ChannelAddress addr;
    unsigned int val;
};

struct TDCV1190Data
{
    ChannelAddress addr;
    unsigned int val;
};

struct EPICSRawData
{
    const char *buf;
};

struct GEMRawData
{
    APVAddress addr;
    const uint32_t *buf;
    uint32_t size;
};

struct GEMZeroSupData
{
    APVAddress addr;
    unsigned char channel;
    unsigned char time_sample;
    unsigned short adc_value;
};

#endif

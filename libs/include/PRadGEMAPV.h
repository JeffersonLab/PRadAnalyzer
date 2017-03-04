#ifndef PRAD_GEM_APV_H
#define PRAD_GEM_APV_H

#include <vector>
#include <fstream>
#include <iostream>
#include "PRadEventStruct.h"
#include "datastruct.h"

//1 time sample data have 128 channel
#define APV_CHANNEL_SIZE 128

// 12 words before real time sample data
#define TIME_SAMPLE_DIFF 140 // 12 + APV_CHANNEL_SIZE

// arbitrary number, additional buffer add on time sample data
// this depends on the configuration in the readout list
#define APV_EXTEND_SIZE 130


class PRadGEMFEC;
class PRadGEMPlane;
class TH1I;

class PRadGEMAPV
{
public:
    struct Pedestal
    {
        float offset;
        float noise;

        // initialize with large noise level so there will be no hits instead
        // of maximum hits when gem is not correctly initialized
        Pedestal() : offset(0.), noise(5000.)
        {};
        Pedestal(const float &o, const float &n)
        : offset(o), noise(n)
        {};
    };

    struct StripNb
    {
        unsigned char local;
        int plane;
    };

public:
    // constrcutor
    PRadGEMAPV(const int &orient,
               const int &header_level,
               const std::string &status,
               const uint32_t &time_sample = 3,
               const float &common_threshold = 20.,
               const float &zero_threshold = 5.,
               const float &cross_threshold = 8.);

    // copy/move constructors
    PRadGEMAPV(const PRadGEMAPV &p);
    PRadGEMAPV(PRadGEMAPV &&p);

    // destructor
    virtual ~PRadGEMAPV();

    // copy/move assignment operators
    PRadGEMAPV &operator =(const PRadGEMAPV &p);
    PRadGEMAPV &operator =(PRadGEMAPV &&p);

    // member functions
    void ClearData();
    void ClearPedestal();
    void CreatePedHist();
    void ReleasePedHist();
    void FillPedHist();
    void ResetPedHist();
    void FitPedestal();
    void FillRawData(const uint32_t *buf, const uint32_t &siz);
    void FillZeroSupData(const uint32_t &ch, const uint32_t &ts, const unsigned short &val);
    void FillZeroSupData(const uint32_t &ch, const std::vector<float> &vals);
    void SplitData(const uint32_t &buf, float &word1, float &word2);
    void UpdatePedestal(std::vector<Pedestal> &ped);
    void UpdatePedestal(const Pedestal &ped, const uint32_t &index);
    void UpdatePedestal(const float &offset, const float &noise, const uint32_t &index);
    void ZeroSuppression();
    void CommonModeCorrection(float *buf, const uint32_t &size);
    void CommonModeCorrection_Split(float *buf, const uint32_t &size);
    void CollectZeroSupHits(std::vector<GEM_Data> &hits);
    void CollectZeroSupHits();
    void ResetHitPos();
    void PrintOutPedestal(std::ofstream &out);
    StripNb MapStrip(int ch);
    bool IsCrossTalkStrip(const uint32_t &strip) const;

    // get parameters
    int GetFECID() const {return fec_id;};
    int GetADCChannel() const {return adc_ch;};
    APVAddress GetAddress() const {return APVAddress(fec_id, adc_ch);};
    uint32_t GetNTimeSamples() const {return time_samples;};
    uint32_t GetTimeSampleSize() const {return APV_CHANNEL_SIZE;};
    int GetOrientation() const {return orient;};
    int GetPlaneIndex() const {return plane_index;};
    int GetHeaderLevel() const {return header_level;};
    bool GetSplitStatus() const {return split;};
    float GetCommonModeThresLevel() const {return common_thres;};
    float GetZeroSupThresLevel() const {return zerosup_thres;};
    float GetCrossTalkThresLevel() const {return crosstalk_thres;};
    uint32_t GetBufferSize() const {return buffer_size;};
    int GetLocalStripNb(const uint32_t &ch) const;
    int GetPlaneStripNb(const uint32_t &ch) const;
    PRadGEMFEC *GetFEC() const {return fec;};
    PRadGEMPlane *GetPlane() const {return plane;};
    std::vector<TH1I *> GetHistList() const;
    std::vector<Pedestal> GetPedestalList() const;
    float GetMaxCharge(const uint32_t &ch) const;
    float GetAveragedCharge(const uint32_t &ch) const;
    float GetIntegratedCharge(const uint32_t &ch) const;

    // set parameters
    void SetFEC(PRadGEMFEC *f, int adc_ch, bool force_set = false);
    void UnsetFEC(bool force_unset = false);
    void SetDetectorPlane(PRadGEMPlane *p, int pl_idx, bool force_set = false);
    void UnsetDetectorPlane(bool force_unset = false);
    void SetTimeSample(const uint32_t &t);
    void SetOrientation(const int &o) {orient = o;};
    void SetHeaderLevel(const int &h) {header_level = h;};
    void SetCommonModeThresLevel(const float &t) {common_thres = t;};
    void SetZeroSupThresLevel(const float &t) {zerosup_thres = t;};
    void SetCrossTalkThresLevel(const float &t) {crosstalk_thres = t;};

private:
    void initialize();
    void getAverage(float &ave, const float *buf, const uint32_t &set = 0);
    uint32_t getTimeSampleStart();
    void buildStripMap();

private:
    PRadGEMFEC *fec;
    PRadGEMPlane *plane;
    int fec_id;
    int adc_ch;
    int plane_index;

    uint32_t time_samples;
    int orient;
    int header_level;
    bool split;
    float common_thres;
    float zerosup_thres;
    float crosstalk_thres;
    uint32_t buffer_size;
    uint32_t ts_begin;
    float *raw_data;
    Pedestal pedestal[APV_CHANNEL_SIZE];
    StripNb strip_map[APV_CHANNEL_SIZE];
    bool hit_pos[APV_CHANNEL_SIZE];
    TH1I *offset_hist[APV_CHANNEL_SIZE];
    TH1I *noise_hist[APV_CHANNEL_SIZE];
};

std::ostream &operator <<(std::ostream &os, const APVAddress &ad);

#endif

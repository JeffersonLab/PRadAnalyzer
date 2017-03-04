#ifndef PRAD_GEM_PLANE_H
#define PRAD_GEM_PLANE_H

#include <vector>
#include <string>
#include "PRadGEMAPV.h"

// unit is mm
// X plane did a shift to move the X origin to center hole
// origin shift: 550.4/2 - (44-pitch)/2 = 253.2 mm
#define STRIP_PITCH 0.4
#define X_SHIFT 253.2


class PRadGEMDetector;
class PRadGEMCluster;
// these two structure will be used for cluster reconstruction, defined at the end
struct StripHit;
struct StripCluster;

class PRadGEMPlane
{
public:
    enum PlaneType
    {
        Plane_X,
        Plane_Y,
        Plane_Max
    };

    static const char *GetPlaneTypeName(int enumVal);
    static int GetPlaneTypeID(const char *name);

public:
    // constructors
    PRadGEMPlane(PRadGEMDetector *det = nullptr);
    PRadGEMPlane(const std::string &n, const int &t, const float &s, const int &c,
                 const int &o, const int &d, PRadGEMDetector *det = nullptr);

    // copy/move constructors
    PRadGEMPlane(const PRadGEMPlane &that);
    PRadGEMPlane(PRadGEMPlane &&that);

    // destructor
    virtual ~PRadGEMPlane();

    // copy/move assignment operators
    PRadGEMPlane &operator= (const PRadGEMPlane &rhs);
    PRadGEMPlane &operator= (PRadGEMPlane &&rhs);

    // public member functions
    void ConnectAPV(PRadGEMAPV *apv, const int &index);
    void DisconnectAPV(const uint32_t &plane_index, bool force_disconn);
    void DisconnectAPVs();
    void AddStripHit(int strip, float charge, bool xtalk, int fec, int adc);
    void ClearStripHits();
    void CollectAPVHits();
    float GetStripPosition(const int &plane_strip) const;
    void FormClusters(PRadGEMCluster *method);

    // set parameter
    void SetDetector(PRadGEMDetector *det, bool force_set = false);
    void UnsetDetector(bool force_unset = false);
    void SetName(const std::string &n) {name = n;};
    void SetType(const PlaneType &t) {type = t;};
    void SetSize(const float &s) {size = s;};
    void SetOrientation(const int &o) {orient = o;};
    void SetCapacity(int c);

    // get parameter
    PRadGEMDetector *GetDetector() const {return detector;};
    const std::string &GetName() const {return name;};
    PlaneType GetType() const {return type;};
    float GetSize() const {return size;};
    int GetCapacity() const {return apv_list.size();};
    int GetOrientation() const {return orient;};
    std::vector<PRadGEMAPV*> GetAPVList() const;
    std::vector<StripHit> &GetStripHits() {return strip_hits;};
    const std::vector<StripHit> &GetStripHits() const {return strip_hits;};
    std::vector<StripCluster> &GetStripClusters() {return strip_clusters;};
    const std::vector<StripCluster> &GetStripClusters() const {return strip_clusters;};

private:
    PRadGEMDetector *detector;
    std::string name;
    PlaneType type;
    float size;
    int connector;
    int orient;
    int direction;
    std::vector<PRadGEMAPV*> apv_list;

    // plane raw hits and clusters
    std::vector<StripHit> strip_hits;
    std::vector<StripCluster> strip_clusters;
};

struct StripHit
{
    int strip;
    float charge;
    float position;
    bool cross_talk;
    APVAddress apv_addr;

    StripHit()
    : strip(0), charge(0.), position(0.), cross_talk(false), apv_addr(-1, -1)
    {};
    StripHit(int s, float c, float p, bool f = false, int fec = -1, int adc = -1)
    : strip(s), charge(c), position(p), cross_talk(f), apv_addr(fec, adc)
    {};
};

struct StripCluster
{
    float position;
    float peak_charge;
    float total_charge;
    std::vector<StripHit> hits;

    StripCluster()
    : position(0.), peak_charge(0.), total_charge(0.)
    {};

    StripCluster(const std::vector<StripHit> &p)
    : position(0.), peak_charge(0.), total_charge(0.), hits(p)
    {};

    StripCluster(std::vector<StripHit> &&p)
    : position(0.), peak_charge(0.), total_charge(0.), hits(std::move(p))
    {};

    bool IsCrossTalk()
    const
    {
        for(auto &hit : hits)
        {
            if(!hit.cross_talk)
                return false;
        }
        return true;
    }
};

#endif

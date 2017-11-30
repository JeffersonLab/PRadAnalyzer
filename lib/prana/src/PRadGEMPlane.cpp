//============================================================================//
// GEM plane class                                                            //
// A GEM plane is a component of GEM detector, it should never exist without  //
// a detector, thus the memory of plane will be managed in GEM detector       //
//                                                                            //
// A GEM plane can be connected to several APV units                          //
// GEM hits are collected and grouped on plane level                          //
//                                                                            //
// Chao Peng, Frame work of this class                                        //
// Xinzhan Bai, position and charge calculation method                        //
// 10/07/2016                                                                 //
//============================================================================//

#include "PRadGEMPlane.h"
#include "PRadGEMDetector.h"
#include "PRadGEMCluster.h"
#include "PRadGEMAPV.h"
#include <iostream>
#include <iterator>
#include <algorithm>



//============================================================================//
// constructor, assigment operator, destructor                                //
//============================================================================//

// constructor
PRadGEMPlane::PRadGEMPlane(PRadGEMDetector *det)
: detector(det), name("Undefined"), type(Plane_X), size(0.), orient(0)
{
    // place holder
}

// constructor
PRadGEMPlane::PRadGEMPlane(const std::string &n, const int &t, const float &s,
                           const int &c, const int &o, const int &d, PRadGEMDetector *det)
: detector(det), name(n), size(s), orient(o), direction(d)
{
    type = (PlaneType)t;

    apv_list.resize(c, nullptr);
}

// copy constructor
// connections between it and apv/detector won't be copied
PRadGEMPlane::PRadGEMPlane(const PRadGEMPlane &that)
: detector(nullptr), name(that.name), type(that.type), size(that.size), orient(that.orient),
  direction(that.direction), strip_hits(that.strip_hits), strip_clusters(that.strip_clusters)
{
    apv_list.resize(that.apv_list.size(), nullptr);
}

// move constructor
PRadGEMPlane::PRadGEMPlane(PRadGEMPlane &&that)
: detector(nullptr), name(std::move(that.name)), type(that.type), size(that.size),
  orient(that.orient), direction(that.direction), strip_hits(std::move(strip_hits)),
  strip_clusters(std::move(that.strip_clusters))
{
    apv_list.resize(that.apv_list.size(), nullptr);
}

// destructor
PRadGEMPlane::~PRadGEMPlane()
{
    UnsetDetector();
    DisconnectAPVs();
}

// copy assignment
PRadGEMPlane &PRadGEMPlane::operator =(const PRadGEMPlane &rhs)
{
    if(this == &rhs)
        return *this;

    PRadGEMPlane that(rhs);
    *this = std::move(that);
    return *this;
}

PRadGEMPlane &PRadGEMPlane::operator =(PRadGEMPlane &&rhs)
{
    if(this == &rhs)
        return *this;

    name = std::move(rhs.name);
    type = rhs.type;
    size = rhs.size;
    orient = rhs.orient;
    direction = rhs.direction;

    strip_hits = std::move(rhs.strip_hits);
    strip_clusters = std::move(rhs.strip_clusters);
    return *this;
}



//============================================================================//
// Public Member Functions                                                    //
//============================================================================//

// set detector to the plane
void PRadGEMPlane::SetDetector(PRadGEMDetector *det, bool force_set)
{
    if(det == detector)
        return;

    if(!force_set)
        UnsetDetector();

    detector = det;
}

// disconnect the detector
void PRadGEMPlane::UnsetDetector(bool force_unset)
{
    if(!detector)
        return;

    if(!force_unset)
        detector->DisconnectPlane(type, true);

    detector = nullptr;
}

// change the capacity
void PRadGEMPlane::SetCapacity(int c)
{
    // capacity cannot be negative
    if(c < 0) c = 0;

    if((uint32_t)c < apv_list.size())
    {
        std::cout << "PRad GEM Plane Warning: Reduce the connectors on plane "
                  << name << " from " << apv_list.size() << " to " << c
                  << ". Thus it will lose the connection between APVs that beyond "
                  << c
                  << std::endl;

        for(uint32_t i = c; i < apv_list.size(); ++i)
        {
            if(apv_list[i] != nullptr)
                apv_list[i]->UnsetDetectorPlane(true);
        }
    }

    apv_list.resize(c, nullptr);
}

// connect an APV to the plane
void PRadGEMPlane::ConnectAPV(PRadGEMAPV *apv, const int &index)
{
    if(apv == nullptr)
        return;

    if((uint32_t)index >= apv_list.size()) {
        std::cout << "PRad GEM Plane Warning: Failed to connect plane " << name
                  << " with APV " << apv->GetAddress()
                  << ". Plane connectors are not enough, have " << apv_list.size()
                  << ", this APV is to be connected at " << index
                  << std::endl;
        return;
    }

    if(apv_list[index] != nullptr) {
        std::cout << "PRad GEM Plane Warning: The connector " << index
                  << " of plane " << name << " is connected to APV " << apv->GetAddress()
                  << ", replace the connection."
                  << std::endl;
        return;
    }

    apv_list[index] = apv;
    apv->SetDetectorPlane(this, index);
}

// disconnect an APV
void PRadGEMPlane::DisconnectAPV(const uint32_t &plane_index, bool force_disconn)
{
    if(plane_index >= apv_list.size())
        return;

    auto &apv = apv_list[plane_index];

    if(!apv)
        return;

    if(!force_disconn)
        apv->UnsetDetectorPlane(true);

    apv = nullptr;
}

// reset all APV connections
void PRadGEMPlane::DisconnectAPVs()
{
    for(auto &apv : apv_list)
    {
        if(apv) {
            apv->UnsetDetectorPlane(true);
            apv = nullptr;
        }
    }
}

// get existing APV list
std::vector<PRadGEMAPV*> PRadGEMPlane::GetAPVList()
const
{
    // since the apv list may contain nullptr,
    // only pack connected APVs and return
    std::vector<PRadGEMAPV*> result;

    for(const auto &apv : apv_list)
    {
        if(apv != nullptr)
            result.push_back(apv);
    }

    return result;
}

// calculate strip position by plane strip index
float PRadGEMPlane::GetStripPosition(const int &plane_strip)
const
{
    float position;

    if(type == Plane_X) {
        position = -0.5*(size + 31*STRIP_PITCH) + STRIP_PITCH*plane_strip - X_SHIFT;
    } else {
        position = -0.5*(size - STRIP_PITCH) + STRIP_PITCH*plane_strip;
    }

    return direction*position;
}

// clear the stored plane hits
void PRadGEMPlane::ClearStripHits()
{
    strip_hits.clear();
}

// add a plane hit
void PRadGEMPlane::AddStripHit(int strip, float charge, bool xtalk, int fec, int adc)
{
    // floating strip removal
    // hard coded because it is specific to PRad GEM setting
    if((type == Plane_X) &&
       ((strip < 16) || (strip > 1391)))
       return;

    strip_hits.emplace_back(strip, charge, GetStripPosition(strip), xtalk, fec, adc);
}

// collect hits from the connected APVs
void PRadGEMPlane::CollectAPVHits()
{
    ClearStripHits();

    for(auto &apv : apv_list)
    {
        if(apv != nullptr)
            apv->CollectZeroSupHits();
    }
}

// form clusters by the clustering method
void PRadGEMPlane::FormClusters(PRadGEMCluster *method)
{
    method->FormClusters(strip_hits, strip_clusters);
}

//============================================================================//
// Plane Type Enum Related                                                    //
//============================================================================//
static const char *__plane_type_list[] = {"X", "Y", "Undefined"};

const char *PRadGEMPlane::GetPlaneTypeName(int enumVal)
{
    if(enumVal < 0 || enumVal > (int)Plane_Max)
        return "";

    return __plane_type_list[enumVal];
}

int PRadGEMPlane::GetPlaneTypeID(const char *name)
{
    for(int i = 0; i < (int)Plane_Max; ++i)
        if(strcmp(name, __plane_type_list[i]) == 0)
            return i;

    std::cerr << "PRad GEM Plane Error: Undefined plane type "
              << name << ", please check GetPlaneTypeID() in PRadGEMPlane class."
              << std::endl;
    // not found
    return -1;
}


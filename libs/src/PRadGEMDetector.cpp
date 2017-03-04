//============================================================================//
// GEM detector class                                                         //
// A detector has several planes (currently they are X, Y)                    //
// The planes are managed by detector, thus copy a detector will also copy    //
// the planes that consist of it                                              //
//                                                                            //
// Chao Peng                                                                  //
// 10/07/2016                                                                 //
//============================================================================//

#include "PRadGEMDetector.h"
#include "PRadGEMSystem.h"
#include "PRadGEMCluster.h"
#include "PRadGEMAPV.h"
#include <algorithm>
#include <iostream>



//============================================================================//
// constructor, assigment operator, destructor                                //
//============================================================================//

// constructor
PRadGEMDetector::PRadGEMDetector(const std::string &readoutBoard,
                                 const std::string &detectorType,
                                 const std::string &detector,
                                 PRadGEMSystem *g)
: PRadDetector(detector),
  gem_srs(g), type(detectorType), readout_board(readoutBoard)
{
    planes.resize(PRadGEMPlane::Plane_Max, nullptr);

    gem_hits.reserve(GEM_CLUSTERS_BUFFER);
}

// copy and move assignment will copy or move the planes because planes are 
// managed by detectors, but it won't copy the connection to gem system and so
// won't the id assigned by gem system
// copy constructor
PRadGEMDetector::PRadGEMDetector(const PRadGEMDetector &that)
: PRadDetector(that), gem_srs(nullptr), type(that.type),
  readout_board(that.readout_board), gem_hits(that.gem_hits)
{
    for(auto &plane : that.planes)
    {
        if(plane != nullptr)
            planes.push_back(new PRadGEMPlane(*plane));
        else
            planes.push_back(nullptr);
    }

    ConnectPlanes();
}

// move constructor
PRadGEMDetector::PRadGEMDetector(PRadGEMDetector &&that)
: PRadDetector(that), gem_srs(nullptr), type(std::move(that.type)),
  readout_board(std::move(that.readout_board)), planes(std::move(that.planes)),
  gem_hits(std::move(gem_hits))
{
    // reset the planes' detector
    ConnectPlanes();
}

// destructor
PRadGEMDetector::~PRadGEMDetector()
{
    UnsetSystem();

    for(auto &plane : planes)
    {
        if(plane)
            plane->UnsetDetector(true);

        delete plane;
    }
}

// copy assignment operator
PRadGEMDetector &PRadGEMDetector::operator =(const PRadGEMDetector &rhs)
{
    if(this == &rhs)
        return *this;

    PRadGEMDetector that(rhs); // use copy constructor
    *this = std::move(that); // use move assignment
    return *this;
}

// move assignment operator
PRadGEMDetector &PRadGEMDetector::operator =(PRadGEMDetector &&rhs)
{
    if(this == &rhs)
        return *this;

    PRadDetector::operator=(rhs);
    type = std::move(rhs.type);
    readout_board = std::move(rhs.readout_board);
    planes = std::move(rhs.planes);
    gem_hits = std::move(gem_hits);

    ConnectPlanes();

    return *this;
}

//============================================================================//
// Public Member Functions                                                    //
//============================================================================//

// set a new system, disconnect from the previous one
void PRadGEMDetector::SetSystem(PRadGEMSystem *sys, bool force_set)
{
    if(sys == gem_srs)
        return;

    if(!force_set)
        UnsetSystem();

    gem_srs = sys;
}

void PRadGEMDetector::UnsetSystem(bool force_unset)
{
    if(!gem_srs)
        return;

    if(!force_unset)
        gem_srs->DisconnectDetector(det_id, true);

    gem_srs = nullptr;
}

// add plane to the detector
bool PRadGEMDetector::AddPlane(const int &type,
                               const std::string &name,
                               const double &size,
                               const int &conn,
                               const int &ori,
                               const int &dir)
{
    PRadGEMPlane *new_plane = new PRadGEMPlane(name, type, size, conn, ori, dir);

    // successfully added plane in
    if(AddPlane(new_plane))
        return true;

    delete new_plane;
    return false;
}

// add plane to the detector
bool PRadGEMDetector::AddPlane(PRadGEMPlane *plane)
{
    if(plane == nullptr)
        return false;

    int idx = (int)plane->GetType();

    if(plane->GetDetector() != nullptr) {
        std::cerr << "PRad GEM Detector Error: "
                  << "Trying to add plane " << plane->GetName()
                  << " to detector " << det_name
                  << ", but the plane is belong to " << plane->GetDetector()->GetName()
                  << ". Action aborted."
                  << std::endl;
        return false;
    }

    if(planes[idx] != nullptr) {
        std::cerr << "PRad GEM Detector Error: "
                  << "Trying to add multiple planes with the same type " << idx
                  << "to detector " << det_name << ". Action aborted."
                  << std::endl;
        return false;
    }

    plane->SetDetector(this);
    planes[idx] = plane;
    return true;
}

// remove plane
void PRadGEMDetector::RemovePlane(const int &type)
{
    if((uint32_t)type >= planes.size())
        return;

    auto &plane = planes[type];

    if(plane) {
        plane->UnsetDetector(true);
        delete plane, plane = nullptr;
    }
}

// remove plane
void PRadGEMDetector::DisconnectPlane(const int &type, bool force_disconn)
{
    if((uint32_t)type >= planes.size())
        return;

    auto &plane = planes[type];

    if(!plane)
        return;

    if(!force_disconn)
        plane->UnsetDetector(true);

    plane = nullptr;
}

// Asure the connection to all the planes
void PRadGEMDetector::ConnectPlanes()
{
    for(auto &plane : planes)
    {
        if(plane != nullptr)
            plane->SetDetector(this);
    }
}

// reconstruct hits on planes, need PRadGEMCluster as an input
void PRadGEMDetector::Reconstruct(PRadGEMCluster *gem_recon)
{
    // group strip hits into clusters
    for(auto &plane : planes)
    {
        if(plane == nullptr)
            continue;
        plane->FormClusters(gem_recon);
    }

    // Cartesian reconstruction method
    // reconstruct event hits from clusters
    PRadGEMPlane *plane_x = GetPlane(PRadGEMPlane::Plane_X);
    PRadGEMPlane *plane_y = GetPlane(PRadGEMPlane::Plane_Y);
    // do not have these two planes, cannot reconstruct
    if(!plane_x || !plane_y)
        return;
    gem_recon->CartesianReconstruct(plane_x->GetStripClusters(),
                                    plane_y->GetStripClusters(),
                                    gem_hits,
                                    det_id);
}

// collect all the hits from APVs
void PRadGEMDetector::CollectHits()
{
    for(auto &plane : planes)
    {
        if(plane != nullptr)
            plane->CollectAPVHits();
    }
}

// clear all the hits on plane
void PRadGEMDetector::ClearHits()
{
    for(auto &plane : planes)
    {
        if(plane != nullptr)
            plane->ClearStripHits();
    }
}

// clear all the hits and clusters
void PRadGEMDetector::Reset()
{
    ClearHits();
    gem_hits.clear();
}

PRadGEMPlane *PRadGEMDetector::GetPlane(const std::string &type)
const
{
    uint32_t idx = (uint32_t) PRadGEMPlane::GetPlaneTypeID(type.c_str());

    if(idx >= planes.size())
        return nullptr;

    return planes.at(idx);
}

// get the plane list
std::vector<PRadGEMPlane*> PRadGEMDetector::GetPlaneList()
const
{
    // since it allows nullptr in planes
    // for safety issue, only pack existing planes and return
    std::vector<PRadGEMPlane*> result;

    for(auto &plane : planes)
    {
        if(plane != nullptr)
            result.push_back(plane);
    }

    return result;
}

// get plane by type
PRadGEMPlane *PRadGEMDetector::GetPlane(const int &type)
const
{
    return planes[type];
}

// get apv lists from a plane
std::vector<PRadGEMAPV*> PRadGEMDetector::GetAPVList(const int &type)
const
{
    if(planes[type] == nullptr)
        return std::vector<PRadGEMAPV*>();

    return planes[type]->GetAPVList();
}


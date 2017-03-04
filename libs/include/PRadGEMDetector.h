#ifndef PRAD_GEM_DETECTOR_H
#define PRAD_GEM_DETECTOR_H

#include <vector>
#include <string>
#include "PRadException.h"
#include "PRadGEMPlane.h"
#include "PRadEventStruct.h"
#include "PRadDetector.h"


// reserve space for faster filling of clusters
#define GEM_CLUSTERS_BUFFER 500

class PRadGEMSystem;
class PRadGEMCluster;
class PRadGEMAPV;

class PRadGEMDetector : public PRadDetector
{
public:
    friend class PRadGEMCluster;

public:
    // constructor
    PRadGEMDetector(const std::string &readoutBoard,
                    const std::string &detectorType,
                    const std::string &detector,
                    PRadGEMSystem *g = nullptr);

    // copy/move constructors
    PRadGEMDetector(const PRadGEMDetector &that);
    PRadGEMDetector(PRadGEMDetector &&that);

    // desctructor
    virtual ~PRadGEMDetector();

    // copy/move assignment operators
    PRadGEMDetector &operator =(const PRadGEMDetector &rhs);
    PRadGEMDetector &operator =(PRadGEMDetector &&rhs);

    // public member functions
    void SetSystem(PRadGEMSystem *sys, bool false_set = false);
    void UnsetSystem(bool false_unset = false);
    bool AddPlane(PRadGEMPlane *plane);
    bool AddPlane(const int &type, const std::string &name, const double &size,
                  const int &conn, const int &ori, const int &dir);
    void RemovePlane(const int &type);
    void DisconnectPlane(const int &type, bool false_disconn = false);
    void ConnectPlanes();
    void Reconstruct(PRadGEMCluster *c);
    void CollectHits();
    void ClearHits();
    void Reset();

    // get parameters
    PRadGEMSystem *GetSystem() const {return gem_srs;};
    const std::string &GetType() const {return type;};
    const std::string &GetReadoutBoard() const {return readout_board;};
    PRadGEMPlane *GetPlane(const int &type) const;
    PRadGEMPlane *GetPlane(const std::string &type) const;
    std::vector<PRadGEMPlane*> GetPlaneList() const;
    std::vector<PRadGEMAPV*> GetAPVList(const int &type) const;
    std::vector<GEMHit> &GetHits() {return gem_hits;};
    const std::vector<GEMHit> &GetHits() const {return gem_hits;};

private:
    PRadGEMSystem *gem_srs;
    std::string type;
    std::string readout_board;
    std::vector<PRadGEMPlane*> planes;
    std::vector<GEMHit> gem_hits;
};

#endif

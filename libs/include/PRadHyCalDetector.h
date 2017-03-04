#ifndef PRAD_HYCAL_DETECTOR_H
#define PRAD_HYCAL_DETECTOR_H

#include <vector>
#include <string>
#include <iostream>
#include <unordered_map>
#include "PRadException.h"
#include "PRadHyCalModule.h"
#include "PRadDetector.h"
#include "PRadEventStruct.h"


// value to judge if two modules are connected at corner, quantized to module size
#define CORNER_ADJACENT 1.5
// value to judge if two modules are sharing a side line
#define SIDE_ADJACENT 1.3

class PRadHyCalSystem;
class PRadHyCalCluster;
// these two structure will be used for cluster reconstruction, defined at the end
struct ModuleHit;
struct ModuleCluster;

class PRadHyCalDetector : public PRadDetector
{
public:
    friend class PRadHyCalSystem;
    enum Sector
    {
        // undefined
        Undefined_Sector = -1,
        // normal sectors
        Center = 0,
        Top = 1,
        Right = 2,
        Bottom = 3,
        Left = 4,
        // max number of sectors
        Max_Sector,
    };

public:
    // constructor
    PRadHyCalDetector(const std::string &name = "HyCal", PRadHyCalSystem *sys = nullptr);

    // copy/move constructors
    PRadHyCalDetector(const PRadHyCalDetector &that);
    PRadHyCalDetector(PRadHyCalDetector &&that);

    // desctructor
    virtual ~PRadHyCalDetector();

    // copy/move assignment operators
    PRadHyCalDetector &operator =(const PRadHyCalDetector &rhs);
    PRadHyCalDetector &operator =(PRadHyCalDetector &&rhs);

    // public member functions
    void SetSystem(PRadHyCalSystem *sys, bool force_set = false);
    void UnsetSystem(bool force_unset = false);
    virtual void ReadModuleList(const std::string &path);
    void ReadVirtualModuleList(const std::string &path);
    void ReadCalibrationFile(const std::string &path);
    void SaveModuleList(const std::string &path) const;
    void SaveCalibrationFile(const std::string &path) const;
    bool AddModule(PRadHyCalModule *module);
    void RemoveModule(int id);
    void RemoveModule(const std::string &name);
    void RemoveModule(PRadHyCalModule *module);
    void DisconnectModule(int id, bool force_disconn = false);
    void DisconnectModule(const std::string &name, bool force_disconn = false);
    void DisconnectModule(PRadHyCalModule *module, bool force_disconn = false);
    void SortModuleList();
    void ClearModuleList();
    void OutputModuleList(std::ostream &os) const;
    void Reset();

    // hits/clusters reconstruction
    void Reconstruct(PRadHyCalCluster *method);
    void CreateDeadHits();
    void CollectHits();
    void ClearHits();

    // get parameters
    PRadHyCalSystem *GetSystem() const {return system;};
    PRadHyCalModule *GetModule(const int &primex_id) const;
    PRadHyCalModule *GetModule(const std::string &module_name) const;
    PRadHyCalModule *GetModule(const float &x, const float &y) const;
    double GetEnergy() const;
    const std::vector<PRadHyCalModule*> &GetModuleList() const {return module_list;};
    const std::vector<ModuleHit> &GetModuleHits() const {return module_hits;};
    const std::vector<ModuleCluster> &GetModuleClusters() const {return module_clusters;};
    std::vector<HyCalHit> &GetHits() {return hycal_hits;};
    const std::vector<HyCalHit> &GetHits() const {return hycal_hits;};

public:
    static int get_sector_id(const char *name);
    static const char *get_sector_name(int sec);
    static float hit_distance(const ModuleHit &m1, const ModuleHit &m2);

protected:
    virtual void setLayout(PRadHyCalModule &module) const;

protected:
    PRadHyCalSystem *system;
    std::vector<PRadHyCalModule*> module_list;
    std::unordered_map<int, PRadHyCalModule*> id_map;
    std::unordered_map<std::string, PRadHyCalModule*> name_map;
    std::vector<ModuleHit> module_hits;
    std::vector<ModuleHit> dead_hits;
    std::vector<ModuleCluster> module_clusters;
    std::vector<HyCalHit> hycal_hits;
};

struct ModuleHit
{
    int id;                         // module id
    unsigned int flag;              // module flag
    int sector;                     // hycal sector
    PRadHyCalModule::Geometry geo;  // geometry
    float energy;                   // participated energy, may be splitted
    bool real;                      // false for virtual hit to correct leakage

    ModuleHit(bool r = true)
    : id(0), flag(0), sector(0), energy(0), real(r)
    {};

    ModuleHit(PRadHyCalModule *m, float e, bool r = true)
    : energy(e), real(r)
    {
        id = m->GetID();
        flag = m->GetLayoutFlag();
        sector = m->GetSectorID();
        geo = m->GetGeometry();
    };

    bool operator ==(const ModuleHit &rhs) const {return id == rhs.id;};
};

struct ModuleCluster
{
    ModuleHit center;               // center hit
    std::vector<ModuleHit> hits;    // hits group
    float energy;                   // cluster energy
    float leakage;                  // energy leakage

    ModuleCluster()
    : energy(0), leakage(0)
    {
        hits.reserve(100);
    }

    ModuleCluster(const ModuleHit &hit)
    : center(hit), energy(0), leakage(0)
    {
        hits.reserve(100);
    }

    void AddHit(const ModuleHit &hit)
    {
        hits.emplace_back(hit);
        energy += hit.energy;
    }

    void Merge(const ModuleCluster &that)
    {
        hits.reserve(hits.size() + that.hits.size());
        hits.insert(hits.end(), that.hits.begin(), that.hits.end());

        energy += that.energy;
        leakage += that.leakage;

        if(center.energy < that.center.energy)
            center = that.center;
    }

    void FindCenter()
    {
        if(hits.empty())
            return;

        float max_e = center.energy;
        ModuleHit *cptr = nullptr;
        for(auto &hit : hits)
        {
            if(hit.energy > max_e) {
                max_e = hit.energy;
                cptr = &hit;
            }
        }

        if(cptr)
            center = *cptr;
    }
};

#endif

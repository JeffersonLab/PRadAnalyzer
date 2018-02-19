#ifndef PRAD_HYCAL_DETECTOR_H
#define PRAD_HYCAL_DETECTOR_H

#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <unordered_map>
#include "ConfigParser.h"
#include "PRadException.h"
#include "PRadHyCalModule.h"
#include "PRadDetector.h"
#include "PRadEventStruct.h"



class PRadHyCalSystem;
class PRadHyCalCluster;

class PRadHyCalDetector : public PRadDetector
{
public:
    friend class PRadHyCalSystem;

    enum SectorType
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
    // macro in ConfigParser.h
    ENUM_MAP(SectorType, "Center|Top|Right|Bottom|Left");

    struct SectorInfo
    {
        int id, mtype;
        double msize_x, msize_y;
        std::vector<Point2D<double>> boundpts;

        void Init(int i, PRadHyCalModule *m)
        {
            id = i;
            mtype = m->GetType();
            msize_x = m->GetSizeX();
            msize_y = m->GetSizeY();
        }

        void SetBoundary(double x1, double y1, double x2, double y2)
        {
            boundpts.clear();
            boundpts.emplace_back(x1, y2);
            boundpts.emplace_back(x1, y1);
            boundpts.emplace_back(x2, y1);
            boundpts.emplace_back(x2, y2);
        }

        void GetBoundary(double &x1, double &y1, double &x2, double &y2)
        const
        {
            if(boundpts.size() < 4)
                x1 = 0., y1 = 0., x2 = 0., y2 = 0.;

            auto &min = boundpts[1], &max = boundpts[3];
            x1 = min.x, y1 = min.y, x2 = max.x, y2 = max.y;
        }
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
    virtual bool ReadModuleList(const std::string &path);
    bool ReadVModuleList(const std::string &path);
    bool ReadCalibrationFile(const std::string &path);
    void InitLayout();
    void UpdateSectorInfo();
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
    void ClearVModuleList();
    void OutputModuleList(std::ostream &os) const;
    void Reset();

    // hits related
    void UpdateDeadModules();
    void AddHit(const HyCalHit &hit) {hycal_hits.emplace_back(hit);}
    void AddHit(HyCalHit &&hit) {hycal_hits.emplace_back(hit);}
    void ClearHits() {hycal_hits.clear();}

    // get parameters
    PRadHyCalSystem *GetSystem() const {return system;}
    PRadHyCalModule *GetModule(int primex_id) const;
    PRadHyCalModule *GetModule(double x, double y) const;
    PRadHyCalModule *GetModule(const std::string &module_name) const;
    double GetEnergy() const;
    int GetSectorID(double x, double y) const;
    const std::vector<PRadHyCalModule*> &GetModuleList() const {return module_list;}
    std::vector<HyCalHit> &GetHits() {return hycal_hits;}
    const std::vector<HyCalHit> &GetHits() const {return hycal_hits;}
    const std::vector<SectorInfo> &GetSectorInfo() const {return sector_info;}

    // quantized distance between modules
    double QuantizedDist(const PRadHyCalModule *m1, const PRadHyCalModule *m2) const;
    double QuantizedDist(double x1, double y1, double x2, double y2) const;
    double QuantizedDist(double x1, double y1, int s1, double x2, double y2, int s2) const;
    void QuantizedDist(const PRadHyCalModule *m1, const PRadHyCalModule *m2,
                       double &dx, double &dy) const;
    void QuantizedDist(double x1, double y1, double x2, double y2,
                       double &dx, double &dy) const;
    void QuantizedDist(double x1, double y1, int s1, double x2, double y2, int s2,
                       double &dx, double &dy) const;

public:
    // resolution formula, E, a, b, c should have consistent unit
    static inline double resolution(double E, double a, double b, double c)
    {
        // a/E ++ b/sqrt(E) ++ c, ++ means quadratic sum
        return sqrt(a/E*a/E + b*b/E + c*c);
    }

protected:
    virtual void setLayout(PRadHyCalModule &module) const;

protected:
    PRadHyCalSystem *system;
    std::vector<PRadHyCalModule*> module_list, vmodule_list;
    std::unordered_map<int, PRadHyCalModule*> id_map;
    std::unordered_map<std::string, PRadHyCalModule*> name_map;
    std::vector<HyCalHit> hycal_hits;
    std::vector<SectorInfo> sector_info;
};

#endif


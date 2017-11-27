#ifndef PRAD_HYCAL_MODULE_H
#define PRAD_HYCAL_MODULE_H

#include <cmath>
#include <string>
#include <iostream>
#include "PRadCalibConst.h"
#include "PRadTriggerConst.h"
// Geometry and Layout definition
#include "generalstruct.h"

// crystal module starts from id 1000
#define PWO_ID0 1000


class PRadHyCalDetector;
class PRadADCChannel;
class PRadTDCChannel;

class PRadHyCalModule
{
public:
    friend class PRadHyCalDetector;

    enum ModuleType
    {
        // undefined
        Undefined_Type = -1,
        // normal types
        PbGlass = 0,
        PbWO4 = 1,
        // max number of types
        Max_Type,
    };

    struct Neighbor
    {
        PRadHyCalModule *ptr;
        double dx, dy, dist;

        Neighbor(PRadHyCalModule *p = nullptr, double x = 0., double y = 0.)
        : ptr(p), dx(x), dy(y)
        {
            dist = std::sqrt(dx*dx + dy*dy);
        }
    };

public:
    // constructors
    PRadHyCalModule(int pid,
                    const Geometry &geo = Geometry(),
                    PRadHyCalDetector *det = nullptr);
    PRadHyCalModule(const std::string &name,
                    const Geometry &geo = Geometry(),
                    PRadHyCalDetector *det = nullptr);

    // copy/move constructors
    PRadHyCalModule(const PRadHyCalModule &that);
    PRadHyCalModule(PRadHyCalModule &&that);

    // destructor
    virtual ~PRadHyCalModule();

    // asignment operators
    PRadHyCalModule &operator =(const PRadHyCalModule &rhs);
    PRadHyCalModule &operator =(PRadHyCalModule &&rhs);

    // set members
    void SetDetector(PRadHyCalDetector *det, bool force_set = false);
    void UnsetDetector(bool force_unset = false);
    void SetChannel(PRadADCChannel *ch, bool force_set = false);
    void UnsetChannel(bool force_unset = false);
    void SetGeometry(const Geometry &geo) {geometry = geo;}
    void SetLayout(const Layout &lay) {layout = lay;}
    void SetLayoutFlag(unsigned int &flag) {layout.flag = flag;}
    void SetCalibConst(const PRadCalibConst &c) {cal_const = c;}
    void SetTriggerConst(const PRadTriggerConst &c) {trg_const = c;}
    void GainCorrection(const double &g, const int &ref) {cal_const.GainCorrection(g, ref);}
    void AddNeighbor(PRadHyCalModule *m, double dx, double dy) {neighbors.emplace_back(m, dx, dy);}
    void ClearNeighbors() {neighbors.clear();}
    void RemoveNeighbor(PRadHyCalModule *m);

    // energy related
    double Calibration(const unsigned short &adcVal) const;
    double GetEnergy() const;
    double GetEnergy(const double &value) const;

    // check
    bool IsHyCalModule() const {return (geometry.type == PbGlass) || (geometry.type == PbWO4);}
    bool IsLeadTungstate() const {return geometry.type == PbWO4;}
    bool IsLeadGlass() const {return geometry.type == PbGlass;}
    bool IsNeighbor(int id, bool square_or_circle = true) const;

    // get members
    PRadHyCalDetector *GetDetector() const {return detector;};
    unsigned short GetID() const {return id;}
    const std::string &GetName() const {return name;}
    const Geometry &GetGeometry() const {return geometry;}
    const Layout &GetLayout() const {return layout;}
    const PRadCalibConst &GetCalibConst() const {return cal_const;}
    PRadADCChannel *GetChannel() const {return daq_ch;}

    // get specific information
    std::string GetTypeName() const;
    int GetType() const {return geometry.type;}
    double GetX() const {return geometry.x;}
    double GetY() const {return geometry.y;}
    double GetZ() const {return geometry.z;}
    double GetSizeX() const {return geometry.size_x;}
    double GetSizeY() const {return geometry.size_y;}
    double GetSizeZ() const {return geometry.size_z;}
    void GetBoundary(double &xmin, double &ymin, double &zmin, double &xmax, double &ymax, double &zmax) const;
    int GetSectorID() const {return layout.sector;}
    std::string GetSectorName() const;
    int GetRow() const {return layout.row;}
    int GetColumn() const {return layout.column;}
    unsigned int GetLayoutFlag() const {return layout.flag;}
    double GetCalibrationFactor() const {return cal_const.factor;}
    double GetNonLinearConst() const {return cal_const.non_linear;}
    double GetCalibrationEnergy() const {return cal_const.base_energy;}
    double GetReferenceGain(int ref) const {return cal_const.GetRefGain(ref);}
    double GetTriggerEfficiency(double energy) const {return trg_const.GetTriggerEfficiency(energy);}
    PRadTDCChannel *GetTDC() const;
    const std::vector<Neighbor> &GetNeighbors() {return neighbors;}

    // compare operator
    bool operator < (const PRadHyCalModule &rhs) const
    {
        return id < rhs.id;
    }

public:
    // static functions
    static int name_to_primex_id(const std::string &name);
    static int get_module_type(const char *name);
    static const char *get_module_type_name(int type);
    static double distance(const PRadHyCalModule &m1, const PRadHyCalModule &m2);

protected:
    PRadHyCalDetector *detector;
    PRadADCChannel *daq_ch;
    std::string name;
    int id;
    Geometry geometry;
    Layout layout;
    PRadCalibConst cal_const;
    PRadTriggerConst trg_const;
    std::vector<Neighbor> neighbors;
};

std::ostream &operator <<(std::ostream &os, const PRadHyCalModule &m);
#endif

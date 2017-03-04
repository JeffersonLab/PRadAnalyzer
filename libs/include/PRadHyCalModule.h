#ifndef PRAD_HYCAL_MODULE_H
#define PRAD_HYCAL_MODULE_H

#include <string>
#include <iostream>
#include "PRadCalibConst.h"


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

    struct Geometry
    {
        int type;
        double size_x;
        double size_y;
        double size_z;
        double x;
        double y;
        double z;

        Geometry()
        : type(-1), size_x(0), size_y(0), size_z(0), x(0), y(0), z(0)
        {};

        Geometry(int t, double sx, double sy, double sz,
                double pos_x, double pos_y, double pos_z)
        : type(t), size_x(sx), size_y(sy), size_z(sz), x(pos_x), y(pos_y), z(pos_z)
        {};
    };

    struct Layout
    {
        unsigned int flag;
        int sector;
        int row;
        int column;

        Layout() : flag(0), sector(-1), row(0), column(0)
        {};

        Layout(unsigned int f, int s, int r, int c)
        : flag(f), sector(s), row(r), column(c)
        {};
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
    void SetGeometry(const Geometry &geo) {geometry = geo;};
    void SetLayout(const Layout &lay) {layout = lay;};
    void SetLayoutFlag(unsigned int &flag) {layout.flag = flag;};
    void SetCalibConst(const PRadCalibConst &c) {cal_const = c;};
    void SetTriggerEfficiency(const double &eff) {trigger_eff = eff;};
    void GainCorrection(const double &g, const int &ref) {cal_const.GainCorrection(g, ref);};

    // energy related
    double Calibration(const unsigned short &adcVal) const;
    double GetEnergy() const;
    double GetEnergy(const double &value) const;

    // check type
    bool IsHyCalModule() const {return (geometry.type == PbGlass) || (geometry.type == PbWO4);};
    bool IsLeadTungstate() const {return geometry.type == PbWO4;};
    bool IsLeadGlass() const {return geometry.type == PbGlass;};

    // get members
    unsigned short GetID() const {return id;};
    const std::string &GetName() const {return name;};
    const Geometry &GetGeometry() const {return geometry;};
    const Layout &GetLayout() const {return layout;};
    const PRadCalibConst &GetCalibConst() const {return cal_const;};
    PRadADCChannel *GetChannel() const {return daq_ch;};

    // get specific information
    std::string GetTypeName() const;
    int GetType() const {return geometry.type;};
    double GetX() const {return geometry.x;};
    double GetY() const {return geometry.y;};
    double GetZ() const {return geometry.z;};
    double GetSizeX() const {return geometry.size_x;};
    double GetSizeY() const {return geometry.size_y;};
    double GetSizeZ() const {return geometry.size_z;};
    int GetSectorID() const {return layout.sector;};
    std::string GetSectorName() const;
    int GetRow() const {return layout.row;};
    int GetColumn() const {return layout.column;};
    unsigned int GetLayoutFlag() const {return layout.flag;};
    double GetCalibrationFactor() const {return cal_const.factor;};
    double GetNonLinearConst() const {return cal_const.non_linear;};
    double GetCalibrationEnergy() const {return cal_const.base_energy;};
    double GetReferenceGain(int ref) const {return cal_const.GetRefGain(ref);};
    double GetTriggerEfficiency() const {return trigger_eff;};
    PRadTDCChannel *GetTDC() const;

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
    double trigger_eff;
};

std::ostream &operator <<(std::ostream &os, const PRadHyCalModule &m);
#endif

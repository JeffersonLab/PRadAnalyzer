//============================================================================//
// Basic DAQ channel unit                                                     //
// It could be HyCal Module, LMS PMT or Scintillator                          //
// Chao Peng                                                                  //
// 02/17/2016                                                                 //
//============================================================================//

#include "PRadHyCalModule.h"
#include "PRadHyCalDetector.h"
#include "PRadADCChannel.h"
#include "PRadEventStruct.h"
#include <exception>
#include <iomanip>
#include <cstring>

// enum name lists
static const char *__module_type_list[] = {"PbGlass", "PbWO4"};



//============================================================================//
// Constructors, Destructor, Assignment Operators                             //
//============================================================================//

// constructors
PRadHyCalModule::PRadHyCalModule(const std::string &n,
                                 const Geometry &geo,
                                 PRadHyCalDetector *det)
: detector(det), daq_ch(nullptr), name(n), geometry(geo), trigger_eff(1.0)
{
    id = name_to_primex_id(n);
}

PRadHyCalModule::PRadHyCalModule(int pid, const Geometry &geo, PRadHyCalDetector *det)
: detector(det), daq_ch(nullptr), id(pid), geometry(geo), trigger_eff(1.0)
{
    if(geo.type == PbGlass)
        name = "G";
    if(geo.type == PbWO4)
        name = "W";

    name += std::to_string(id);
}

// copy constructor
PRadHyCalModule::PRadHyCalModule(const PRadHyCalModule &that)
: detector(nullptr), daq_ch(nullptr), name(that.name), id(that.id),
  geometry(that.geometry), layout(that.layout), cal_const(that.cal_const),
  trigger_eff(that.trigger_eff)
{
    // place holder
}

// move constructor
PRadHyCalModule::PRadHyCalModule(PRadHyCalModule &&that)
: detector(nullptr), daq_ch(nullptr), name(std::move(that.name)), id(that.id),
  geometry(that.geometry), layout(that.layout), cal_const(that.cal_const),
  trigger_eff(that.trigger_eff)
{
    // place holder
}

// destructor
PRadHyCalModule::~PRadHyCalModule()
{
    UnsetDetector();
    UnsetChannel();
}

// copy assignment operator
PRadHyCalModule &PRadHyCalModule::operator =(const PRadHyCalModule &rhs)
{
    if(this == &rhs)
        return *this;

    name = rhs.name;
    id = rhs.id;
    geometry = rhs.geometry;
    cal_const = rhs.cal_const;
    trigger_eff = rhs.trigger_eff;
    return *this;
}

PRadHyCalModule &PRadHyCalModule::operator =(PRadHyCalModule &&rhs)
{
    if(this == &rhs)
        return *this;

    name = std::move(rhs.name);
    id = rhs.id;
    geometry = rhs.geometry;
    cal_const = std::move(rhs.cal_const);
    trigger_eff = rhs.trigger_eff;
    return *this;
}

//============================================================================//
// Public Member Functions                                                    //
//============================================================================//

// set detector
void PRadHyCalModule::SetDetector(PRadHyCalDetector *det, bool force_set)
{
    if(det == detector)
        return;

    if(!force_set)
        UnsetDetector();

    detector = det;
}

// disconnect the detector
void PRadHyCalModule::UnsetDetector(bool force_unset)
{
    if(!detector)
        return;

    if(!force_unset)
        detector->DisconnectModule(this, true);

    detector = nullptr;
}

// set daq channel
void PRadHyCalModule::SetChannel(PRadADCChannel *ch, bool force_set)
{
    if(ch == daq_ch)
        return;

    if(!force_set)
        UnsetChannel();

    daq_ch = ch;
}

// disconnect the daq channel
void PRadHyCalModule::UnsetChannel(bool force_unset)
{
    if(!daq_ch)
        return;

    if(!force_unset)
        daq_ch->UnsetModule(true);

    daq_ch = nullptr;
}

// get module type name
std::string PRadHyCalModule::GetTypeName()
const
{
    return std::string(get_module_type_name(geometry.type));
}

// get sector name
std::string PRadHyCalModule::GetSectorName()
const
{
    return std::string(PRadHyCalDetector::get_sector_name(layout.sector));
}

PRadTDCChannel *PRadHyCalModule::GetTDC()
const
{
    if(daq_ch)
        return daq_ch->GetTDC();
    return nullptr;
}

double PRadHyCalModule::GetEnergy()
const
{
    if(daq_ch)
        return cal_const.Calibration(daq_ch->GetReducedValue());

    return 0.;
}

double PRadHyCalModule::GetEnergy(const double &value)
const
{
    if(value > 0.)
        return cal_const.Calibration(value);

    return 0.;
}

//============================================================================//
// Public Static Member Functions                                             //
//============================================================================//

// convert name to primex id, it is highly specific for HyCal setup
int PRadHyCalModule::name_to_primex_id(const std::string &name)
{
    try {
        // lead tungstate module
        if(name.at(0) == 'W' || name.at(0) == 'w')
            return std::stoi(name.substr(1)) + 1000;

        // lead glass module
        if(name.at(0) == 'G' || name.at(0) == 'g')
            return std::stoi(name.substr(1));

    } catch(std::exception &e) {
        std::cerr << e.what() << std::endl;
    }

    // unknown module
    std::cerr << "Cannot auto determine id from mdoule" << name << std::endl;
    return -1;
}

// get enum ModuleType by its name
int PRadHyCalModule::get_module_type(const char *name)
{
    for(int i = 0; i < (int)Max_Type; ++i)
        if(strcmp(name, __module_type_list[i]) == 0)
            return i;

    std::cerr << "PRad HyCal Module Error: Cannot find type " << name
              << ", please check the definition in PRadHyCalModule."
              << std::endl;
    // not found
    return -1;
}

// get name of ModuleType
const char *PRadHyCalModule::get_module_type_name(int type)
{
    if(type < 0 || type >= (int)Max_Type)
        return "Undefined";
    else
        return __module_type_list[type];
}

double PRadHyCalModule::distance(const PRadHyCalModule &m1, const PRadHyCalModule &m2)
{
    double x_dis = m1.geometry.x - m2.geometry.x;
    double y_dis = m1.geometry.y - m2.geometry.y;
    return sqrt(x_dis*x_dis + y_dis*y_dis);
}



//============================================================================//
// Other Functions                                                            //
//============================================================================//

// output hycal module information
std::ostream &operator <<(std::ostream &os, const PRadHyCalModule &m)
{
    // name
    os << std::setw(8) << m.GetName()
    // type
       << std::setw(10) << m.GetTypeName()
    // geometry
       << std::setw(8) << m.GetSizeX()
       << std::setw(8) << m.GetSizeY()
       << std::setw(8) << m.GetSizeZ()
       << std::setw(12) << m.GetX()
       << std::setw(12) << m.GetY()
       << std::setw(12) << m.GetZ();

    return os;
}


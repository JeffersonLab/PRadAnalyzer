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




//============================================================================//
// Constructors, Destructor, Assignment Operators                             //
//============================================================================//

// constructors
PRadHyCalModule::PRadHyCalModule(const std::string &n,
                                 const Geometry &geo,
                                 PRadHyCalDetector *det)
: detector(det), daq_ch(nullptr), name(n), geometry(geo)
{
    id = name_to_primex_id(n);
}

PRadHyCalModule::PRadHyCalModule(int pid, const Geometry &geo, PRadHyCalDetector *det)
: detector(det), daq_ch(nullptr), id(pid), geometry(geo)
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
  trg_const(that.trg_const)
{
    // place holder
}

// move constructor
PRadHyCalModule::PRadHyCalModule(PRadHyCalModule &&that)
: detector(nullptr), daq_ch(nullptr), name(std::move(that.name)), id(that.id),
  geometry(that.geometry), layout(that.layout), cal_const(that.cal_const),
  trg_const(that.trg_const)
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
    trg_const = rhs.trg_const;
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
    trg_const = rhs.trg_const;
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
    return Type2str(geometry.type);
}

// get sector name
std::string PRadHyCalModule::GetSectorName()
const
{
    return std::string(PRadHyCalDetector::SectorType2str(layout.sector));
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
    if(daq_ch && !daq_ch->IsDead())
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

void PRadHyCalModule::GetBoundary(double &xmin, double &ymin, double &zmin,
                                  double &xmax, double &ymax, double &zmax)
const
{
    xmin = geometry.x - geometry.size_x/2.;
    xmax = geometry.x + geometry.size_x/2.;
    ymin = geometry.y - geometry.size_y/2.;
    ymax = geometry.y + geometry.size_y/2.;
    zmin = geometry.z - geometry.size_z/2.;
    zmax = geometry.z + geometry.size_z/2.;
}

void PRadHyCalModule::RemoveNeighbor(PRadHyCalModule *m)
{
    for(auto it = neighbors.begin(); it != neighbors.end(); ++it)
    {
        if(it->ptr->GetID() == m->GetID()) {
            neighbors.erase(it);
            return;
        }
    }
}

bool PRadHyCalModule::IsNeighbor(int id, bool square_or_circle)
const
{
    for(auto &m : neighbors)
    {
        if(id == m->GetID()) {
            // square range
            if(square_or_circle) {
                return (std::abs(m.dx) < 1.01 && std::abs(m.dy) < 1.01);
            // circle range, use 1.2 for the transition region
            } else {
                return m.dist < 1.2;
            }
        }
    }

    return false;
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
            return std::stoi(name.substr(1)) + PWO_ID0;

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


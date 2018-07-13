//============================================================================//
// HyCal detector class                                                       //
//                                                                            //
// Chao Peng                                                                  //
// 11/11/2016                                                                 //
//============================================================================//

#include "PRadHyCalDetector.h"
#include "PRadHyCalSystem.h"
#include "TH1.h"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include "canalib.h"



// a helper function to determine the quantized distance between modules
inline void qdist(double x1, double y1, int s1, double x2, double y2, int s2,
                  const std::vector<PRadHyCalDetector::SectorInfo> &secs,
                  double &dx, double &dy)
{
    const auto &sec1 = secs[s1], &sec2 = secs[s2];

    // in different sections
    if(s1 != s2) {
        // NOTICE highly specific for the current HyCal layout
        // the center sector is for crystal modules
        const auto &center = secs[static_cast<int>(PRadHyCalDetector::Center)];
        const auto &boundary = center.boundpts;

        // in different sectors, and with different types of module
        if(sec1.mtype != sec2.mtype) {
            for(size_t i = 0; i < boundary.size(); ++i)
            {
                // boundary line from two points
                size_t ip = (i == 0) ? boundary.size() - 1 : i - 1;
                auto &p1 = boundary[ip], &p2 = boundary[i];

                double xc = 0., yc = 0.;
                int inter = cana::intersection(p1.x, p1.y, p2.x, p2.y, x1, y1, x2, y2, xc, yc);

                if(inter == 0) {
                    dx = (x2 - xc)/sec2.msize_x + (xc - x1)/sec1.msize_x;
                    dy = (y2 - yc)/sec2.msize_y + (yc - y1)/sec1.msize_y;
                    // there will be only one boundary satisfies all the conditions
                    return;
                }
            }
        // in different sectors but with the same type of module
        // 2 possibilities, pass through the center part or not
        } else {
            double xc[2], yc[2];
            size_t ic = 0;
            for(size_t i = 0; i < boundary.size(); ++i)
            {
                size_t ip = (i == 0) ? boundary.size() - 1 : i - 1;
                auto &p1 = boundary[ip], &p2 = boundary[i];
                int inter = cana::intersection(p1.x, p1.y, p2.x, p2.y, x1, y1, x2, y2,
                                               xc[ic], yc[ic]);

                // found two points
                if(inter == 0 && ic++ > 0) break;
            }

            // across over the centeral part
            if(ic > 1) {
                double dxt = x2 - x1, dyt = y2 - y1;
                double dxc = xc[0] - xc[1], dyc = yc[0] - yc[1];

                // the two segments should have the same sign
                double sign = (dxt*dxc > 0.) ? 1. : -1.;
                dxc *= sign;
                dyc *= sign;

                dx = (dxt - dxc)/sec1.msize_x + dxc/center.msize_x;
                dy = (dyt - dyc)/sec1.msize_y + dyc/center.msize_y;
                return;
            }
        }
    }

    dx = (x2 - x1)/sec1.msize_x;
    dy = (y2 - y1)/sec1.msize_y;
}




//============================================================================//
// constructor, assigment operator, destructor                                //
//============================================================================//

// constructor
PRadHyCalDetector::PRadHyCalDetector(const std::string &det, PRadHyCalSystem *sys)
: PRadDetector(det), system(sys)
{
    // place holder
}

// copy and move assignment will copy or move the modules, but the connection to
// HyCal system and the connections between modules and DAQ units won't be copied
// copy constructor
PRadHyCalDetector::PRadHyCalDetector(const PRadHyCalDetector &that)
: PRadDetector(that), system(nullptr), hycal_hits(that.hycal_hits),
  sector_info(that.sector_info), res_pars(that.res_pars)
{
    for(auto module : that.module_list)
    {
        AddModule(new PRadHyCalModule(*module));
    }

    for(auto vmodule : that.vmodule_list)
    {
        PRadHyCalModule *vm = new PRadHyCalModule(*vmodule);
        vm->SetDetector(this);
        vmodule_list.push_back(vm);
    }
}

// move constructor
PRadHyCalDetector::PRadHyCalDetector(PRadHyCalDetector &&that)
: PRadDetector(that), system(nullptr), module_list(std::move(that.module_list)),
  vmodule_list(std::move(vmodule_list)), id_map(std::move(that.id_map)),
  name_map(std::move(that.name_map)), hycal_hits(std::move(that.hycal_hits)),
  sector_info(std::move(that.sector_info)), res_pars(std::move(that.res_pars))
{
    // reset the connections between module and HyCal
    for(auto module : module_list)
        module->SetDetector(this, true);
    that.module_list.clear();

    for(auto vmodule : vmodule_list)
        vmodule->SetDetector(this, true);
    that.vmodule_list.clear();
}

// destructor
PRadHyCalDetector::~PRadHyCalDetector()
{
    UnsetSystem();

    ClearModuleList();
    ClearVModuleList();
}

// copy assignment operator
PRadHyCalDetector &PRadHyCalDetector::operator =(const PRadHyCalDetector &rhs)
{
    if(this == &rhs)
        return *this;

    PRadHyCalDetector that(rhs); // use copy constructor
    *this = std::move(that);     // use move assignment
    return *this;
}

// move assignment operator
PRadHyCalDetector &PRadHyCalDetector::operator =(PRadHyCalDetector &&rhs)
{
    if(this == &rhs)
        return *this;

    PRadDetector::operator =(rhs);
    module_list = std::move(rhs.module_list);
    vmodule_list = std::move(rhs.vmodule_list);
    id_map = std::move(rhs.id_map);
    name_map = std::move(rhs.name_map);
    hycal_hits = std::move(rhs.hycal_hits);
    sector_info = std::move(rhs.sector_info);
    res_pars = std::move(rhs.res_pars);

    for(auto module : module_list)
        module->SetDetector(this, true);
    rhs.module_list.clear();

    for(auto vmodule : vmodule_list)
        vmodule->SetDetector(this, true);
    rhs.vmodule_list.clear();

    return *this;
}



//============================================================================//
// Public Member Functions                                                    //
//============================================================================//

// set a new system, disconnect from the previous one
void PRadHyCalDetector::SetSystem(PRadHyCalSystem *sys, bool force_set)
{
    if(sys == system)
        return;

    // force to set system without unset original system
    // it will be useful for destruction
    if(!force_set)
        UnsetSystem();

    system = sys;
}

void PRadHyCalDetector::UnsetSystem(bool force_unset)
{
    // only do this if system exists
    if(!system)
        return;

    if(!force_unset)
        system->DisconnectDetector();

    system = nullptr;
}

// read module list
bool PRadHyCalDetector::ReadModuleList(const std::string &path)
{
    if(path.empty())
        return false;

    ConfigParser c_parser;
    if(!c_parser.ReadFile(path)) {
        std::cerr << "PRad HyCal Detector Error: Failed to read module list file "
                  << "\"" << path << "\"."
                  << std::endl;
        return false;
    }

    // clear all modules
    ClearModuleList();

    std::string name;
    std::string type, sector;
    Geometry geo;

    // some info that is not read from list
    while (c_parser.ParseLine())
    {
        if(!c_parser.CheckElements(8))
            continue;

        c_parser >> name >> type
                 >> geo.size_x >> geo.size_y >> geo.size_z
                 >> geo.x >> geo.y >> geo.z;

        geo.type = PRadHyCalModule::str2Type(type.c_str());

        PRadHyCalModule *module = new PRadHyCalModule(name, geo);

        // failed to add module to detector
        if(!AddModule(module))
            delete module;
    }

    // update sector information and build neighbors
    InitLayout();

    return true;
}

// read virtual module list, this helps leakage correction
bool PRadHyCalDetector::ReadVModuleList(const std::string &path)
{
    if(path.empty())
        return false;

    ConfigParser c_parser;
    if(!c_parser.ReadFile(path)) {
        std::cerr << "PRad HyCal Detector Error: Failed to read virtual module"
                  << " list file \"" << path << "\"."
                  << std::endl;
        return false;
    }

    ClearVModuleList();

    std::string name, type;
    Geometry geo;

    // some info that is not read from list
    while (c_parser.ParseLine())
    {
        if(!c_parser.CheckElements(8))
            continue;

        c_parser >> name >> type
                 >> geo.size_x >> geo.size_y >> geo.size_z
                 >> geo.x >> geo.y >> geo.z;

        geo.type = PRadHyCalModule::str2Type(type.c_str());

        PRadHyCalModule *vmodule = new PRadHyCalModule(-1, geo, this);
        vmodule->name = name;
        vmodule->layout.sector = GetSectorID(vmodule->GetX(), vmodule->GetY());
        vmodule_list.push_back(vmodule);
    }

    return true;
}

// read calibration constants file
bool PRadHyCalDetector::ReadCalibrationFile(const std::string &path)
{
    if(path.empty())
        return false;

    ConfigParser c_parser;
    if(!c_parser.ReadFile(path)) {
        std::cerr << "PRad HyCal Detector Error: Failed to read calibration file "
                  << " \"" << path << "\""
                  << std::endl;
        return false;
    }

    std::string name;
    double factor, Ecal, nl;

    while(c_parser.ParseLine())
    {
        if(!c_parser.CheckElements(4, -1)) // more than 4 elements
            continue;

        c_parser >> name >> factor >> Ecal >> nl;
        std::vector<double> gains;
        while(c_parser.NbofElements())
        {
            gains.push_back(c_parser.TakeFirst<double>());
        }

        PRadCalibConst cal_const(factor, Ecal, nl, gains);

        PRadHyCalModule *module = GetModule(name);
        if(module) {
            module->SetCalibConst(cal_const);
        } else {
            std::cout << "PRad HyCal Detector Warning: Cannot find HyCal module "
                      << name << ", skipped its update for calibration constant."
                      << std::endl;
        }
    }

    return true;
}

void PRadHyCalDetector::SaveModuleList(const std::string &path)
const
{
    std::ofstream outf(path);

    OutputModuleList(outf);

    outf.close();
}

void PRadHyCalDetector::SaveCalibrationFile(const std::string &path)
const
{
    std::ofstream outf(path);

    for(auto &module : module_list)
    {
        outf << std::setw(8) << module->GetName()
             << module->GetCalibConst()
             << std::endl;
    }

    outf.close();
}

// add a HyCal module to the detector
bool PRadHyCalDetector::AddModule(PRadHyCalModule *module)
{
    if(module == nullptr)
        return false;

    int id = module->GetID();

    if(id_map.find(id) != id_map.end()) {
        std::cerr << "PRad HyCal Detector Error: "
                  << "Module " << id << " exists, abort adding module "
                  << "with the same id."
                  << std::endl;
        return false;
    }
    const std::string &name = module->GetName();

    if(name_map.find(name) != name_map.end()) {
        std::cerr << "PRad HyCal Detector Error: "
                  << "Module " << name << " exists, abort adding module "
                  << "with the same name."
                  << std::endl;
        return false;
    }

    module->SetDetector(this);
    setLayout(*module);
    module_list.push_back(module);
    name_map[name] = module;
    id_map[id] = module;

    return true;
}

// remove module
void PRadHyCalDetector::RemoveModule(int id)
{
    RemoveModule(GetModule(id));
}

void PRadHyCalDetector::RemoveModule(const std::string &name)
{
    RemoveModule(GetModule(name));
}

void PRadHyCalDetector::RemoveModule(PRadHyCalModule *module)
{
    if(!module)
        return;

    id_map.erase(module->GetID());
    name_map.erase(module->GetName());

    module->UnsetDetector(true);
    for(auto &neighbor : module->neighbors)
    {
        neighbor->RemoveNeighbor(module);
    }

    delete module;

    // rebuild module list
    module_list.clear();
    for(auto &it : id_map)
        module_list.push_back(it.second);
}

// disconnect module
void PRadHyCalDetector::DisconnectModule(int id, bool force_disconn)
{
    DisconnectModule(GetModule(id), force_disconn);
}

void PRadHyCalDetector::DisconnectModule(const std::string &name, bool force_disconn)
{
    DisconnectModule(GetModule(name), force_disconn);
}

void PRadHyCalDetector::DisconnectModule(PRadHyCalModule *module, bool force_disconn)
{
    if(!module)
        return;

    id_map.erase(module->GetID());
    name_map.erase(module->GetName());

    if(!force_disconn)
        module->UnsetDetector(true);

    // rebuild module list
    module_list.clear();
    for(auto &it : id_map)
        module_list.push_back(it.second);
}

void PRadHyCalDetector::SortModuleList()
{
    std::sort(module_list.begin(), module_list.end(),
             [](const PRadHyCalModule *m1, const PRadHyCalModule *m2)
             {
                return *m1 < *m2;
             });
}

void PRadHyCalDetector::ClearModuleList()
{
    for(auto module : module_list)
    {
        // prevent module calling RemoveModule upon destruction
        module->UnsetDetector(true);
        delete module;
    }

    module_list.clear();
    id_map.clear();
    name_map.clear();
}

void PRadHyCalDetector::ClearVModuleList()
{
    for(auto vmodule : vmodule_list)
    {
        vmodule->UnsetDetector(true);
        delete vmodule;
    }

    vmodule_list.clear();
}

void PRadHyCalDetector::OutputModuleList(std::ostream &os)
const
{
    for(auto module : module_list)
    {
        os << *module << std::endl;
    }
}

void PRadHyCalDetector::Reset()
{
    hycal_hits.clear();
}

// prepare dead and virtual hits for the leakage correction in reconstruction
void PRadHyCalDetector::UpdateDeadModules()
{
    // initialize
    for(auto module: module_list)
    {
        module->ClearVirtNeighbors();
        // clear the dead module flag
        CLEAR_BIT(module->layout.flag, kDeadModule);
    }

    // set flag for dead modules
    for(auto module : module_list)
    {
        // module is not connected to a adc channel or the channel is dead
        if(!module->GetChannel() || module->GetChannel()->IsDead()) {
            SET_BIT(module->layout.flag, kDeadModule);
            // set bit for the dead module neighbors
            for(auto &m : module->neighbors)
            {
                SET_BIT(m->layout.flag, kDeadNeighbor);
                m->AddVirtNeighbor(module, -m.dx, -m.dy);
            }
        }
    }

    // check if the module is a neighbor of a dead module, and set a bit for
    // the future correction
    for(auto module : module_list)
    {
        // not boundary modules
        if(!TEST_BIT(module->layout.flag, kInnerBound) &&
           !TEST_BIT(module->layout.flag, kOuterBound))
            continue;

        for(auto &vm : vmodule_list) {
            double dx, dy;
            qdist(module->GetX(), module->GetY(), module->GetSectorID(),
                  vm->GetX(), vm->GetY(), vm->GetSectorID(),
                  sector_info, dx, dy);
            if(std::abs(dx) < 1.01 && std::abs(dy) < 1.01) {
                module->AddVirtNeighbor(vm, dx, dy);
            }
        }
    }
}

PRadHyCalModule *PRadHyCalDetector::GetModule(int id)
const
{
    auto it = id_map.find(id);
    if(it == id_map.end())
        return nullptr;
    return it->second;
}

PRadHyCalModule *PRadHyCalDetector::GetModule(const std::string &name)
const
{
    auto it = name_map.find(name);
    if(it == name_map.end())
        return nullptr;
    return it->second;
}

PRadHyCalModule *PRadHyCalDetector::GetModule(double x, double y)
const
{
    for(auto &module : module_list)
    {
        double pos_x = module->GetX();
        double size_x = module->GetSizeX();
        if((x > pos_x + size_x/2.) || (x < pos_x - size_x/2.))
            continue;
        double pos_y = module->GetY();
        double size_y = module->GetSizeY();
        if((y > pos_y + size_y/2.) || (y < pos_y - size_x/2.))
            continue;

        return module;
    }
    return nullptr;
}

double PRadHyCalDetector::GetEnergy()
const
{
    double energy = 0.;
    for(auto &module : module_list)
    {
        energy += module->GetEnergy();
    }
    return energy;
}

PRadHyCalDetector::ResRegion PRadHyCalDetector::GetResRegion(const PRadHyCalModule *c)
const
{
    if(TEST_BIT(c->GetLayoutFlag(), kTransition)) {
        return Transition;
    } else {
        switch(c->GetType())
        {
        case PRadHyCalModule::PbWO4:
            return PbWO4;
        case PRadHyCalModule::PbGlass:
            return PbGlass;
        default:
            return Undefined_ResRegion;
        }
    }
}

// input/output in the unit of MeV
double PRadHyCalDetector::GetEneRes(const PRadHyCalModule *c, double E)
const
{
    return GetEneRes(GetResRegion(c), E);
}

// input/output in the unit of MeV
double PRadHyCalDetector::GetEneRes(ResRegion r, double E)
const
{
    int rid = static_cast<int>(r);
    if(rid < 0) return 0.;
    return E*resolution(E/1000., res_pars.ene[rid])/100.;
}

// input/output in the unit of MeV/mm
double PRadHyCalDetector::GetPosRes(const PRadHyCalModule *c, double E)
const
{
    return GetPosRes(GetResRegion(c), E);
}

// input/output in the unit of MeV/mm
double PRadHyCalDetector::GetPosRes(ResRegion r, double E)
const
{
    int rid = static_cast<int>(r);
    if(rid < 0) return 0.;
    return resolution(E/1000., res_pars.pos[rid]);
}

// set energy resolution parameters
bool PRadHyCalDetector::SetEneRes(ResRegion r, double a, double b, double c)
{
    int rid = static_cast<int>(r);
    if(rid < 0) return false;

    res_pars.ene[rid][0] = a;
    res_pars.ene[rid][1] = b;
    res_pars.ene[rid][2] = c;
    return true;
}

// set position resolution parameters
bool PRadHyCalDetector::SetPosRes(ResRegion r, double a, double b, double c)
{
    int rid = static_cast<int>(r);
    if(rid < 0) return false;

    res_pars.pos[rid][0] = a;
    res_pars.pos[rid][1] = b;
    res_pars.pos[rid][2] = c;
    return true;
}

// get the sector id for quantized distance, highly specific for HyCal layout
// Notice that out of hycal will also be given a valid id, corresponding to the
// closest lead glass sector, this is intended
int PRadHyCalDetector::GetSectorID(double x, double y)
const
{
    // get the central area
    auto &center = sector_info[static_cast<int>(Center)];
    double x1, y1, x2, y2;
    center.GetBoundary(x1, y1, x2, y2);

    if(x > x2) return (y < y1)? static_cast<int>(Bottom): static_cast<int>(Right);
    if(x < x1) return (y > y2)? static_cast<int>(Top): static_cast<int>(Left);
    if(y > y2) return (x > x2)? static_cast<int>(Right) : static_cast<int>(Top);
    if(y < y1) return (x < x1)? static_cast<int>(Left) : static_cast<int>(Bottom);

    // inside
    return static_cast<int>(Center);
}

void PRadHyCalDetector::InitLayout()
{
    // update sector information first, this is essential for QuantizedDist
    UpdateSectorInfo();

    // update neighbors for each module
    for(auto &module : module_list)
        module->ClearNeighbors();

    for(auto it = module_list.begin(); it != module_list.end(); ++it)
    {
        for(auto itn = std::next(it); itn != module_list.end(); ++itn)
        {
            double dx, dy;
            qdist((*it)->GetX(), (*it)->GetY(), (*it)->GetSectorID(),
                  (*itn)->GetX(), (*itn)->GetY(), (*itn)->GetSectorID(),
                  sector_info, dx, dy);
            if(std::abs(dx) < 1.01 && std::abs(dy) < 1.01) {
                (*it)->AddNeighbor(*itn, dx, dy);
                (*itn)->AddNeighbor(*it, -dx, -dy);
            }
        }
    }
}

void PRadHyCalDetector::UpdateSectorInfo()
{
    int Ns = static_cast<int>(Max_Sector);
    sector_info.clear();
    sector_info.resize(Ns);
    std::vector<bool> init(Ns, false);
    std::vector<double> xmin(Ns), ymin(Ns), xmax(Ns), ymax(Ns);

    double x1, y1, z1, x2, y2, z2;
    for(auto &module : module_list)
    {
        int i = module->GetSectorID();
        if(i < 0 || i >= Ns)
            continue;

        module->GetBoundary(x1, y1, z1, x2, y2, z2);
        // update boundary
        if(init[i]) {
            xmin[i] = std::min(xmin[i], x1);
            ymin[i] = std::min(ymin[i], y1);
            xmax[i] = std::max(xmax[i], x2);
            ymax[i] = std::max(ymax[i], y2);
        // initialize
        } else {
            sector_info[i].Init(i, module);
            xmin[i] = x1;
            ymin[i] = y1;
            xmax[i] = x2;
            ymax[i] = y2;
            init[i] = true;
        }
    }

    // form boundary points
    for(int i = 0; i < Ns; ++i)
        sector_info[i].SetBoundary(xmin[i], ymin[i], xmax[i], ymax[i]);
}

// the distance quantized by the module size
// module size is dependent on the Moliere radius
double PRadHyCalDetector::QuantizedDist(const PRadHyCalModule *m1, const PRadHyCalModule *m2)
const
{
    double dx, dy;
    qdist(m1->GetX(), m1->GetY(), m1->GetSectorID(),
          m2->GetX(), m2->GetY(), m2->GetSectorID(),
          sector_info, dx, dy);
    return std::sqrt(dx*dx + dy*dy);
}

double PRadHyCalDetector::QuantizedDist(double x1, double x2, double y1, double y2)
const
{
    double dx, dy;
    qdist(x1, y1, GetSectorID(x1, y1), x2, y2, GetSectorID(x2, y2), sector_info, dx, dy);
    return std::sqrt(dx*dx + dy*dy);
}

double PRadHyCalDetector::QuantizedDist(double x1, double y1, int s1,
                                        double x2, double y2, int s2)
const
{
    double dx, dy;
    qdist(x1, y1, s1, x2, y2, s2, sector_info, dx, dy);
    return std::sqrt(dx*dx + dy*dy);
}

void PRadHyCalDetector::QuantizedDist(const PRadHyCalModule *m1, const PRadHyCalModule *m2,
                                      double &dx, double &dy)
const
{
    qdist(m1->GetX(), m1->GetY(), m1->GetSectorID(),
          m2->GetX(), m2->GetY(), m2->GetSectorID(),
          sector_info, dx, dy);
}

void PRadHyCalDetector::QuantizedDist(double x1, double x2, double y1, double y2,
                                      double &dx, double &dy)
const
{
    qdist(x1, y1, GetSectorID(x1, y1), x2, y2, GetSectorID(x2, y2), sector_info, dx, dy);
}

void PRadHyCalDetector::QuantizedDist(double x1, double y1, int s1,
                                      double x2, double y2, int s2,
                                      double &dx, double &dy)
const
{
    qdist(x1, y1, s1, x2, y2, s2, sector_info, dx, dy);
}

// using primex id to get layout information
// TODO now it is highly specific to the current HyCal layout, make it configurable
void PRadHyCalDetector::setLayout(PRadHyCalModule &module)
const
{
    int pid = module.GetID(), sector = -1, col = 0, row = 0;
    unsigned int flag;

    // calculate geometry information
    flag = 0;
    if(pid > PWO_ID0) {
        // crystal module
        pid -= PWO_ID0 + 1;
        sector = (int)Center;
        row = pid/34 + 1;
        col = pid%34 + 1;

        // set flag
        SET_BIT(flag, kPbWO4);
        if(row <= 19 && row >= 16 && col <= 19 && row >= 16)
            SET_BIT(flag, kInnerBound);
        if(row == 1 || row == 34 || col == 1 || col == 34)
            SET_BIT(flag, kTransition);

    } else {
        // lead glass module
        pid -= 1;
        int g_row = pid/30 + 1;
        int g_col = pid%30 + 1;

        // set flag
        SET_BIT(flag, kPbGlass);
        if(g_row == 1 || g_row == 30 || g_col == 1 || g_col == 30)
            SET_BIT(flag, kOuterBound);

        // there are 4 sectors for lead glass
        // top sector
        if(g_col <= 24 && g_row <= 6) {
            sector = (int)Top;
            row = g_row;
            col = g_col;
            if(row == 6 && col >= 6)
                SET_BIT(flag, kTransition);
        }
        // right sector
        if(g_col > 24 && g_row <= 24) {
            sector = (int)Right;
            row = g_row;
            col = g_col - 24;
            if(col == 1 && row >= 6)
                SET_BIT(flag, kTransition);
        }
        // bottom sector
        if(g_col > 6 && g_row > 24) {
            sector = (int)Bottom;
            row = g_row - 24;
            col = g_col - 6;
            if(row == 1 && col < 20)
                SET_BIT(flag, kTransition);
        }
        // left sector
        if(g_col <= 6 && g_row > 6) {
            sector = (int)Left;
            row = g_row - 6;
            col = g_col;
            if(col == 6 && row < 20)
                SET_BIT(flag, kTransition);
        }
    }

    module.SetLayout(Layout(flag, sector, row-1, col-1));
}

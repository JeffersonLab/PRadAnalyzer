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


// enum name lists
static const char *__hycal_sector_list[] = {"Center", "Top", "Right", "Bottom", "Left"};



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
: PRadDetector(that), system(nullptr), module_hits(that.module_hits),
  dead_hits(that.dead_hits), module_clusters(that.module_clusters),
  hycal_hits(that.hycal_hits)
{
    for(auto module : that.module_list)
    {
        AddModule(new PRadHyCalModule(*module));
    }
}

// move constructor
PRadHyCalDetector::PRadHyCalDetector(PRadHyCalDetector &&that)
: PRadDetector(that), system(nullptr), module_list(std::move(that.module_list)),
  id_map(std::move(that.id_map)), name_map(std::move(that.name_map)),
  module_hits(std::move(that.module_hits)), dead_hits(std::move(that.dead_hits)),
  module_clusters(std::move(that.module_clusters)), hycal_hits(std::move(that.hycal_hits))
{
    // reset the connections between module and HyCal
    for(auto module : module_list)
        module->SetDetector(this);
}

// destructor
PRadHyCalDetector::~PRadHyCalDetector()
{
    UnsetSystem();

    // release modules
    for(auto module : module_list)
    {
        // prevent module calling RemoveModule upon destruction
        module->UnsetDetector(true);
        delete module;
    }
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
    id_map = std::move(rhs.id_map);
    name_map = std::move(rhs.name_map);
    module_hits = std::move(rhs.module_hits);
    dead_hits = std::move(rhs.dead_hits);
    module_clusters = std::move(rhs.module_clusters);
    hycal_hits = std::move(rhs.hycal_hits);

    for(auto module : module_list)
        module->SetDetector(this);

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

        geo.type = PRadHyCalModule::get_module_type(type.c_str());

        PRadHyCalModule *module = new PRadHyCalModule(name, geo);

        // failed to add module to detector
        if(!AddModule(module))
            delete module;
    }

    // sort the module by id
    SortModuleList();

    // update sector information
    UpdateSectorInfo();

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

// prepare dead hits for the leakage correction in reconstruction
void PRadHyCalDetector::CreateDeadHits()
{
    // clear current dead hits
    dead_hits.clear();

    // create dead hits
    for(auto module : module_list)
    {
        // clear the dead module flag
        CLEAR_BIT(module->layout.flag, kDeadModule);

        // module is not connected to a adc channel or the channel is dead
        if(!module->GetChannel() || module->GetChannel()->IsDead()) {
            SET_BIT(module->layout.flag, kDeadModule);
            dead_hits.emplace_back(module->GetID(),         // id
                                   module->GetGeometry(),   // geometry
                                   module->GetLayout(),     // layout
                                   0.,                      // energy
                                   false);                  // virtual
        }
    }

    // check if any modules are very close to this dead module, and set a bit
    // for the future correction
    for(auto module : module_list)
    {
        const auto &geo = module->GetGeometry();

        for(auto &dead : dead_hits)
        {
            float dx = (geo.x - dead.geo.x)/(geo.size_x + dead.geo.size_x);
            float dy = (geo.y - dead.geo.y)/(geo.size_y + dead.geo.size_y);

            if(sqrt(dx*dx + dy*dy)*2. < CORNER_ADJACENT)
            {
               SET_BIT(module->layout.flag, kDeadNeighbor);
               break;
            }
        }
    }
}

// hits/clusters reconstruction
void PRadHyCalDetector::Reconstruct(PRadHyCalCluster *method)
{
    // clear containers
    hycal_hits.clear();

    // group module hits into clusters
    method->FormCluster(module_hits, module_clusters);

    for(auto &cluster : module_clusters)
    {
        // discard cluster that does not satisfy certain conditions
        if(!method->CheckCluster(cluster))
            continue;

        // leakage correction for dead modules
        method->LeakCorr(cluster, dead_hits);

        // the center module does not exist should be a fatal problem, thus no
        // safety check here
        PRadHyCalModule *center = GetModule(cluster.center.id);

        // get non-linear correction factor
        float lin_corr = center->GetCalibConst().NonLinearCorr(cluster.energy);

        // reconstruct hit the position based on the cluster
        HyCalHit hit = method->Reconstruct(cluster, lin_corr);

        // add timing information
        PRadTDCChannel *tdc = center->GetTDC();
        if(tdc)
            hit.set_time(tdc->GetTimeMeasure());

        // final hit reconstructed
        hycal_hits.emplace_back(std::move(hit));
    }
}

// collect hits from modules
void PRadHyCalDetector::CollectHits()
{
    module_hits.clear();

    for(auto &module : module_list)
    {
        float energy = module->GetEnergy();
        if(energy > 0)
            module_hits.emplace_back(module->GetID(),           // id
                                     module->GetGeometry(),     // geometry
                                     module->GetLayout(),       // layout
                                     energy);                   // energy
    }
}

// clear existing hits
void PRadHyCalDetector::ClearHits()
{
    module_hits.clear();
}

PRadHyCalModule *PRadHyCalDetector::GetModule(const int &id)
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

PRadHyCalModule *PRadHyCalDetector::GetModule(const float &x, const float &y)
const
{
    for(auto &module : module_list)
    {
        float pos_x = module->GetX();
        float size_x = module->GetSizeX();
        if((x > pos_x + size_x/2.) || (x < pos_x - size_x/2.))
            continue;
        float pos_y = module->GetY();
        float size_y = module->GetSizeY();
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

        module->GetBoundaries(x1, y1, z1, x2, y2, z2);
        // update boundary
        if(init[i]) {
            xmin[i] = std::min(xmin[i], x1);
            ymin[i] = std::min(ymin[i], y1);
            xmax[i] = std::max(xmax[i], x2);
            ymax[i] = std::max(ymax[i], y2);
        // initialize
        } else {
            xmin[i] = x1;
            ymin[i] = y1;
            xmax[i] = x2;
            ymax[i] = y2;
            init[i] = true;
            sector_info[i].msize_x = module->GetSizeX();
            sector_info[i].msize_y = module->GetSizeY();
        }
    }

    // form boundary points
    for(int i = 0; i < Ns; ++i)
    {
        sector_info[i].boundpts.clear();
        // align in counter-clockwise
        sector_info[i].boundpts.emplace_back(xmin[i], ymax[i]);
        sector_info[i].boundpts.emplace_back(xmin[i], ymin[i]);
        sector_info[i].boundpts.emplace_back(xmax[i], ymin[i]);
        sector_info[i].boundpts.emplace_back(xmax[i], ymax[i]);
    }
}

// the distance quantized by the module size
// module size is dependent on the Moliere radius
double PRadHyCalDetector::QuantizedDist(PRadHyCalModule *m1, PRadHyCalModule *m2)
const
{
    double x1 = m1->GetX(), y1 = m1->GetY(), x2 = m2->GetX(), y2 = m2->GetY();
    double dx = 0., dy = 0.;
    // in the same sector
    if(m1->GetSectorID() == m2->GetSectorID()) {
        dx = (x1 - x2)/m1->GetSizeX();
        dy = (y1 - y2)/m1->GetSizeY();

    // in different sectors, and with different types of module
    // NOTICE highly specific for the current HyCal layout
    // the center sector is for crystal modules
    } else if (m1->GetType() != m2->GetType()) {
        auto &boundary = sector_info[static_cast<int>(Center)].boundpts;
        for(size_t i = 0; i < boundary.size(); ++i)
        {
            // boundary line from two points
            size_t ip = (i == 0) ? boundary.size() - 1 : i - 1;
            auto &p1 = boundary[ip], &p2 = boundary[i];

            double xc = 0., yc = 0.;
            int inter = cana::intersection(p1.x, p1.y, p2.x, p2.y, x1, y1, x2, y2, xc, yc);

            if(inter == 0) {
                dx = (x1 - xc)/m1->GetSizeX() + (xc - x2)/m2->GetSizeX();
                dy = (y1 - yc)/m1->GetSizeX() + (yc - y2)/m2->GetSizeY();
                // there will be only one boundary satisfies all the conditions
                break;
            }
        }
    // in different sectors but with the same type of module
    } else {
        // 2 possibilities, pass through the center part or not
        auto &spwo = sector_info[static_cast<int>(Center)];

        double xc[2], yc[2];
        int ic = 0;
        for(size_t i = 0; i < spwo.boundpts.size(); ++i)
        {
            size_t ip = (i == 0) ? spwo.boundpts.size() - 1 : i - 1;
            auto &p1 = spwo.boundpts[ip], &p2 = spwo.boundpts[i];
            int inter = cana::intersection(p1.x, p1.y, p2.x, p2.y, x1, y1, x2, y2, xc[ic], yc[ic]);

            // find two points
            if(inter == 0 && ic++ > 0) break;
        }

        if(ic > 1) {
            double dxc = std::abs(xc[0] - xc[1]);
            double dyc = std::abs(yc[0] - yc[1]);
            double dxt = std::abs(x1 - x2) - dxc;
            double dyt = std::abs(y1 - y2) - dyc;
            dx = dxt/m1->GetSizeX() + dxc/spwo.msize_x;
            dy = dyt/m1->GetSizeY() + dyc/spwo.msize_y;
        } else {
            dx = (x1 - x2)/(m1->GetSizeX());
            dy = (y1 - y2)/(m1->GetSizeY());
        }
    }

    return sqrt(dx*dx + dy*dy);
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

// quantize the distance between two modules by there sizes
// by this way we can indiscriminately check modules with different size
// only useful for adjacent module checking
float PRadHyCalDetector::hit_distance(const ModuleHit &m1, const ModuleHit &m2)
{
    float dx = (m1.geo.x - m2.geo.x)/(m1.geo.size_x + m2.geo.size_x);
    float dy = (m1.geo.y - m2.geo.y)/(m1.geo.size_y + m2.geo.size_y);

    return sqrt(dx*dx + dy*dy)*2.;
}

// get enum HyCalSector by its name
int PRadHyCalDetector::get_sector_id(const char *name)
{
    for(int i = 0; i < (int)Max_Sector; ++i)
        if(strcmp(name, __hycal_sector_list[i]) == 0)
            return i;

    std::cerr << "PRad HyCal Detector Error: Cannot find sector " << name
              << ", please check the definition in PRadHyCalModule."
              << std::endl;
    // not found
    return -1;
}

// get name of HyCalSector
const char *PRadHyCalDetector::get_sector_name(int sec)
{
    if(sec < 0 || sec >= (int)Max_Sector)
        return "Undefined";
    else
        return __hycal_sector_list[sec];
}

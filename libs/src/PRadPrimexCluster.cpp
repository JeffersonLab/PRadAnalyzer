//============================================================================//
// A C++ wrapper for island reconstruction method from PrimEx code            //
//                                                                            //
// Weizhi Xiong, Chao Peng                                                    //
// 09/28/2016                                                                 //
//============================================================================//

#include "PRadPrimexCluster.h"




#define ICH(M,N) __prcl_ich[N-1][M-1]
// a temporary storage to convert row and col back to id
// which is needed in fetching result from island.F
static int __prcl_ich[MROW][MCOL];

PRadPrimexCluster::PRadPrimexCluster(const std::string &path)
{
    // configuration
    Configure(path);
}

PRadPrimexCluster::~PRadPrimexCluster()
{
    // place holder
}

PRadHyCalCluster *PRadPrimexCluster::Clone()
const
{
    return new PRadPrimexCluster(*this);
}

void PRadPrimexCluster::Configure(const std::string &path)
{
    PRadHyCalCluster::Configure(path);

    bool verbose = !path.empty();

    // adj_dist is used for merging clusters separated by sectors
    bool corner = getDefConfig<bool>("Corner Connection", false, verbose);
    if(corner)
        adj_dist = CORNER_ADJACENT;
    else
        adj_dist = SIDE_ADJACENT;

    // set the min module energy for all the module type
    float univ_min_energy = getDefConfig<float>("Min Module Energy", 0., false);
    min_module_energy.resize(PRadHyCalModule::Max_Type, univ_min_energy);
    // update the min module energy if some type is specified
    // the key is "Min Module Energy [typename]"
    for(unsigned int i = 0; i < min_module_energy.size(); ++i)
    {
        // determine key name
        std::string type = PRadHyCalModule::get_module_type_name(i);
        std::string key = "Min Module Energy [" + type + "]";
        auto value = GetConfigValue(key);
        if(!value.IsEmpty())
            min_module_energy[i] = value.Float();
    }

    // pass parameters to fortran code
    // MeV to GeV
    SET_EMIN  = min_cluster_energy*0.001;   // banks->CONFIG->config->CLUSTER_ENERGY_MIN;
    SET_EMAX  = 9.9;                        // banks->CONFIG->config->CLUSTER_ENERGY_MAX;
    SET_HMIN  = min_cluster_size;           // banks->CONFIG->config->CLUSTER_MIN_HITS_NUMBER;
    SET_MINM  = min_center_energy*0.001;    // banks->CONFIG->config->CLUSTER_MAX_CELL_MIN_ENERGY;
}

void PRadPrimexCluster::LoadCrystalProfile(const std::string &path)
{
    if(path.empty())
        return;

    char c_path[path.size()];
    strcpy(c_path, path.c_str());
    load_pwo_prof_(c_path, strlen(c_path));
}

void PRadPrimexCluster::LoadLeadGlassProfile(const std::string &path)
{
    if(path.empty())
        return;

    char c_path[path.size()];
    strcpy(c_path, path.c_str());
    load_lg_prof_(c_path, strlen(c_path));
}

void PRadPrimexCluster::UpdateModuleStatus(const std::vector<PRadHyCalModule*> &mlist)
{
    // initialize module status table
    for(int k = 0; k < MSECT; ++k)
        for(int i = 0; i < MCOL; ++i)
            for(int j = 0; j < MROW; ++j)
                module_status[k][i][j] = -1;

    for(auto &module : mlist)
    {
        int sect = module->GetSectorID();
        int col = module->GetColumn();
        int row = module->GetRow();

        if(TEST_BIT(module->GetLayoutFlag(), kDeadModule))
            module_status[sect][col][row] = 1;
        else
            module_status[sect][col][row] = 0;

    }
}

void PRadPrimexCluster::FormCluster(std::vector<ModuleHit> &hits,
                                    std::vector<ModuleCluster> &clusters)
const
{
    // clear container first
    clusters.clear();

    // build a hit map
    std::map<int, ModuleHit*> hit_map;
    for(auto &hit : hits)
    {
        hit_map[hit.id] = &hit;
    }

    // call island reconstruction of each sectors
    // HyCal has 5 sectors, 4 for lead glass one for crystal
    std::vector<std::vector<ModuleCluster>> sect_clusters;
    sect_clusters.resize(MSECT);
    for(int isect = 0; isect < MSECT; ++isect)
    {
        callIsland(hits, isect);
        sect_clusters[isect] = getIslandResult(hit_map);
    }

    // glue clusters separated by the sector
    for(size_t i = 0; i < MSECT; ++i)
    {
        for(int j = i + 1; j < MSECT; ++j)
        {
            glueClusters(sect_clusters.at(i), sect_clusters.at(j));
        }
    }

    // add clusters to the container
    for(auto &sect : sect_clusters)
    {
        for(auto &cluster : sect)
        {
            clusters.emplace_back(std::move(cluster));
        }
    }
}

void PRadPrimexCluster::callIsland(const std::vector<ModuleHit> &hits, int isect)
const
{
    // determine sector's column and row ranges
    ISECT = isect;
    int coloffset = 0, rowoffset = 0;
    switch(isect)
    {
    case 0:
        NCOL = 34; NROW = 34;
        SET_XSIZE = CRYS_SIZE_X; SET_YSIZE = CRYS_SIZE_Y;
        break;
    case 1:
        NCOL = 24; NROW =  6;
        SET_XSIZE = GLASS_SIZE; SET_YSIZE = GLASS_SIZE;
        coloffset = 0, rowoffset = 0;
        break;
    case 2:
        NCOL =  6; NROW =  24;
        SET_XSIZE = GLASS_SIZE; SET_YSIZE = GLASS_SIZE;
        coloffset = 24, rowoffset = 0;
        break;
    case 3:
        NCOL = 24; NROW =  6;
        SET_XSIZE = GLASS_SIZE; SET_YSIZE = GLASS_SIZE;
        coloffset =  6, rowoffset = 24;
        break;
    case 4:
        NCOL =  6; NROW = 24;
        SET_XSIZE = GLASS_SIZE; SET_YSIZE = GLASS_SIZE;
        coloffset = 0, rowoffset = 6;
        break;
    default:
        printf("call_island bad sector given : %i\n",isect);
        exit(1);
    }

    // load module status and reset energies
    for(int icol = 1; icol <= NCOL; ++icol)
    {
        for(int irow = 1; irow <= NROW; ++irow)
        {
            ECH(icol,irow) = 0;
            ICH(icol,irow) = -1;
            STAT_CH(icol,irow) = module_status[isect][icol-1][irow-1];
        }
    }

    // load hits in this sector
    for(auto &hit : hits)
    {
        // not belong to this sector or energy is too low
        if((hit.sector != isect) ||
           (hit.energy < min_module_energy.at(hit.geo.type)))
            continue;

        const int &id  = hit.id;
        int column, row;
        if(id > 1000) {
            column = (id-1001)%NCOL+1;
            row    = (id-1001)/NROW+1;
        } else {
            column = (id-1)%(NCOL+NROW)+1-coloffset;
            row    = (id-1)/(NCOL+NROW)+1-rowoffset;
        }

        // discretize to 0.1 MeV
        ECH(column,row) = int(hit.energy*10. + 0.5);
        ICH(column,row) = id;
    }

    // call fortran code to reconstruct clusters
    main_island_();
}

// get result from fortran island code
std::vector<ModuleCluster> PRadPrimexCluster::getIslandResult(const std::map<int, ModuleHit*> &hmap)
const
{
    std::vector<ModuleCluster> res;
    res.reserve(10);

    for(int k = 0; k < adcgam_cbk_.nadcgam; ++k)
    {
        int ncl = std::min(adcgam_cbk_.u.iadcgam[k][8], MAX_CC);
        ModuleCluster cluster;
        // GeV to MeV
        cluster.energy = adcgam_cbk_.u.fadcgam[k][0]*1000.;
        cluster.leakage = cluster.energy;

        for(int j = 0; j < ncl; ++j)
        {
            int add = ICL_INDEX(k,j);
            int kx = (add/100), ky = add%100;
            int id = ICH(kx, ky);
            // convert back from 0.1 MeV
            float ecell = 0.1*(float)ICL_IENER(k,j);
            cluster.leakage -= ecell;
            auto it = hmap.find(id);
            if(it != hmap.end())
            {
                ModuleHit hit(*it->second);
                hit.energy = ecell;
                cluster.hits.push_back(hit);
            }
        }
        cluster.FindCenter();
        res.push_back(cluster);
    }

    return res;
}

void PRadPrimexCluster::glueClusters(std::vector<ModuleCluster> &base,
                                     std::vector<ModuleCluster> &sector)
const
{
    // merge clusters between sectors
    for(auto it = sector.begin(); it != sector.end(); ++it)
    {
        for(auto it2 = base.begin(); it2 != base.end(); ++it2)
        {
            if(checkTransAdj(*it, *it2)) {
                it2->Merge(*it);
                sector.erase(it--);
                break;
            }
        }
    }
}

inline bool PRadPrimexCluster::checkTransAdj(const ModuleCluster &c1,
                                             const ModuleCluster &c2)
const
{
    // don't merge the clusters in the same sector
    if(c1.center.sector == c2.center.sector)
        return false;

    for(auto &m1 : c1.hits)
    {
        for(auto &m2 : c2.hits)
        {
            if(PRadHyCalDetector::hit_distance(m1, m2) < adj_dist) {
                return true;
            }
        }
    }

    return false;
}

void PRadPrimexCluster::LeakCorr(ModuleCluster &, const std::vector<ModuleHit> &)
const
{
    // place holder
    // TODO island.F has corrected the leakage
    // here we probably can add some inforamtion about the leakage correction
}

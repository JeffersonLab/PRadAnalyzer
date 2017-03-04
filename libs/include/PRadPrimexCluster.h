#ifndef PRAD_PRIMEX_CLUSTER_H
#define PRAD_PRIMEX_CLUSTER_H

#include <map>
#include <string>
#include <vector>
#include "PRadHyCalCluster.h"

//this is a c++ wrapper around the primex island algorithm
//used for HyCal cluster reconstruction

//define some global constants
#define MAX_CC 60
#define MSECT 5
#define MCOL 34
#define MROW 34

#define CRYS_SIZE_X 2.077   // real X-size of crystal
#define CRYS_SIZE_Y 2.075   // real Y-size of crystal
#define GLASS_SIZE 3.815    // real size of glass

extern "C"
{
    void load_pwo_prof_(char* config_dir, int str_len);
    void load_lg_prof_(char* config_dir, int str_len);
    void main_island_();

    extern struct
    {
        int ech[MROW][MCOL];
    } ech_common_;

    extern struct
    {
        int stat_ch[MROW][MCOL];
    } stat_ch_common_;

    extern struct
    {
        int icl_index[MAX_CC][200], icl_iener[MAX_CC][200];
    } icl_common_;

    extern struct
    {
        float xsize, ysize, mine, maxe;
        int min_dime;
        float minm;
        int ncol, nrow;
        float zhycal;
        int isect;
    } set_common_;

    extern struct
    {
        int nadcgam;
        union
        {
            int iadcgam[50][11];
            float fadcgam[50][11];
        } u;
    } adcgam_cbk_;

    extern struct
    {
        float fa[100];
    } hbk_common_;

    #define ECH(M,N) ech_common_.ech[N-1][M-1]
    #define STAT_CH(M,N) stat_ch_common_.stat_ch[N-1][M-1]
    #define ICL_INDEX(M,N) icl_common_.icl_index[N][M]
    #define ICL_IENER(M,N) icl_common_.icl_iener[N][M]
    #define HEGEN(N) read_mcfile_com_.hegen[N-1]
    #define SET_XSIZE set_common_.xsize
    #define SET_YSIZE set_common_.ysize
    #define SET_EMIN  set_common_.mine
    #define SET_EMAX  set_common_.maxe
    #define SET_HMIN  set_common_.min_dime
    #define SET_MINM  set_common_.minm
    #define NCOL      set_common_.ncol
    #define NROW      set_common_.nrow
    #define ZHYCAL    set_common_.zhycal
    #define ISECT     set_common_.isect
    #define FA(N) hbk_common_.fa[N-1]
}

class PRadPrimexCluster : public PRadHyCalCluster
{
public:
    PRadPrimexCluster(const std::string &path = "");
    virtual ~PRadPrimexCluster();
    PRadHyCalCluster *Clone() const;

    void Configure(const std::string &path);
    void LoadCrystalProfile(const std::string &path);
    void LoadLeadGlassProfile(const std::string &path);
    void UpdateModuleStatus(const std::vector<PRadHyCalModule*> &mlist);
    void FormCluster(std::vector<ModuleHit> &hits,
                     std::vector<ModuleCluster> &clusters) const;
    void LeakCorr(ModuleCluster &c, const std::vector<ModuleHit> &dead) const;

private:
    void callIsland(const std::vector<ModuleHit> &hits, int isect) const;
    std::vector<ModuleCluster> getIslandResult(const std::map<int, ModuleHit*> &hmap) const;
    void glueClusters(std::vector<ModuleCluster> &b, std::vector<ModuleCluster> &s) const;
    bool checkTransAdj(const ModuleCluster &c1, const ModuleCluster &c2) const;

private:
    float adj_dist;
    std::vector<float> min_module_energy;
    int module_status[MSECT][MCOL][MROW];
};

#endif

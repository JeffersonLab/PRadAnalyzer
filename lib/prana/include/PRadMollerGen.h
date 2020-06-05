#ifndef PRAD_MOLLER_GEN_H
#define PRAD_MOLLER_GEN_H

#include <vector>
#include "canalib.h"

#define MAX_T_BINS 30000
#define MAX_V_BINS 30000

extern "C"
{
    // interface to visit meradgen
    void merad_init(double Elab);
    double merad_sig(double t, double pl, int type);
    double merad_sig2(double t, int type);
    double merad_sigir2(double vmin, double t);
    double merad_sigir(double vmin, double t, double pl);
    double merad_sigfs(double v, double t, double pl, double sig0);
    double merad_sigfh(double v, double t, double pl);
    double merad_sample_t1(double t, double vgen, double plin, double rnd);
    double merad_sample_z(double t, double t1gen, double vgen, double plin, double rnd);

    // the grid used in meradgen declared in
    // fortran/include/merad_grid.inc
    // dimension must be changed accordingly
#define MERAD_NT1 30
#define MERAD_NZ 60

    extern struct
    {
        double grt1[MERAD_NT1], grz[MERAD_NZ];
    } merad_grid_;

    extern struct
    {
        double distsit1[4*MERAD_NT1], distart1[4*MERAD_NT1];
        double distsiz[MERAD_NZ], distarz[MERAD_NZ];
    } merad_dist_;
};

// v distribution structure
struct VDist : public cana::val_cdf
{
    double sig;

    VDist(double v, double s)
    : val_cdf(v, 0.), sig(s)
    {}
};

// t(Q^2) distribution structure
struct TDist : public cana::val_cdf
{
    // point information
    double sig_nrad, sig_rad;
    // v related
    std::vector<VDist> v_dist;

    // constructor
    TDist(double t, double nrad, std::vector<VDist> &&dist)
    : val_cdf(t, 0.)
    {
        sig_nrad = nrad;
        v_dist = dist;
        sig_rad = dist.back().cdf;
    }
};

// unit MeV, degree and nb
class PRadMollerGen
{
public:
    PRadMollerGen(double vmin = 5, double vmax = 400, int mbins = 100,
                  double t_res = 1e-6, double v_res = 1e-4);
    virtual ~PRadMollerGen();

    double Generate(double Es, double min_angle, double max_angle, int nevents,
                    const char *path, bool verbose = true) const;
    void GetXS(double Es, double angle,
               double &sig_born, double &sig_nrad, double &sig_rad) const;
    void GetXSdQsq(double s, double t,
                   double &sig_born, double &sig_nrad, double &sig_rad) const;

    double GetNonRadXSdQsq(double s, double t) const;
    std::vector<VDist> GetRadVDist(double s, double t) const;

    void SetMinBins(unsigned int mbins) {min_bins = mbins;}
    void SetPrecision(double t_res, double v_res) {t_prec = t_res; v_prec = v_res;}
    unsigned int GetMinBins() const {return min_bins;}
    double GetTDistPrecision() const {return t_prec;}
    double GetVDistPrecision() const {return v_prec;}

    // static functions
    static double SigmaBorn(double s, double t);
    static void SigmaVph(double s, double t,
                         double &sig_0, double &sig_S, double &sig_vert, double &sig_B);
    static double SigmaIR(double s, double t, double v_max);
    static double SigmaFs(double s, double t, double v_min, double res = 1e-3);
    static double SigmaRad(double s, double t, double v_min, double v_max, double res = 1e-4);
    static bool MomentumRec(double *k2, double *p2, double *k,
                            double s, double t, double t1, double v, double z,
                            double rnd1 = 0., double rnd2 = 0.);

private:
    std::vector<TDist> init_grids(double s, double t_min, double t_max, bool verbose) const;

private:
    // v_min defines the minimum photon energy that to be generated (hard photons)
    // v_cut defines the integration range to the highest photon energy
    double v_min, v_cut;
    // minimum theta bins
    unsigned int min_bins;
    // required theta interpolation precision
    double t_prec, v_prec;
};

#endif


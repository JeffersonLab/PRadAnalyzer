#ifndef PRAD_EP_ELAS_GEN_H
#define PRAD_EP_ELAS_GEN_H

#include <vector>
#include "canalib.h"

// unit MeV, degree and nb
class PRadEpElasGen
{
public:
    PRadEpElasGen(double vmin = 5, double vmax = 400, int mbins = 100,
                  double t_res = 1e-6, double v_res = 1e-4);
    virtual ~PRadEpElasGen();


    void GetXSdQsq(double S, double Q2,
                   double &sig_born, double &sig_nrad, double &sig_rad) const;
    void GetEMFF(double Q2, double &GE, double &GM) const;
    void GetHadStrFunc(double Q2, double &F1, double &F2) const;

    void CalcVphIR(double S, double Q2, double v_max,
                   double &sig_born, double &sig_AMM,
                   double &delta_VR, double &delta_vac, double &delta_inf) const;
    double SigmaBorn(double S, double Q2) const;
    double SigmaBrem(double v, double t, double phik, double S, double Q2) const;
    double SigmaBrem_phik(double v, double t, double S, double Q2, bool finite) const;
    double SigmaBrem_phik_v(double t, double v_min, double v_max, double S, double Q2, bool finite) const;
    double SigmaFh(double t1, double t2, double v1, double v2, double S, double Q2) const;
    double SigmaFs(double t1, double t2, double v1, double v2, double S, double Q2) const;

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


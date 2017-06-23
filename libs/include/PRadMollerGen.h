#ifndef PRAD_MOLLER_GEN_H
#define PRAD_MOLLER_GEN_H

extern "C"
{
    // interface to visit meradgen
    void merad_init(double Elab);
    double merad_sig(double t, double pl, int type);
    double merad_sig2(double t, int type);
    double merad_sigir2(double vmin, double t);
    double merad_sigir(double vmin, double t, double pl);
    double merad_sigfs(double vmin, double t, double pl);
    double merad_sigfh(double vmin, double vmax, double t, double pl);
    double merad_sample_t1(double t, double vgen, double plin, double rnd);
    double merad_sample_z(double t, double t1gen, double vgen, double plin, double rnd);

    // the grid used in meradgen declared in
    // fortran/include/merad_grid.inc
    // dimension must be changed accordingly
#define MERAD_NV 3000
#define MERAD_NT1 30
#define MERAD_NZ 60

    extern struct
    {
        double grv[MERAD_NV], grt1[MERAD_NT1], grz[MERAD_NZ];
    } merad_grid_;

    extern struct
    {
        double distsiv[MERAD_NV], distarv[MERAD_NV];
        double distsit1[4*MERAD_NT1], distart1[4*MERAD_NT1];
        double distsiz[MERAD_NZ], distarz[MERAD_NZ];
    } merad_dist_;
};

// unit MeV, degree and nb
class PRadMollerGen
{
public:
    PRadMollerGen(double vmin = 5, double vmax = 400, int mbins = 100, double prec = 1e-6);
    virtual ~PRadMollerGen();

    double Generate(double Es, double min_angle, double max_angle, int nevents,
                    const char *path, bool verbose = true) const;
    void GetXS(double Es, double angle,
               double &sig_born, double &sig_nrad, double &sig_rad) const;
    void GetXSdQsq(double s, double t,
                   double &sig_born, double &sig_nrad, double &sig_rad) const;

    void SetMinBins(unsigned int mbins) {min_bins = mbins;}
    void SetPrecision(double prec) {req_prec = prec;}
    unsigned int GetMinBins() const {return min_bins;}
    double GetPrecision() const {return req_prec;}

    // static functions
    static void SigmaVph(double s, double t,
                         double &sig_0, double &sig_S, double &sig_vert, double &sig_B);
    static void SigmaIR(double s, double t, double v_max,
                        double &delta_1H, double &delta_1S, double &delta_1inf);
    static void SigmaF(double s, double t, double v_min, double v_max,
                       double &sig_Fs, double &sig_Fh);
    static void MomentumRec(double *k2, double *p2, double *k,
                            double s, double t, double t1, double v, double z,
                            double rnd1 = 0., double rnd2 = 0.);

private:
    // v_min defines the minimum photon energy that to be generated (hard photons)
    // v_cut defines the integration range to the highest photon energy
    double v_min, v_cut;
    // minimum theta bins
    unsigned int min_bins;
    // required theta interpolation precision
    double req_prec;
};

#endif


#ifndef PRAD_MOLLER_GEN_H
#define PRAD_MOLLER_GEN_H

extern "C"
{
    // interface to visit meradgen
    void merad_init(double Elab);
    double merad_sigfs(double vmin, double t, double pl);
    double merad_sigfh(double vmin, double vmax, double t, double pl);

    // the grid used in meradgen declared in
    // fortran/include/merad_grid.inc
    // dimension must be changed accordingly
#define MERAD_NV 60
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
    PRadMollerGen(double vmin = 5, double vmax = 400);
    virtual ~PRadMollerGen();

    void Generate(double Es, double min_angle, double max_angle, int nevents) const;
    void GetXS(double Es, double angle, double &sig_born, double &sig_nrad, double &sig_rad) const;

private:
    void moller_Vph(double s, double t, double u0,
                    double &sig_0, double &sig_S, double &sig_vert, double &sig_B) const;
    void moller_IR(double s, double t, double u0,
                   double &delta_1H, double &delta_1S, double &delta_1inf) const;

private:
    // v_min defines the minimum photon energy that to be generated (hard photons)
    // v_cut defines the integration range to the highest photon energy
    double v_min, v_cut;
};

#endif


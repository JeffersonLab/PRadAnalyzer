#ifndef PRAD_MOLLER_GEN_H
#define PRAD_MOLLER_GEN_H

extern "C"
{
    // interface to visit meradgen
    void merad_init(double Elab);
    double merad_sigfs(double vmin, double t, double pl, double born);
    double merad_sigfh(double vmin, double vmax, double t, double pl);

    // the grid used in meradgen declared in
    // fortran/include/merad_grid.inc
    // dimension must be changed accordingly
    extern struct
    {
        double grv[60], grt1[30], grz[60];
    } merad_grid_;
};

// unit MeV, degree and nb
class PRadMollerGen
{
public:
    PRadMollerGen(double vmin = 1, double vmax = 400);
    virtual ~PRadMollerGen();

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


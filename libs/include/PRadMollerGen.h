#ifndef PRAD_MOLLER_GEN_H
#define PRAD_MOLLER_GEN_H

extern "C"
{
    void merad_init(double Elab);
    double merad_fsir(double t, double t1, double v, double z, double sig0, int *nn, int ikey);
};

// unit MeV, degree and nb
class PRadMollerGen
{
public:
    PRadMollerGen(double ph_min = 1e-8, double ph_cut = 10.);
    virtual ~PRadMollerGen();

    void GetXS(double Es, double angle, double &sig_born, double &sig_nrad, double &sig_rad) const;

private:
    void moller_Vph(double s, double t, double u0,
                    double &sig_0, double &sig_S, double &sig_vert, double &sig_B) const;
    void moller_IR(double s, double t, double u0,
                   double &delta_1H, double &delta_1S, double &delta_1inf) const;
    double merad_fsirv(double v, double t, double sig0) const;

private:
    // v_min defines the minimum photon energy that to be generated (hard photons)
    // v_cut defines the integration range to the highest photon energy
    double v_min, v_cut;
};

#endif


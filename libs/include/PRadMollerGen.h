#ifndef PRAD_MOLLER_GEN_H
#define PRAD_MOLLER_GEN_H

// unit MeV, degree and nb
class PRadMollerGen
{
public:
    PRadMollerGen(double ph_min = 0., double ph_cut = 10.);
    virtual ~PRadMollerGen();

    double GetBornXS(const double &Es, const double &angle);
    double GetNonRadXS(const double &Es, const double &angle);

private:
    double moller_nonrad(double *p1, double *k1, double *k2, int type = 1);

private:
    // v_min defines the minimum photon energy that to be generated (hard photons)
    // v_cut defines the integration range to the highest photon energy
    double v_min, v_cut;
};

#endif


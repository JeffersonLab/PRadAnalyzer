#ifndef PRAD_MOLLER_GEN_H
#define PRAD_MOLLER_GEN_H

class PRadMollerGen
{
public:
    PRadMollerGen();
    virtual ~PRadMollerGen();

    double GetBornXS(const double &Es, const double &angle);
    double GetNonRadXS(const double &Es, const double &angle);

private:
    double moller_nonrad(double *p1, double *k1, double *k2, int type = 1);

};

#endif


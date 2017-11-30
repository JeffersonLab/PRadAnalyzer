#ifndef PRAD_TRIGGER_CONST_H
#define PRAD_TRIGGER_CONST_H

#include <vector>
#include <iostream>
#include <initializer_list>

// we need 3 parameters to calculate the trigger efficiency
#define TRGEFF_NPAR 3

class PRadTriggerConst
{
public:
    PRadTriggerConst();
    PRadTriggerConst(std::initializer_list<double> pars);
    PRadTriggerConst(const double *pars);

    virtual ~PRadTriggerConst();

    void SetTriggerParams(const std::initializer_list<double> &pars);
    void SetTriggerParams(const double *pars);

    // energy needs to be in MeV
    double GetTriggerEfficiency(double energy) const;

private:
    double params[TRGEFF_NPAR];
};

#endif

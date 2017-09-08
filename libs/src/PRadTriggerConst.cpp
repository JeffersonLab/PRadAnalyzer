//============================================================================//
// Trigger constant class, store the trigger efficiency related parameters    //
//                                                                            //
// Chao Peng                                                                  //
// 09/08/2017                                                                 //
//============================================================================//

#include "PRadTriggerConst.h"
#include <cmath>

//============================================================================//
// Constructors, Destructor                                                   //
//============================================================================//

// constructors
PRadTriggerConst::PRadTriggerConst()
{
    for(auto &par : params)
        par = 0.;

    // make sure it returns 1 if not initialized
    params[0] = 1.0;
}

PRadTriggerConst::PRadTriggerConst(std::initializer_list<double> pars)
{
    SetTriggerParams(pars);
}

PRadTriggerConst::PRadTriggerConst(const double *pars)
{
    SetTriggerParams(pars);
}

// destructor
PRadTriggerConst::~PRadTriggerConst()
{
    // place holder
}

void PRadTriggerConst::SetTriggerParams(const std::initializer_list<double> &pars)
{
    if(pars.size() != TRGEFF_NPAR) {
        std::cerr << "Warning: wrong number of parameters to initialize class "
                  << "TriggerConst, received " << pars.size() << ", expected "
                  << TRGEFF_NPAR << std::endl;
    }

    size_t i = 0;
    for(auto &par : pars)
    {
        if(i++ < TRGEFF_NPAR)
            params[i] = par;
    }
}

void PRadTriggerConst::SetTriggerParams(const double *pars)
{
    for(size_t i = 0; i < TRGEFF_NPAR; ++i)
        params[i] = pars[i];
}

double PRadTriggerConst::GetTriggerEfficiency(double energy)
const
{
    return params[0]*(1. - std::exp(params[1]*energy + params[2]));
}


//============================================================================//
// Calibration constant class, store the calibration related parameters       //
// base_factor, base_gains, base_energy directly from the calibration runs    //
// Accept LMS gain from any runs to correct the calibration factor            //
//                                                                            //
// Chao Peng                                                                  //
// 11/12/2016                                                                 //
//============================================================================//

#include "PRadCalibConst.h"
#include <iomanip>
#include <cmath>

//============================================================================//
// Constructors, Destructor                                                   //
//============================================================================//

// constructors
PRadCalibConst::PRadCalibConst(int ref_num)
: factor(0.), base_factor(0.), base_energy(0.), non_linear(0.)
{
    base_gains.resize(ref_num, 0);
}

PRadCalibConst::PRadCalibConst(double f, double e, double nl, const std::vector<double> &g)
: factor(f), base_factor(f), base_energy(e), base_gains(g), non_linear(nl)
{
    // place holder
}

PRadCalibConst::PRadCalibConst(double f, double e, double nl, double *g, int num)
: factor(f), base_factor(f), base_energy(e), non_linear(nl)
{
    for(int i = 0; i < num; ++i)
        base_gains.push_back(g[i]);
}

// destructor
PRadCalibConst::~PRadCalibConst()
{
    // place holder
}

void PRadCalibConst::SetRefGain(double gain, int ref)
{
    if((size_t)ref >= base_gains.size()) {
        std::cerr << "PRadCalibConst Error: Cannot update gain for reference PMT "
                  << ref + 1
                  << ", only has "
                  << base_gains.size()
                  << " reference PMTs"
                  << std::endl;
        return;
    }

    base_gains[ref] = gain;
}

void PRadCalibConst::ClearRefGains()
{
    for(auto &gain : base_gains)
        gain = 0;
}

double PRadCalibConst::GetRefGain(int ref)
const
{
    if((size_t)ref >= base_gains.size())
        return 0.;

    return base_gains.at(ref);
}

void PRadCalibConst::GainCorrection(double gain, int ref)
{
    double base = GetRefGain(ref);

    if((gain > 0.) && (base > 0.)) {
        factor = base_factor * base/gain;
    }
}

// transfer adc value to energy
double PRadCalibConst::Calibration(const double &adc_val)
const
{
    if(adc_val > 0.)
        return adc_val*factor;

    return 0.;
}

// correct the non linear energy response in HyCal, unit is in MeV
double PRadCalibConst::NonLinearCorr(const double &E)
const
{
    return non_linear*(E - base_energy)/1000.;
}

float PRadCalibConst::NonLinearCorr(const float &E)
const
{
    return non_linear*(E - base_energy)/1000.;
}

std::ostream &operator <<(std::ostream &os, const PRadCalibConst &cal)
{
    os << std::setw(12) << cal.GetBaseConst()
       << std::setw(10) << cal.GetCalibEnergy()
       << std::setw(15) << cal.GetNonLinearFactor();
    for(auto gain : cal.GetRefGains())
        os << std::setw(10) << gain;

    return os;
}

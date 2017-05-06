//============================================================================//
// Calculate UNPOLARIZED Moller cross section with radiation effect           //
//                                                                            //
// Main Reference: Eur. Phys. J. A (2015) 51: 1                               //
// Radiative corrections beyond the ultra relativistic limit in unpolarized   //
// ep elastic and Moller scatterings for the PRad Experiment at Jefferson     //
// Laboratory                                                                 //
// I. Akushevich, H. Gao, A. Ilyichev, and M. Meziane                         //
//                                                                            //
// Code Developer: Chao Peng                                                  //
//============================================================================//

#include "PRadMollerGen.h"
#include "canalib.h"
#include <cmath>
#include <iostream>

// TODO define lamda for now, but this part will be cancelled finally
const double lamda = 1.;

// constructor
PRadMollerGen::PRadMollerGen()
{
    // place holder
}

// destructor
PRadMollerGen::~PRadMollerGen()
{
    // place holder
}

// square
inline double pow2(double val) {return val*val;}
// cubic
inline double pow3(double val) {return val*val*val;}
// power of 4
inline double pow4(double val) {double val2 = pow2(val); return val2*val2;}

// Non-rad or Born level cross section for Moller scattering
// input three arrays with fixed size - 4 elements each
// p1, initial target four-momentum
// k1, initial beam electron four-momentum
// k2, final beam electron four-momentum
// type, 0. born, others. non-rad
// output cross section with the same unit from input
double PRadMollerGen::moller_nonrad(double *p1, double *k1, double *k2, int type)
{
    // Mandelstam variables
    double u0, s, t;
    // s = (k1 + p1)^2
    s = -(pow2(k1[1] + p1[1]) + pow2(k1[2] + p1[2]) + pow2(k1[3] + p1[3]))
        + pow2(k1[0] + p1[0]);

    // t = (k1 - k2)^2
    t = -(pow2(k1[1] - k2[1]) + pow2(k1[2] - k2[2]) + pow2(k1[3] - k2[3]))
        + pow2(k1[0] - k2[0]);

    // u0 = (k2 - p1)^2
    u0 = -(pow2(k2[1] - p1[1]) + pow2(k2[2] - p1[2]) + pow2(k2[3] - p1[3]))
         + pow2(k2[0] - p1[0]);

    double s2t = s*s*t, st2 = s*t*t, s3 = pow3(s);
    double u02t = u0*u0*t, u0t2 = u0*t*t, u03 = pow3(u0), t3 = pow3(t);

    // frequently used variables
    double m = cana::ele_mass, m2 = m*m;
    double alp_pi = cana::alpha/cana::pi;
    double pi2 = cana::pi*cana::pi;
    double alp2 = cana::alpha*cana::alpha;
    double alp3 = alp2*cana::alpha;

    double xi_s = sqrt(1. - 4.*m2/s);
    double xi_t = sqrt(1 - 4.*m2/t);
    double xi_u0 = sqrt(1. - 4.*m2/u0);
    double xi_s2 = xi_s*xi_s, xi_s4 = xi_s2*xi_s2;
    double xi_t2 = xi_t*xi_t, xi_t4 = xi_t2*xi_t2;
    double xi_u02 = xi_u0*xi_u0, xi_u04 = xi_u02*xi_u02;

    // equation (49), Born Level
    double sig0 = u0*u0/xi_s2/4./s*(4.*xi_u04 - pow2(1. - xi_u02)*(2. + t/u0)) - s*s*xi_s4/u0;
    sig0 *= 2.*cana::pi*alp2/st2;

    if(type == 0)
        return sig0;

    // other frequently used variables
    // Q^2 related (t)
    double log_m = log(lamda/m);
    double Q2_m = -t + 2.*m2;
    double lamda_m = t*t - 4.*m2*t;
    double slamda_m = sqrt(lamda_m);
    double L_m = 1./slamda_m*log((slamda_m - t)/(slamda_m + t));

    // s related
    double log_s = log((1. + xi_s)/(1. - xi_s));
    double log_2s = log((1. + xi_s)/2./xi_s);
    double Li2_s = cana::spence((xi_s - 1.)/2./xi_s);
    double Li2_sp = cana::spence(2.*xi_s/(xi_s + 1.));

    // t related
    double log_t = log((1. + xi_t)/(xi_t - 1.));
    double log_2t = log((1. + xi_t)/2.);
    double log_4t = log((xi_t2 - 1.)/4.);
    double Li2_t = cana::spence((1. - xi_t)/2.);

    // u0 related
    double log_u0 = log((1. + xi_u0)/(xi_u0 - 1.));
    double log_2u0 = log((xi_u0 - 1.)/2./xi_u0);
    double log_2u0p = log((xi_u0 + 1.)/2./xi_u0);
    double Li2_u0 = cana::spence((1. + xi_u0)/2./xi_u0);


    // vacuum polarization for all leptons, factorized part
    // equation (41) with Q^2 -> -t
    double delta_vac = 0.;
    double lepton_mass[3] = {cana::ele_mass, cana::mu_mass, cana::tau_mass};
    for(auto &vac_m : lepton_mass)
    {
        double vac_m2 = vac_m*vac_m;
        double vac_slamda_m = sqrt(t*t - 4.*vac_m2*t);
        double vac_L_m = 1./vac_slamda_m*log((vac_slamda_m - t)/(vac_slamda_m + t));
        delta_vac += 2./3.*(-t + 2*vac_m2)*vac_L_m - 10./9. - 8./3.*vac_m2/t*(1. - 2.*vac_m2*vac_L_m);
    }

    // vertex correction, factorized part
    // equation (36) with Q^2 -> -t
    double delta_ver = 2.*(Q2_m*L_m - 1.)*log_m + (4.*m2 - 3./2.*t)*L_m - 2.
                       - Q2_m/slamda_m*(lamda_m*L_m*L_m/2. + 2.*cana::spence(2.*slamda_m/(slamda_m - t)) - pi2/2.);

    // vertex correction, non-factorized part, anomalous magnetic moment
    // euqation (52)
    double sig_AMM = 4.*alp3/st2/xi_t*m2*log_t*(3.*(s - 2.*m2)/u0 + (10.*m2 - 3.*u0)/(s - 4.*m2));

    // box diagram, factorized part
    // equation (54)
    double delta_box = (1. + xi_s2)/xi_s*(-4.*log_s*log_m + log_s*log_s - 2*pi2 + 4.*Li2_sp)
                       + (1 + xi_u02)/xi_u0*(4.*log_u0*log_m - log_u0*log_u0 + 2.*log_2u0p*log_2u0p - pi2/3. + 4*Li2_u0);

    // box diagram, non-factorized part
    // equation (A.1) and (A.2)
    double sig_b1, sig_b2;

    sig_b1 = 1./12./xi_s/t * ((xi_s2 + 1.)*(xi_s4 - 6.*xi_s2 - 3.)*s2t
                              - 2.*xi_s2*pow3(xi_s2 + 1.)*s3 - 12.*xi_s2*st2 - 4.*t3)
                           * (4.*pi2 + 3.*log_s*log_s - 6.*log_2s*log_2s - 12.*Li2_s - 6.*log_s*log(-xi_s2*s/t))
             + 1./12./xi_t2/xi_t/t * (xi_t2*(xi_t2 - 3.)*(3.*xi_t2 + 1.)*t3
                                      - 2.*u0t2*(3.*xi_t4*xi_t2 + 2.*xi_t4 + 10.*xi_t2 - 1.)
                                      - 4.*u02t*(5.*xi_t4 + 4.*xi_t2 - 1.) - 16.*xi_t2*u03)
                                   * (4.*pi2 - 6.*log_2t*log_2t + 3.*log_t*log_t - 12.*Li2_t)
             + 1./xi_s*log_s*(xi_s2*s + t)*((xi_s4 + xi_s2)*s - 2.*(xi_s2 - 2.)*t)
             + log_4t*(2.*t*t - (xi_s4 + xi_s2)*s*s + (3.*xi_s2 - 1.)*s*t - 2.*s*(t + 2.*u0)/xi_t2);

    sig_b2 = 1./12./xi_u0/t * (4.*t3 - 2.*u0t2*(xi_u04 - 6.*xi_u02 - 1.)
                               + u02t*(-xi_u04*xi_u02 + xi_u04 + 9.*xi_u02 + 7.) + 2.*pow3(xi_u02 + 1.)*u03)
             + (3.*log_u0*log_u0 - 6.*log_2u0*log_2u0 - 12.*Li2_u0 + 6.*log_u0*log(xi_u02*u0/t) + pi2)
             + 1./12./xi_t2/xi_t/t * (xi_t2*(-xi_t4 + 2.*xi_t2 + 3.)*t3
                                      + 2.*(xi_t4*xi_t2 - 4.*xi_t4 + 8.*xi_t2 + 1.)*u0t2
                                      + 4.*(3.*xi_t4 + 1.)*u02t + 16.*xi_t2*u03)
                                   * (-6.*log_2t*log_2t + 3.*log_t*log_t - 12.*Li2_t + 4.*pi2)
             + log_4t*(2.*u0/xi_t2*(xi_t2*t + t + 2.*u0) + (t - u0)*(2.*t + xi_u02*u0 + u0))
             - 1./xi_u0*log_u0*(xi_u02*(t - u0) - 2.*t)*(2.*t + xi_u02*u0 + u0);

    double sig_box = alp3/xi_s2/s2t/u0*(sig_b1 + sig_b2);

    return sig0*(1. + alp_pi*delta_vac + 2.*alp_pi*delta_ver + alp_pi/2.*delta_box)
           + sig_AMM + sig_box;
}


// get cros section
// input beam energy (MeV), angle (deg)
// output Born level cross section (nb)
double PRadMollerGen::GetBornXS(const double &Es, const double &angle)
{
    double m = cana::ele_mass;
    double theta = angle*cana::deg2rad;
    double p_tot = sqrt(pow2(Es) - pow2(m));

    // incident electron kinematics
    double k1[4] = {Es, 0., 0., p_tot};
    double p1[4] = {m, 0., 0., 0.};


    // Energy of the scattered electron k1
    double cos_E = (Es - m)*pow2(cos(theta));
    double E1 = m*(Es + m + cos_E)/(Es + m - cos_E);
    double k1_tot = sqrt(pow2(E1) - pow2(m));

    double k2[4] = {E1, k1_tot*sin(theta), 0., k1_tot*cos(theta)};

    // convert MeV^-2 to nbarn
    return cana::hbarc2*1e7*(moller_nonrad(p1, k1, k2, 0) + moller_nonrad(k1, p1, k2, 0));
}

double PRadMollerGen::GetNonRadXS(const double &Es, const double &angle)
{
    double m = cana::ele_mass;
    double theta = angle*cana::deg2rad;
    double p_tot = sqrt(pow2(Es) - pow2(m));

    // incident electron kinematics
    double k1[4] = {Es, 0., 0., p_tot};
    double p1[4] = {m, 0., 0., 0.};


    // Energy of the scattered electron k1
    double cos_E = (Es - m)*pow2(cos(theta));
    double E1 = m*(Es + m + cos_E)/(Es + m - cos_E);
    double k1_tot = sqrt(pow2(E1) - pow2(m));

    double k2[4] = {E1, k1_tot*sin(theta), 0., k1_tot*cos(theta)};

    // convert MeV^-2 to nbarn
    return cana::hbarc2*1e7*(moller_nonrad(p1, k1, k2) + moller_nonrad(k1, p1, k2));
}

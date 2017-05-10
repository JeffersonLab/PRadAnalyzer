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

// TODO not finished, now only has the non-radiative part
// v_min is not used because it is to separate the radiative part

#include "PRadMollerGen.h"
#include "canalib.h"
#include <cmath>
#include <iostream>
#include <iomanip>

//#define MOLLER_TEST_URA

// some constant values to be used
const double m = cana::ele_mass;
const double m2 = m*m;
const double alp_pi = cana::alpha/cana::pi;
const double pi2 = cana::pi*cana::pi;
const double alp2 = cana::alpha*cana::alpha;
const double alp3 = alp2*cana::alpha;
// convert MeV^-2 to nbarn
const double unit = cana::hbarc2*1e7;
// square
inline double pow2(double val) {return val*val;}
// cubic
inline double pow3(double val) {return val*val*val;}
// power of 4
inline double pow4(double val) {double val2 = pow2(val); return val2*val2;}




// constructor
PRadMollerGen::PRadMollerGen(double vmin, double vmax)
: v_min(vmin), v_cut(vmax)
{
    // place holder
}

// destructor
PRadMollerGen::~PRadMollerGen()
{
    // place holder
}

// get cross section
// input beam energy (MeV), angle (deg)
// output Born, non-radiative, radiative cross sections (nb)
void PRadMollerGen::GetXS(double Es, double angle, double &sig_born, double &sig_nrad, double &sig_rad)
const
{
    double theta = angle*cana::deg2rad;
    // conversion from dsigma/dy (y = -t/s) to dsigma/dOmega
    double jacob = 4.*m*(Es - m)/pow2(Es + m -(Es - m)*pow2(cos(theta)))/2./cana::pi;
    double p_tot = sqrt(pow2(Es) - pow2(m));

    // incident electron kinematics
    double k1[4] = {Es, 0., 0., p_tot};
    double p1[4] = {m, 0., 0., 0.};


    // Energy of the scattered electron k1
    double cos_E = (Es - m)*pow2(cos(theta));
    double E1 = m*(Es + m + cos_E)/(Es + m - cos_E);
    double k1_tot = sqrt(pow2(E1) - pow2(m));

    double k2[4] = {E1, k1_tot*sin(theta), 0., k1_tot*cos(theta)};

    // Mandelstam variables
    double s, t, u0;
    // s = (k1 + p1)^2
    s = pow2(k1[0] + p1[0]) - pow2(k1[1] + p1[1]) - pow2(k1[2] + p1[2])
        - pow2(k1[3] + p1[3]);

    // t = (k2 - k1)^2
    t = pow2(k2[0] - k1[0]) - pow2(k2[1] - k1[1]) - pow2(k2[2] - k1[2])
        - pow2(k2[3] - k1[3]);

    // u0 = (k2 - p1)^2
    u0 = pow2(k2[0] - p1[0]) - pow2(k2[1] - p1[1]) - pow2(k2[2] - p1[2])
         - pow2(k2[3] - p1[3]);

    // virtual photon part of Moller cross section
    // the t and u channels should be calculated together
    double sig_0t, sig_0u, sig_St, sig_Su, sig_vertt, sig_vertu, sig_Bt, sig_Bu;
    // t channel
    moller_Vph(s, t, u0, sig_0t, sig_St, sig_vertt, sig_Bt);
    // u0 channel
    moller_Vph(s, u0, t, sig_0u, sig_Su, sig_vertu, sig_Bu);

    // infrared divergent part of real photon emission
    // the t and u channels are not separated
    double delta_1H, delta_1S, delta_1inf;
    moller_IR(s, t, u0, delta_1H, delta_1S, delta_1inf);

    // initialize MERADGEN
    merad_init(Es);
    double v_limit = (s + t - 4.*m2)*0.999;
    double v_max = (v_limit > v_cut) ? v_cut : v_limit;
    // the "soft" Bremsstrahlung part of the radiative cross section
    // blow v_min, photon emission is not detectable
    double sig_Fs = merad_sigfs(v_min, t, 0., sig_0t + sig_0u);
    // the "hard" Bremsstrahlung part of the radiative cross section
    double sig_Fh = merad_sigfh(v_min, v_max, t, 0.);

    // t and u channels together
    // born level cross section
    sig_born = (sig_0t + sig_0u)*jacob*unit;
    sig_nrad = ((1. + alp_pi*(delta_1H + delta_1S))*exp(alp_pi*delta_1inf)*(sig_0t + sig_0u)
                + sig_St + sig_Su + sig_vertt + sig_vertu + sig_Bt + sig_Bu + sig_Fs)
               * jacob*unit;
    sig_rad = sig_Fh*jacob*unit;
}

// Cross section including virtual photon part for Moller scattering
// input Mandelstam variables s, t, u0
// output cross section including virtual photon effects
void PRadMollerGen::moller_Vph(double s, double t, double u0,
                               double &sig_0, double &sig_S, double &sig_vert, double &sig_B)
const
{
    double s2t = s*s*t, st2 = s*t*t, s3 = pow3(s);
    double u02t = u0*u0*t, u0t2 = u0*t*t, u03 = pow3(u0), t3 = pow3(t);

    // frequently used variables
    double xi_s = sqrt(1. - 4.*m2/s);
    double xi_t = sqrt(1 - 4.*m2/t);
    double xi_u0 = sqrt(1. - 4.*m2/u0);
    double xi_s2 = xi_s*xi_s, xi_s4 = xi_s2*xi_s2;
    double xi_t2 = xi_t*xi_t, xi_t4 = xi_t2*xi_t2;
    double xi_u02 = xi_u0*xi_u0, xi_u04 = xi_u02*xi_u02;

    // equation (49), Born Level
    // NOTE that there is an additional s in the denominator, but this s leads
    // to a wrong dimension, and disagrees with the URA form
    sig_0 = (u0*u0/xi_s2/4./s*(4.*xi_u04 - pow2(1. - xi_u02)*(2. + t/u0)) - s*s*xi_s4/u0)
            * 2.*cana::pi*alp2/t/t;

    // singularity term, appears in delta_ver and delta_box
    // we only need the divergence free part of the sigma_ver and simga_box,
    // which are obtained by substituting lamda = m so log(lamda/m) = 0
    double log_m = 0.; // log(lamda/m), where lamda is the infinitesimal photon mass

    // other frequently used variables
    // Q^2 (-t) related, equation (27) - (29)
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

    // equation (50)
    sig_S = alp_pi*delta_vac*sig_0;

    // vertex correction, factorized part
    // equation (36) with Q^2 -> -t
    double delta_vert = 2.*(Q2_m*L_m - 1.)*log_m + (4.*m2 - 3./2.*t)*L_m - 2.
                        - Q2_m/slamda_m*(lamda_m*L_m*L_m/2. + 2.*cana::spence(2.*slamda_m/(slamda_m - t)) - pi2/2.);

    // vertex correction, non-factorized part, anomalous magnetic moment
    // euqation (52)
    double sig_AMM = 4.*alp3/st2/xi_t*m2*log_t*(3.*(s - 2.*m2)/u0 + (10.*m2 - 3.*u0)/(s - 4.*m2));

    // equation (51)
    sig_vert = 2.*alp_pi*delta_vert*sig_0 + sig_AMM;

    // box diagram, factorized part
    // equation (54)
    double delta_box = (1. + xi_s2)/xi_s*(-4.*log_s*log_m + log_s*log_s - 2*pi2 + 4.*Li2_sp)
                       + (1 + xi_u02)/xi_u0*(4.*log_u0*log_m - log_u0*log_u0 + 2.*log_2u0p*log_2u0p - pi2/3. + 4*Li2_u0);

    // box diagram, non-factorized part
    // equation (A.1) and (A.2)
    double sig_B1, sig_B2;

    sig_B1 = 1./12./xi_s/t * ((xi_s2 + 1.)*(xi_s4 - 6.*xi_s2 - 3.)*s2t
                              - 2.*xi_s2*pow3(xi_s2 + 1.)*s3 - 12.*xi_s2*st2 - 4.*t3)
                           * (4.*pi2 + 3.*log_s*log_s - 6.*log_2s*log_2s - 12.*Li2_s - 6.*log_s*log(-xi_s2*s/t))
             + 1./12./xi_t2/xi_t/t * (xi_t2*(xi_t2 - 3.)*(3.*xi_t2 + 1.)*t3
                                      - 2.*u0t2*(3.*xi_t4*xi_t2 + 2.*xi_t4 + 10.*xi_t2 - 1.)
                                      - 4.*u02t*(5.*xi_t4 + 4.*xi_t2 - 1.) - 16.*xi_t2*u03)
                                   * (4.*pi2 - 6.*log_2t*log_2t + 3.*log_t*log_t - 12.*Li2_t)
             + 1./xi_s*log_s*(xi_s2*s + t)*((xi_s4 + xi_s2)*s - 2.*(xi_s2 - 2.)*t)
             + log_4t*(2.*t*t - (xi_s4 + xi_s2)*s*s + (3.*xi_s2 - 1.)*s*t - 2.*s*(t + 2.*u0)/xi_t2);

    sig_B2 = 1./12./xi_u0/t * (4.*t3 - 2.*u0t2*(xi_u04 - 6.*xi_u02 - 1.)
                               + u02t*(-xi_u04*xi_u02 + xi_u04 + 9.*xi_u02 + 7.) + 2.*pow3(xi_u02 + 1.)*u03)
             + (3.*log_u0*log_u0 - 6.*log_2u0*log_2u0 - 12.*Li2_u0 + 6.*log_u0*log(xi_u02*u0/t) + pi2)
             + 1./12./xi_t2/xi_t/t * (xi_t2*(-xi_t4 + 2.*xi_t2 + 3.)*t3
                                      + 2.*(xi_t4*xi_t2 - 4.*xi_t4 + 8.*xi_t2 + 1.)*u0t2
                                      + 4.*(3.*xi_t4 + 1.)*u02t + 16.*xi_t2*u03)
                                   * (-6.*log_2t*log_2t + 3.*log_t*log_t - 12.*Li2_t + 4.*pi2)
             + log_4t*(2.*u0/xi_t2*(xi_t2*t + t + 2.*u0) + (t - u0)*(2.*t + xi_u02*u0 + u0))
             - 1./xi_u0*log_u0*(xi_u02*(t - u0) - 2.*t)*(2.*t + xi_u02*u0 + u0);

    // equation (53)
    sig_B = alp_pi/2.*delta_box*sig_0 + alp3/xi_s2/s2t/u0*(sig_B1 + sig_B2);
}

// Infrared part of the Moller cross sections with real photon emission
// input Mandelstam variables
// NOTE that the t and u0 channel are not separated
// output 3 factorized part for the infrared part of the radiative cross section
void PRadMollerGen::moller_IR(double s, double t, double u0,
                              double &delta_1H, double &delta_1S, double &delta_1inf)
const
{
    // frequently used variables
    double xi_s = sqrt(1. - 4.*m2/s);
    double xi_t = sqrt(1 - 4.*m2/t);
    double xi_u0 = sqrt(1. - 4.*m2/u0);
    double xi_s2 = xi_s*xi_s;
    double xi_t2 = xi_t*xi_t;
    double xi_u02 = xi_u0*xi_u0;
    double log_s = log((1. + xi_s)/(1. - xi_s));
    double log_t = log((1. + xi_t)/(xi_t - 1.));
    double log_u0 = log((1. + xi_u0)/(xi_u0 - 1.));

    // equation (A.5) - (A.13)
    double v_limit = (s*t + sqrt(s*(s - 4.*m2)*t*(t - 4.*m2)))/2./m2;
    double v_max = (v_min > v_limit) ? v_limit : v_min;
    double z_u1 = sqrt((xi_u02*(v_max + u0) - v_max)/u0)/xi_u0;
    double z_u2 = sqrt((v_max + xi_u02*u0)/(v_max + u0))/xi_u0;
    auto H = [v_max](const double &ch)
             {
                 double xi_ch = sqrt(1. - 4.*m2/ch), xi_ch2 = xi_ch*xi_ch;
                 double z_ch = xi_ch/v_max*(sqrt(xi_ch2*ch*ch - 2.*ch*v_max + v_max*v_max) - xi_ch*ch) + 1.;
                 double z_1 = 1. + xi_ch;
                 double z_2 = pow2(1. + xi_ch)/(1. - xi_ch);
                 double z_3 = 1. - xi_ch;
                 double z_4 = pow2(1. - xi_ch)/(1. + xi_ch);
                 double Li2_z1 = cana::spence(z_ch/z_1);
                 double Li2_z2 = cana::spence(z_ch/z_2);
                 double Li2_z3 = cana::spence(z_ch/z_3);
                 double Li2_z4 = cana::spence(z_ch/z_4);
                 // NOTICE: fabs is not in the original function
                 // however, the term in log can be negative and thus result in
                 // undefined behavior for real numbers
                 double log_term = log(pow2((xi_ch + 1.)/(xi_ch - 1.)))*log(fabs((pow2(z_ch - 1.) - xi_ch2)/(1. - xi_ch2)));

                 return (xi_ch2 + 1.)/2./xi_ch * (Li2_z1 + Li2_z2 - Li2_z3 - Li2_z4 - log_term);
             };

    // equation (A.3)
    delta_1H = log(1. + v_max/m2) + H(s) - H(t)
               + (xi_u02 + 1.)/2./xi_u0
                  * (cana::spence(4.* xi_u0/pow2(xi_u0 + 1.))
                     - cana::spence(-4.*xi_u0/pow2(xi_u0 - 1.))
                     - 2.*cana::spence(2.*xi_u0/(xi_u0 - 1.))
                     + 2.*cana::spence(2.*xi_u0/(xi_u0 + 1.))
                     + cana::spence(2.*(z_u1 - 1.)*xi_u0/pow2(xi_u0 - 1.))
                     + cana::spence(-2.*(z_u1 + 1.)*xi_u0/pow2(xi_u0 - 1.))
                     - cana::spence(-2.*(z_u1 - 1.)*xi_u0/pow2(xi_u0 + 1.))
                     - cana::spence(2.*(z_u1 + 1.)*xi_u0/pow2(xi_u0 + 1.))
                     + 2.*cana::spence(-(z_u2 - 1.)*xi_u0/(xi_u0 - 1.))
                     + 2.*cana::spence((z_u2 + 1.)*xi_u0/(xi_u0 - 1.))
                     - 2.*cana::spence((z_u2 + 1.)*xi_u0/(xi_u0 + 1.))
                     - 2.*cana::spence((1. - z_u2)*xi_u0/(xi_u0 + 1.))
                     + 2.*log_u0*log((xi_u02*z_u2*z_u2 - 1.)/(xi_u02 - 1.)));


    // equation (A.14) - (A.15)
    auto S_phi = [](const double &s_1, const double &s_2, const double &s_3)
                 {
                     double lamda_1 = s_1*s_1 - 16.*m2*m2, slamda_1 = sqrt(lamda_1);
                     double lamda_2 = s_2*s_2 - 16.*m2*m2, slamda_2 = sqrt(lamda_2);
                     double lamda_3 = s_3*s_3 - 16.*m2*m2, slamda_3 = sqrt(lamda_3);
                     // z_u and z_d
                     double z_ud[2] = {slamda_1/slamda_2 - 1., (s_1*s_2 - 4.*m2*s_3)/lamda_2 - 1.};
                     // z_1, z_2, z_3, z_4
                     double z[4] = {1./slamda_2*(4.*m2*(s_3 - slamda_3)/(s_2 - slamda_2) - s_1 - slamda_2),
                                    1./slamda_2*(4.*m2*(s_3 + slamda_3)/(s_2 - slamda_2) - s_1 - slamda_2),
                                    1./slamda_2*(s_1 - slamda_2 - 4.*m2*(s_3 + slamda_3)/(s_2 + slamda_2)),
                                    1./slamda_2*(s_1 - slamda_2 - 4.*m2*(s_3 - slamda_3)/(s_2 + slamda_2))};
                     // Sj
                     double Sj[4] = {1, 1, -1, -1};
                     // (-1)^(i + 1), i from 1 to 4 but index is from 0 to 3
                     double Si[4] = {1, -1, 1, -1};
                     // z_u term - z_d term
                     double Sk[2] = {1, -1};

                     double result = 0.;
                     for(int k = 0; k < 2; ++k)
                     {
                         // TODO, check with authors
                         // it is noted in the reference that
                         // S_phi(s1, s2, s3) == S_phi(s2, s1, s3)
                         // but this part could not satisfy the relation
                         double term = log((s_2 - slamda_2)/(s_2 + slamda_2))
                                       * log((z_ud[k] - z[0])*(z_ud[k] - z[2])/(z_ud[k] - z[1])/(z_ud[k] - z[3]));

                         double sum_term = 0.;
                         for(int i = 0; i < 4; ++i)
                         {
                             for(int j = 0; j < 4; ++j)
                             {
                                 double sign = Sj[j]*Si[i];
                                 if(i == j) {
                                     sum_term += sign*log(pow2(z_ud[k] - z[i]))/2.;
                                 } else {
                                     // the input to log may be negative and thus
                                     // the result may be a complex number
                                     // TODO check with authors for this part
                                     double spence_term = cana::spence((z_ud[k] - z[i])/(z[j] - z[i]));
                                     sum_term += sign*(log(fabs(z_ud[k] - z[i]))*log(fabs(z[i] - z[j])) - spence_term);
                                 }
                             }
                         }

                         result += s_3/2./slamda_3*(/*term +*/ sum_term)*Sk[k];
                     }
                     return result;
                 };

    // equation (A.4)
    delta_1S = (xi_s2 + 1.)/2./xi_s*(log_s*log_s + log_s + cana::spence(4.*xi_s/pow2(xi_s + 1.)))
               - (xi_t2 + 1.)/2./xi_t*(log_t*log_t - log_t + cana::spence(4.*xi_t/pow2(xi_t + 1.)))
               - (xi_u02 + 1.)/2./xi_u0*(log_u0*log_u0 - log_u0 + cana::spence(4.*xi_u0/pow2(xi_u0 + 1.)))
               - S_phi(-(xi_u02 + 1.)*u0, (xi_s2 + 1.)*s, -(xi_t2 + 1.)*t)
               + S_phi(-(xi_u02 + 1.)*u0, -(xi_t2 + 1.)*t, (xi_s2 + 1.)*s)
               - S_phi(-(xi_t2 + 1.)*t, (xi_s2 + 1.)*s, -(xi_u02 + 1.)*u0) + 1.;


    // equation (60)
    double J_0 = -2.*((xi_s2 + 1.)/xi_s*log_s - (xi_t2 + 1.)/xi_t*log_t
                      - (xi_u02 + 1.)/xi_u0*log_u0 + 2.);
    // equation (66)
    delta_1inf = J_0*log(v_max/m2);

// test difference with URA
#ifdef MOLLER_TEST_URA
    // equation (61)
    double delta_1H_URA = log(-t/m2)*(log(pow2(t*(s + t))*(s - v_max)/s/v_max/(v_max - t)/pow2(s + t - v_max)) + 1.)
                          - pow2(log(-t/m2))/2.
                          + 2.*(-cana::spence(v_max/(s + t)) + cana::spence(v_max/s) - cana::spence(v_max/t))
                          + cana::spence((s - v_max)/s)
                          - cana::spence((t - v_max)/t)
                          + log((s + t)/(s + t - v_max))*log((s + t)*(s + t - v_max)/t/t)
                          + log((s - v_max)/s)*log((v_max - s)/t)
                          - pow2(log(-v_max/t))/2.
                          - pow2(log(1. - v_max/t))
                          + log(-v_max/t) - pi2/6.;

    // equation (62)
    double delta_1S_URA = 1. - (log(-t/m2) - 1.)*log(s*(s + t)/t/t)
                          + log(-t/m2)*(3. - 2.*log((s + t)/s))
                          - 5./2.*pow2(log(-t/m2)) - pow2(log((s + t)/s))/2. - pi2/3.;

    // equation (63)
    double J_0_URA = -4.*(1. + log(m2*s/t/u0));

    std::cout << S_phi(-(xi_u02 + 1.)*u0, (xi_s2 + 1.)*s, -(xi_t2 + 1.)*t) << ", "
              << S_phi((xi_s2 + 1.)*s, -(xi_u02 + 1.)*u0, -(xi_t2 + 1.)*t)
              << std::endl;

    std::cout << " | d_1H: "
              << std::setw(10) << delta_1H << ", "
              << std::setw(10) << delta_1H_URA
              << " | d_1S: "
              << std::setw(10) << delta_1S << ", "
              << std::setw(10) << delta_1S_URA
              << " | J_0:  "
              << std::setw(10) << J_0 << ", "
              << std::setw(10) << J_0_URA << " |"
              << std::endl;
#endif // MOLLER_TEST_URA
// end test
}

// real photon emission part of the Moller scattering
// input four momentum of the particles - arrays with size 4, energy is at 0
// p1, initial target four-momentum
// k1, initial beam electron four-momentum
// k2, final beam electron four-momentum
// type, 0. born, others. non-rad
// output cross section with the same unit from input
/*
double PRadMollerGen::moller_rad(double *p1, double *k1, double *k2)
{
    return 0.;
}
*/

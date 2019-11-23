//============================================================================//
// Monte Carlo event generator for UNPOLARIZED Elastic e-p scattering beyond  //
// ultra-relativistic approximation (URA)                                     //
//----------------------------------------------------------------------------//
// References:                                                                //
// [1] Eur. Phys. J. A51, 1 (2015)                                            //
//     Radiative corrections beyond the ultra relativistic limit in           //
//     unpolarized ep elastic and Moller scatterings for the PRad Experiment  //
//     at Jefferson Laboratory                                                //
//     I. Akushevich, H. Gao, A. Ilyichev, and M. Meziane                     //
// [2] newep event generator                                                  //
//     https://github.com/JeffersonLab/PRadSim/tree/master/evgen/newep        //
//     C. Gu                                                                  //
//----------------------------------------------------------------------------//
// Code Developer: Chao Peng  10/06/2017                                      //
//============================================================================//


#include "PRadEpElasGen.h"
#include <cmath>
#include <random>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include "PRadBenchMark.h"

#define PROGRESS_EVENT_COUNT 1000
#define PROGRESS_BIN_COUNT 10

//============================================================================//
// Helper constants, structures and functions                                 //
//============================================================================//

// some constant values to be used
const double m = cana::ele_mass;
const double m2 = m*m;
const double M = cana::proton_mass;
const double M2 = M*M;
const double twopi = 2.*cana::pi;
const double alp_pi = cana::alpha/cana::pi;
const double pi2 = cana::pi*cana::pi;
const double alp2 = cana::alpha*cana::alpha;
const double alp3 = alp2*cana::alpha;
// convert MeV^-2 to nbarn
// const double unit = cana::hbarc2*1e7;

// random number generator
static cana::rand_gen<> rng;

// some inlines
// square
inline double pow2(double val) {return val*val;}
// cubic
inline double pow3(double val) {return val*val*val;}
// power of 4
inline double pow4(double val) {double val2 = pow2(val); return val2*val2;}

inline double tau_ext1(double v, double S, double Q2)
{
    return 1./(8.*m2*M2 + 2.*S*S)*(4.*m2*(Q2 + 1.) - 2.*S*Q2);
}

inline double tau_ext2(double v, double S, double Q2)
{
    double X = S - Q2 - v;
    return 1./(8.*m2*M2 + 2.*X*X)*(4.*m2*(Q2 + 1.) + 2.*X*Q2);
}

//============================================================================//
// Constructor, destructor                                                    //
//============================================================================//

// constructor
PRadEpElasGen::PRadEpElasGen(double vmin, double vmax, int nbins, double q2_res, double v_res)
: v_min(vmin), v_cut(vmax), /* min_bins(nbins), q2_prec(q2_res), */ v_prec(v_res)
{
    // place holder
}

// destructor
PRadEpElasGen::~PRadEpElasGen()
{
    // place holder
}



//============================================================================//
// Public functions                                                           //
//============================================================================//

// get the electricmagnetic form factors for proton
void PRadEpElasGen::GetEMFF(double Q2, double &GE, double &GM)
const
{
    static double GEp1[6] = {1., 2.90966, -1.11542229, 3.866171e-2, 0., 0.};
    static double GEp2[6] = {1., 14.5187212, 40.88333, 99.999998, 4.579e-5, 10.3580447};
    static double GMp1[6] = {1., -1.43573, 1.19052066, 2.5455841e-1, 0., 0.};
    static double GMp2[6] = {1., 9.70703681, 3.7357e-4, 6.0e-8, 9.9527277, 12.7977739};

    // convert MeV^2 to GeV^2
    double tau = Q2/4./M2;
    double x[6];
    x[0] = 1.;
    for(int i = 1; i < 6; ++i) x[i] = x[i - 1]*tau;

    double GEp_nu = 0., GEp_de = 0., GMp_nu = 0., GMp_de = 0.;
    for(int i = 0; i < 6; ++i)
    {
        GEp_nu += GEp1[i]*x[i];
        GEp_de += GEp2[i]*x[i];
        GMp_nu += GMp1[i]*x[i];
        GMp_de += GMp2[i]*x[i];
    }

    GE = GEp_nu/GEp_de;
    GM = 2.792782*GMp_nu/GMp_de;
}

// translate EM form factors to hadronic structure functions
void PRadEpElasGen::GetHadStrFunc(double Q2, double &F1, double &F2)
const
{
    double GE, GM;
    GetEMFF(Q2, GE, GM);
    double tau = Q2/4./M2;

    F1 = 4.*tau*M2*GM*GM;
    F2 = 4.*M2*(GE*GE + tau*GM*GM)/(1. + tau);
}

// get differential cross section dsigma/dQ2
// input variables S, Q^2 in MeV^2
// output Born, non-radiative, radiative cross sections (MeV^-4)
void PRadEpElasGen::GetXSdQsq(double S, double Q2, double &sig_born, double &sig_nrad, double &sig_rad)
const
{
    // equation (12), limitation of photonic variable v, scaled by 0.99
    // substitute Q2*(Q2 + 4m2) with lambda_m in equation (29)
    double lambda_S = S*S - 4.*m2*M2, lambda_m = Q2*(Q2 + 4.*m2);
    double v_limit = 0.99*2.*Q2*(lambda_S - Q2*(S + m2 + M2))/(Q2*(S + 2.*m2) + sqrt(lambda_S*lambda_m));

    double v2 = cana::clamp(v_cut, v_cut, v_limit);
    double v1 = cana::clamp(v_min, v_min, v2);

    sig_born = SigmaBorn(S, Q2);

    // equation (39) without hard photon emission part of sigmaF
    sig_nrad = SigmaVphIR(S, Q2, v1) + SigmaFs(1e-6, v1, S, Q2);

    // hard photon emission part
    sig_rad = SigmaFh(v1, v2, S, Q2);
}

// Differential cross section dsigma/dQ^2 at Born level
// input variable S, Q^2 in MeV^2
double PRadEpElasGen::SigmaBorn(double S, double Q2)
const
{
    double F01, F02;
    GetHadStrFunc(Q2, F01, F02);

    // equation (3) ~ (8)
    double X = S - Q2;
    double lambda_S = S*S - 4.*m2*M2;
    double theta_B1 = Q2 - 2.*m2;
    double theta_B2 = 1./(2.*M2)*(S*X - M2*Q2);

    return twopi*alp2/lambda_S/Q2/Q2*(F01*theta_B1 + F02*theta_B2);
}

// S_phi function defined in equation (35)
inline double S_phi(double s, double l, double a, double b)
{
    static double delta[4] = {1, 1, -1, -1};
    static double sign[4] = {-1, 1, -1, 1};
    double sqrt_l = sqrt(l), sqrt_b = sqrt(b);
    double D = (s + a)*(l*a - s*b) + pow2(l + b)/4.;
    double gamma_u = (sqrt(b + l) - sqrt_b)/sqrt_l;
    // NOTICE
    // gamma_1,2 is different with paper [1], but it has to be dimension-less
    // so the paper probably has a misprint on this term
    double gamma[2] = { -pow2(sqrt_b - sqrt_l)/(b - l),
                        pow2(sqrt_b + sqrt_l)/(b - l) };
    double i_sign, a_j, tau_j, delta_j, gamma_i, gamma_jk;

    auto S_term = [&i_sign, &gamma_i, &gamma_jk] (double g)
                  { return cana::spence((gamma_i - g)/(gamma_i - gamma_jk))
                           + cana::spence((g + i_sign)/(gamma_jk + i_sign)); };

    double res = 0.;
    for(int i = 0; i < 2; ++i)
    {
        i_sign = sign[i];
        gamma_i = gamma[i];
        for(int j = 0; j < 4; ++j)
        {
            delta_j = delta[j];
            a_j = s - delta_j*sqrt_l;
            tau_j = -a*sqrt_l + delta_j*(b - l)/2. + sign[j]*sqrt(D);
            for(int k = 0; k < 2; ++k)
            {
                gamma_jk = -(a_j*sqrt_b - sign[k]*sqrt(b*a_j*a_j + tau_j*tau_j))/tau_j;
                res += i_sign*delta_j*(S_term(gamma_u) - S_term(gamma[0]));
            }
        }
    }

    return res*s/2./sqrt_l;
}

// Cross section including virtual photon part and the Infrared part of the
// photon emission of the ep elastic cross section
// input variables S, Q^2 in MeV^2, v_max for the photonic variable in MeV^2
double PRadEpElasGen::SigmaVphIR(double S, double Q2, double v1)
const
{
    // substitute S - Q2 with X
    double X = S - Q2;

    // equation (27) ~ (34)
    double Q2_m = Q2 + 2.*m2;
    double lambda_m = Q2*(Q2 + 4.*m2), sqrt_lm = sqrt(lambda_m);
    double L_m = log((sqrt_lm + Q2)/(sqrt_lm - Q2))/sqrt_lm;
    double lambda_S = S*S - 4.*m2*M2, sqrt_lS = sqrt(lambda_S);  // equation (4)
    double L_S = log((S + sqrt_lS)/(S - sqrt_lS))/sqrt_lS;
    double lambda_X0 = X*X - 4.*m2*M2, sqrt_lX0 = sqrt(lambda_X0);
    double L_X0 = log((X + sqrt_lX0)/(X - sqrt_lX0))/sqrt_lX0;
    double a = (S*X - 2*M2*(Q2 - 2.*m2))/2./M2;
    double b = (Q2*(S*X - M2*Q2) - m2*Q2*(Q2 + 4.*M2))/M2;

    double v_max = 2.*Q2*(lambda_S - Q2*(S + m2 + M2))/(Q2*(S + 2.*m2) + sqrt(lambda_S*lambda_m));

    double F01, F02;
    GetHadStrFunc(Q2, F01, F02);
    // equation (7), (8)
    double theta_B1 = Q2 - 2.*m2;
    double theta_B2 = 1./(2.*M2)*(S*X - M2*Q2);

    double sig_born, sig_AMM, delta_VR, delta_vac, delta_inf, delta_add;

    // equation (3)
    sig_born = twopi*alp2/lambda_S/Q2/Q2*(F01*theta_B1 + F02*theta_B2);

    // equation (40)
    delta_VR = 2.*(Q2_m*L_m - 1.)*log(v_max/m/M)
               + (S*L_S + X*L_X0)/2. + S_phi(Q2_m, lambda_m, a, b)
               + (3./2.*Q2 + 4.*m2)*L_m - 2.
               - Q2_m/sqrt_lm*(lambda_m*L_m*L_m/2.
                               + 2.*cana::spence((2.*sqrt_lm)/(Q2 + sqrt_lm))
                               - pi2/2.);

    // equation (41)
    delta_vac = 0.;
    double lepton_mass[3] = {cana::ele_mass, cana::mu_mass, cana::tau_mass};
    for(auto &vm : lepton_mass)
    {
        double vm2 = vm*vm;
        double vsqrt_lm = sqrt(Q2*(Q2 + 4.*vm2));
        double vL_m = log((vsqrt_lm + Q2)/(vsqrt_lm - Q2))/vsqrt_lm;
        delta_vac += 2./3.*(Q2 + 2.*vm2)*vL_m - 10./9. + 8./3.*vm2/Q2*(1. - 2.*vm2*vL_m);
    }


    // equation (42)
    delta_inf = (Q2_m*L_m - 1.)*log(v_max*v_max/S/X);

    // factorized part from sigma_Fs
    // adding this term just change the v_max in delta_VR to v_min in the first order
    delta_add = -2.*(Q2_m*L_m - 1.)*log(v_max/v1);

    // equation (38)
    sig_AMM = alp3*m2*L_m*(12.*M2*F01 - (Q2 + 4.*M2)*F02)/(2.*M2*Q2*lambda_S);

    return sig_born*(1. + alp_pi*(delta_VR + delta_vac - delta_inf))*exp(alp_pi*delta_inf)
           + alp_pi*delta_add*sig_born*exp(alp_pi*delta_inf) + sig_AMM;
}

// The Bremsstrahlung differential cross section with hard photon emission
// radiative cross section dsig/dQ2/dtau/dv/dphik
double PRadEpElasGen::SigmaBrem(double v, double tau, double phik, double S, double Q2, bool finite)
const
{
    // varibles to simplify equations
    double lambda_S = S*S - 4.*m2*M2;
    double R = v/(1. + tau);
    double t = Q2 + tau*R;
    double X = S - R - t;
    double SmX = S - X; // R + t = Q2 + v
    double SpX = S + X;

    // equation (13)
    double lambda_Y = SmX*SmX + 4.*M2*Q2;
    double sqrt_lY = sqrt(lambda_Y);
    double tau_min = (SmX - sqrt_lY)/2./M2;
    double tau_max = (SmX + sqrt_lY)/2./M2;

    // equation (22), (23)
    double lambda_Z = (tau - tau_min)*(tau_max - tau)*(S*X*Q2 - M2*Q2*Q2 - m2*lambda_Y);
    double sqrt_lZ = sqrt(lambda_Z);
    double z1 = (Q2*SpX + tau*(S*SmX + 2.*M2*Q2) - 2.*M*cos(phik)*sqrt_lZ)/lambda_Y;
    // NOTICE, probably typo in the original paper
    double z2 = (Q2*SpX + tau*(X*SmX - 2.*M2*Q2) - 2.*M*cos(phik)*sqrt_lZ)/lambda_Y;

    double F = 1./(sqrt_lY*twopi);
    double F_d = F/(z1*z2);
    double F_1p = F*(1./z1 + 1./z2);
    double F_2p = F*m2*(1./z2/z2 + 1./z1/z1);
    double F_2m = F*m2*(1./z2/z2 - 1./z1/z1);

    double F_IR = F_2p - (Q2 + 2.*m2)*F_d;

    // equation (43)
    double res;
    // second term
    if(finite) {
        res = alp_pi*F_IR/R*SigmaBorn(S, Q2)/(1. + tau);
    // first term
    } else {
        // equation (16) ~ (21)
        double theta_11 = 4.*(Q2 - 2.*m2)*F_IR;
        double theta_12 = 4.*tau*F_IR;
        double theta_13 = -4.*F - 2.*tau*tau*F_d;
        double theta_21 = 2.*(S*X - M2*Q2)*F_IR/M2;
        double theta_22 = (2.*SpX*F_2m + SpX*SmX*F_1p + 2.*(SmX - 2.*M2*tau)*F_IR - tau*SpX*SpX*F_d)/2./M2;
        double theta_23 = (4.*M2*F + (4.*m2 + 2.*M2*tau*tau - SmX*tau)*F_d - SpX*F_1p)/2./M2;

        // R^(j - 2)*theta_ij
        double theta_1j = theta_11/R + theta_12 + theta_13*R;
        double theta_2j = theta_21/R + theta_22 + theta_23*R;

        double F01, F02;
        GetHadStrFunc(t, F01, F02);

        res = -alp3/2./lambda_S*(theta_1j*F01 + theta_2j*F02)/t/t/(1. + tau);
    }

    return res;
}

// Analytical integration of SigmaBrem over phik dsig/dQ2/dtau/dv
// finite = true means the finite part of this differential cross section
double PRadEpElasGen::SigmaBrem_phik(double v, double tau, double S, double Q2, bool finite)
const
{
    // varibles to simplify equations
    double lambda_S = S*S - 4.*m2*M2;
    double R = v/(1. + tau);
    double t = Q2 + tau*R;
    double X = S - R - t;
    double SmX = S - X; // R + t = Q2 + v
    double SpX = S + X;

    double lambda_Y = SmX*SmX + 4.*M2*Q2;
    double sqrt_lY = sqrt(lambda_Y);

    // according to ELRADGEN code
    double b2 = (-lambda_Y*tau + SpX*SmX*tau + 2.*SpX*Q2)/2.;
    double b1 = (-lambda_Y*tau - SpX*SmX*tau - 2.*SpX*Q2)/2.;
    double c1 = -(4.*(M2*tau*tau - SmX*tau - Q2)*m2 - pow2(S*tau + Q2));
    double c2 = -(4.*(M2*tau*tau - SmX*tau - Q2)*m2 - pow2(tau*X - Q2));
    double sc1 = sqrt(c1);
    double sc2 = sqrt(c2);

    double F = 1./sqrt_lY;
    double F_d = (SpX*(SmX*tau + 2.*Q2))/(sc1*sc2*(sc1 + sc2));
    double F_1p = 1./sc1 + 1./sc2;
    double F_2p = m2*(b2/sc2/c2 - b1/sc1/c1);
    double F_2m = m2*(b2/sc2/c2 + b1/sc1/c1);

    double F_IR = F_2p - (Q2 + 2.*m2)*F_d;

    // equation (43)
    double res;
    // second term
    if(finite) {
        res = alp_pi*F_IR/R*SigmaBorn(S, Q2)/(1. + tau);
    // first term
    } else {
        // equation (16) ~ (21)
        double theta_11 = 4.*(Q2 - 2.*m2)*F_IR;
        double theta_12 = 4.*tau*F_IR;
        double theta_13 = -4.*F - 2.*tau*tau*F_d;
        double theta_21 = 2.*(S*X - M2*Q2)*F_IR/M2;
        double theta_22 = (2.*SpX*F_2m + SpX*SmX*F_1p + 2.*(SmX - 2.*M2*tau)*F_IR - tau*SpX*SpX*F_d)/2./M2;
        double theta_23 = (4.*M2*F + (4.*m2 + 2.*M2*tau*tau - SmX*tau)*F_d - SpX*F_1p)/2./M2;

        // R^(j - 2)*theta_ij
        double theta_1j = theta_11/R + theta_12 + theta_13*R;
        double theta_2j = theta_21/R + theta_22 + theta_23*R;

        double F01, F02;
        GetHadStrFunc(t, F01, F02);

        res = -alp3/2./lambda_S*(theta_1j*F01 + theta_2j*F02)/t/t/(1. + tau);
    }

    return res;
}

// numerical integration of SigmaBrem_phik over tau, dsig/dQ2/dv
// finite = true means the integration of the second term of equation (43)
double PRadEpElasGen::SigmaBrem_phik_tau(double v, double S, double Q2, bool finite, double prec)
const
{
    // wrapper of member function
    auto fn = [this] (double tau, double v, double S, double Q2, bool finite)
              { return SigmaBrem_phik(v, tau, S, Q2, finite); };

    // equation (13), tau range
    double sqrt_lY = sqrt(pow2(Q2 + v) + 4.*M2*Q2);
    double tau_min = (Q2 + v - sqrt_lY)/2./M2;
    double tau_max = (Q2 + v + sqrt_lY)/2./M2;

    // high peak at the center
    double tau_step = (tau_max - tau_min)*0.01;
    double tau_left = tau_ext1(v, S, Q2) - tau_step, tau_right = tau_ext2(v, S, Q2) + tau_step;
    double res1 = cana::simpson(fn, tau_min, tau_left, 200, v, S, Q2, finite);
    double res2 = cana::simpson_prec(fn, tau_left, tau_right, prec, v, S, Q2, finite);
    double res3 = cana::simpson(fn, tau_right, tau_max, 200, v, S, Q2, finite);
    return res1 + res2 + res3;
}

// hard Bremsstrahlung part
double PRadEpElasGen::SigmaFh(double v1, double v2, double S, double Q2)
const
{
    auto fn = [this] (double v, double S, double Q2)
              { return SigmaBrem_phik_tau(v, S, Q2, false, 1e-4); };

    return cana::simpson_prec(fn, v1, v2, v_prec, S, Q2);
}

// soft Bremsstrahlung part
double PRadEpElasGen::SigmaFs(double v1, double v2, double S, double Q2)
const
{
    // two terms are separated because it helps the integration converge
    auto fn = [this] (double v, double S, double Q2)
              { return SigmaBrem_phik_tau(v, S, Q2, true, 1e-4)
                       + SigmaBrem_phik_tau(v, S, Q2, false, 1e-4);};

    return cana::simpson_prec(fn, v1, v2, v_prec, S, Q2);
}

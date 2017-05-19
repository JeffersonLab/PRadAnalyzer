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

// TODO not finished, now only has the cross sections part
// need to implement the sampling part

#include "PRadMollerGen.h"
#include "canalib.h"
#include <cmath>
#include <random>
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

// some inlines
// square
inline double pow2(double val) {return val*val;}
// cubic
inline double pow3(double val) {return val*val*val;}
// power of 4
inline double pow4(double val) {double val2 = pow2(val); return val2*val2;}

// get the symmetric Moller pair angle in Lab frame
inline double get_sym_angle(double Es)
{
    return acos(sqrt((Es + m)/(Es + 3.*m)))*cana::rad2deg;
}

// get the polar angle of the other Moller electron in Lab frame
inline double get_other_polar(double Es, double angle)
{
    double theta = angle*cana::deg2rad;
    double cos2 = pow2(cos(theta));
    double rm = (Es - m)/(Es + m);
    double E1 = m*(1. + rm*cos2)/(1. - rm*cos2);
    double k2p = sqrt(pow2(E1) - pow2(m));
    double k1p = sqrt(pow2(Es) - pow2(m));
    return atan(k2p*sin(theta)/(k1p - k2p*cos(theta)))*cana::rad2deg;
}

// get the beta for CM frame
inline double get_CM_beta(double Es)
{
    return sqrt((Es - cana::ele_mass)/(Es + cana::ele_mass));
}

// boost the four momentum along z axis
inline void four_momentum_boost_z(double *p, double *pp, double beta)
{
    double gamma = sqrt(1./(1. - beta*beta));
    p[0] = gamma*(pp[0] - beta*pp[3]);
    p[1] = pp[1];
    p[2] = pp[2];
    p[3] = gamma*(pp[3] - beta*pp[0]);
}

// linear interpolation from two arrays (one for distribution, one for value)
inline double interp_dist(double rnd, double *dist, double *val, int np)
{
    double rdist = rnd*dist[np - 1];
    auto itv = cana::binary_search_interval(&dist[0], &dist[np], rdist);
    // should not happen
    if(itv.first == &dist[np] || itv.second == &dist[np])
    {
        std::cerr << "Interpolation error for random value = " << rdist
                  << ", distribution CDF starts at " << dist[0]
                  << ", ends at " << dist[np - 1]
                  << std::endl;
        return 0.;
    }

    int i1 = itv.first - dist;
    int i2 = itv.second - dist;

    if(i1 == i2) {
        return val[i1];
    } else {
        return cana::linear_interp(dist[i1], val[i1], dist[i2], val[i2], rdist);
    }
}

// Mandelstam variables for Moller process
inline void get_moller_stu(double Es, double angle, double &s, double &t, double &u)
{

    double theta = angle*cana::deg2rad;
    double k1p = sqrt(Es*Es - m2);
    double cosE = (Es - m)*pow2(cos(theta));
    double Ep = m*(Es + m + cosE)/(Es + m - cosE);
    double k2p = sqrt(Ep*Ep - m*m);

    s = 2.*m*(Es + m);
    t = pow2(Es - Ep) - pow2(k2p*sin(theta)) - pow2(k2p*cos(theta) - k1p);
    u = 4.*m2 - s - t;
}


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

// generate Moller events, unit is in MeV, degree
void PRadMollerGen::Generate(double Es, double min_angle, double max_angle, int nevents)
const
{
    // sanity check
    if(min_angle > max_angle || nevents < 1 || Es < 0.) {
        std::cerr << "Invalid input(s), please make sure these inputs are correct:\n"
                  << "Es = " << Es << " MeV\n"
                  << "Angle range = " << min_angle << " ~ " << max_angle << "\n"
                  << "Number of events: " << nevents
                  << std::endl;
        return;
    }

    // set up random number generator
    // a "true" random number engine, good for seed
    std::random_device rd;
    // Mersenne Twister pseudo-randome number generator
    std::mt19937 rng{rd()};
    // uniform distribution
    std::uniform_real_distribution<double> uni_dist(0., 1.);

    // since the cross section has already covered both t and u channels
    // we should calculate the if the angle ranges include the symmetric point
    // if it does, then the angle beyond symmetric point should not be generated
    // to avoid duplicating Moller pairs
    double sym_angle = get_sym_angle(Es);

    // choose lower angle side to generate
    if(sym_angle > min_angle && sym_angle < max_angle) {
        // update min angle so it can cover the input max_angle
        double min_angle2 = get_other_polar(Es, max_angle);
        if(min_angle > min_angle2) min_angle = min_angle2;
        max_angle = sym_angle;
    }

    // prepare grid for interpolation of angle
    int abins = 100;
    double angle_step = (max_angle - min_angle)/(double)abins;
    struct CDF_Angle {double cdf, angle, sig_born, sig_nrad, sig_rad, v_cdf[MERAD_NV], v_val[MERAD_NV];};
    std::vector<CDF_Angle> angle_dist(abins + 1);

    // first point
    CDF_Angle &fp = angle_dist[0];
    fp.angle = min_angle;
    fp.cdf = 0.;
    GetXS(Es, fp.angle, fp.sig_born, fp.sig_nrad, fp.sig_rad);
    for(int i = 1; i <= abins; ++i)
    {
        CDF_Angle &point = angle_dist[i];
        CDF_Angle &prev = angle_dist[i - 1];
        point.angle = min_angle + angle_step*i;
        // calculate cross sections
        GetXS(Es, point.angle, point.sig_born, point.sig_nrad, point.sig_rad);

        point.cdf = prev.cdf + angle_step*(point.sig_nrad + point.sig_rad + prev.sig_nrad + prev.sig_rad)/2.;

        // copy v distribution that calculated in MERADGEN
        for(int j = 0; j < MERAD_NV; ++j)
        {
            point.v_cdf[j] = merad_dist_.distsiv[j];
            point.v_val[j] = merad_dist_.distarv[j];
        }
    }

    // convert uniform distribution to cross-section vs. angle distribution
    auto comp = [](const CDF_Angle &point, const double &val)
                {
                    return point.cdf - val;
                };

    double beta_CM = get_CM_beta(Es);
    double angle, s, t, u, v, t1, z;
    double k2[4], p2[4], k[4], k2_CM[4], p2_CM[4], k_CM[4];

    for(int i = 0; i < nevents; ++i)
    {
        double rnd = uni_dist(rng)*angle_dist.back().cdf;
        auto interval = cana::binary_search_interval(angle_dist.begin(), angle_dist.end(), rnd, comp);

        // should not happen
        if(interval.first == angle_dist.end() || interval.second == angle_dist.end()) {
            std::cerr << "Could not find CDF value at " << rnd << std::endl;
            i--;
            continue;
        }

        // check if this event has a hard photon emission
        double rnd2 = uni_dist(rng);
        double sig_rad, sig_nrad, rnd_rad;

        // exact matched one point
        if(interval.first == interval.second) {
            angle = interval.first->angle;
            get_moller_stu(Es, angle, s, t, u);

            sig_rad = interval.first->sig_rad;
            sig_nrad = interval.first->sig_nrad;
            rnd_rad = (rnd2*(sig_nrad + sig_rad) - sig_nrad)/sig_rad;
            if(rnd_rad > 0.) {
                v = interp_dist(rnd_rad, interval.first->v_cdf, interval.first->v_val, MERAD_NV);
                t1 = merad_sample_t1(t, v, 0., uni_dist(rng));
                z = merad_sample_z(t, t1, v, 0., uni_dist(rng));
            } else {
                v = 0., t1 = t, z = 0.;
            }
        // in an interval, interpolate everything between two points
        } else {
            angle = cana::linear_interp(interval.first->cdf, interval.first->angle,
                                             interval.second->cdf, interval.second->angle,
                                             rnd);
            get_moller_stu(Es, angle, s, t, u);

            sig_nrad = cana::linear_interp(interval.first->cdf, interval.first->sig_nrad,
                                                  interval.second->cdf, interval.second->sig_nrad,
                                                  rnd);
            sig_rad = cana::linear_interp(interval.first->cdf, interval.first->sig_rad,
                                                 interval.second->cdf, interval.second->sig_rad,
                                                 rnd);
            rnd_rad = (rnd2*(sig_nrad + sig_rad) - sig_nrad)/sig_rad;
            if(rnd_rad > 0.) {
                double v1 = interp_dist(rnd_rad, interval.first->v_cdf, interval.first->v_val, MERAD_NV);
                double v2 = interp_dist(rnd_rad, interval.second->v_cdf, interval.second->v_val, MERAD_NV);

                v = cana::linear_interp(interval.first->cdf, v1,
                                        interval.second->cdf, v2,
                                        rnd);
                t1 = merad_sample_t1(t, v, 0., uni_dist(rng));
                z = merad_sample_z(t, t1, v, 0., uni_dist(rng));
            } else {
                v = 0., t1 = t, z = 0.;
            }
        }

        MomentumRec(k2_CM, p2_CM, k_CM, s, t, t1, v, z, uni_dist(rng)*2.*cana::pi);
        four_momentum_boost_z(k2, k2_CM, -beta_CM);
        four_momentum_boost_z(p2, p2_CM, -beta_CM);
        four_momentum_boost_z(k, k_CM, -beta_CM);

        std::cout << i << ", " << angle << ", " << k2[0] << ", " << p2[0] << ", " << k[0] << std::endl;
    }
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

    // Mandelstam variables
    double s, t, u0;
    get_moller_stu(Es, angle, s, t, u0);

    // virtual photon part of Moller cross section
    // the t and u channels should be calculated separately
    double sig_0t, sig_0u, sig_St, sig_Su, sig_vertt, sig_vertu, sig_Bt, sig_Bu;
    // t channel
    moller_Vph(s, t, sig_0t, sig_St, sig_vertt, sig_Bt);
    // u0 channel
    moller_Vph(s, u0, sig_0u, sig_Su, sig_vertu, sig_Bu);

    // infrared divergent part of real photon emission
    // NOTE that the t and u channels are not separated
    double delta_1H, delta_1S, delta_1inf;
    moller_IR(s, t, delta_1H, delta_1S, delta_1inf);

    // infrared free part of real photon emission
    double sig_Fs, sig_Fh;
    double v_limit = (s*t + sqrt(s*(s - 4.*m2)*t*(t - 4.*m2)))/2./m2;
    double v_max = (v_limit > v_cut) ? v_cut : v_limit;
    moller_IRF(s, t, v_min, v_max, sig_Fs, sig_Fh);

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
void PRadMollerGen::moller_Vph(double s, double t,
                               double &sig_0, double &sig_S, double &sig_vert, double &sig_B)
const
{
    double u0 = 4.*m2 - s - t;
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
    // we found that the equation in the paper is wrong in dimension
    // By comparing with the equation (8) in Ref.
    // N.M. Shumeiko and J.G. Suarez,
    // Journal of Physics G: Nuclear and Particle Physics 26, 2, 113 (2000).
    // There is a factor of (s - 2.0*m^2) missing due to misprint
    sig_0 = (u0*u0/xi_s2/4./s*(4.*xi_u04 - pow2(1. - xi_u02)*(2. + t/u0)) - s*s*xi_s4/u0)
            * 2.*cana::pi*alp2/t/t/s*(s - 2.*m2);

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
// input Mandelstam variables s, t
// NOTE that the t and u0 channel are not separated
// output 3 factorized part for the infrared part of the radiative cross section
void PRadMollerGen::moller_IR(double s, double t,
                              double &delta_1H, double &delta_1S, double &delta_1inf)
const
{
    double u0 = 4.*m2 - s - t;
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
    double v_max = (v_cut > v_limit) ? v_limit : v_cut;
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
                 double log_term = log(pow2((xi_ch + 1.)/(xi_ch - 1.)))
                                   * log(fabs((pow2(z_ch - 1.) - xi_ch2)/(1. - xi_ch2)));

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
                         // TODO, check with the authors
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

                         result += s_3/2./slamda_3*(/*term*/ + sum_term)*Sk[k];
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

    double s1 = -(xi_u02 + 1.)*u0;
    double s2 = (xi_s2 + 1.)*s;
    double s3 = -(xi_t2 + 1.)*t;

    std::cout << "S_phi test: "
              << s1 << ", " << s2 << ", " << s3 << ", "
              << S_phi(s1, s2, s3) << ", "
              << S_phi(s2, s1, s3)
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

// real photon emission (infrared free part) of the Moller scattering
// s, t: input, Mandelstam variables s, t in MeV^2
// v_min: input, the separation of "soft" and "hard" Bremsstrahlung
// v_max: input, the upper limit v of the Bremsstrahlung integration
// sig_Fs: output, "soft" Bremsstrahlung cross section
// sig_Fh: output, "hard" Bremsstrahlung cross section
void PRadMollerGen::moller_IRF(double s, double t, double v_min, double v_max,
                               double &sig_Fs, double &sig_Fh)
const
{
    // initialize MERADGEN
    merad_init(s);

    // the "soft" Bremsstrahlung part of the radiative cross section
    // blow v_min, photon emission is not detectable
    sig_Fs = merad_sigfs(v_min, t, 0.);
    // the "hard" Bremsstrahlung part of the radiative cross section
    sig_Fh = merad_sigfh(v_min, v_max, t, 0.);
}


// reconstruct the four momentum of outgoing particles by invariants and azimuthal angle
// units are in MeV and rad
// incident particles are in CM frame
// k1[0] = p1[0] = sqrt(s)/2, k1p = p1p = sqrt(lamda_s/s)/2
void PRadMollerGen::MomentumRec(double *k2, double *p2, double *k,
                                double s, double t, double t1, double v, double z, double phi)
{
    // frequently used variables
    double lamda_s = s*(s - 4.*m2);
    double lamda_1 = pow2(s - v) - 4.*s*m2;
    double lamda_2 = 2.*t + s - v - 4.*m2;
    double lamda_3 = -s*t*(s + t - v - 4.*m2) - pow2(m*v);
    double lamda_4 = s*(s - v - 4.*m2) - (s + v)*z;
    double lamda_5 = v*z*(s - v - z) - pow2(m*(v + z));
    double lamda_6 = s*(v - z) - v*(v + z);
    double lamda_7 = (s + 2.*t1 - z - 4.*m2)*lamda_1 - lamda_2*lamda_4;
    double lamda_8 = 16.*lamda_3*lamda_5 - lamda_7*lamda_7;

    double lamda_34 = lamda_3*lamda_4;
    double lamda_27 = lamda_2*lamda_7;
    double lamda_36 = lamda_3*lamda_6;
    double lamda_s18 = lamda_s*lamda_1*lamda_8;
    double lamda_1_s3 = lamda_1*sqrt(lamda_s*lamda_3);

    // outgoing electron 1
    k2[0] = (s - v)/sqrt(s)/2.;
    k2[1] = sqrt(lamda_3/lamda_s)*cos(phi);
    k2[2] = sqrt(lamda_3/lamda_s)*sin(phi);
    k2[3] = sqrt(s/lamda_s)*lamda_2/2.;

    // outgoing electron 2
    p2[0] = (s - z)/sqrt(s)/2.;
    p2[1] = -(sqrt(lamda_s18)*sin(phi) + (4.*lamda_34 + s*lamda_27)*cos(phi))/4./lamda_1_s3;
    p2[2] = (sqrt(lamda_s18)*cos(phi) - (4.*lamda_34 + s*lamda_27)*sin(phi))/4./lamda_1_s3;
    p2[3] = sqrt(s/lamda_s)*(lamda_7 - lamda_2*lamda_4)/lamda_1/2.;

    // outgoing photon
    k[0] = (v + z)/sqrt(s)/2.;
    k[1] = (sqrt(lamda_s18)*sin(phi) + (4.*lamda_36 - s*lamda_27)*cos(phi))/4./lamda_1_s3;
    k[2] = (-sqrt(lamda_s18)*cos(phi) + (4.*lamda_36 - s*lamda_27)*sin(phi))/4./lamda_1_s3;
    k[3] = sqrt(s/lamda_s)*(lamda_7 + lamda_2*lamda_6)/lamda_1/2.;
}

//============================================================================//
// Monte Carlo event generator for UNPOLARIZED Moller scattering beyond       //
// ultra-relativistic approximation (URA)                                     //
//----------------------------------------------------------------------------//
// References:                                                                //
// [1] Eur. Phys. J. A51, 1 (2015)                                            //
//     Radiative corrections beyond the ultra relativistic limit in           //
//     unpolarized ep elastic and Moller scatterings for the PRad Experiment  //
//     at Jefferson Laboratory                                                //
//     I. Akushevich, H. Gao, A. Ilyichev, and M. Meziane                     //
//                                                                            //
// [2] Comput. Phys. Commun. 176, 218 (2007)                                  //
//     MERADGEN 1.0: Monte Carlo generator for the simulation of radiative    //
//     events in parity conserving doubly-polarized Møller scattering         //
//     A. Afanasev, E. Chudakov, A. Ilyichev, V. Zykunov                      //
//                                                                            //
// [3] Journal of Physics G: Nuclear and Particle Physics 26, 2, 113 (2000)   //
//     The QED lowest-order radiative corrections to the two polarized        //
//     identical fermion scattering                                           //
//     N.M. Shumeiko and J.G. Suarez                                          //
//----------------------------------------------------------------------------//
// Code Developer: Chao Peng                                                  //
// Thanks to A. Ilyichev and C. Gu who helped review the formula              //
//============================================================================//


#include "PRadMollerGen.h"
#include "canalib.h"
#include <cmath>
#include <random>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "PRadBenchMark.h"

#define PROGRESS_EVENT_COUNT 1000
#define PROGRESS_BIN_COUNT 10

//#define MOLLER_TEST_MERA
//#define MOLLER_TEST_URA
//#define MOLLER_TEST_KIN

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

// calculate momentum square from four momentum
inline double calc_p2(double *p)
{
    return p[1]*p[1] + p[2]*p[2] + p[3]*p[3];
}

// calculate mass square from four momentum
inline double calc_mass2(double *p)
{
    return p[0]*p[0] - calc_p2(p);
}

// calculate the polar angle from four momentum
inline double calc_polar(double *p)
{
    return atan(sqrt((p[1]*p[1] + p[2]*p[2])/(p[3]*p[3])));
}

// calculate the azimuthal angle from four momentum
inline double calc_azimuthal(double *p)
{
    return atan2(p[1], p[2]);
}



// constructor
PRadMollerGen::PRadMollerGen(double vmin, double vmax, int bins)
: v_min(vmin), v_cut(vmax), theta_bins(bins)
{
    // place holder
}

// destructor
PRadMollerGen::~PRadMollerGen()
{
    // place holder
}

// generate Moller events, unit is in MeV, degree
// return integrated luminosity in the unit of nb^-1
double PRadMollerGen::Generate(double Es, double min_angle, double max_angle,
                               int nevents, const char *save_path, bool verbose)
const
{
    // sanity check
    if(min_angle > max_angle || nevents < 1 || Es < 0.) {
        std::cerr << "Invalid input(s), please make sure these inputs are correct:\n"
                  << "Es = " << Es << " MeV\n"
                  << "Angle range = " << min_angle << " ~ " << max_angle << "\n"
                  << "Number of events: " << nevents
                  << std::endl;
        return 0.;
    }

    PRadBenchMark timer;
    if(verbose) {
        std::cout << "Preparing distributions in theta bins to sample events..."
                  << std::endl;
    }

    // set up random number generator
    // a "true" random number engine, good for seed
    std::random_device rd;
    // Mersenne Twister pseudo-randome number generator
    std::mt19937 rng{rd()};
    // uniform distribution
    std::uniform_real_distribution<double> uni_dist(0., 1.);

    /* TODO, this method was trying to avoid duplicating events, but it fails
     * with radiation effects
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
    */

    // prepare grid for interpolation of angle
    double angle_step = (max_angle - min_angle)/(double)theta_bins;
    struct CDF_Angle {double cdf, angle, sig_born, sig_nrad, sig_rad, v_cdf[MERAD_NV], v_val[MERAD_NV];};
    std::vector<CDF_Angle> angle_dist(theta_bins + 1);

    // first point
    CDF_Angle &fp = angle_dist[0];
    fp.angle = min_angle;
    fp.cdf = 0.;
    GetXS(Es, fp.angle, fp.sig_born, fp.sig_nrad, fp.sig_rad);
    for(int i = 1; i <= theta_bins; ++i)
    {
        // show progress
        if(verbose && i%PROGRESS_BIN_COUNT == 0) {
            std::cout <<"------[ bin " << i << "/" << theta_bins << " ]---"
                      << "---[ " << timer.GetElapsedTimeStr() << " ]---"
                      << "---[ " << timer.GetElapsedTime()/(double)i
                      << " ms/bin ]------\r"
                      << std::flush;
        }

        CDF_Angle &point = angle_dist[i];
        CDF_Angle &prev = angle_dist[i - 1];
        point.angle = min_angle + angle_step*i;
        // calculate cross sections
        GetXS(Es, point.angle, point.sig_born, point.sig_nrad, point.sig_rad);

        // solid angle d_Omega = sin(theta)*d_theta*d_phi
        double curr_xs = 2.*cana::pi*sin(point.angle*cana::deg2rad)*(point.sig_nrad + point.sig_rad);
        double prev_xs = 2.*cana::pi*sin(prev.angle*cana::deg2rad)*(prev.sig_nrad + prev.sig_rad);

        // trapezoid rule for integration
        point.cdf = angle_step*cana::deg2rad*(curr_xs + prev_xs)/2. + prev.cdf;

        // copy v distribution that calculated in MERADGEN
        for(int j = 0; j < MERAD_NV; ++j)
        {
            point.v_cdf[j] = merad_dist_.distsiv[j];
            point.v_val[j] = merad_dist_.distarv[j];
        }
    }

    if(verbose) {
        std::cout <<"------[ bin " << theta_bins << "/" << theta_bins << " ]---"
                  << "---[ " << timer.GetElapsedTimeStr() << " ]---"
                  << "---[ " << timer.GetElapsedTime()/(double)theta_bins
                  << " ms/bin ]------\n"
                  << "Preparation done! Now start events generation..."
                  << std::endl;
    }
    timer.Reset();

    // prepare variables
    double beta_CM = get_CM_beta(Es);
    double angle, s, t, u, v, t1, z;
    double k2[4], p2[4], k[4], k2_CM[4], p2_CM[4], k_CM[4];
    std::ofstream fout(save_path);

    // convert uniform distribution to cross-section vs. angle distribution
    // lamda function for binary search
    auto comp = [](CDF_Angle point, double val) { return point.cdf - val;};

    // event sampling
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

        // exactly matched one point
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

        MomentumRec(k2_CM, p2_CM, k_CM, s, t, t1, v, z, uni_dist(rng), uni_dist(rng));
        four_momentum_boost_z(k2, k2_CM, -beta_CM);
        four_momentum_boost_z(p2, p2_CM, -beta_CM);
        four_momentum_boost_z(k, k_CM, -beta_CM);

        fout << sqrt(calc_p2(k2)) << "   " << calc_polar(k2) << "   " << calc_azimuthal(k2) << "   "
             << sqrt(calc_p2(p2)) << "   " << calc_polar(p2) << "   " << calc_azimuthal(p2) << "   ";
        if(rnd_rad > 0.) {
            fout << sqrt(calc_p2(k)) << "   " << calc_polar(k) << "   " << calc_azimuthal(k) << std::endl;
        } else {
            fout << 0. << "   " << 0. << "   " << 0. << std::endl;
        }

        // show progress
        if(verbose && i%PROGRESS_EVENT_COUNT == 0) {
            std::cout <<"------[ ev " << i << "/" << nevents << " ]---"
                      << "---[ " << timer.GetElapsedTimeStr() << " ]---"
                      << "---[ " << timer.GetElapsedTime()/(double)i
                      << " ms/ev ]------\r"
                      << std::flush;
        }

#ifdef MOLLER_TEST_KIN
        std::cout << i << ", " << angle << ", "
                  << k2[0] << ", " << calc_mass2(k2) << ", "
                  << p2[0] << ", " << calc_mass2(p2) << ", "
                  << k[0] << ", " << calc_mass2(k) << ", "
                  << z << ", " << 2.*(k[0]*k2[0] - k[1]*k2[1] - k[2]*k2[2] - k[3]*k2[3]) << ", "
                  << t1 << ", " << pow2(p2[0] - m) - p2[1]*p2[1] - p2[2]*p2[2] - p2[3]*p2[3] << ", "
                  << pow2(Es - k2[0] - k[0]) - pow2(k2[1] + k[1]) - pow2(k2[2] + k[2]) - pow2(sqrt(Es*Es - m*m) - k2[3] - k[3])
                  << std::endl;
        std::cout << sqrt(s) << ", "
                  << k2_CM[0] + p2_CM[0] + k_CM[0] << ", "
                  << k2_CM[1] + p2_CM[1] + k_CM[1] << ", "
                  << k2_CM[2] + p2_CM[2] + k_CM[2] << ", "
                  << k2_CM[3] + p2_CM[3] + k_CM[3]
                  << std::endl;
        std::cout << Es + m << ", "
                  << k2[0] + p2[0] + k[0] << ", "
                  << k2[1] + p2[1] + k[1] << ", "
                  << k2[2] + p2[2] + k[2] << ", "
                  << k2[3] + p2[3] + k[3]
                  << std::endl;
#endif //MOLLER_TEST_KIN
    }

    if(verbose) {
        std::cout <<"------[ ev " << nevents << "/" << nevents << " ]---"
                  << "---[ " << timer.GetElapsedTimeStr() << " ]---"
                  << "---[ " << timer.GetElapsedTime()/(double)nevents
                  << " ms/ev ]------\n"
                  << "Events generation done! Saved in file \"" << save_path << "\".\n"
                  << "Integrated luminosity = " << (double)nevents/angle_dist.back().cdf
                  << " nb^-1."
                  << std::endl;
    }

    return (double)nevents/angle_dist.back().cdf;
}

// get differential cross section dsigma/dOmega
// input beam energy (MeV), angle (deg)
// output Born, non-radiative, radiative cross sections (nb)
void PRadMollerGen::GetXS(double Es, double angle, double &sig_born, double &sig_nrad, double &sig_rad)
const
{
    double theta = angle*cana::deg2rad;

    // Mandelstam variables
    double s, t, u0;
    get_moller_stu(Es, angle, s, t, u0);
    // conversion from dsigma/dQsq (Q^2 = -t) to dsigma/dOmega
    double jacob = (s - 2.*m2)*4.*m*(Es - m)/pow2(Es + m -(Es - m)*pow2(cos(theta)))/2./cana::pi;

    GetXSdQsq(s, t, sig_born, sig_nrad, sig_rad);

    // convert to dsigma/dOmega, using unit nb
    sig_born *= jacob*unit;
    sig_nrad *= jacob*unit;
    sig_rad *= jacob*unit;
}

// get differential cross section dsigma/dQ2
// input Mandelstam variables s, t (MeV^2)
// output Born, non-radiative, radiative cross sections (MeV^-4)
void PRadMollerGen::GetXSdQsq(double s, double t, double &sig_born, double &sig_nrad, double &sig_rad)
const
{
    double u0 = 4.*m2 - s -t;
    // virtual photon part of Moller cross section
    // the t and u channels should be calculated separately
    double sig_0t, sig_0u, sig_St, sig_Su, sig_vertt, sig_vertu, sig_Bt, sig_Bu;
    // t channel
    SigmaVph(s, t, sig_0t, sig_St, sig_vertt, sig_Bt);
    // u0 channel
    SigmaVph(s, u0, sig_0u, sig_Su, sig_vertu, sig_Bu);

    // limitation on variable v, 99% of its allowed kinematic value
    double v_limit = 0.99*(s*t + sqrt(s*(s - 4.*m2)*t*(t - 4.*m2)))/2./m2;

    // infrared divergent part of real photon emission
    // NOTE that the t and u channels are not separated
    double v_ir = (v_min > v_limit) ? v_limit : v_min;
    double delta_1H, delta_1S, delta_1inf;
    SigmaIR(s, t, v_ir, delta_1H, delta_1S, delta_1inf);

    // infrared free part of real photon emission
    double sig_Fs, sig_Fh;
    double v_f = (v_cut > v_limit) ? v_limit : v_cut;
    SigmaF(s, t, v_ir, v_f, sig_Fs, sig_Fh);

    // t and u channels together
    // born level cross section
    sig_born = sig_0t + sig_0u;


    // non-radiative cross section
    double sig_S = sig_St + sig_Su;
    double sig_vert = sig_vertt + sig_vertu;
    double sig_B = sig_Bt + sig_Bu;
    double sig_IR = alp_pi*(delta_1H + delta_1S + delta_1inf)*sig_born;

    sig_nrad = sig_born + sig_IR + sig_S + sig_vert + sig_B + sig_Fs;

#ifdef MOLLER_TEST_MERA
    merad_init(s);
    double sig_vr = merad_sig(t, 0., 1);
    double sig_B2 = merad_sig(t, 0., 2);
    double sig_IR2 = merad_sigir(v_ir, t, 0.);

    double sig_nrad2 = sig_born + sig_IR2 + sig_vr + sig_B2 + sig_Fs;

    std::cout << "PRADMOLL: " << s << ", " << t << ", " << (sig_vert + sig_S)/sig_born << ", "
              << sig_B/sig_born << ", " << sig_IR/sig_born << ", "
              << sig_nrad/sig_born << std::endl;
    std::cout << "MERADGEN: " << s << ", " << t << ", " << sig_vr/sig_born << ", "
              << sig_B2/sig_born << ", " << sig_IR2/sig_born << ", "
              << sig_nrad2/sig_born << std::endl;

    sig_nrad = sig_nrad2;

#endif //MOLLER_TEST_MERA

    // radiative cross section
    sig_rad = sig_Fh;
}

// Cross section including virtual photon part for Moller scattering
// input Mandelstam variables s, t, u0
// output cross section including virtual photon effects
void PRadMollerGen::SigmaVph(double s, double t,
                             double &sig_0, double &sig_S, double &sig_vert, double &sig_B)
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

    // equation (49) in [1], Born Level
    sig_0 = (u0*u0/xi_s2/4./s*(4.*xi_u04 - pow2(1. - xi_u02)*(2. + t/u0)) - s*s*xi_s4/u0)
            * 2.*cana::pi*alp2/t/t/s;

    // singularity term, appears in delta_ver and delta_box
    // we only need the divergence free part of the sigma_ver and simga_box,
    // which are obtained by substituting lamda = m so log(lamda/m) = 0
    double log_m = 0.; // log(lamda/m), where lamda is the infinitesimal photon mass

    // other frequently used variables
    // Q^2 (-t) related, equation (27) - (29) in [1]
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
    // equation (41) in [1] with Q^2 -> -t
    double delta_vac = 0.;
    double lepton_mass[3] = {cana::ele_mass, cana::mu_mass, cana::tau_mass};
    for(auto &vac_m : lepton_mass)
    {
        double vac_m2 = vac_m*vac_m;
        double vac_slamda_m = sqrt(t*t - 4.*vac_m2*t);
        double vac_L_m = log((vac_slamda_m - t)/(vac_slamda_m + t))/vac_slamda_m;
        delta_vac += 2./3.*(-t + 2*vac_m2)*vac_L_m - 10./9. - 8./3.*vac_m2/t*(1. - 2.*vac_m2*vac_L_m);
    }

    // equation (50) in [1]
    sig_S = alp_pi*delta_vac*sig_0;

    // vertex correction, factorized part
    // equation (36) in [1] with Q^2 -> -t
    double delta_vert = 2.*(Q2_m*L_m - 1.)*log_m + (4.*m2 - 3./2.*t)*L_m - 2.
                        - Q2_m/slamda_m*(lamda_m*L_m*L_m/2. + 2.*cana::spence(2.*slamda_m/(slamda_m - t)) - pi2/2.);

    // vertex correction, non-factorized part, anomalous magnetic moment
    // equation (52) in [1]
    double sig_AMM = -4.*alp3/st2/xi_t*m2*log_t*(3.*(s - 2.*m2)/u0 + (10.*m2 - 3.*u0)/(s - 4.*m2));

    // equation (51) in [1]
    sig_vert = 2.*alp_pi*delta_vert*sig_0 + sig_AMM;

    // box diagram, factorized part
    // equation (54) in [1]
    double delta_box = (1. + xi_s2)/xi_s*(-4.*log_s*log_m + log_s*log_s - 2.*pi2 + 4.*Li2_sp)
                       + (1. + xi_u02)/xi_u0*(4.*log_u0*log_m - log_u0*log_u0 + 2.*log_2u0p*log_2u0p - pi2/3. + 4.*Li2_u0);

    // box diagram, non-factorized part
    // equation (A.1) and (A.2) in [1]
    // NOTICE there is a misprint in the (A.2) and we fixed it here
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
                            * (3.*log_u0*log_u0 - 6.*log_2u0*log_2u0 - 12.*Li2_u0 + 6.*log_u0*log(xi_u02*u0/t) + pi2)
             + 1./12./xi_t2/xi_t/t * (xi_t2*(-xi_t4 + 2.*xi_t2 + 3.)*t3
                                      + 2.*(xi_t4*xi_t2 - 4.*xi_t4 + 8.*xi_t2 + 1.)*u0t2
                                      + 4.*(3.*xi_t4 + 1.)*u02t + 16.*xi_t2*u03)
                                   * (-6.*log_2t*log_2t + 3.*log_t*log_t - 12.*Li2_t + 4.*pi2)
             + log_4t*(2.*u0/xi_t2*(xi_t2*t + t + 2.*u0) + (t - u0)*(2.*t + xi_u02*u0 + u0))
             - 1./xi_u0*log_u0*(xi_u02*(t - u0) - 2.*t)*(2.*t + xi_u02*u0 + u0);

    // equation (53) in [1]
    sig_B = alp_pi/2.*delta_box*sig_0 + alp3/xi_s2/s2t/u0*(sig_B1 + sig_B2);
}

// Infrared part of the Moller cross sections with real photon emission
// input Mandelstam variables s, t
// NOTE that the t and u0 channel are not separated
// output 3 factorized part for the infrared part of the radiative cross section
void PRadMollerGen::SigmaIR(double s, double t, double v_max,
                            double &delta_1H, double &delta_1S, double &delta_1inf)
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

    // equation (A.5) - (A.13) in [1]
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

    // equation (A.3) in [1]
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


    // equation (A.14) - (A.15) in [1]
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
                                     sum_term += sign*pow2(log(fabs(z_ud[k] - z[i])))/2.;
                                 } else {
                                     // the input to log may be negative and thus
                                     // the result may be a complex number
                                     // TODO check with authors for this part
                                     double spence_term = cana::spence((z_ud[k] - z[i])/(z[j] - z[i]));
                                     sum_term += sign*(log(fabs(z_ud[k] - z[i]))*log(fabs(z[i] - z[j])) - spence_term);
                                 }
                             }
                         }

                         result += s_3/2./slamda_3*(term + sum_term)*Sk[k];
                     }
                     return result;
                 };

    // equation (A.4) in [1]
    delta_1S = (xi_s2 + 1.)/2./xi_s*(log_s*log_s + log_s + cana::spence(4.*xi_s/pow2(xi_s + 1.)))
               - (xi_t2 + 1.)/2./xi_t*(log_t*log_t - log_t + cana::spence(4.*xi_t/pow2(xi_t + 1.)))
               - (xi_u02 + 1.)/2./xi_u0*(log_u0*log_u0 - log_u0 + cana::spence(4.*xi_u0/pow2(xi_u0 + 1.)))
               - S_phi(-(xi_u02 + 1.)*u0, (xi_s2 + 1.)*s, -(xi_t2 + 1.)*t)
               + S_phi(-(xi_u02 + 1.)*u0, -(xi_t2 + 1.)*t, (xi_s2 + 1.)*s)
               - S_phi(-(xi_t2 + 1.)*t, (xi_s2 + 1.)*s, -(xi_u02 + 1.)*u0) + 1.;


    // equation (60) in [1]
    double J_0 = -2.*((xi_s2 + 1.)/xi_s*log_s - (xi_t2 + 1.)/xi_t*log_t
                      - (xi_u02 + 1.)/xi_u0*log_u0 + 2.);
    // equation (66) in [1]
    delta_1inf = J_0*log(v_max/m2);

// test difference with URA
#ifdef MOLLER_TEST_URA
    // equation (61) in [1]
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

    // equation (62) in [1]
    double delta_1S_URA = 1. - (log(-t/m2) - 1.)*log(s*(s + t)/t/t)
                          + log(-t/m2)*(3. - 2.*log((s + t)/s))
                          - 5./2.*pow2(log(-t/m2)) - pow2(log((s + t)/s))/2. - pi2/3.;

    // equation (63) in [1]
    double J_0_URA = -4.*(1. + log(m2*s/t/u0));

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
void PRadMollerGen::SigmaF(double s, double t, double v_min, double v_max,
                           double &sig_Fs, double &sig_Fh)
{
    // initialize MERADGEN
    merad_init(s);

    // the "soft" Bremsstrahlung part of the radiative cross section
    // blow v_min, photon emission is not detectable
    sig_Fs = merad_sigfs(v_min, t, 0.);
    // the "hard" Bremsstrahlung part of the radiative cross section
    sig_Fh = merad_sigfh(v_min, v_max, t, 0.);
}


// reconstruct the four momentum of outgoing particles by invariants and
// azimuthal angle based on the Appendix A in ref. [2]
// units are in MeV
// incident particles are in CM frame
// k1[0] = p1[0] = sqrt(s)/2, k1p = p1p = sqrt(lamda_s/s)/2
// rnd1 is a random number (0, 1) to sample azimuthal angle phi
// rnd2 is a random number (0, 1) to switch sign
void PRadMollerGen::MomentumRec(double *k2, double *p2, double *k,
                                double s, double t, double t1, double v, double z,
                                double rnd1, double rnd2)
{
    double phi = rnd1*cana::pi*2.;

    // frequently used variables
    double lamda_1, lamda_2, lamda_3, lamda_4, lamda_5, lamda_6, lamda_7, lamda_8;

    double lamda_s = s*(s - 4.*m2);

    // when it is non-radiative, some of the lamda is exactly 0
    // but due to the precision of real number, it may get a very small negative
    // number in the end, which destroy the momentum reconstruction
    // thus separate two cases here
    if(v == 0.) {
        lamda_1 = s*s - 4.*s*m2;
        lamda_2 = 2.*t + s - 4.*m2;
        lamda_3 = -s*t*(s + t - 4.*m2);
        lamda_4 = s*(s - 4.*m2);
        lamda_5 = 0.;
        lamda_6 = 0.;
        lamda_7 = 0.;
        lamda_8 = 0.;
    } else {
        lamda_1 = pow2(s - v) - 4.*s*m2;
        lamda_2 = 2.*t + s - v - 4.*m2;
        lamda_3 = -s*t*(s + t - v - 4.*m2) - pow2(m*v);
        lamda_4 = s*(s - v - 4.*m2) - (s + v)*z;
        lamda_5 = v*z*(s - v - z) - pow2(m*(v + z));
        lamda_6 = s*(v - z) - v*(v + z);
        lamda_7 = (s + 2.*t1 - z - 4.*m2)*lamda_1 - lamda_2*lamda_4;
        lamda_8 = 16.*lamda_3*lamda_5 - lamda_7*lamda_7;
    }

    double lamda_34 = lamda_3*lamda_4;
    double lamda_27 = lamda_2*lamda_7;
    double lamda_36 = lamda_3*lamda_6;
    double lamda_1_s3 = lamda_1*sqrt(lamda_s*lamda_3);

    // NOTICE, this sign change is from MERADGEN's code but it is not mentioned
    // in ref. [2]
    // it probably is related to t/u channel photon emission
    double slamda_s18 = sqrt(lamda_s*lamda_1*lamda_8);
    if(rnd2 > 0.5)
        slamda_s18 *= -1.;
/*
    // outgoing electron 1
    k2[0] = (s - v)/sqrt(s)/2.;
    k2[1] = sqrt(lamda_3/lamda_s)*cos(phi);
    k2[2] = sqrt(lamda_3/lamda_s)*sin(phi);
    k2[3] = sqrt(s/lamda_s)*lamda_2/2.;

    // outgoing electron 2
    // NOTE that the formula in ref. [2] cannot satisfy energy and momentum
    // conservation, changing the sign before terms have lamda_7 will recover
    // the conservation laws, thanks to C. Gu who pointed out this
    p2[0] = (s - z)/sqrt(s)/2.;
    p2[1] = -(slamda_s18*sin(phi) + (4.*lamda_34 - s*lamda_27)*cos(phi))/4./lamda_1_s3;
    p2[2] = (slamda_s18*cos(phi) - (4.*lamda_34 - s*lamda_27)*sin(phi))/4./lamda_1_s3;
    p2[3] = sqrt(s/lamda_s)*(-lamda_7 - lamda_2*lamda_4)/lamda_1/2.;
*/
    // outgoing photon
    k[0] = (v + z)/sqrt(s)/2.;
    k[1] = (slamda_s18*sin(phi) + (4.*lamda_36 - s*lamda_27)*cos(phi))/4./lamda_1_s3;
    k[2] = (-slamda_s18*cos(phi) + (4.*lamda_36 - s*lamda_27)*sin(phi))/4./lamda_1_s3;
    k[3] = sqrt(s/lamda_s)*(lamda_7 + lamda_2*lamda_6)/lamda_1/2.;

    // p2 - p1
    double vp[4];
    vp[0] = -z/sqrt(s)/2.;
    vp[1] = -(slamda_s18*sin(phi) + (4.*lamda_34 - s*lamda_27)*cos(phi))/4./lamda_1_s3;
    vp[2] = (slamda_s18*cos(phi) - (4.*lamda_34 - s*lamda_27)*sin(phi))/4./lamda_1_s3;
    vp[3] = (lamda_s*lamda_1 - s*(lamda_7 + lamda_2*lamda_4))/2./sqrt(s*lamda_s)/lamda_1;

    // outgoing electron 1
    // p2 - p1 = k1 - k2 - k
    // k2 = k1 - k - vp
    k2[0] = sqrt(s)/2. - vp[0] - k[0];
    k2[1] = -vp[1] - k[1];
    k2[2] = -vp[2] - k[2];
    k2[3] = sqrt(lamda_s/s)/2. - vp[3] - k[3];

    // outgoing electron 2
    // p2 = vp + p1
    p2[0] = vp[0] + sqrt(s)/2.;
    p2[1] = vp[1];
    p2[2] = vp[2];
    p2[3] = vp[3] - sqrt(lamda_s/s)/2.;

}

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
#include <cmath>
#include <random>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include "PRadBenchMark.h"

#define PROGRESS_EVENT_COUNT 1000
#define PROGRESS_BIN_COUNT 10

//#define MOLLER_TEST_MERA
//#define MOLLER_TEST_URA
//#define MOLLER_TEST_KIN

//============================================================================//
// Helper constants, structures and functions                                 //
//============================================================================//

// some constant values to be used
const double m = cana::ele_mass;
const double m2 = m*m;
const double alp_pi = cana::alpha/cana::pi;
const double pi2 = cana::pi*cana::pi;
const double alp2 = cana::alpha*cana::alpha;
const double alp3 = alp2*cana::alpha;
// convert MeV^-2 to nbarn
const double unit = cana::hbarc2*1e7;

// random number generator
static cana::rand_gen<> rng;

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
    // (p1 + p2)/(gamma1*m1 + gamma2*m2) = sqrt(Es^2 - m^2)/(Es + m)
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

// progress output
void show_progress(const PRadBenchMark &timer, int count, int max, const char *str, bool end = false)
{
    std::cout << "------[ " << str << " " << count << "/" << max << " ]---"
              << "---[ " << timer.GetElapsedTimeStr() << " ]---";

    double avg = (count > 0) ? (double)timer.GetElapsedTime()/(double)count : 0.;
    std::cout << "---[ " << avg << " ms/" << str << " ]------";

    if(!end)
        std::cout << "\r" << std::flush;
    else
        std::cout << std::endl;
}

// refine t bin til the interpolation precision reaches required value
inline void refine_t_bin(const PRadMollerGen &model,
                         std::vector<TDist> &container, size_t i, size_t f,
                         double prec, double s)
{
    // safety check
    if(container.size() >= MAX_T_BINS)
        return;

    const TDist &beg = container.at(i), &end = container.at(f);
    double t = (beg.val + end.val)/2.;

    container.emplace_back(t, model.GetNonRadXSdQsq(s, t), model.GetRadVDist(s, t));
    const TDist &center = container.back();

    if(std::abs(1. - 2.*center.sig_nrad/(beg.sig_nrad + end.sig_nrad)) > prec ||
       std::abs(1. - 2.*center.sig_rad/(beg.sig_rad + end.sig_rad)) > prec)
    {
        refine_t_bin(model, container, i, container.size() - 1, prec, s);
        refine_t_bin(model, container, container.size() - 1, f, prec, s);
    }

}

// refine v bin til the interpolation precision reaches required value
inline void refine_v_bin(std::vector<VDist> &container, size_t i, size_t f,
                         double prec, double t)
{
    if(container.size() >= MAX_V_BINS)
        return;

    auto &beg = container.at(i), &end = container.at(f);
    double v = (beg.val + end.val)/2.;

    container.emplace_back(v, merad_sigfh(v, t, 0.));
    auto &center = container.back();

    if(std::abs(1. - 2.*center.sig/(beg.sig + end.sig)) > prec) {
        refine_v_bin(container, i, container.size() - 1, prec, t);
        refine_v_bin(container, container.size() - 1, f, prec, t);
    }
}



//============================================================================//
// Constructor, destructor                                                    //
//============================================================================//

// constructor
PRadMollerGen::PRadMollerGen(double vmin, double vmax, int nbins, double t_res, double v_res)
: v_min(vmin), v_cut(vmax), min_bins(nbins), t_prec(t_res), v_prec(v_res)
{
    // place holder
}

// destructor
PRadMollerGen::~PRadMollerGen()
{
    // place holder
}



//============================================================================//
// Public functions                                                           //
//============================================================================//

// generate Moller events, unit is in MeV, degree
// return integrated luminosity in the unit of nb^-1
double PRadMollerGen::Generate(double Es, double min_angle, double max_angle,
                               int nevents, const char *save_path, bool verbose)
const
{
    // sanity check
    if(min_angle > max_angle || nevents < 1 || Es < 0.) {
        std::cerr << "Invalid inputs, please make sure these inputs are correct:\n"
                  << "Es = " << Es << " MeV\n"
                  << "Angle range = " << min_angle << " ~ " << max_angle << "\n"
                  << "Number of events: " << nevents
                  << std::endl;
        return 0.;
    }

    // kinematic variables
    double s, t, u, v, t1, z, t_min, t_max;
    double beta_CM = get_CM_beta(Es);

    // determine t range
    get_moller_stu(Es, min_angle, s, t_min, u);
    get_moller_stu(Es, max_angle, s, t_max, u);

    // prepare grid for interpolation of angle
    std::vector<TDist> t_dist = init_grids(s, t_min, t_max, verbose);

    // prepare variables
    double k2[4], p2[4], k[4], k2_CM[4], p2_CM[4], k_CM[4];
    std::ofstream fout(save_path);

    // event sampling
    PRadBenchMark timer;
    int iev = 0;
    while(iev < nevents)
    {
        double rnd = rng()*t_dist.back().cdf;
        auto interval = cana::binary_search_interval(t_dist.begin(), t_dist.end(), rnd);
        const auto &fp = interval.first, &sp = interval.second;

        // should not happen
        if(fp == t_dist.end() || sp == t_dist.end()) {
            std::cerr << "Could not find CDF value at " << rnd << std::endl;
            continue;
        }

        // check if this event has a hard photon emission
        double rnd2 = rng();
        double sig_rad, sig_nrad, rnd_rad;

        // exactly matched one point
        if(fp == sp) {
            t = fp->val;

            sig_rad = fp->sig_rad;
            sig_nrad = fp->sig_nrad;
            rnd_rad = (rnd2*(sig_nrad + sig_rad) - sig_nrad)/sig_rad;
            if(rnd_rad > 0.) {
                v = cana::uni2dist(fp->v_dist.begin(), fp->v_dist.end(), rnd*fp->v_dist.back().cdf);
                t1 = merad_sample_t1(t, v, 0., rng());
                z = merad_sample_z(t, t1, v, 0., rng());
            } else {
                v = 0., t1 = t, z = 0.;
            }
        // in an interval, interpolate everything between two points
        } else {
            t = cana::linear_interp(fp->cdf, fp->val, sp->cdf, sp->val, rnd);
            sig_nrad = cana::linear_interp(fp->cdf, fp->sig_nrad, sp->cdf, sp->sig_nrad, rnd);
            sig_rad = cana::linear_interp(fp->cdf, fp->sig_rad, sp->cdf, sp->sig_rad, rnd);

            rnd_rad = (rnd2*(sig_nrad + sig_rad) - sig_nrad)/sig_rad;
            if(rnd_rad > 0.) {
                double v1 = cana::uni2dist(fp->v_dist.begin(),
                                           fp->v_dist.end(),
                                           rnd_rad*fp->v_dist.back().cdf);
                double v2 = cana::uni2dist(sp->v_dist.begin(),
                                           sp->v_dist.end(),
                                           rnd_rad*sp->v_dist.back().cdf);

                v = cana::linear_interp(fp->cdf, v1, sp->cdf, v2, rnd);
                t1 = merad_sample_t1(t, v, 0., rng());
                z = merad_sample_z(t, t1, v, 0., rng());
            } else {
                v = 0., t1 = t, z = 0.;
            }
        }

        if(!MomentumRec(k2_CM, p2_CM, k_CM, s, t, t1, v, z, rng(), rng())) {
            continue;
        }
        iev ++;
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
        if(verbose && (iev % PROGRESS_EVENT_COUNT == 0))
            show_progress(timer, iev, nevents, "ev");

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
        show_progress(timer, nevents, nevents, "ev", true);
        std::cout << "Events generation done! Saved in file \"" << save_path << "\".\n"
                  << "Integrated luminosity = " << (double)nevents/(unit*t_dist.back().cdf)
                  << " nb^-1."
                  << std::endl;
    }

    // return integrated luminosity
    return (double)nevents/(t_dist.back().cdf*unit);
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
    double u0 = 4.*m2 - s - t;
    // virtual photon part of Moller cross section
    // the t and u channels should be calculated separately
    double sig_0t, sig_0u, sig_St, sig_Su, sig_vertt, sig_vertu, sig_Bt, sig_Bu;
    // t channel
    SigmaVph(s, t, sig_0t, sig_St, sig_vertt, sig_Bt);
    // u0 channel
    SigmaVph(s, u0, sig_0u, sig_Su, sig_vertu, sig_Bu);

    // limitation on variable v, 99% of its allowed kinematic value
    double v_limit = 0.99*(s*t + sqrt(s*(s - 4.*m2)*t*(t - 4.*m2)))/2./m2;

    // t and u channels together
    // born level cross section
    sig_born = sig_0t + sig_0u;

    // non-radiative cross section
    double sig_S = sig_St + sig_Su;
    double sig_vert = sig_vertt + sig_vertu;
    double sig_B = sig_Bt + sig_Bu;

    // v1 = vmin separates soft and hard photon
    // v2 = vcut is the cut end of photons
    double v2 = cana::clamp(v_cut, v_cut, v_limit);
    double v1 = cana::clamp(v_min, v_min, v2);

    // infrared divergent part of real photon emission
    // NOTE that the t and u channels are not separated
    double sig_IR = SigmaIR(s, t, v1);

    // infrared free part of real photon emission, "soft" part
    double sig_Fs = SigmaFs(s, t, v1);

    sig_nrad = sig_born + sig_IR + sig_S + sig_vert + sig_B + sig_Fs;

    // radiative cross section, infrared free part of real photon emission, "hard" part
    sig_rad = SigmaRad(s, t, v1, v2);

#ifdef MOLLER_TEST_MERA
    merad_init(s);
    double sig_vr2 = merad_sig(t, 0., 1);
    double sig_B2 = merad_sig(t, 0., 2);
    double sig_IR2 = merad_sigir(v_ir, t, 0.);

    double sig_nrad2 = sig_born + sig_IR2 + sig_vr2 + sig_B2 + sig_Fs;

    double sig_vr3 = merad_sig2(t, 1);
    double sig_B3 = merad_sig2(t, 2);
    double sig_IR3 = merad_sigir2(v_ir, t);

    double sig_nrad3 = sig_born + sig_IR3 + sig_vr3 + sig_B3 + sig_Fs;

    std::cout << "KINEMATICS: " << s << ", " << t << std::endl;

    std::cout << "PRADMOLL: " << (sig_vert + sig_S)/sig_born << ", "
              << sig_B/sig_born << ", " << sig_IR/sig_born << ", "
              << sig_nrad/sig_born << std::endl;

    std::cout << "MERADGEN: " << sig_vr2/sig_born << ", "
              << sig_B2/sig_born << ", " << sig_IR2/sig_born << ", "
              << sig_nrad2/sig_born << std::endl;

    std::cout << "MEHDIGEN: " << sig_vr3/sig_born << ", "
              << sig_B3/sig_born << ", " << sig_IR3/sig_born << ", "
              << sig_nrad3/sig_born << std::endl;

#endif //MOLLER_TEST_MERA
}

// similar to GetXSdQsq, but only calculate the non-radiative part
double PRadMollerGen::GetNonRadXSdQsq(double s, double t)
const
{
    double u0 = 4.*m2 - s - t;
    // virtual photon part of Moller cross section
    // the t and u channels should be calculated separately
    double sig_0t, sig_0u, sig_St, sig_Su, sig_vertt, sig_vertu, sig_Bt, sig_Bu;
    // t channel
    SigmaVph(s, t, sig_0t, sig_St, sig_vertt, sig_Bt);
    // u0 channel
    SigmaVph(s, u0, sig_0u, sig_Su, sig_vertu, sig_Bu);

    // t and u channels together
    // born level cross section
    double sig_born = sig_0t + sig_0u;

    // non-radiative cross section
    double sig_S = sig_St + sig_Su;
    double sig_vert = sig_vertt + sig_vertu;
    double sig_B = sig_Bt + sig_Bu;

    double v1 = cana::clamp(v_min, v_min, v_cut);

    // infrared divergent part of real photon emission
    // NOTE that the t and u channels are not separated
    double sig_IR = SigmaIR(s, t, v1);

    // infrared free part of real photon emission, "soft" part
    double sig_Fs = SigmaFs(s, t, v1);

    return sig_born + sig_IR + sig_S + sig_vert + sig_B + sig_Fs;
}

// returns the whole v distribution of radiative part in val-cdf way
std::vector<VDist> PRadMollerGen::GetRadVDist(double s, double t)
const
{
    std::vector<VDist> res;
    res.reserve(MAX_V_BINS);

    // initialize MERADGEN
    merad_init(s);

    // limitation on variable v, 99% of its allowed kinematic value
    double v_limit = 0.99*(s*t + sqrt(s*(s - 4.*m2)*t*(t - 4.*m2)))/2./m2;

    // integration range
    double v2 = cana::clamp(v_cut, v_cut, v_limit);
    double v1 = cana::clamp(v_min, v_min, v2);

    // set initial v bins
    double v_step = (v2 - v1)/(double)min_bins;
    for(size_t i = 0; i <= min_bins; ++i)
    {
        double v = v1 + i*v_step;
        res.emplace_back(v, merad_sigfh(v, t, 0.));
    }

    // refine v bins to reach required precision
    // the content in the loop will change container's size, so a fix number needs
    // to be used in for loop
    for(size_t i = 1; i <= min_bins; ++i)
    {
        refine_v_bin(res, i-1, i, v_prec, t);
    }

    // sort in v transcendent
    std::sort(res.begin(), res.end(), [] (const VDist &b1, const VDist &b2)
                                         {
                                             return b1.val < b2.val;
                                         });

    for(size_t i = 1; i < res.size(); ++i)
    {
        auto &prev = res.at(i - 1);
        auto &curr = res.at(i);
        curr.cdf = prev.cdf + (curr.val - prev.val)*(prev.sig + curr.sig)/2.;
    }

    res.shrink_to_fit();

    return res;
}



//============================================================================//
// Static functions                                                           //
//============================================================================//

// Cross section at Born level in t and u channel
// input Mandelstam variables s, t
double PRadMollerGen::SigmaBorn(double s, double t)
{
    double u0 = 4.*m2 - s - t;
    // frequently used variables
    double xi_s = sqrt(1. - 4.*m2/s);
    double xi_t = sqrt(1 - 4.*m2/t);
    double xi_u0 = sqrt(1. - 4.*m2/u0);
    double xi_s2 = xi_s*xi_s, xi_s4 = xi_s2*xi_s2;
    double xi_t2 = xi_t*xi_t, xi_t4 = xi_t2*xi_t2;
    double xi_u02 = xi_u0*xi_u0, xi_u04 = xi_u02*xi_u02;

    // equation (49) in [1]
    // t channel
    double sig_0t = (u0*u0/xi_s2/4./s*(4.*xi_u04 - pow2(1. - xi_u02)*(2. + t/u0)) - s*s*xi_s4/u0)
                    * 2.*cana::pi*alp2/t/t/s;
    // u channel t <-> u0
    double sig_0u = (t*t/xi_s2/4./s*(4.*xi_t4 - pow2(1. - xi_t2)*(2. + u0/t)) - s*s*xi_s4/t)
                    * 2.*cana::pi*alp2/u0/u0/s;

    return sig_0t + sig_0u;
}

// Cross section including virtual photon part for Moller scattering
// t channel only
// input Mandelstam variables s, t
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
    // which are obtained by substituting lambda = m so log(lambda/m) = 0
    double log_m = 0.; // log(lambda/m), where lambda is the infinitesimal photon mass

    // other frequently used variables
    // Q^2 (-t) related, equation (27) - (29) in [1]
    double Q2_m = -t + 2.*m2;
    double lambda_m = t*t - 4.*m2*t;
    double slambda_m = sqrt(lambda_m);
    double L_m = 1./slambda_m*log((slambda_m - t)/(slambda_m + t));

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
    for(auto &vm : lepton_mass)
    {
        double vm2 = vm*vm;
        double vslambda_m = sqrt(t*t - 4.*vm2*t);
        double vL_m = log((vslambda_m - t)/(vslambda_m + t))/vslambda_m;
        delta_vac += 2./3.*(-t + 2*vm2)*vL_m - 10./9. - 8./3.*vm2/t*(1. - 2.*vm2*vL_m);
    }

    // equation (50) in [1]
    sig_S = alp_pi*delta_vac*sig_0;

    // vertex correction, factorized part
    // equation (36) in [1] with Q^2 -> -t
    double delta_vert = 2.*(Q2_m*L_m - 1.)*log_m + (4.*m2 - 3./2.*t)*L_m - 2.
                        - Q2_m/slambda_m*(lambda_m*L_m*L_m/2. + 2.*cana::spence(2.*slambda_m/(slambda_m - t)) - pi2/2.);

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
// t and u channel
// input Mandelstam variables s, t
// output the cross section, including  3 factorized parts
double PRadMollerGen::SigmaIR(double s, double t, double v_max)
{
    // 3 factors;
    double delta_1H, delta_1S, delta_1inf;

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
                     double lambda_1 = s_1*s_1 - 16.*m2*m2, slambda_1 = sqrt(lambda_1);
                     double lambda_2 = s_2*s_2 - 16.*m2*m2, slambda_2 = sqrt(lambda_2);
                     double lambda_3 = s_3*s_3 - 16.*m2*m2, slambda_3 = sqrt(lambda_3);
                     // z_u and z_d
                     double z_ud[2] = {slambda_1/slambda_2 - 1., (s_1*s_2 - 4.*m2*s_3)/lambda_2 - 1.};
                     // z_1, z_2, z_3, z_4
                     double z[4] = {1./slambda_2*(4.*m2*(s_3 - slambda_3)/(s_2 - slambda_2) - s_1 - slambda_2),
                                    1./slambda_2*(4.*m2*(s_3 + slambda_3)/(s_2 - slambda_2) - s_1 - slambda_2),
                                    1./slambda_2*(s_1 - slambda_2 - 4.*m2*(s_3 + slambda_3)/(s_2 + slambda_2)),
                                    1./slambda_2*(s_1 - slambda_2 - 4.*m2*(s_3 - slambda_3)/(s_2 + slambda_2))};
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
                         double term = log((s_2 - slambda_2)/(s_2 + slambda_2))
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

                         result += s_3/2./slambda_3*(term + sum_term)*Sk[k];
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

    return alp_pi*(delta_1H + delta_1S + delta_1inf)*SigmaBorn(s, t);
}

// real photon emission (infrared free part) of the Moller scattering
// t and u channel
// s, t: input, Mandelstam variables s, t in MeV^2
// sig_born : input, Born level cross section, must be both t and u channel
// v_min: input, the separation of "soft" and "hard" Bremsstrahlung
// v_max: input, the upper limit v of the Bremsstrahlung integration
// res: relatively precision in integration
// SigmaFs: output, "soft" Bremsstrahlung cross section
// SigmaRad: output, "hard" Bremsstrahlung cross section
double PRadMollerGen::SigmaFs(double s, double t, double v_min, double res)
{
    // initialize MERADGEN
    merad_init(s);

    // the "soft" Bremsstrahlung part of the radiative cross section
    // blow v_min, photon emission is not detectable
    return cana::simpson_prec(merad_sigfs, 1e-30, v_min, res, t, 0., SigmaBorn(s, t));
}

double PRadMollerGen::SigmaRad(double s, double t, double v_min, double v_max, double res)
{
    // initialize MERADGEN
    merad_init(s);

    // the "hard" Bremsstrahlung part of the radiative cross section
    return cana::simpson_prec(merad_sigfh, v_min, v_max, res, t, 0.);
}


// reconstruct the four momentum of outgoing particles by invariants and
// azimuthal angle based on the Appendix A in ref. [2]
// units are in MeV
// incident particles are in CM frame
// k1[0] = p1[0] = sqrt(s)/2, k1p = p1p = sqrt(lambda_s/s)/2
// rnd1 is a random number (0, 1) to sample azimuthal angle phi
// rnd2 is a random number (0, 1) to switch sign
bool PRadMollerGen::MomentumRec(double *k2, double *p2, double *k,
                                double s, double t, double t1, double v, double z,
                                double rnd1, double rnd2)
{
    double phi = rnd1*cana::pi*2.;

    // frequently used variables
    double lambda_1, lambda_2, lambda_3, lambda_4, lambda_5, lambda_6, lambda_7, lambda_8;

    double lambda_s = s*(s - 4.*m2);

    // when it is non-radiative, some of the lambda is exactly 0
    // but due to the precision of real number, it may get a very small negative
    // number in the end, which destroy the momentum reconstruction
    // thus separate two cases here
    if(v == 0.) {
        lambda_1 = s*s - 4.*s*m2;
        lambda_2 = 2.*t + s - 4.*m2;
        lambda_3 = -s*t*(s + t - 4.*m2);
        lambda_4 = s*(s - 4.*m2);
        lambda_5 = 0.;
        lambda_6 = 0.;
        lambda_7 = 0.;
        lambda_8 = 0.;
    } else {
        lambda_1 = pow2(s - v) - 4.*s*m2;
        lambda_2 = 2.*t + s - v - 4.*m2;
        lambda_3 = -s*t*(s + t - v - 4.*m2) - pow2(m*v);
        lambda_4 = s*(s - v - 4.*m2) - (s + v)*z;
        lambda_5 = v*z*(s - v - z) - pow2(m*(v + z));
        lambda_6 = s*(v - z) - v*(v + z);
        lambda_7 = (s + 2.*t1 - z - 4.*m2)*lambda_1 - lambda_2*lambda_4;
        lambda_8 = 16.*lambda_3*lambda_5 - lambda_7*lambda_7;
    }

    double lambda_34 = lambda_3*lambda_4;
    double lambda_27 = lambda_2*lambda_7;
    double lambda_36 = lambda_3*lambda_6;
    double lambda_1_s3 = lambda_1*sqrt(lambda_s*lambda_3);

    double slambda_s18 = sqrt(lambda_s*lambda_1*lambda_8);
    // TODO, sometimes lambda_8 < 0 (given a large lambda_7) and results in undefined slambda_s18
    // Need to investigate the reason
    if (std::isnan(slambda_s18)) {
        return false;
    }

    // NOTICE, this sign change is from MERADGEN's code but it is not mentioned
    // in ref. [2]
    // it probably is related to t/u channel photon emission
    if (rnd2 > 0.5)
        slambda_s18 *= -1.;
/*
    // outgoing electron 1
    k2[0] = (s - v)/sqrt(s)/2.;
    k2[1] = sqrt(lambda_3/lambda_s)*cos(phi);
    k2[2] = sqrt(lambda_3/lambda_s)*sin(phi);
    k2[3] = sqrt(s/lambda_s)*lambda_2/2.;

    // outgoing electron 2
    // NOTE that the formula in ref. [2] cannot satisfy energy and momentum
    // conservation, changing the sign before terms have lambda_7 will recover
    // the conservation laws, thanks to C. Gu who pointed out this
    p2[0] = (s - z)/sqrt(s)/2.;
    p2[1] = -(slambda_s18*sin(phi) + (4.*lambda_34 - s*lambda_27)*cos(phi))/4./lambda_1_s3;
    p2[2] = (slambda_s18*cos(phi) - (4.*lambda_34 - s*lambda_27)*sin(phi))/4./lambda_1_s3;
    p2[3] = sqrt(s/lambda_s)*(-lambda_7 - lambda_2*lambda_4)/lambda_1/2.;
*/
    // outgoing photon
    k[0] = (v + z)/sqrt(s)/2.;
    k[1] = (slambda_s18*sin(phi) + (4.*lambda_36 - s*lambda_27)*cos(phi))/4./lambda_1_s3;
    k[2] = (-slambda_s18*cos(phi) + (4.*lambda_36 - s*lambda_27)*sin(phi))/4./lambda_1_s3;
    k[3] = sqrt(s/lambda_s)*(lambda_7 + lambda_2*lambda_6)/lambda_1/2.;

    // p2 - p1
    double vp[4];
    vp[0] = -z/sqrt(s)/2.;
    vp[1] = -(slambda_s18*sin(phi) + (4.*lambda_34 - s*lambda_27)*cos(phi))/4./lambda_1_s3;
    vp[2] = (slambda_s18*cos(phi) - (4.*lambda_34 - s*lambda_27)*sin(phi))/4./lambda_1_s3;
    vp[3] = (lambda_s*lambda_1 - s*(lambda_7 + lambda_2*lambda_4))/2./sqrt(s*lambda_s)/lambda_1;

    // outgoing electron 1
    // p2 - p1 = k1 - k2 - k
    // k2 = k1 - k - vp
    k2[0] = sqrt(s)/2. - vp[0] - k[0];
    k2[1] = -vp[1] - k[1];
    k2[2] = -vp[2] - k[2];
    k2[3] = sqrt(lambda_s/s)/2. - vp[3] - k[3];

    // outgoing electron 2
    // p2 = vp + p1
    p2[0] = vp[0] + sqrt(s)/2.;
    p2[1] = vp[1];
    p2[2] = vp[2];
    p2[3] = vp[3] - sqrt(lambda_s/s)/2.;

    return true;
}



//============================================================================//
// Private Function                                                           //
//============================================================================//

// initialize theta grids for events generation
std::vector<TDist> PRadMollerGen::init_grids(double s, double t_min, double t_max, bool verbose)
const
{
    std::vector<TDist> res;
    res.reserve(MAX_T_BINS);
    PRadBenchMark timer;

    if(verbose) {
        std::cout << "Initializing grids in theta to sample events..."
                  << std::endl;
    }

    double t_step = (t_max - t_min)/(double)min_bins;

    for(unsigned int i = 0; i <= min_bins; ++i)
    {
        // new point
        double t = t_min + t_step*i;
        res.emplace_back(t, GetNonRadXSdQsq(s, t), GetRadVDist(s, t));
        if(verbose) show_progress(timer, i, min_bins, "bin");
    }

    if(verbose) {
        show_progress(timer, min_bins, min_bins, "bin", true);
        std::cout << "Initialization done! Now refine binning to reach precision "
                  << "t: " << t_prec << ", v: " << v_prec <<  std::endl;
    }

    timer.Reset();

    // refine t bin
    // the content in the loop will change container's size, so a fix number needs
    // to be used in for loop
    for(unsigned int i = 1; i <= min_bins; ++i)
    {
        refine_t_bin(*this, res, i - 1, i, t_prec, s);
        if(verbose) show_progress(timer, i, min_bins, "bin");
    }

    // sort in Q2 transcendent (t descendant) order
    std::sort(res.begin(), res.end(), [] (const TDist &b1, const TDist &b2)
                                         {
                                             return b1.val > b2.val;
                                         });

    size_t tot_vbins = res.front().v_dist.size();
    // calculate cdf for each bin
    for(auto curr = res.begin(), prev = curr++;
        curr != res.end();
        curr++, prev++)
    {
        tot_vbins += curr->v_dist.size();
        double curr_xs = (curr->sig_nrad + curr->sig_rad);
        double prev_xs = (prev->sig_nrad + prev->sig_rad);

        // trapezoid rule for integration, Q2 = -t
        curr->cdf = prev->cdf + (prev->val - curr->val)*(curr_xs + prev_xs)/2.;
    }

    if(verbose) {
        show_progress(timer, min_bins, min_bins, "bin", true);
        std::cout << "Interpolation grids finalized! \n"
                  << "Total number of grids t: " << res.size() - 1 << ", v: "
                  << tot_vbins << std::endl;
    }

    res.shrink_to_fit();

    return res;
}

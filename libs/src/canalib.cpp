//============================================================================//
// Some functions that might be useful for data analysis                      //
//                                                                            //
// Chao Peng                                                                  //
// 11/18/2016                                                                 //
//============================================================================//

#include "canalib.h"

#define __GAMMA_NUM 9
#define __GAMMA_G 7.0
static double __gamma_c[] = {0.9999999999998099,
                             6.7652036812188510E2,
                            -1.2591392167224028E3,
                             7.7132342877765313E2,
                            -1.7661502916214059E2,
                             1.2507343278686905E1,
                            -1.3857109526572012E-1,
                             9.9843695780195716E-6,
                             1.5056327351493116E-7};

#define __SPENCE_NUM 8
#define __SPENCE_NMAX 50
#define __SPENCE_NMAX_PREC 1000
static double __spence_c[] = {-1.1741940560772957946600E-1,
                              -2.7618966846029390643791E-2,
                              -8.0493987190845793511240E-3,
                              -2.7095568666150792944136E-3,
                              -7.1455906877666711465857E-4,
                               4.1757495974272487715417E-2,
                              -4.9028996486663818655E-2,
                              -9.08073640732783360};


// sigmoid function
double cana::sigmoid(const double &a, const double &p)
{
	return 1./(1. + std::exp(-a/p));
}

// gamma function
double cana::gamma(const double &z)
{

    if(z < 1) {
        return pi*(1 - z)/sin(pi*(1 - z))/cana::gamma(2 - z);
    } else if(z == 1) {
        return 1;
    } else {
        double ag = __gamma_c[0];
        for(int k = 1; k < __GAMMA_NUM; ++k)
            ag += __gamma_c[k]/(z - 1. + k);

        double output = 0.5*log(2*pi)
                        + (z - 0.5)*log(z - 0.5 + __GAMMA_G)
                        - (z - 0.5 + __GAMMA_G)
                        + log(ag);
        return exp(output);
    }
}

// spence function
double cana::spence(const double &z, const double &res)
{
    if(z > 1) {
        return 2*pi*pi/6 - log(z)*log(z)/2 - cana::spence(1/z, res);
    } else if (z == 1) {
        return pi*pi/6;
    } else if (z > 0.5) {
        return pi*pi/6 - log(1-z)*log(z) - cana::spence(1-z, res);
    } else if (z == 0.5) {
        return pi*pi/6/2 - log(0.5)*log(0.5)/2;
    } else if (z > 0) {
        return cana::spence_tr(z, res, __SPENCE_NMAX); // do nothing, fall into the bottom session
    } else if (z == 0) {
        return 0;
    } else if (z > -0.95) {
        return cana::spence_tr(z, res, __SPENCE_NMAX_PREC);
    } else if (z > -1.05) {
        // poly fit
        double output = 0;
        double dz = z + 1;
        for(int i = 0; i < __SPENCE_NUM; ++i)
            output += __spence_c[i]*std::pow(dz, i);
        return -(1 + output*dz*dz)*pi*pi/6/2 + dz*log(2);
    } else {
        return -pi*pi/6 - log(-z)*log(-z)/2 - cana::spence(1/z, res);
    }
}

// truncation of spence function
double cana::spence_tr(const double &z, const double &res, const int &nmax)
{
    // calculate spence until res is reached
    double output = 0.;
    int n = 0;
    while(++n <= nmax)
    {
        double nth = std::pow(z, n)/n/n;
        output += nth;

        if(std::abs(nth) < res)
            break;
    }

    return output;
}

// simpson integration
double cana::simpson(double begin, double end,
                     double (*f)(const double&), double step, int Nmin)
{
    int Nsteps = (end - begin)/step;
    int Nbins = std::max(Nmin, Nsteps)/2;
    double s = (end - begin)/(double)(2.*Nbins);

    // first bin
    double result = (*f)(begin) + 4.*(*f)(begin + s) + (*f)(end);
    double x = begin + 2.*s;
    int i = 1;
    while(i++ < Nbins)
    {
        result += 2.*(*f)(x) + 4.*(*f)(x + s);
        x += 2.*s;
    }
    return result*s/3.;
}


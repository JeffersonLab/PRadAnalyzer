//============================================================================//
// Some functions that might be useful for data analysis                      //
//                                                                            //
// Chao Peng                                                                  //
// 11/18/2016                                                                 //
//============================================================================//

#include "canalib.h"
#include "legendre_table.h"
#include <limits>



// sigmoid function
double cana::sigmoid(double a, double p)
{
	return 1./(1. + std::exp(-a/p));
}

// gamma function
double cana::gamma(double z)
{
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
        return std::exp(output);
    }
}

// Universal Landau distribution function
double cana::landau(double x)
{
    // value cut-off
    // below this value, the probability is known to be very small,
    // and the integration range and precision is not appropriate
    if(x < -3.6)
        return 0.;

    auto expr = [x](double t)
                {
                    return std::sin(pi*t)*std::exp(-t*(x + std::log(t)));
                };

    return simpson(expr, 1e-22, 500, 500000)/pi;
}

// a fit method to calculate landau distribution
// much faster than doing integration
// algorithm from Root(CERN), which is based on CERNLIB G110 denlan
double cana::landau_fit(double x)
{
    static double p1[5] = {0.4259894875,-0.1249762550, 0.03984243700, -0.006298287635,   0.001511162253};
    static double q1[5] = {1.0         ,-0.3388260629, 0.09594393323, -0.01608042283,    0.003778942063};

    static double p2[5] = {0.1788541609, 0.1173957403, 0.01488850518, -0.001394989411,   0.0001283617211};
    static double q2[5] = {1.0         , 0.7428795082, 0.3153932961,   0.06694219548,    0.008790609714};

    static double p3[5] = {0.1788544503, 0.09359161662,0.006325387654, 0.00006611667319,-0.000002031049101};
    static double q3[5] = {1.0         , 0.6097809921, 0.2560616665,   0.04746722384,    0.006957301675};

    static double p4[5] = {0.9874054407, 118.6723273,  849.2794360,   -743.7792444,      427.0262186};
    static double q4[5] = {1.0         , 106.8615961,  337.6496214,    2016.712389,      1597.063511};

    static double p5[5] = {1.003675074,  167.5702434,  4789.711289,    21217.86767,     -22324.94910};
    static double q5[5] = {1.0         , 156.9424537,  3745.310488,    9834.698876,      66924.28357};

    static double p6[5] = {1.000827619,  664.9143136,  62972.92665,    475554.6998,     -5743609.109};
    static double q6[5] = {1.0         , 651.4101098,  56974.73333,    165917.4725,     -2815759.939};

    static double a1[3] = {0.04166666667,-0.01996527778, 0.02709538966};

    static double a2[2] = {-1.845568670,-4.284640743};

    double u, ue, us, denlan;

    if (x < -5.5) {
        u = std::exp(x + 1.0);
        // cut-off
        if(u < 1e-10)
            return 0.0;
        ue = std::exp(-1./u);
        us = std::sqrt(u);
        denlan = 0.3989422803*(ue/us)*(1 + (a1[0] + (a1[1] + a1[2]*u)*u)*u);
    } else if(x < -1.) {
        u = std::exp(-x - 1.);
        denlan = std::exp(-u)*std::sqrt(u)*
                 (p1[0]+(p1[1]+(p1[2]+(p1[3]+p1[4]*x)*x)*x)*x)/
                 (q1[0]+(q1[1]+(q1[2]+(q1[3]+q1[4]*x)*x)*x)*x);
    } else if(x < 1.) {
        denlan = (p2[0]+(p2[1]+(p2[2]+(p2[3]+p2[4]*x)*x)*x)*x)/
                 (q2[0]+(q2[1]+(q2[2]+(q2[3]+q2[4]*x)*x)*x)*x);
    } else if(x < 5.) {
        denlan = (p3[0]+(p3[1]+(p3[2]+(p3[3]+p3[4]*x)*x)*x)*x)/
                 (q3[0]+(q3[1]+(q3[2]+(q3[3]+q3[4]*x)*x)*x)*x);
    } else if(x < 12.) {
        u = 1./x;
        denlan = u*u*(p4[0]+(p4[1]+(p4[2]+(p4[3]+p4[4]*u)*u)*u)*u)/
                 (q4[0]+(q4[1]+(q4[2]+(q4[3]+q4[4]*u)*u)*u)*u);
    } else if(x < 50.) {
        u = 1./x;
        denlan = u*u*(p5[0]+(p5[1]+(p5[2]+(p5[3]+p5[4]*u)*u)*u)*u)/
                 (q5[0]+(q5[1]+(q5[2]+(q5[3]+q5[4]*u)*u)*u)*u);
    } else if(x < 300.) {
        u = 1./x;
        denlan = u*u*(p6[0]+(p6[1]+(p6[2]+(p6[3]+p6[4]*u)*u)*u)*u)/
                 (q6[0]+(q6[1]+(q6[2]+(q6[3]+q6[4]*u)*u)*u)*u);
    } else {
        u = 1./(x - x*std::log(x)/(x + 1.));
        denlan = u*u*(1+(a2[0]+a2[1]*u)*u);
    }

    return denlan;
}

// Landau straggling probability function
double cana::landau_straggle(double x, double xi, double x0, bool fit)
{
    if(xi < 0.)
        return 0.;

    double lamda = (x - x0)/xi;

    if(fit)
        return landau_fit(lamda)/xi;
    else
        return landau(lamda)/xi;
}

// spence function
double cana::spence(double z, double res)
{
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

    if(z > 1) {
        return 2*pi*pi/6. - log(z)*log(z)/2. - cana::spence(1/z, res);
    } else if (z == 1) {
        return pi*pi/6;
    } else if (z > 0.5) {
        return pi*pi/6. - log(1-z)*log(z) - cana::spence(1-z, res);
    } else if (z == 0.5) {
        return pi*pi/6./2. - log(0.5)*log(0.5)/2.;
    } else if (z > 0) {
        return cana::spence_tr(z, res, __SPENCE_NMAX);
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
        return -(1 + output*dz*dz)*pi*pi/6./2. + dz*log(2);
    } else {
        return -pi*pi/6. - log(-z)*log(-z)/2. - cana::spence(1/z, res);
    }
}

// truncation of spence function
double cana::spence_tr(double z, double res, int nmax)
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

// calculate legendre polynomial from look-up table
// based on the code from http://www.holoborodko.com/pavel/?page_id=679
// Copyright (c)2007-2010 Pavel Holoborodko.
cana::legendre_nodes cana::calc_legendre_nodes(int n, double prec)
{
    cana::legendre_nodes res;
    res.order = n;
    // not valid
    if(n < 2) return res;

    // reserve space for solutions
	int m = (n + 1) >> 1;
    res.weights.resize(m);

	double x0,  x1,  dx;	// Abscissas
	double w0,  w1,  dw;	// Weights
	double P0, P_1, P_2;	// Legendre polynomial values
	double dpdx;			// Legendre polynomial derivative
	double t0, t1, t2, t3;

	t0 = 1. - (1. - 1./(double)n)/(8.*(double)(n*n));
	t1 = 1./(4.*(double)n + 2.);

	for (int i = 1; i <= m; i++)
	{
		// Find i-th root of Legendre polynomial

		// Initial guess
		x0 = cos(pi*(double)((i << 2) - 1)*t1)*t0;

		// Newton iterations, at least one
		int j = 0;
        w0 = 0.;
		dx = dw = std::numeric_limits<double>::max();

        // loop to get required precision
		do
        {
			// Compute Legendre polynomial value at x0
			P_1 = 1.0;
			P0  = x0;

			for (int k = 2; k <= n; k++)
			{
				P_2 = P_1;
				P_1 = P0;
				t2  = x0*P_1;

			    // Use look-up table for small n
                if(k < LEGENDRE_TABLE_SIZE) {
				    P0 = t2 + legendre_table[k]*(t2 - P_2);
                } else {
                    t3 = (double)(k - 1)/(double)k;
                    P0 = t2 + t3*(t2 - P_2);
                }
			}

			// Compute Legendre polynomial derivative at x0
			dpdx = ((x0*P0 - P_1)*(double)n)/(x0*x0 - 1.);

			// Newton step
			x1 = x0 - P0/dpdx;

			// Weight computing
			w1 = 2./((1. - x1*x1)*dpdx*dpdx);

			// Compute weight w0 on first iteration, needed for dw
			if (j == 0) w0 = 2./((1. - x0*x0)*dpdx*dpdx);

			dx = x0 - x1;
			dw = w0 - w1;

			x0 = x1;
			w0 = w1;
			j++;
		} while((std::abs(dx) > prec || std::abs(dw) > prec) && j < 1000);

        int index = (m - 1) - (i - 1);
        res.weights[index].x = x1;
        res.weights[index].w = w1;
	}

	return res;
}


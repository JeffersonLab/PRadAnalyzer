#include "cana_integrate.h"
#include "cana_utils.h"
#include "data/legendre_table.h"

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



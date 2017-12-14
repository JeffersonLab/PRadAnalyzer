#ifndef CANA_INTEGRATE_H
#define CANA_INTEGRATE_H

#include <cmath>
#include <vector>
#include <cstddef>

// limit of bins for simpson integration by precision
#define MAX_SIMPSON_BINS 50000
// limit of bin size for simpson method (relative diff. b/a - 1.)
#define MIN_SIMPSON_SIZE 1e-15

namespace cana
{
    // weight point for gauss legendre integration
    struct weight_point
    {
        double x, w;

        weight_point() : x(0.), w(0.) {}
        weight_point(double xi, double wi) : x(xi), w(wi) {}
    };
    struct legendre_nodes { int order; std::vector<weight_point> weights; };

    // legendre polynomial calculation to get the weight table
    legendre_nodes calc_legendre_nodes(int n, double prec = 1e-10);
    // gauss-legendre quadrature
    template<typename F, typename... Args>
    double gauss_quad(const legendre_nodes &ln, F &&f, double a, double b, Args&&... args)
    {
        // invalid legendre nodes
        if(ln.order < 2) return 0.;

        // convert integration range to [-1, 1]
        double A = (b - a)/2., B = (b + a)/2.;

        // different first point handling for even or odd case
        const auto &fp = ln.weights.front();
        double s = (ln.order&1)?
                   (fp.w*f(B, args...)) :
                   (fp.w*f(B + A*fp.x, args...) + fp.w*f(B - A*fp.x, args...));

        for(size_t i = 1; i < ln.weights.size(); ++i)
        {
            const auto &p = ln.weights.at(i);
            s += p.w*(f(B + A*p.x, args...) + f(B - A*p.x, args...));
        }

        return A*s;
    }

    template<class T, typename F, typename... Args>
    double gauss_quad(const legendre_nodes &ln, F (T::*f), T *t, double a, double b, Args&&... args)
    {
        // wrapper of a member function
        auto fn = [t, f] (double val, Args&&... args2)
                  {
                      return (t->*f)(val, args2...);
                  };

        return gauss_quad(ln, fn, a, b, args...);
    }


    // simpson integration
    template<typename F, typename... Args>
    double simpson(F &&f, double a, double b, int steps, Args&&... args)
    {
        double s = (b - a)/(double)steps;
        double res = f(a, args...) + f(b, args...) + 4.*f(a + s/2., args...);
        for(int i = 1; i < steps; ++i)
        {
            res += 4.*f(a + s*i + s/2., args...) + 2.*f(a + s*i, args...);
        }

        return s/6.*res;
    }

    // simpson integration for class member function
    template<class T, typename F, typename... Args>
    double simpson(F (T::*f), T *t, double a, double b, int steps, Args&&... args)
    {
        double s = (b - a)/(double)steps;
        double res = (t->*f)(a, args...)
                     + (t->*f)(b, args...)
                     + 4.*(t->*f)(a + s/2., args...);
        for(int i = 1; i < steps; ++i)
        {
            res += 4.*(t->*f)(a + s*i + s/2., args...)
                   + 2.*(t->*f)(a + s*i, args...);
        }

        return s/6.*res;
    }

    // helper function for simpson integration, keep refine binning if the
    // precision is not reached
    template<typename F, typename... Args>
    inline double simpson_prec_helper(F &&f, double a, double f_a, double b, double f_b, double prec, int &count, Args&&... args)
    {
        double c = (a + b)/2.;
        double f_c = f(c, args...);
        if((++count < MAX_SIMPSON_BINS) &&
           (f_c != 0.) &&
           (std::abs(a/c - 1.) > MIN_SIMPSON_SIZE) &&
           (std::abs(1. - (f_a + f_b)/2./f_c) > prec)) {

            return simpson_prec_helper(f, a, f_a, c, f_c, prec, count, args...)
                   + simpson_prec_helper(f, c, f_c, b, f_b, prec, count, args...);
        } else {
            return (b - a)*(f_a + f_b + f_c*4.)/6.;
        }
    }

    // simpson integration according to required precision
    template<typename F, typename... Args>
    double simpson_prec(F &&f, double a, double b, double prec, Args&&... args)
    {
        double f_a = f(a, args...);
        double f_b = f(b, args...);
        int count = 0;
        return simpson_prec_helper(f, a, f_a, b, f_b, prec, count, args...);
    }

    template<class T, typename F, typename... Args>
    double simpson_prec(F (T::*f), T *t, double a, double b, double prec, Args&&... args)
    {
        double f_a = (t->*f)(a, args...);
        double f_b = (t->*f)(b, args...);
        int count = 0;

        // wrapper member function
        auto fn = [t, f] (double val, Args&&... args2)
                  {
                      return (t->*f)(val, args2...);
                  };

        return simpson_prec_helper(fn, a, f_a, b, f_b, prec, count, args...);
    }
} // namespace cana

#endif // CANA_INTEGRATE_H

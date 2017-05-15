#ifndef C_ANA_LIB_H
#define C_ANA_LIB_H
// some useful tools for the whole library

#include <iterator>
#include <algorithm>
#include <utility>
#include <cmath>
#include <cstdlib>

namespace cana
{
    const static double alpha = 7.297352568E-3;     // 1./137.03599911
    const static double pi = 3.1415926535897932;    // pi
    const static double rad2deg = 57.2957795131;    // rad to degree
    const static double deg2rad = 0.01745329252;    // degree to rad
    const static double ele_mass = 0.510998918;     // MeV
    const static double mu_mass = 105.6583745;      // MeV
    const static double tau_mass = 1776.82;         // MeV
    const static double proton_mass = 938.272046;   // MeV
    const static double neutron_mass = 939.5654133; // MeV
    const static double hbarc = 197.326968;         // hbar*c (MeV*fm)
    const static double hbarc2 = 38937.9323;        // (hbar*c)^2 (MeV*fm)^2
    const static double amu = 931.494043;           // MeV per amu

    inline double sigmoid(double a, double p);
    double gamma(double z);
    double landau(double x);
    double landau_fit(double x);
    double landau_straggle(double x, double xi = 1., double x0 = 0., bool fit = true);
    double spence(double z, double res = 1e-15);
    inline double spence_tr(double z, double res, int nmax);

    template<typename T>
    inline T clamp(T val, T min, T max)
    {
        if(val < min) return min;
        if(val > max) return max;
        return val;
    }

    // simpson integration
    double simpson(double begin, double end, double (*f)(double), double step, int Nmin);

    // simpson for lamda expression
    template<class Lamda_func>
    double simpson(double begin, double end, int Nbins, Lamda_func func)
    {
        double s = (end - begin)/(double)(2.*Nbins);

        double result = func(begin) + 4.*func(begin + s) + func(end);
        double x = begin + 2.*s;
        int i = 1;
        while(i++ < Nbins)
        {
            result += 2.*func(x) + 4.*func(x + s);
            x += 2.*s;
        }

        return result*s/3.;
    }

    template<class T>
    double simpson(double begin, double end,
                   double (T::*f)(const double&), T *t, double step, int Nmin)
    {
        int Nsteps = (end - begin)/step;
        int Nbins = std::max(Nmin, Nsteps)/2;
        double s = (end - begin)/(double)(2.*Nbins);

        // first bin
        double result = (t->*f)(begin) + 4.*(t->*f)(begin + s) + (t->*f)(end);
        double x = begin + 2.*s;
        int i = 1;
        while(i++ < Nbins)
        {
            result += 2.*(t->*f)(x) + 4.*(t->*f)(x + s);
            x += 2.*s;
        }

        return result*s/3.;
    }

    template<class T, typename... Args>
    double simpson(double begin, double end, double step, int Nmin,
                   double (T::*f)(const double&, const Args& ...), T *t, const Args&... args)
    {
        int Nsteps = (end - begin)/step;
        int Nbins = std::max(Nmin, Nsteps)/2;
        double s = (end - begin)/(double)(2.*Nbins);

        // first bin
        double result =  (t->*f)(begin, args...)
                       + 4.*(t->*f)(begin + s, args...)
                       + (t->*f)(end, args...);
        double x = begin + 2.*s;
        int i = 1;
        while(i++ < Nbins)
        {
            result += 2.*(t->*f)(x, args...) + 4.*(t->*f)(x + s, args...);
            x += 2.*s;
        }

        return result*s/3.;
    }

    template<class T, typename... Args>
    double simpson(double begin, double end, double step, int Nmin,
                   double (T::*f) (double, Args ...) const, const T *t, Args... args)
    {
        int Nsteps = (end - begin)/step;
        int Nbins = std::max(Nmin, Nsteps)/2;
        double s = (end - begin)/(double)(2.*Nbins);

        // first bin
        double result =  (t->*f)(begin, args...)
                       + 4.*(t->*f)(begin + s, args...)
                       + (t->*f)(end, args...);
        double x = begin + 2.*s;
        int i = 1;
        while(i++ < Nbins)
        {
            result += 2.*(t->*f)(x, args...) + 4.*(t->*f)(x + s, args...);
            x += 2.*s;
        }

        return result*s/3.;
    }

    // the function is based on c++ source code
    // it adds permutation parity track
    template<class BidirIt>
    bool permutate(BidirIt first, BidirIt last, int &parity)
    {
        if (first == last) return false;
        BidirIt i = last;
        if (first == --i) return false;

        while (true) {
            BidirIt i1, i2;

            i1 = i;
            if (*--i < *i1) {
                i2 = last;
                while (!(*i < *--i2))
                    ;
                std::iter_swap(i, i2);
                std::reverse(i1, last);
                size_t swap = std::distance(i1, last)/2 + 1;
                // odd number of swaps
                if(swap&1)  parity *= -1;
                // even number of swaps, no change needed
                return true;
            }
            if (i == first) {
                std::reverse(first, last);
                size_t swap = std::distance(first, last)/2;
                // odd number of swaps
                if(swap&1)  parity *= -1;
                // even number of swaps, no change needed
                return false;
            }
        }
    }

    template<class RdmaccIt, typename T>
    RdmaccIt binary_search(RdmaccIt beg, RdmaccIt end, const T &val)
    {
        RdmaccIt not_found = end;

        RdmaccIt mid = beg + (end - beg)/2;
        while(mid != end && *mid != val)
        {
            if(*mid > val)
                end = mid;
            else
                beg = mid + 1;
            mid = beg + (end - beg)/2;
        }

        if(mid == end)
            return not_found;

        return mid;
    }

    // Comp(*it, val) output should be defined as
    // =0 : ==
    // >0 : >
    // <0 : <
    template<class RdmaccIt, typename T, class Comp>
    RdmaccIt binary_search(RdmaccIt beg, RdmaccIt end, const T &val, Comp comp)
    {
        RdmaccIt not_found = end;

        RdmaccIt mid = beg + (end - beg)/2;
        while(mid != end && comp(*mid,val) != 0)
        {
            if(comp(*mid, val) > 0)
                end = mid;
            else
                beg = mid + 1;
            mid = beg + (end - beg)/2;
        }

        if(mid == end)
            return not_found;

        return mid;
    }

    template<class RdmaccIt, typename T>
    std::pair<RdmaccIt, RdmaccIt> binary_search_interval(RdmaccIt beg,
                                                         RdmaccIt end,
                                                         const T &val)
    {
        RdmaccIt first = beg, last = end;
        RdmaccIt mid = beg + (end - beg)/2;
        while(mid != end)
        {
            if(*mid == val)
                return std::make_pair(mid, mid);
            else if(*mid > val)
                end = mid;
            else
                beg = mid + 1;
            mid = beg + (end - beg)/2;
        }

        if(mid == last)
            return std::make_pair(last, last);

        if(*mid < val) {
            return std::make_pair(mid, mid + 1);
        } else {
            if(mid == first)
                return std::make_pair(last, last);
            return std::make_pair(mid - 1, mid);
        }

    }

    template<class RdmaccIt, typename T, class Compare>
    std::pair<RdmaccIt, RdmaccIt> binary_search_interval(RdmaccIt beg,
                                                         RdmaccIt end,
                                                         const T &val,
                                                         Compare comp)
    {
        RdmaccIt first = beg, last = end;
        RdmaccIt mid = beg + (end - beg)/2;
        while(mid != end)
        {
            if(comp(*mid, val) == 0)
                return std::make_pair(mid, mid);
            else if(comp(*mid, val) > 0)
                end = mid;
            else
                beg = mid + 1;
            mid = beg + (end - beg)/2;
        }

        if(mid == last)
            return std::make_pair(last, last);

        if(comp(*mid, val) < 0) {
            return std::make_pair(mid, mid + 1);
        } else {
            if(mid == first)
                return std::make_pair(last, last);
            return std::make_pair(mid - 1, mid);
        }
    }

    template<class RdmaccIt, typename T>
    RdmaccIt binary_search_close_less(RdmaccIt beg, RdmaccIt end, const T &val)
    {
        if(beg == end)
            return end;
        if(*(end - 1) <= val)
            return end - 1;

        RdmaccIt first = beg, last = end;
        RdmaccIt mid = beg + (end - beg)/2;
        while(mid != end)
        {
            if(*mid == val)
                return mid;
            if(*mid > val)
                end = mid;
            else
                beg = mid + 1;
            mid = beg + (end - beg)/2;
        }

        if(*mid < val) {
            return mid;
        } else {
            if(mid == first)
                return last;
            return mid - 1;
        }
    }

    template<class Iter, typename T>
    bool is_in(const T &val, Iter beg, Iter end)
    {
        for(Iter it = beg; it != end; ++it)
        {
            if(*it == val)
                return true;
        }

        return false;
    }

    // a fast double to int conversion from Lua
    // magical number is 2^51 + 2^52, it pushes double to the range of 2^52 ~ 2^53,
    // within this range, the mantissa part is exactly the same as an integer
    // it will be safe to simply cast it to a 32-bit integer
    // Assumed little endian here, change i[0] to i[1] for big endian
    union i_cast {double d; int i[2];};
    inline int double2int(double d)
    {
        volatile union i_cast u;
        u.d = d + 6755399441055744.0;
        return u.i[0];
    }

    // Check if a point is in a polygon
    // Implementation of the method described by Dan Sunday
    // Check http://geomalgorithms.com/a03-_inclusion.html for details
    // helper function to determine if the point is at the left side of a line
    // > 0 : (x,y) is at left of the line by (x0, y0) and (x1, y1)
    // = 0 : on the line
    // < 0 : right of the line
    template<typename T>
    inline T is_left(const T &x, const T &y,
                     const T &x0, const T &y0, const T &x1, const T &y1)
    {
        return (x1 - x0) * (y - y0) - (x - x0) * (y1 - y0);
    }

    // input: A 2d point and the iterators of polygon vertices
    //        Vertices order determine how the polygon is reconstructed
    //        The polygon is assumed to be closed, the last boundary connects
    //        its first and last vertices
    // output: winding number, 0 means outside
    template<class Iter, typename T>
    int inside_polygon_2d(const T& point, Iter pv_beg, Iter pv_end)
    {
        int wn = 0;    // the  winding number counter

        // loop through all edges of the polygon
        Iter it_next = pv_beg;
        it_next++;
        for(Iter it = pv_beg; it != pv_end; it++, it_next++)
        {
            // last vertex, the boundary is formed between it and the first one
            if(it_next == pv_end)
                it_next = pv_beg;

            // start y <= point.y
            if(it->y <= point.y)
            {
                // upward crossing
                if(it_next->y > point.y)
                    // left of this boundary
                    if(is_left(point.x, point.y, it->x, it->y, it_next->x, it_next->y) > 0)
                        ++wn;
            }
            else
            {
                // downward crossing
                if(it_next->y <= point.y)
                    // right of this boundary
                    if(is_left(point.x, point.y, it->x, it->y, it_next->x, it_next->y) < 0)
                        --wn;
            }
        }

        // wn is positivie for counter clockwise and negative for clockwise
        return abs(wn);
    }
};

#endif

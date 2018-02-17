#ifndef CANA_UTILS_H
#define CANA_UTILS_H

#include <iterator>
#include <utility>
#include <cmath>
#include <algorithm>

namespace cana
{
    const static double alpha = 7.297352568E-3;     // 1./137.03599911
    const static double pi = 3.14159265358979323846;// pi
    const static double euler = 0.5772156649;       // Euler constant gamma
    const static double rad2deg = 57.2957795131;    // rad to degree
    const static double deg2rad = 0.01745329252;    // degree to rad
    const static double ele_mass = 0.510998918;     // MeV
    const static double mu_mass = 105.6583745;      // MeV
    const static double tau_mass = 1776.82;         // MeV
    const static double proton_mass = 938.272046;   // MeV
    const static double neutron_mass = 939.5654133; // MeV
    const static double hbarc = 197.32696979;       // hbar*c (MeV*fm)
    const static double hbarc2 = 38937.93664;       // (hbar*c)^2 (MeV*fm)^2
    const static double amu = 931.494043;           // MeV per amu

    // commonly used functions
    inline double sigmoid(double a, double p);
    double gamma(double z);
    double landau(double x);
    double landau_fit(double x);
    double landau_straggle(double x, double xi = 1., double x0 = 0., bool fit = true);
    double spence(double z, double res = 1e-15);
    inline double spence_tr(double z, double res, int nmax);

    // utils
    template<typename T> inline T pow2(T val) {return val*val;}
    template<typename T> inline T pow3(T val) {return val*val*val;}
    inline bool is_odd(int i) {return i&1;}

    // clamp values to be restricted inside [min, max]
    template<typename T>
    inline T clamp(T val, T min, T max)
    {
        if(val < min) return min;
        if(val > max) return max;
        return val;
    }

    // iteration to check if the value is included
    template<class Iter, typename T>
    inline bool is_in(Iter beg, Iter end, T val)
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

    // check the range
    template<typename T> inline void update_max(T &max, T val) {if(val > max) max = val;}
    template<typename T> inline void update_min(T &min, T val) {if(val < min) min = val;}
    template<class Iter, typename T>
    inline void calc_range(Iter beg, Iter end, T &min, T &max)
    {
        if(beg == end) return;
        min = *beg, max = *beg;

        for(Iter it = std::next(beg); it != end; ++it)
        {
            update_min(min, *it);
            update_max(max, *it);
        }
    }

    // iteratively solving function f(x) = 0
    // it is based on Newton's method, convergence rate is quadratic, but
    // convergence is not guaranteed, good initial value helps convergence
    template<typename T, class func>
    T solve_func(func &&f, T initial, T res, unsigned int max_iter = 500)
    {
        while(max_iter--)
        {
            T fval = f(initial);

            if(std::abs(fval) <= res) break;
            // Newton's method, x1 = x0 - f(x0)/f'(x0)
            // Approximately f'(x0) = (f(x0 + delta) - f(x0 - delta))/2delta
            initial += -fval/((f(initial + res) - f(initial - res))/2./res);
        }

        return initial;
    }

    // Bisection method, it is for monotonic function
    // linear rate to converge, convergence is guaranteed if the correct range is provided
    // the helper function
    template<typename T, class func>
    T solve_func2_helper(func &&f, T fa, T a, T fb, T b, T res, unsigned int iter)
    {
        T c = (a + b)/2., fc = f(c);
        if(std::abs(fc) <= res || iter == 0) return c;
        if(fc*fa < 0.)
            return solve_func2_helper(f, fa, a, fc, c, res, iter - 1);
        else
            return solve_func2_helper(f, fc, c, fb, b, res, iter - 1);
    }

    // main function for the Bisection method
    template<typename T, class func>
    T solve_func2(func &&f, T a, T b, T res, unsigned int max_iter = 500)
    {
        T fa = f(a), fb = f(b);
        if(std::abs(fa) <= res) return a;
        if(std::abs(fb) <= res) return b;
        return solve_func2_helper(f, fa, a, fb, b, res, max_iter);
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

    //  binary search, iterator can be randomly accessed
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
    // == : 0,  > : >0, < : <0
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

    // binary search in an interval
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

    // binary search in an interval
    // out put of Comp must be defined as
    // == : 0,  > : >0, < : <0
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

    // binary search, return the closest smaller value of the input if the same
    // value is not found
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

} // namespace cana

#endif // CANA_UTILS_H

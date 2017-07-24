#ifndef C_ANA_LIB_H
#define C_ANA_LIB_H
// some useful tools for the whole library

#include <iterator>
#include <algorithm>
#include <utility>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>
#include <functional>



// limit of bins for simpson integration by precision
#define MAX_SIMPSON_BINS 60000

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

    // clamp values to be restricted inside [min, max]
    template<typename T>
    inline T clamp(T val, T min, T max)
    {
        if(val < min) return min;
        if(val > max) return max;
        return val;
    }

    // linear interpolation of two points (x1, y1), (x2, y2)
    template<typename T>
    inline T linear_interp(T x1, T y1, T x2, T y2, T val)
    {
        if(x1 == x2)
            return y1;
        return ((val - x1)*y2 + (x2 - val)*y1)/(x2 - x1);
    }

    // linear interpolation of three points (x1, y1), (x2, y2), (x3, y3)
    template<typename T>
    inline T parabolic_interp(T x1, T y1, T x2, T y2, T x3, T y3, T val)
    {
        if(x1 == x2)
            return linear_interp(x1, y1, x3, y3, val);
        if(x2 == x3)
            return linear_interp(x1, y1, x2, y2, val);

        double x = val - x2, xp = x3 - x2, xm = x1 - x2;
        double a = y2;
        double c = ((y3 - y1)*(xp + xm)/(xp - xm) - (y3 + y1 - 2.*y2))/(xp*xm*2.);
        double b = (y3 - y1)/(xp - xm) - c*(xp + xm);
        return a + b*x + c*x*x;
    }

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
        // wrapper member function
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
        if(++count < MAX_SIMPSON_BINS && std::abs(1. - (f_a + f_b)/2./f_c) > prec) {
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

    // iteration to check if the value is included
    template<class Iter, typename T>
    inline bool is_in(const T &val, Iter beg, Iter end)
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

    // seeding random engine
    template<class T = std::mt19937, std::size_t N = T::state_size>
    auto seeded_random_engine() -> typename std::enable_if<!!N, T>::type
    {
        typename T::result_type random_data[N];
        std::random_device rd;
        std::generate(std::begin(random_data), std::end(random_data), std::ref(rd));
        std::seed_seq seeds(std::begin(random_data), std::end(random_data));
        T re(seeds);
        return re;
    }

    template<typename T = double, class Engine = std::mt19937>
    class rand_gen
    {
    public:
        // constructor
        rand_gen()
        : engine(seeded_random_engine<Engine>())
        {
            typename Engine::result_type range = engine.max() - engine.min();
            divisor = static_cast<T>(range) + 1;

            if(divisor <= 0) {
                std::cerr << "CRandom: Output type cannot handle engine range, "
                          << "undefined value may be returned."
                          << std::endl;
            }
        }

        // get a random number in (0, 1)
        T operator() ()
        {
            return static_cast<T>(engine() - engine.min())/divisor;
        }

        // get a random number in (min, max)
        T operator() (T min, T max)
        {
            return static_cast<T>(engine() - engine.min())/divisor*(max - min) + min;
        }

    private:
        Engine engine;
        T divisor;
    };

    // simple structure for converting uniform distribution to any distribution
    struct val_cdf
    {
        double val, cdf;

        // constructors
        val_cdf() {}
        val_cdf(double v, double c) : val(v), cdf(c) {}

        // for binary search
        bool operator ==(double v) const {return cdf == v;}
        bool operator <(double v) const {return cdf < v;}
        bool operator >(double v) const {return cdf > v;}
        bool operator != (double v) const {return cdf != v;}
    };

    // convert uniform distribution to another distribution according to the
    // val-cdf distribution provided
    template<class RdmIt, typename T>
    inline T uni2dist(RdmIt begin, RdmIt end, T val)
    {
        auto itv = cana::binary_search_interval(begin, end, val);

        // should not happen
        if(itv.first == end || itv.second == end)
        {
            std::cerr << "Fail to convert distribution from value = " << val << std::endl;
            return 0.;
        }

        if(itv.first == itv.second) {
            return itv.first->val;
        } else {
            return cana::linear_interp(itv.first->cdf, itv.first->val,
                                       itv.second->cdf, itv.second->val,
                                       val);
        }
    }

} // namespace cana

#endif

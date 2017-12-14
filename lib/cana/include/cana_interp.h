#ifndef CANA_INTERP_H
#define CANA_INTERP_H

#include "cana_utils.h"
#include <cmath>

namespace cana
{
    struct PointErr
    {
        double x, y, stat, syst;

        PointErr() {}
        PointErr(double xi, double yi, double st, double sy)
        : x(xi), y(yi), stat(st), syst(sy) {}

        // for binary search
        bool operator ==(double xi) const {return x == xi;}
        bool operator <(double xi) const {return x < xi;}
        bool operator >(double xi) const {return x > xi;}
        bool operator != (double xi) const {return x != xi;}
    };

    // interpolation for 1 dimensional array, cubic interpolation will be used
    // for most of the points, in the case of lacking points, quadratic or linear
    // interpolation will be used
    // data should be ordered in x
    bool interp_1d(const std::vector<PointErr> &data, PointErr &p, double res = 1e-3);


    // interpolations, no safety checks, need to make sure all the points are
    // different and ordered
    //
    // linear interpolation of two points (x1, y1), (x2, y2)
    //
    template<typename T>
    inline void linear_coeff(T x1, T x2, T val, T *c)
    {
        c[0] = (x2 - val)/(x2 - x1);
        c[1] = (val - x1)/(x2 - x1);
    }

    template<typename T>
    inline T linear_interp(T x1, T y1, T x2, T y2, T val)
    {
        T c[2];
        linear_coeff(x1, x2, val, c);
        return c[0]*y1 + c[1]*y2;
    }

    template<class Point>
    inline void linear_interp(Point p1, Point p2, Point &p)
    {
        decltype(p.x) c[2];
        linear_coeff(p1.x, p2.x, p.x, c);
        p.y = c[0]*p1.y + c[1]*p2.y;
    }

    // linear interpolation with error propagation
    template<class Point>
    inline void linear_interp_err(Point p1, Point p2, Point &p, bool err_cor = false)
    {
        decltype(p.x) c[2];
        linear_coeff(p1.x, p2.x, p.x, c);

        p.y = c[0]*p1.y + c[1]*p2.y;
        if(err_cor)
            p.err = c[0]*p1.err + c[1]*p2.err;
        else
            p.err = std::sqrt(pow2(c[0]*p1.err) + pow2(c[1]*p2.err));
    }

    //
    // quadratic interpolation of three points (x1, y1), (x2, y2), (x3, y3)
    //
    template<typename T>
    inline void quad_coeff(T x1, T x2, T x3, T val, T *c)
    {
        T dx1 = val - x1, dx2 = val - x2, dx3 = val - x3;
        T a12 = x1 - x2, a13 = x1 - x3, a23 = x2 - x3;
        c[0] = dx2*dx3/(a12*a13);
        c[1] = dx1*dx3/(-a23*a12);
        c[2] = dx1*dx2/(a13*a23);
    }

    template<typename T>
    inline T quad_interp(T x1, T y1, T x2, T y2, T x3, T y3, T val)
    {
        T c[3];
        quad_coeff(x1, x2, x3, val, c);
        return c[0]*y1 + c[1]*y2 + c[2]*y3;
    }

    template<class Point>
    inline void quad_interp(Point p1, Point p2, Point p3, Point &p)
    {
        decltype(p.x) c[3];
        quad_coeff(p1.x, p2.x, p3.x, p.x, c);
        p.y = c[0]*p1.y + c[1]*p2.y + c[2]*p3.y;
    }

    // quadratic interpolation with error propagation
    template<class Point>
    inline void quad_interp_err(Point p1, Point p2, Point p3, Point &p, bool err_corr = false)
    {
        decltype(p.x) c[3];
        quad_coeff(p1.x, p2.x, p3.x, p.x, c);
        p.y = c[0]*p1.y + c[1]*p2.y + c[2]*p3.y;
        if(err_corr)
            p.err = c[0]*p1.err + c[1]*p2.err + c[2]*p3.err;
        else
            p.err = std::sqrt(pow2(c[0]*p1.err) + pow2(c[1]*p2.err) + pow2(c[2]*p3.err));
    }

    //
    // cubic interpolation of four points
    //
    template<typename T>
    inline void cubic_coeff(T x1, T x2, T x3, T x4, T val, T *c)
    {
        T dx1 = val - x1, dx2 = val - x2, dx3 = val - x3, dx4 = val - x4;
        T a12 = x1 - x2, a13 = x1 - x3, a14 = x1 - x4;
        T a23 = x2 - x3, a24 = x2 - x4, a34 = x3 - x4;
        c[0] = dx2*dx3*dx4/(a12*a13*a14);
        c[1] = dx1*dx3*dx4/(-a12*a23*a24);
        c[2] = dx1*dx2*dx4/(a13*a23*a34);
        c[3] = dx1*dx2*dx3/(-a14*a24*a34);
    }

    template<typename T>
    inline T cubic_interp(T x1, T y1, T x2, T y2, T x3, T y3, T x4, T y4, T val)
    {
        T c[4];
        cubic_coeff(x1, x2, x3, x4, val, c);
        return c[0]*y1 + c[1]*y2 + c[2]*y3 + c[3]*y4;
    }

    template<class Point>
    inline void cubic_interp(Point p1, Point p2, Point p3, Point p4, Point &p)
    {
        decltype(p.x) c[4];
        cubic_coeff(p1.x, p2.x, p3.x, p4.x, p.x, c);
        p.y = c[0]*p1.y + c[1]*p2.y + c[2]*p3.y + c[3]*p4.y;
    }

    template<class Point>
    inline void cubic_interp_err(Point p1, Point p2, Point p3, Point p4, Point &p, bool err_corr = false)
    {
        decltype(p.x) c[4];
        cubic_coeff(p1.x, p2.x, p3.x, p4.x, p.x, c);
        p.y = c[0]*p1.y + c[1]*p2.y + c[2]*p3.y + c[3]*p4.y;
        if(err_corr)
            p.err = c[0]*p1.err + c[1]*p2.err + c[3]*p3.err + c[4]*p4.err;
        else
            p.err = std::sqrt(pow2(c[0]*p1.err) + pow2(c[1]*p2.err) + pow2(c[2]*p3.err) + pow2(c[3]*p4.err));
    }

    // nth order polynomial fit, n = number of points - 1
    template<typename T>
    inline T poly_interp(int Np, T *x, T *y, T val)
    {
        if(Np < 2) return 0.;

        T *dx = new T[Np];

        for(int i = 0; i < Np; ++i)
        {
            dx[i] = val - x[i];
        }

        T res = 0.;

        for(int i = 1; i < Np; ++i)
        {
            T ci = 1.;
            for(int j = 1; j < Np; ++j)
            {
                if(j != i) ci *= dx[j]/(x[i] - x[j]);
            }

            res += ci*y[i];
        }

        delete[] dx;

        return res;
    }

} // namespace cana

#endif // CANA_INTERP_H

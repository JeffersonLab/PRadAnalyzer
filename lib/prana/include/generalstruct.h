#ifndef GENERAL_STRUCT_H
#define GENERAL_STRUCT_H

#include <cmath>
#include <iostream>
// general structures that will be used among several classes

// geometry for a HyCal module
struct Geometry
{
    int type;
    double size_x, size_y, size_z;
    double x, y, z;

    Geometry()
    : type(-1), size_x(0), size_y(0), size_z(0), x(0), y(0), z(0)
    {};

    Geometry(int t, double sx, double sy, double sz,
             double pos_x, double pos_y, double pos_z)
    : type(t), size_x(sx), size_y(sy), size_z(sz), x(pos_x), y(pos_y), z(pos_z)
    {};
};

// layout for a HyCal module
struct Layout
{
    unsigned int flag;
    int sector, row, column;

    Layout() : flag(0), sector(-1), row(0), column(0)
    {};

    Layout(unsigned int f, int s, int r, int c)
    : flag(f), sector(s), row(r), column(c)
    {};
};

// 2D point
template<typename T>
struct Point2D
{
    T x, y;

    Point2D() : x(0.), y(0.) {}
    Point2D(T xi, T yi) : x(xi), y(yi) {}
    Point2D(const T *v) : x(v[0]), y(v[1]) {}

    T norm() const
    {
        return std::sqrt(x*x + y*y);
    }

    template<typename T2>
    T dist(const Point2D<T2> &rhs) const
    {
        return (*this - rhs).norm();
    }

    template<typename T2>
    T dot(const Point2D<T2> &rhs) const
    {
        return x*rhs.x + y*rhs.y;
    }

    inline Point2D<T> rotate(double a) const
    {
        return Point2D<T>(x*std::cos(a) + y*std::sin(a), -x*std::sin(a) + y*std::cos(a));
    }

    inline Point2D<T> rotate_inv(double a) const
    {
        return Point2D<T>(x*std::cos(a) - y*std::sin(a), x*std::sin(a) + y*std::cos(a));
    }

    Point2D<T> operator -() const
    {
        return Point2D<T>(-x, -y);
    }

    template<typename T2>
    bool operator ==(const Point2D<T2> &rhs) const
    {
        return (x == rhs.x) && (y == rhs.y);
    }

    template<typename T2>
    bool operator !=(const Point2D<T2> &rhs) const
    {
        return (x != rhs.x) || (y != rhs.y);
    }

    template<typename T2>
    Point2D<T> operator ()(const Point2D<T2> &v) const
    {
        return Point2D<T>(v.x, v.y);
    }

    template<typename T2>
    Point2D<T> operator +(const Point2D<T2> &rhs) const
    {
        return Point2D<T>(x + rhs.x, y + rhs.y);
    }
    template<typename T2>
    Point2D<T> &operator +=(const Point2D<T2> &rhs)
    {
        *this = *this + rhs;
        return *this;
    }

    template<typename T2>
    Point2D<T> operator -(const Point2D<T2> &rhs) const
    {
        return Point2D<T>(x - rhs.x, y - rhs.y);
    }
    template<typename T2>
    Point2D<T> &operator -=(const Point2D<T2> &rhs)
    {
        *this = *this - rhs;
        return *this;
    }

    template<typename T2>
    Point2D<T> operator *(const T2 &rhs) const
    {
        return Point2D<T>(x*rhs, y*rhs);
    }
    template<typename T2>
    Point2D<T> &operator *=(const T2 &rhs)
    {
        *this = *this * rhs;
        return *this;
    }

    template<typename T2>
    Point2D<T> operator /(const T2 &rhs) const
    {
        return Point2D<T>(x/rhs, y/rhs);
    }
    template<typename T2>
    Point2D<T> &operator /=(const T2 &rhs)
    {
        *this = *this / rhs;
        return *this;
    }
};

template<typename T1, typename T2>
Point2D<T2> operator *(const T1 &lhs, const Point2D<T2> &rhs)
{
    return rhs*lhs;
}

template<typename T>
std::ostream &operator <<(std::ostream &os, const Point2D<T> &p)
{
    os << "(" <<  p.x << ", " << p.y << ")";
    return os;
}


// 3D point
template<typename T>
struct Point3D
{
    T x, y, z;

    Point3D() : x(0.), y(0.), z(0.) {}
    Point3D(T xi, T yi, T zi) : x(xi), y(yi), z(zi) {}
    Point3D(const T *v) : x(v[0]), y(v[1]), z(v[2]) {}

    T norm() const
    {
        return std::sqrt(x*x + y*y + z*z);
    }

    template<typename T2>
    T dist(const Point3D<T2> &rhs) const
    {
        return (*this - rhs).norm();
    }

    template<typename T2>
    T dot(const Point3D<T2> &rhs) const
    {
        return x*rhs.x + y*rhs.y + z*rhs.z;
    }

    template<typename T2>
    Point3D<T> cross(const Point3D<T2> &rhs) const
    {
        return Point3D<T>(y*rhs.z - z*rhs.y, z*rhs.x - x*rhs.z, x*rhs.y - y*rhs.x);
    }

    template<typename T2>
    inline Point3D<T> rotate(const Point3D<T2> &rot) const
    {
        // Rxyz = RxRyRz
        double cx = std::cos(rot.x), sx = std::sin(rot.x);
        double cy = std::cos(rot.y), sy = std::sin(rot.y);
        double cz = std::cos(rot.z), sz = std::sin(rot.z);

        return Point3D<T>(cy*cz*x - cy*sz*y + sy*z,
                          (cx*sz + sx*sy*cz)*x + (cx*cz - sx*sy*sz)*y - sx*cy*z,
                          (sx*sz - cx*sy*cz)*x + (sx*cz + cx*sy*sz)*y + cx*cy*z);
    }

    template<typename T2>
    inline Point3D<T> rotate_inv(const Point3D<T2> &rot) const
    {
        // Rxyz = RxRyRz
        double cx = std::cos(rot.x), sx = std::sin(rot.x);
        double cy = std::cos(rot.y), sy = std::sin(rot.y);
        double cz = std::cos(rot.z), sz = std::sin(rot.z);

        return Point3D<T>(cy*cz*x + (cx*sz + sx*sy*cz)*y + (sx*sz - cx*sy*cz)*z,
                          -cy*sz*x + (cx*cz - sx*sy*sz)*y + (sx*cz + cx*sy*sz)*z,
                          sy*x - sx*cy*y + cx*cy*z);
    }

    inline Point3D<T> rot_x(double a) const
    {
        // Rx(a) = ( 1           0         0  )
        //         ( 0       cos(a)    sin(a) )
        //         ( 0      -sin(a)    cos(a) )
        return Point3D<T>(x, y*std::cos(a) + z*std::sin(a), -y*std::sin(a) + z*std::cos(a));
    }

    inline Point3D<T> rot_y(double a) const
    {
        // Ry(a) = ( cos(a)      0    -sin(a) )
        //         ( 0           1         0  )
        //         ( sin(a)      0     cos(a) )
        return Point3D<T>(x*std::cos(a) - z*std::sin(a), y, x*std::sin(a) + z*std::cos(a));
    }

    inline Point3D<T> rot_z(double a) const
    {
        // Rz(a) = ( cos(a)  sin(a)        0  )
        //         (-sin(a)  cos(a)        0  )
        //         ( 0           0         1  )
        return Point3D<T>(x*std::cos(a) + y*std::sin(a), -x*std::sin(a) + y*std::cos(a), z);
    }

    // find the intersect point of a line and a plane
    // (this_point, p2) forms the line and (p3, normal) forms the plane
    template<typename T2>
    Point3D<T> intersect_plane(const Point3D<T2> &p2, const Point3D<T2> &p3,
                               const Point3D<T2> &normal) const
    {
        T alpha = normal.dot(p3 - *this)/normal.dot(p2 - *this);
        return *this + alpha*(p2 - *this);
    }

    Point3D<T> operator -() const
    {
        return Point3D<T>(-x, -y, -z);
    }

    template<typename T2>
    bool operator ==(const Point3D<T2> &rhs) const
    {
        return (x == rhs.x) && (y == rhs.y) && (z == rhs.z);
    }

    template<typename T2>
    bool operator !=(const Point3D<T2> &rhs) const
    {
        return (x != rhs.x) || (y != rhs.y) || (z != rhs.z);
    }

    template<typename T2>
    Point3D<T> operator ()(const Point3D<T2> &v) const
    {
        return Point3D<T>(v.x, v.y, v.z);
    }

    template<typename T2>
    Point3D<T> operator +(const Point3D<T2> &rhs) const
    {
        return Point3D<T>(x + rhs.x, y + rhs.y, z + rhs.z);
    }
    template<typename T2>
    Point3D<T> &operator +=(const Point3D<T2> &rhs)
    {
        *this = *this + rhs;
        return *this;
    }

    template<typename T2>
    Point3D<T> operator -(const Point3D<T2> &rhs) const
    {
        return Point3D<T>(x - rhs.x, y - rhs.y, z - rhs.z);
    }
    template<typename T2>
    Point3D<T> &operator -=(const Point3D<T2> &rhs)
    {
        *this = *this - rhs;
        return *this;
    }

    template<typename T2>
    Point3D<T> operator *(const T2 &rhs) const
    {
        return Point3D<T>(x*rhs, y*rhs, z*rhs);
    }
    template<typename T2>
    Point3D<T> &operator *=(const T2 &rhs)
    {
        *this = *this * rhs;
        return *this;
    }

    template<typename T2>
    Point3D<T> operator /(const T2 &rhs) const
    {
        return Point3D<T>(x/rhs, y/rhs, z/rhs);
    }
    template<typename T2>
    Point3D<T> &operator /=(const T2 &rhs)
    {
        *this = *this / rhs;
        return *this;
    }
};

template<typename T1, typename T2>
Point3D<T2> operator *(const T1 &lhs, const Point3D<T2> &rhs)
{
    return rhs*lhs;
}

template<typename T>
std::ostream &operator <<(std::ostream &os, const Point3D<T> &p)
{
    os << "(" <<  p.x << ", " << p.y << ", " << p.z << ")";
    return os;
}


// 3D transformation
template<typename T>
struct Transform3D
{
    Point3D<T> trans, rot;

    Transform3D() : trans(0., 0., 0.), rot(0. ,0. ,0.) {}
    Transform3D(T x, T y, T z) : trans(x, y, z), rot(0., 0., 0.) {}
    Transform3D(T x, T y, T z, T rx, T ry, T rz) : trans(x, y, z), rot(rx, ry, rz) {}

    // these functions help to retrieve values in array or set values in array
    T GetCoord(int i)
    {
        if(i == 0) return trans.x;
        if(i == 1) return trans.y;
        if(i == 2) return trans.z;
        if(i == 3) return rot.x;
        if(i == 4) return rot.y;
        if(i == 5) return rot.z;
        return 0.;
    }

    void SetCoord(int i, T val)
    {
        if(i == 0) trans.x = val;
        if(i == 1) trans.y = val;
        if(i == 2) trans.z = val;
        if(i == 3) rot.x = val;
        if(i == 4) rot.y = val;
        if(i == 5) rot.z = val;
    }
};

#endif

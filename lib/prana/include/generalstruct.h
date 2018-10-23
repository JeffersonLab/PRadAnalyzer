#ifndef GENERAL_STRUCT_H
#define GENERAL_STRUCT_H

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

    T dot(const Point2D<T> &rhs) const
    {
        return x*rhs.x + y*rhs.y;
    }

    Point2D<T> operator +(const Point2D<T> &rhs) const
    {
        return Point2D<T>(x + rhs.x, y + rhs.y);
    }
    Point2D<T> operator -(const Point2D<T> &rhs) const
    {
        return Point2D<T>(x - rhs.x, y - rhs.y);
    }
    Point2D<T> operator *(const T &rhs) const
    {
        return Point2D<T>(x*rhs, y*rhs);
    }
    Point2D<T> operator /(const T &rhs) const
    {
        return Point2D<T>(x/rhs, y/rhs);
    }
};

template<typename T>
Point2D<T> operator *(double lhs, const Point2D<T> &rhs)
{
    return Point2D<T>(rhs.x*lhs, rhs.y*lhs);
}

// 3D point
template<typename T>
struct Point3D
{
    T x, y, z;

    Point3D() : x(0.), y(0.), z(0.) {}
    Point3D(T xi, T yi, T zi) : x(xi), y(yi), z(zi) {}

    T dot(const Point3D<T> &rhs) const
    {
        return x*rhs.x + y*rhs.y + z*rhs.z;
    }

    // find the intersect point of a line and a plane
    // (this_point, p2) forms the line and (p3, normal) forms the plane
    Point3D<T> intersect_plane(const Point3D<T> &p2, const Point3D<T> &p3,
                               const Point3D<T> &normal)
    {
        T alpha = normal.dot(p3 - *this)/normal.dot(p2 - *this);
        return *this + alpha*(p2 - *this);
    }

    Point3D<T> operator +(const Point3D<T> &rhs) const
    {
        return Point3D<T>(x + rhs.x, y + rhs.y, z + rhs.z);
    }
    Point3D<T> operator -(const Point3D<T> &rhs) const
    {
        return Point3D<T>(x - rhs.x, y - rhs.y, z - rhs.z);
    }
    Point3D<T> operator *(const T &rhs) const
    {
        return Point3D<T>(x*rhs, y*rhs, z*rhs);
    }
    Point3D<T> operator /(const T &rhs) const
    {
        return Point3D<T>(x/rhs, y/rhs, z/rhs);
    }
};

template<typename T>
Point3D<T> operator *(double lhs, const Point3D<T> &rhs)
{
    return Point3D<T>(rhs.x*lhs, rhs.y*lhs, rhs.z*lhs);
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

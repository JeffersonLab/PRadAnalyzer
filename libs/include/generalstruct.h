#ifndef GENERAL_STRUCT_H
#define GENERAL_STRUCT_H

// general structures that will be used among several classes

// geometry for a HyCal module
struct Geometry
{
    int type;
    double size_x;
    double size_y;
    double size_z;
    double x;
    double y;
    double z;

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
    int sector;
    int row;
    int column;

    Layout() : flag(0), sector(-1), row(0), column(0)
    {};

    Layout(unsigned int f, int s, int r, int c)
    : flag(f), sector(s), row(r), column(c)
    {};
};

// 3D point
struct Point
{
    float x;
    float y;
    float z;

    Point() : x(0), y(0), z(0)
    {};
    Point(float xi, float yi, float zi)
    : x(xi), y(yi), z(zi)
    {};
};

// detector coordinates in a certain frame
struct DetCoord
{
    int run_number; // associated run
    int det_enum;   // detector index
    float x_ori;    // origin x
    float y_ori;    // origin y
    float z_ori;    // origin z
    float theta_x;  // tilting angle on x axis
    float theta_y;  // tilting angle on y axis
    float theta_z;  // tilting angle on z axis

    DetCoord()
    : run_number(0), det_enum(0),
      x_ori(0), y_ori(0), z_ori(0),
      theta_x(0), theta_y(0), theta_z(0)
    {};

    DetCoord(int r, int i, float x, float y, float z)
    : run_number(r), det_enum(i),
      x_ori(x), y_ori(y), z_ori(z),
      theta_x(0), theta_y(0), theta_z(0)
    {};

    DetCoord(int r, int i, float x, float y, float z, float tx, float ty, float tz)
    : run_number(r), det_enum(i),
      x_ori(x), y_ori(y), z_ori(z),
      theta_x(tx), theta_y(ty), theta_z(tz)
    {};

    // these functions help to retrieve values in array or set values in array
    float get_dim_coord(int i)
    {
        if(i == 0) return x_ori;
        if(i == 1) return y_ori;
        if(i == 2) return z_ori;
        if(i == 3) return theta_x;
        if(i == 4) return theta_y;
        if(i == 5) return theta_z;
        return 0.;
    }

    void set_dim_coord(int i, float val)
    {
        if(i == 0) x_ori = val;
        if(i == 1) y_ori = val;
        if(i == 2) z_ori = val;
        if(i == 3) theta_x = val;
        if(i == 4) theta_y = val;
        if(i == 5) theta_z = val;
    }
};

#endif

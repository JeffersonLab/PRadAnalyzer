#ifndef PRAD_DET_COOR_H
#define PRAD_DET_COOR_H

#include <string>
#include <vector>
#include <iostream>
#include <map>
#include <cmath>
#include "canalib.h"
#include "PRadDetector.h"
#include "PRadEventStruct.h"


class PRadHyCalDetector;
class PRadGEMDetector;

class PRadCoordSystem
{
public:
public:
    PRadCoordSystem(const std::string &path = "", const int &run = 0);
    virtual ~PRadCoordSystem();

    // manipulate coordinates database
    void LoadCoordData(const std::string &path, const int &run = 0);
    void SaveCoordData(const std::string &path);
    void ChooseCoord(int run_number);

    // set members
    void SetCurrentCoord(const std::vector<DetCoord> &coords);

    // get members
    const std::map<int ,std::vector<DetCoord>> &GetCoordsData() const {return coords_data;}
    std::vector<DetCoord> GetCurrentCoords() const {return current_coord;}

    // basic transform functions
    void Transform(int det_id, float &x, float &y, float &z) const;

public:
    // template functions
    // transform for clusters with det_id
    template<class T>
    void Transform(T &t)
    const
    {
        Transform(t.det_id, t.x, t.y, t.z);
    }

    // transform for clusters with specified det_id
    template<class T>
    void Transform(int det_id, T &t)
    const
    {
        Transform(det_id, t.x, t.y, t.z);
    }

    template<class T>
    void Transform(int det_id, T *t, int NCluster)
    const
    {
        for(int i = 0; i < NCluster; ++i)
        {
            Transform(det_id, t[i].x, t[i].y, t[i].z);
        }
    }

    // transform for clusters, accepts iterator
    template<class T_it>
    void Transform(int det_id, T_it first, T_it last)
    const
    {
        for(T_it it = first; it != last; ++it)
        {
            Transform(det_id, (*it).x, (*it).y, (*it).z);
        }
    }


    // z-projection
    // by default it projects to HyCal surface from origin
    template<class T>
    void Projection(T &t, const Point &pi = target(),
                    int det_id = (int)PRadDetector::HyCal)
    const
    {
        float zf = current_coord[det_id].z_ori;
        Projection(t.x, t.y, t.z, pi.x, pi.y, pi.z, zf);
    }

    // projection for clusters, accepts array
    template<class T>
    void Projection(T *t, int NCluster, const Point &pi, const float &zf)
    const
    {
        for(int i = 0; i < NCluster; ++i)
        {
            Projection(t[i].x, t[i].y, t[i].z, pi.x, pi.y, pi.z, zf);
        }
    }

    template<class T>
    void Projection(T *t, int NCluster, const Point &pi = target(),
                    int det_id = (int)PRadDetector::HyCal)
    const
    {
        float zf = current_coord[det_id].z_ori;
        for(int i = 0; i < NCluster; ++i)
        {
            Projection(t[i].x, t[i].y, t[i].z, pi.x, pi.y, pi.z, zf);
        }
    }

    // projection for clusters, accepts iterator
    template<class T_it>
    void Projection(T_it first, T_it last, const Point &pi, const float &zf)
    const
    {
        for(T_it it = first; it != last; ++it)
        {
            Projection((*it).x, (*it).y, (*it).z, pi.x, pi.y, pi.z, zf);
        }
    }

    template<class T_it>
    void Projection(T_it first, T_it last, const Point &pi = target(),
                    int det_id = (int)PRadDetector::HyCal)
    const
    {
        float zf = current_coord[det_id].z_ori;
        for(T_it it = first; it != last; ++it)
        {
            Projection((*it).x, (*it).y, (*it).z, pi.x, pi.y, pi.z, zf);
        }
    }


public:
    //static public members
    static Point origin() {return Point(0., 0., 0.);}
    // target cell center is at z = 88.9 mm from survey data
    static Point target() {return Point(0., 0., 88.9);}
    static Point beamline(const float &z) {return Point(0., 0., z);}

    // basic projection functions
    static void Projection(float &x, float &y, float &z,
                           const float &xi, const float &yi, const float &zi,
                           const float &zf);
    static void Projection(Point &p, const Point &pi, const float &zf);
    static void Projection(float &x, float &y, float &z, const float &zf);
    static void Projection(float &x, float &y, float &z, const Point &pi, const float &zf);
    static float ProjectionDistance(Point p1, Point p2, double proj_z = 5725);

    template<class T1, class T2>
    static float ProjectionDistance(const T1 &t1, const T2 &t2, double proj_z = 5725)
    {
        return ProjectionDistance(Point(t1.x, t1.y, t1.z), Point(t2.x, t2.y, t2.z), proj_z);
    }

    template<class T>
    inline float GetPolarAngle(const T &hit, const Point &O = target())
    {
        float z = hit.z - O.z;
        float r = sqrt(pow(z, 2) + pow(hit.x - O.x, 2) + pow(hit.y - O.y, 2));
        return acos(z/r)*cana::rad2deg;
    }

    template<class T>
    inline float GetAzimuthalAngle(const T &hit, const Point &O = target())
    {
        return atan2(hit.y - O.y, hit.x - O.x)*cana::rad2deg;
    }


protected:
    // offsets data, run number as key, order is important, thus use map instead of hash map
    std::map<int, std::vector<DetCoord>> coords_data;
    std::vector<DetCoord> current_coord;
};

std::ostream &operator <<(std::ostream &os, const DetCoord &coord);
#endif

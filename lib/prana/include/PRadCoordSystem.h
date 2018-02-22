#ifndef PRAD_DET_COOR_H
#define PRAD_DET_COOR_H

#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include "canalib.h"
#include "PRadDetector.h"
#include "PRadEventStruct.h"



typedef Point3D<float> Point;
typedef Transform3D<float> DetCoord;

// coordinates information in a run
struct RunCoord
{
    int run_number;
    std::vector<DetCoord> dets;

    RunCoord(int r = 0) : run_number(r)
    {
        dets.resize(static_cast<size_t>(PRadDetector::Max_Dets));
    }

    void Clear()
    {
        run_number = 0;
        dets.clear();
        dets.resize(static_cast<size_t>(PRadDetector::Max_Dets));
    }

    // define operators for searching run number
    bool operator == (const int &run) const {return run_number == run;}
    bool operator != (const int &run) const {return run_number != run;}
    bool operator > (const int &run) const {return run_number > run;}
    bool operator < (const int &run) const {return run_number < run;}
};

std::ostream &operator <<(std::ostream &os, const RunCoord &coord);

class PRadCoordSystem
{
public:
    PRadCoordSystem(const std::string &path = "", const int &run = 0);
    virtual ~PRadCoordSystem();

    // manipulate coordinates database
    void LoadCoordData(const std::string &path, const int &run = 0);
    void SaveCoordData(const std::string &path);
    void ChooseCoord(int run_number, bool warn_not_found = true);
    void ChooseCoordAt(int index);

    // set members
    bool SetCurrentCoord(const RunCoord &coords);

    // get members
    const std::vector<RunCoord> &GetCoordsData() const {return coords_data;}
    RunCoord GetCurrentCoords() const {return current_coord;}

    // basic transform functions
    void Transform(int det_id, float &x, float &y, float &z) const;

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

    // transform the hits on the detector
    template<class DetPtr>
    void TransformHits(DetPtr det)
    const
    {
        for(auto it = det->GetHits().begin(); it != det->GetHits().end(); ++it)
        {
            Transform(det->GetDetID(), it->x, it->y, it->z);
        }
    }

    // z-projection
    // by default it projects to HyCal surface from origin
    template<class T>
    void Projection(T &t, const Point &pi = target(),
                    int det_id = (int)PRadDetector::HyCal)
    const
    {
        float zf = current_coord.dets[det_id].trans.z;
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
        float zf = current_coord.dets[det_id].trans.z;
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
        float zf = current_coord.dets[det_id].trans.z;
        for(T_it it = first; it != last; ++it)
        {
            Projection((*it).x, (*it).y, (*it).z, pi.x, pi.y, pi.z, zf);
        }
    }

public:
    //static public members
    static inline float hycal_z() {return 5731.22;}
    static inline Point origin() {return Point(0., 0., 0.);}
    // target cell center is at z = 88.9 mm from survey data
    static inline Point target() {return Point(0., 0., 88.9);}
    static inline Point beamline(const float &z) {return Point(0., 0., z);}

    // basic projection functions
    static void Projection(float &x, float &y, float &z,
                           const float &xi, const float &yi, const float &zi,
                           const float &zf);
    static void Projection(Point &p, const Point &pi, const float &zf);
    static void Projection(float &x, float &y, float &z, const float &zf);
    static void Projection(float &x, float &y, float &z, const Point &pi, const float &zf);
    static float ProjectionDistance(Point p1, Point p2, Point ori, float proj_z);
    static Point ProjectionCoordDiff(Point p1, Point p2, Point ori, float proj_z);

    template<class T1, class T2>
    static inline float ProjectionDistance(const T1 &t1, const T2 &t2, Point ori = target(), float proj_z = hycal_z())
    {
        return ProjectionDistance(Point(t1.x, t1.y, t1.z), Point(t2.x, t2.y, t2.z), ori, proj_z);
    }

    template<class T1, class T2>
    static inline Point ProjectionCoordDiff(const T1 &t1, const T2 &t2, Point ori = target(), float proj_z = hycal_z())
    {
        return ProjectionCoordDiff(Point(t1.x, t1.y, t1.z), Point(t2.x, t2.y, t2.z), ori, proj_z);
    }

    static inline float GetPolarAngle(const Point &p, const Point &O = target())
    {
        float z = p.z - O.z, x = p.x - O.x, y = p.y - O.y;
        float r = sqrt(z*z + x*x + y*y);
        return acos(z/r)*cana::rad2deg;
    }

    template<class T>
    static inline float GetPolarAngle(const T &hit, const Point &O = target())
    {
        return GetPolarAngle(Point(hit.x, hit.y, hit.z), O);
    }

    static inline float GetAzimuthalAngle(const Point &p, const Point &O = target())
    {
        return atan2(p.y - O.y, p.x - O.x)*cana::rad2deg;
    }

    template<class T>
    static inline float GetAzimuthalAngle(const T &hit, const Point &O = target())
    {
        return atan2(hit.y - O.y, hit.x - O.x)*cana::rad2deg;
    }


protected:
    RunCoord current_coord;
    std::vector<RunCoord> coords_data;
};

#endif

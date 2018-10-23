//============================================================================//
// An example to search the best shift and tilting angles                     //
//                                                                            //
// Chao Peng                                                                  //
// 12/06/2017                                                                 //
//============================================================================//

#include "ConfigParser.h"
#include "generalstruct.h"
#include "canalib.h"
#include "PRadCoordSystem.h"
#include <unordered_map>
#include <cmath>
#include <iostream>

using namespace std;

#define HYCAL_Z (5732.27 - 88.9)
#define GEM1_Z (5314.55 - 88.9)
#define GEM2_Z (5274.91 - 88.9)
#define GEM_Z ((GEM1_Z + GEM2_Z)/2.)
#define PROGRESS_COUNT 1000


Point Transform(Point trans, Point rot, Point coord)
{
    float x = coord.x;
    float y = coord.y;
    float z = coord.z;

    float xt, yt, zt;
    // firstly do the angle tilting
    // basic rotation matrix
    // Rx(a) = ( 1           0         0  )
    //         ( 0       cos(a)    sin(a) )
    //         ( 0      -sin(a)    cos(a) )
    xt = x, yt = y, zt = z;
    y = yt*cos(rot.x*0.001) + zt*sin(rot.x*0.001);
    z = -yt*sin(rot.x*0.001) + zt*cos(rot.x*0.001);

    // Ry(a) = ( cos(a)      0    -sin(a) )
    //         ( 0           1         0  )
    //         ( sin(a)      0     cos(a) )
    xt = x, yt = y, zt = z;
    x = xt*cos(rot.y*0.001) - zt*sin(rot.y*0.001);
    z = xt*sin(rot.y*0.001) + zt*cos(rot.y*0.001);

    // Rz(a) = ( cos(a)  sin(a)        0  )
    //         (-sin(a)  cos(a)        0  )
    //         ( 0           0         1  )
    xt = x, yt = y, zt = z;
    x = xt*cos(rot.z*0.001) + yt*sin(rot.z*0.001);
    y = -xt*sin(rot.z*0.001) + yt*cos(rot.z*0.001);

    // then correct the origin
    x += trans.x;
    y += trans.y;
    z += trans.z;

    return Point(x, y, z);
}

Point RevTransform(Point trans, Point rot, Point coord)
{
    float x = coord.x;
    float y = coord.y;
    float z = coord.z;

    trans = -1.*trans;
    rot = -1.*rot;

    // then correct the origin
    x += trans.x;
    y += trans.y;
    z += trans.z;

    float xt, yt, zt;
    // Rz(a) = ( cos(a)  sin(a)        0  )
    //         (-sin(a)  cos(a)        0  )
    //         ( 0           0         1  )
    xt = x, yt = y, zt = z;
    x = xt*cos(rot.z*0.001) + yt*sin(rot.z*0.001);
    y = -xt*sin(rot.z*0.001) + yt*cos(rot.z*0.001);

    // Ry(a) = ( cos(a)      0    -sin(a) )
    //         ( 0           1         0  )
    //         ( sin(a)      0     cos(a) )
    xt = x, yt = y, zt = z;
    x = xt*cos(rot.y*0.001) - zt*sin(rot.y*0.001);
    z = xt*sin(rot.y*0.001) + zt*cos(rot.y*0.001);

    // firstly do the angle tilting
    // basic rotation matrix
    // Rx(a) = ( 1           0         0  )
    //         ( 0       cos(a)    sin(a) )
    //         ( 0      -sin(a)    cos(a) )
    xt = x, yt = y, zt = z;
    y = yt*cos(rot.x*0.001) + zt*sin(rot.x*0.001);
    z = -yt*sin(rot.x*0.001) + zt*cos(rot.x*0.001);

    return Point(x, y, z);
}

unordered_map<string, Point> TakeDifference(const unordered_map<string, Point> &geometry,
                                            const Point &shift, const Point &tilt)
{
    unordered_map<string, Point> res = geometry;
    // target center in the general coordinate system
    auto target = Point(0., 0., 0.);

    // a point on GEM plane and its normal
    auto gemp = Point(0., 0., GEM_Z);
    auto norm = Point(0., 0., 1.);

    // apply shift and tilt to the plane
    gemp = Transform(shift, tilt, gemp);
    norm = Transform(shift, tilt, norm);

    for (auto &it : res)
    {
        // find geometry
        auto geo = it.second;

        // find intersection on the gem plane
        auto p = geo.intersect_plane(target, gemp, norm);

        // transform back to GEM coordinate
        auto p2 = RevTransform(gemp, tilt, p);

        // assume GEM is well positioned, then transform
        auto p3 = Transform(Point(0., 0., GEM_Z), Point(0., 0., 0.), p2);

        // difference gem - hycal on the hycal plane
        auto diff = PRadCoordSystem::ProjectionCoordDiff(p3, geo, target, HYCAL_Z);
        // cout << it.first << ": " << diff.x << ", " << diff.y << ", " << diff.z << endl;

        it.second = diff;
    }

    return res;
}

unordered_map<string, Point> ReadGeometry(const string &path)
{
    unordered_map<string, Point> res;
    ConfigParser parser;
    parser.OpenFile(path);
    string name, type;
    float sx, sy, sz, x, y, z;
    // hycal coordiate in the general coordinate system
    Point hycal_center = Point(0., 0., HYCAL_Z);
    Point hycal_tilt = Point(0., 0., 0.);
    while(parser.ParseLine())
    {
        parser >> name >> type >> sx >> sy >> sz >> x >> y >> z;
        // transform the geometry to the general coordinate
        res[name] = Transform(hycal_center, hycal_tilt, Point(x, y, z));
    }

    return res;
}

struct PosDiff
{
    float x, y, res_x, res_y;
};

unordered_map<string, PosDiff> ReadPositionDiffs(const string &path)
{
    unordered_map<string, PosDiff> res;
    ConfigParser parser;
    parser.OpenFile(path);
    string name;
    PosDiff pdiff;
    while(parser.ParseLine())
    {
        parser >> name >> pdiff.x >> pdiff.y >> pdiff.res_x >> pdiff.res_y;
        res[name] = pdiff;
    }

    return res;
}

float TiltEstimator(const unordered_map<string, PosDiff> observed_diffs,
                    const unordered_map<string, Point> test_diffs)
{
    float total_sig = 0.;
    for (auto &it : observed_diffs)
    {
        float sig_x = (test_diffs.find(it.first)->second.x - it.second.x)/it.second.res_x;
        float sig_y = (test_diffs.find(it.first)->second.y - it.second.y)/it.second.res_y;
        total_sig += sig_x*sig_x + sig_y*sig_y;
    }

    return std::sqrt(total_sig/(float)observed_diffs.size());
}

struct Ranger
{
    vector<float> min, max, step;
    vector<int> counts;
    int64_t size;
    size_t dim;

    Ranger(vector<float> min, vector<float> max, vector<int> counts)
    : min(min), max(max), counts(counts)
    {
        size = 1;
        dim = counts.size();
        for (size_t i = 0; i < dim; ++i)
        {
            size *= counts[i] + 1;
            if (counts[i] == 0) {
                step.push_back(0.);
            } else {
                step.push_back((max[i] - min[i])/(float)counts[i]);
            }
        }
        size += 1;
    }

    vector<float> At(int64_t count) const
    {
        vector<float> res;
        for(size_t i = 0; i < dim; ++i)
        {
            int curr = count%(counts[i] + 1);
            count /= counts[i] + 1;
            res.push_back(min[i] + step[i]*curr);
        }

        return res;
    }
};

int main(int argc, char *argv[])
{
    // geometry
    auto geometry = ReadGeometry("database/hycal_module.txt");
    auto obdiffs = ReadPositionDiffs("database/transform/delta_pos_table_2GeV.txt");

    // search range
    // shift_x, shift_y, tilt_x, tilt_y, tilt_z
    vector<float> min = {-5, -5, 0, 0, 0.8};
    vector<float> max = {5, 5, 0, 0, 1.2};
    vector<int> counts= {20, 20, 0, 0, 100};

    Ranger range(min, max, counts);

    vector<float> best = {0., 0., 0., 0., 0.};
    auto diffs = TakeDifference(geometry, Point(0., 0., 0.), Point(0., 0., 0.));
    float best_est = TiltEstimator(obdiffs, diffs);
    cout << "Initial estimator = " << best_est << endl;
    unordered_map<string, Point> best_diffs = diffs;

    for(int64_t i = 0; i < range.size; ++i)
    {
        auto vals = range.At(i);
        if(i%PROGRESS_COUNT == 0) {
            cout <<"------[ " << i << "/" << range.size << " ]---"
                 << "---[ Best estimator = "  << best_est << " ]---"
                 << "---[ shift and tilt = (" << best[0] << ", " << best[1]
                 << ", " << best[2] << ", " << best[3] << ", " << best[4]
                 << ") ]------"
                 << "\r" << flush;
        }
        auto shift = Point(vals[0], vals[1], 0.);
        auto tilt = Point(vals[2], vals[3], vals[4]);
        auto diffs = TakeDifference(geometry, shift, tilt);
        auto est = TiltEstimator(obdiffs, diffs);
        if (est < best_est) {
            best = vals;
            best_est = est;
            best_diffs = diffs;
        }
    }

    cout <<"------[ " << range.size << " points tested ]---"
         << "---[ Best estimator = "  << best_est << " ]---"
         << "---[ shift and tilt = (" << best[0] << ", " << best[1]
         << ", " << best[2] << ", " << best[3] << ", " << best[4]
         << ") ]------"
         << endl;

    std::ofstream ofsx("best_diff_x.txt");
    std::ofstream ofsy("best_diff_y.txt");
    for (auto &it : best_diffs)
    {
        ofsx << it.first << "  " << it.second.x << endl;
        ofsy << it.first << "  " << it.second.y << endl;
    }
    ofsx.close();
    ofsy.close();

    return 0;
}


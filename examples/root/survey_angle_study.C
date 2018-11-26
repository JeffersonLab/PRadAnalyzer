#include "PRadCoordSystem.h"
#include "canalib.h"

typedef Point3D<double> Pt;

inline Pt get_rotation(Pt x, Pt y, Pt z)
{
    Pt ux = x/x.norm(), uy = y/y.norm(), uz = z/z.norm();
    return Pt(atan2(-uz.y, uz.z), asin(uz.x), atan2(-uy.x, ux.x));
}

void survey_angle_study()
{
    Pt gem2[4] = {Pt(-80.6308, 104.0049, -390.9068),
                  Pt(-80.0651, 104.0043, -390.9078),
                  Pt(-80.0663, 102.7043, -390.9088),
                  Pt(-80.6321, 102.7046, -390.9079)};
    Pt center(-80.59993, 103.35537, -385.64404);

    cout << "surveyed data:" << endl;
    for(int i = 0; i < 4; ++i)
    {
        cout << 1000.*(gem2[i] - center) << endl;
    }

//    Pt vx = gem2[1] - gem2[0];
//    Pt vy = gem2[0] - gem2[3];
    Pt vx = (gem2[1] + gem2[2])/2. - (gem2[0] + gem2[3])/2.;
    Pt vy = (gem2[0] + gem2[1])/2. - (gem2[2] + gem2[3])/2.;

    // rotation angles
    Pt rot = get_rotation(vx, vy, vx.cross(vy));
    cout << "rotation angles: \n" << rot*cana::rad2deg << endl;

    cout << "backward rotated data:" << endl;
    for(int i = 0; i < 4; ++i)
    {
        //cout << 1000.*(gem2p[i] - hycalp) << endl;
        cout << 1000.*(gem2[i] - center).rotate_inv(rot) << endl;
    }
}

// test if the transform can reconstruct the detector plane
// all the coordinates input are relative to the detector center (0., 0., 0.)
// anchor: the anchor point to do rotation
// system_origin: the origin of the system frame we will need to the reconstruct detector plane
// rotation: rotation angles of the detector plane, unit is degree
void anchor_choice_test(Pt anchor = Pt(10., 10., 0.),
                        Pt system_origin = Pt(-100., -100., -100.),
                        Pt rotation = Pt(5., 10., 15.))
{
    // test four points on the plane, relative to the detector center
    vector<Pt> test_pts{Pt(-10., 10., 0.),
                        Pt(10., 10., 0.),
                        Pt(10., -10., 0.),
                        Pt(-10., -10., 0.)};

    // transform to system frame
    // rotate around detector origin
    // translate according to system origin
    cout << "Real detector positions:" << endl;
    vector<Pt> real_pts;
    for(auto &pt : test_pts)
    {
        real_pts.push_back(pt.rotate(rotation*cana::deg2rad) - system_origin);
        cout << real_pts.back() << endl;
    }
    // we also need to know where is the anchor point
    Pt real_anchor = anchor.rotate(rotation*cana::deg2rad) - system_origin;

    // start with the real detector positions
    // get rotation angles back
    Pt vx = real_pts[1] - real_pts[0];
    Pt vy = real_pts[1] - real_pts[2];
    Pt rot = get_rotation(vx, vy, vx.cross(vy));
    cout << "Rotation angles: " << rot*cana::rad2deg << endl;

    // reconstruct
    // firstly translate to reference frame (anchor as origin)
    // then rotation
    // then translate back to system frame (with real anchor position)
    cout << "Reconstructed detector positions:" << endl;
    for(auto &pt : test_pts)
    {
        pt -= anchor;
        pt = pt.rotate(rot);
        pt += real_anchor;
        cout << pt << endl;
    }

}


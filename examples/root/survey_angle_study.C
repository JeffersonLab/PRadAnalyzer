#include "PRadCoordSystem.h"
#include "canalib.h"

typedef Point3D<double> Pt;

inline Pt get_rotation(Pt x, Pt y, Pt z)
{
    Pt ux = x/x.norm(), uy = y/y.norm(), uz = z/z.norm();
    return Pt(atan2(-uz.y, uz.z), asin(uz.x), atan2(-uy.x, ux.x));
}

inline Pt switch_handedness(Pt rot)
{
    return Pt(-rot.x, -rot.y, rot.z);
}

// this study demonstrates the procedure to reconstruct the detector plane
// in our analysis frame based on the survey results from CEBAF system
void reconstruct_study()
{
    // survey data
    Pt survey[4] = {Pt(-80.6308, 104.0049, -390.9068),
                    Pt(-80.0651, 104.0043, -390.9078),
                    Pt(-80.0663, 102.7043, -390.9088),
                    Pt(-80.6321, 102.7046, -390.9079)};
    Pt scenter(-80.59993, 103.35537, -385.64404);

    // what we would like to see the points in our system
    cout << "surveyed points in analysis frame:" << endl;
    for(auto &p : survey)
    {
        auto ap = 1000.*(p - scenter);
        ap.z = -ap.z;
        cout << ap << endl;
    }

    // now get what we can observe
    // choose an anchor point and use it as the origin of detector frame
    Pt anchor = (survey[0] + survey[1] + survey[2] + survey[3])/4.;

    // get rotation angles
    Pt vx = (survey[1] + survey[2])/2. - (survey[0] + survey[3])/2.;
    Pt vy = (survey[0] + survey[1])/2. - (survey[2] + survey[3])/2.;
    Pt vz = vx.cross(vy);
    Pt rot = get_rotation(vy.cross(vz), vy, vz);

    // rotate back and get the points coordinates in detector frame
    // these will be what we observed
    Pt observe[4];
    cout << "observed points in detector frame:" << endl;
    for(int i = 0; i < 4; ++i)
    {
        observe[i] = 1000.*(survey[i] - anchor).rotate_inv(rot);
        observe[i].z = -observe[i].z;
        cout << observe[i] << endl;
    }
    cout << "observed anchor point in analysis frame:" << endl;
    anchor = 1000.*(anchor - scenter);
    anchor.z = -anchor.z;
    cout << anchor << endl;
    cout << "reconstruction rotation angles:" << endl;
    cout << switch_handedness(rot)*cana::rad2deg << endl;

    // reconstruct from the observed points
    // 1. rotate around anchor
    // 2. translate according to anchor position
    cout << "reconstruct position" << endl;
    for(auto &p : observe)
    {
        cout << p.rotate(switch_handedness(rot)) + anchor << endl;
    }
}

// get gem1 rotation from survey data
void gem1_study()
{
    // unfortunately we only have two surveyed points for gem1
    Pt sgem1[2] = {Pt(-81.1393, 102.7057, -390.9489),
                   Pt(-81.1374, 104.0043, -390.9473)};
    Pt scenter(-80.59993, 103.35537, -385.64404);
    // surveyed rotations, it does not provide roll angle
    Pt srot = Pt(0.0718, -0.1566, 0.)*cana::deg2rad;

    Pt gem1[2];
    cout << "surveyed data:" << endl;
    for(int i = 0; i < 2; ++i)
    {
        cout << 1000.*(sgem1[i] - scenter) << endl;
        gem1[i] = (sgem1[i] - scenter).rotate_inv(srot);
    }

    // we need calculate roll angle, however the y axis are not well defined on
    // GEM plane, using these two points to define y
    Pt vy = gem1[1] - gem1[0];
    Pt vx(1, 0, 0);
    Pt vz = vx.cross(vy);

    // roll angle
    srot.z = get_rotation(vy.cross(vz), vy, vz).z;
    cout << "extracted rotation angles: \n" << srot*cana::rad2deg << endl;

    cout << "backward rotated data:" << endl;
    for(int i = 0; i < 2; ++i)
    {
        cout << 1000.*(sgem1[i] - scenter).rotate_inv(srot) << endl;
    }

    cout << "rotation angles for analysis frame: \n"
         << switch_handedness(srot)*cana::rad2deg << endl;
}

// get gem2 rotation from survey data
void gem2_study()
{
    Pt sgem2[4] = {Pt(-80.6308, 104.0049, -390.9068),
                   Pt(-80.0651, 104.0043, -390.9078),
                   Pt(-80.0663, 102.7043, -390.9088),
                   Pt(-80.6321, 102.7046, -390.9079)};
    Pt scenter(-80.59993, 103.35537, -385.64404);
    // surveyed rotations, it does not provide roll angle
    Pt srot = Pt(0.0453, 0.0939, 0.)*cana::deg2rad;

    Pt gem2[4];
    cout << "surveyed data:" << endl;
    for(int i = 0; i < 4; ++i)
    {
        cout << 1000.*(sgem2[i] - scenter) << endl;
        gem2[i] = (sgem2[i] - scenter).rotate_inv(srot);
    }

    // we need calculate roll angle, however the y axis are not well defined on
    // GEM plane, using the mid points line as y
    Pt vx = (gem2[1] + gem2[2])/2. - (gem2[0] + gem2[3])/2.;
    Pt vy = (gem2[0] + gem2[1])/2. - (gem2[2] + gem2[3])/2.;
    Pt vz = vx.cross(vy);

    // roll angle
    srot.z = get_rotation(vy.cross(vz), vy, vz).z;
    cout << "reconstructed rotation angles: \n" << srot*cana::rad2deg << endl;

    cout << "backward rotated data:" << endl;
    for(int i = 0; i < 4; ++i)
    {
        cout << 1000.*(sgem2[i] - scenter).rotate_inv(srot) << endl;
    }

    cout << "rotation angles for analysis frame: \n"
         << switch_handedness(srot)*cana::rad2deg << endl;
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


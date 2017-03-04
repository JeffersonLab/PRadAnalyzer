//============================================================================//
// Transform the HyCal and GEM to beam center frame                           //
//                                                                            //
// Weizhi Xiong, first version, framework of coordinates system               //
// Xinzhan Bai, studied coordinates offsets, generated database               //
// Chao Peng, changed the coordinates transform method, now all detectors     //
//            are transformed directly to beam frame                          //
// 10/21/2016                                                                 //
//============================================================================//

#include "PRadCoordSystem.h"
#include "PRadHyCalDetector.h"
#include "PRadGEMDetector.h"
#include "ConfigParser.h"
#include <fstream>
#include <iomanip>



//============================================================================//
// Constructor, Destructor                                                    //
//============================================================================//

// constructor
PRadCoordSystem::PRadCoordSystem(const std::string &path, const int &run)
{
    if(!path.empty())
        LoadCoordData(path, run);
}

// destructor
PRadCoordSystem::~PRadCoordSystem()
{
    // place holder
}



//============================================================================//
// Public Member Functions                                                    //
//============================================================================//

// load coordinates data, the format should be
// run_number, det_name, origin_x, origin_y, origin_z, theta_x, theta_y, theta_z
void PRadCoordSystem::LoadCoordData(const std::string &path, const int &chosen_run)
{
    ConfigParser c_parser;

    if(!c_parser.ReadFile(path)) {
        std::cerr << "PRad Coord System Error: Cannot open data file "
                  << "\"" << path << "\"."
                  << std::endl;
        return;
    }

    coords_data.clear();
    current_coord.clear();

    while(c_parser.ParseLine())
    {
        if(c_parser.NbofElements() < 8)
            continue;

        int run;
        std::string det_name;
        float x, y, z, theta_x, theta_y, theta_z;

        c_parser >> run >> det_name
                 >> x >> y >> z >> theta_x >> theta_y >> theta_z;

        size_t index = 0;
        for(; index < (size_t) PRadDetector::Max_Dets; ++index)
        {
            if(det_name.find(PRadDetector::getName(index)) != std::string::npos)
                break;
        }

        if(index >= (size_t) PRadDetector::Max_Dets) {
            std::cout << "PRad Coord System Warning: Unrecognized detector "
                      << det_name << ", skipped reading its offsets."
                      << std::endl;
            continue;
        }

        DetCoord new_off(run, index, x, y, z, theta_x, theta_y, theta_z);

        auto it = coords_data.find(run);
        if(it == coords_data.end()) {
            // create a new entry
            std::vector<DetCoord> new_entry((size_t)PRadDetector::Max_Dets);
            new_entry[index] = new_off;

            coords_data[run] = new_entry;
        } else {
            it->second[index] = new_off;
        }
    }

    ChooseCoord(chosen_run);
}

// save the data
void PRadCoordSystem::SaveCoordData(const std::string &path)
{
    std::ofstream output(path);

    // output headers
    output << "# The least run will be chosen as the default offset" << std::endl
           << "# Lists the origin offsets of detectors to beam center "
           << "and the tilting angles" << std::endl
           << "# Units are in mm and mradian" << std::endl
           << "#" << std::setw(7) << "run"
           << std::setw(10) << "detector"
           << std::setw(12) << "x_origin"
           << std::setw(12) << "y_origin"
           << std::setw(12) << "z_origin"
           << std::setw(12) << "x_tilt"
           << std::setw(12) << "y_tilt"
           << std::setw(12) << "z_tilt"
           << std::endl;

    for(auto &it : coords_data)
    {
        for(auto &coord : it.second)
        {
            output << coord << std::endl;
        }
    }

    output.close();
}

// choose the coordinate offsets from the database
void PRadCoordSystem::ChooseCoord(int run_number)
{
    if(coords_data.empty()) {
        std::cerr << "PRad Coord System Error: Database is empty, make sure you "
                  << "have loaded the correct coordinates data."
                  << std::endl;
        return;
    }

    // choose the default run
    if(run_number == 0) {
        current_coord = coords_data.begin()->second;
        return;
    }

    auto it = coords_data.find(run_number);
    if((it == coords_data.end())) {
        current_coord = coords_data.begin()->second;
        std::cout << "PRad Coord System Warning: Cannot find run " << run_number
                  << " in the current database, choose the default run "
                  << coords_data.begin()->first
                  << std::endl;
    } else {
        current_coord = it->second;
    }
}

// choose coordinates
void PRadCoordSystem::SetCurrentCoord(const std::vector<PRadCoordSystem::DetCoord> &coords)
{
    current_coord = coords;

    current_coord.resize((int)PRadDetector::Max_Dets);

    coords_data[current_coord.begin()->run_number] = current_coord;
}

// Transform the detector frame to beam frame
// it corrects the tilting angle first, and then correct origin
void PRadCoordSystem::Transform(int det_id, float &x, float &y, float &z)
const
{
    const DetCoord &coord = current_coord.at(det_id);

    // firstly do the angle tilting
    // basic rotation matrix
    // Rx(a) = ( 1           0         0  )
    //         ( 0       cos(a)   -sin(a) )
    //         ( 0       sin(a)    cos(a) )
    y = y*cos(coord.theta_x*0.001) + z*sin(coord.theta_x*0.001);
    z = -y*sin(coord.theta_x*0.001) + z*cos(coord.theta_x*0.001);

    // Ry(a) = ( cos(a)      0     sin(a) )
    //         ( 0           1         0  )
    //         (-sin(a)      0     cos(a) )
    x = x*cos(coord.theta_y*0.001) - z*sin(coord.theta_y*0.001);
    z = x*sin(coord.theta_y*0.001) + z*cos(coord.theta_x*0.001);

    // Rz(a) = ( cos(a) -sin(a)        0  )
    //         ( sin(a)  cos(a)        0  )
    //         ( 0           0         1  )
    x = x*cos(coord.theta_z*0.001) + y*sin(coord.theta_z*0.001);
    y = -x*sin(coord.theta_z*0.001) + y*cos(coord.theta_z*0.001);

    // then correct the origin
    x += coord.x_ori;
    y += coord.y_ori;
    z += coord.z_ori;
}

// projection from (xi, yi, zi) to zf
void PRadCoordSystem::Projection(float &x, float &y, float &z,
                                 const float &xi, const float &yi, const float &zi,
                                 const float &zf)
{
    // no need for projection
    if(z == zf)
        return;

    float kx = (xi - x)/(zi - z);
    float ky = (yi - y)/(zi - z);
    x += kx*(zf - z);
    y += ky*(zf - z);
    z = zf;
}

void PRadCoordSystem::Projection(Point &p, const Point &pi, const float &zf)
{
    Projection(p.x, p.y, p.z, pi.x, pi.y, pi.z, zf);
}

// by default project from origin (0, 0, 0)
void PRadCoordSystem::Projection(float &x, float &y, float &z, const float &zf)
{
    Projection(x, y, z, 0, 0, 0, zf);
}

void PRadCoordSystem::Projection(float &x, float &y, float &z, const Point &pi, const float &zf)
{
    Projection(x, y, z, pi.x, pi.y, pi.z, zf);
}

float PRadCoordSystem::ProjectionDistance(PRadCoordSystem::Point p1, PRadCoordSystem::Point p2)
{
    // on the same plane
    if(p1.z == p2.z)
        return sqrt((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y));

    // project to further z
    if(p1.z < p2.z) {
        Projection(p1, target(), p2.z);
    } else {
        Projection(p2, target(), p1.z);
    }

    return sqrt((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y));
}

// origin of the beam center frame
PRadCoordSystem::Point PRadCoordSystem::origin()
{
    return Point(0., 0., 0.);
}

// target center at the beam center frame
PRadCoordSystem::Point PRadCoordSystem::target()
{
    return Point(0., 0., 88.9);
}

// beam line point
PRadCoordSystem::Point PRadCoordSystem::BeamLine(const float &z)
{
    return Point(0., 0., z);
}

std::ostream &operator <<(std::ostream &os, const PRadCoordSystem::DetCoord &det)
{
    return os << std::setw(8)  << det.run_number
              << std::setw(12) << PRadDetector::getName(det.det_enum)
              << std::setw(12) << det.x_ori
              << std::setw(12) << det.y_ori
              << std::setw(12) << det.z_ori
              << std::setw(8) << det.theta_x
              << std::setw(8) << det.theta_y
              << std::setw(8) << det.theta_z;
}

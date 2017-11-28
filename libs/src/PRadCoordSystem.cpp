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
#include "canalib.h"
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
    current_coord.Clear();

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

        DetCoord new_off(x, y, z, theta_x, theta_y, theta_z);

        auto it_pair = cana::binary_search_interval(coords_data.begin(), coords_data.end(), run);
        // not find the exact matched entry
        if(it_pair.second == coords_data.end() || it_pair.first != it_pair.second) {
            // create a new entry
            RunCoord new_entry(run);
            new_entry.dets[index] = new_off;
            // insert the entry before the upper bound to keep the order
            coords_data.insert(it_pair.second, new_entry);
        } else {
            it_pair.first->dets[index] = new_off;
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

    for(auto &coord : coords_data)
    {
        output << coord << std::endl;
    }

    output.close();
}

// choose the coordinate offsets from the database
void PRadCoordSystem::ChooseCoord(int run, bool warn_not_found)
{
    if(coords_data.empty()) {
        std::cerr << "PRad Coord System Error: Database is empty, make sure you "
                  << "have loaded the correct coordinates data."
                  << std::endl;
        return;
    }

    // choose default run
    if(run <= 0) {
        current_coord = coords_data.front();
        return;
    }

    // exceeds the lower bound
    if(run <= coords_data.front().run_number) {
        current_coord = coords_data.front();
    // exceeds the upper bound
    } else if (run >= coords_data.back().run_number) {
        current_coord = coords_data.back();
    // find interval that the run sits in
    } else {
        auto it_pair = cana::binary_search_interval(coords_data.begin(), coords_data.end(), run);
        // always choose the nearest previous run
        current_coord = *it_pair.first;
    }

    // warn not exact
    if(warn_not_found && run != current_coord.run_number) {
        std::cout << "PRad Coord System: Cannot find run <" << run
                  << "> in the current database, choose run <"
                  << current_coord.run_number << "> instead."
                  << std::endl;
    }
}

// choose coordinates by index in the data container
void PRadCoordSystem::ChooseCoordAt(int idx)
{
    if(idx < 0 || idx >= (int)coords_data.size())
        return;

    current_coord = coords_data.at(idx);
}

// set and update the current coordinates
bool PRadCoordSystem::SetCurrentCoord(const RunCoord &coords)
{
    if(coords.dets.size() != static_cast<size_t>(PRadDetector::Max_Dets))
        return false;

    current_coord = coords;

    auto itp = cana::binary_search_interval(coords_data.begin(), coords_data.end(), coords.run_number);
    // add a new entry
    if(itp.second == coords_data.end() || itp.first != itp.second) {
        coords_data.insert(itp.second, coords);
    // exact match
    } else {
        *itp.first = coords;
    }

    return true;
}

// Transform the detector frame to beam frame
// it corrects the tilting angle first, and then correct origin
void PRadCoordSystem::Transform(int det_id, float &x, float &y, float &z)
const
{
    const DetCoord &coord = current_coord.dets.at(det_id);

    // firstly do the angle tilting
    // basic rotation matrix
    // Rx(a) = ( 1           0         0  )
    //         ( 0       cos(a)   -sin(a) )
    //         ( 0       sin(a)    cos(a) )
    y = y*cos(coord.rot.x*0.001) + z*sin(coord.rot.x*0.001);
    z = -y*sin(coord.rot.x*0.001) + z*cos(coord.rot.x*0.001);

    // Ry(a) = ( cos(a)      0     sin(a) )
    //         ( 0           1         0  )
    //         (-sin(a)      0     cos(a) )
    x = x*cos(coord.rot.y*0.001) - z*sin(coord.rot.y*0.001);
    z = x*sin(coord.rot.y*0.001) + z*cos(coord.rot.y*0.001);

    // Rz(a) = ( cos(a) -sin(a)        0  )
    //         ( sin(a)  cos(a)        0  )
    //         ( 0           0         1  )
    x = x*cos(coord.rot.z*0.001) + y*sin(coord.rot.z*0.001);
    y = -x*sin(coord.rot.z*0.001) + y*cos(coord.rot.z*0.001);

    // then correct the origin
    x += coord.trans.x;
    y += coord.trans.y;
    z += coord.trans.z;
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

float PRadCoordSystem::ProjectionDistance(Point p1, Point p2, Point ori, float proj_z)
{
    Projection(p1, ori, proj_z);
    Projection(p2, ori, proj_z);

    return sqrt((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y));
}

Point PRadCoordSystem::ProjectionCoordDiff(Point p1, Point p2, Point ori, float proj_z)
{
    Projection(p1, ori, proj_z);
    Projection(p2, ori, proj_z);

    return Point(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z);
}

std::ostream &operator <<(std::ostream &os, const RunCoord &coord)
{
    for(size_t i = 0; i < coord.dets.size(); ++i)
    {
        const auto &det = coord.dets.at(i);
        os << std::setw(8)  << coord.run_number
           << std::setw(12) << PRadDetector::getName(i)
           << std::setw(12) << det.trans.x
           << std::setw(12) << det.trans.y
           << std::setw(12) << det.trans.z
           << std::setw(8) << det.rot.x
           << std::setw(8) << det.rot.y
           << std::setw(8) << det.rot.z
           << std::endl;
    }

    return os;
}

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
#include <utility>
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

// data structure to store detector info
struct det_setup
{
    std::vector<std::pair<int, int>> run_ranges;
    std::vector<DetCoord> dets;
};

// a helper function to process texts to get detector setup
det_setup process_detector_setup(const std::string &str)
{
    det_setup result;
    result.dets.resize(static_cast<size_t>(PRadDetector::Max_Dets));

    auto text = ConfigParser::break_into_blocks(str, "{", "}");
    ConfigParser c_parser;
    c_parser.SetSplitters("=:");
    std::string var_name, var_value;

    // setup detectors
    for(auto &block : text.blocks)
    {
        if(!ConfigParser::case_ins_equal(block.label, "Detector"))
            continue;
        // find detector block
        c_parser.ReadBuffer(block.content.c_str());
        while(c_parser.ParseLine())
        {
            c_parser >> var_name >> var_value;
            // name to detector id
            int id = PRadDetector::str2DetEnum(var_name.c_str());
            if(id < 0) {
                std::cout << "PRad Coord. System Warning: Unrecognized detector "
                          << var_name << std::endl;
                continue;
            }
            // convert to numbers
            auto vals = ConfigParser::stods(var_value, ",", " \t");
            // set detectors
            for(int i = 0; (i < (int)vals.size()) && (i < 6); ++i)
            {
                result.dets[id].SetCoord(i, vals[i]);
            }
        }
    }

    // setup run ranges
    c_parser.ReadBuffer(text.residual.c_str());
    while(c_parser.ParseLine())
    {
        c_parser >> var_name >> var_value;
        if(!ConfigParser::case_ins_equal(var_name, "Run"))
            continue;

        // convert numbers to set run range
        auto vals = ConfigParser::split(var_value, ",");
        for(auto &val : vals)
        {
            // it could be split by -, i.e., a - b
            auto runs = ConfigParser::split(val, "-");
            if(runs.empty()) {
                std::cout << "PRad Coord. System Warning: Cannot process run range "
                          << val << std::endl;
                continue;
            }
            // begin and end
            int begin = std::stoi(ConfigParser::trim(runs.front(), " \t"));
            int end = std::stoi(ConfigParser::trim(runs.back(), " \t"));
            result.run_ranges.emplace_back(begin, end);
        }
    }

    return result;
}

// a helper function to process texts and add coordinates
void process_beam_position(const std::string &str,
                           const std::vector<det_setup> &setups,
                           std::vector<RunCoord> &container)
{
    ConfigParser c_parser;
    int run;
    double beam_x, beam_y, target_z;
    c_parser.ReadBuffer(str.c_str());

    while(c_parser.ParseLine())
    {
        // read-in run data
        c_parser >> run >> beam_x >> beam_y >> target_z;

        // new run coordinates
        RunCoord coord;
        coord.run_number = run;

        // no setups
        if(setups.empty()) {
            std::cout << "PRad Coord. System Warning: No detector setup found, "
                      << "assume every detector is at the origin."
                      << std::endl;
            coord.dets.resize(static_cast<size_t>(PRadDetector::Max_Dets));
        } else {
            // find the setup
            size_t is = 0;
            for(size_t i = 0; i < setups.size(); ++i)
            {
                for(auto &range : setups[i].run_ranges)
                {
                    if((run >= range.first) && (run <= range.second)) {
                        is = i;
                        break;
                    }
                }
            }

            // convert reference coordinate system to the beam coordinate system
            coord.dets = setups[is].dets;
        }

        // apply beam and target position
        for(auto &det : coord.dets)
        {
            det.trans += Point(beam_x, beam_y, target_z);
        }
        container.emplace_back(coord);
    }
}

// load coordinates data, the format should be
// run_number, det_name, origin_x, origin_y, origin_z, theta_x, theta_y, theta_z
void PRadCoordSystem::LoadCoordData(const std::string &path, const int &chosen_run)
{
    ConfigParser c_parser;
    // read in file
    std::string buffer = ConfigParser::file_to_string(path);

    if(!c_parser.ReadFile(path)) {
        std::cerr << "PRad Coord System Error: Invalid coordinate data file "
                  << "\"" << path << "\"."
                  << std::endl;
        return;
    }

    // remove comments
    ConfigParser::comment_between(buffer, "/*", "*/");
    ConfigParser::comment_line(buffer, "//", "\n");
    ConfigParser::comment_line(buffer, "#", "\n");

    // get content blocks
    auto text = ConfigParser::break_into_blocks(buffer, "{", "}");

    // clear containers
    coords_data.clear();
    current_coord.Clear();

    // retrieve detector setups
    std::vector<det_setup> setups;
    for(auto &block : text.blocks)
    {
        if(ConfigParser::case_ins_equal(block.label, "DetectorSetup")) {
            setups.emplace_back(process_detector_setup(block.content));
        }
    }

    // process beam positions, this should be done after detector setup
    for(auto &block : text.blocks)
    {
        if(ConfigParser::case_ins_equal(block.label, "BeamPosition")) {
            process_beam_position(block.content, setups, coords_data);
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
    // we use mrad as the unit
    auto p =Point3D<float>(x, y, z).transform(coord.trans, coord.rot*0.001);
    x = p.x, y = p.y, z= p.z;
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
           << std::setw(12) << PRadDetector::DetEnum2str(i)
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

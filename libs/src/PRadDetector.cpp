//============================================================================//
// A global list for all the detectors, change name here                      //
//                                                                            //
//============================================================================//

#include "PRadDetector.h"
#include <iostream>
#include <cstring>

static const char *__detector_list[] = {"HyCal", "PRadGEM1", "PRadGEM2", "Undefined"};

const char *PRadDetector::getName(int enumVal)
{
    if(enumVal < 0 || enumVal > (int)Max_Dets) {
        std::cerr << "PRad Detector Error: Does not support Detector ID "
                  << enumVal
                  << ", please check if the Detector ID exists in the PRadDetector list."
                  << std::endl;
        return "";
    }

    return __detector_list[enumVal];
}

int PRadDetector::getID(const char *name)
{
    for(int i = 0; i < (int)Max_Dets; ++i)
        if(strcmp(name, __detector_list[i]) == 0)
            return i;

    std::cerr << "PRad Detectors Error: Cannot find " << name
              << ", please check if the detector name exists in the PRadDetector list."
              << std::endl;
    // not found
    return -1;
}

PRadDetector::PRadDetector(int id)
: det_id(id), det_name(getName(id))
{
    // place holder
}

PRadDetector::PRadDetector(const std::string &n)
: det_id(getID(n.c_str())), det_name(n)
{

}

//============================================================================//
// A global list for all the detectors, change name here                      //
//                                                                            //
//============================================================================//

#include "PRadDetector.h"
#include <iostream>
#include <cstring>

PRadDetector::PRadDetector(int id)
: det_id(id)
{
    det_name = DetEnum2str(static_cast<DetEnum>(id));
}

PRadDetector::PRadDetector(const std::string &n)
: det_name(n)
{
    det_id = str2DetEnum(n.c_str());
}

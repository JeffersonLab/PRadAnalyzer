#ifndef PRAD_DETECTOR_H
#define PRAD_DETECTOR_H

#include <string>
#include "ConfigParser.h"

class PRadDetector
{
public:
    enum DetEnum
    {
        HyCal = 0,
        PRadGEM1,
        PRadGEM2,
        Max_Dets
    };
    // macro in ConfigParser.h
    ENUM_MAP(DetEnum, "HyCal|PRadGEM1|PRadGEM2|Undefined")

public:
    PRadDetector(int id);
    PRadDetector(const std::string &name);

    int GetDetID() const {return det_id;}
    const std::string &GetName() const {return det_name;}

protected:
    int det_id;
    std::string det_name;
};

#endif

#ifndef PRAD_DETECTOR_H
#define PRAD_DETECTOR_H

#include <string>

class PRadDetector
{
public:
    enum DetectorEnum
    {
        HyCal = 0,
        PRadGEM1,
        PRadGEM2,
        Max_Dets
    };

    static const char *getName(int enumVal);
    static int getID(const char *);

public:
    PRadDetector(int id);
    PRadDetector(const std::string &name);

    int GetDetID() const {return det_id;};
    const std::string &GetName() const {return det_name;};

protected:
    int det_id;
    std::string det_name;
};

#endif

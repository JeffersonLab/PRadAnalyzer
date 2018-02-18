#ifndef PRAD_HYCAL_RECONSTRUCTOR_H
#define PRAD_HYCAL_RECONSTRUCTOR_H

#include <vector>
#include <string>
#include "PRadEventStruct.h"
#include "PRadHyCalCluster.h"
#include "PRadSquareCluster.h"
#include "PRadIslandCluster.h"
#include "PRadClusterProfile.h"
#include "ConfigParser.h"



// class to manage the various hycal clustering methods
class PRadHyCalReconstructor
{
public:
    // clustering method enum
    enum MethodEnum
    {
        Undefined = -1,
        Island = 0,
        Square,
        Max_Methods,
    };
    // macro in ConfigParser.h
    ENUM_MAP(MethodEnum, "Island|Square");

public:
    PRadHyCalReconstructor();
    PRadHyCalReconstructor(const std::string &name, const std::string &conf_path);

    // copy/move constructors
    PRadHyCalReconstructor(const PRadHyCalReconstructor &that);
    PRadHyCalReconstructor(PRadHyCalReconstructor &&that);

    // destructor
    virtual ~PRadHyCalReconstructor();

    // copy/move assignment operators
    PRadHyCalReconstructor &operator =(const PRadHyCalReconstructor &rhs);
    PRadHyCalReconstructor &operator =(PRadHyCalReconstructor &&rhs);

    // core functions
    void Reconstruct(class PRadHyCalDetector *det);
    void Reconstruct(class PRadHyCalDetector *det, const EventData &event);

    // profile related
    class PRadClusterProfile *GetProfile() {return profile;}

    // methods information
    bool SetMethod(const std::string &name, const std::string &conf_path);
    PRadHyCalCluster *GetMethod() const {return method;}
    MethodEnum GetMethodType() const {return type;}
    std::string GetMethodName() const {return MethodEnum2str(type);}
    std::vector<std::string> GetMethodNames() const;



private:
    MethodEnum type;
    PRadHyCalCluster *method;
    class PRadClusterProfile *profile;
};

#endif // PRAD_HYCAL_RECONSTRUCTOR

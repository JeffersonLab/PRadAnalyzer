//============================================================================//
// Class to manage HyCal clustering method and reconstruct hits for HyCal     //
// Different clustering methods can be switched                               //
//                                                                            //
// Chao Peng                                                                  //
// 02/18/2018                                                                 //
//============================================================================//

#include "PRadHyCalReconstructor.h"
#include "PRadHyCalDetector.h"
#include "PRadHyCalSystem.h"
#include "PRadClusterProfile.h"
#include <iostream>



// constructors
PRadHyCalReconstructor::PRadHyCalReconstructor()
: type(Undefined), method(nullptr), profile(new PRadClusterProfile())
{
    // place holder
}

PRadHyCalReconstructor::PRadHyCalReconstructor(const std::string &name, const std::string &conf_path)
: profile(new PRadClusterProfile())
{
    SetMethod(name, conf_path);
}

// copy/move constructor
PRadHyCalReconstructor::PRadHyCalReconstructor(const PRadHyCalReconstructor &that)
: type(that.type)
{
    method = that.method->Clone();
    profile = new PRadClusterProfile(*that.profile);
}

PRadHyCalReconstructor::PRadHyCalReconstructor(PRadHyCalReconstructor &&that)
: type(that.type)
{
    method = that.method;
    that.method = nullptr;
    profile = that.profile;
    that.profile = nullptr;
}

// destructor
PRadHyCalReconstructor::~PRadHyCalReconstructor()
{
    delete method, method = nullptr;
    delete profile, profile = nullptr;
}

// copy/move assignment operator
PRadHyCalReconstructor &PRadHyCalReconstructor::operator =(const PRadHyCalReconstructor &rhs)
{
    type = rhs.type;
    method = rhs.method->Clone();
    profile = new PRadClusterProfile(*rhs.profile);
    return *this;
}

PRadHyCalReconstructor &PRadHyCalReconstructor::operator =(PRadHyCalReconstructor &&rhs)
{
    type = rhs.type;
    method = rhs.method;
    rhs.method = nullptr;
    profile = rhs.profile;
    rhs.profile = nullptr;
    return *this;
}



// reconstruct the event to clusters
void PRadHyCalReconstructor::Reconstruct(PRadHyCalDetector *hycal, const EventData &event)
{
    // cannot reconstruct without necessary objects
    if(!hycal || !method) {
        std::cerr << "PRad HyCal Reconstructor Error: undefined method or null "
                  << "detector pointer. Abort event reconstruction."
                  << std::endl;
        return;
    }

    // no need to reconstruct non-physics event
    if(!event.is_physics_event())
        return;

    auto sys = hycal->GetSystem();

    if(!sys) {
        std::cerr << "PRad HyCal Reconstructor Error: cannot reconstruct a given "
                  << "event without a DAQ system connected to the detector. "
                  << "Abort event reconstruction."
                  << std::endl;
        return;
    }

    // collect hits
    method->ClearHits();
    for(auto adc : event.get_adc_data())
    {
        auto channel = sys->GetADCChannel(adc.channel_id);
        if(!channel) continue;

        auto module = channel->GetModule();
        if(!module) continue;

        double val = (double)adc.value - channel->GetPedestal().mean;

        method->AddHit(ModuleHit(module, module->GetID(), module->GetEnergy(val)));
    }

    method->Reconstruct(hycal, profile);

    // add timing information
    for(auto &hit : hycal->GetHits())
    {
        auto center = hycal->GetModule(hit.cid);
        if(!center) continue;

        auto tdc = center->GetTDC();
        if(!tdc) continue;

        auto id = tdc->GetID();
        std::vector<uint16_t> time;
        for(auto &tdc : event.tdc_data)
        {
            if(tdc.channel_id == id)
                time.push_back(tdc.value);
        }
        hit.set_time(time);
    }

}

void PRadHyCalReconstructor::Reconstruct(PRadHyCalDetector *hycal)
{
    // cannot reconstruct without necessary objects
    if(!hycal || !method) {
        std::cerr << "PRad HyCal Reconstructor Error: undefined method or null detector pointer. "
                  << "Abort event reconstruction."
                  << std::endl;
        return;
    }

    // collect hits
    method->CollectHits(hycal);

    // reconstruct
    method->Reconstruct(hycal, profile);

    // add timing information
    for(auto &hit : hycal->GetHits())
    {
        auto center = hycal->GetModule(hit.cid);
        if(!center) continue;

        auto tdc = center->GetTDC();
        if(tdc) hit.set_time(tdc->GetTimeMeasure());
    }
}

bool PRadHyCalReconstructor::SetMethod(const std::string &name, const std::string &conf)
{
    PRadHyCalCluster *newone = nullptr;

    auto newtype = str2MethodEnum(name.c_str());
    if(newtype == type && method && conf.size()) {
        method->Configure(conf);
        return true;
    }

    // create corresponding method
    switch(newtype)
    {
    case Island: newone = new PRadIslandCluster(conf); break;
    case Square: newone = new PRadSquareCluster(conf); break;
    default: break;
    }

    // warn the failure of set method
    if(newone == nullptr) {
        auto method_names = GetMethodNames();
        std::cout << "PRad HyCal System Warning: Cannot find clustering method "
                  << name << ", skip setting method.\n"
                  << "Available methods are: \n";

        for(auto &n : method_names)
        {
            std::cout << "\t" << n << "\n";
        }
        std::cout << std::endl;

        return false;
    }

    // free previous method
    delete method;
    method = newone;
    type = newtype;
    return true;
}

std::vector<std::string> PRadHyCalReconstructor::GetMethodNames()
const
{
    std::vector<std::string> res;
    for(int i = 0; i < static_cast<int>(Max_Methods); ++i)
    {
        res.emplace_back(MethodEnum2str(i));
    }

    return res;
}

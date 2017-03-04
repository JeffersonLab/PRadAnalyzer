//============================================================================//
// GEM FEC class                                                              //
// FEC contains several APVs                                                  //
// The memory of included APV classes will be managed by FEC                  //
//                                                                            //
// Chao Peng                                                                  //
// 10/07/2016                                                                 //
//============================================================================//

#include "PRadGEMFEC.h"
#include "PRadGEMSystem.h"
#include "PRadGEMAPV.h"
#include <iostream>
#include <algorithm>



//============================================================================//
// constructor, assigment operator, destructor                                //
//============================================================================//

// constructor
PRadGEMFEC::PRadGEMFEC(const int &i,
                       const std::string &p,
                       const int &slots,
                       PRadGEMSystem *g)
: gem_srs(g), id(i), ip(p)
{
    // open slots for inserting APVs
    adc_list.resize(slots, nullptr);
}

// copy constructor
PRadGEMFEC::PRadGEMFEC(const PRadGEMFEC &that)
: gem_srs(nullptr), id(that.id), ip(that.ip)
{
    // open slots for inserting APVs
    adc_list.resize(that.adc_list.size(), nullptr);

    for(uint32_t i = 0; i < that.adc_list.size(); ++i)
    {
        // add the copied apv to the same slot
        if(that.adc_list.at(i) != nullptr)
            AddAPV(new PRadGEMAPV(*that.adc_list.at(i)), i);
    }
}

// move constructor
PRadGEMFEC::PRadGEMFEC(PRadGEMFEC &&that)
: gem_srs(nullptr), id(that.id), ip(std::move(that.ip)), adc_list(std::move(that.adc_list))
{
    // reset the connection
    for(uint32_t i = 0; i < adc_list.size(); ++i)
        adc_list[i]->SetFEC(this, i);
}

// desctructor
PRadGEMFEC::~PRadGEMFEC()
{
    UnsetSystem();

    Clear();
}

// copy constructor
PRadGEMFEC &PRadGEMFEC::operator =(const PRadGEMFEC &rhs)
{
    Clear();

    id = rhs.id;
    ip = rhs.ip;
    adc_list.resize(rhs.adc_list.size(), nullptr);

    for(uint32_t i = 0; i < rhs.adc_list.size(); ++i)
    {
        // add the copied apv to the same slot
        if(rhs.adc_list.at(i) != nullptr)
            AddAPV(new PRadGEMAPV(*rhs.adc_list.at(i)), i);
    }
    return *this;
}

// move constructor
PRadGEMFEC &PRadGEMFEC::operator =(PRadGEMFEC &&rhs)
{
    Clear();

    id = rhs.id;
    ip = std::move(rhs.ip);
    adc_list = std::move(rhs.adc_list);

    // reset the connection
    for(uint32_t i = 0; i < adc_list.size(); ++i)
        adc_list[i]->SetFEC(this, i);
    return *this;
}



//============================================================================//
// Public Member Functions                                                    //
//============================================================================//

// set the GEM System for FEC, and disconnect it from the previous GEM System
void PRadGEMFEC::SetSystem(PRadGEMSystem *g, bool force_set)
{
    if(g == gem_srs)
        return;

    if(!force_set)
        UnsetSystem();

    gem_srs = g;
}

void PRadGEMFEC::UnsetSystem(bool force_unset)
{
    if(!gem_srs)
        return;

    if(!force_unset)
        gem_srs->DisconnectFEC(id, true);

    gem_srs = nullptr;
}

// change the capacity
void PRadGEMFEC::SetCapacity(int slots)
{
    // capacity cannot be negative
    if(slots < 0) slots = 0;

    if((uint32_t)slots < adc_list.size())
    {
        std::cout << "PRad GEM FEC Warning: Reduce the slots in FEC "
                  << id << " from " << adc_list.size() << " to " << slots
                  << ". All APVs beyond " << slots << " will be removed. "
                  << std::endl;

        for(uint32_t i = slots; i < adc_list.size(); ++i)
        {
            if(adc_list[i] != nullptr)
                adc_list[i]->UnsetFEC(true);
        }
    }

    adc_list.resize(slots, nullptr);
}

// add an apv to fec, return false if failed
bool PRadGEMFEC::AddAPV(PRadGEMAPV *apv, const int &slot)
{
    if(apv == nullptr)
        return false;

    if((uint32_t)slot >= adc_list.size()) {
        std::cerr << "GEM FEC " << id
                  << ": Abort to add an apv to adc channel "
                  << apv->GetADCChannel()
                  << ", this FEC only has " << adc_list.size()
                  << " channels. (defined in PRadGEMFEC.h)"
                  << std::endl;
        return false;
    }

    if(adc_list.at(slot) != nullptr) {
        std::cerr << "GEM FEC " << id
                  << ": Abort to add an apv to adc channel "
                  << apv->GetADCChannel()
                  << ", channel is occupied"
                  << std::endl;
        return false;
    }

    adc_list[slot] = apv;
    apv->SetFEC(this, slot);
    return true;
}

// remove apv in the slot
void PRadGEMFEC::RemoveAPV(const int &slot)
{
    if((uint32_t)slot >= adc_list.size())
        return;

    auto &apv = adc_list[slot];
    if(apv) {
        apv->UnsetFEC(true);
        delete apv, apv = nullptr;
    }
}

// disconnect apv in the slot
void PRadGEMFEC::DisconnectAPV(const int &slot, bool force_disconn)
{
    if((uint32_t)slot >= adc_list.size())
        return;

    auto &apv = adc_list[slot];
    if(!apv)
        return;

    if(!force_disconn)
        apv->UnsetFEC(true);

    apv = nullptr;
}

// clear all the apvs in fec
void PRadGEMFEC::Clear()
{
    for(auto &adc : adc_list)
    {
        // prevent calling remove apv
        if(adc)
            adc->UnsetFEC(true);
        delete adc;
    }
}

// get the apv in the slot
PRadGEMAPV *PRadGEMFEC::GetAPV(const int &slot)
const
{
    if((uint32_t)slot >= adc_list.size())
        return nullptr;

    return adc_list[slot];
}

std::vector<PRadGEMAPV*> PRadGEMFEC::GetAPVList()
const
{
    std::vector<PRadGEMAPV*> result;
    for(auto &adc : adc_list)
    {
        if(adc != nullptr)
            result.push_back(adc);
    }
    return result;
}


#ifndef PRAD_GEM_FEC_H
#define PRAD_GEM_FEC_H

#include <string>
#include <vector>
#include <unordered_map>
#include "PRadEventStruct.h"
#include "PRadGEMAPV.h"

// maximum channels in a FEC
#define FEC_CAPACITY 9

class PRadGEMSystem;

class PRadGEMFEC
{
public:
    // constructor
    PRadGEMFEC(const int &i,
               const std::string &p,
               const int &slots = FEC_CAPACITY,
               PRadGEMSystem *g = nullptr);

    // copy/move constructors
    PRadGEMFEC(const PRadGEMFEC &that);
    PRadGEMFEC(PRadGEMFEC &&that);

    // descructor
    virtual ~PRadGEMFEC();

    // copy/move assignment operators
    PRadGEMFEC &operator =(const PRadGEMFEC &rhs);
    PRadGEMFEC &operator =(PRadGEMFEC &&rhs);

    // public member functions
    void SetSystem(PRadGEMSystem *g, bool force_set = false);
    void UnsetSystem(bool force_unset = false);
    void SetCapacity(int slots);
    bool AddAPV(PRadGEMAPV *apv, const int &slot);
    void RemoveAPV(const int &slot);
    void DisconnectAPV(const int &slot, bool force_disconn = false);
    void Clear();

    // get parameters
    PRadGEMSystem *GetSystem() const {return gem_srs;};
    int GetID() const {return id;};
    const std::string &GetIP() const {return ip;};
    uint32_t GetCapacity() const {return adc_list.size();};
    PRadGEMAPV *GetAPV(const int &slot) const;
    std::vector<PRadGEMAPV*> GetAPVList() const;

    // functions apply to all apv members
    // functions apply to all apv members
    template<typename... Args>
    void APVControl(void (PRadGEMAPV::*act)(Args...), Args&&... args)
    {
        for(auto apv : adc_list)
        {
            if(apv != nullptr)
                (apv->*act)(std::forward<Args>(args)...);
        }
    }


private:
    PRadGEMSystem *gem_srs;
    int id;
    std::string ip;
    std::vector<PRadGEMAPV*> adc_list;
};

#endif

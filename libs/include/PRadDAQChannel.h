#ifndef PRAD_DAQ_CHANNEL_H
#define PRAD_DAQ_CHANNEL_H

#include <string>
#include <iostream>
#include "datastruct.h"

class PRadDAQChannel
{
public:
    // constructors
    PRadDAQChannel(ChannelAddress addr, bool dead = false);
    PRadDAQChannel(int id, ChannelAddress addr, bool dead = false);
    PRadDAQChannel(std::string name, ChannelAddress addr, bool dead = false);
    PRadDAQChannel(int id, std::string name, ChannelAddress addr, bool dead = false);

    // destructor
    virtual ~PRadDAQChannel();

    // set members
    void SetID(const unsigned short &id) {ch_id = id;};
    void SetName(const std::string &name) {ch_name = name;};
    void SetDead(bool d) {ch_dead = d;};
    void SetAddress(const ChannelAddress &addr) {ch_address = addr;};

    // check if it is a dead channel
    bool IsDead() const {return ch_dead;};
    // get members
    unsigned short GetID() const {return ch_id;};
    const std::string &GetName() const {return ch_name;};
    ChannelAddress GetAddress() const {return ch_address;};

    bool operator <(const PRadDAQChannel &rhs)
    const
    {
        return ch_id < rhs.ch_id;
    }

protected:
    int ch_id;
    std::string ch_name;
    ChannelAddress ch_address;
    bool ch_dead;
};

std::ostream &operator <<(std::ostream &os, const ChannelAddress &addr);
std::ostream &operator <<(std::ostream &os, const PRadDAQChannel &ch);
#endif

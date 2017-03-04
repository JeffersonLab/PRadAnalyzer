//============================================================================//
// Basic DAQ channel unit                                                     //
//                                                                            //
// Chao Peng                                                                  //
// 12/11/2016                                                                 //
//============================================================================//

#include "PRadDAQChannel.h"
#include <iomanip>


//============================================================================//
// Constructors, Destructors, Assignment Operators                            //
//============================================================================//

// constructors
PRadDAQChannel::PRadDAQChannel(ChannelAddress addr, bool d)
: ch_id(-1), ch_name("Undefined"), ch_address(addr), ch_dead(d)
{
    // place holder
}

PRadDAQChannel::PRadDAQChannel(int id, ChannelAddress addr, bool d)
: ch_id(id), ch_name("Undefined"), ch_address(addr), ch_dead(d)
{
    // place holder
}

PRadDAQChannel::PRadDAQChannel(std::string name, ChannelAddress addr, bool d)
: ch_id(-1), ch_name(name), ch_address(addr), ch_dead(d)
{
    // place holder
}

PRadDAQChannel::PRadDAQChannel(int id, std::string name, ChannelAddress addr, bool d)
: ch_id(id), ch_name(name), ch_address(addr), ch_dead(d)
{
    // place holder
}

// destructor
PRadDAQChannel::~PRadDAQChannel()
{
    // place holder
}



//============================================================================//
// Other Functions                                                            //
//============================================================================//

// show channel address
std::ostream &operator <<(std::ostream &os, const ChannelAddress &addr)
{
    return os << std::setw(6) << addr.crate
              << std::setw(6) << addr.slot
              << std::setw(6) << addr.channel;
}

// show channel
std::ostream &operator <<(std::ostream &os, const PRadDAQChannel &ch)
{
    return os << std::setw(8) << ch.GetName() << ch.GetAddress();
}

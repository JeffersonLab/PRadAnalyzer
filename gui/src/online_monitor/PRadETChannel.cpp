//============================================================================//
// A C++ wrapper class for C based ET                                         //
//                                                                            //
// Chao Peng                                                                  //
// 02/27/2016                                                                 //
//============================================================================//

#include "PRadETChannel.h"
#include "PRadETStation.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <sys/types.h>
#include <unistd.h>
#include <time.h>
#include <pthread.h>

using namespace std;

PRadETChannel::PRadETChannel(size_t size)
: curr_stat(nullptr), et_id(nullptr), bufferSize(size)
{
    buffer = new uint32_t[bufferSize];
}

PRadETChannel::~PRadETChannel()
{
    if(buffer != nullptr)
        delete[](buffer), buffer = nullptr;

    // force close ET
    ForceClose();
}

// Close ET connection
void PRadETChannel::ForceClose()
{
    if((et_id != nullptr) && et_alive(et_id)) {
        et_forcedclose(et_id);
        et_id = nullptr;
    }
}

// Open ET
void PRadETChannel::Open(const char* ipAddr, int tcpPort, const char* etFile)
{
    // Use a direct connection to the ET system
    config.SetCast(ET_DIRECT);

    // Set the ip address and tcp port
    config.SetHost(ipAddr);
    config.SetServerPort(tcpPort);

    int charSize = strlen(etFile)+1;
    char *fileName = new char[charSize];
    strncpy(fileName, etFile, charSize);

    // Open et client
    int status = et_open(&et_id, fileName, config.Get());
    delete fileName;

    if(status != ET_OK) {
        throw(PRadException(PRadException::ET_CONNECT_ERROR, "et_client: cannot open et client!"));
    }

    /* set level of debug output */
    et_system_setdebug(et_id, ET_DEBUG_INFO);
}

void PRadETChannel::NewStation(const string &name, int mode)
{
    auto it = stations.find(name);
    if(it == stations.end()) {
        curr_stat = new PRadETStation(this, name, mode);
        stations[string(name)] = curr_stat;
    }
}

void PRadETChannel::SwitchStation(const string &name)
{
    auto it = stations.find(name);
    if(it != stations.end()) {
        curr_stat = it->second;
    } else {
        cout << "ET Channel Warning: station " << name << " does not exist!" << endl;
    }
}

void PRadETChannel::RemoveStation(const string &name)
{
    try {
        if(et_id != nullptr && et_alive(et_id)) {
            auto it = stations.find(name);
            if(it != stations.end()) {
                it->second->Remove();
                stations.erase(it);
            } else {
                cout << "ET Channel Warning: station " << name << " does not exist!" << endl;
            }
        } else {
            cout << "ET Channel Warning: cannot remove station while disconnected from ET!" << endl;
        }
    } catch (PRadException e) {
        throw e;
    }
}

PRadETStation* PRadETChannel::GetStation(const string &name)
{
    auto it = stations.find(name);
    if(it != stations.end()) {
        return it->second;
    } else {
        return nullptr;
    }
}

// Attach station
void PRadETChannel::AttachStation()
{
    try {
        curr_stat->Create();
        curr_stat->Attach();
    } catch (PRadException e) {
        throw(e);
    }
    cout << "Successfully attached to ET!" << endl;
}


// Read one event from ET station, return true if success
bool PRadETChannel::Read()
{
    // check if et is opened or alive
    if(et_id == nullptr || !et_alive(et_id))
        throw(PRadException(PRadException::ET_READ_ERROR,"et_client: et is not opened or dead!"));

    et_att_id att = curr_stat->GetAttachID();

    // get the event
    int status = et_event_get(et_id, att, &etEvent, ET_ASYNC, nullptr);

    switch(status)
    {
    case ET_OK:
        break;
    case ET_ERROR_EMPTY:
        return false;
    case ET_ERROR_DEAD:
        throw(PRadException(PRadException::ET_READ_ERROR,"et_client: et is dead!"));
    case ET_ERROR_TIMEOUT:
        throw(PRadException(PRadException::ET_READ_ERROR,"et_client: got timeout!!"));
    case ET_ERROR_BUSY:
        throw(PRadException(PRadException::ET_READ_ERROR,"et_client: station is busy!"));
    case ET_ERROR_WAKEUP:
        throw(PRadException(PRadException::ET_READ_ERROR,"et_client: someone told me to wake up."));
    default:
        throw(PRadException(PRadException::ET_READ_ERROR,"et_client: unkown error!"));
    }

    // copy the data buffer
    copyEvent();

    // put back the event
    status = et_event_put(et_id, att, etEvent);

    switch(status)
    {
    case ET_OK:
        break;
    case ET_ERROR_DEAD:
        throw(PRadException(PRadException::ET_READ_ERROR,"et_client: et is dead!"));
    default:
        throw(PRadException(PRadException::ET_READ_ERROR,"et_client: unkown error!"));
    }

    return true;
}


bool PRadETChannel::Write(void *buf, int nbytes)
{
    // check if et is opened or alive
    if(et_id == nullptr || !et_alive(et_id))
        throw(PRadException(PRadException::ET_READ_ERROR,"et_client: et is not opened or dead!"));

    et_att_id att = curr_stat->GetAttachID();

    int status = et_event_new(et_id, att, &etEvent, ET_SLEEP, nullptr, nbytes);

    switch(status) {
    case ET_OK:
        break;
    case ET_ERROR_EMPTY:
        return false;
    case ET_ERROR_DEAD:
        throw(PRadException(PRadException::ET_READ_ERROR,"et_client: et is dead!"));
    case ET_ERROR_TIMEOUT:
        throw(PRadException(PRadException::ET_READ_ERROR,"et_client: got timeout!!"));
    case ET_ERROR_BUSY:
        throw(PRadException(PRadException::ET_READ_ERROR,"et_client: station is busy!"));
    case ET_ERROR_WAKEUP:
        throw(PRadException(PRadException::ET_READ_ERROR,"et_client: someone told me to wake up."));
    default:
        throw(PRadException(PRadException::ET_READ_ERROR,"et_client: unkown error!"));
    }

    // build et event
    void *data;
    et_event_getdata(etEvent, &data);
    memcpy((void *) data, (const void *) buf, nbytes);
    et_event_setlength(etEvent, nbytes);

    // put back the event
    status = et_event_put(et_id, att, etEvent);

    switch(status)
    {
    case ET_OK:
        break;
    case ET_ERROR_DEAD:
        throw(PRadException(PRadException::ET_READ_ERROR,"et_client: et is dead!"));
    default:
        throw(PRadException(PRadException::ET_READ_ERROR,"et_client: unkown error!"));
    }

    return true;
}

void PRadETChannel::copyEvent()
{
    void *data;
    et_event_getdata(etEvent, &data);
    et_event_getlength(etEvent, &bufferSize);
    bufferSize /= 4; // from byte to int32 words

    uint32_t *data_buffer = (uint32_t*) data;
    size_t index = 0;
    // check if it is a block header
    if(bufferSize >= 8 && data_buffer[7] == 0xc0da0100) {
        index += 8;
        bufferSize -= 8;
    }

    for(size_t i = 0; i < bufferSize; ++i)
    {
        buffer[i] = data_buffer[index+i];
    }
}


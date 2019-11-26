//============================================================================//
// A C++ wrapper class for C based ET station                                 //
//                                                                            //
// Chao Peng                                                                  //
// 04/01/2016                                                                 //
//============================================================================//

#include "PRadETStation.h"
#include "PRadETChannel.h"


PRadETStation::PRadETStation(PRadETChannel *p, std::string n, int mode)
: et_system(p), name(n)
{
    PreSetting(mode);
}

PRadETStation::~PRadETStation()
{
    Remove();
}

void PRadETStation::PreSetting(int mode)
{
    // Generic settings
    config.SetUser(ET_STATION_USER_MULTI);
    config.SetRestore(ET_STATION_RESTORE_OUT);
    config.SetPrescale(1);
    config.SetCUE(ET_CHUNK_SIZE);

    // TODO, change to meaningful settings
    std::vector<int> selections{17,15,-1,-1};
    char fName[] = "et_my_function";
    char libName[] = "libet_user.so";

    // some pre-defined settings
    // TODO, make these settings selectable
    switch(mode)
    {
    case 1:
        config.SetSelect(ET_STATION_SELECT_ALL);
        config.SetBlock(ET_STATION_BLOCKING);
        break;
    case 2:
        config.SetSelect(ET_STATION_SELECT_ALL);
        config.SetBlock(ET_STATION_NONBLOCKING);
        break;
    case 3:
        config.SetSelect(ET_STATION_SELECT_MATCH);
        config.SetBlock(ET_STATION_BLOCKING);
        config.SetSelectWords(selections);
        break;
    case 4:
        config.SetSelect(ET_STATION_SELECT_MATCH);
        config.SetBlock(ET_STATION_NONBLOCKING);
        config.SetSelectWords(selections);
        break;
    case 5:
        config.SetSelect(ET_STATION_SELECT_USER);
        config.SetBlock(ET_STATION_BLOCKING);
        config.SetSelectWords(selections);
        config.SetFunction(fName);
        config.SetLib(libName);
        break;
    case 6:
        config.SetSelect(ET_STATION_SELECT_USER);
        config.SetBlock(ET_STATION_NONBLOCKING);
        config.SetSelectWords(selections);
        config.SetFunction(fName);
        config.SetLib(libName);
        break;
    }
}

// Create station
void PRadETStation::Create()
{
    char s_name[256];
    strcpy(s_name, name.c_str());

    /* create the station */
    int status = et_station_create(et_system->GetID(), &station_id, s_name, config.Get());

    if(status < ET_OK) {
        if(status == ET_ERROR_EXISTS) {
            /* station_id contains pointer to existing station */;
            throw(PRadException(PRadException::ET_STATION_CREATE_ERROR, "et_client: station already exists!"));
        } else if(status == ET_ERROR_TOOMANY) {
            throw(PRadException(PRadException::ET_STATION_CREATE_ERROR, "et_client: too many stations created!"));
        } else {
            throw(PRadException(PRadException::ET_STATION_CREATE_ERROR, "et_client: error in station creation!"));
        }
    }
}

void PRadETStation::Attach()
{
    if(et_station_attach(et_system->GetID(), station_id, &attach_id) < ET_OK) {
        throw(PRadException(PRadException::ET_STATION_ATTACH_ERROR, "et_client: error in station attach!"));
    }
}

void PRadETStation::Detach()
{
    if(et_station_detach(et_system->GetID(), attach_id) < ET_OK) {
        throw(PRadException(PRadException::ET_STATION_ATTACH_ERROR, "et_client: error in station dettach!"));
    }
}

void PRadETStation::Remove()
{
    if(et_station_remove(et_system->GetID(), station_id) < ET_OK) {
        throw(PRadException(PRadException::ET_STATION_ATTACH_ERROR, "et_client: error in station remove!"));
    }
}


//============================================================================//
// A C++ wrapper class to use CAENHVWrapper library                           //
//                                                                            //
// Chao Peng                                                                  //
// 05/17/2016                                                                 //
//============================================================================//

#include "CAENHVSystem.h"
#include <cstring>

using namespace std;
using namespace CAENHV;

//============================================================================//
// CAEN HV Channel                                                            //
//============================================================================//
CAEN_Channel::~CAEN_Channel()
{
}

void CAEN_Channel::SetPower(const bool& on)
{
    int err;
    unsigned int val = on ? 1 : 0;
    int handle = mother->GetHandle();
    unsigned short slot = mother->GetSlot();

    if(on && mother->GetPrimaryChannelNumber() >= 0) {
        unsigned short list[2];
        list[0] = mother->GetPrimaryChannelNumber();
        list[1] = channel;
        unsigned int vallist[2] = {val, val};
        err = CAENHV_SetChParam(handle, slot, "Pw", 2, list, vallist);
    } else {
        err = CAENHV_SetChParam(handle, slot, "Pw", 1, &channel, &val);
    }

    CAEN_ShowError("HV Channel Set Power", err);
}

void CAEN_Channel::SetVoltage(const float& v)
{
    float val = v;
    if(v > limit) {
        cerr << "HV Channel ERROR: Trying to set voltage " << v
             << " V, which exceeds limit " << limit
             << " V for channel " << name
             << ". Change it to the limit. " << endl;
        val = v;
    }

    int handle = mother->GetHandle();
    unsigned short slot = mother->GetSlot();
    int err = CAENHV_SetChParam(handle, slot, "V0Set", 1, &channel, &val);

    CAEN_ShowError("HV Channel Set Voltage", err);

    if(err == CAENHV_OK)
        Vset = val;
}

void CAEN_Channel::SetName(const string &n)
{
    int handle = mother->GetHandle();
    unsigned short slot = mother->GetSlot();
    int err = CAENHV_SetChName(handle, slot, 1, &channel, n.c_str());

    CAEN_ShowError("HV Channel Set Name", err);

    if(err == CAENHV_OK)
        name = n;
}

void CAEN_Channel::SetLimit(const float &l)
{
    limit = l;
}

void CAEN_Channel::CheckStatus()
{
    int handle = mother->GetHandle();
    unsigned short slot = mother->GetSlot();
    unsigned int status;
    int err = CAENHV_GetChParam(handle, slot, "Status", 1, &channel, &status);
    CAEN_ShowError("HV Channel Read Status", err);
    CAEN_ShowChError(name, status);
}

void CAEN_Channel::ReadVoltage()
{
    int handle = mother->GetHandle();
    unsigned short slot = mother->GetSlot();
    unsigned int pw;
    int err;

    err = CAENHV_GetChParam(handle, slot, "Pw", 1, &channel, &pw);
    CAEN_ShowError("HV Channel Read Power", err);
    on_off = (bool)pw;

    err = CAENHV_GetChParam(handle, slot, "VMon", 1, &channel, &Vmon);
    CAEN_ShowError("HV Channel Read Voltage", err);

    err = CAENHV_GetChParam(handle, slot, "V0Set", 1, &channel, &Vset);
    CAEN_ShowError("HV Channel Read Voltage Set", err);

    if(Vset > limit) {
        cerr << "HV Channel ERROR: Current voltage for channel " << name
             << " is " << Vset
             << " V, which exceeds the limit " << limit << " V"
             << endl;
    }
}

void CAEN_Channel::UpdateVoltage(const bool &pw, const float &vm, const float &vs)
{
    on_off = pw;
    Vmon = vm;
    Vset = vs;

    if(Vset > limit) {
        cerr << "HV Channel ERROR: Current voltage for channel " << name
             << " is " << Vset
             << " V, which exceeds the limit " << limit << " V"
             << endl;
    }
}

//============================================================================//
// CAEN HV Board                                                              //
//============================================================================//
CAEN_Board::~CAEN_Board()
{
    for(auto &channel : channelList)
        delete channel;
}

int CAEN_Board::GetHandle()
{
    return mother->GetHandle();
}

CAEN_Channel *CAEN_Board::GetPrimaryChannel()
{
    return GetChannel(primary);
}

CAEN_Channel *CAEN_Board::GetChannel(int i)
{
    unsigned int index  = (unsigned int) i;
    if(index >= channelList.size())
        return nullptr;

    return channelList[index];
}

void CAEN_Board::ReadBoardMap()
{
    channelList.clear();
    if(nChan < 1)
        return;

    int err;
    float monVals[nChan], setVals[nChan];
    unsigned int pwON[nChan];
    unsigned short list[nChan];
    char nameList[nChan][MAX_CH_NAME];

    for(int k = 0; k < nChan; ++k)
        list[k] = k;

    err = CAENHV_GetChName(mother->GetHandle(), slot, nChan, list, nameList);
    CAEN_ShowError("HV Board Read Name", err);

    err = CAENHV_GetChParam(mother->GetHandle(), slot, "Pw", nChan, list, pwON);
    CAEN_ShowError("HV Board Read Power", err);

    err = CAENHV_GetChParam(mother->GetHandle(), slot, "VMon", nChan, list, monVals);
    CAEN_ShowError("HV Board Read Voltage", err);

    err = CAENHV_GetChParam(mother->GetHandle(), slot, "V0Set", nChan, list, setVals);
    CAEN_ShowError("HV Board Read Voltage Set", err);

    for(int k = 0; k < nChan; ++k)
    {
        string ch_name = nameList[k];
        channelList.push_back(new CAEN_Channel(this, k, ch_name, pwON[k], monVals[k], setVals[k]));
    }

    if(model.find("1932") != string::npos)
        primary = 0;
}

void CAEN_Board::ReadVoltage()
{
    int err;
    int handle = mother->GetHandle();
    float monVals[nChan], setVals[nChan];
    unsigned int pw[nChan];
    unsigned short list[nChan];
    char nameList[nChan][MAX_CH_NAME];

    for(int k = 0; k < nChan; ++k)
        list[k] = k;

    err = CAENHV_GetChName(handle, slot, nChan, list, nameList);
    CAEN_ShowError("HV Board Read Name", err);

    err = CAENHV_GetChParam(handle, slot, "Pw", nChan, list, pw);
    CAEN_ShowError("HV Board Read Power", err);

    err = CAENHV_GetChParam(handle, slot, "VMon", nChan, list, monVals);
    CAEN_ShowError("HV Board Read Voltage", err);

    err = CAENHV_GetChParam(handle, slot, "V0Set", nChan, list, setVals);
    CAEN_ShowError("HV Board Read Voltage Set", err);

    for(int k = 0; k < nChan; ++k)
    {
        if(channelList.at(k)->GetName() != nameList[k]) {
            cerr << "Different channel name, need to update HV system map!"
                 << " Read: " << nameList[k]
                 << " Was: " << channelList.at(k)->GetName()
                 << endl;
            continue;
        }
        channelList.at(k)->UpdateVoltage(pw[k], monVals[k], setVals[k]);
    }
}

void CAEN_Board::CheckStatus()
{
    int err;
    int handle = mother->GetHandle();
    unsigned int status[nChan];
    unsigned short list[nChan];

    for(int k = 0; k < nChan; ++k)
        list[k] = k;

    err = CAENHV_GetChParam(handle, slot, "Status", nChan, list, status);
    CAEN_ShowError("HV Board Read Voltage Set", err);

    for(int k = 0; k < nChan; ++k)
    {
        CAEN_ShowChError(channelList.at(k)->GetName(), status[k]);
    }
}

unsigned short CAEN_Board::GetFirmware()
{
    unsigned short lsb = fmwLSB;
    unsigned short msb = fmwMSB;

    return ((lsb << 8) | msb);
}

void CAEN_Board::SetVoltage(const vector<float> &Vset)
{
    if(Vset.size() != nChan) {
        cout << "HV Board Set Voltage Warning: Mismatched size in value list, "
             << "input size is " << Vset.size()
             << ", should be " << nChan
             << ". Operation terminated."
             << endl;
        return;
    }

    int err;
    int handle = mother->GetHandle();
    unsigned short list[nChan];
    float val[nChan];

    for(int k = 0; k < nChan; ++k)
    {
        list[k] = k;
        val[k] = Vset[k];
    }

    err = CAENHV_SetChParam(handle, slot, "V0Set", nChan, list, val);
    CAEN_ShowError("HV Board Set Voltage", err);
}

void CAEN_Board::SetPower(const bool &on_off)
{
    int err;
    int handle = mother->GetHandle();
    unsigned short list[nChan];
    unsigned int val[nChan];

    for(int k = 0; k < nChan; ++k)
    {
        list[k] = k;
        val[k] = (on_off)? 1: 0;
    }

    err = CAENHV_SetChParam(handle, slot, "Pw", nChan, list, val);
    CAEN_ShowError("HV Board Set Power", err);
}

void CAEN_Board::SetPower(const vector<unsigned int> &on_off)
{
    if(on_off.size() != nChan) {
        cout << "HV Board Set Power Warning: Mismatched size in value list, "
             << "input size is " << on_off.size()
             << ", should be " << nChan
             << ". Operation terminated."
             << endl;
        return;
    }

    int err;
    int handle = mother->GetHandle();
    unsigned short list[nChan];
    unsigned int val[nChan];

    for(int k = 0; k < nChan; ++k)
    {
        list[k] = k;
        val[k] = on_off[k];
    }

    err = CAENHV_SetChParam(handle, slot, "Pw", nChan, list, val);
    CAEN_ShowError("HV Board Set Power", err);
}

//============================================================================//
// CAEN HV Crate                                                              //
//============================================================================//
CAEN_Crate::~CAEN_Crate()
{
    for(auto &board : boardList)
        delete board;

    DeInitialize();
}

bool CAEN_Crate::Initialize()
{
    char arg[32];
    strcpy(arg, ip.c_str());
    int err = CAENHV_InitSystem(sys_type,
                                link_type,
                                arg,
                                username.c_str(),
                                password.c_str(),
                                &handle);
    if(err != CAENHV_OK) {
        CAEN_ShowError("HV Crate Initialize", err);
        return false;
    }

    if(!mapped) {// fist time initialize, does not have crate map
        ReadCrateMap();
    }

    return true;
}

bool CAEN_Crate::DeInitialize()
{
    int err = CAENHV_DeinitSystem(handle);

    if(err != CAENHV_OK) {
        CAEN_ShowError("HV Crate DeInitialize", err);
        return false;
    }

    Clear();

    return true;
}

void CAEN_Crate::Clear()
{
    handle = -1;
    boardList.clear();
    mapped = false;
    for(auto &i : slot_map)
    {
        i = 0;
    }
}

void CAEN_Crate::ReadCrateMap()
{
    if(handle < 0) {
        cerr << "HV Crate Read Map Error: crate "
             << name << " is not initialized!"
             << endl;
        return;
    }

    boardList.clear();
    for(auto &i : slot_map)
    {
        i = 0;
    }

    unsigned short NbofSlot;
    unsigned short *NbofChList;
    char *modelList, *descList;
    unsigned short *serNumList;
    unsigned char *fmwMinList, *fmwMaxList;

    int err = CAENHV_GetCrateMap(handle, &NbofSlot, &NbofChList, &modelList, &descList, &serNumList, &fmwMinList, &fmwMaxList);

    char *m = modelList, *d = descList;
    if(err == CAENHV_OK) {
        for(int slot = 0; slot < NbofSlot; ++slot, m += strlen(m) + 1, d += strlen(d) + 1) {

            if(!NbofChList[slot])
                continue;

            // TODO, get rid of this hard coded exception
            if((id == 5 && slot == 14) ||
               (id == 4 && slot == 12) ||
               (id == 4 && slot == 14)
              )
                continue;

            CAEN_Board *newBoard = new CAEN_Board(this, m, d, slot, NbofChList[slot], serNumList[slot], fmwMinList[slot], fmwMaxList[slot]);
            newBoard->ReadBoardMap();

            slot_map[slot] = boardList.size();
            boardList.push_back(newBoard);

        }
        mapped = true;
    } else {
        CAEN_ShowError("HV Crate Read Map", err);
    }

    free(NbofChList);
    free(modelList);
    free(descList);
    free(serNumList);
    free(fmwMinList);
    free(fmwMaxList);
}

CAEN_Board *CAEN_Crate::GetBoard(const unsigned short &slot)
{
    if(slot >= 50) {
        cerr << "Crate does not have slot "  << slot << endl;
        return nullptr;
    }
    size_t index = slot_map[slot];
    if(index >= boardList.size())
        return nullptr;

    return boardList[index];
}

void CAEN_Crate::HeartBeat()
{
    char sw[30];
    int err = CAENHV_GetSysProp(handle, "SwRelease", sw);
    CAEN_ShowError("HV Crate Heartbeat", err);
}

void CAEN_Crate::CheckStatus()
{
    for(auto &board : boardList)
        board->CheckStatus();
}

void CAEN_Crate::ReadVoltage()
{
    for(auto &board : boardList)
        board->ReadVoltage();
}

void CAEN_Crate::SetPower(const bool &on_off)
{
    for(auto &board : boardList)
        board->SetPower(on_off);
}

//============================================================================//
// CAEN HV General                                                            //
//============================================================================//

ostream &operator<<(ostream &os, CAEN_Board &b)
{
    return os << b.GetModel() << ", " << b.GetDescription() << ", "
              << b.GetSlot() << ", " << b.GetSize() << ", "
              << b.GetSerialNum() << ", "
              << b.GetFirmware() << ".";
}

void CAEN_ShowError(const string &prefix, const int &err, const bool &ShowSuccess)
{
    if(err == CAENHV_OK && !ShowSuccess)
        return;

    string result = prefix + " ERROR: ";

    switch(err)
    {
    case 0: cout << prefix << ": Command is successfully executed," << endl; return;
    case 1: result += "Error of operatived system"; break;
    case 2: result += "Write error in communication channel"; break;
    case 3: result += "Read error in communication channel"; break;
    case 4: result += "Time out in server communication"; break;
    case 5: result += "Command Front End application is down"; break;
    case 6: result += "Communication with system not yet connected by a Login command"; break;
    case 7: result += "Communication with a not present board/slot"; break;
    case 8: result += "Communication with RS232 not yet implemented"; break;
    case 9: result += "User memory not sufficient"; break;
    case 10: result += "Value out of range"; break;
    case 11: result += "Execute command not yet implemented"; break;
    case 12: result += "Get Property not yet implemented"; break;
    case 13: result += "Set Property not yet implemented"; break;
    case 14: result += "Property not found"; break;
    case 15: result += "Execute command not found"; break;
    case 16: result += "No System property"; break;
    case 17: result += "No get property"; break;
    case 18: result += "No set property"; break;
    case 19: result += "No execute command"; break;
    case 20: result += "Device configuration changed"; break;
    case 21: result += "Property of param not found"; break;
    case 22: result += "Param not found"; break;
    case 23: result += "No data present"; break;
    case 24: result += "Device already open"; break;
    case 25: result += "To Many devices opened"; break;
    case 26: result += "Function Parameter not valid"; break;
    case 27: result += "Function not available for the connected device"; break;
    case 0x1001: result += "Device already connected"; break;
    case 0x1002: result += "Device not connected"; break;
    case 0x1003: result += "Operating system error"; break;
    case 0x1004: result += "Login failed"; break;
    case 0x1005: result += "Logout failed"; break;
    case 0x1006: result += "Link type not supported"; break;
    case 0x1007: result += "Incorrect username/password"; break;
    default: result += "Unknown error code"; break;
    }

    cerr << result << endl;
}


float CAEN_VoltageLimit(const string &name)
{
    // Lead Glass
    if(name[0] == 'G') return 1950;
    // Lead Tungstate
    if(name[0] == 'W') return 1450;
    // LMS Reference PMT
    if(name[0] == 'L') return 2000;
    // Scintillator
    if(name[0] == 'S') return 2000;
    // Primary Channel
    if(name[0] == 'P') return 3000;
    // PrimEx Veto Counter Channels
    if(name[0] == 'H') return 2000;

    // not configured
    return 1500;
}

void CAEN_ShowChError(const string &n, const unsigned int &err_bit)
{
    // known problematic channels
    if((n == "W305") || n == "G900") return;

    if(err_bit&(1 << 3)) cerr << "Channel " << n << " is in overcurrent!" << endl;
    if(err_bit&(1 << 4)) cerr << "Channel " << n << " is in overvoltage!" << endl;
    if(err_bit&(1 << 5)) cerr << "Channel " << n << " is in undervoltage!" << endl;
    if(err_bit&(1 << 6)) cerr << "Channel " << n << " is in external trip!" << endl;
    if(err_bit&(1 << 7)) cerr << "Channel " << n << " is in max voltage!" << endl;
    if(err_bit&(1 << 8)) cerr << "Channel " << n << " is in external disable!" << endl;
    if(err_bit&(1 << 9)) cerr << "Channel " << n << " is in internal trip!" << endl;
    if(err_bit&(1 << 10)) cerr << "Channel " << n << " is in calibration error!" << endl;
    if(err_bit&(1 << 11)) cerr << "Channel " << n << " is unplugged!" << endl;
    if(err_bit&(1 << 13)) cerr << "Channel " << n << " is in overvoltage protection!" << endl;
    if(err_bit&(1 << 14)) cerr << "Channel " << n << " is in power fail!" << endl;
    if(err_bit&(1 << 15)) cerr << "Channel " << n << " is in temperature error!" << endl;
}

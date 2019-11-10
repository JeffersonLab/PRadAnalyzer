//============================================================================//
// PRad HyCal System, it contains both HyCal detector and its DAQ system      //
// The connections between HyCal and DAQ system are managed by the system     //
//                                                                            //
// Chao Peng                                                                  //
// 11/12/2016                                                                 //
//============================================================================//

#include "PRadHyCalSystem.h"
#include "PRadHyCalCluster.h"
#include "PRadInfoCenter.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include "canalib.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1D.h"



//============================================================================//
// Constructor, Destructor, Assignment Operators                              //
//============================================================================//

// constructor
PRadHyCalSystem::PRadHyCalSystem(const std::string &path)
: hycal(new PRadHyCalDetector("HyCal", this))
{
    // reserve enough buckets for the adc maps
    adc_addr_map.reserve(ADC_BUCKETS);
    adc_name_map.reserve(ADC_BUCKETS);

    // initialize energy histogram
    energy_hist = new TH1D("HyCal Energy", "Total Energy (MeV)", 2000, 0, 2500);

    if(!path.empty())
        Configure(path);
}

// copy constructor
// it does not only copy the members, but also copy the connections between the
// members
PRadHyCalSystem::PRadHyCalSystem(const PRadHyCalSystem &that)
: ConfigObject(that), hycal(nullptr), recon(that.recon), cal_period(that.cal_period)
{
    // copy detector
    if(that.hycal) {
        hycal = new PRadHyCalDetector(*that.hycal);
        hycal->SetSystem(this, true);
    }

    // copy histogram
    energy_hist = new TH1D(*that.energy_hist);

    // copy tdc
    for(auto tdc : that.tdc_list)
    {
        AddTDCChannel(new PRadTDCChannel(*tdc));
    }

    // copy adc
    for(auto adc : that.adc_list)
    {
        PRadADCChannel *new_adc = new PRadADCChannel(*adc);
        AddADCChannel(new_adc);
        // copy the connections between adc and tdc
        if(!adc->GetTDC())
            continue;
        PRadTDCChannel *tdc = GetTDCChannel(adc->GetTDC()->GetName());
        tdc->ConnectChannel(new_adc);
    }

    // build connections between adc channels and modules
    BuildConnections();
}

// move constructor
PRadHyCalSystem::PRadHyCalSystem(PRadHyCalSystem &&that)
: ConfigObject(that), recon(std::move(that.recon)), cal_period(std::move(that.cal_period)),
  adc_list(std::move(that.adc_list)), tdc_list(std::move(that.tdc_list)),
  adc_addr_map(std::move(that.adc_addr_map)), adc_name_map(std::move(that.adc_name_map)),
  tdc_addr_map(std::move(that.tdc_addr_map)), tdc_name_map(std::move(that.tdc_name_map))
{
    hycal = that.hycal;
    that.hycal = nullptr;
    hycal->SetSystem(this, true);

    energy_hist = that.energy_hist;
    that.energy_hist = nullptr;
}

// destructor
PRadHyCalSystem::~PRadHyCalSystem()
{
    delete hycal;
    delete energy_hist;
    ClearADCChannel();
    ClearTDCChannel();
}

// copy assignment operator
PRadHyCalSystem &PRadHyCalSystem::operator =(const PRadHyCalSystem &rhs)
{
    if(this == &rhs)
        return *this;

    PRadHyCalSystem that(rhs); // copy constructor
    *this = std::move(that); // move assignment operator
    return *this;
}

// move assignment operator
PRadHyCalSystem &PRadHyCalSystem::operator =(PRadHyCalSystem &&rhs)
{
    if(this == &rhs)
        return *this;

    ConfigObject::operator =(rhs);

    // release memories
    delete hycal;
    delete energy_hist;
    ClearADCChannel();
    ClearTDCChannel();

    hycal = rhs.hycal;
    rhs.hycal = nullptr;
    hycal->SetSystem(this, true);
    recon = std::move(rhs.recon);
    energy_hist = rhs.energy_hist;
    rhs.energy_hist = nullptr;
    cal_period = std::move(rhs.cal_period);

    adc_list = std::move(rhs.adc_list);
    tdc_list = std::move(rhs.tdc_list);
    adc_addr_map = std::move(rhs.adc_addr_map);
    adc_name_map = std::move(rhs.adc_name_map);
    tdc_addr_map = std::move(rhs.tdc_addr_map);
    tdc_name_map = std::move(rhs.tdc_name_map);

    return *this;
}



//============================================================================//
// Public Member Functions                                                    //
//============================================================================//

// configure HyCal System
void PRadHyCalSystem::Configure(const std::string &path)
{
    ConfigObject::Configure(path);

    if(hycal) {
        hycal->ReadModuleList(GetConfigValue<std::string>("Module List"));
        hycal->ReadVModuleList(GetConfigValue<std::string>("Virtual Module List"));
    }

    // channel, pedestal and gain factors
    ReadChannelList(GetConfigValue<std::string>("DAQ Channel List"));

    // trigger efficiency
    ReadTriggerEffFile(GetConfigValue<std::string>("Trigger Efficiency Map"));

    // choose clustering method
    recon.SetClusterMethod(GetConfigValue<std::string>("Cluster Method"));
    recon.SetPositionMethod(GetConfigValue<std::string>("Position Method"));

    // configurate reconstructor
    recon.Configure(GetConfigValue<std::string>("Reconstructor Configuration"));

    // lamda to help find strings with "key [type]" in configuration
    auto findstr = [&](const std::string &key, const std::string &type)
                   {
                       return GetConfigValue(key + " [" + type + "]");
                   };

    // load cluster profile
    for(int i = 0; i < static_cast<int>(PRadHyCalModule::Max_Types); ++i)
    {
        auto value = findstr("Cluster Profile", PRadHyCalModule::Type2str(i));
        if(!value.IsEmpty())
            recon.LoadProfile(i, value.String());
    }

    // load density parameters
    for(int i = 0; i < static_cast<int>(PRadClusterDensity::Max_SetEnums); ++i)
    {
        std::string set_name = PRadClusterDensity::SetEnum2str(i);
        auto value1 = findstr("Density Profile", set_name);
        auto value2 = findstr("S-Shape Energy Profile", set_name);
        recon.LoadDensityParams(i, value1.String(), value2.String());
    }

    // read calibration period
    std::string file_path = ConfigParser::form_path(
                            GetConfigValue<std::string>("Calibration Folder"),
                            GetConfigValue<std::string>("Calibration Period File"));
    ReadCalPeriodFile(file_path);

    // set run number
    int run_number;
    CONF_CONN(run_number, "Run Number", 0, false);
    ChooseRun(run_number, false);

    // set resolution for detector
    if(hycal) {
        auto warn = [](int size, const std::string &str)
                    {
                        std::cout << "PRad HyCal System Warning: Expected 3 parameters "
                                  << "for resolution, but received " << size << " in "
                                  << "\"" << str << "\". Skip setting." << std::endl;
                    };

        for(int i = 0; i < static_cast<int>(PRadHyCalDetector::Max_ResRegions); ++i)
        {
            auto type = static_cast<PRadHyCalDetector::ResRegion>(i);
            auto valstr = findstr("Energy Resolution", PRadHyCalDetector::ResRegion2str(i));
            auto vals = ConfigParser::stofs(valstr, ",", " \t");
            if(vals.size() == 3) {
                hycal->SetEneRes(type, vals[0], vals[1], vals[2]);
            } else {
                warn(vals.size(), valstr);
            }

            valstr = findstr("Position Resolution", PRadHyCalDetector::ResRegion2str(i));
            vals = ConfigParser::stofs(valstr, ",", " \t");
            if(vals.size() == 3) {
                hycal->SetPosRes(type, vals[0], vals[1], vals[2]);
            } else {
                warn(vals.size(), valstr);
            }
        }
    }
}

// read DAQ channel list
bool PRadHyCalSystem::ReadChannelList(const std::string &path)
{
    if(path.empty())
        return false;

    ConfigParser c_parser;
    // set special splitter
    c_parser.SetSplitters(",: \t");
    if(!c_parser.ReadFile(path)) {
        std::cerr << "PRad HyCal System Error: Failed to read channel list file "
                  << "\"" << path << "\"."
                  << std::endl;
        return false;
    }

    // we accept 2 types of channels
    // tdc, adc
    std::vector<std::string> types = {"TDC", "ADC"};
    // tdc args: name crate slot channel
    // adc args: name crate slot channel tdc
    std::vector<int> expect_args = {4, 4};
    std::vector<int> option_args = {0, 1};

    // this vector is to store all the following arguments
    std::vector<std::vector<std::vector<ConfigValue>>> ch_args(types.size());

    // read all the elements in
    while(c_parser.ParseLine())
    {
        std::string type = c_parser.TakeFirst();
        size_t i = 0;
        for(; i < types.size(); ++i)
        {
            if(ConfigParser::case_ins_equal(type, types.at(i))) {
                // only save elements from expected format
                if(c_parser.CheckElements(expect_args.at(i), option_args.at(i)))
                    ch_args[i].push_back(c_parser.TakeAll<std::vector>());
                break;
            }
        }

        if(i >= types.size()) { // did not find any type
            std::cout << "PRad HyCal System Warning: Undefined channel type "
                      << type << " in channel list file "
                      << "\"" << path << "\""
                      << std::endl;
        }
    }

    // create TDC first, since it will be needed by ADC channels
    for(auto &args : ch_args[0])
    {
        std::string name(args[0]);
        ChannelAddress addr(args[1].UInt(), args[2].UInt(), args[3].UInt());
        PRadTDCChannel *new_tdc = new PRadTDCChannel(name, addr);
        if(!AddTDCChannel(new_tdc)) // failed to add tdc
            delete new_tdc;
    }

    // create ADC channels, and add them to TDC groups
    for(auto &args : ch_args[1])
    {
        std::string name(args[0]);
        ChannelAddress addr(args[1].UInt(), args[2].UInt(), args[3].UInt());
        PRadADCChannel *new_adc = new PRadADCChannel(name, addr);
        if(!AddADCChannel(new_adc)) { // failed to add adc
            delete new_adc;
            continue;
        }

        // no tdc group specified
        if(args.size() < 5)
            continue;

        // add this adc to tdc group
        std::string tdc_name(args[4]);
        // this adc has no tdc connection
        if(ConfigParser::case_ins_equal(tdc_name, "NONE") ||
           ConfigParser::case_ins_equal(tdc_name, "N/A"))
            continue;

        PRadTDCChannel *tdc = GetTDCChannel(tdc_name);
        if(tdc) {
            tdc->ConnectChannel(new_adc);
        } else {
            std::cout << "PRad HyCal System Warning: ADC Channel " << name
                      << " belongs to TDC Group " << tdc_name
                      << ", but the TDC Channel does not exist in the system."
                      << std::endl;
        }
    }

    // build connection between modules and channels
    BuildConnections();
    return true;
}

// build connections between ADC channels and HyCal modules
void PRadHyCalSystem::BuildConnections()
{
    if(!hycal) {
        std::cout << "PRad HyCal System Warning: HyCal detector does not exist "
                  << "in the system, abort building connections between ADCs "
                  << "and modules"
                  << std::endl;
        return;
    }

    // connect based on name
    for(auto &module : hycal->GetModuleList())
    {
        PRadADCChannel *adc = GetADCChannel(module->GetName());

        if(!adc) { // did not find adc
            std::cout << "PRad HyCal System Warning: Module "
                      << module->GetName()
                      << " has no corresponding ADC channel in the system."
                      << std::endl;
            continue;
        }

        adc->SetModule(module);
        module->SetChannel(adc);
    }
}

// read module status file
bool PRadHyCalSystem::ReadRunInfoFile(const std::string &path)
{
    if(path.empty())
        return false;

    ConfigParser c_parser;
    if(!c_parser.ReadFile(path)) {
        std::cerr << "PRad HyCal System Error: Failed to read status file "
                  << "\"" << path << "\""
                  << std::endl;
        return false;
    }

    std::string name;
    std::vector<double> ref_gain;
    unsigned int ref = 0;

    // first line will be gains for 3 reference PMTs
    if(c_parser.ParseLine()) {
        c_parser >> name;

        if(!ConfigParser::case_ins_equal(name, "REF_GAIN")) {
            std::cerr << "PRad HyCal System Error: Expected Reference PMT info "
                      << "(started by REF_GAIN) as the first input. Aborted status "
                      << "file reading from "
                      << "\"" << path << "\""
                      << std::endl;
            return false;
        }

        // fill in reference PMT gains
        while(c_parser.NbofElements() > 1)
            ref_gain.push_back(c_parser.TakeFirst().Double());

        // get suggested reference number
        c_parser >> ref;
        ref--;

        if(ref >= ref_gain.size()) {
            std::cerr << "PRad HyCal System Error: Unknown Reference PMT "
                      << ref + 1
                      << ", only has " << ref_gain.size()
                      << " Ref. PMTs"
                      << std::endl;
            return false;
        }
    }

    double lms_mean, lms_sig, ped_mean, ped_sig;
    unsigned int status;
    // following lines will be information about modules
    while(c_parser.ParseLine())
    {
        if(!c_parser.CheckElements(6))
            continue;

        c_parser >> name >> ped_mean >> ped_sig >> lms_mean >> lms_sig >> status;

        PRadADCChannel *ch = GetADCChannel(name);
        if(ch) {
            ch->SetPedestal(ped_mean, ped_sig);
            ch->SetDead(status&1);
            PRadHyCalModule *module = ch->GetModule();
            if(module)
                module->GainCorrection((lms_mean - ped_mean)/ref_gain[ref], ref);
        } else {
            std::cout << "PRad HyCal System Warning: Cannot find ADC Channel "
                      << name << ", skip status update and gain correction."
                      << std::endl;
        }
    }

    // finished reading, inform detector to update virtual and dead module neighbors
    if(hycal)
        hycal->UpdateDeadModules();

    return true;
}

// update the trigger efficiency
bool PRadHyCalSystem::ReadTriggerEffFile(const std::string &path)
{
    if(path.empty())
        return false;

    ConfigParser c_parser;
    if(!c_parser.ReadFile(path)) {
        std::cerr << "PRad HyCal System Error: Failed to read trigger efficiency file "
                  << "\"" << path << "\""
                  << std::endl;
        return false;
    }

    std::string name;
    // TRGEFF_NPAR defined in PRadTriggerConst.h
    double pars[TRGEFF_NPAR];

    while(c_parser.ParseLine())
    {
        if(!c_parser.CheckElements(1 + TRGEFF_NPAR))
            continue;

        c_parser >> name;
        for(auto &par : pars)
            par = c_parser.TakeFirst().Double();

        PRadHyCalModule *module = GetModule(name);

        if(module) {
            module->SetTriggerConst(PRadTriggerConst(pars));
        }
    }

    return true;
}

// read file that contains the information about calibration period
bool PRadHyCalSystem::ReadCalPeriodFile(const std::string &path)
{
    if(path.empty())
        return false;

    ConfigParser c_parser;
    if(!c_parser.ReadFile(path)) {
        std::cerr << "PRad HyCal System Error: Failed to read calibration period file "
                  << "\"" << path << "\""
                  << std::endl;
        return false;
    }

    cal_period.clear();

    int period, sub_period, begin, end;

    while(c_parser.ParseLine())
    {
        if(!c_parser.CheckElements(4))
            continue;

        c_parser >> period >> sub_period >> begin >> end;

        // no need to update
        cal_period.emplace_back(begin, end, period, sub_period);
    }

    return true;
}

// set run number from data file path and update related file
void PRadHyCalSystem::ChooseRun(const std::string &path, bool verbose)
{
    if(PRadInfoCenter::SetRunNumber(path)) {
        if(verbose) {
            std::cout << "PRad HyCal System: Set Run Number "
                      << PRadInfoCenter::GetRunNumber()
                      << std::endl;
        }
        UpdateRunFiles(verbose);
    }
}

// set run number and update related file
void PRadHyCalSystem::ChooseRun(int run, bool verbose)
{
    if(PRadInfoCenter::SetRunNumber(run)) {
        if(verbose) {
            std::cout << "PRad HyCal System: Set Run Number "
                      << PRadInfoCenter::GetRunNumber()
                      << std::endl;
        }
        UpdateRunFiles(verbose);
    }
}

// read run number from information center and update related file
void PRadHyCalSystem::UpdateRunFiles(bool verbose)
{
    int run = PRadInfoCenter::GetRunNumber();

    // update config value first since file path will need them
    SetConfigValue("Run Number", run);

    auto it = cana::binary_search(cal_period.begin(), cal_period.end(), run);
    if(it == cal_period.end()) {
        std::cout << "PRad HyCal System Warning: Cannot find calibration period "
                  << "for run " << run << ", assuming period 1-1."
                  << std::endl;
        SetConfigValue("Period", 1);
        SetConfigValue("Sub-period", 1);
    } else {
        SetConfigValue("Period", it->main);
        SetConfigValue("Sub-period", it->sub);
    }

    std::string file_path;
    // calibration file
    file_path = ConfigParser::form_path(GetConfigValue<std::string>("Calibration Folder"),
                                        GetConfigValue<std::string>("Calibration File"));

    if(hycal && hycal->ReadCalibrationFile(file_path) && verbose) {
        std::cout << "PRad HyCal System: Read Calibration File "
                  << "\"" << file_path << "\""
                  << std::endl;
    }

    // run info file
    // calibration file should be read first, since the gain will be corrected
    // based on the read calibration constants
    file_path = ConfigParser::form_path(GetConfigValue<std::string>("Run Info Folder"),
                                        GetConfigValue<std::string>("Run Info File"));

    if(ReadRunInfoFile(file_path) && verbose) {
        std::cout << "PRad HyCal System: Read Run Info File "
                  << "\"" << file_path << "\""
                  << std::endl;
    }

    // choose density profile set for reconstructor
    if(run < 1362) {
        recon.ChooseDensitySet(PRadClusterDensity::Set_1GeV);
    } else {
        recon.ChooseDensitySet(PRadClusterDensity::Set_2GeV);
    }
}

// update the event info to DAQ system
void PRadHyCalSystem::ChooseEvent(const EventData &event)
{
    // clear all the channels
    for(auto &ch : adc_list)
    {
        ch->SetValue(0);
    }

    for(auto &ch : tdc_list)
    {
        ch->ClearTimeMeasure();
    }

    for(auto &adc : event.adc_data)
    {
        if(adc.channel_id >= adc_list.size())
            continue;

        adc_list[adc.channel_id]->SetValue(adc.value);
    }

    for(auto &tdc : event.tdc_data)
    {
        if(tdc.channel_id >= tdc_list.size())
            continue;

        tdc_list[tdc.channel_id]->AddTimeMeasure(tdc.value);
    }
}

// reset current histograms and detector status
void PRadHyCalSystem::Reset()
{
    if(hycal)
        hycal->Reset();

    for(auto &adc : adc_list)
        adc->Reset();
    for(auto &tdc : tdc_list)
        tdc->Reset();
    ResetEnergyHist();
}

// add detector, remove the original detector
void PRadHyCalSystem::SetDetector(PRadHyCalDetector *h)
{
    RemoveDetector();

    hycal = h;

    if(hycal)
        hycal->SetSystem(this);
}

// remove current detector
void PRadHyCalSystem::RemoveDetector()
{
    if(hycal) {
        hycal->UnsetSystem(true);
        delete hycal, hycal = nullptr;
    }
}

void PRadHyCalSystem::DisconnectDetector(bool force_disconn)
{
    if(hycal) {
        if(!force_disconn)
            hycal->UnsetSystem(true);
        hycal = nullptr;
    }
}

// add adc channel
bool PRadHyCalSystem::AddADCChannel(PRadADCChannel *adc)
{
    if(!adc)
        return false;

    if(GetADCChannel(adc->GetAddress())) {
        std::cerr << "PRad HyCal System Error: Failed to add ADC channel "
                  << adc->GetAddress() << ", a channel with the same address exists."
                  << std::endl;
        return false;
    }

    if(GetADCChannel(adc->GetName())) {
        std::cerr << "PRad HyCal System Error: Failed to add ADC channel "
                  << adc->GetName() << ", a channel with the same name exists."
                  << std::endl;
        return false;
    }

    adc->SetID(adc_list.size());
    adc_list.push_back(adc);
    adc_name_map[adc->GetName()] = adc;
    adc_addr_map[adc->GetAddress()] = adc;
    return true;
}

// add tdc channel
bool PRadHyCalSystem::AddTDCChannel(PRadTDCChannel *tdc)
{
    if(!tdc)
        return false;

    if(GetTDCChannel(tdc->GetAddress())) {
        std::cerr << "PRad HyCal System Error: Failed to add TDC channel "
                  << tdc->GetAddress() << ", a channel with the same address exists."
                  << std::endl;
        return false;
    }

    if(GetTDCChannel(tdc->GetName())) {
        std::cerr << "PRad HyCal System Error: Failed to add ADC channel "
                  << tdc->GetName() << ", a channel with the same name exists."
                  << std::endl;
        return false;
    }

    tdc->SetID(tdc_list.size());
    tdc_list.push_back(tdc);
    tdc_name_map[tdc->GetName()] = tdc;
    tdc_addr_map[tdc->GetAddress()] = tdc;
    return true;
}

void PRadHyCalSystem::ClearADCChannel()
{
    for(auto &adc : adc_list)
        delete adc;
    adc_list.clear();
    adc_name_map.clear();
    adc_addr_map.clear();
}

void PRadHyCalSystem::ClearTDCChannel()
{
    for(auto &tdc : tdc_list)
        delete tdc;
    tdc_list.clear();
    tdc_name_map.clear();
    tdc_addr_map.clear();
}

PRadHyCalModule *PRadHyCalSystem::GetModule(const int &id)
const
{
    if(hycal)
        return hycal->GetModule(id);
    return nullptr;
}

PRadHyCalModule *PRadHyCalSystem::GetModule(const std::string &name)
const
{
    if(hycal)
        return hycal->GetModule(name);
    return nullptr;
}

std::vector<PRadHyCalModule*> PRadHyCalSystem::GetModuleList()
const
{
    if(hycal)
        return hycal->GetModuleList();

    return std::vector<PRadHyCalModule*>();
}

PRadADCChannel *PRadHyCalSystem::GetADCChannel(const int &id)
const
{
    if((size_t)id >= adc_list.size())
        return nullptr;
    return adc_list.at(id);
}

PRadADCChannel *PRadHyCalSystem::GetADCChannel(const std::string &name)
const
{
    auto it = adc_name_map.find(name);
    if(it != adc_name_map.end())
        return it->second;
    return nullptr;
}

PRadADCChannel *PRadHyCalSystem::GetADCChannel(const ChannelAddress &addr)
const
{
    auto it = adc_addr_map.find(addr);
    if(it != adc_addr_map.end())
        return it->second;
    return nullptr;
}

PRadTDCChannel *PRadHyCalSystem::GetTDCChannel(const int &id)
const
{
    if((size_t)id >= tdc_list.size())
        return nullptr;
    return tdc_list.at(id);
}

PRadTDCChannel *PRadHyCalSystem::GetTDCChannel(const std::string &name)
const
{
    auto it = tdc_name_map.find(name);
    if(it != tdc_name_map.end())
        return it->second;
    return nullptr;
}

PRadTDCChannel *PRadHyCalSystem::GetTDCChannel(const ChannelAddress &addr)
const
{
    auto it = tdc_addr_map.find(addr);
    if(it != tdc_addr_map.end())
        return it->second;
    return nullptr;
}

void PRadHyCalSystem::Sparsify(const EventData &event)
{
    for(auto &adc : event.adc_data)
    {
       if(adc.channel_id < adc_list.size())
           adc_list[adc.channel_id]->Sparsify();
    }
}

double PRadHyCalSystem::GetEnergy(const EventData &event)
const
{
    double energy = 0.;
    for(auto &adc : event.adc_data)
    {
        if(adc.channel_id >= adc_list.size())
            continue;

        energy += adc_list[adc.channel_id]->GetEnergy(adc.value);
    }

    return energy;
}

// histogram manipulation
void PRadHyCalSystem::FillHists(const EventData &event)
{
    double energy = 0.;

    // adc hists for all types of events
    for(auto &adc : event.get_adc_data())
    {
        if(adc.channel_id >= adc_list.size())
            continue;

        PRadADCChannel *channel = adc_list.at(adc.channel_id);
        channel->FillHist(adc.value, event.get_trigger());
        energy += channel->GetEnergy(adc.value);
    }

    // energy and tdc for only physics events
    if(!event.is_physics_event())
        return;

    energy_hist->Fill(energy);

    for(auto &tdc : event.get_tdc_data())
    {
        if(tdc.channel_id >= tdc_list.size())
            continue;

        PRadTDCChannel *channel = tdc_list.at(tdc.channel_id);
        channel->FillHist(tdc.value);
    }
}

void PRadHyCalSystem::FillEnergyHist()
{
    if(!hycal)
        return;

    energy_hist->Fill(hycal->GetEnergy());
}

void PRadHyCalSystem::FillEnergyHist(const double &e)
{
    energy_hist->Fill(e);
}

void PRadHyCalSystem::FillEnergyHist(const EventData &event)
{
    energy_hist->Fill(GetEnergy(event));
}

void PRadHyCalSystem::ResetEnergyHist()
{
    energy_hist->Reset();
}

void PRadHyCalSystem::SaveHists(const std::string &path)
const
{
    TFile f(path.c_str(), "recreate");

    energy_hist->Write();

    // tdc hists
    TDirectory *cur_dir = f.mkdir("TDC Histograms");
    cur_dir->cd();

    auto ch_id = [] (const PRadDAQChannel &ch)
                 {
                     int id = ch.GetName().at(0)*10000;

                     size_t i = 1;
                     for(; i < ch.GetName().size(); ++i)
                     {
                        if(isdigit(ch.GetName().at(i)))
                            break;
                     }
                     if(i < ch.GetName().size())
                         id += std::stoi(ch.GetName().substr(i));

                     return id;
                 };
    auto ch_order = [ch_id] (const PRadDAQChannel *a, const PRadDAQChannel *b)
                    {return ch_id(*a) < ch_id(*b);};

    // copy the tdc list
    auto tlist = tdc_list;
    std::sort(tlist.begin(), tlist.end(), ch_order);

    for(auto tdc : tlist)
    {
        tdc->GetHist()->Write();
    }

    // adc hists
    f.cd();
    cur_dir = f.mkdir("ADC Histograms");
    cur_dir->cd();

    // copy the adc list
    auto ch_list = adc_list;
    std::sort(ch_list.begin(), ch_list.end(), ch_order);

    TDirectory *mod_dir[PRadHyCalModule::Max_Types];
    for(int i = 0; i < (int) PRadHyCalModule::Max_Types; ++i)
    {
        mod_dir[i] = cur_dir->mkdir(PRadHyCalModule::Type2str(i).c_str());
    }
    TDirectory *other_dir = cur_dir->mkdir("Others");

    for(auto channel : ch_list)
    {
        TDirectory *ch_dir;

        if(!channel->GetModule() || channel->GetModule()->GetType() < 0) {
            ch_dir = other_dir->mkdir(channel->GetName().c_str());
        } else {
            int type = channel->GetModule()->GetType();
            ch_dir = mod_dir[type]->mkdir(channel->GetName().c_str());
        }
        ch_dir->cd();

        std::vector<TH1*> hists = channel->GetHistList();
        for(auto hist : hists)
        {
            hist->Write();
        }
    }

    f.Close();
}

std::vector<double> PRadHyCalSystem::FitHist(const std::string &channel,
                                             const std::string &hist_name,
                                             const std::string &fit_function,
                                             const double &range_min,
                                             const double &range_max,
                                             const bool &verbose)
const
throw(PRadException)
{
    PRadADCChannel *ch = GetADCChannel(channel);
    if(!ch) {
        throw PRadException("Fit Histogram Failure", "Channel " + channel + " does not exist!");
    }

    TH1 *hist = ch->GetHist(hist_name);
    if(hist == nullptr) {
        throw PRadException("Fit Histogram Failure", "Histogram " + hist_name + " does not exist!");
    }

    int beg_bin = hist->GetXaxis()->FindBin(range_min);
    int end_bin = hist->GetXaxis()->FindBin(range_max) - 1;

    if(!hist->Integral(beg_bin, end_bin)) {
        throw PRadException("Fit Histogram Failure", "Histogram does not have entries in specified range!");
    }

    TF1 *fit = new TF1("newfit", fit_function.c_str(), range_min, range_max);

    hist->Fit(fit, "qR");

    TF1 *myfit = (TF1*) hist->GetFunction("newfit");

    // pack parameters, print out result if verbose is true
    std::vector<double> result;

    if(verbose)
        std::cout << "Fit histogram " << hist->GetTitle()
             //<< " with expression " << myfit->GetFormula()->GetExpFormula().Data()
             << std::endl;

    for(int i = 0; i < myfit->GetNpar(); ++i)
    {
        result.push_back(myfit->GetParameter(i));
        if(verbose)
            std::cout << "Parameter " << i << ", "
                      << myfit->GetParameter(i)
                      << std::endl;
    }

    delete fit;

    return result;
}

void PRadHyCalSystem::FitPedestal()
{
    for(auto &channel : adc_list)
    {
        TH1 *ped_hist = channel->GetHist("Pedestal");

        if(ped_hist == nullptr || ped_hist->Integral() < 1000)
            continue;

        ped_hist->Fit("gaus", "qww");

        TF1 *myfit = (TF1*) ped_hist->GetFunction("gaus");
        double p0 = myfit->GetParameter(1);
        double p1 = myfit->GetParameter(2);

        channel->SetPedestal(p0, p1);
    }
}

void PRadHyCalSystem::CorrectGainFactor(int ref)
{
// We had some reference PMT shifts, one happened at run 1228
// If try to do the gain correction that the PMT shifts happened between, it
// won't give the correct calibration constant
#define PED_LED_REF 1000  // separation value for led signal and pedestal signal of reference PMT
#define PED_LED_HYC 30 // separation value for led signal and pedestal signal of all HyCal Modules


    std::string reference = "LMS" + std::to_string(ref);

    // firstly, get the reference factor from LMS PMT
    // LMS 2 seems to be the best one for fitting
    PRadADCChannel *ref_ch = GetADCChannel(reference);
    if(ref_ch == nullptr) {
        std::cerr << "PRad HyCal System Error:"
                  << " Cannot find the reference PMT channel " << reference
                  << " for gain factor correction, abort gain correction."
                  << std::endl;
        return;
    }

    // reference pmt has both pedestal and alpha source signals in this histogram
    TH1* ref_alpha = ref_ch->GetHist("Physics");
    TH1* ref_led = ref_ch->GetHist("LMS");
    if(ref_alpha == nullptr || ref_led == nullptr) {
        std::cerr << "PRad HyCal System Error: "
                  << "Cannot find the histograms of reference PMT, "
                  << "abort gain correction."
                  << std::endl;
        return;
    }

    int sep_bin = ref_alpha->GetXaxis()->FindBin(PED_LED_REF);
    int end_bin = ref_alpha->GetNbinsX(); // 1 for overflow bin and 1 for underflow bin

    if((ref_alpha->Integral(0, sep_bin) < 1000) ||
       (ref_alpha->Integral(sep_bin, end_bin) < 1000)) {
        std::cerr << "PRad HyCal System Error: "
                  << "Not enough entries in pedestal histogram of reference PMT, "
                  << "abort gain correction."
                  << std::endl;
        return;
    }

    // lamda expression, fit a gaussian and return the mean value
    auto fit_gaussian = [] (TH1* hist,
                            const int &range_min = 0,
                            const int &range_max = 8191,
                            const double &warn_ratio = 0.06)
                        {
                            int beg_bin = hist->GetXaxis()->FindBin(range_min);
                            int end_bin = hist->GetXaxis()->FindBin(range_max) - 1;

                            if(hist->Integral(beg_bin, end_bin) < 1000) {
                                std::cout << "PRad HyCal System Warning: "
                                          << "Not enough entries in histogram "
                                          << hist->GetName()
                                          << ". Abort fitting!"
                                          << std::endl;
                                return 0.;
                            }

                            TF1 *fit = new TF1("tmpfit", "gaus", range_min, range_max);

                            hist->Fit(fit, "qR");
                            TF1 *hist_fit = hist->GetFunction("tmpfit");
                            double mean = hist_fit->GetParameter(1);
                            double sigma = hist_fit->GetParameter(2);
                            if(sigma/mean > warn_ratio) {
                                std::cout << "PRad HyCal System Warning: "
                                          << "Bad fit for " << hist->GetTitle()
                                          << ". Mean: " << mean
                                          << ", sigma: " << sigma
                                          << std::endl;
                            }
                            delete fit;
                            return mean;
                        };

    double ped_mean = fit_gaussian(ref_alpha, 0, PED_LED_REF, 0.02);
    double alpha_mean = fit_gaussian(ref_alpha, PED_LED_REF + 1, 8191, 0.05);
    double led_mean = fit_gaussian(ref_led);

    if(ped_mean == 0. || alpha_mean == 0. || led_mean == 0.) {
        std::cerr << "PRad HyCal System Error: Failed to get gain factor from "
                  << reference << ", abort gain correction."
                  << std::endl;
        return;
    }

    double ref_factor = (led_mean - ped_mean)/(alpha_mean - ped_mean);

    for(auto channel : adc_list)
    {
        PRadHyCalModule *module = channel->GetModule();
        if(!module)
            continue;

        TH1 *hist = channel->GetHist("LMS");
        if(!hist)
            continue;

        double ch_led = fit_gaussian(hist) - channel->GetPedestal().mean;

        if(ch_led > PED_LED_HYC) {// meaningful led signal
            module->GainCorrection(ch_led/ref_factor, ref);
        } else {
            std::cout << "PRad HyCal System Error: Gain factor of "
                      << module->GetName()
                      << " is not updated due to bad fit of LED signal."
                      << std::endl;
        }
    }
}


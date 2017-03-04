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
: hycal(new PRadHyCalDetector("HyCal", this)), recon(nullptr)
{
    // reserve enough buckets for the adc maps
    adc_addr_map.reserve(ADC_BUCKETS);
    adc_name_map.reserve(ADC_BUCKETS);

    // initialize energy histogram
    energy_hist = new TH1D("HyCal Energy", "Total Energy (MeV)", 2000, 0, 2500);

    // hycal clustering methods
    AddClusterMethod("Square", new PRadSquareCluster());
    AddClusterMethod("Island", new PRadIslandCluster());
#ifdef USE_PRIMEX_METHOD
    AddClusterMethod("Primex", new PRadPrimexCluster());
#endif
    if(!path.empty())
        Configure(path);
}

// copy constructor
// it does not only copy the members, but also copy the connections between the
// members
PRadHyCalSystem::PRadHyCalSystem(const PRadHyCalSystem &that)
: ConfigObject(that), hycal(nullptr), cal_period(that.cal_period)
{
    // copy detector
    if(that.hycal) {
        hycal = new PRadHyCalDetector(*that.hycal);
        hycal->SetSystem(this, true);
    }

    // copy histogram
    energy_hist = new TH1D(*that.energy_hist);

    // copy reconstruction method
    for(auto &it : that.recon_map)
    {
        recon_map[it.first] = it.second->Clone();
    }
    SetClusterMethod(that.GetClusterMethodName());

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
: ConfigObject(that), cal_period(std::move(that.cal_period)),
  adc_list(std::move(that.adc_list)), tdc_list(std::move(that.tdc_list)),
  adc_addr_map(std::move(that.adc_addr_map)), adc_name_map(std::move(that.adc_name_map)),
  tdc_addr_map(std::move(that.tdc_addr_map)), tdc_name_map(std::move(that.tdc_name_map)),
  recon_map(std::move(that.recon_map))
{
    hycal = that.hycal;
    that.hycal = nullptr;
    hycal->SetSystem(this, true);

    recon = that.recon;
    that.recon = nullptr;

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
    ClearClusterMethods();
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
    ClearClusterMethods();

    hycal = rhs.hycal;
    rhs.hycal = nullptr;
    hycal->SetSystem(this, true);
    recon = rhs.recon;
    rhs.recon = nullptr;
    energy_hist = rhs.energy_hist;
    rhs.energy_hist = nullptr;
    cal_period = std::move(rhs.cal_period);

    adc_list = std::move(rhs.adc_list);
    tdc_list = std::move(rhs.tdc_list);
    adc_addr_map = std::move(rhs.adc_addr_map);
    adc_name_map = std::move(rhs.adc_name_map);
    tdc_addr_map = std::move(rhs.tdc_addr_map);
    tdc_name_map = std::move(rhs.tdc_name_map);
    recon_map = std::move(rhs.recon_map);

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
        hycal->ReadModuleList(GetConfig<std::string>("Module List"));
    }

    // channel, pedestal and gain factors
    ReadChannelList(GetConfig<std::string>("DAQ Channel List"));

    // trigger efficiency
    ReadTriggerEffFile(GetConfig<std::string>("Trigger Efficiency Map"));

    // reconstruction configuration
    SetClusterMethod(GetConfig<std::string>("Cluster Method"));
    if(recon)
        recon->Configure(GetConfig<std::string>("Cluster Configuration"));

    // load profile
    std::string pwo_prof, lg_prof;
    pwo_prof = GetConfig<std::string>("Lead Tungstate Profile");
    PRadClusterProfile::Instance().LoadProfile((int)PRadHyCalModule::PbWO4, pwo_prof);
    lg_prof = GetConfig<std::string>("Lead Glass Profile");
    PRadClusterProfile::Instance().LoadProfile((int)PRadHyCalModule::PbGlass, lg_prof);

#ifdef USE_PRIMEX_METHOD
    // original primex method needs to load the profile into fortran code
    PRadPrimexCluster *method = static_cast<PRadPrimexCluster*>(GetClusterMethod("Primex"));
    if(method) {
        method->LoadCrystalProfile(pwo_prof);
        method->LoadLeadGlassProfile(lg_prof);
    }
#endif

    // read calibration period
    std::string file_path = ConfigParser::form_path(
                            GetConfig<std::string>("Calibration Folder"),
                            GetConfig<std::string>("Calibration Period File"));
    ReadCalPeriodFile(file_path);

    int run_number = getDefConfig<int>("Run Number", 1291, false);
    ChooseRun(run_number, false);
}

// read DAQ channel list
void PRadHyCalSystem::ReadChannelList(const std::string &path)
{
    if(path.empty())
        return;

    ConfigParser c_parser;
    // set special splitter
    c_parser.SetSplitters(",: \t");
    if(!c_parser.ReadFile(path)) {
        std::cerr << "PRad HyCal System Error: Failed to read channel list file "
                  << "\"" << path << "\"."
                  << std::endl;
        return;
    }

    // we accept 2 types of channels
    // tdc, adc
    std::vector<std::string> types = {"TDC", "ADC"};
    // tdc args: name crate slot channel
    // adc args: name crate slot channel tdc
    std::vector<int> expect_args = {4, 4};
    std::vector<int> option_args = {0, 1};

    // this vector is to store all the following arguments
    std::vector<std::vector<ConfigValue>> ch_args[types.size()];

    // read all the elements in
    while(c_parser.ParseLine())
    {
        std::string type = c_parser.TakeFirst();
        size_t i = 0;
        for(; i < types.size(); ++i)
        {
            if(ConfigParser::strcmp_case_insensitive(type, types.at(i))) {
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
        if(ConfigParser::strcmp_case_insensitive(tdc_name, "NONE") ||
           ConfigParser::strcmp_case_insensitive(tdc_name, "N/A"))
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
void PRadHyCalSystem::ReadRunInfoFile(const std::string &path)
{
    if(path.empty())
        return;

    ConfigParser c_parser;
    if(!c_parser.ReadFile(path)) {
        std::cerr << "PRad HyCal System Error: Failed to read status file "
                  << "\"" << path << "\""
                  << std::endl;
        return;
    }

    std::string name;
    std::vector<double> ref_gain;
    unsigned int ref = 0;

    // first line will be gains for 3 reference PMTs
    if(c_parser.ParseLine()) {
        c_parser >> name;

        if(!ConfigParser::strcmp_case_insensitive(name, "REF_GAIN")) {
            std::cerr << "PRad HyCal System Error: Expected Reference PMT info "
                      << "(started by REF_GAIN) as the first input. Aborted status "
                      << "file reading from "
                      << "\"" << path << "\""
                      << std::endl;
            return;
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
            return;
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

    // finished reading, inform detector to update dead module listg
    if(hycal)
        hycal->CreateDeadHits();

#ifdef USE_PRIMEX_METHOD
    // original primex method needs to load the profile into fortran coe
    PRadPrimexCluster *method = static_cast<PRadPrimexCluster*>(GetClusterMethod("Primex"));
    if(method && hycal) {
        method->UpdateModuleStatus(hycal->GetModuleList());
    }
#endif
}

// update the trigger efficiency
void PRadHyCalSystem::ReadTriggerEffFile(const std::string &path)
{
    if(path.empty())
        return;

    ConfigParser c_parser;
    if(!c_parser.ReadFile(path)) {
        std::cerr << "PRad HyCal System Error: Failed to read trigger efficiency file "
                  << "\"" << path << "\""
                  << std::endl;
        return;
    }

    std::string name;
    double eff;
    while(c_parser.ParseLine())
    {
        if(!c_parser.CheckElements(2))
            continue;

        c_parser >> name >> eff;

        PRadHyCalModule *module = GetModule(name);

        if(module) {
            module->SetTriggerEfficiency(eff);
        }
    }
}

// read file that contains the information about calibration period
void PRadHyCalSystem::ReadCalPeriodFile(const std::string &path)
{
    if(path.empty())
        return;

    ConfigParser c_parser;
    if(!c_parser.ReadFile(path)) {
        std::cerr << "PRad HyCal System Error: Failed to read calibration period file "
                  << "\"" << path << "\""
                  << std::endl;
        return;
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
                  << "for run " << run << ", assuming period 5-1."
                  << std::endl;
        SetConfigValue("Period", 5);
        SetConfigValue("Sub-period", 1);
    } else {
        SetConfigValue("Period", it->main);
        SetConfigValue("Sub-period", it->sub);
    }

    // run info file
    std::string file_path;
    file_path = ConfigParser::form_path(GetConfig<std::string>("Run Info Folder"),
                                        GetConfig<std::string>("Run Info File"));

    ReadRunInfoFile(file_path);

    if(verbose) {
        std::cout << "PRad HyCal System: Read Run Info File "
                  << "\"" << file_path << "\""
                  << std::endl;
    }

    // calibration file
    file_path = ConfigParser::form_path(GetConfig<std::string>("Calibration Folder"),
                                        GetConfig<std::string>("Calibration File"));
    if(hycal) {
        hycal->ReadCalibrationFile(file_path);

        if(verbose) {
            std::cout << "PRad HyCal System: Read Calibration File "
                      << "\"" << file_path << "\""
                      << std::endl;
        }
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

// reconstruct the event to clusters
void PRadHyCalSystem::Reconstruct(const EventData &event)
{
    // cannot reconstruct without necessary objects
    if(!hycal || !recon)
        return;

    // no need to reconstruct non-physics event
    if(!event.is_physics_event())
        return;

    // collect hits from eventdata
    auto &hits = hycal->module_hits;

    hits.clear();

    for(auto adc : event.get_adc_data())
    {
        if(adc.channel_id >= adc_list.size())
            continue;

        PRadHyCalModule *module = adc_list.at(adc.channel_id)->GetModule();
        if(!module)
            continue;

        double val = (double)adc.value - adc_list.at(adc.channel_id)->GetPedestal().mean;

        hits.emplace_back(module, module->GetEnergy(val));
    }

    // reoncsturct
    hycal->Reconstruct(recon);
}

void PRadHyCalSystem::Reconstruct()
{
    // cannot reconstruct without necessary objects
    if(!hycal || !recon)
        return;

    // collect current hits
    hycal->CollectHits();

    // reconstruct
    hycal->Reconstruct(recon);
}

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

bool PRadHyCalSystem::AddClusterMethod(const std::string &name, PRadHyCalCluster *c)
{
    // automatically set the first method as the default one
    if(recon_map.empty())
        recon = c;

    std::string key = ConfigParser::str_upper(name);
    auto it = recon_map.find(key);
    // exists, skip
    if(it != recon_map.end()) {
        std::cerr << "PRad HyCal System Error: Clustering method " << name
                  << " exits in the system, abort adding duplicated method."
                  << std::endl;
        return false;
    }

    recon_map[key] = c;
    return true;
}

void PRadHyCalSystem::RemoveClusterMethod(const std::string &name)
{
    std::string key = ConfigParser::str_upper(name);
    auto it = recon_map.find(key);

    if(it != recon_map.end()) {
        if(it->second == recon)
            recon = nullptr;
        delete it->second;
        recon_map.erase(it);
    } else {
        std::cout << "PRad HyCal System Warning: Cannot find clustering method "
                  << name << ", skip removing method."
                  << std::endl;
    }

}

void PRadHyCalSystem::ClearClusterMethods()
{
    for(auto &it : recon_map)
    {
        delete it.second;
    }

    recon_map.clear();
    recon = nullptr;
}

void PRadHyCalSystem::SetClusterMethod(const std::string &name)
{
    if(name.empty())
        return;

    std::string key = ConfigParser::str_upper(name);
    auto it = recon_map.find(key);

    if(it != recon_map.end()) {
        recon = it->second;
    } else {
        std::cout << "PRad HyCal System Warning: Cannot find clustering method "
                  << name << ", skip setting method."
                  << std::endl;
    }
}


PRadHyCalCluster *PRadHyCalSystem::GetClusterMethod(const std::string &name)
const
{
    if(name.empty())
        return recon;

    std::string key = ConfigParser::str_upper(name);
    auto it = recon_map.find(key);

    if(it != recon_map.end())
        return it->second;

    return nullptr;
}

std::string PRadHyCalSystem::GetClusterMethodName()
const
{
    if(recon == nullptr)
        return "";

    for(auto &it : recon_map)
    {
        if(it.second == recon)
            return it.first;
    }

    return "";
}

std::vector<std::string> PRadHyCalSystem::GetClusterMethodNames()
const
{
    std::vector<std::string> result;

    for(auto &it : recon_map)
    {
        if(it.second != nullptr)
            result.push_back(it.first);
    }

    return result;
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

    TDirectory *mod_dir[PRadHyCalModule::Max_Type];
    for(int i = 0; i < (int) PRadHyCalModule::Max_Type; ++i)
    {
        mod_dir[i] = cur_dir->mkdir(PRadHyCalModule::get_module_type_name(i));
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


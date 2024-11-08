//============================================================================//
// An example showing how to replay the dst files and extract all the phyiscs //
// events information and save to root files                                  //
//                                                                            //
// Weizhi Xiong                                                               //
// 08/02/2017                                                                 //
//============================================================================//

#include "PRadDataHandler.h"
#include "PRadDSTParser.h"
#include "PRadInfoCenter.h"
#include "PRadBenchMark.h"
#include "PRadEPICSystem.h"
#include "PRadTaggerSystem.h"
#include "PRadHyCalSystem.h"
#include "PRadGEMSystem.h"
#include "PRadCoordSystem.h"
#include "PRadDetMatch.h"
#include "canalib.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TTree.h"
#include <iostream>
#include <string>
#include <vector>
#include <cassert>
#include <algorithm>
#include "TSystem.h"
#include "TROOT.h"
#include "canalib.h"
#include "TClass.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TMath.h"

using namespace std;
#define MAXC 100
#define PROGRESS_COUNT 10000
#define T_BLOCKS 2156
const Float_t HYCAL_Z = 5646.12;
const float targetMass[2] = {938.272046, 0.510998928};

struct HyCalModule{
    int    cid;
    float  x;
    float  y;
    float  sizex;
    float  sizey;
    HyCalModule() : cid(0), x(0), y(0), sizex(0), sizey(0) {}
    HyCalModule(int s, float x1, float y1, float size1, float size2) :
    cid(s), x(x1), y(y1), sizex(size1), sizey(size2) {}
    void SetInfo(int s, float x1, float y1, float size1, float size2) {
        cid    = s;
        x      = x1;
        y      = y1;
        sizex  = size1;
        sizey  = size2;
    }
    bool Contain(double x1, double y1){
        if (fabs(x1-x) <= sizex/2. && fabs(y1-y) <= sizey/2.) return true;
        else return false;
    }
};

string inputDir = "../../../../dst_files/sel_pro";
bool DoSShapeCorr = false;

PRadEPICSystem *epics;
PRadHyCalSystem *hycal;
PRadGEMSystem *gem;
PRadCoordSystem *coord_sys;
PRadDetMatch *det_match;

void physReplay(const string &path);
void InitTree(TTree* t);
ostream &operator <<(ostream &os, const PRadBenchMark &timer);

Int_t eventNumber;
Float_t Ebeam;
//for HyCal
Int_t clusterN;
Float_t totalE;
Float_t liveCharge;
UShort_t trgType;
UShort_t trgTime;
Float_t clusterE[MAXC], clusterX[MAXC], clusterY[MAXC], clusterZ[MAXC];
Float_t clusterX1[MAXC], clusterY1[MAXC], clusterZ1[MAXC];
UInt_t clusterFlag[MAXC];
UInt_t clusterMatch[MAXC];
Short_t clusterNHit[MAXC];
Short_t clusterCID[MAXC];
UShort_t clusterTime[3][MAXC];

//for GEM
Float_t GEMX[MAXC], GEMY[MAXC], GEMZ[MAXC], GEMChargeX[MAXC], GEMChargeY[MAXC], GEMPeakX[MAXC], GEMPeakY[MAXC];
Short_t GEMSizeX[MAXC], GEMSizeY[MAXC], GEMID[MAXC];
Float_t ONGEMX[MAXC], ONGEMY[MAXC], ONGEMZ[MAXC];

vector<string> inputFiles;
vector<HyCalModule> moduleList;
TH2D* moduleECorrHist[2][T_BLOCKS];
TFile* moduleECorrFile[2];


int startRun = 0;
int endRun = 9999;
void LoadECorrHist();
float GetModuleECorr(float x, float y, int t);

float GetExpectedEnergy(int type, float Eb, float t)
{
    float theta = t/180.*TMath::Pi();
    float expectE = 0.;

    if (type == 0){
        expectE = Eb*targetMass[type] / ( Eb*(1.-cos(theta)) + targetMass[type] );
    }else{
        double C = cos(theta);
        expectE = targetMass[1]*(Eb + C*C*Eb + targetMass[1] - C*C*targetMass[1])/(Eb - C*C*Eb + targetMass[1] + C*C*targetMass[1]);
    }
    return expectE;
}
//___________________________________________________________________________________
int GetPrimExID(string & s){
    if (s[0] == 'W' || s[0] == 'G'){
        int id = stoi(s.substr(1, s.length() - 1));
        if (s[0] == 'W') id += 1000;
        return id;
    }else{
        return -1;
    }
}
//___________________________________________________________________________________
int GetRunNumber(string run)
{
    string sub = run.substr(strlen(run.c_str())-12 , 4);
    if (sub.at(0) == '0'){
        return stoi(sub.substr(1,4));
    }else{
        return stoi(sub);
    }
}
//__________________________________________________________________________________
bool SortFile(string name1, string name2) {
    return (GetRunNumber(name1) < GetRunNumber(name2)) ;
}
//__________________________________________________________________________________
bool SortForR(HyCalModule const & a, HyCalModule const & b){
    return (a.x*a.x + a.y*a.y) < (b.x*b.x + b.y*b.y);
}
//__________________________________________________________________________________
void FindInputFiles()
{
    inputFiles.clear();
    void* dirp = gSystem->OpenDirectory(inputDir.c_str());
    if (!dirp) {
        cerr << "Error: Can't open input file directory "<< inputDir << endl;
        return;
    }
    cout<<"Scanning for input files in "<< inputDir <<" "<<flush;
    const char* dir_item;
    while( (dir_item = gSystem->GetDirEntry(dirp)) ){
        if (!strncmp("prad_", dir_item, 5) && !strncmp("_sel.dst", dir_item+strlen(dir_item)-8, 8)){
        string thisName = inputDir;
        thisName.append("/");
        thisName.append(dir_item);
        int runNumber = GetRunNumber(thisName);
        if (runNumber <= endRun && runNumber >= startRun)
            inputFiles.push_back(thisName);
        }
    }
    cout<<endl<< "Found "<<inputFiles.size() <<" input files "<<endl;
    std::sort(inputFiles.begin(), inputFiles.end(), SortFile);
    gSystem->FreeDirectory(dirp);
}
//____________________________________________________________________________________

int main(int argc, char *argv[])
{
    if(argc != 3) {
        cout << "usage: physReplay start_run end_run ..." << endl;
        return 0;
    }

    // get the prad root directory
    string prad_root = getenv("PRAD_PATH");
    if(prad_root.size() && prad_root.back() != '/') prad_root += "/";

    startRun = stoi(argv[1]);
    endRun   = stoi(argv[2]);
    FindInputFiles();

    // initialize objects
    epics = new PRadEPICSystem(prad_root + "config/epics_channels.conf");
    gem = new PRadGEMSystem(prad_root + "config/gem.conf");
    hycal = new PRadHyCalSystem(prad_root + "config/hycal.conf");
    coord_sys = new PRadCoordSystem(prad_root + "database/coordinates.dat");
    det_match = new PRadDetMatch(prad_root + "config/det_match.conf");


    if (DoSShapeCorr){
        ConfigParser parser;
        if (!parser.OpenFile("./database/hycal_module.txt")){
            cout<<"cannot find offset database"<<endl;
            exit(0);
        }
        while(parser.ParseLine()){
            string name = parser.TakeFirst();
            string type = parser.TakeFirst();
            float input[6];
            for (int i=0; i<6; i++) input[i] = parser.TakeFirst().Float();

                HyCalModule thisModule(GetPrimExID(name), input[3], input[4], input[0], input[1]);
                moduleList.push_back(thisModule);
        }
        std::sort(moduleList.begin(), moduleList.end(), SortForR);
        parser.CloseFile();
        LoadECorrHist();
    }

    for(unsigned int i = 0; i < inputFiles.size(); ++i)
    {
        cout<<"analyzing file "<<inputFiles[i]<<endl;
        physReplay(inputFiles[i]);
    }

    return 0;
}
//___________________________________________________________________________________________
void physReplay(const string &path)
{
    hycal->ChooseRun(path);
    int runNumber = PRadInfoCenter::GetRunNumber();
    coord_sys->ChooseCoord(PRadInfoCenter::GetRunNumber());

    PRadDSTParser *dst_parser = new PRadDSTParser();
    dst_parser->OpenInput(path);

    string outPrefix = "physReplay";
    outPrefix += "_newAna";

    TFile* outFile = new TFile(Form("%s_%d.root", outPrefix.c_str(), runNumber), "RECREATE");
    TTree* t = new TTree("T", "T");
    InitTree(t);

    int count = 0;
    PRadBenchMark timer;
    PRadHyCalDetector *hycal_det = hycal->GetDetector();
    PRadGEMDetector *gem_det1 = gem->GetDetector("PRadGEM1");
    PRadGEMDetector *gem_det2 = gem->GetDetector("PRadGEM2");
    int beam_energy_ch = epics->GetChannel("MBSY2C_energy");


    while(dst_parser->Read())
    {
        if(dst_parser->EventType() == PRadDSTParser::Type::event) {
            auto event = dst_parser->GetEvent();

	        //if (event.trigger != LMS_Alpha) continue;

            // only interested in physics event
            if(!event.is_physics_event())
                continue;

            if((++count)%PROGRESS_COUNT == 0) {
                cout <<"------[ ev " << count << " ]---"
                     << "---[ " << timer.GetElapsedTimeStr() << " ]---"
                     << "---[ " << timer.GetElapsedTime()/(double)count << " ms/ev ]------"
                     << "\r" << flush;
            }

            //if (count == 10000) break;

            // update run information
            PRadInfoCenter::Instance().UpdateInfo(event);

            liveCharge = PRadInfoCenter::GetLiveBeamCharge();



            // reconstruct
            hycal->ChooseEvent(event);
            hycal->Reconstruct();
            gem->Reconstruct(event);

            // get reconstructed clusters
            auto &hycal_hit = hycal_det->GetHits();
            auto &gem1_hit = gem_det1->GetHits();
            auto &gem2_hit = gem_det2->GetHits();

            clusterN = hycal_hit.size() > MAXC ? MAXC : hycal_hit.size();

            //save the original HyCal recon position before projection and transform, used for energy s-shape correction
            vector<float> x_save;  x_save.clear();
            vector<float> y_save;  y_save.clear();
            vector<float> z_save;  z_save.clear();
            vector<int>  id_save; id_save.clear();

            for (int i=0; i<clusterN; i++){
                x_save.push_back(hycal_hit[i].x);
                y_save.push_back(hycal_hit[i].y);
                z_save.push_back(hycal_hit[i].z);
                id_save.push_back(hycal_hit[i].cid);
            }

            // coordinates transform, projection
            coord_sys->Transform(PRadDetector::HyCal, hycal_hit.begin(), hycal_hit.end());
            coord_sys->Transform(PRadDetector::PRadGEM1, gem1_hit.begin(), gem1_hit.end());
            coord_sys->Transform(PRadDetector::PRadGEM2, gem2_hit.begin(), gem2_hit.end());

            // hits matching, return matched index
            coord_sys->Projection(hycal_hit.begin(), hycal_hit.end());
            auto matched = det_match->Match(hycal_hit, gem1_hit, gem2_hit);

            // project to HyCal surface

            coord_sys->Projection(matched.begin(), matched.end());

            eventNumber = event.event_number;

            totalE = 0.;
            trgType = 0;
            trgTime = 65535;

            PRadTDCChannel * thisChannel = NULL;
            if (event.trigger == PHYS_LeadGlassSum){
            	//cout<<"LG trigger"<<endl;
            	thisChannel = hycal->GetTDCChannel("T1");
            	trgType = 1;
           	}
            else if (event.trigger == PHYS_TotalSum){
            	//cout<<"total sum trigger"<<endl;
            	thisChannel = hycal->GetTDCChannel("T2");
            	trgType = 2;
            }
            //assert (thisChannel != NULL);
            vector<unsigned short> time = thisChannel->GetTimeMeasure();
            if (time.size()!=0) trgTime = time[0];


            //save info for HyCal
            for (int i=0; i<clusterN; i++){
                clusterE[i]    = hycal_hit[i].E;
                clusterX[i]    = hycal_hit[i].x;
                totalE        += clusterE[i];
                clusterY[i]    = hycal_hit[i].y;
                clusterZ[i]    = hycal_hit[i].z;
                clusterFlag[i] = hycal_hit[i].flag;
                clusterNHit[i] = hycal_hit[i].nblocks;
                clusterCID[i]  = hycal_hit[i].cid;
                clusterTime[0][i] = hycal_hit[i].time[0];
                clusterTime[1][i] = hycal_hit[i].time[1];
                clusterTime[2][i] = hycal_hit[i].time[2];
                //init array for GEM
                GEMX[i] = 0;
                GEMY[i] = 0;
                GEMZ[i] = 0;
                ONGEMX[i] = 0;
                ONGEMY[i] = 0;
                ONGEMZ[i] = 0;
                GEMChargeX[i] = 0;
                GEMChargeY[i] = 0;
                GEMPeakX[i] = 0;
                GEMPeakY[i] = 0;
                GEMSizeX[i] = 0;
                GEMSizeY[i] = 0;
                GEMID[i] = -1;
                // matching flag
                clusterMatch[i] = 0;


                int index = -1;
                for (unsigned int j=0; j<id_save.size(); j++){
                    if (id_save[j] == hycal_hit[i].cid) {
                        index = j;
                        break;
                    }
                }

                if (index < 0) cout<<"warning: did not find the right index !!!!"<<endl;

                clusterX1[i] = x_save[index];
                clusterY1[i] = y_save[index];
                clusterZ1[i] = z_save[index];
            }

            for(auto &hit : matched)
            {
                GEMX[hit.hycal_idx]       = hit.x;
                GEMY[hit.hycal_idx]       = hit.y;
                GEMZ[hit.hycal_idx]       = hit.z;
                clusterMatch[hit.hycal_idx]= hit.mflag;

                GEMHit* thisHit = 0;
                if (TEST_BIT(hit.mflag, kGEM1Match)) thisHit = &hit.gem1.front();
                else if (TEST_BIT(hit.mflag, kGEM2Match))thisHit = &hit.gem2.front();

                GEMChargeX[hit.hycal_idx] = thisHit->x_charge;
                GEMChargeY[hit.hycal_idx] = thisHit->y_charge;
                GEMPeakX[hit.hycal_idx]   = thisHit->x_peak;
                GEMPeakY[hit.hycal_idx]   = thisHit->y_peak;
                GEMSizeX[hit.hycal_idx]   = thisHit->x_size;
                GEMSizeY[hit.hycal_idx]   = thisHit->y_size;
                GEMID[hit.hycal_idx]      = thisHit->det_id;
                ONGEMX[hit.hycal_idx]     = thisHit->x;
                ONGEMY[hit.hycal_idx]     = thisHit->y;
                ONGEMZ[hit.hycal_idx]     = thisHit->z;
            }

            if (DoSShapeCorr){
                for (int i=0; i<clusterN; i++){
                    float theta = atan(sqrt(clusterX[i]*clusterX[i] + clusterY[i]*clusterY[i])/HYCAL_Z)*180./TMath::Pi();
                    float expectE = GetExpectedEnergy(0, Ebeam, theta);
                    int regionType = 0;
                    if (TEST_BIT(clusterFlag[i], kTransition)) regionType = 1;
                    if (TEST_BIT(clusterFlag[i], kPbGlass))    regionType = 2;

                    float reso = 0.024;
                    if (regionType != 0) reso = 0.062;

                    int correctType = 0;

                    if (theta < 2.5){
                        float eCut = reso/sqrt(expectE/1000.)*expectE;

                        int cutSize = 4.;

                        if (regionType == 0) {
                            cutSize += 1.5;
                        }else{
                            cutSize += 4;
                        }


                        if (clusterE[i] <= expectE - cutSize*eCut) correctType = 1;
                    }else{
                        if (startRun > 1345 && clusterE[i] < 1500) correctType = 1;
                        if (endRun <= 1345 && clusterE[i] < 600) correctType = 1;
                    }


                    float factor = GetModuleECorr(clusterX1[i], clusterY1[i], correctType);
                    clusterE[i] /= factor;
                }
            }

            t->Fill();

        } else if(dst_parser->EventType() == PRadDSTParser::Type::epics) {
            auto epics_ev = dst_parser->GetEPICS();
	        // save epics into handler, otherwise get epicsvalue won't work
	        epics->AddEvent(epics_ev);
            // only update beam energy when there is an epics event
            Ebeam = epics_ev.values.at(beam_energy_ch);
        }
    }

    dst_parser->CloseInput();

    outFile->cd();
    t->Write();
    outFile->Close();

    delete dst_parser;
}
//___________________________________________________________________________________________
void InitTree(TTree* t)
{
    t->Branch("EventNumber",          &eventNumber,          "EventNumber/I"           );
    t->Branch("nHit",                 &clusterN,             "nHit/I"                  );
    t->Branch("EBeam",                &Ebeam,                "EBeam/F"                 );
    t->Branch("TotalE",               &totalE,               "TotalE/F"                );
    t->Branch("liveCharge",           &liveCharge,           "liveCharge/F"            );
    t->Branch("TriggerType",          &trgType,              "trgType/s"               );
    t->Branch("TriggerTime",          &trgTime,              "trgTime/s"               );
    t->Branch("Hit.Flag",             &clusterFlag[0],       "Hit.Flag[nHit]/i"        );
    t->Branch("Hit.Match",            &clusterMatch[0],      "Hit.Match[nHit]/i"       );
    t->Branch("Hit.X",                &clusterX[0],          "Hit.X[nHit]/F"           );
    t->Branch("Hit.Y",                &clusterY[0],          "Hit.Y[nHit]/F"           );
    t->Branch("Hit.Z",                &clusterZ[0],          "Hit.Z[nHit]/F"           );
    t->Branch("Hit.X1",               &clusterX1[0],         "Hit.X1[nHit]/F"          );
    t->Branch("Hit.Y1",               &clusterY1[0],         "Hit.Y1[nHit]/F"          );
    t->Branch("Hit.Z1",               &clusterZ1[0],         "Hit.Z1[nHit]/F"          );
    t->Branch("Hit.E",                &clusterE[0],          "Hit.E[nHit]/F"           );
    t->Branch("Hit.NModule",          &clusterNHit[0],       "Hit.NModule[nHit]/S"     );
    t->Branch("Hit.CID",              &clusterCID[0],        "Hit.CID[nHit]/S"         );
    t->Branch("Hit.Time.1",           &clusterTime[0][0],    "Hit.Time.1[nHit]/s"      );
    t->Branch("Hit.Time.2",           &clusterTime[1][0],    "Hit.Time.2[nHit]/s"      );
    t->Branch("Hit.Time.3",           &clusterTime[2][0],    "Hit.Time.3[nHit]/s"      );
    t->Branch("Hit.GEM.X",            &GEMX[0],              "Hit.GEM.X[nHit]/F"       );
    t->Branch("Hit.GEM.Y",            &GEMY[0],              "Hit.GEM.Y[nHit]/F"       );
    t->Branch("Hit.GEM.Z",            &GEMZ[0],              "Hit.GEM.Z[nHit]/F"       );

    t->Branch("Hit.GEM.Charge.X",     &GEMChargeX[0],        "Hit.GEM.Charge.X[nHit]/F");
    t->Branch("Hit.GEM.Charge.Y",     &GEMChargeY[0],        "Hit.GEM.Charge.Y[nHit]/F");
    t->Branch("Hit.GEM.Peak.X",       &GEMPeakX[0],          "Hit.GEM.Peak.X[nHit]/F"  );
    t->Branch("Hit.GEM.Peak.Y",       &GEMPeakY[0],          "Hit.GEM.Peak.Y[nHit]/F"  );
    t->Branch("Hit.GEM.Size.X",       &GEMSizeX[0],          "Hit.GEM.Size.X[nHit]/S"  );
    t->Branch("Hit.GEM.Size.Y",       &GEMSizeY[0],          "Hit.GEM.Size.Y[nHit]/S"  );
    t->Branch("Hit.GEM.ID",           &GEMID[0],             "Hit.GEM.ID[nHit]/S"      );
    t->Branch("Hit.ONGEM.X",          &ONGEMX[0],            "Hit.ONGEM.X[nHit]/F"     );
    t->Branch("Hit.ONGEM.Y",          &ONGEMY[0],            "Hit.ONGEM.Y[nHit]/F"     );
    t->Branch("Hit.ONGEM.Z",          &ONGEMZ[0],            "Hit.ONGEM.Z[nHit]/F"     );
}
//________________________________________________________________________________________
void LoadECorrHist()
{
    for (int i=0; i<T_BLOCKS; i++){
        moduleECorrHist[0][i] = 0;
        moduleECorrHist[1][i] = 0;
    }

    if (startRun > 1345){
        moduleECorrFile[0] = new TFile("./database/2GeV_ep_cali_total_25blocks.root", "READ");
        moduleECorrFile[1] = new TFile("./database/2GeV_ee_cali_total_25blocks.root", "READ");
    }else{
        moduleECorrFile[0] = new TFile("./database/1GeV_ep_cali_total_25blocks.root", "READ");
        moduleECorrFile[1] = new TFile("./database/1GeV_ee_cali_total_25blocks.root", "READ");
    }
    TIter next(moduleECorrFile[0]->GetListOfKeys());
    TKey *thisKey;
    while ( (thisKey = (TKey*)next()) ){
        TClass *thisClass = gROOT->GetClass(thisKey->GetClassName());
        if (!thisClass->InheritsFrom("TDirectory")) continue;
        TDirectory *thisDir = (TDirectory*)thisKey->ReadObj();
        string name = thisKey->GetName();
        int id = stoi(name);

        moduleECorrHist[0][id-1] = (TH2D*)thisDir->Get(Form("avg_ep_ratio_%04d", id));
    }

    TIter next1(moduleECorrFile[1]->GetListOfKeys());
    TKey *thisKey1;
    while ( (thisKey1 = (TKey*)next1()) ){
        TClass *thisClass = gROOT->GetClass(thisKey1->GetClassName());
        if (!thisClass->InheritsFrom("TDirectory")) continue;
        TDirectory *thisDir = (TDirectory*)thisKey1->ReadObj();
        string name = thisKey1->GetName();
        int id = stoi(name);

        moduleECorrHist[1][id-1] = (TH2D*)thisDir->Get(Form("avg_ee_ratio_%04d", id));
    }
}
//________________________________________________________________________________________
float GetModuleECorr(float x, float y, int t)
{
    float tx = 0;
    float ty = 0;
    int id = -1;

    for (unsigned int i=0; i<moduleList.size(); i++){
        if (moduleList[i].Contain(x, y)) {
            tx = (x - moduleList[i].x)/moduleList[i].sizex;
            ty = (y - moduleList[i].y)/moduleList[i].sizey;
            id = moduleList[i].cid;
            break;
        }
    }
    if (id <= 0) return 1.;

    if (moduleECorrHist[t][id-1] == 0) {
        cout<<"module ECorr hist does not exist for "<<id<<" "<<t<<endl;
        return 1.;
    }

    int binx = moduleECorrHist[t][id-1]->GetXaxis()->FindBin(tx);
    int biny = moduleECorrHist[t][id-1]->GetYaxis()->FindBin(ty);

    double factor = moduleECorrHist[t][id-1]->GetBinContent(binx, biny);

    if ( ! (fabs(tx) <0.3 && fabs(ty) < 0.3 ) ){
        double tdx = tx > 0. ? 0.5 - tx : tx - -0.5;
        double tdy = ty > 0. ? 0.5 - ty : ty - -0.5;
        double tdd = tdx > tdy ? tdy : tdx;

        if (tdd > 0.2) cout<<"something not right"<<" "<<tdd<<endl;

        double c1 = 1.004;
        double c2 = 0.996;

        if (id < 950) { c1 = 1.008; c2 = 0.992; }

        double addFactor = (c1-c2)/(-0.1)*(tdd-0.05) + c1;



        factor *= addFactor;


    }

    if (factor < 0.7 || factor > 1.3) {cout<<"unusual ECorr factor for "<<id<<" "<<factor<<" "<<binx<<" "<<biny<<endl; return 1.; }

    return factor;
}


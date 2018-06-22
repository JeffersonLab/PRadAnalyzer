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
#include "TSystem.h"
#include "TROOT.h"
#include "canalib.h"
using namespace std;

PRadEPICSystem *epics;
PRadHyCalSystem *hycal;
PRadGEMSystem *gem;
PRadCoordSystem *coord_sys;
PRadDetMatch *det_match;
void physReplay(const string &path);
void InitTree(TTree* t);
ostream &operator <<(ostream &os, const PRadBenchMark &timer);

#define MAXC 100
#define PROGRESS_COUNT 10000
Int_t eventNumber;
Float_t Ebeam;
//for HyCal
Int_t clusterN;
Float_t totalE;
Float_t clusterE[MAXC], clusterX[MAXC], clusterY[MAXC], clusterZ[MAXC];
UInt_t clusterFlag[MAXC];
UInt_t clusterMatch[MAXC];
Short_t clusterNHit[MAXC];
Short_t clusterCID[MAXC];
UShort_t clusterTime[3][MAXC];

//for GEM
Float_t GEMX[MAXC], GEMY[MAXC], GEMZ[MAXC], GEMChargeX[MAXC], GEMChargeY[MAXC], GEMPeakX[MAXC], GEMPeakY[MAXC];
Short_t GEMSizeX[MAXC], GEMSizeY[MAXC], GEMID[MAXC];
Float_t ONGEMX[MAXC], ONGEMY[MAXC], ONGEMZ[MAXC];

string inputDir = "../../../dst_files/background";
vector<string> inputFiles;
int startRun = 0;
int endRun = 9999;

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
        cout << "usage: physReplay <start_run> <end_run> ..." << endl;
        return -1;
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

    TFile* outFile = new TFile(Form("physReplay_%d.root", runNumber), "RECREATE");
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

            // only interested in physics event
            if(!event.is_physics_event())
                continue;

            if((++count)%PROGRESS_COUNT == 0) {
                cout <<"------[ ev " << count << " ]---"
                     << "---[ " << timer.GetElapsedTimeStr() << " ]---"
                     << "---[ " << timer.GetElapsedTime()/(double)count << " ms/ev ]------"
                     << "\r" << flush;
            }

            //if (count == 100000) break;

            // update run information
            PRadInfoCenter::Instance().UpdateInfo(event);

            // reconstruct
            hycal->ChooseEvent(event);
            hycal->Reconstruct();
            gem->Reconstruct(event);

            // get reconstructed clusters
            auto &hycal_hit = hycal_det->GetHits();
            auto &gem1_hit = gem_det1->GetHits();
            auto &gem2_hit = gem_det2->GetHits();

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

            clusterN = hycal_hit.size() > MAXC ? MAXC : hycal_hit.size();

            totalE = 0.;
            //save info for HyCal
            for (int i=0; i<clusterN; i++){
                clusterE[i]    = hycal_hit[i].E;
                totalE += clusterE[i];
                clusterX[i]    = hycal_hit[i].x;
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
    t->Branch("Hit.Flag",             &clusterFlag[0],       "Hit.Flag[nHit]/i"        );
    t->Branch("Hit.Match",            &clusterMatch[0],      "Hit.Match[nHit]/i"       );
    t->Branch("Hit.X",                &clusterX[0],          "Hit.X[nHit]/F"           );
    t->Branch("Hit.Y",                &clusterY[0],          "Hit.Y[nHit]/F"           );
    t->Branch("Hit.Z",                &clusterZ[0],          "Hit.Z[nHit]/F"           );
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


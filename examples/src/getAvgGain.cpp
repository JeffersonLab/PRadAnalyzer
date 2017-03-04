//============================================================================//
// Calculate the mean value of gain factor for each channel, for a given      //
// period based on the LMS root file produced during zero-suppression replay, //
// these files can be found on ifarm /work/hallb/prad/replay_LMS              //
//                                                                            //
// Weizhi Xiong                                                               //
// 10/23/2016                                                                 //
//============================================================================//

#include "PRadDataHandler.h"
#include "TFile.h"
#include "TH1.h"
#include "TMath.h"
#include "TAxis.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TKey.h"
#include "TClass.h"
#include "TF1.h"
#include <math.h>
#include <iostream>
#include <string>
#include <vector>
#include <cassert>
#include <iomanip>

#define BLOCKS 1728
#define NREFERENCE 3
#define TOLERANCE 0.06

using namespace std;
//global variables
string inputDir = "./calibration/LMS"; //define input root file here
string outDir   = "./calibration/LMS"; //define output dat file here
vector<string> inputFiles;
TFile* rootFile;
bool doAverage;
int startRun;
int endRun;
string channelName[BLOCKS];
double  pedMean[BLOCKS];
double  pedSigma[BLOCKS];
double  LMSMean[BLOCKS];
double  LMSSigma[BLOCKS];
double  refLMSMean[NREFERENCE];
double  refLMSSigma[NREFERENCE];
double  refPedMean[NREFERENCE];
double  refPedSigma[NREFERENCE];
double  refAlphaMean[NREFERENCE];
double  refGain[NREFERENCE];

//functions
void  InitArray();
int   GetRunNumber(string run);
bool  SortFile(string name1, string name2);
void  FindInputFiles();
void  FitHistograms();
void  FitHyCalHist(TKey* k, int & i);
void  FitReferenceHist(TKey* k);
double GetMaxBinCenterInRange(TH1* h, double xmin, double xmax);
void  WriteOutput(double nTotal, int rNumber = -1);

int main(int argc, char * argv [])
{
    char *ptr;
    doAverage = false;
    startRun = 0;
    endRun = 0;
    for(int i = 0; i < argc; ++i){
        ptr = argv[i];
        if(*(ptr++) == '-') {
        switch(*(ptr++)){
            case 'a':
                doAverage = true;
                break;
            case 's':
                startRun = stoi(argv[++i]);
                break;
            case 'e':
                endRun = stoi(argv[++i]);
                break;
            default:
                printf("Unkown option!\n");
                exit(1);
            }
        }
    }
    if (endRun == 0) endRun = startRun;
    if (endRun < startRun){
        cout<<"end Run less than start Run"<<endl;
        exit(1);
    }
    if (startRun == 0){
        cout<<"need to specify start run at least"<<endl;
        exit(1);
    }

    InitArray();

    FindInputFiles();//search for root files that within the directory

    FitHistograms();



}
//_________________________________________________________________________________
void InitArray()
{

    for (int i=0; i<BLOCKS; i++){
        pedMean[i]  = 0.;
        pedSigma[i] = 0.;
        LMSMean[i]  = 0.;
        LMSSigma[i] = 0.;
    }
    for (int i=0; i<NREFERENCE; i++){
        refLMSMean[i]   = 0.;
        refPedMean[i]   = 0.;
        refAlphaMean[i] = 0.;
        refLMSSigma[i]  = 0.;
        refPedSigma[i]  = 0.;
        refGain[i]      = 0.;
    }
}
//_________________________________________________________________________________
void FitHistograms()
{
    for (unsigned int i=0; i<inputFiles.size(); i++){
        rootFile = new TFile(inputFiles.at(i).c_str(), "READ");
        cout<<"analyzing file "<<inputFiles[i]<<endl;
        TIter next(rootFile->GetListOfKeys());
        TKey *thisKey;
        int ich = 0;
        while ( (thisKey = (TKey*)next()) ){
            TClass *thisClass = gROOT->GetClass(thisKey->GetClassName());
            if (!thisClass->InheritsFrom("TList")) continue;
            if ( (thisKey->GetName())[0] == 'W' || (thisKey->GetName())[0] == 'G' ){
                if (i==0) channelName[ich] = string(thisKey->GetName());

                FitHyCalHist(thisKey, ich);
                ich++;
            }
            if ( !strncmp("LMS", thisKey->GetName(), 3))
            FitReferenceHist(thisKey);
        }
        if (!doAverage) WriteOutput(1., GetRunNumber(inputFiles[i]) );
        delete rootFile; rootFile = 0;
    }
    if (doAverage) WriteOutput((double)inputFiles.size());
}
//_________________________________________________________________________________
void FitHyCalHist(TKey* k, int & i)
{
    if (!doAverage){
        pedMean[i]  = 0.; pedSigma[i] = 0.; LMSMean[i]  = 0.; LMSSigma[i] = 0.;
    }
    TList *thisList = (TList*)k->ReadObj();
    string LMSHistName(k->GetName());
    LMSHistName+="_LMS";
    string PedHistName(k->GetName());
    PedHistName+="_PED";
    TH1 *theLMSHist = (TH1*)thisList->FindObject(LMSHistName.c_str());
    TH1 *thePedHist = (TH1*)thisList->FindObject(PedHistName.c_str());
    //fit HyCal Module LMS
    TAxis *LMSAxis  = theLMSHist->GetXaxis();
    double integral  = theLMSHist->Integral(LMSAxis->FindBin(1),
                                           LMSAxis->FindBin(8190));
    if (integral > 1000){
        theLMSHist->Fit("gaus", "Qww", "", 1, 8190);
        TF1*  thisFunc = theLMSHist->GetFunction("gaus");
        double mean     = thisFunc->GetParameter(1);
        double sigma    = thisFunc->GetParameter(2);
        LMSMean[i]     = doAverage ? LMSMean[i]  + mean  : mean;
        LMSSigma[i]    = doAverage ? LMSSigma[i] + sigma : sigma;
        if (sigma/mean > TOLERANCE) cout<<"bad fitting for channel "<<channelName[i]<<endl;
    }
    //fit HyCal Module pedestal
    TAxis *pedAxis  = thePedHist->GetXaxis();
    int   binmax    = thePedHist->GetMaximumBin();
    double x         = pedAxis->GetBinCenter(binmax);
    integral        = thePedHist->Integral(pedAxis->FindBin(x - 30),
                                           pedAxis->FindBin(x + 30));
    if (integral > 1000){
        thePedHist->Fit("gaus", "Qww", "", x - 30, x + 30);
        TF1* thisFunc = thePedHist->GetFunction("gaus");
        double mean    = thisFunc->GetParameter(1);
        double sigma   = thisFunc->GetParameter(2);
        pedMean[i]    = doAverage ? pedMean[i]  + mean  : mean;
        pedSigma[i]   = doAverage ? pedSigma[i] + sigma : sigma;
        if (sigma/mean > TOLERANCE) cout<<"bad fitting for channel "<<channelName[i]<<endl;
    }

}
//_________________________________________________________________________________
void FitReferenceHist(TKey* k)
{
    TList *thisList = (TList*)k->ReadObj();
    int iref = ( (k->GetName())[3] - '0' ) - 1;
    assert(iref >=0 && iref < NREFERENCE);

    //if not averaging things over many runs, clean arrays before filling
    if (!doAverage){
        refLMSMean[iref] = 0.; refPedMean[iref] = 0.; refAlphaMean[iref] = 0.; refGain[iref] = 0.;
        refLMSSigma[iref] = 0.; refPedSigma[iref] = 0.;
    }
    string LMSHistName(k->GetName());
    LMSHistName+="_LMS";
    string PhysHistName(k->GetName());
    PhysHistName+="_PHYS";

    TH1 *theLMSHist = (TH1*)thisList->FindObject(LMSHistName.c_str());
    TH1 *thePhysHist = (TH1*)thisList->FindObject(PhysHistName.c_str());

    TAxis *physAxis = thePhysHist->GetXaxis();
    double integral1 = thePhysHist->Integral(physAxis->FindBin(1), physAxis->FindBin(1000));
    double integral2 = thePhysHist->Integral(physAxis->FindBin(1000), physAxis->FindBin(8190));

    TAxis *LMSAxis  = theLMSHist->GetXaxis();
    double integral3 = theLMSHist->Integral(LMSAxis->FindBin(1), LMSAxis->FindBin(8190));

    double up = 0, down = 0;
    if (integral1 > 1000){
        double center = GetMaxBinCenterInRange(thePhysHist, 1, 1000);
        TF1 *g1 = new TF1("g1", "gaus", center - 20, center + 20);
        thePhysHist->Fit(g1, "RQww");
        double mean = g1->GetParameter(1);
        double sigma = g1->GetParameter(2);
        refPedMean[iref]  += mean;
        refPedSigma[iref] += sigma;
        up = -1.*mean; down = -1.*mean;
        if (sigma/mean > TOLERANCE) cout<<"bad fitting for pedestal of reference PMT "<<iref + 1<<endl;
        delete g1;
    }

    if (integral2 > 1000){
        double center = GetMaxBinCenterInRange(thePhysHist, 1000, 8190);
        TF1 *g1 = new TF1("g1", "gaus", center - 300, center + 300);
        thePhysHist->Fit(g1, "RQww");
        double mean = g1->GetParameter(1);
        double sigma = g1->GetParameter(2);
        refAlphaMean[iref] += mean;
        down += mean;
        if (sigma/mean > TOLERANCE) cout<<"bad fitting for alpha of reference PMT "<<iref + 1<<endl;
        delete g1;
    }

    if (integral3 > 1000){
        theLMSHist->Fit("gaus", "Qww", "", 1, 8190);
        TF1*  thisFunc = theLMSHist->GetFunction("gaus");
        double mean     = thisFunc->GetParameter(1);
        double sigma    = thisFunc->GetParameter(2);
        refLMSMean[iref]  += mean;
        refLMSSigma[iref] += sigma;
        up += mean;
        if (sigma/mean > TOLERANCE) cout<<"bad fitting for LMS of reference PMT "<<iref + 1<<endl;
    }

    if (down > 0){
        //if trying to estimate the average, the following two ways are actually difference,
        //one is the many of all reference gains, the other is using the sum of all sparsified LMS
        //divided by all sparsified alpha, not sure which one is better though

        //refGain[iref] += ((refLMSMean[iref] - refPedMean[iref]) / (refAlphaMean[iref] - refPedMean[iref]));
        refGain[iref] += up/down;
    }
}
//_________________________________________________________________________________
void WriteOutput(double nTotal, int rNumber)
{
    if (rNumber < 0) assert(doAverage);
    ofstream outFile;
    string outFileName(outDir);
    outFileName += "/db_prad_baseinfo";
    if (doAverage){
        outFileName += "_"+to_string(startRun)+"_"+to_string(endRun)+".dat";
    }else{
        outFileName += "_"+to_string(rNumber)+".dat";
    }
    outFile.open(outFileName);

    for (int i=0; i<BLOCKS; i++){
        outFile<<setw(12)<<channelName[i]<<setw(12)<<pedMean[i]/nTotal
               <<setw(12)<<pedSigma[i]/nTotal<<setw(12)<<LMSMean[i]/nTotal
               <<setw(12)<<LMSSigma[i]/nTotal<<endl;
    }
    for (int i=0; i<NREFERENCE; i++){
        string thisName(Form("LMS%d", i+1));
        outFile<<setw(12)<<thisName<<setw(12)<<refPedMean[i]/nTotal
               <<setw(12)<<refPedSigma[i]/nTotal<<setw(12)<<refLMSMean[i]/nTotal
               <<setw(12)<<refLMSSigma[i]/nTotal<<endl;
    }
    outFile.close();

    ofstream gainFile;
    string gainFileName(outDir);
    gainFileName += "/prad_gain.dat";
    gainFile.open(gainFileName, std::ofstream::out | std::ofstream::app);

    string period;
    if (doAverage) period = Form("%d_%d", startRun, endRun);
    else period = Form("%d", rNumber);

    gainFile<<setw(12)<<period<<setw(12)<<refGain[0]/nTotal<<setw(12)<<
    refGain[1]/nTotal<<setw(12)<<refGain[2]/nTotal<<endl;
    gainFile.close();
}
//_________________________________________________________________________________
int GetRunNumber(string run)
{
    string sub = run.substr(strlen(run.c_str())-13 , 4);
    if (sub.at(0) == '0'){
        return stoi(sub.substr(1,4));
    }else{
        return stoi(sub);
    }
}
//__________________________________________________________________________________
bool SortFile(string name1, string name2)
{
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
        if (!strncmp("prad_", dir_item, 5) && !strncmp("LMS.root", dir_item+strlen(dir_item)-8, 8)){
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
//__________________________________________________________________________________
double GetMaxBinCenterInRange(TH1* h, double xmin, double xmax)
{
    double max = 0;
    double save = 0;
    for (int i=1; i<h->GetXaxis()->GetNbins(); i++){
        if (h->GetBinCenter(i) < xmin || h->GetBinCenter(i) > xmax) continue;
        if (h->GetBinContent(i) > save) {
            save = h->GetBinContent(i);
            max  = h->GetBinCenter(i);
        }
    }
    return max;
}

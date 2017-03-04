#include "PRadDSTParser.h"
#include "PRadEvioParser.h"
#include "PRadBenchMark.h"
#include "PRadHyCalSystem.h"
#include "PRadGEMSystem.h"
#include "PRadEPICSystem.h"
#include "PRadCoordSystem.h"
#include "PRadDetMatch.h"
#include "PRadInfoCenter.h"
#include "PRadBenchMark.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"

#include <iostream>
#include <string>
#include <vector>

#define T_BLOCKS 2156
#define PROGRESS_COUNT 10000
#define MAXCLUSTER 6

using namespace std;

const bool isLowE = true; //true if 1.1 GeV

struct PairCoor
{
  float x1;
  float y1;
  float x2;
  float y2;
  PairCoor(float &xx1, float& yy1, float& xx2, float& yy2)
  :x1(xx1), y1(yy1), x2(xx2), y2(yy2) {}
};

enum ParticleType
{
    kProton = 0,
    kElectron
};

// combine gem hit info to hycal hit
class CombinedHit : public HyCalHit
{
public:
    float x_gem;
    float y_gem;

    // be abled to assgin HyCalHit to this class
    CombinedHit &operator =(const HyCalHit &hit)
    {
        // assign hit to the base class HyCalHit
        HyCalHit::operator =(hit);
        return *this;
    }

    // be abled to assgin GEMHit to this class
    CombinedHit &operator =(const GEMHit &hit)
    {
        x_gem = hit.x;
        y_gem = hit.y;
        return *this;
    }
};

//global variables
int clusterN, eventNumber;
int currentRunNumber;
float Ebeam, HyCalZ;
float offsetGEM[4];
const float targetMass[2] = {938.272046, 0.510998928};
const float innerBoundary[2] = {83.079,     83};
vector< PairCoor > pairHyCal;
vector< PairCoor > pairGEM;

PRadHyCalSystem *hycal_sys;
PRadGEMSystem *gem_sys;
PRadEPICSystem *epics;
PRadCoordSystem *coord_sys;
PRadDetMatch *det_match;
PRadDSTParser *dst_parser;
PRadHyCalDetector *hycal;
PRadGEMDetector *gem1;
PRadGEMDetector *gem2;
const EventData *current_event;

vector<int> innerModule;
vector<int> boundModule;
vector<int> countFill;
vector<TH2F*> profile;
float moduleEnergy[64][12][12] = {0.};
int innerModuleList[12] = {1526, 1527, 1528, 1529, 1560, 1563, 1594, 1597, 1628, 1629, 1630, 1631};
//int innerModuleList[12] = {1458, 1459, 1460, 1461, 1696, 1697, 1698, 1699, 1531, 1565, 1599, 1633};
float pwo_profile[501][501] = {0.};
//functions
void Helper();
void InitHistogram();
void InitInnerModule();
void InitProfileData();
vector<string> FindInputFiles(const string &in_dir, int start, int end);
int  GetRunNumber(string run);
bool SortFile(string name1, string name2);
bool SortForE(CombinedHit const & a, CombinedHit const & b);
void InnerModuleAnalyzer(CombinedHit* h);
void FillInnerHist(int& id, float& E, float& expectE, ParticleType type);
void EPAnalyzer(CombinedHit * h);
void MollerAnalyzer(CombinedHit * h, int index1 = 0, int index2 = 1);
void MultiHitAnalyzer(CombinedHit* h, int &n);
void OffsetAnalyzer(CombinedHit* h);
void LinesIntersect(float &xsect, float &ysect, float &xsect_GEM, float & ysect_GEM);
void GetIntersect(float& x1, float& x2, float& y1,
                  float& y2, float& x, float& y);//for monitoring beam spot
float GetExpectedEnergy(ParticleType type, float& x, float &y);
float GetExpectedEFromProfile(float& dx, float& dy);
float GetElossIonElectron(float &theta, float& E);
int   IsNeighborToInner(int& id);
bool MatchedGEM(const CombinedHit &hit);

//histograms
TH1F *h_beam_e;
TH2F *ep_x_y;
TH2F *GEM_x_y;
TH2F *GEM_ep_x_y;
TH1F *ep_ratio_all;
TH2F *ep_angle_E;

TH2F *ee1_x_y;
TH2F *ee2_x_y;
TH2F *sym_ee_x_y;
TH1F *sym_ee_E;
TH1F *ee_ratio_all;
TH2F *ee_angle_E;
TH2F *eloss;
TH2F *total_angle_E;

TH1F *deltaCoor[2];
TH1F *deltaCoor_GEM[2];
TH1F *coorIntersect[2];
TH1F *coorIntersect_GEM[2];
TH1F *sym_ee_r[2];

TH1F *ep_ratio[T_BLOCKS];
TH1F *ee_ratio[T_BLOCKS];
TH1F *ep_energy[T_BLOCKS];
TH1F *ee_energy[T_BLOCKS];

TH1F *inner_ratio[12];
TH1F *inner_deltax[12];
TH1F *inner_deltay[12];
TH2F *ratio_x[12];
TH2F *ratio_y[12];


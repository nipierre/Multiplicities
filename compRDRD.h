#include <vector>
#include <string>
#include <set>
#include <map>
#include <utility>
#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TMatrixD.h>
#include <TMatrixTUtils.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TLine.h>
#include <TGaxis.h>
#include <TLatex.h>

// COLORS

#define RST  "\x1B[0m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KYEL  "\x1B[33m"
#define KBLU  "\x1B[34m"
#define KMAG  "\x1B[35m"
#define KCYN  "\x1B[36m"
#define KWHT  "\x1B[37m"

#define FRED(x) KRED x RST
#define FGRN(x) KGRN x RST
#define FYEL(x) KYEL x RST
#define FBLU(x) KBLU x RST
#define FMAG(x) KMAG x RST
#define FCYN(x) KCYN x RST
#define FWHT(x) KWHT x RST

#define BOLD(x) "\x1B[1m" x RST
#define UNDL(x) "\x1B[4m" x RST

using namespace std;

//Qty containers

vector<Double_t> fQ2;
vector<Double_t> fQ2kin[5];
vector<Double_t> fQ2local;
vector<Double_t> fXBj;
vector<Double_t> fXBjkin[5];
vector<Double_t> fYBj;
vector<Double_t> fYBjkin[5];
vector<Double_t> fWBj;
vector<Double_t> fWBjkin[5];
vector<Double_t> fNu;
vector<Double_t> fNukin[5];
vector<Double_t> fMu[5];
vector<Double_t> fMup[5];
vector<Double_t> fThetaMu[3];
vector<Double_t> fTheta[5];
vector<Double_t> fPhi[5];
vector<Double_t> fVertex[5];

vector<Double_t> fQ2RD2;
vector<Double_t> fQ2kinRD2[5];
vector<Double_t> fXBjRD2;
vector<Double_t> fXBjkinRD2[5];
vector<Double_t> fYBjRD2;
vector<Double_t> fYBjkinRD2[5];
vector<Double_t> fWBjRD2;
vector<Double_t> fWBjkinRD2[5];
vector<Double_t> fNuRD2;
vector<Double_t> fNukinRD2[5];
vector<Double_t> fMuRD2[5];
vector<Double_t> fMupRD2[5];
vector<Double_t> fThetaRD2Mu[3];
vector<Double_t> fThetaRD2[5];
vector<Double_t> fPhiRD2[5];
vector<Double_t> fVertexRD2[5];

vector<Double_t> fX;
vector<Double_t> fY;
vector<Double_t> fZ;
vector<Double_t> fXRD2;
vector<Double_t> fYRD2;
vector<Double_t> fZRD2;

//Misc

set<Double_t> fLHsec_set;
Double_t* fLHsec_tab;
Double_t fLHsec;
int fId;
int fId_loose;
int fId_severe;
Double_t fNu_max[3][12];
Double_t fNu_min[3][12];
vector<Double_t> fXv;
vector<Double_t> fYv;
vector<Double_t> fZv;
vector<Double_t> fRv;

bool fAllDISflag;
bool fAllDISflag_RD2;


int xbin, ybin, zbin;

// Cuts

Double_t fXmin;
Double_t fXmax;
Double_t fYmin;
Double_t fYmax;
Double_t fWmin;
Double_t fWmax;
Double_t fPmin;
Double_t fPmax;


//Draw

TH1F* fKinematicsRD[5][13];
TH1F* fKinematicsRD2[5][13];
TH1F* fKinematicsRatio[5][13];
TCanvas c1("c1","",3200,1600);
TCanvas c2("c2","",3200,1600);
TCanvas c3("c3","",3200,1600);
TCanvas c4("c4","",3200,1600);
TCanvas c5("c5","",3200,1600);
TCanvas c6("c6","",3200,1600);
TCanvas c8("c8","",1600,1600);
TCanvas c9("c9","",1600,1600);
TCanvas c10("c10","",1600,1600);
TCanvas c11("c11","",1600,1600);
TCanvas c12("c12","",1600,1600);
TCanvas c13("c13","",1600,1600);
TCanvas c14("c14","",3200,1600);
TCanvas c15("c15","",3200,1600);
TCanvas c16("c16","",3200,1600);
TCanvas c17("c17","",3200,1600);
TCanvas c18("c18","",3200,1600);
TCanvas c19("c19","",3200,1600);
TCanvas c20("c20","",1600,1600);
TCanvas c21("c21","",3200,1600);
TCanvas c22("c22","",1600,1600);

vector<double> fError, fErrorRD, fErrorRD2;
Int_t fLineStyle[7] = {3,3,3,1,3,3,3};

TLine* l1[17][7];

//Graphic Style

Int_t fMarkerColor[2] = {4,2};
Int_t fMarkerStyle[2] = {24,20};

//Constants

static const Double_t fM_p = 938.272046/(1e3);
static const Double_t fM_mu = 105.6583715/(1e3);
static const Double_t fM_K = 493.677/(1e3);
static const Double_t fM_pi = 139.57018/(1e3);

string trigname[5] = {"MT","LT","OT","LAST",""};

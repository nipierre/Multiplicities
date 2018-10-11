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
TH2F* fThetaRDp[3];
TH1F* fKinematicsRD2[5][13];
TH2F* fThetaRD2p[3];
TH1F* fKinematicsRatio[5][13];
TCanvas c1("Kin_Q2 Trigger","Kin_Q2 Trigger",3200,1600);
TCanvas c2("Kin_xBj Trigger","Kin_xBj Trigger",3200,1600);
TCanvas c3("Kin_y Trigger","Kin_y Trigger",3200,1600);
TCanvas c4("Kin_z Trigger","Kin_z Trigger",3200,1600);
TCanvas c5("Kin_w Trigger","Kin_w Trigger",3200,1600);
TCanvas c6("Kin_nu Trigger","Kin_nu Trigger",3200,1600);
TCanvas c8("Kin_Q2","Kin_Q2",1600,1600);
TCanvas c9("Kin_xBj","Kin_xBj",1600,1600);
TCanvas c10("Kin_y","Kin_y",1600,1600);
TCanvas c11("Kin_z","Kin_z",1600,1600);
TCanvas c12("Kin_w","Kin_w",1600,1600);
TCanvas c13("Kin_nu","Kin_nu",1600,1600);
TCanvas c14("Kin_mu","Kin_mu",3200,1600);
TCanvas c15("Kin_mup","Kin_mup",3200,1600);
TCanvas c16("Kin_theta","Kin_theta",3200,1600);
TCanvas c17("Kin_phi","Kin_phi",3200,1600);
TCanvas c18("Kin_vertex","Kin_vertex",3200,1600);
TCanvas c19("Kin_hadron_p Trigger","Kin_hadron_p Trigger",3200,1600);
TCanvas c20("Kin_hadron_p","Kin_hadron_p",1600,1600);
TCanvas c21("Kin_p_T Trigger","Kin_hadron_p_T Trigger",3200,1600);
TCanvas c22("Kin_p_T","Kin_hadron_p_T",1600,1600);

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
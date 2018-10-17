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

vector<Double_t> fQ2_MC;
vector<Double_t> fQ2kinMC[5];
vector<Double_t> fXBj_MC;
vector<Double_t> fXBjkinMC[5];
vector<Double_t> fYBj_MC;
vector<Double_t> fYBjkinMC[5];
vector<Double_t> fWBj_MC;
vector<Double_t> fWBjkinMC[5];
vector<Double_t> fNu_MC;
vector<Double_t> fNukinMC[5];
vector<Double_t> fMuMC[5];
vector<Double_t> fMupMC[5];
vector<Double_t> fThetaMCMu[3];
vector<Double_t> fThetaMC[5];
vector<Double_t> fPhiMC[5];
vector<Double_t> fVertexMC[5];

vector<Double_t> fX;
vector<Double_t> fY;
vector<Double_t> fZ;
vector<Double_t> fXMC;
vector<Double_t> fYMC;
vector<Double_t> fZMC;

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
bool fAllDISflag_MC;


int xbin, ybin, zbin, xbin_MC, ybin_MC, zbin_u;

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

TH1F* fKinematicsRD[5][20];
TH2F* fThetaRDp[3];
TH1F* fKinematicsMC[5][20];
TH2F* fThetaMCp[3];
TH1F* fKinematicsRatio[5][20];
Int_t fCountingRD[5][20];
Int_t fCountingMC[5][20];
TCanvas c1("Kin_{Q2} Trigger","Kin_{Q2} Trigger",3200,1600);
TCanvas c2("Kin_x^{Bj} Trigger","Kin_x^{Bj} Trigger",3200,1600);
TCanvas c3("Kin_y Trigger","Kin_y Trigger",3200,1600);
TCanvas c4("Kin_z Trigger","Kin_z Trigger",3200,1600);
TCanvas c5("Kin_w Trigger","Kin_w Trigger",3200,1600);
TCanvas c6("Kin_nu Trigger","Kin_nu Trigger",3200,1600);
TCanvas c7("Kin_Phi Trigger","Kin_Phi Trigger",3200,1600);
TCanvas c8("Kin_Q^{2}","Kin_Q^{2}",3200,1600);
TCanvas c9("Kin_x^{Bj}","Kin_x^{Bj}",3200,1600);
TCanvas c10("Kin_y","Kin_y",3200,1600);
TCanvas c11("Kin_z","Kin_z",3200,1600);
TCanvas c12("Kin_w","Kin_w",3200,1600);
TCanvas c13("Kin_nu","Kin_nu",3200,1600);
TCanvas c14("Kin_mu","Kin_mu",3200,1600);
TCanvas c15("Kin_mup","Kin_mup",3200,1600);
TCanvas c16("Kin_theta","Kin_theta",3200,1600);
TCanvas c17("Kin_phi","Kin_phi",3200,1600);
TCanvas c18("Kin_vertex","Kin_vertex",3200,1600);
TCanvas c19("Kin_hadron_p Trigger","Kin_hadron_p Trigger",3200,1600);
TCanvas c20("Kin_hadron_theta Trigger","Kin_hadron_theta Trigger",3200,1600);
TCanvas c21("Kin_hadron_phi Trigger","Kin_hadron_phi Trigger",3200,1600);
TCanvas c22("Kin_hadron_p","Kin_hadron_p",3200,1600);
TCanvas c23("Kin_hadron_theta","Kin_hadron_theta",3200,1600);
TCanvas c24("Kin_hadron_phi","Kin_hadron_phi",3200,1600);
TCanvas c25("Kin_hadron_phi_pl Trigger","Kin_hadron_phi_pl Trigger",3200,1600);
TCanvas c26("Kin_hadron_phi_pl","Kin_hadron_phi_pl",3200,1600);
TCanvas c27("Kin_p_T Trigger","Kin_hadron_p_T Trigger",3200,1600);
TCanvas c28("Kin_p_T","Kin_hadron_p_T",3200,1600);
TCanvas c29("photoelec","photoelec",2000,2000);
TCanvas c30("Q^2","Q^2",2000,2000);
TCanvas c31("nu","nu",2000,2000);
TCanvas c32("z","z",2000,2000);
TCanvas c33("E_mu","E_mu",3200,1600);
TCanvas c34("Thetay_mu","Thetay_mu",3200,1600);
TCanvas c35("Thetax_mu","Thetax_mu",3200,1600);
TCanvas c36("Thetaxy_mu","Thetaxy_mu",3200,1600);
TCanvas c37("Vertex (Hadron+e) Trig","Vertex (Hadron+e) Trig",3200,1600);
TCanvas c38("Vertex (Hadron+e)","Vertex (Hadron+e)",3200,1600);
TCanvas c39("Vertex (e) Trig","Vertex (e) Trig",3200,1600);
TCanvas c40("Vertex (e)","Vertex (e)",3200,1600);
TCanvas c41("Vertex (Hadron) Trig","Vertex (Hadron) Trig",3200,1600);
TCanvas c42("Vertex (Hadron)","Vertex (Hadron)",3200,1600);

vector<double> fError, fErrorRD, fErrorMC;
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

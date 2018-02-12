#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <set>
#include <map>
#include <utility>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TMatrixD.h>
#include <TMatrixTUtils.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TGraphErrors.h>
#include <TLine.h>

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

//Structs

struct Wrapper { Double_t tab[2][2][4]; };
struct Pvsz { vector<Double_t> vec[2][5]; };
struct hadiden { vector<Double_t> vec; };
struct studyxy { vector<Double_t> vec[2]; };
struct Recovery { vector<Double_t> vec[2][4][4]; };
struct Recovery_tank { Double_t tab[2][4][4]; };


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
vector<Double_t> fThetaMC[5];
vector<Double_t> fPhiMC[5];
vector<Double_t> fVertexMC[5];

vector<Double_t> fX;
vector<Double_t> fY;
vector<Double_t> fZ;
vector<Double_t> fXMC;
vector<Double_t> fYMC;
vector<Double_t> fZMC;

vector<Pvsz> fPvsz;
vector<Pvsz> fPvsz_err;
vector<hadiden> fHadiden;

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
hadiden fRstudy[24];
studyxy fRstudy_xy[24];

bool fAllDISflag;
bool fAllDISflag_MC;

//Counting

Double_t fBP = 0;
Double_t fRmu = 0;
Double_t fBMS = 0;
Double_t fBEC = 0;
Double_t fTarg = 0;
Double_t fCell = 0;
Double_t fTrig = 0;
Double_t fQ2test = 0;
Double_t fYBjtest = 0;
Double_t fXBjtest = 0;
Double_t fWBjtest = 0;
Double_t fXX0test = 0;
Double_t fMom = 0;
Double_t fTRICH = 0;
Double_t fPosRICH = 0;
Double_t fHplus = 0;
Double_t fHminus = 0;
Double_t fPiplus = 0;
Double_t fPiminus = 0;
Double_t fKplus = 0;
Double_t fKminus = 0;
Double_t fPplus = 0;
Double_t fPminus = 0;
Double_t fMCHplus = 0;
Double_t fMCHminus = 0;
Double_t fMCPiplus = 0;
Double_t fMCPiminus = 0;
Double_t fMCKplus = 0;
Double_t fMCKminus = 0;
Double_t fMCPplus = 0;
Double_t fMCPminus = 0;
Double_t fPiplus_true = 0;
Double_t fPiminus_true = 0;
Double_t fKplus_true = 0;
Double_t fKminus_true = 0;
Double_t fPplus_true = 0;
Double_t fPminus_true = 0;
Double_t fPiplus_err = 0;
Double_t fPiminus_err = 0;
Double_t fKplus_err = 0;
Double_t fKminus_err = 0;
Double_t fPplus_err = 0;
Double_t fPminus_err = 0;

//Binning
Wrapper fMC1[9][5][12];
Wrapper fMC2[9][5][12];
Wrapper fMultRatio[9][5][12];
Double_t fNDIS_evt_MC1[9][5][12];
Double_t fNDIS_evt_MC2[9][5][12];
int xbin, ybin, zbin, xbin_MC, ybin_MC, zbin_u;
Double_t fZrange[13] = {.20,.25,.30,.35,.40,.45,.50,.55,.60,.65,.70,.75,.85};
Double_t fXrange[10] = {.004,.01,.02,.03,.04,.06,.1,.14,.18,.4};
Double_t fYrange[6] = {.1,.15,.2,.3,.5,.7};
Double_t fRcutval[24] = {1.68812,1.6915,1.69572,1.69733,1.71178,1.74735,1.74682,1.7846,1.80058,1.81382,1.83367,1.84183,1.84587,1.8423,1.8376,1.8368,1.84023,1.84309,1.85645,1.86316,1.85021,1.84775,1.84463,1.84185};
int fFlag[3][9][5][12];
int fFlag_MC[3][9][5][12];
Double_t fZ_bin_width[12] = {.05,.05,.05,.05,.05,.05,.05,.05,.05,.05,.05,.1};

Recovery fMeanvalues[9][5][12]; // tab[][][i], iC[0,3] : x,y,Q2,z
Recovery_tank fMeanvalues_size[9][5][12];
Recovery_tank fMeanvalues_data[9][5][12];

Wrapper fBinning[9][5][12];
Wrapper fBinning_loose[9][5][12];
Wrapper fBinning_severe[9][5][12];

studyxy fR_xy[24];

//Draw

TH1F* fKinematicsMC1[5][12];
TH1F* fKinematicsMC2[5][12];
TH1F* fKinematicsRatio[5][12];
TCanvas c1("Kin_Q^{2} Trigger","Kin_Q^{2} Trigger",3200,1600);
TCanvas c2("Kin_x^{Bj} Trigger","Kin_x^{Bj} Trigger",3200,1600);
TCanvas c3("Kin_y Trigger","Kin_y Trigger",3200,1600);
TCanvas c4("Kin_z Trigger","Kin_z Trigger",3200,1600);
TCanvas c5("Kin_w Trigger","Kin_w Trigger",3200,1600);
TCanvas c6("Kin_#nu Trigger","Kin_#nu Trigger",3200,1600);
TCanvas c7("Kin_#Phi Trigger","Kin_#Phi Trigger",3200,1600);
TCanvas c8("Kin_Q^{2}","Kin_Q^{2}",3200,1600);
TCanvas c9("Kin_x^{Bj}","Kin_x^{Bj}",3200,1600);
TCanvas c10("Kin_y","Kin_y",3200,1600);
TCanvas c11("Kin_z","Kin_z",3200,1600);
TCanvas c12("Kin_w","Kin_w",3200,1600);
TCanvas c13("Kin_#nu","Kin_#nu",3200,1600);
TCanvas c14("Kin_#mu","Kin_#mu",3200,1600);
TCanvas c15("Kin_#mup","Kin_#mup",3200,1600);
TCanvas c16("Kin_#theta","Kin_#theta",3200,1600);
TCanvas c17("Kin_#phi","Kin_#phi",3200,1600);
TCanvas c18("Kin_vertex","Kin_vertex",3200,1600);
TCanvas c19("MultRatio Hadron","MultRatio Hadron",3200,1600);
TCanvas c20("MultRatio Pion","MultRatio Hadron",3200,1600);
TCanvas c21("MultRatio Kaon","MultRatio Hadron",3200,1600);

TLine* l1[5];

//Graphic Style

Int_t fMarkerColor[2] = {4,2};
Int_t fMarkerStyle[2] = {24,20};
Int_t fMarkerColorAlt[5] = {2,95,209,226,221};
Int_t fMarkerStyleAlt[5][2] = {{24,20},{26,22},{25,21},{27,33},{28,34}};

//Constants

static const Double_t fM_p = 938.272046/(1e3);
static const Double_t fM_mu = 105.6583715/(1e3);
static const Double_t fM_K = 493.677/(1e3);
static const Double_t fM_pi = 139.57018/(1e3);

string trigname[5] = {"MT","LT","OT","LAST","All Trig"};

//----------------------------------------------------------------------------
//--------- RICH unfolding ---------------------------------------------------
//----------------------------------------------------------------------------

TMatrixD rich_mat_p[2][10];
TMatrixD rich_mat_m[2][10];
TMatrixD inv_rich_p[2][10];
TMatrixD inv_rich_m[2][10];
TMatrixD err_rich_p[3][3];
TMatrixD err_rich_m[3][3];
Double_t pi_unfolding_err_p[2][10][3];
Double_t pi_unfolding_err_m[2][10][3];
Double_t k_unfolding_err_p[2][10][3];
Double_t k_unfolding_err_m[2][10][3];
Double_t p_unfolding_err_p[2][10][3];
Double_t p_unfolding_err_m[2][10][3];
Int_t mat_bin[2][10];
Int_t err_bin[2][10];
Double_t cov1_pi[2][3];
Double_t cov1_k[2][3];
Double_t cov1_p[2][3];
Double_t cov2[2][3];
TMatrixD pi_sigma_uni(3,3);
TMatrixD k_sigma_uni(3,3);
TMatrixD p_sigma_uni(3,3);
TMatrixD pi_vect(3,1);
TMatrixD k_vect(3,1);
TMatrixD p_vect(3,1);
Double_t dummy;

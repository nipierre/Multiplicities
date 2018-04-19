#ifndef __CUTMAKER_SPLIT_H_INCLUDED__
#define __CUTMAKER_SPLIT_H_INCLUDED__

#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TMatrixD.h>
#include <TMatrixTUtils.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TGraphErrors.h>
#include <TGaxis.h>
#include <vector>
#include <set>
#include <map>
#include <utility>

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

struct Wrapper { Double_t tab[2][2][4]; };
struct Multiplicities { Double_t tab[2][3][4]; };
struct Recovery { vector<Double_t> vec[2][4][4]; };
struct Recovery_tank { Double_t tab[2][4][4]; };
struct Pvsz { vector<Double_t> vec[2][5]; };
struct hadiden { vector<Double_t> vec; };
struct studyxy { vector<Double_t> vec[2]; };

//DISevt
Int_t            fRunNo;
vector<Int_t>    fSpillNo;
vector<Int_t>    fEvtInSpill;
vector<bool>     fTrigMask;
vector<Int_t>    fEvNo;
vector<Double_t> fX;
vector<Double_t> fY;
vector<Double_t> fZ;
vector<Double_t> fP0x;
vector<Double_t> fP0y;
vector<Double_t> fP0z;
vector<Double_t> fP1x;
vector<Double_t> fP1y;
vector<Double_t> fP1z;
vector<Double_t> fE_beam;
vector<Double_t> fE_mu_prim;
vector<Double_t> fXX0;
vector<Double_t> fHM04x;
vector<Double_t> fHM04y;
vector<Double_t> fHM05x;
vector<Double_t> fHM05y;
vector<Double_t> fHO03x;
vector<Double_t> fHO03y;
vector<Double_t> fHO04x;
vector<Double_t> fHO04y;
vector<bool>     fSaved;
vector<bool>     fCellCrossed;
vector<bool>     fBackPropFlag;

//Hadrons
vector<Double_t> fP;
vector<Double_t> fTh;
vector<Double_t> fPh;
vector<Double_t> fHXX0;
vector<bool>     fInHCALacc;
vector<Double_t> fHCAL;
vector<Int_t>    fCharge;
vector<Double_t> fThRICH;
vector<Double_t> fLH[6];
vector<Int_t>    fMCpid;
vector<Double_t> fMM01x;
vector<Double_t> fMM01y;
vector<Double_t> fMM02x;
vector<Double_t> fMM02y;
vector<Double_t> fMM03x;
vector<Double_t> fMM03y;
vector<Double_t> fZ2Ax;
vector<Double_t> fZ2Ay;
vector<Double_t> fZ2Bx;
vector<Double_t> fZ2By;
vector<Double_t> fRICHx;
vector<Double_t> fRICHy;
vector<Double_t> fAtsize;

vector<Double_t> fQ2;
vector<Double_t> fQ2kin;
vector<Double_t> fQ2local;
vector<Double_t> fXBj;
vector<Double_t> fXBjkin;
vector<Double_t> fYBj;
vector<Double_t> fYBjkin;
vector<Double_t> fWBj;
vector<Double_t> fWBjkin;
vector<Double_t> fNu;
vector<Double_t> fNukin;
vector<Double_t> fLHpi;
vector<Double_t> fLHK;

vector<Pvsz> fPvsz;
vector<Pvsz> fPvsz_err;
vector<hadiden> fHadiden;

//Misc

set<Double_t> fLHsec_set;
Double_t* fLHsec_tab;
Double_t fLHsec;
Int_t fId;
Int_t fId_loose;
Int_t fId_severe;
Double_t fNu_max[3][12];
Double_t fNu_min[3][12];
vector<Double_t> fXv;
vector<Double_t> fYv;
vector<Double_t> fZv;

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

Double_t fInclusiveRCproton[30][19];
Double_t fSemiInclusiveRCproton[9][6][14];

//Binning
Wrapper fBinning[9][5][12];
Wrapper fBinning_loose[9][5][12];
Wrapper fBinning_severe[9][5][12];
Multiplicities fMultiplicities[9][5][12];
Wrapper fAcceptance[9][5][12];
Wrapper fRich_sys_err[9][5][12]; // tab[][0][] : stat, tab[][1][] : sys
Recovery fMeanvalues[9][5][12]; // tab[][][i], iC[0,3] : x,y,Q2,z
Recovery_tank fMeanvalues_size[9][5][12];
Recovery_tank fMeanvalues_data[9][5][12];
Double_t fNDIS_evt[3][9][5][12];
Double_t fNDIS_evt_err[3][9][5][12];
Int_t xbin, ybin, zbin;
Double_t fZrange[13] = {.20,.25,.30,.35,.40,.45,.50,.55,.60,.65,.70,.75,.85};
Double_t fXrange[10] = {.004,.01,.02,.03,.04,.06,.1,.14,.18,.4};
Double_t fYrange[6] = {.1,.15,.2,.3,.5,.7};
Double_t fRcutval[24] = {1.68812,1.6915,1.69572,1.69733,1.71178,1.74735,1.74682,1.7846,1.80058,1.81382,1.83367,1.84183,1.84587,1.8423,1.8376,1.8368,1.84023,1.84309,1.85645,1.86316,1.85021,1.84775,1.84463,1.84185};
Int_t fFlag[3][9][5][12];
Double_t fZ_bin_width[12] = {.05,.05,.05,.05,.05,.05,.05,.05,.05,.05,.05,.1};
hadiden fRstudy[24];
studyxy fRstudy_xy[24];
studyxy fR_xy[24];

//Draw

TH1F* fKinematics[6];
TH2F* fKinematics2D;
TH2F* fHO03;
TH2F* fHO04;
TH2F* fTarget2D;
TH2F* fRICHLH;
TCanvas c1("Kin_Q^{2}","Kin_Q^{2}",3200,1600);
TCanvas c2("Kin_x^{Bj}","Kin_x^{Bj}",3200,1600);
TCanvas c3("Kin_y","Kin_y",3200,1600);
TCanvas c4("Kin_z","Kin_z",3200,1600);
TCanvas c5("Kin_w","Kin_w",3200,1600);
TCanvas c6("Kin_#nu","Kin_#nu",3200,1600);
TCanvas c7("Kin_xy","Kin_xy",3200,1600);
TCanvas c8("target_xy","target_xy",3200,1600);
TCanvas c9("RICH LH","RICH LH",3200,1600);
TCanvas c10("HO03","HO03",3200,1600);
TCanvas c11("HO04","HO04",3200,1600);

//Graphic Style

Int_t fMarkerColor[2] = {4,2};
Int_t fMarkerStyle[2] = {24,20};

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


//Constants

static const Double_t fM_p = 938.272046/(1e3);
static const Double_t fM_mu = 105.6583715/(1e3);
static const Double_t fM_K = 493.677/(1e3);
static const Double_t fM_pi = 139.57018/(1e3);

#endif

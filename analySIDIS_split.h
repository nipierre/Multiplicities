#ifndef __CUTMAKER_SPLIT_H_INCLUDED__
#define __CUTMAKER_SPLIT_H_INCLUDED__

#include <iostream>
#include <iomanip>
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
#include <math.h>
#include <TColor.h>
#include <TROOT.h>
#include <TStyle.h>

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

struct Wrapper { Float_t tab[2][2][2][4]; };
struct Multiplicities { Float_t tab[2][3][4]; };
struct Recovery { vector<Float_t> vec[2][4][4]; };
struct Recovery_tank { Float_t tab[2][4][4]; };
struct Pvsz { vector<Float_t> vec[2][6]; };
struct hadiden { vector<Float_t> vec; };
struct studyxy { vector<Float_t> vec[2]; };

//DISevt
Int_t            fRunNo;
vector<Int_t>    fSpillNo;
vector<Int_t>    fEvtInSpill;
vector<bool>     fTrigMask;
vector<Int_t>    fEvNo;
vector<Float_t> fX;
vector<Float_t> fY;
vector<Float_t> fZ;
vector<Float_t> fP0x;
vector<Float_t> fP0y;
vector<Float_t> fP0z;
vector<Float_t> fP1x;
vector<Float_t> fP1y;
vector<Float_t> fP1z;
vector<Float_t> fE_beam;
vector<Float_t> fE_mu_prim;
vector<Float_t> fXX0;
vector<Float_t> fHM04x;
vector<Float_t> fHM04y;
vector<Float_t> fHM05x;
vector<Float_t> fHM05y;
vector<Float_t> fHO03x;
vector<Float_t> fHO03y;
vector<Float_t> fHO04x;
vector<Float_t> fHO04y;
vector<bool>     fSaved;
vector<bool>     fCellCrossed;

//Hadrons
vector<Float_t> fP;
vector<Float_t> fTh;
vector<Float_t> fPh;
vector<Float_t> fHXX0;
vector<bool>     fInHCALacc;
vector<Float_t> fHCAL;
vector<Int_t>    fCharge;
vector<Float_t> fThRICH;
vector<Float_t> fLH[6];
vector<Int_t>    fMCpid;
vector<Float_t> fMM01x;
vector<Float_t> fMM01y;
vector<Float_t> fMM02x;
vector<Float_t> fMM02y;
vector<Float_t> fMM03x;
vector<Float_t> fMM03y;
vector<Float_t> fZ2Ax;
vector<Float_t> fZ2Ay;
vector<Float_t> fZ2Bx;
vector<Float_t> fZ2By;
vector<Float_t> fRICHx;
vector<Float_t> fRICHy;
vector<Float_t> fAtsize;

vector<Float_t> fQ2;
vector<Float_t> fQ2kin;
vector<Float_t> fQ2local;
vector<Float_t> fXBj;
vector<Float_t> fXBjkin;
vector<Float_t> fYBj;
vector<Float_t> fYBjkin;
vector<Float_t> fWBj;
vector<Float_t> fWBjkin;
vector<Float_t> fNu;
vector<Float_t> fNukin;
vector<Float_t> fLHpi;
vector<Float_t> fLHK;
vector<vector<Float_t>> fTheta;
vector<vector<Float_t>> fpT;

vector<Pvsz> fPvsz;
vector<Pvsz> fPvsz_err;
vector<hadiden> fHadiden;

//Misc

set<Float_t> fLHsec_set;
Float_t* fLHsec_tab;
Float_t fLHsec;
Int_t fId;
Int_t fId_loose;
Int_t fId_severe;
Float_t fNu_max[3][12];
Float_t fNu_min[3][12];
vector<Float_t> fXv;
vector<Float_t> fYv;
vector<Float_t> fZv;
vector<Float_t> fRv;
Int_t fMuCharge;

//Counting

Double_t fBP = 0;
Double_t fRmu = 0;
Double_t fBMS = 0;
Double_t fMuchi2 = 0;
Double_t fBEC = 0;
Double_t fVtx = 0;
Double_t fTarg = 0;
Double_t fCell = 0;
Double_t fTrig = 0;
Double_t fMupchi2 = 0;
Double_t fMZfirst = 0;
Double_t fQ2test = 0;
Double_t fYBjtest = 0;
Double_t fXBjtest = 0;
Double_t fWBjtest = 0;
Double_t fHadrons = 0;
Double_t fXX0test = 0;
Double_t fChi2Hadron = 0;
Double_t fHZfirst = 0;
Double_t fHZlast = 0;
Double_t fMom = 0;
Double_t fTRICH = 0;
Double_t fPosRICH = 0;
Double_t fZtest = 0;
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
Double_t fFilesNumber = 0;
vector<Int_t> fPeriodBit;
vector<string> fPeriodName;

Float_t fInclusiveRCproton[30][19];
Float_t fSemiInclusiveRCproton[9][6][14];

Float_t fCepi[2][2][9][6][12];
Float_t fCepiVtx[2][2][9][6][12][4];
Float_t fCepiTh[2][2][9][6][12][8];
Float_t fCepipT[2][2][9][6][12][10];


//Binning
Wrapper fBinning[9][6][12];
Wrapper fBinning_zvtx[9][6][12][4];
Wrapper fBinning_theta[9][6][12][8];
Wrapper fBinning_pt[9][6][12][10];
Wrapper fBinning_loose[9][6][12];
Wrapper fBinning_severe[9][6][12];
Multiplicities fMultiplicities[9][6][12];
Wrapper fAcceptance[9][6][12];
Wrapper fRich_sys_err[9][6][12]; // tab[][0][] : stat, tab[][1][] : sys
Recovery fMeanvalues[9][6][12]; // tab[][][i], iC[0,3] : x,y,Q2,z
Recovery_tank fMeanvalues_size[9][6][12];
Recovery_tank fMeanvalues_data[9][6][12];
Float_t fNDIS_evt[3][2][9][6][12];
Float_t fNDIS_evt_err[3][2][9][6][12];
Float_t fNDIS_evt_zvtx[3][2][9][6][12][4];
Float_t fNDIS_evt_err_zvtx[3][2][9][6][12][4];
Int_t xbin, ybin, zbin, zlabbin, thbin, ptbin;
Float_t fZrange[13] = {.20,.25,.30,.35,.40,.45,.50,.55,.60,.65,.70,.75,.85};
Float_t fXrange[10] = {.004,.01,.02,.03,.04,.06,.1,.14,.18,.4};
Float_t fYrange[7] = {.1,.15,.2,.3,.5,.7,.9};
Float_t fRcutval[24] = {1.68812,1.6915,1.69572,1.69733,1.71178,1.74735,1.74682,1.7846,1.80058,1.81382,1.83367,1.84183,1.84587,1.8423,1.8376,1.8368,1.84023,1.84309,1.85645,1.86316,1.85021,1.84775,1.84463,1.84185};
Int_t fFlag[3][9][5][12];
Float_t fZ_bin_width[12] = {.05,.05,.05,.05,.05,.05,.05,.05,.05,.05,.05,.1};
hadiden fRstudy[24];
studyxy fRstudy_xy[24];
studyxy fR_xy[24];

//Draw

TH1F* fKinematics[6];
TH2F* fKinematics2D[2];
TH2F* fKinematicsRICH;
TH2F* fHO03;
TH2F* fHO04;
TH2F* fTarget2D;
TH2F* fRICHLH;
TH2F* fAllTarget[12];
TH2F* fInTarget[12];
TH1F* fZvtx[2];
TH1F* fQ2k[2];
TH1F* fYk[2];
TH1F* fThRich[2];
TH1F* fPk[2];
TH1F* fZk[2];
TCanvas c1("Kin_Q^{2}","Kin_Q^{2}");
TCanvas c2("Kin_x^{Bj}","Kin_x^{Bj}");
TCanvas c3("Kin_y","Kin_y");
TCanvas c4("Kin_z","Kin_z");
TCanvas c5("Kin_w","Kin_w");
TCanvas c6("Kin_#nu","Kin_#nu");
TCanvas c7("Kin_xy","Kin_xy");
TCanvas c71("Kin_xQ2","Kin_xQ2");
TCanvas c8("target_xy","target_xy");
TCanvas c9("RICH LH","RICH LH");
TCanvas c10("HO03","HO03");
TCanvas c11("HO04","HO04");
TCanvas c12("RICH_spec","RICH_spec");
TCanvas c13("Target_cut","Target_cut");
TCanvas c14("Target_cut2","Target_cut2");
TCanvas c15("Cuts","Cuts");

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
Float_t pi_unfolding_err_p[2][10][3];
Float_t pi_unfolding_err_m[2][10][3];
Float_t k_unfolding_err_p[2][10][3];
Float_t k_unfolding_err_m[2][10][3];
Float_t p_unfolding_err_p[2][10][3];
Float_t p_unfolding_err_m[2][10][3];
Int_t mat_bin[2][10];
Int_t err_bin[2][10];
Float_t cov1_pi[2][3];
Float_t cov1_k[2][3];
Float_t cov1_p[2][3];
Float_t cov2[2][3];
TMatrixD pi_sigma_uni(3,3);
TMatrixD k_sigma_uni(3,3);
TMatrixD p_sigma_uni(3,3);
TMatrixD pi_vect(3,1);
TMatrixD k_vect(3,1);
TMatrixD p_vect(3,1);
Float_t dummy;


//Constants

static const Float_t fM_p = 938.272046/(1e3);
static const Float_t fM_mu = 105.6583715/(1e3);
static const Float_t fM_K = 493.677/(1e3);
static const Float_t fM_pi = 139.57018/(1e3);

#endif

#ifndef __DVM_H_INCLUDED__
#define __DVM_H_INCLUDED__

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

//Structs

struct Wrapper { Double_t tab[2][2][4]; };
struct Pvsz { vector<Double_t> vec[2][5]; };
struct hadiden { vector<Double_t> vec; };
struct studyxy { vector<Double_t> vec[2]; };
struct Recovery { vector<Double_t> vec[2][4][4]; };
struct Recovery_tank { Double_t tab[2][4][4]; };

//Misc

int fId;
Double_t fNu_max[3][12];
Double_t fNu_min[3][12];
vector<Double_t> fXv;
vector<Double_t> fYv;
vector<Double_t> fZv;
vector<Double_t> fRv;

bool fAllDISflag;
bool fAllDISflag_MC;

//Weighted events
Double_t SIDIS_WEIGHT;
Double_t SIDIS_EVENTS;
Double_t RHO_WEIGHT;
Double_t RHO_EVENTS;
Double_t PHI_WEIGHT;
Double_t PHI_EVENTS;


//Binning
Wrapper fSIDIS[9][6][12];
Double_t fSIDIS_tot[2][4];
Wrapper fRho[9][6][12];
Double_t fRho_tot[2][4];
Wrapper fPhi[9][6][12];
Double_t fPhi_tot[2][4];
Wrapper fDVM_h[9][6][12];
Double_t fNDIS_evt_SIDIS[9][6];
Double_t fNDIS_SIDIS_tot;
Double_t fNDIS_evt_rho[9][6];
Double_t fNDIS_rho_tot;
Double_t fNDIS_evt_phi[9][6];
Double_t fNDIS_phi_tot;
Double_t fDVM_DIS_pi[9][6];
Double_t fDVM_DIS_K[9][6];
Double_t fDVM_pi[2][9][6][12];
Double_t fDVM_K[2][9][6][12];
Double_t fDVM_pi_err[2][9][6][12];
Double_t fDVM_K_err[2][9][6][12];
int xbin, ybin, zbin, xbin_MC, ybin_MC, zbin_u;
Double_t fZrange[13] = {.20,.25,.30,.35,.40,.45,.50,.55,.60,.65,.70,.75,.85};
Double_t fXrange[10] = {.004,.01,.02,.03,.04,.06,.1,.14,.18,.4};
Double_t fYrange[6] = {.1,.15,.2,.3,.5,.7};
int fFlag[3][9][5][12];
Double_t fZ_bin_width[12] = {.05,.05,.05,.05,.05,.05,.05,.05,.05,.05,.05,.1};

Wrapper fBinning[9][5][12];

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

#endif

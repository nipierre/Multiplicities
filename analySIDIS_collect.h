#ifndef __CUTMAKER_COLLECT_H_INCLUDED__
#define __CUTMAKER_COLLECT_H_INCLUDED__

#include <vector>
#include <set>
#include <map>
#include <utility>
#include <TLatex.h>
#include <TLegend.h>

using namespace std;

struct Wrapper { Double_t tab[2][2][2][4]; };
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
vector<Double_t> fQ2local;
vector<Double_t> fXBj;
vector<Double_t> fYBj;
vector<Double_t> fWBj;
vector<Double_t> fNu;

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

int fNumberPeriod;
vector<int> fPeriods;

Double_t fSemiInclusiveRCproton[2][9][6][14];
Double_t fDiffVectorMeson[2][9][6][12][4];
Double_t fQelCorr[9][6];

Double_t PeriodFlux[2][11] = {{0,0,0,0,0,0,8.20114,4.40163,6.43954,6.06487,0},
                              {0,0,0,0,0,0,9.14961,6.67131,6.86646,7.27168,2.08672}};

Double_t PeriodFluxTot;

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
Double_t fChi2Hadron = 0;
Double_t fZfirst = 0;
Double_t fZlast = 0;
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

//Binning
Wrapper fBinning[9][6][12];
Wrapper fBinning_zvtx[9][6][12][4];
Wrapper fBinning_period[11][9][6][12];
Wrapper fBinning_period_zvtx[11][9][6][12][4];
Wrapper fBinning_period_theta[11][9][6][12][8];
Wrapper fBinning_period_pt[11][9][6][12][10];
Wrapper fBinning_loose_period[11][9][6][12];
Wrapper fBinning_severe_period[11][9][6][12];
Wrapper fBinning_yavg[9][12];
Wrapper fBinning_zavg[9];
Wrapper fBinning_yavg_period[11][9][12];
Wrapper fBinning_zavg_period[11][9];
Multiplicities fMultiplicities[9][6][12];
Multiplicities fMultiplicities_zvtx[9][6][12][4];
Multiplicities fMultiplicities_theta[9][6][12][8];
Multiplicities fMultiplicities_thetaint[9][12];
Multiplicities fMultiplicities_pt[9][6][12][10];
Multiplicities fMultiplicities_ptint[9][6][12];
Multiplicities fMultiplicities_yavg[9][12];
Multiplicities fMultiplicities_theta_yavg[9][12][8];
Multiplicities fMultiplicities_pt_yavg[9][12];
Multiplicities fMultiplicities_zvtx_yavg[9][12][4];
Multiplicities fMultiplicities_zavg[9];
Wrapper fDiffractiveVectorMeson[9][6][12];
Wrapper fAcceptance[11][9][6][12];
Wrapper fAcceptance_zvtx[11][9][6][12][4];
Wrapper fAcceptance_theta[11][9][6][12][8];
Wrapper fAcceptance_pt[11][9][6][12][10];
Wrapper fAcceptance_yavg[11][9][12];
Wrapper fAcceptance_yavg_weighted[9][12];
Wrapper fRich_sys_err_period[11][9][6][12]; // tab[][0][] : stat, tab[][1][] : sys
Wrapper fRich_sys_err_yavg[9][12];
Wrapper fRich_sys_err_zavg[9];
Recovery fMeanvalues[9][6][12]; // tab[][][i], iC[0,3] : x,y,Q2,z
Recovery_tank fMeanvalues_size[9][6][12];
Recovery_tank fMeanvalues_size_periods[11][9][6][12];
Recovery_tank fMeanvalues_data[9][6][12];
Recovery_tank fMeanvalues_yavg[9][12];
Recovery_tank fMeanvalues_data_periods[11][9][6][12];
Recovery_tank fMeanvalues_temp[9][6][12];
Double_t fNDIS_evt_period[11][3][2][9][6][12];
Double_t fNDIS_evt_zvtx_period[11][3][2][9][6][12][4];
Double_t fNDIS_evt_yavg_period[11][3][2][9][12];
Double_t fNDIS_evt_zavg[3][2][9];
Double_t fNDIS_evt_err_period[11][3][2][9][6][12];
Double_t fNDIS_evt_err_zvtx_period[11][3][2][9][6][12][4];
Double_t fNDIS_evt_err_yavg_period[11][3][2][9][12];
Int_t xbin, ybin, zbin;
Double_t fZrange[13] = {.20,.25,.30,.35,.40,.45,.50,.55,.60,.65,.70,.75,.85};
Double_t fThrange[9] = {0.,.015,.025,.035,.045,.058,.072,.088,.2};
Double_t fpTrange[11] = {.02,.08,.14,.23,.35,.52,.76,1.12,1.52,2.05,3.};
Double_t fXrange[10] = {.004,.01,.02,.03,.04,.06,.1,.14,.18,.4};
Double_t fYrange[7] = {.1,.15,.2,.3,.5,.7,.9};
Double_t fRcutval[24] = {1.68812,1.6915,1.69572,1.69733,1.71178,1.74735,1.74682,1.7846,1.80058,1.81382,1.83367,1.84183,1.84587,1.8423,1.8376,1.8368,1.84023,1.84309,1.85645,1.86316,1.85021,1.84775,1.84463,1.84185};
Int_t fFlag[3][9][5][12];
Double_t fZ_bin_width[12] = {.05,.05,.05,.05,.05,.05,.05,.05,.05,.05,.05,.1};
Double_t fTh_bin_width[8] = {.015,.01,.01,.01,.013,.014,.016,.112};
Double_t fpT_bin_width[10] = {.06,.06,.09,.12,.17,.24,.36,.40,.53,.95};
hadiden fRstudy[24];
studyxy fRstudy_xy[24];
studyxy fR_xy[24];

//Graphic Style

int fMarkerColor[6] = {2,95,209,226,4,221};
int fMarkerColorZvtx[4][2] = {{kBlue,kRed},{kBlue-7,kRed-7},{kBlue-9,kRed-9},{kBlue-10,kRed-10}};
int fMarkerStyle[6][2] = {{20,20},{26,22},{25,21},{27,33},{28,34},{30,29}};

//Constants

static const Double_t fM_p = 938.272046/(1e3);
static const Double_t fM_mu = 105.6583715/(1e3);
static const Double_t fM_K = 493.677/(1e3);
static const Double_t fM_pi = 139.57018/(1e3);

#endif

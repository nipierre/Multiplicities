#include <iostream>
#include <iomanip>
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

struct Wrapper { Double_t tab[2][2][5]; };
struct Pvsz { vector<Double_t> vec[2][5]; };
struct hadiden { vector<Double_t> vec; };
struct studyxy { vector<Double_t> vec[2]; };


//Qty containers

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
vector<Double_t> fMu;
vector<Double_t> fMu0Cov;
vector<Double_t> fNukin;

vector<Double_t> fQ2_MC;
vector<Double_t> fQ2kinMC;
vector<Double_t> fXBj_MC;
vector<Double_t> fXBjkinMC;
vector<Double_t> fYBj_MC;
vector<Double_t> fYBjkinMC;
vector<Double_t> fWBj_MC;
vector<Double_t> fWBjkinMC;
vector<Double_t> fNu_MC;
vector<Double_t> fMuMC;
vector<Double_t> fNukinMC;

vector<Double_t> fX;
vector<Double_t> fY;
vector<Double_t> fZ;
vector<Double_t> fXMC;
vector<Double_t> fYMC;
vector<Double_t> fZMC;

vector<Double_t> fHM04x;
vector<Double_t> fHM04y;
vector<Double_t> fHM05x;
vector<Double_t> fHM05y;
vector<Double_t> fHL04x;
vector<Double_t> fHL04y;
vector<Double_t> fHL05x;
vector<Double_t> fHL05y;
vector<Double_t> fHO03x;
vector<Double_t> fHO03y;
vector<Double_t> fHO04x;
vector<Double_t> fHO04y;
vector<Double_t> fHG01x;
vector<Double_t> fHG01y;
vector<Double_t> fHG021x;
vector<Double_t> fHG022x;
vector<Double_t> fHG021y;
vector<Double_t> fHG022y;

vector<Double_t> fHM04MCx;
vector<Double_t> fHM04MCy;
vector<Double_t> fHM05MCx;
vector<Double_t> fHM05MCy;
vector<Double_t> fHL04MCx;
vector<Double_t> fHL04MCy;
vector<Double_t> fHL05MCx;
vector<Double_t> fHL05MCy;
vector<Double_t> fHO03MCx;
vector<Double_t> fHO03MCy;
vector<Double_t> fHO04MCx;
vector<Double_t> fHO04MCy;
vector<Double_t> fHG01MCx;
vector<Double_t> fHG01MCy;
vector<Double_t> fHG021MCx;
vector<Double_t> fHG022MCx;
vector<Double_t> fHG021MCy;
vector<Double_t> fHG022MCy;

vector<Double_t> fTCx;
vector<Double_t> fTCy;

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
Double_t fMCDIS = 0;
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
Wrapper fRcstr[9][6][12];
Wrapper fGnrt[9][6][12];
Wrapper fRcstr_c[9][6][12];
Wrapper fRcstr_zvtx[9][6][12][4];
Wrapper fGnrt_zvtx[9][6][12][4];
Double_t fNDIS_evt[3][9][6][12];
Double_t fNDIS_evt_c[3][9][6][12];
Double_t fNDIS_evt_MC[3][9][6][12];
Double_t fNDIS_evt_zvtx[3][9][6][12][4];
Double_t fNDIS_evt_MC_zvtx[3][9][6][12][4];
int xbin, ybin, zbin, xbin_MC, ybin_MC, zbin_u;
Double_t fZrange[13] = {.20,.25,.30,.35,.40,.45,.50,.55,.60,.65,.70,.75,.85};
Double_t fXrange[10] = {.004,.01,.02,.03,.04,.06,.1,.14,.18,.4};
Double_t fYrange[7] = {.1,.15,.2,.3,.5,.7,.9};
Double_t fRcutval[24] = {1.68812,1.6915,1.69572,1.69733,1.71178,1.74735,1.74682,1.7846,1.80058,1.81382,1.83367,1.84183,1.84587,1.8423,1.8376,1.8368,1.84023,1.84309,1.85645,1.86316,1.85021,1.84775,1.84463,1.84185};
int fFlag[3][9][6][12];
int fFlag_MC[3][9][6][12];
Double_t fZ_bin_width[12] = {.05,.05,.05,.05,.05,.05,.05,.05,.05,.05,.05,.1};

//Draw

TH1F* fKinematics[7];
TH2F* fKinematics2D;
TH2F* fBeamCovariance;
TH2F* fHM04;
TH2F* fHM05;
TH2F* fHL04;
TH1F* fHL04x1D;
TH1F* fHL04MCx1D;
TH2F* fHL05;
TH1F* fHL05x1D;
TH1F* fHL05MCx1D;
TH2F* fHLx2D;
TH2F* fHO03;
TH2F* fHO04;
TH2F* fHG01;
TH2F* fHG021;
TH2F* fHG022;
TH2F* fHM04MC;
TH2F* fHM05MC;
TH2F* fHL04MC;
TH2F* fHL04MCx2D;
TH2F* fHL05MC;
TH2F* fHL05MCx2D;
TH2F* fHO03MC;
TH2F* fHO04MC;
TH2F* fHG01MC;
TH2F* fHG021MC;
TH2F* fHG022MC;
TH2F* fTarget2D;
TH2F* fTrigCov;
TH1F* fVertexHadron[3];
TH1F* fVertexHadronMC[3];
TCanvas c1("Kin_Q^{2}","Kin_Q^{2}",3200,1600);
TCanvas c2("Kin_x^{Bj}","Kin_x^{Bj}",3200,1600);
TCanvas c3("Kin_y","Kin_y",3200,1600);
TCanvas c4("Kin_z","Kin_z",3200,1600);
TCanvas c5("Kin_w","Kin_w",3200,1600);
TCanvas c6("Kin_#nu","Kin_#nu",3200,1600);
TCanvas c7("Kin_xy","Kin_xy",3200,1600);
TCanvas c8("target_xy","target_xy",3200,1600);
TCanvas c9("HM04X1","HM04X1",3200,1600);
TCanvas c10("HM05X1","HM05X1",3200,1600);
TCanvas c11("HL04Y1","HL04Y1",3200,1600);
TCanvas c12("HL05Y1","HL05Y1",3200,1600);
TCanvas c13("HO03Y1","HO03Y1",3200,1600);
TCanvas c14("HO04Y1","HO04Y1",3200,1600);
TCanvas c15("HG01Y1","HG01Y1",3200,1600);
TCanvas c16("HG02Y1","HG02Y1",3200,1600);
TCanvas c17("HG02Y2","HG02Y2",3200,1600);
TCanvas c26("HL04Y1y","HL04Y1y",3200,1600);
TCanvas c27("HL05Y1y","HL05Y1y",3200,1600);
TCanvas c28("HLx","HLx",3200,1600);
TCanvas c29("HM04MCX1","HM04MCX1",3200,1600);
TCanvas c30("HM05MCX1","HM05MCX1",3200,1600);
TCanvas c31("HL04MCY1","HL04MCY1",3200,1600);
TCanvas c32("HL05MCY1","HL05MCY1",3200,1600);
TCanvas c33("HO03MCY1","HO03MCY1",3200,1600);
TCanvas c34("HO04MCY1","HO04MCY1",3200,1600);
TCanvas c35("HG01MCY1","HG01MCY1",3200,1600);
TCanvas c36("HG02MCY1","HG02MCY1",3200,1600);
TCanvas c37("HG02MCY2","HG02MCY2",3200,1600);
TCanvas c38("HL04MCY1x","HL04MCY1x",3200,1600);
TCanvas c39("HL05MCY1x","HL05MCY1x",3200,1600);
TCanvas c40("E_{#mu}","E_{#mu}",3200,1600);
TCanvas c41("Trigger_Coverage","Trigger_Coverage",3200,1600);
TCanvas c42("Vertex Hadron+e","Vertex Hadron+e",3200,1600);
TCanvas c43("Vertex MC","Vertex MC",1200,1200);
TCanvas c44("Beam Covariance","Beam Covariance",1200,1200);

TH1F* fKinematicsMC[7];
TH2F* fKinematics2DMC;
TH2F* fTarget2DMC;
TH1F* fVertexStudyMC[4];
TH2F* fVertexStudyMC2D[3];
TCanvas c18("KinMC_Q^{2}","KinMC_Q^{2}",3200,1600);
TCanvas c19("KinMC_x^{Bj}","KinMC_x^{Bj}",3200,1600);
TCanvas c20("KinMC_y","KinMC_y",3200,1600);
TCanvas c21("KinMC_z","KinMC_z",3200,1600);
TCanvas c22("KinMC_w","KinMC_w",3200,1600);
TCanvas c23("KinMC_#nu","KinMC_#nu",3200,1600);
TCanvas c24("KinMC_xy","KinMC_xy",3200,1600);
TCanvas c25("targetMC_xy","targetMC_xy",3200,1600);

//Graphic Style

Int_t fMarkerColor[2] = {4,2};
Int_t fMarkerStyle[2] = {24,20};

//Constants

static const Double_t fM_p = 938.272046/(1e3);
static const Double_t fM_mu = 105.6583715/(1e3);
static const Double_t fM_K = 493.677/(1e3);
static const Double_t fM_pi = 139.57018/(1e3);

// Reweighting
Double_t Cth[2][31] = {{1.70, 1.45, 1.32, 1.25, 1.23, 1.20, 1.12, 1.01, 0.99, 0.96, 0.97, 0.88, 0.87, 0.83, 0.82, 0.73, 0.70, 0.69, 0.69, 0.65, 0.68, 0.62, 0.61, 0.64, 0.62, 0.65, 0.60, 0.58, 0.70, 0.68, 0.69},
                       {1.66, 1.40, 1.32, 1.32, 1.31, 1.23, 1.15, 1.08, 1.05, 1.03, 0.98, 0.95, 0.89, 0.85, 0.80, 0.77, 0.72, 0.71, 0.64, 0.62, 0.64, 0.55, 0.55, 0.56, 0.50, 0.56, 0.58, 0.51, 0.51, 0.54, 0.53}};

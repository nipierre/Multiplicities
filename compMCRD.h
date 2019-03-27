#include <vector>
#include <string>
#include <set>
#include <map>
#include <utility>
#include <iostream>
#include <iomanip>
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
#include <TRatioPlot.h>
#include <TLegend.h>

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

Double_t fNEventsRD = 0;
Double_t fNEventsMC = 0;
Double_t fHadronRD = 0;
Double_t fHadronMC = 0;
Double_t fHadronpRD = 0;
Double_t fHadronpMC = 0;
Double_t fHadronmRD = 0;
Double_t fHadronmMC = 0;

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
TH2F* fECAL0RD;
TH1F* fVertexRD[5];
TH1F* fKinematicsMC[5][20];
TH2F* fThetaMCp[3];
TH2F* fECAL0MC;
TH1F* fVertexMCb[5];
TRatioPlot* fKinematicsRatio[5][20];
Int_t fCountingRD[5][20];
Int_t fCountingMC[5][20];
TH1F* fTarget[2];
TH1F* fTargetMC[2];
TH2F* fTarget2D[2];
TH2F* fTarget2DMC[2];
string unit[20] = {"Q^{2}","x","y","z","W","#nu","#mu","#mu'","#theta_{#mu}","#Phi_{#mu}",
                   "z_{vtx}","#Phi_{e,prod.pl}","p_{h}","#theta_{h}","#Phi_{h}","#Phi_{h,prod.pl}","p_{T}","z_{vtx,h}","z_{vtx,e}","z_{vtx,h}"};
TCanvas c1("Kin_{Q2} Trigger","Kin_{Q2} Trigger");
TCanvas c2("Kin_x^{Bj} Trigger","Kin_x^{Bj} Trigger");
TCanvas c3("Kin_y Trigger","Kin_y Trigger");
TCanvas c4("Kin_z Trigger","Kin_z Trigger");
TCanvas c5("Kin_w Trigger","Kin_w Trigger");
TCanvas c6("Kin_nu Trigger","Kin_nu Trigger");
TCanvas c7("Kin_Phi Trigger","Kin_Phi Trigger");
TCanvas c8("Kin_Q^{2}","Kin_Q^{2}");
TCanvas c9("Kin_x^{Bj}","Kin_x^{Bj}");
TCanvas c10("Kin_y","Kin_y");
TCanvas c11("Kin_z","Kin_z");
TCanvas c12("Kin_w","Kin_w");
TCanvas c13("Kin_nu","Kin_nu");
TCanvas c14("Kin_mu","Kin_mu");
TCanvas c15("Kin_mup","Kin_mup");
TCanvas c16("Kin_theta","Kin_theta");
TCanvas c17("Kin_phi","Kin_phi");
TCanvas c18("Kin_vertex","Kin_vertex");
TCanvas c19("Kin_hadron_p Trigger","Kin_hadron_p Trigger");
TCanvas c20("Kin_hadron_theta Trigger","Kin_hadron_theta Trigger");
TCanvas c21("Kin_hadron_phi Trigger","Kin_hadron_phi Trigger");
TCanvas c22("Kin_hadron_p","Kin_hadron_p");
TCanvas c23("Kin_hadron_theta","Kin_hadron_theta");
TCanvas c24("Kin_hadron_phi","Kin_hadron_phi");
TCanvas c25("Kin_hadron_phi_pl Trigger","Kin_hadron_phi_pl Trigger");
TCanvas c26("Kin_hadron_phi_pl","Kin_hadron_phi_pl");
TCanvas c27("Kin_p_T Trigger","Kin_hadron_p_T Trigger");
TCanvas c28("Kin_p_T","Kin_hadron_p_T");
TCanvas c29("photoelec","photoelec",2000,2000);
TCanvas c30("Q^2","Q^2",2000,2000);
TCanvas c31("nu","nu",2000,2000);
TCanvas c32("z","z",2000,2000);
TCanvas c33("E_mu","E_mu");
TCanvas c34("Thetay_mu","Thetay_mu");
TCanvas c35("Thetax_mu","Thetax_mu");
TCanvas c36("Thetaxy_mu","Thetaxy_mu");
TCanvas c37("Vertex (Hadron+e) Trig","Vertex (Hadron+e) Trig");
TCanvas c38("Vertex (Hadron+e)","Vertex (Hadron+e)");
TCanvas c39("Vertex (e) Trig","Vertex (e) Trig");
TCanvas c40("Vertex (e)","Vertex (e)");
TCanvas c41("Vertex (Hadron) Trig","Vertex (Hadron) Trig");
TCanvas c42("Vertex (Hadron)","Vertex (Hadron)");
TCanvas c43("ECAL0","ECAL0");
TCanvas c44("Vertex Endpoint","Vertex Endpoint");
TCanvas c45("Target","Target");

vector<double> fError, fErrorRD, fErrorMC;
Int_t fLineStyle[7] = {3,3,3,1,3,3,3};

TLine* l1[20][7];

//Graphic Style

Int_t fMarkerColor[2] = {4,2};
Int_t fMarkerStyle[2] = {24,20};

//Constants

static const Double_t fM_p = 938.272046/(1e3);
static const Double_t fM_mu = 105.6583715/(1e3);
static const Double_t fM_K = 493.677/(1e3);
static const Double_t fM_pi = 139.57018/(1e3);

string trigname[5] = {"MT","LT","OT","LAST",""};

// ROOT tree

TFile *mf;
TTree* DIS;
TTree* Hadron;
TTree* DISMC;
TTree* HadronMC;

Double_t xd, xh, xd_MC, xh_MC;
Double_t yd, yh, yd_MC, yh_MC;
Double_t zh, zh_MC;
Double_t xVTXd, xVTXh, xVTXd_MC, xVTXh_MC;
Double_t yVTXd, yVTXh, yVTXd_MC, yVTXh_MC;
Double_t zVTXd, zVTXh, zVTXd_MC, zVTXh_MC;
Int_t trigd, trigh, trigd_MC, trigh_MC;
Double_t nud, nuh, nud_MC, nuh_MC;
Int_t mu_charged, mu_chargeh, mu_charged_MC, mu_chargeh_MC;
Double_t eVTX_MC;
Double_t Wd, Wh, Wd_MC, Wh_MC;
Double_t thh, th_MC;
Double_t thRICH, thRICH_MC;
Double_t richpipe, richpipe_MC;
Double_t phh, ph_MC;
Double_t phpl, phpl_MC;
Double_t phad, thChad;
Double_t phad_MC, thChad_MC;
Double_t EECAL1, EECAL2, EHCAL1, EHCAL2;
Double_t EECAL1_MC, EECAL2_MC, EHCAL1_MC, EHCAL2_MC;
Int_t isinECAL, isinHCAL;
Int_t isinECAL_MC, isinHCAL_MC;
Int_t PID, PID_MC;
Double_t XX0h, XX0h_MC;
vector<Double_t> xdv, xhv, xd_MCv, xh_MCv;
vector<Double_t> ydv, yhv, yd_MCv, yh_MCv;
vector<Double_t> zhv, zh_MCv;
vector<Double_t> xVTXdv, xVTXhv, xVTXd_MCv, xVTXh_MCv;
vector<Double_t> yVTXdv, yVTXhv, yVTXd_MCv, yVTXh_MCv;
vector<Double_t> zVTXdv, zVTXhv, zVTXd_MCv, zVTXh_MCv;
vector<Int_t> trigdv, trighv, trigd_MCv, trigh_MCv;
vector<Double_t> nudv, nuhv, nud_MCv, nuh_MCv;
vector<Int_t> mu_chargedv, mu_chargehv, mu_charged_MCv, mu_chargeh_MCv;
vector<Double_t> eVTX_MCv;
vector<Double_t> Wdv, Whv, Wd_MCv, Wh_MCv;
vector<Double_t> thhv, th_MCv;
vector<Double_t> thRICHv, thRICH_MCv;
vector<Double_t> richpipev, richpipe_MCv;
vector<Double_t> phhv, ph_MCv;
vector<Double_t> phplv, phpl_MCv;
vector<Int_t> PIDv, PID_MCv;
vector<Double_t> XX0hv, XX0h_MCv;
vector<Double_t> pv, thCv;
vector<Double_t> p_MCv, thC_MCv;
vector<Double_t> EECAL1v, EECAL2v, EHCAL1v, EHCAL2v;
vector<Double_t> EECAL1_MCv, EECAL2_MCv, EHCAL1_MCv, EHCAL2_MCv;
vector<Int_t> isinECALv, isinHCALv;
vector<Int_t> isinECAL_MCv, isinHCAL_MCv;

Double_t Cth[2][2][31] = {0};

Double_t Cthr[2][31] = {{1.70, 1.45, 1.32, 1.25, 1.23, 1.20, 1.12, 1.01, 0.99, 0.96, 0.97, 0.88, 0.87, 0.83, 0.82, 0.73, 0.70, 0.69, 0.69, 0.65, 0.68, 0.62, 0.61, 0.64, 0.62, 0.65, 0.60, 0.58, 0.70, 0.68, 0.69},
                        {1.66, 1.40, 1.32, 1.32, 1.31, 1.23, 1.15, 1.08, 1.05, 1.03, 0.98, 0.95, 0.89, 0.85, 0.80, 0.77, 0.72, 0.71, 0.64, 0.62, 0.64, 0.55, 0.55, 0.56, 0.50, 0.56, 0.58, 0.51, 0.51, 0.54, 0.53}};

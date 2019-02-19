#include <vector>
#include <string>
#include <set>
#include <map>
#include <utility>
#include <TLatex.h>

using namespace std;

//Structs

struct Wrapper { double tab[2][2][2][4]; };
struct Pvsz { vector<double> vec[2][5]; };
struct hadiden { vector<double> vec; };
struct studyxy { vector<double> vec[2]; };


//Qty containers

vector<double> fQ2;
vector<double> fQ2local;
vector<double> fXBj;
vector<double> fYBj;
vector<double> fWBj;
vector<double> fNu;

vector<double> fQ2_MC;
vector<double> fXBj_MC;
vector<double> fYBj_MC;
vector<double> fWBj_MC;
vector<double> fNu_MC;

//Misc

set<double> fLHsec_set;
double* fLHsec_tab;
double fLHsec;
int fId;
int fId_loose;
int fId_severe;
double fNu_max[3][12];
double fNu_min[3][12];
vector<double> fXv;
vector<double> fYv;
vector<double> fZv;

bool fAllDISflag;
bool fAllDISflag_MC;

//Counting

double fBP = 0;
double fRmu = 0;
double fBMS = 0;
double fBEC = 0;
double fTarg = 0;
double fCell = 0;
double fTrig = 0;
double fQ2test = 0;
double fYBjtest = 0;
double fXBjtest = 0;
double fWBjtest = 0;
double fXX0test = 0;
double fMom = 0;
double fTRICH = 0;
double fPosRICH = 0;
double fHplus = 0;
double fHminus = 0;
double fPiplus = 0;
double fPiminus = 0;
double fKplus = 0;
double fKminus = 0;
double fPplus = 0;
double fPminus = 0;
double fMCHplus = 0;
double fMCHminus = 0;
double fMCPiplus = 0;
double fMCPiminus = 0;
double fMCKplus = 0;
double fMCKminus = 0;
double fMCPplus = 0;
double fMCPminus = 0;
double fPiplus_true = 0;
double fPiminus_true = 0;
double fKplus_true = 0;
double fKminus_true = 0;
double fPplus_true = 0;
double fPminus_true = 0;
double fPiplus_err = 0;
double fPiminus_err = 0;
double fKplus_err = 0;
double fKminus_err = 0;
double fPplus_err = 0;
double fPminus_err = 0;

//Binning
Wrapper fRcstr[9][6][12];
Wrapper fGnrt[9][6][12];
Wrapper fRcstr_c[9][6][12];
Double_t fRcstr_yavg[2][4];
Double_t fGnrt_yavg[2][4];
Double_t fRcstr_c_yavg[2][4];
Wrapper fRcstr_zvtx[9][6][12][4];
Wrapper fGnrt_zvtx[9][6][12][4];
Wrapper fAcceptance[9][6][12];
Wrapper fAcceptance_yavg[9][12];
Wrapper fAcceptance_zvtx[9][6][12][4];
Double_t fNDIS_evt[3][2][9][6][12];
Double_t fNDIS_evt_c[3][2][9][6][12];
Double_t fNDIS_evt_MC[3][2][9][6][12];
Double_t fNDIS_evt_yavg[2][4];
Double_t fNDIS_evt_c_yavg[2][4];
Double_t fNDIS_evt_MC_yavg[2][4];
Double_t fNDIS_evt_zvtx[3][2][9][6][12][4];
Double_t fNDIS_evt_MC_zvtx[3][2][9][6][12][4];
int xbin, ybin, zbin, xbin_MC, ybin_MC;
Double_t fZrange[13] = {.20,.25,.30,.35,.40,.45,.50,.55,.60,.65,.70,.75,.85};
Double_t fXrange[10] = {.004,.01,.02,.03,.04,.06,.1,.14,.18,.4};
Double_t fYrange[7] = {.1,.15,.2,.3,.5,.7,.9};
Double_t fRcutval[24] = {1.68812,1.6915,1.69572,1.69733,1.71178,1.74735,1.74682,1.7846,1.80058,1.81382,1.83367,1.84183,1.84587,1.8423,1.8376,1.8368,1.84023,1.84309,1.85645,1.86316,1.85021,1.84775,1.84463,1.84185};
int fFlag[3][9][6][12];
int fFlag_MC[3][9][6][12];
Double_t fZ_bin_width[12] = {.05,.05,.05,.05,.05,.05,.05,.05,.05,.05,.05,.1};

//Graphic Style

int fMarkerColor[6] = {2,95,209,226,4,221};
int fMarkerStyle[6][2] = {{24,20},{26,22},{25,21},{27,33},{28,34},{30,29}};
int fMarkerStyleb[2][2] = {{24,20},{27,33}};

//Constants

static const double fM_p = 938.272046/(1e3);
static const double fM_mu = 105.6583715/(1e3);
static const double fM_K = 493.677/(1e3);
static const double fM_pi = 139.57018/(1e3);

#include <vector>
#include <string>
#include <set>
#include <map>
#include <utility>

using namespace std;

//Structs

struct Wrapper { double tab[2][2][4]; };
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
Wrapper fRcstr[9][5][12];
Wrapper fGnrt[9][5][12];
Wrapper fRcstr_c[9][5][12];
Wrapper fAcceptance[9][5][12];
double fNDIS_evt[3][9][5][12];
double fNDIS_evt_c[3][9][5][12];
double fNDIS_evt_MC[3][9][5][12];
int xbin, ybin, zbin, xbin_MC, ybin_MC;
double fZrange[13] = {.20,.25,.30,.35,.40,.45,.50,.55,.60,.65,.70,.75,.85};
double fXrange[10] = {.004,.01,.02,.03,.04,.06,.1,.14,.18,.4};
double fYrange[6] = {.1,.15,.2,.3,.5,.7};
double fRcutval[24] = {1.68812,1.6915,1.69572,1.69733,1.71178,1.74735,1.74682,1.7846,1.80058,1.81382,1.83367,1.84183,1.84587,1.8423,1.8376,1.8368,1.84023,1.84309,1.85645,1.86316,1.85021,1.84775,1.84463,1.84185};
int fFlag[3][9][5][12];
int fFlag_MC[3][9][5][12];
double fZ_bin_width[12] = {.05,.05,.05,.05,.05,.05,.05,.05,.05,.05,.05,.1};

TLine* l1[5];

//Graphic Style

int fMarkerColor[5] = {2,95,209,226,221};
int fMarkerStyle[5][2] = {{24,20},{26,22},{25,21},{27,33},{28,34}};

//Constants

static const double fM_p = 938.272046/(1e3);
static const double fM_mu = 105.6583715/(1e3);
static const double fM_K = 493.677/(1e3);
static const double fM_pi = 139.57018/(1e3);

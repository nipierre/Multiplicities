#include <vector>
#include <string>
#include <set>
#include <map>
#include <utility>
#include <TLatex.h>

using namespace std;

//Structs

struct Wrapper { double tab[2][2][2][5]; };

//Binning
Wrapper fRcstr[9][6][12];
Wrapper fGnrt[9][6][12];
Wrapper fRcstr_c[9][6][12];
Double_t fRcstr_yavg[4];
Double_t fGnrt_yavg[4];
Double_t fRcstr_c_yavg[4];
Wrapper fRcstr_zvtx[9][6][12][4];
Wrapper fGnrt_zvtx[9][6][12][4];
Wrapper fAcceptance[9][6][12];
Wrapper fAcceptance_yavg[9][12];
Wrapper fAcceptance_zvtx[9][6][12][4];
Double_t fNDIS_evt[3][2][9][6][12];
Double_t fNDIS_evt_c[3][2][9][6][12];
Double_t fNDIS_evt_MC[3][2][9][6][12];
Double_t fNDIS_evt_yavg[4];
Double_t fNDIS_evt_c_yavg[4];
Double_t fNDIS_evt_MC_yavg[4];
Double_t fNDIS_evt_zvtx[3][2][9][6][12][4];
Double_t fNDIS_evt_MC_zvtx[3][2][9][6][12][4];

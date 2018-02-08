#include <iostream>
#include <fstream>
#include <string>
#include <TFile.h>
#include <TGraph.h>
#include <TCanvas.h>

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

vector<double> fZ[9][5];
vector<double> fQ2mean[9][5];

vector<double> fDfav[9][5];
vector<double> fDstr[9][5];
vector<double> fDunf[9][5];
vector<double> fDunf1[9][5];
vector<double> fDunf2[9][5];

TGraph* fDfavG[9][5];
TGraph* fDstrG[9][5];
TGraph* fDunfG[9][5];
TGraph* fDunf1G[9][5];
TGraph* fDunf2G[9][5];

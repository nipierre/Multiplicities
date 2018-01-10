#include <iostream>
#include <fstream>
#include <string>
#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TMatrixD.h>
#include "LHAPDF/GridPDF.h"

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

double fX[9][5][12];
double fY[9][5][12];
double fZ[9][5][12];
double fQ2[9][5][12];

double fKp_p[9][5][12];
double fKm_p[9][5][12];
double fKp_d[9][5][12];
double fKm_d[9][5][12];
double fKpm_d[9][5][12];

double fPip_p[9][5][12];
double fPim_p[9][5][12];
double fPip_d[9][5][12];
double fPim_d[9][5][12];

double fDfav[9][5][12];
double fDstr[9][5][12];
double fDunf[9][5][12];
double fDunf1[9][5][12];
double fDunf2[9][5][12];

string fLHGrid;
TMatrixD fCoeff2E(2,2);
TMatrixD fCoeff(3,3);
TMatrixD fCoeff4E(4,4);
TMatrixD fMult2E(2,1);
TMatrixD fMult(3,1);
TMatrixD fMult4E(4,1);

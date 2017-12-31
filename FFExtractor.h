#include <iostream>
#include <fstream>
#include <string>
#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TMatrixD.h>
#include "LHAPDF/LHAPDF.h"

double fX[9][5][10];
double fY[9][5][10];
double fZ[9][5][10];
double fQ2[9][5][10];

double fKp_p[9][5][10];
double fKm_p[9][5][10];
double fKp_d[9][5][10];
double fKm_d[9][5][10];
double fKpm_d[9][5][10];

double fDfav[9][5][10];
double fDstr[9][5][10];
double fDunf[9][5][10];
double fDunf1[9][5][10];
double fDunf2[9][5][10];

string fLHGrid;
TMatrixD fCoeff(3,3);
TMatrixD fCoeff4E(4,4);
TMatrixD fMult(3,1);
TMatrixD fMult4E(4,1);

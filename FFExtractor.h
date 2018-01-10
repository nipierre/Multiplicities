#include <iostream>
#include <fstream>
#include <string>
#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TMatrixD.h>
#include "LHAPDF/GridPDF.h"

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

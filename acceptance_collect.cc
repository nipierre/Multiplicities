/*
LOG
16/03/16 : apparently working version for [3,40] range
22/03/16 : changing storage form Lyon to CASTOR
0X/04/16 : minor fixes
20/04/16 : added parallelized version
22/04/16 : lifted version
*/

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
#include <TGraphErrors.h>
#include <TArrow.h>

#include "acceptance_collect.h"

// Flags
#define Y2006 0
#define Y2012 0
#define Y2016 1
#define SPREAD 1

// Outputs
#if BW
#define dirroot "/u/sciteam/pierre1/production/software/Multiplicities/acceptance"
#else
#define dirroot "/sps/compass/npierre/Multiplicities/acceptance"
#endif

using namespace std;

Double_t RelDiff(int c, int x, int y, int z, int m, int had)
{
  Double_t min=fAcceptance_zvtx[x][y][z][0].tab[c][m][0][had];
  Double_t max=fAcceptance_zvtx[x][y][z][0].tab[c][m][0][had];

  for(int i=1; i<4; i++)
  {
    if(fAcceptance_zvtx[x][y][z][i].tab[c][m][0][had]>max)
      max=fAcceptance_zvtx[x][y][z][i].tab[c][m][0][had];
    if(fAcceptance_zvtx[x][y][z][i].tab[c][m][0][had]<max)
      min=fAcceptance_zvtx[x][y][z][i].tab[c][m][0][had];
  }

  return (min ? Double_t((max-min)/min) : 0);

}

void yavg(int c, int x, int z)
{
  int rec[2][4] = {{0,0,0,0},{0,0,0,0}};
  int gen[2][4] = {{0,0,0,0},{0,0,0,0}};
  int dis[2][3] = {{0,0,0},{0,0,0}};
  int dis_MC[2][3] = {{0,0,0},{0,0,0}};
  int rec_c[2][4] = {{0,0,0,0},{0,0,0,0}};
  int dis_c[2][3] = {{0,0,0},{0,0,0}};

  for(int i=0; i<4; i++)
  {
    for(int j=0; j<2; j++)
    {
      fRcstr_yavg[j][i]=0;
      fNDIS_evt_yavg[j][i]=0;
      fGnrt_yavg[j][i]=0;
      fNDIS_evt_MC_yavg[j][i]=0;
      fRcstr_c_yavg[j][i]=0;
      fNDIS_evt_c_yavg[j][i]=0;
    }
  }
  for(int i=0; i<6; i++)
  {
    for(int j=0; j<2; j++)
    {
      fRcstr_yavg[j][0]+=fRcstr[x][i][z].tab[c][j][0][0];
      if(fRcstr[x][i][z].tab[c][j][0][0]) rec[j][0]++;
      fRcstr_yavg[j][1]+=fRcstr[x][i][z].tab[c][j][0][1];
      if(fRcstr[x][i][z].tab[c][j][0][1]) rec[j][1]++;
      fRcstr_yavg[j][2]+=fRcstr[x][i][z].tab[c][j][0][2];
      if(fRcstr[x][i][z].tab[c][j][0][2]) rec[j][2]++;
      fRcstr_yavg[j][3]+=fRcstr[x][i][z].tab[c][j][0][3];
      if(fRcstr[x][i][z].tab[c][j][0][3]) rec[j][3]++;
      fNDIS_evt_yavg[j][0]+=fNDIS_evt[0][j][x][i][z];
      if(fNDIS_evt[0][j][x][i][z]) dis[j][0]++;
      fNDIS_evt_yavg[j][1]+=fNDIS_evt[1][j][x][i][z];
      if(fNDIS_evt[1][j][x][i][z]) dis[j][1]++;
      fNDIS_evt_yavg[j][2]+=fNDIS_evt[2][j][x][i][z];
      if(fNDIS_evt[2][j][x][i][z]) dis[j][2]++;
      fGnrt_yavg[j][0]+=fGnrt[x][i][z].tab[c][j][0][0];
      if(fGnrt[x][i][z].tab[c][j][0][0]) gen[j][0]++;
      fGnrt_yavg[j][1]+=fGnrt[x][i][z].tab[c][j][0][1];
      if(fGnrt[x][i][z].tab[c][j][0][1]) gen[j][1]++;
      fGnrt_yavg[j][2]+=fGnrt[x][i][z].tab[c][j][0][2];
      if(fGnrt[x][i][z].tab[c][j][0][2]) gen[j][2]++;
      fGnrt_yavg[j][3]+=fGnrt[x][i][z].tab[c][j][0][3];
      if(fGnrt[x][i][z].tab[c][j][0][3]) gen[j][3]++;
      fNDIS_evt_MC_yavg[j][0]+=fNDIS_evt_MC[0][j][x][i][z];
      if(fNDIS_evt_MC[0][j][x][i][z]) dis_MC[j][0]++;
      fNDIS_evt_MC_yavg[j][1]+=fNDIS_evt_MC[1][j][x][i][z];
      if(fNDIS_evt_MC[1][j][x][i][z]) dis_MC[j][1]++;
      fNDIS_evt_MC_yavg[j][2]+=fNDIS_evt_MC[2][j][x][i][z];
      if(fNDIS_evt_MC[2][j][x][i][z]) dis_MC[j][2]++;
      fRcstr_c_yavg[j][0]+=fRcstr_c[x][i][z].tab[c][j][0][0];
      if(fRcstr_c[x][i][z].tab[c][j][0][0]) rec_c[j][0]++;
      fRcstr_c_yavg[j][1]+=fRcstr_c[x][i][z].tab[c][j][0][1];
      if(fRcstr_c[x][i][z].tab[c][j][0][1]) rec_c[j][1]++;
      fRcstr_c_yavg[j][2]+=fRcstr_c[x][i][z].tab[c][j][0][2];
      if(fRcstr_c[x][i][z].tab[c][j][0][2]) rec_c[j][2]++;
      fRcstr_c_yavg[j][3]+=fRcstr_c[x][i][z].tab[c][j][0][3];
      if(fRcstr_c[x][i][z].tab[c][j][0][3]) rec_c[j][3]++;
      fNDIS_evt_c_yavg[j][0]+=fNDIS_evt_c[0][j][x][i][z];
      if(fNDIS_evt_c[0][j][x][i][z]) dis_c[j][0]++;
      fNDIS_evt_c_yavg[j][1]+=fNDIS_evt_c[1][j][x][i][z];
      if(fNDIS_evt_c[1][j][x][i][z]) dis_c[j][1]++;
      fNDIS_evt_c_yavg[j][2]+=fNDIS_evt_c[2][j][x][i][z];
      if(fNDIS_evt_c[2][j][x][i][z]) dis_c[j][2]++;
    }
  }
  // (rec[0] ? fRcstr_yavg[0] /= rec[0] : 0);
  // (rec[1] ? fRcstr_yavg[1] /= rec[1] : 0);
  // (rec[2] ? fRcstr_yavg[2] /= rec[2] : 0);
  // (rec[3] ? fRcstr_yavg[3] /= rec[3] : 0);
  // (dis[0] ? fNDIS_evt_yavg[0] /= dis[0] : 0);
  // (dis[1] ? fNDIS_evt_yavg[1] /= dis[1] : 0);
  // (dis[2] ? fNDIS_evt_yavg[2] /= dis[2] : 0);
  // (gen[0] ? fGnrt_yavg[0] /= gen[0] : 0);
  // (gen[1] ? fGnrt_yavg[1] /= gen[1] : 0);
  // (gen[2] ? fGnrt_yavg[2] /= gen[2] : 0);
  // (gen[3] ? fGnrt_yavg[3] /= gen[3] : 0);
  // (dis_MC[0] ? fNDIS_evt_MC_yavg[0] /= dis_MC[0] : 0);
  // (dis_MC[1] ? fNDIS_evt_MC_yavg[1] /= dis_MC[1] : 0);
  // (dis_MC[2] ? fNDIS_evt_MC_yavg[2] /= dis_MC[2] : 0);
  // (rec_c[0] ? fRcstr_c_yavg[0] /= rec_c[0] : 0);
  // (rec_c[1] ? fRcstr_c_yavg[1] /= rec_c[1] : 0);
  // (rec_c[2] ? fRcstr_c_yavg[2] /= rec_c[2] : 0);
  // (rec_c[3] ? fRcstr_c_yavg[3] /= rec_c[3] : 0);
  // (dis_c[0] ? fNDIS_evt_c_yavg[0] /= dis_c[0] : 0);
  // (dis_c[1] ? fNDIS_evt_c_yavg[1] /= dis_c[1] : 0);
  // (dis_c[2] ? fNDIS_evt_c_yavg[2] /= dis_c[2] : 0);
}

void resetValues()
{
  for(int c=0; c<2; c++)
  {
    for(int i=0; i<9; i++)
    {
      for(int j=0; j<6; j++)
      {
        for(int k=0; k<12; k++)
        {
          if(!c)
          {
            for(int n=0; n<3; n++)
            {
              for(int m=0; m<2; m++)
              {
                fNDIS_evt[n][m][i][j][k]=0; fNDIS_evt_c[n][m][i][j][k]=0; fNDIS_evt_MC[n][m][i][j][k]=0;
                for(int l=0; l<4; l++)
                {
                  fNDIS_evt_zvtx[n][m][i][j][k][l]=0; fNDIS_evt_MC_zvtx[n][m][i][j][k][l]=0;
                }
              }
            }
          }
          for(int n=0; n<3; n++)
          {
            for(int m=0; m<2; m++)
            {
              fRcstr[i][j][k].tab[c][m][0][n]=0; fRcstr_c[i][j][k].tab[c][m][0][n]=0; fGnrt[i][j][k].tab[c][m][0][n]=0;
              for(int l=0; l<4; l++)
              {
                fRcstr_zvtx[i][j][k][l].tab[c][m][0][n]=0; fGnrt_zvtx[i][j][k][l].tab[c][m][0][n]=0;
              }
            }
          }
        }
      }
    }
  }
}

int main(int argc, char **argv)
{

  if(argc < 2)
  {
    cout << "ERROR : Not enough arguments." << endl;
    cout << "Asked : at least 1 *** Received : " << argc-1 << endl;
    cout << "./acccollect periodFile" << endl;

    return 1;
  }

  int year=0;

  if(Y2006) year=2006;
  else if(Y2012) year=2012;
  else if(Y2016) year=2016;

  // Files

  ifstream periods(argv[1]);
  string filelist, periodName;
  int periodBit;
  while(periods >> periodName)
  {
    periods >> periodBit;
    if(!periodBit) continue;

    double dummyd;

    ifstream DIS_file(Form("acceptance/%d/DIS/DIS_%s.txt",year,periodName.c_str()));
    ifstream DIS_zvtx_file(Form("acceptance/%d/DIS/DIS_zvtx_%s.txt",year,periodName.c_str()));
    ifstream had_file(Form("acceptance/%d/hadron/hadron_%s.txt",year,periodName.c_str()));
    ifstream had_zvtx_file(Form("acceptance/%d/hadron/hadron_zvtx_%s.txt",year,periodName.c_str()));
    ifstream had_th_file(Form("acceptance/%d/hadron/hadron_theta_%s.txt",year,periodName.c_str()));
    ifstream had_pt_file(Form("acceptance/%d/hadron/hadron_pt_%s.txt",year,periodName.c_str()));

    for(int c=0; c<2; c++)
    {
      for(int i=0; i<9; i++)
      {
        for(int j=0; j<6; j++)
        {
          for(int k=0; k<12; k++)
          {
            if(!c)
            {
              DIS_file >> dummyd;
              fNDIS_evt[0][1][i][j][k] += dummyd;
              DIS_file >> dummyd;
              fNDIS_evt_c[0][1][i][j][k] += dummyd;
              DIS_file >> dummyd;
              fNDIS_evt_MC[0][1][i][j][k] += dummyd;
              DIS_file >> dummyd;
              fNDIS_evt[1][1][i][j][k] += dummyd;
              DIS_file >> dummyd;
              fNDIS_evt_c[1][1][i][j][k] += dummyd;
              DIS_file >> dummyd;
              fNDIS_evt_MC[1][1][i][j][k] += dummyd;
              DIS_file >> dummyd;
              fNDIS_evt[2][1][i][j][k] += dummyd;
              DIS_file >> dummyd;
              fNDIS_evt_c[2][1][i][j][k] += dummyd;
              DIS_file >> dummyd;
              fNDIS_evt_MC[2][1][i][j][k] += dummyd;

              DIS_file >> dummyd;
              fNDIS_evt[0][0][i][j][k] += dummyd;
              DIS_file >> dummyd;
              fNDIS_evt_c[0][0][i][j][k] += dummyd;
              DIS_file >> dummyd;
              fNDIS_evt_MC[0][0][i][j][k] += dummyd;
              DIS_file >> dummyd;
              fNDIS_evt[1][0][i][j][k] += dummyd;
              DIS_file >> dummyd;
              fNDIS_evt_c[1][0][i][j][k] += dummyd;
              DIS_file >> dummyd;
              fNDIS_evt_MC[1][0][i][j][k] += dummyd;
              DIS_file >> dummyd;
              fNDIS_evt[2][0][i][j][k] += dummyd;
              DIS_file >> dummyd;
              fNDIS_evt_c[2][0][i][j][k] += dummyd;
              DIS_file >> dummyd;
              fNDIS_evt_MC[2][0][i][j][k] += dummyd;

              DIS_zvtx_file >> dummyd;
              fNDIS_evt_zvtx[0][1][i][j][k][0] += dummyd;
              // cout << dummyd << " ";
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_zvtx[0][1][i][j][k][1] += dummyd;
              // cout << dummyd << " ";
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_zvtx[0][1][i][j][k][2] += dummyd;
              // cout << dummyd << " ";
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_zvtx[0][1][i][j][k][3] += dummyd;
              // cout << dummyd << " ";

              DIS_zvtx_file >> dummyd;
              fNDIS_evt_MC_zvtx[0][1][i][j][k][0] += dummyd;
              // cout << dummyd << " ";
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_MC_zvtx[0][1][i][j][k][1] += dummyd;
              // cout << dummyd << " ";
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_MC_zvtx[0][1][i][j][k][2] += dummyd;
              // cout << dummyd << " ";
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_MC_zvtx[0][1][i][j][k][3] += dummyd;
              // cout << dummyd << " ";

              DIS_zvtx_file >> dummyd;
              fNDIS_evt_zvtx[1][1][i][j][k][0] += dummyd;
              // cout << dummyd << " ";
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_zvtx[1][1][i][j][k][1] += dummyd;
              // cout << dummyd << " ";
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_zvtx[1][1][i][j][k][2] += dummyd;
              // cout << dummyd << " ";
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_zvtx[1][1][i][j][k][3] += dummyd;
              // cout << dummyd << " ";

              DIS_zvtx_file >> dummyd;
              fNDIS_evt_MC_zvtx[1][1][i][j][k][0] += dummyd;
              // cout << dummyd << " ";
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_MC_zvtx[1][1][i][j][k][1] += dummyd;
              // cout << dummyd << " ";
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_MC_zvtx[1][1][i][j][k][2] += dummyd;
              // cout << dummyd << " ";
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_MC_zvtx[1][1][i][j][k][3] += dummyd;
              // cout << dummyd << " ";

              DIS_zvtx_file >> dummyd;
              fNDIS_evt_zvtx[2][1][i][j][k][0] += dummyd;
              // cout << dummyd << " ";
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_zvtx[2][1][i][j][k][1] += dummyd;
              // cout << dummyd << " ";
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_zvtx[2][1][i][j][k][2] += dummyd;
              // cout << dummyd << " ";
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_zvtx[2][1][i][j][k][3] += dummyd;
              // cout << dummyd << " ";

              DIS_zvtx_file >> dummyd;
              fNDIS_evt_MC_zvtx[2][1][i][j][k][0] += dummyd;
              // cout << dummyd << " ";
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_MC_zvtx[2][1][i][j][k][1] += dummyd;
              // cout << dummyd << " ";
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_MC_zvtx[2][1][i][j][k][2] += dummyd;
              // cout << dummyd << " ";
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_MC_zvtx[2][1][i][j][k][3] += dummyd;
              // cout << dummyd << " ";

              DIS_zvtx_file >> dummyd;
              fNDIS_evt_zvtx[0][0][i][j][k][0] += dummyd;
              // cout << dummyd << " ";
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_zvtx[0][0][i][j][k][1] += dummyd;
              // cout << dummyd << " ";
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_zvtx[0][0][i][j][k][2] += dummyd;
              // cout << dummyd << " ";
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_zvtx[0][0][i][j][k][3] += dummyd;
              // cout << dummyd << " ";

              DIS_zvtx_file >> dummyd;
              fNDIS_evt_MC_zvtx[0][0][i][j][k][0] += dummyd;
              // cout << dummyd << " ";
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_MC_zvtx[0][0][i][j][k][1] += dummyd;
              // cout << dummyd << " ";
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_MC_zvtx[0][0][i][j][k][2] += dummyd;
              // cout << dummyd << " ";
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_MC_zvtx[0][0][i][j][k][3] += dummyd;
              // cout << dummyd << " ";

              DIS_zvtx_file >> dummyd;
              fNDIS_evt_zvtx[1][0][i][j][k][0] += dummyd;
              // cout << dummyd << " ";
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_zvtx[1][0][i][j][k][1] += dummyd;
              // cout << dummyd << " ";
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_zvtx[1][0][i][j][k][2] += dummyd;
              // cout << dummyd << " ";
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_zvtx[1][0][i][j][k][3] += dummyd;
              // cout << dummyd << " ";

              DIS_zvtx_file >> dummyd;
              fNDIS_evt_MC_zvtx[1][0][i][j][k][0] += dummyd;
              // cout << dummyd << " ";
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_MC_zvtx[1][0][i][j][k][1] += dummyd;
              // cout << dummyd << " ";
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_MC_zvtx[1][0][i][j][k][2] += dummyd;
              // cout << dummyd << " ";
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_MC_zvtx[1][0][i][j][k][3] += dummyd;
              // cout << dummyd << " ";

              DIS_zvtx_file >> dummyd;
              fNDIS_evt_zvtx[2][0][i][j][k][0] += dummyd;
              // cout << dummyd << " ";
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_zvtx[2][0][i][j][k][1] += dummyd;
              // cout << dummyd << " ";
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_zvtx[2][0][i][j][k][2] += dummyd;
              // cout << dummyd << " ";
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_zvtx[2][0][i][j][k][3] += dummyd;
              // cout << dummyd << " ";

              DIS_zvtx_file >> dummyd;
              fNDIS_evt_MC_zvtx[2][0][i][j][k][0] += dummyd;
              // cout << dummyd << " ";
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_MC_zvtx[2][0][i][j][k][1] += dummyd;
              // cout << dummyd << " ";
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_MC_zvtx[2][0][i][j][k][2] += dummyd;
              // cout << dummyd << " ";
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_MC_zvtx[2][0][i][j][k][3] += dummyd;
              // cout << dummyd << " ";

              // cout <<endl;

            }

            for(int ll=0; ll<4; ll++)
            {
              had_file >> dummyd;
              fRcstr[i][j][k].tab[c][1][0][ll] += dummyd;
              had_file >> dummyd;
              fRcstr_c[i][j][k].tab[c][1][0][ll] += dummyd;
              had_file >> dummyd;
              fGnrt[i][j][k].tab[c][1][0][ll] += dummyd;

              had_zvtx_file >> dummyd;
              fRcstr_zvtx[i][j][k][0].tab[c][1][0][ll] += dummyd;
              // cout << dummyd << " ";
              had_zvtx_file >> dummyd;
              fRcstr_zvtx[i][j][k][1].tab[c][1][0][ll] += dummyd;
              // cout << dummyd << " ";
              had_zvtx_file >> dummyd;
              fRcstr_zvtx[i][j][k][2].tab[c][1][0][ll] += dummyd;
              // cout << dummyd << " ";
              had_zvtx_file >> dummyd;
              fRcstr_zvtx[i][j][k][3].tab[c][1][0][ll] += dummyd;
              // cout << dummyd << " ";
              had_zvtx_file >> dummyd;
              fGnrt_zvtx[i][j][k][0].tab[c][1][0][ll] += dummyd;
              // cout << dummyd << " ";
              had_zvtx_file >> dummyd;
              fGnrt_zvtx[i][j][k][1].tab[c][1][0][ll] += dummyd;
              // cout << dummyd << " ";
              had_zvtx_file >> dummyd;
              fGnrt_zvtx[i][j][k][2].tab[c][1][0][ll] += dummyd;
              // cout << dummyd << " ";
              had_zvtx_file >> dummyd;
              fGnrt_zvtx[i][j][k][3].tab[c][1][0][ll] += dummyd;
              // cout << dummyd << " ";

              for(int th=0; th<8; th++)
              {
                had_th_file >> dummyd;
                fRcstr_theta[i][j][k][th].tab[c][1][0][ll] += dummyd;
              }
              for(int th=0; th<8; th++)
              {
                had_th_file >> dummyd;
                fGnrt_theta[i][j][k][th].tab[c][1][0][ll] += dummyd;
              }

              for(int pt=0; pt<10; pt++)
              {
                had_pt_file >> dummyd;
                fRcstr_pt[i][j][k][pt].tab[c][1][0][ll] += dummyd;
              }
              for(int pt=0; pt<10; pt++)
              {
                had_pt_file >> dummyd;
                fGnrt_pt[i][j][k][pt].tab[c][1][0][ll] += dummyd;
              }
            }

            for(int ll=0; ll<4; ll++)
            {
              had_file >> dummyd;
              fRcstr[i][j][k].tab[c][0][0][ll] += dummyd;
              had_file >> dummyd;
              fRcstr_c[i][j][k].tab[c][0][0][ll] += dummyd;
              had_file >> dummyd;
              fGnrt[i][j][k].tab[c][0][0][ll] += dummyd;

              had_zvtx_file >> dummyd;
              fRcstr_zvtx[i][j][k][0].tab[c][0][0][ll] += dummyd;
              // cout << dummyd << " ";
              had_zvtx_file >> dummyd;
              fRcstr_zvtx[i][j][k][1].tab[c][0][0][ll] += dummyd;
              // cout << dummyd << " ";
              had_zvtx_file >> dummyd;
              fRcstr_zvtx[i][j][k][2].tab[c][0][0][ll] += dummyd;
              // cout << dummyd << " ";
              had_zvtx_file >> dummyd;
              fRcstr_zvtx[i][j][k][3].tab[c][0][0][ll] += dummyd;
              // cout << dummyd << " ";
              had_zvtx_file >> dummyd;
              fGnrt_zvtx[i][j][k][0].tab[c][0][0][ll] += dummyd;
              // cout << dummyd << " ";
              had_zvtx_file >> dummyd;
              fGnrt_zvtx[i][j][k][1].tab[c][0][0][ll] += dummyd;
              // cout << dummyd << " ";
              had_zvtx_file >> dummyd;
              fGnrt_zvtx[i][j][k][2].tab[c][0][0][ll] += dummyd;
              // cout << dummyd << " ";
              had_zvtx_file >> dummyd;
              fGnrt_zvtx[i][j][k][3].tab[c][0][0][ll] += dummyd;
              // cout << dummyd << " ";
              // cout << endl;

              for(int th=0; th<8; th++)
              {
                had_th_file >> dummyd;
                fRcstr_theta[i][j][k][th].tab[c][0][0][ll] += dummyd;
              }
              for(int th=0; th<8; th++)
              {
                had_th_file >> dummyd;
                fGnrt_theta[i][j][k][th].tab[c][0][0][ll] += dummyd;
              }

              for(int pt=0; pt<10; pt++)
              {
                had_pt_file >> dummyd;
                fRcstr_pt[i][j][k][pt].tab[c][0][0][ll] += dummyd;
              }
              for(int pt=0; pt<10; pt++)
              {
                had_pt_file >> dummyd;
                fGnrt_pt[i][j][k][pt].tab[c][0][0][ll] += dummyd;
              }
            }
          }
        }
      }
    }

    DIS_file.close();
    had_file.close();
    DIS_zvtx_file.close();
    had_zvtx_file.close();
    had_th_file.close();
    had_pt_file.close();

    TCanvas c5("Hadron_Acceptance","Hadron_Acceptance",3200,1600);
    TCanvas c6("Pion_Acceptance","Pion_Acceptance",3200,1600);
    TCanvas c7("Kaon_Acceptance","Kaon_Acceptance",3200,1600);
    TCanvas c88("Proton_Acceptance","Proton_Acceptance",3200,1600);
    TCanvas* c8[12];
    TCanvas c9("Hadron_Acceptance_yavg","Hadron_Acceptance_yavg",3200,1600);
    TCanvas c10("Pion_Acceptance_yavg","Pion_Acceptance_yavg",3200,1600);
    TCanvas c11("Kaon_Acceptance_yavg","Kaon_Acceptance_yavg",3200,1600);

    c5.SetFillColor(0);
    c6.SetFillColor(0);
    c7.SetFillColor(0);
    c88.SetFillColor(0);
    c9.SetFillColor(0);
    c10.SetFillColor(0);
    c11.SetFillColor(0);

    if(SPREAD)
    {
      c5.Divide(9,5,0,0);
      c6.Divide(9,5,0,0);
      c7.Divide(9,5,0,0);
      c88.Divide(9,5,0,0);
    }
    else
    {
      c5.Divide(5,2,0,0);
      c6.Divide(5,2,0,0);
      c7.Divide(5,2,0,0);
      c88.Divide(5,2,0,0);
    }
    c9.Divide(5,2,0,0);
    c10.Divide(5,2,0,0);
    c11.Divide(5,2,0,0);

    TGraphErrors* H_acc[2][9][6][2];
    TGraphErrors* P_acc[2][9][6][2];
    TGraphErrors* K_acc[2][9][6][2];
    TGraphErrors* PR_acc[2][9][6][2];
    TGraphErrors* H_y[2][9][2];
    TGraphErrors* P_y[2][9][2];
    TGraphErrors* K_y[2][9][2];
    TGraphErrors* H_corr_zvtx[2][9][6][12][2];
    TGraphErrors* P_corr_zvtx[2][9][6][12][2];
    TGraphErrors* K_corr_zvtx[2][9][6][12][2];

    for(int i=0; i<12; i++)
    {
      c8[i] = new TCanvas(Form("Hadron_Acceptance_zvtx_%d",i),Form("Hadron_Acceptance_zvtx_%d",i),3200,1600);
      c8[i]->SetFillColor(0);
      c8[i]->Divide(5,2,0,0);
    }

    double z_range[12] = {.225,.275,.325,.375,.425,.475,.525,.575,.625,.675,.725,.8};
    double zvtx_range[4] = {-281.19,-221.19,-161.19,-101.19};

    ofstream ofs(Form("%s/%d/acceptance_%s.txt",dirroot,year,periodName.c_str()), std::ofstream::out | std::ofstream::trunc);
    ofstream ofs_yavg(Form("%s/%d/acceptance_yavg_%s.txt",dirroot,year,periodName.c_str()), std::ofstream::out | std::ofstream::trunc);
    ofstream ofs_zvtx(Form("%s/%d/acceptance_vtx_%s.txt",dirroot,year,periodName.c_str()), std::ofstream::out | std::ofstream::trunc);
    ofstream ofs_theta(Form("%s/%d/acceptance_theta_%s.txt",dirroot,year,periodName.c_str()), std::ofstream::out | std::ofstream::trunc);
    ofstream ofs_pt(Form("%s/%d/acceptance_pt_%s.txt",dirroot,year,periodName.c_str()), std::ofstream::out | std::ofstream::trunc);
    ofstream ofs_reld(Form("%s/%d/reldiff_vtx_%s.txt",dirroot,year,periodName.c_str()), std::ofstream::out | std::ofstream::trunc);
    ofstream ofs_test("test.txt", std::ofstream::out | std::ofstream::trunc);
    ofstream lepto(Form("%s/%d/lepto_%s.txt",dirroot,year,periodName.c_str()), std::ofstream::out | std::ofstream::trunc);

    for(int c=0; c<2; c++)
    {
      for(int i=0; i<9; i++)
      {
        for(int m=0; m<2; m++)
        {
          std::vector<double> p_y;
          std::vector<double> k_y;
          std::vector<double> h_y;
          std::vector<double> p_y_err;
          std::vector<double> k_y_err;
          std::vector<double> h_y_err;

          std::vector<double> z_range_p_y;
          std::vector<double> z_range_k_y;
          std::vector<double> z_range_h_y;

          for(int l=0; l<12; l++)
          {
            z_range_p_y.push_back(z_range[l]);
            z_range_k_y.push_back(z_range[l]);
            z_range_h_y.push_back(z_range[l]);
          }

          for(int j=0; j<6; j++)
          {
              std::vector<double> p_a;
              std::vector<double> k_a;
              std::vector<double> h_a;
              std::vector<double> pr_a;
              std::vector<double> p_err;
              std::vector<double> k_err;
              std::vector<double> h_err;
              std::vector<double> pr_err;

              std::vector<double> z_range_p;
              std::vector<double> z_range_k;
              std::vector<double> z_range_h;
              std::vector<double> z_range_pr;

              for(int l=0; l<12; l++)
              {
                z_range_p.push_back(z_range[l]);
                z_range_k.push_back(z_range[l]);
                z_range_h.push_back(z_range[l]);
                z_range_pr.push_back(z_range[l]);
              }

              for(int k=0; k<12; k++)
              {
                std::vector<double> p_corr;
                std::vector<double> k_corr;
                std::vector<double> h_corr;
                std::vector<double> p_cerr;
                std::vector<double> k_cerr;
                std::vector<double> h_cerr;

                std::vector<double> zvtx_range_p;
                std::vector<double> zvtx_range_k;
                std::vector<double> zvtx_range_h;

                for(int l=0; l<4; l++)
                {
                  zvtx_range_p.push_back(zvtx_range[l]);
                  zvtx_range_k.push_back(zvtx_range[l]);
                  zvtx_range_h.push_back(zvtx_range[l]);
                }

                /*cout << "data count : " << fNDIS_evt_MC[0][i][j][k] << endl;
                cout << "MC count : " << fNDIS_evt_MC[0][i][j][k] << endl;
                cout << "Rcst count : " << fRcstr[i][j][k].tab[c][0][3] << endl;
                cout << "Gnrt count : " << fGnrt[i][j][k].tab[c][0][3] << endl;*/

                //cout << fRcstr[i][j][k].tab[c][0][2] << " " << fGnrt[i][j][k].tab[c][0][2] << endl;

                fAcceptance[i][j][k].tab[c][m][0][0] = ((fNDIS_evt[0][m][i][j][k] && fNDIS_evt_MC[0][m][i][j][k] && fGnrt[i][j][k].tab[c][m][0][0]) ? double((fRcstr[i][j][k].tab[c][m][0][0]/fNDIS_evt[0][m][i][j][k])/(fGnrt[i][j][k].tab[c][m][0][0]/fNDIS_evt_MC[0][m][i][j][k])) : 0);
                fAcceptance[i][j][k].tab[c][m][0][1] = ((fNDIS_evt[1][m][i][j][k] && fNDIS_evt_MC[1][m][i][j][k] && fGnrt[i][j][k].tab[c][m][0][1]) ? double((fRcstr[i][j][k].tab[c][m][0][1]/fNDIS_evt[1][m][i][j][k])/(fGnrt[i][j][k].tab[c][m][0][1]/fNDIS_evt_MC[1][m][i][j][k])) : 0);
                fAcceptance[i][j][k].tab[c][m][0][2] = ((fNDIS_evt[2][m][i][j][k] && fNDIS_evt_MC[2][m][i][j][k] && fGnrt[i][j][k].tab[c][m][0][2]) ? double((fRcstr[i][j][k].tab[c][m][0][2]/fNDIS_evt[2][m][i][j][k])/(fGnrt[i][j][k].tab[c][m][0][2]/fNDIS_evt_MC[2][m][i][j][k])) : 0);
                fAcceptance[i][j][k].tab[c][m][0][3] = ((fNDIS_evt[0][m][i][j][k] && fNDIS_evt_MC[0][m][i][j][k] && fGnrt[i][j][k].tab[c][m][0][3]) ? double((fRcstr[i][j][k].tab[c][m][0][3]/fNDIS_evt[0][m][i][j][k])/(fGnrt[i][j][k].tab[c][m][0][3]/fNDIS_evt_MC[0][m][i][j][k])) : 0);

                if(fAcceptance[i][j][k].tab[c][m][0][0]<0) fAcceptance[i][j][k].tab[c][m][0][0]=0;
                if(fAcceptance[i][j][k].tab[c][m][0][1]<0) fAcceptance[i][j][k].tab[c][m][0][1]=0;
                if(fAcceptance[i][j][k].tab[c][m][0][2]<0) fAcceptance[i][j][k].tab[c][m][0][2]=0;
                if(fAcceptance[i][j][k].tab[c][m][0][3]<0) fAcceptance[i][j][k].tab[c][m][0][3]=0;

                //cout << fAcceptance[i][j][k].tab[c][0][0] << " " << fAcceptance[i][j][k].tab[c][0][1] << " " << fAcceptance[i][j][k].tab[c][0][2] << " " << fAcceptance[i][j][k].tab[c][0][3] << endl;

                double rd_prime[3];
                double rd[3];
                double gd[3];
                double rh_prime[4];
                double rh[4];
                double gh[4];

                rd_prime[0] = fNDIS_evt[0][m][i][j][k] - fNDIS_evt_c[0][m][i][j][k];
                rd_prime[1] = fNDIS_evt[1][m][i][j][k] - fNDIS_evt_c[1][m][i][j][k];
                rd_prime[2] = fNDIS_evt[2][m][i][j][k] - fNDIS_evt_c[2][m][i][j][k];

                rd[0] = fNDIS_evt_c[0][m][i][j][k];
                rd[1] = fNDIS_evt_c[1][m][i][j][k];
                rd[2] = fNDIS_evt_c[2][m][i][j][k];

                gd[0] = fNDIS_evt_MC[0][m][i][j][k];
                gd[1] = fNDIS_evt_MC[1][m][i][j][k];
                gd[2] = fNDIS_evt_MC[2][m][i][j][k];

                rh_prime[0] = fRcstr[i][j][k].tab[c][m][0][0] - fRcstr_c[i][j][k].tab[c][m][0][0];
                rh_prime[1] = fRcstr[i][j][k].tab[c][m][0][1] - fRcstr_c[i][j][k].tab[c][m][0][1];
                rh_prime[2] = fRcstr[i][j][k].tab[c][m][0][2] - fRcstr_c[i][j][k].tab[c][m][0][2];
                rh_prime[3] = fRcstr[i][j][k].tab[c][m][0][3] - fRcstr_c[i][j][k].tab[c][m][0][3];

                rh[0] = fRcstr_c[i][j][k].tab[c][m][0][0];
                rh[1] = fRcstr_c[i][j][k].tab[c][m][0][1];
                rh[2] = fRcstr_c[i][j][k].tab[c][m][0][2];
                rh[3] = fRcstr_c[i][j][k].tab[c][m][0][3];

                gh[0] = fGnrt[i][j][k].tab[c][m][0][0];
                gh[1] = fGnrt[i][j][k].tab[c][m][0][1];
                gh[2] = fGnrt[i][j][k].tab[c][m][0][2];
                gh[3] = fGnrt[i][j][k].tab[c][m][0][3];

                fAcceptance[i][j][k].tab[c][m][1][0] = ((rh[0]>0 && gh[0]>0 && gd[0]>0 && ((rd[0]+rd_prime[0])!=0)) ? double(pow(gd[0]/(rd_prime[0]+rd[0]),2)
                                                                                    *((rh[0]+1)*(gh[0]-rh[0]+1)/(pow(gh[0]+2,2)*(gh[0]+3))+rh_prime[0]/pow(gh[0],2)+pow(rh_prime[0],2)/pow(gh[0],3))
                                                                                    +pow(gd[0]/(rd[0]+rd_prime[0]),4)*pow((rh[0]+rh_prime[0])/gh[0],2)*((rd[0]+1)*(gd[0]-rd[0]+1)/(pow(gd[0]+2,2)*(gd[0]+3))+rd_prime[0]/pow(gd[0],2)+pow(rd_prime[0],2)/pow(gd[0],3))) : 0);
                fAcceptance[i][j][k].tab[c][m][1][1] = ((rh[1]>0 && gh[1]>0 && gd[1]>0 && ((rd[1]+rd_prime[1])!=0)) ? double(pow(gd[1]/(rd_prime[1]+rd[1]),2)
                                                                                    *((rh[1]+1)*(gh[1]-rh[1]+1)/(pow(gh[1]+2,2)*(gh[1]+3))+rh_prime[1]/pow(gh[1],2)+pow(rh_prime[1],2)/pow(gh[1],3))
                                                                                    +pow(gd[1]/(rd[1]+rd_prime[1]),4)*pow((rh[1]+rh_prime[1])/gh[1],2)*((rd[1]+1)*(gd[1]-rd[1]+1)/(pow(gd[1]+2,2)*(gd[1]+3))+rd_prime[1]/pow(gd[1],2)+pow(rd_prime[1],2)/pow(gd[1],3))) : 0);
                fAcceptance[i][j][k].tab[c][m][1][2] = ((rh[2]>0 && gh[2]>0 && gd[2]>0 && ((rd[2]+rd_prime[2])!=0)) ? double(pow(gd[2]/(rd_prime[2]+rd[2]),2)
                                                                                    *((rh[2]+1)*(gh[2]-rh[2]+1)/(pow(gh[2]+2,2)*(gh[2]+3))+rh_prime[2]/pow(gh[2],2)+pow(rh_prime[2],2)/pow(gh[2],3))
                                                                                    +pow(gd[2]/(rd[2]+rd_prime[2]),4)*pow((rh[2]+rh_prime[2])/gh[2],2)*((rd[2]+1)*(gd[2]-rd[2]+1)/(pow(gd[2]+2,2)*(gd[2]+3))+rd_prime[2]/pow(gd[2],2)+pow(rd_prime[2],2)/pow(gd[2],3))) : 0);
                fAcceptance[i][j][k].tab[c][m][1][3] = ((rh[3]>0 && gh[3]>0 && gd[0]>0 && ((rd[0]+rd_prime[0])!=0)) ? double(pow(gd[0]/(rd_prime[0]+rd[0]),2)
                                                                                    *((rh[3]+1)*(gh[3]-rh[3]+1)/(pow(gh[3]+2,2)*(gh[3]+3))+rh_prime[3]/pow(gh[3],2)+pow(rh_prime[3],2)/pow(gh[3],3))
                                                                                    +pow(gd[0]/(rd[0]+rd_prime[0]),4)*pow((rh[3]+rh_prime[3])/gh[3],2)*((rd[0]+1)*(gd[0]-rd[0]+1)/(pow(gd[0]+2,2)*(gd[0]+3))+rd_prime[0]/pow(gd[0],2)+pow(rd_prime[0],2)/pow(gd[0],3))) : 0);

                if(fAcceptance[i][j][k].tab[c][m][0][0]==0) fAcceptance[i][j][k].tab[c][m][1][0]=0;
                if(fAcceptance[i][j][k].tab[c][m][0][1]==0) fAcceptance[i][j][k].tab[c][m][1][1]=0;
                if(fAcceptance[i][j][k].tab[c][m][0][2]==0) fAcceptance[i][j][k].tab[c][m][1][2]=0;
                if(fAcceptance[i][j][k].tab[c][m][0][3]==0) fAcceptance[i][j][k].tab[c][m][1][3]=0;

                if((fAcceptance[i][j][k].tab[c][m][0][3] != 0) && (fAcceptance[i][j][k].tab[c][m][1][3] > 0.2*fAcceptance[i][j][k].tab[c][m][0][3]))
                {
                  fAcceptance[i][j][k].tab[c][m][0][0] = 0;
                  fAcceptance[i][j][k].tab[c][m][0][1] = 0;
                  fAcceptance[i][j][k].tab[c][m][0][2] = 0;
                  fAcceptance[i][j][k].tab[c][m][0][3] = 0;
                  fAcceptance[i][j][k].tab[c][m][1][0] = 0;
                  fAcceptance[i][j][k].tab[c][m][1][1] = 0;
                  fAcceptance[i][j][k].tab[c][m][1][2] = 0;
                  fAcceptance[i][j][k].tab[c][m][1][3] = 0;
                }

                if((i==0 && j==4 && k>4)
                  ||(i==0 && j==3 && k>9)
                  ||(i==1 && j==4 && k>4)
                  ||(i==1 && j==2 && k<2)
                  ||(i==1 && j==1 && k<4)
                  ||(i==2 && j==4 && k>4)
                  ||(i==2 && j==2 && k<2)
                  ||(i==2 && j==1 && k<4)
                  ||(i==2 && j==0 && k<7)
                  ||(i==3 && j==4 && k>4)
                  ||(i==3 && j==2 && k<2)
                  ||(i==3 && j==1 && k<4)
                  ||(i==3 && j==0 && k<7)
                  ||(i==4 && j==4 && k>4)
                  ||(i==4 && j==2 && k<2)
                  ||(i==4 && j==1 && k<4)
                  ||(i==4 && j==0 && k<7)
                  ||(i==5 && j==4 && k>4)
                  ||(i==5 && j==2 && k<2)
                  ||(i==5 && j==1 && k<4)
                  ||(i==5 && j==0 && k<7)
                  ||(i==6 && j==4 && k>4)
                  ||(i==6 && j==3 && k>10)
                  ||(i==6 && j==2 && k<2)
                  ||(i==6 && j==1 && k<4)
                  ||(i==6 && j==0 && k<7)
                  ||(i==7 && j==4)
                  ||(i==7 && j==3 && k>10)
                  ||(i==7 && j==2 && k<2)
                  ||(i==7 && j==1 && k<4)
                  ||(i==7 && j==0 && k<7)
                  ||(i==8 && j==4)
                  ||(i==8 && j==3 && k>8)
                  ||(i==8 && j==2 && k<2)
                  ||(i==8 && j==1 && k<4)
                  ||(i==8 && j==0)
                )
                {
                   fAcceptance[i][j][k].tab[c][m][0][0] = 0;
                   fAcceptance[i][j][k].tab[c][m][0][1] = 0;
                   fAcceptance[i][j][k].tab[c][m][0][2] = 0;
                   fAcceptance[i][j][k].tab[c][m][0][3] = 0;
                   fAcceptance[i][j][k].tab[c][m][1][0] = 0;
                   fAcceptance[i][j][k].tab[c][m][1][1] = 0;
                   fAcceptance[i][j][k].tab[c][m][1][2] = 0;
                   fAcceptance[i][j][k].tab[c][m][1][3] = 0;
                }

                fAcceptance[i][j][k].tab[c][m][0][0] = ((fAcceptance[i][j][k].tab[c][m][0][0]) ? (fAcceptance[i][j][k].tab[c][m][0][0]) : 0);
                fAcceptance[i][j][k].tab[c][m][0][1] = ((fAcceptance[i][j][k].tab[c][m][0][1]) ? (fAcceptance[i][j][k].tab[c][m][0][1]) : 0);
                fAcceptance[i][j][k].tab[c][m][0][2] = ((fAcceptance[i][j][k].tab[c][m][0][2]) ? (fAcceptance[i][j][k].tab[c][m][0][2]) : 0);
                fAcceptance[i][j][k].tab[c][m][0][3] = ((fAcceptance[i][j][k].tab[c][m][0][3]) ? (fAcceptance[i][j][k].tab[c][m][0][3]) : 0);

                p_a.push_back(fAcceptance[i][j][k].tab[c][m][0][0]);
                k_a.push_back(fAcceptance[i][j][k].tab[c][m][0][1]);
                pr_a.push_back(fAcceptance[i][j][k].tab[c][m][0][2]);
                h_a.push_back(fAcceptance[i][j][k].tab[c][m][0][3]);

                p_err.push_back(fAcceptance[i][j][k].tab[c][m][1][0]);
                k_err.push_back(fAcceptance[i][j][k].tab[c][m][1][1]);
                pr_err.push_back(fAcceptance[i][j][k].tab[c][m][1][2]);
                h_err.push_back(fAcceptance[i][j][k].tab[c][m][1][3]);

                for(int l=0; l<4; l++)
                {
                  fAcceptance_zvtx[i][j][k][l].tab[c][m][0][0] = ((fNDIS_evt_zvtx[0][m][i][j][k][l] && fNDIS_evt_MC_zvtx[0][m][i][j][k][l] && fGnrt_zvtx[i][j][k][l].tab[c][m][0][0]) ? double((fRcstr_zvtx[i][j][k][l].tab[c][m][0][0]/fNDIS_evt_zvtx[0][m][i][j][k][l])/(fGnrt_zvtx[i][j][k][l].tab[c][m][0][0]/fNDIS_evt_MC_zvtx[0][m][i][j][k][l])) : 0);
                  fAcceptance_zvtx[i][j][k][l].tab[c][m][0][1] = ((fNDIS_evt_zvtx[1][m][i][j][k][l] && fNDIS_evt_MC_zvtx[1][m][i][j][k][l] && fGnrt_zvtx[i][j][k][l].tab[c][m][0][1]) ? double((fRcstr_zvtx[i][j][k][l].tab[c][m][0][1]/fNDIS_evt_zvtx[1][m][i][j][k][l])/(fGnrt_zvtx[i][j][k][l].tab[c][m][0][1]/fNDIS_evt_MC_zvtx[1][m][i][j][k][l])) : 0);
                  fAcceptance_zvtx[i][j][k][l].tab[c][m][0][2] = ((fNDIS_evt_zvtx[2][m][i][j][k][l] && fNDIS_evt_MC_zvtx[2][m][i][j][k][l] && fGnrt_zvtx[i][j][k][l].tab[c][m][0][2]) ? double((fRcstr_zvtx[i][j][k][l].tab[c][m][0][2]/fNDIS_evt_zvtx[2][m][i][j][k][l])/(fGnrt_zvtx[i][j][k][l].tab[c][m][0][2]/fNDIS_evt_MC_zvtx[2][m][i][j][k][l])) : 0);
                  fAcceptance_zvtx[i][j][k][l].tab[c][m][0][3] = ((fNDIS_evt_zvtx[0][m][i][j][k][l] && fNDIS_evt_MC_zvtx[0][m][i][j][k][l] && fGnrt_zvtx[i][j][k][l].tab[c][m][0][3]) ? double((fRcstr_zvtx[i][j][k][l].tab[c][m][0][3]/fNDIS_evt_zvtx[0][m][i][j][k][l])/(fGnrt_zvtx[i][j][k][l].tab[c][m][0][3]/fNDIS_evt_MC_zvtx[0][m][i][j][k][l])) : 0);

                  if(fAcceptance_zvtx[i][j][k][l].tab[c][m][0][0]<0) fAcceptance_zvtx[i][j][k][l].tab[c][m][0][0]=0;
                  if(fAcceptance_zvtx[i][j][k][l].tab[c][m][0][1]<0) fAcceptance_zvtx[i][j][k][l].tab[c][m][0][1]=0;
                  if(fAcceptance_zvtx[i][j][k][l].tab[c][m][0][2]<0) fAcceptance_zvtx[i][j][k][l].tab[c][m][0][2]=0;
                  if(fAcceptance_zvtx[i][j][k][l].tab[c][m][0][3]<0) fAcceptance_zvtx[i][j][k][l].tab[c][m][0][3]=0;

                  fAcceptance_zvtx[i][j][k][l].tab[c][m][1][0] = fAcceptance[i][j][k].tab[c][m][1][0];
                  fAcceptance_zvtx[i][j][k][l].tab[c][m][1][1] = fAcceptance[i][j][k].tab[c][m][1][1];
                  fAcceptance_zvtx[i][j][k][l].tab[c][m][1][2] = fAcceptance[i][j][k].tab[c][m][1][2];
                  fAcceptance_zvtx[i][j][k][l].tab[c][m][1][3] = fAcceptance[i][j][k].tab[c][m][1][3];


                  // if((i==7 && j==4) || (i==8 && j==0) || (i==8 && j==4))
                  // {
                  //   fAcceptance_zvtx[i][j][k][l].tab[c][m][0][0] = 0;
                  //   fAcceptance_zvtx[i][j][k][l].tab[c][m][0][1] = 0;
                  //   fAcceptance_zvtx[i][j][k][l].tab[c][m][0][2] = 0;
                  //   fAcceptance_zvtx[i][j][k][l].tab[c][m][0][3] = 0;
                  //   fAcceptance_zvtx[i][j][k][l].tab[c][m][1][0] = 0;
                  //   fAcceptance_zvtx[i][j][k][l].tab[c][m][1][1] = 0;
                  //   fAcceptance_zvtx[i][j][k][l].tab[c][m][1][2] = 0;
                  //   fAcceptance_zvtx[i][j][k][l].tab[c][m][1][3] = 0;
                  // }

                  p_corr.push_back(fAcceptance_zvtx[i][j][k][l].tab[c][m][0][0] ? fAcceptance_zvtx[i][j][k][l].tab[c][m][0][0] : 0);
                  k_corr.push_back(fAcceptance_zvtx[i][j][k][l].tab[c][m][0][1] ? fAcceptance_zvtx[i][j][k][l].tab[c][m][0][1] : 0);
                  h_corr.push_back(fAcceptance_zvtx[i][j][k][l].tab[c][m][0][3] ? fAcceptance_zvtx[i][j][k][l].tab[c][m][0][3] : 0);

                  p_cerr.push_back(sqrt(fAcceptance_zvtx[i][j][k][l].tab[c][m][1][0]));
                  k_cerr.push_back(sqrt(fAcceptance_zvtx[i][j][k][l].tab[c][m][1][1]));
                  h_cerr.push_back(sqrt(fAcceptance_zvtx[i][j][k][l].tab[c][m][1][3]));
                }

                for(int th=0; th<8; th++)
                {
                  fAcceptance_theta[i][j][k][th].tab[c][m][0][0] = ((fNDIS_evt[0][m][i][j][k] && fNDIS_evt_MC[0][m][i][j][k] && fGnrt_theta[i][j][k][th].tab[c][m][0][0]) ? double((fRcstr_theta[i][j][k][th].tab[c][m][0][0]/fNDIS_evt[0][m][i][j][k])/(fGnrt_theta[i][j][k][th].tab[c][m][0][0]/fNDIS_evt_MC[0][m][i][j][k])) : 0);
                  fAcceptance_theta[i][j][k][th].tab[c][m][0][1] = ((fNDIS_evt[1][m][i][j][k] && fNDIS_evt_MC[1][m][i][j][k] && fGnrt_theta[i][j][k][th].tab[c][m][0][1]) ? double((fRcstr_theta[i][j][k][th].tab[c][m][0][1]/fNDIS_evt[1][m][i][j][k])/(fGnrt_theta[i][j][k][th].tab[c][m][0][1]/fNDIS_evt_MC[1][m][i][j][k])) : 0);
                  fAcceptance_theta[i][j][k][th].tab[c][m][0][2] = ((fNDIS_evt[2][m][i][j][k] && fNDIS_evt_MC[2][m][i][j][k] && fGnrt_theta[i][j][k][th].tab[c][m][0][2]) ? double((fRcstr_theta[i][j][k][th].tab[c][m][0][2]/fNDIS_evt[2][m][i][j][k])/(fGnrt_theta[i][j][k][th].tab[c][m][0][2]/fNDIS_evt_MC[2][m][i][j][k])) : 0);
                  fAcceptance_theta[i][j][k][th].tab[c][m][0][3] = ((fNDIS_evt[0][m][i][j][k] && fNDIS_evt_MC[0][m][i][j][k] && fGnrt_theta[i][j][k][th].tab[c][m][0][3]) ? double((fRcstr_theta[i][j][k][th].tab[c][m][0][3]/fNDIS_evt[0][m][i][j][k])/(fGnrt_theta[i][j][k][th].tab[c][m][0][3]/fNDIS_evt_MC[0][m][i][j][k])) : 0);

                  if(fAcceptance_theta[i][j][k][th].tab[c][m][0][0]<0) fAcceptance_theta[i][j][k][th].tab[c][m][0][0]=0;
                  if(fAcceptance_theta[i][j][k][th].tab[c][m][0][1]<0) fAcceptance_theta[i][j][k][th].tab[c][m][0][1]=0;
                  if(fAcceptance_theta[i][j][k][th].tab[c][m][0][2]<0) fAcceptance_theta[i][j][k][th].tab[c][m][0][2]=0;
                  if(fAcceptance_theta[i][j][k][th].tab[c][m][0][3]<0) fAcceptance_theta[i][j][k][th].tab[c][m][0][3]=0;

                  fAcceptance_theta[i][j][k][th].tab[c][m][1][0] = fAcceptance[i][j][k].tab[c][m][1][0];
                  fAcceptance_theta[i][j][k][th].tab[c][m][1][1] = fAcceptance[i][j][k].tab[c][m][1][1];
                  fAcceptance_theta[i][j][k][th].tab[c][m][1][2] = fAcceptance[i][j][k].tab[c][m][1][2];
                  fAcceptance_theta[i][j][k][th].tab[c][m][1][3] = fAcceptance[i][j][k].tab[c][m][1][3];


                  // if((i==7 && j==4) || (i==8 && j==0) || (i==8 && j==4))
                  // {
                  //   fAcceptance_theta[i][j][k][th].tab[c][m][0][0] = 0;
                  //   fAcceptance_theta[i][j][k][th].tab[c][m][0][1] = 0;
                  //   fAcceptance_theta[i][j][k][th].tab[c][m][0][2] = 0;
                  //   fAcceptance_theta[i][j][k][th].tab[c][m][0][3] = 0;
                  //   fAcceptance_theta[i][j][k][th].tab[c][m][1][0] = 0;
                  //   fAcceptance_theta[i][j][k][th].tab[c][m][1][1] = 0;
                  //   fAcceptance_theta[i][j][k][th].tab[c][m][1][2] = 0;
                  //   fAcceptance_theta[i][j][k][th].tab[c][m][1][3] = 0;
                  // }

                  // p_th.push_back(fAcceptance_theta[i][j][k][th].tab[c][m][0][0] ? fAcceptance_theta[i][j][k][th].tab[c][m][0][0] : 0);
                  // k_th.push_back(fAcceptance_theta[i][j][k][th].tab[c][m][0][1] ? fAcceptance_theta[i][j][k][th].tab[c][m][0][1] : 0);
                  // h_th.push_back(fAcceptance_theta[i][j][k][th].tab[c][m][0][3] ? fAcceptance_theta[i][j][k][th].tab[c][m][0][3] : 0);
                  //
                  // p_therr.push_back(sqrt(fAcceptance_theta[i][j][k][th].tab[c][m][1][0]));
                  // k_therr.push_back(sqrt(fAcceptance_theta[i][j][k][th].tab[c][m][1][1]));
                  // h_therr.push_back(sqrt(fAcceptance_theta[i][j][k][th].tab[c][m][1][3]));
                }

                for(int pt=0; pt<10; pt++)
                {
                  fAcceptance_pt[i][j][k][pt].tab[c][m][0][0] = ((fNDIS_evt[0][m][i][j][k] && fNDIS_evt_MC[0][m][i][j][k] && fGnrt_pt[i][j][k][pt].tab[c][m][0][0]) ? double((fRcstr_pt[i][j][k][pt].tab[c][m][0][0]/fNDIS_evt[0][m][i][j][k])/(fGnrt_pt[i][j][k][pt].tab[c][m][0][0]/fNDIS_evt_MC[0][m][i][j][k])) : 0);
                  fAcceptance_pt[i][j][k][pt].tab[c][m][0][1] = ((fNDIS_evt[1][m][i][j][k] && fNDIS_evt_MC[1][m][i][j][k] && fGnrt_pt[i][j][k][pt].tab[c][m][0][1]) ? double((fRcstr_pt[i][j][k][pt].tab[c][m][0][1]/fNDIS_evt[1][m][i][j][k])/(fGnrt_pt[i][j][k][pt].tab[c][m][0][1]/fNDIS_evt_MC[1][m][i][j][k])) : 0);
                  fAcceptance_pt[i][j][k][pt].tab[c][m][0][2] = ((fNDIS_evt[2][m][i][j][k] && fNDIS_evt_MC[2][m][i][j][k] && fGnrt_pt[i][j][k][pt].tab[c][m][0][2]) ? double((fRcstr_pt[i][j][k][pt].tab[c][m][0][2]/fNDIS_evt[2][m][i][j][k])/(fGnrt_pt[i][j][k][pt].tab[c][m][0][2]/fNDIS_evt_MC[2][m][i][j][k])) : 0);
                  fAcceptance_pt[i][j][k][pt].tab[c][m][0][3] = ((fNDIS_evt[0][m][i][j][k] && fNDIS_evt_MC[0][m][i][j][k] && fGnrt_pt[i][j][k][pt].tab[c][m][0][3]) ? double((fRcstr_pt[i][j][k][pt].tab[c][m][0][3]/fNDIS_evt[0][m][i][j][k])/(fGnrt_pt[i][j][k][pt].tab[c][m][0][3]/fNDIS_evt_MC[0][m][i][j][k])) : 0);

                  if(fAcceptance_pt[i][j][k][pt].tab[c][m][0][0]<0) fAcceptance_pt[i][j][k][pt].tab[c][m][0][0]=0;
                  if(fAcceptance_pt[i][j][k][pt].tab[c][m][0][1]<0) fAcceptance_pt[i][j][k][pt].tab[c][m][0][1]=0;
                  if(fAcceptance_pt[i][j][k][pt].tab[c][m][0][2]<0) fAcceptance_pt[i][j][k][pt].tab[c][m][0][2]=0;
                  if(fAcceptance_pt[i][j][k][pt].tab[c][m][0][3]<0) fAcceptance_pt[i][j][k][pt].tab[c][m][0][3]=0;

                  fAcceptance_pt[i][j][k][pt].tab[c][m][1][0] = fAcceptance[i][j][k].tab[c][m][1][0];
                  fAcceptance_pt[i][j][k][pt].tab[c][m][1][1] = fAcceptance[i][j][k].tab[c][m][1][1];
                  fAcceptance_pt[i][j][k][pt].tab[c][m][1][2] = fAcceptance[i][j][k].tab[c][m][1][2];
                  fAcceptance_pt[i][j][k][pt].tab[c][m][1][3] = fAcceptance[i][j][k].tab[c][m][1][3];


                  // if((i==7 && j==4) || (i==8 && j==0) || (i==8 && j==4))
                  // {
                  //   fAcceptance_theta[i][j][k][th].tab[c][m][0][0] = 0;
                  //   fAcceptance_theta[i][j][k][th].tab[c][m][0][1] = 0;
                  //   fAcceptance_theta[i][j][k][th].tab[c][m][0][2] = 0;
                  //   fAcceptance_theta[i][j][k][th].tab[c][m][0][3] = 0;
                  //   fAcceptance_theta[i][j][k][th].tab[c][m][1][0] = 0;
                  //   fAcceptance_theta[i][j][k][th].tab[c][m][1][1] = 0;
                  //   fAcceptance_theta[i][j][k][th].tab[c][m][1][2] = 0;
                  //   fAcceptance_theta[i][j][k][th].tab[c][m][1][3] = 0;
                  // }

                  // p_th.push_back(fAcceptance_theta[i][j][k][th].tab[c][m][0][0] ? fAcceptance_theta[i][j][k][th].tab[c][m][0][0] : 0);
                  // k_th.push_back(fAcceptance_theta[i][j][k][th].tab[c][m][0][1] ? fAcceptance_theta[i][j][k][th].tab[c][m][0][1] : 0);
                  // h_th.push_back(fAcceptance_theta[i][j][k][th].tab[c][m][0][3] ? fAcceptance_theta[i][j][k][th].tab[c][m][0][3] : 0);
                  //
                  // p_therr.push_back(sqrt(fAcceptance_theta[i][j][k][th].tab[c][m][1][0]));
                  // k_therr.push_back(sqrt(fAcceptance_theta[i][j][k][th].tab[c][m][1][1]));
                  // h_therr.push_back(sqrt(fAcceptance_theta[i][j][k][th].tab[c][m][1][3]));
                }

                for(int l=4; l>0; l--)
                {
                  if(!p_corr[l-1]) {p_corr.erase(p_corr.begin()+l-1); p_cerr.erase(p_cerr.begin()+l-1); zvtx_range_p.erase(zvtx_range_p.begin()+l-1);}
                  if(!k_corr[l-1]) {k_corr.erase(k_corr.begin()+l-1); k_cerr.erase(k_cerr.begin()+l-1); zvtx_range_k.erase(zvtx_range_k.begin()+l-1);}
                  if(!h_corr[l-1]) {h_corr.erase(h_corr.begin()+l-1); h_cerr.erase(h_cerr.begin()+l-1); zvtx_range_h.erase(zvtx_range_h.begin()+l-1);}
                }

                bool h_corr_empty = 0;

                if(!(int(h_corr.size()))) h_corr_empty = 1;

                H_corr_zvtx[c][i][j][k][m] = new TGraphErrors(int(h_corr.size()),&(zvtx_range_h[0]),&(h_corr[0]),0,&(h_cerr[0]));
                P_corr_zvtx[c][i][j][k][m] = new TGraphErrors(int(p_corr.size()),&(zvtx_range_p[0]),&(p_corr[0]),0,&(p_cerr[0]));
                K_corr_zvtx[c][i][j][k][m] = new TGraphErrors(int(k_corr.size()),&(zvtx_range_k[0]),&(k_corr[0]),0,&(k_cerr[0]));

                H_corr_zvtx[c][i][j][k][m]->SetMarkerColor(fMarkerColor[j]);
                P_corr_zvtx[c][i][j][k][m]->SetMarkerColor(fMarkerColor[j]);
                K_corr_zvtx[c][i][j][k][m]->SetMarkerColor(fMarkerColor[j]);

                H_corr_zvtx[c][i][j][k][m]->SetMarkerSize(3);
                P_corr_zvtx[c][i][j][k][m]->SetMarkerSize(3);
                K_corr_zvtx[c][i][j][k][m]->SetMarkerSize(3);

                H_corr_zvtx[c][i][j][k][m]->SetMarkerStyle(fMarkerStyle[j][c]);
                P_corr_zvtx[c][i][j][k][m]->SetMarkerStyle(fMarkerStyle[j][c]);
                K_corr_zvtx[c][i][j][k][m]->SetMarkerStyle(fMarkerStyle[j][c]);

                H_corr_zvtx[c][i][j][k][m]->GetYaxis()->SetTitle("");
                P_corr_zvtx[c][i][j][k][m]->GetYaxis()->SetTitle("");
                K_corr_zvtx[c][i][j][k][m]->GetYaxis()->SetTitle("");

                H_corr_zvtx[c][i][j][k][m]->GetXaxis()->SetTitle("");
                P_corr_zvtx[c][i][j][k][m]->GetXaxis()->SetTitle("");
                K_corr_zvtx[c][i][j][k][m]->GetXaxis()->SetTitle("");

                H_corr_zvtx[c][i][j][k][m]->SetTitle("");
                P_corr_zvtx[c][i][j][k][m]->SetTitle("");
                K_corr_zvtx[c][i][j][k][m]->SetTitle("");

                if(!h_corr_empty)
                {
                  c8[k]->cd(i+1);
                  gPad->SetFillStyle(4000);
                  if(H_corr_zvtx[c][i][j][k][m])
                  {
                    if(!c && !m && j==3)
                    {
                      H_corr_zvtx[c][i][j][k][m]->Draw("SAMEPA");
                      H_corr_zvtx[c][i][j][k][m]->GetXaxis()->SetLimits(-320,-60);
                      H_corr_zvtx[c][i][j][k][m]->SetMinimum(0.);
                      H_corr_zvtx[c][i][j][k][m]->SetMaximum(2.);
                      H_corr_zvtx[c][i][j][k][m]->GetXaxis()->SetLabelSize(0.06);
                      H_corr_zvtx[c][i][j][k][m]->GetYaxis()->SetLabelSize(0.06);
                      H_corr_zvtx[c][i][j][k][m]->SetTitle("");
                      if(i>4) gPad->SetBottomMargin(.15);
                      if(i==0 || i==5) gPad->SetLeftMargin(.22);
                      if(i==8)
                      {
                        H_corr_zvtx[c][i][j][k][m]->GetXaxis()->SetTitle("#font[ 12]{z_{vtx}}");
                        H_corr_zvtx[c][i][j][k][m]->GetXaxis()->SetTitleSize(0.08);
                        H_corr_zvtx[c][i][j][k][m]->GetXaxis()->SetTitleOffset(.8);
                      }
                      H_corr_zvtx[c][i][j][k][m]->GetXaxis()->SetNdivisions(304,kTRUE);
                      H_corr_zvtx[c][i][j][k][m]->GetYaxis()->SetNdivisions(304,kTRUE);
                      if(i==0)
                      {
                        H_corr_zvtx[c][i][j][k][m]->GetYaxis()->SetTitle("#font[12]{acceptance}^{#font[ 12]{h}}");
                        H_corr_zvtx[c][i][j][k][m]->GetYaxis()->SetTitleSize(0.08);
                      }
                      H_corr_zvtx[c][i][0][k][m]->Draw("SAMEP");
                      H_corr_zvtx[c][i][0][k][m]->GetXaxis()->SetLimits(-320,-60);
                      H_corr_zvtx[c][i][0][k][m]->SetMinimum(0.);
                      H_corr_zvtx[c][i][0][k][m]->SetMaximum(2.);
                      H_corr_zvtx[c][i][1][k][m]->Draw("SAMEP");
                      H_corr_zvtx[c][i][1][k][m]->GetXaxis()->SetLimits(-320,-60);
                      H_corr_zvtx[c][i][1][k][m]->SetMinimum(0.);
                      H_corr_zvtx[c][i][1][k][m]->SetMaximum(2.);
                      H_corr_zvtx[c][i][2][k][m]->Draw("SAMEP");
                      H_corr_zvtx[c][i][2][k][m]->GetXaxis()->SetLimits(-320,-60);
                      H_corr_zvtx[c][i][2][k][m]->SetMinimum(0.);
                      H_corr_zvtx[c][i][2][k][m]->SetMaximum(2.);
                      c8[k]->Range(-320,0.,-60,2.);
                    }
                    else
                    {
                      H_corr_zvtx[c][i][j][k][m]->Draw("SAMEP");
                      H_corr_zvtx[c][i][j][k][m]->GetXaxis()->SetLimits(-320,-60);
                      H_corr_zvtx[c][i][j][k][m]->SetMinimum(0.);
                      H_corr_zvtx[c][i][j][k][m]->SetMaximum(2.);
                    }
                  }
                  c8[k]->Update();
                }

              }

              for(int k=12; k>0; k--)
              {
                if(!p_a[k-1]) {p_a.erase(p_a.begin()+k-1); p_err.erase(p_err.begin()+k-1); z_range_p.erase(z_range_p.begin()+k-1);}
                if(!k_a[k-1]) {k_a.erase(k_a.begin()+k-1); k_err.erase(k_err.begin()+k-1); z_range_k.erase(z_range_k.begin()+k-1);}
                if(!pr_a[k-1]) {pr_a.erase(pr_a.begin()+k-1); pr_err.erase(pr_err.begin()+k-1); z_range_pr.erase(z_range_pr.begin()+k-1);}
                if(!h_a[k-1]) {h_a.erase(h_a.begin()+k-1); h_err.erase(h_err.begin()+k-1); z_range_h.erase(z_range_h.begin()+k-1);}
              }

              bool p_a_empty = 0;
              bool k_a_empty = 0;
              bool pr_a_empty = 0;
              bool h_a_empty = 0;

              if(!(int(p_a.size()))) p_a_empty = 1;
              if(!(int(k_a.size()))) k_a_empty = 1;
              if(!(int(pr_a.size()))) pr_a_empty = 1;
              if(!(int(h_a.size()))) h_a_empty = 1;

              H_acc[c][i][j][m] = new TGraphErrors(int(h_a.size()),&(z_range_h[0]),&(h_a[0]),0,&(h_err[0]));
              P_acc[c][i][j][m] = new TGraphErrors(int(p_a.size()),&(z_range_p[0]),&(p_a[0]),0,&(p_err[0]));
              K_acc[c][i][j][m] = new TGraphErrors(int(k_a.size()),&(z_range_k[0]),&(k_a[0]),0,&(k_err[0]));
              PR_acc[c][i][j][m] = new TGraphErrors(int(pr_a.size()),&(z_range_pr[0]),&(pr_a[0]),0,&(pr_err[0]));

              if(SPREAD)
              {
                if(!c)
                {
                  H_acc[c][i][j][m]->SetMarkerColor(fMarkerColor[4]);
                  P_acc[c][i][j][m]->SetMarkerColor(fMarkerColor[4]);
                  K_acc[c][i][j][m]->SetMarkerColor(fMarkerColor[4]);
                  PR_acc[c][i][j][m]->SetMarkerColor(fMarkerColor[4]);
                }
                else
                {
                  H_acc[c][i][j][m]->SetMarkerColor(fMarkerColor[0]);
                  P_acc[c][i][j][m]->SetMarkerColor(fMarkerColor[0]);
                  K_acc[c][i][j][m]->SetMarkerColor(fMarkerColor[0]);
                  PR_acc[c][i][j][m]->SetMarkerColor(fMarkerColor[0]);
                }
              }
              else
              {
                H_acc[c][i][j][m]->SetMarkerColor(fMarkerColor[j]);
                P_acc[c][i][j][m]->SetMarkerColor(fMarkerColor[j]);
                K_acc[c][i][j][m]->SetMarkerColor(fMarkerColor[j]);
                PR_acc[c][i][j][m]->SetMarkerColor(fMarkerColor[j]);
              }

              H_acc[c][i][j][m]->SetMarkerSize(3);
              P_acc[c][i][j][m]->SetMarkerSize(3);
              K_acc[c][i][j][m]->SetMarkerSize(3);
              PR_acc[c][i][j][m]->SetMarkerSize(3);

              if(SPREAD)
              {
                if(!m)
                {
                  H_acc[c][i][j][m]->SetMarkerStyle(fMarkerStyleb[1][c]);
                  P_acc[c][i][j][m]->SetMarkerStyle(fMarkerStyleb[1][c]);
                  K_acc[c][i][j][m]->SetMarkerStyle(fMarkerStyleb[1][c]);
                  PR_acc[c][i][j][m]->SetMarkerStyle(fMarkerStyleb[1][c]);
                }
                else
                {
                  H_acc[c][i][j][m]->SetMarkerStyle(fMarkerStyleb[0][c]);
                  P_acc[c][i][j][m]->SetMarkerStyle(fMarkerStyleb[0][c]);
                  K_acc[c][i][j][m]->SetMarkerStyle(fMarkerStyleb[0][c]);
                  PR_acc[c][i][j][m]->SetMarkerStyle(fMarkerStyleb[0][c]);
                }
              }
              else
              {
                H_acc[c][i][j][m]->SetMarkerStyle(fMarkerStyle[j][c]);
                P_acc[c][i][j][m]->SetMarkerStyle(fMarkerStyle[j][c]);
                K_acc[c][i][j][m]->SetMarkerStyle(fMarkerStyle[j][c]);
                PR_acc[c][i][j][m]->SetMarkerStyle(fMarkerStyle[j][c]);
              }

              H_acc[c][i][j][m]->GetYaxis()->SetTitle("");
              P_acc[c][i][j][m]->GetYaxis()->SetTitle("");
              K_acc[c][i][j][m]->GetYaxis()->SetTitle("");
              PR_acc[c][i][j][m]->GetYaxis()->SetTitle("");

              H_acc[c][i][j][m]->GetXaxis()->SetTitle("");
              P_acc[c][i][j][m]->GetXaxis()->SetTitle("");
              K_acc[c][i][j][m]->GetXaxis()->SetTitle("");
              PR_acc[c][i][j][m]->GetXaxis()->SetTitle("");

              H_acc[c][i][j][m]->SetTitle("");
              P_acc[c][i][j][m]->SetTitle("");
              K_acc[c][i][j][m]->SetTitle("");
              PR_acc[c][i][j][m]->SetTitle("");

              if(SPREAD)
              {
                if(!h_a_empty)
                {
                  c5.cd(i+1+9*j);
                  gPad->SetFillStyle(4000);
                  if(H_acc[c][i][j][m])
                  {
                    if(!c && !m)
                    {
                      H_acc[c][i][j][m]->Draw("SAMEPA");
                      H_acc[c][i][j][m]->GetXaxis()->SetLimits(0.1,0.9);
                      H_acc[c][i][j][m]->SetMinimum(0.25);
                      H_acc[c][i][j][m]->SetMaximum(1.04);
                      H_acc[c][i][j][m]->GetXaxis()->SetLabelSize(0.06);
                      H_acc[c][i][j][m]->GetYaxis()->SetLabelSize(0.06);
                      H_acc[c][i][j][m]->SetTitle("");
                      if(j==5) gPad->SetBottomMargin(.15);
                      if(i==0) gPad->SetLeftMargin(.22);
                      if(i==8 && j==5)
                      {
                        H_acc[c][i][j][m]->GetXaxis()->SetTitle("#font[ 12]{z}");
                        H_acc[c][i][j][m]->GetXaxis()->SetTitleSize(0.08);
                        H_acc[c][i][j][m]->GetXaxis()->SetTitleOffset(.8);
                      }
                      H_acc[c][i][j][m]->GetXaxis()->SetNdivisions(304,kTRUE);
                      H_acc[c][i][j][m]->GetYaxis()->SetNdivisions(304,kTRUE);
                      if(i==1 && j==0)
                      {
                        H_acc[c][i][j][m]->GetYaxis()->SetTitle("#font[12]{acceptance}^{#font[ 12]{h}}");
                        H_acc[c][i][j][m]->GetYaxis()->SetTitleSize(0.08);
                      }
                      H_acc[c][i][j][m]->Draw("SAMEP");
                      H_acc[c][i][j][m]->GetXaxis()->SetLimits(0.1,0.9);
                      H_acc[c][i][j][m]->SetMinimum(0.25);
                      H_acc[c][i][j][m]->SetMaximum(1.04);
                      c5.Range(0.1,0.,0.9,1.2);
                    }
                    else
                    {
                      H_acc[c][i][j][m]->Draw("SAMEP");
                      H_acc[c][i][j][m]->GetXaxis()->SetLimits(0.1,0.9);
                      H_acc[c][i][j][m]->SetMinimum(0.25);
                      H_acc[c][i][j][m]->SetMaximum(1.04);
                    }
                  }
                  c5.Update();
                }

                if(!p_a_empty)
                {
                  c6.cd(i+1+9*j);
                  if(P_acc[c][i][j][m])
                  {
                    if(!c && !m)
                    {
                      P_acc[c][i][j][m]->Draw("SAMEPA");
                      P_acc[c][i][j][m]->GetXaxis()->SetLimits(0.1,0.9);
                      P_acc[c][i][j][m]->SetMinimum(0.25);
                      P_acc[c][i][j][m]->SetMaximum(1.04);
                      P_acc[c][i][j][m]->GetXaxis()->SetLabelSize(0.06);
                      P_acc[c][i][j][m]->GetYaxis()->SetLabelSize(0.06);
                      P_acc[c][i][j][m]->SetTitle("");
                      if(j==5) gPad->SetBottomMargin(.15);
                      if(i==0) gPad->SetLeftMargin(.22);
                      if(i==8 && j==5)
                      {
                        P_acc[c][i][j][m]->GetXaxis()->SetTitle("#font[ 12]{z}");
                        P_acc[c][i][j][m]->GetXaxis()->SetTitleSize(0.08);
                        P_acc[c][i][j][m]->GetXaxis()->SetTitleOffset(.8);
                      }
                      P_acc[c][i][j][m]->GetXaxis()->SetNdivisions(304,kTRUE);
                      P_acc[c][i][j][m]->GetYaxis()->SetNdivisions(304,kTRUE);
                      if(i==2 && j==0)
                      {
                        P_acc[c][i][j][m]->GetYaxis()->SetTitle("#font[12]{acceptance}^{#font[ 12]{#pi}}");
                        P_acc[c][i][j][m]->GetYaxis()->SetTitleSize(0.08);
                      }
                      P_acc[c][i][j][m]->Draw("SAMEP");
                      P_acc[c][i][j][m]->GetXaxis()->SetLimits(0.1,0.9);
                      P_acc[c][i][j][m]->SetMinimum(0.25);
                      P_acc[c][i][j][m]->SetMaximum(1.04);
                      c6.Range(0.,0.,1.,1.2);
                    }
                    else
                    {
                      P_acc[c][i][j][m]->Draw("SAMEP");
                      P_acc[c][i][j][m]->GetXaxis()->SetLimits(0.1,0.9);
                      P_acc[c][i][j][m]->SetMinimum(0.25);
                      P_acc[c][i][j][m]->SetMaximum(1.04);
                    }
                  }
                  c6.Update();
                }

                if(!k_a_empty)
                {
                  c7.cd(i+1+9*j);
                  if(K_acc[c][i][j][m])
                  {
                    if(!c && !m)
                    {
                      K_acc[c][i][j][m]->Draw("SAMEPA");
                      K_acc[c][i][j][m]->GetXaxis()->SetLimits(0.1,0.9);
                      K_acc[c][i][j][m]->SetMinimum(0.25);
                      K_acc[c][i][j][m]->SetMaximum(1.04);
                      K_acc[c][i][j][m]->GetXaxis()->SetLabelSize(0.06);
                      K_acc[c][i][j][m]->GetYaxis()->SetLabelSize(0.06);
                      K_acc[c][i][j][m]->SetTitle("");
                      if(j==5) gPad->SetBottomMargin(.15);
                      if(i==0) gPad->SetLeftMargin(.22);
                      if(i==8 && j==5)
                      {
                        K_acc[c][i][j][m]->GetXaxis()->SetTitle("#font[ 12]{z}");
                        K_acc[c][i][j][m]->GetXaxis()->SetTitleSize(0.08);
                        K_acc[c][i][j][m]->GetXaxis()->SetTitleOffset(.8);
                      }
                      K_acc[c][i][j][m]->GetXaxis()->SetNdivisions(304,kTRUE);
                      K_acc[c][i][j][m]->GetYaxis()->SetNdivisions(304,kTRUE);
                      if(i==2 && j==0)
                      {
                        K_acc[c][i][j][m]->GetYaxis()->SetTitle("#font[12]{acceptance}^{#font[ 12]{K}}");
                        K_acc[c][i][j][m]->GetYaxis()->SetTitleSize(0.08);
                      }
                      K_acc[c][i][j][m]->Draw("SAMEP");
                      K_acc[c][i][j][m]->GetXaxis()->SetLimits(0.1,0.9);
                      K_acc[c][i][j][m]->SetMinimum(0.25);
                      K_acc[c][i][j][m]->SetMaximum(1.04);
                      c7.Range(0.,0.,1.,1.2);
                    }
                    else
                    {
                      K_acc[c][i][j][m]->Draw("SAMEP");
                      K_acc[c][i][j][m]->GetXaxis()->SetLimits(0.1,0.9);
                      K_acc[c][i][j][m]->SetMinimum(0.25);
                      K_acc[c][i][j][m]->SetMaximum(1.04);
                    }
                  }
                  c7.Update();
                }

                if(!pr_a_empty)
                {
                  c88.cd(i+1+9*j);
                  if(PR_acc[c][i][j][m])
                  {
                    if(!c && !m)
                    {
                      PR_acc[c][i][j][m]->Draw("SAMEPA");
                      PR_acc[c][i][j][m]->GetXaxis()->SetLimits(0.1,0.9);
                      PR_acc[c][i][j][m]->SetMinimum(0.25);
                      PR_acc[c][i][j][m]->SetMaximum(1.04);
                      PR_acc[c][i][j][m]->GetXaxis()->SetLabelSize(0.06);
                      PR_acc[c][i][j][m]->GetYaxis()->SetLabelSize(0.06);
                      PR_acc[c][i][j][m]->SetTitle("");
                      if(j==5) gPad->SetBottomMargin(.15);
                      if(i==0) gPad->SetLeftMargin(.22);
                      if(i==8 && j==5)
                      {
                        PR_acc[c][i][j][m]->GetXaxis()->SetTitle("#font[ 12]{z}");
                        PR_acc[c][i][j][m]->GetXaxis()->SetTitleSize(0.08);
                        PR_acc[c][i][j][m]->GetXaxis()->SetTitleOffset(.8);
                      }
                      PR_acc[c][i][j][m]->GetXaxis()->SetNdivisions(304,kTRUE);
                      PR_acc[c][i][j][m]->GetYaxis()->SetNdivisions(304,kTRUE);
                      if(i==2 && j==0)
                      {
                        PR_acc[c][i][j][m]->GetYaxis()->SetTitle("#font[12]{acceptance}^{#font[ 12]{p}}");
                        PR_acc[c][i][j][m]->GetYaxis()->SetTitleSize(0.08);
                      }
                      PR_acc[c][i][j][m]->Draw("SAMEP");
                      PR_acc[c][i][j][m]->GetXaxis()->SetLimits(0.1,0.9);
                      PR_acc[c][i][j][m]->SetMinimum(0.25);
                      PR_acc[c][i][j][m]->SetMaximum(1.04);
                      c88.Range(0.,0.,1.,1.2);
                    }
                    else
                    {
                      PR_acc[c][i][j][m]->Draw("SAMEP");
                      PR_acc[c][i][j][m]->GetXaxis()->SetLimits(0.1,0.9);
                      PR_acc[c][i][j][m]->SetMinimum(0.25);
                      PR_acc[c][i][j][m]->SetMaximum(1.04);
                    }
                  }
                  c88.Update();
                }
              }
              else
              {
                if(!h_a_empty)
                {
                  c5.cd(i+1);
                  gPad->SetFillStyle(4000);
                  if(H_acc[c][i][j][m])
                  {
                    if(!c && !m && j==3)
                    {
                      H_acc[c][i][j][m]->Draw("SAMEPA");
                      H_acc[c][i][j][m]->GetXaxis()->SetLimits(0.1,0.9);
                      H_acc[c][i][j][m]->SetMinimum(0.);
                      H_acc[c][i][j][m]->SetMaximum(1.2);
                      H_acc[c][i][j][m]->GetXaxis()->SetLabelSize(0.06);
                      H_acc[c][i][j][m]->GetYaxis()->SetLabelSize(0.06);
                      H_acc[c][i][j][m]->SetTitle("");
                      if(j==5) gPad->SetBottomMargin(.15);
                      if(i==0 || i==5) gPad->SetLeftMargin(.22);
                      if(i==8)
                      {
                        H_acc[c][i][j][m]->GetXaxis()->SetTitle("#font[ 12]{z}");
                        H_acc[c][i][j][m]->GetXaxis()->SetTitleSize(0.08);
                        H_acc[c][i][j][m]->GetXaxis()->SetTitleOffset(.8);
                      }
                      H_acc[c][i][j][m]->GetXaxis()->SetNdivisions(304,kTRUE);
                      H_acc[c][i][j][m]->GetYaxis()->SetNdivisions(304,kTRUE);
                      if(i==2 && j==0)
                      {
                        H_acc[c][i][j][m]->GetYaxis()->SetTitle("#font[12]{acceptance}^{#font[ 12]{h}}");
                        H_acc[c][i][j][m]->GetYaxis()->SetTitleSize(0.08);
                      }
                      H_acc[c][i][0][m]->Draw("SAMEP");
                      H_acc[c][i][0][m]->GetXaxis()->SetLimits(0.1,0.9);
                      H_acc[c][i][0][m]->SetMinimum(0.);
                      H_acc[c][i][0][m]->SetMaximum(1.2);
                      H_acc[c][i][1][m]->Draw("SAMEP");
                      H_acc[c][i][1][m]->GetXaxis()->SetLimits(0.1,0.9);
                      H_acc[c][i][1][m]->SetMinimum(0.);
                      H_acc[c][i][1][m]->SetMaximum(1.2);
                      H_acc[c][i][2][m]->Draw("SAMEP");
                      H_acc[c][i][2][m]->GetXaxis()->SetLimits(0.1,0.9);
                      H_acc[c][i][2][m]->SetMinimum(0.);
                      H_acc[c][i][2][m]->SetMaximum(1.2);
                      c5.Range(0.1,0.,0.9,1.2);
                    }
                    else
                    {
                      H_acc[c][i][j][m]->Draw("SAMEP");
                      H_acc[c][i][j][m]->GetXaxis()->SetLimits(0.1,0.9);
                      H_acc[c][i][j][m]->SetMinimum(0.);
                      H_acc[c][i][j][m]->SetMaximum(1.2);
                    }
                  }
                  c5.Update();
                }

                if(!p_a_empty)
                {
                  c6.cd(i+1);
                  if(P_acc[c][i][j][m])
                  {
                    if(!c && !m && j==3)
                    {
                      P_acc[c][i][j][m]->Draw("SAMEPA");
                      P_acc[c][i][j][m]->GetXaxis()->SetLimits(0.1,0.9);
                      P_acc[c][i][j][m]->SetMinimum(0.);
                      P_acc[c][i][j][m]->SetMaximum(1.2);
                      P_acc[c][i][j][m]->GetXaxis()->SetLabelSize(0.06);
                      P_acc[c][i][j][m]->GetYaxis()->SetLabelSize(0.06);
                      P_acc[c][i][j][m]->SetTitle("");
                      if(i>4) gPad->SetBottomMargin(.15);
                      if(i==0 || i==5) gPad->SetLeftMargin(.22);
                      if(i==8)
                      {
                        P_acc[c][i][j][m]->GetXaxis()->SetTitle("#font[ 12]{z}");
                        P_acc[c][i][j][m]->GetXaxis()->SetTitleSize(0.08);
                        P_acc[c][i][j][m]->GetXaxis()->SetTitleOffset(.8);
                      }
                      P_acc[c][i][j][m]->GetXaxis()->SetNdivisions(304,kTRUE);
                      P_acc[c][i][j][m]->GetYaxis()->SetNdivisions(304,kTRUE);
                      if(i==0)
                      {
                        P_acc[c][i][j][m]->GetYaxis()->SetTitle("#font[12]{acceptance}^{#font[ 12]{#pi}}");
                        P_acc[c][i][j][m]->GetYaxis()->SetTitleSize(0.08);
                      }
                      P_acc[c][i][0][m]->Draw("SAMEP");
                      P_acc[c][i][0][m]->GetXaxis()->SetLimits(0.1,0.9);
                      P_acc[c][i][0][m]->SetMinimum(0.);
                      P_acc[c][i][0][m]->SetMaximum(1.2);
                      P_acc[c][i][1][m]->Draw("SAMEP");
                      P_acc[c][i][1][m]->GetXaxis()->SetLimits(0.1,0.9);
                      P_acc[c][i][1][m]->SetMinimum(0.);
                      P_acc[c][i][1][m]->SetMaximum(1.2);
                      P_acc[c][i][2][m]->Draw("SAMEP");
                      P_acc[c][i][2][m]->GetXaxis()->SetLimits(0.1,0.9);
                      P_acc[c][i][2][m]->SetMinimum(0.);
                      P_acc[c][i][2][m]->SetMaximum(1.2);
                      c6.Range(0.,0.,1.,1.2);
                    }
                    else
                    {
                      P_acc[c][i][j][m]->Draw("SAMEP");
                      P_acc[c][i][j][m]->GetXaxis()->SetLimits(0.1,0.9);
                      P_acc[c][i][j][m]->SetMinimum(0.);
                      P_acc[c][i][j][m]->SetMaximum(1.2);
                    }
                  }
                  c6.Update();
                }

                if(!k_a_empty)
                {
                  c7.cd(i+1);
                  if(K_acc[c][i][j][m])
                  {
                    if(!c && !m && j==3)
                    {
                      K_acc[c][i][j][m]->Draw("SAMEPA");
                      K_acc[c][i][j][m]->GetXaxis()->SetLimits(0.1,0.9);
                      K_acc[c][i][j][m]->SetMinimum(0.);
                      K_acc[c][i][j][m]->SetMaximum(1.2);
                      K_acc[c][i][j][m]->GetXaxis()->SetLabelSize(0.06);
                      K_acc[c][i][j][m]->GetYaxis()->SetLabelSize(0.06);
                      K_acc[c][i][j][m]->SetTitle("");
                      if(i>4) gPad->SetBottomMargin(.15);
                      if(i==0 || i==5) gPad->SetLeftMargin(.22);
                      if(i==8)
                      {
                        K_acc[c][i][j][m]->GetXaxis()->SetTitle("#font[ 12]{z}");
                        K_acc[c][i][j][m]->GetXaxis()->SetTitleSize(0.08);
                        K_acc[c][i][j][m]->GetXaxis()->SetTitleOffset(.8);
                      }
                      K_acc[c][i][j][m]->GetXaxis()->SetNdivisions(304,kTRUE);
                      K_acc[c][i][j][m]->GetYaxis()->SetNdivisions(304,kTRUE);
                      if(i==0)
                      {
                        K_acc[c][i][j][m]->GetYaxis()->SetTitle("#font[12]{acceptance}^{#font[ 12]{K}}");
                        K_acc[c][i][j][m]->GetYaxis()->SetTitleSize(0.08);
                      }
                      K_acc[c][i][0][m]->Draw("SAMEP");
                      K_acc[c][i][0][m]->GetXaxis()->SetLimits(0.1,0.9);
                      K_acc[c][i][0][m]->SetMinimum(0.);
                      K_acc[c][i][0][m]->SetMaximum(1.2);
                      K_acc[c][i][1][m]->Draw("SAMEP");
                      K_acc[c][i][1][m]->GetXaxis()->SetLimits(0.1,0.9);
                      K_acc[c][i][1][m]->SetMinimum(0.);
                      K_acc[c][i][1][m]->SetMaximum(1.2);
                      K_acc[c][i][2][m]->Draw("SAMEP");
                      K_acc[c][i][2][m]->GetXaxis()->SetLimits(0.1,0.9);
                      K_acc[c][i][2][m]->SetMinimum(0.);
                      K_acc[c][i][2][m]->SetMaximum(1.2);
                      c7.Range(0.,0.,1.,1.2);
                    }
                    else
                    {
                      K_acc[c][i][j][m]->Draw("SAMEP");
                      K_acc[c][i][j][m]->GetXaxis()->SetLimits(0.1,0.9);
                      K_acc[c][i][j][m]->SetMinimum(0.);
                      K_acc[c][i][j][m]->SetMaximum(1.2);
                    }
                  }
                  c7.Update();
                }

                if(!pr_a_empty)
                {
                  c88.cd(i+1);
                  if(PR_acc[c][i][j][m])
                  {
                    if(!c && !m && j==3)
                    {
                      PR_acc[c][i][j][m]->Draw("SAMEPA");
                      PR_acc[c][i][j][m]->GetXaxis()->SetLimits(0.1,0.9);
                      PR_acc[c][i][j][m]->SetMinimum(0.);
                      PR_acc[c][i][j][m]->SetMaximum(1.2);
                      PR_acc[c][i][j][m]->GetXaxis()->SetLabelSize(0.06);
                      PR_acc[c][i][j][m]->GetYaxis()->SetLabelSize(0.06);
                      PR_acc[c][i][j][m]->SetTitle("");
                      if(i>4) gPad->SetBottomMargin(.15);
                      if(i==0 || i==5) gPad->SetLeftMargin(.22);
                      if(i==8)
                      {
                        PR_acc[c][i][j][m]->GetXaxis()->SetTitle("#font[ 12]{z}");
                        PR_acc[c][i][j][m]->GetXaxis()->SetTitleSize(0.08);
                        PR_acc[c][i][j][m]->GetXaxis()->SetTitleOffset(.8);
                      }
                      PR_acc[c][i][j][m]->GetXaxis()->SetNdivisions(304,kTRUE);
                      PR_acc[c][i][j][m]->GetYaxis()->SetNdivisions(304,kTRUE);
                      if(i==0)
                      {
                        PR_acc[c][i][j][m]->GetYaxis()->SetTitle("#font[12]{acceptance}^{#font[ 12]{p}}");
                        PR_acc[c][i][j][m]->GetYaxis()->SetTitleSize(0.08);
                      }
                      PR_acc[c][i][0][m]->Draw("SAMEP");
                      PR_acc[c][i][0][m]->GetXaxis()->SetLimits(0.1,0.9);
                      PR_acc[c][i][0][m]->SetMinimum(0.);
                      PR_acc[c][i][0][m]->SetMaximum(1.2);
                      PR_acc[c][i][1][m]->Draw("SAMEP");
                      PR_acc[c][i][1][m]->GetXaxis()->SetLimits(0.1,0.9);
                      PR_acc[c][i][1][m]->SetMinimum(0.);
                      PR_acc[c][i][1][m]->SetMaximum(1.2);
                      PR_acc[c][i][2][m]->Draw("SAMEP");
                      PR_acc[c][i][2][m]->GetXaxis()->SetLimits(0.1,0.9);
                      PR_acc[c][i][2][m]->SetMinimum(0.);
                      PR_acc[c][i][2][m]->SetMaximum(1.2);
                      c88.Range(0.,0.,1.,1.2);
                    }
                    else
                    {
                      PR_acc[c][i][j][m]->Draw("SAMEP");
                      PR_acc[c][i][j][m]->GetXaxis()->SetLimits(0.1,0.9);
                      PR_acc[c][i][j][m]->SetMinimum(0.);
                      PR_acc[c][i][j][m]->SetMaximum(1.2);
                    }
                  }
                  c88.Update();
                }
              }
            }

            for(int k=0; k<12; k++)
            {
              yavg(c,i,k);

              fAcceptance_yavg[i][k].tab[c][m][0][0] = ((fNDIS_evt_yavg[m][0] && fNDIS_evt_MC_yavg[m][0] && fGnrt_yavg[m][0]) ? double((fRcstr_yavg[m][0]/fNDIS_evt_yavg[m][0])/(fGnrt_yavg[m][0]/fNDIS_evt_MC_yavg[m][0])) : 0);
              fAcceptance_yavg[i][k].tab[c][m][0][1] = ((fNDIS_evt_yavg[m][1] && fNDIS_evt_MC_yavg[m][1] && fGnrt_yavg[m][1]) ? double((fRcstr_yavg[m][1]/fNDIS_evt_yavg[m][1])/(fGnrt_yavg[m][1]/fNDIS_evt_MC_yavg[m][1])) : 0);
              fAcceptance_yavg[i][k].tab[c][m][0][2] = ((fNDIS_evt_yavg[m][2] && fNDIS_evt_MC_yavg[m][2] && fGnrt_yavg[m][2]) ? double((fRcstr_yavg[m][2]/fNDIS_evt_yavg[m][2])/(fGnrt_yavg[m][2]/fNDIS_evt_MC_yavg[m][2])) : 0);
              fAcceptance_yavg[i][k].tab[c][m][0][3] = ((fNDIS_evt_yavg[m][0] && fNDIS_evt_MC_yavg[m][0] && fGnrt_yavg[m][3]) ? double((fRcstr_yavg[m][3]/fNDIS_evt_yavg[m][0])/(fGnrt_yavg[m][3]/fNDIS_evt_MC_yavg[m][0])) : 0);

              if(fAcceptance_yavg[i][k].tab[c][m][0][0]<0) fAcceptance_yavg[i][k].tab[c][m][0][0]=0;
              if(fAcceptance_yavg[i][k].tab[c][m][0][1]<0) fAcceptance_yavg[i][k].tab[c][m][0][1]=0;
              if(fAcceptance_yavg[i][k].tab[c][m][0][2]<0) fAcceptance_yavg[i][k].tab[c][m][0][2]=0;
              if(fAcceptance_yavg[i][k].tab[c][m][0][3]<0) fAcceptance_yavg[i][k].tab[c][m][0][3]=0;

              double rd_prime[3];
              double rd[3];
              double gd[3];
              double rh_prime[4];
              double rh[4];
              double gh[4];

              rd_prime[0] = fNDIS_evt_yavg[m][0] - fNDIS_evt_c_yavg[m][0];
              rd_prime[1] = fNDIS_evt_yavg[m][1] - fNDIS_evt_c_yavg[m][1];
              rd_prime[2] = fNDIS_evt_yavg[m][2] - fNDIS_evt_c_yavg[m][2];

              rd[0] = fNDIS_evt_c_yavg[m][0];
              rd[1] = fNDIS_evt_c_yavg[m][1];
              rd[2] = fNDIS_evt_c_yavg[m][2];

              gd[0] = fNDIS_evt_MC_yavg[m][0];
              gd[1] = fNDIS_evt_MC_yavg[m][1];
              gd[2] = fNDIS_evt_MC_yavg[m][2];

              rh_prime[0] = fRcstr_yavg[m][0] - fRcstr_c_yavg[m][0];
              rh_prime[1] = fRcstr_yavg[m][1] - fRcstr_c_yavg[m][1];
              rh_prime[2] = fRcstr_yavg[m][2] - fRcstr_c_yavg[m][2];
              rh_prime[3] = fRcstr_yavg[m][3] - fRcstr_c_yavg[m][3];

              rh[0] = fRcstr_c_yavg[m][0];
              rh[1] = fRcstr_c_yavg[m][1];
              rh[2] = fRcstr_c_yavg[m][2];
              rh[3] = fRcstr_c_yavg[m][3];

              gh[0] = fGnrt_yavg[m][0];
              gh[1] = fGnrt_yavg[m][1];
              gh[2] = fGnrt_yavg[m][2];
              gh[3] = fGnrt_yavg[m][3];

              fAcceptance_yavg[i][k].tab[c][m][1][0] = sqrt(((rh[0]>0 && gh[0]>0 && gd[0]>0 && ((rd[0]+rd_prime[0])!=0)) ? double(pow(gd[0]/(rd_prime[0]+rd[0]),2)
                                                                                  *((rh[0]+1)*(gh[0]-rh[0]+1)/(pow(gh[0]+2,2)*(gh[0]+3))+rh_prime[0]/pow(gh[0],2)+pow(rh_prime[0],2)/pow(gh[0],3))
                                                                                  +pow(gd[0]/(rd[0]+rd_prime[0]),4)*pow((rh[0]+rh_prime[0])/gh[0],2)*((rd[0]+1)*(gd[0]-rd[0]+1)/(pow(gd[0]+2,2)*(gd[0]+3))+rd_prime[0]/pow(gd[0],2)+pow(rd_prime[0],2)/pow(gd[0],3))) : 0));
              fAcceptance_yavg[i][k].tab[c][m][1][1] = sqrt(((rh[1]>0 && gh[1]>0 && gd[1]>0 && ((rd[1]+rd_prime[1])!=0)) ? double(pow(gd[1]/(rd_prime[1]+rd[1]),2)
                                                                                  *((rh[1]+1)*(gh[1]-rh[1]+1)/(pow(gh[1]+2,2)*(gh[1]+3))+rh_prime[1]/pow(gh[1],2)+pow(rh_prime[1],2)/pow(gh[1],3))
                                                                                  +pow(gd[1]/(rd[1]+rd_prime[1]),4)*pow((rh[1]+rh_prime[1])/gh[1],2)*((rd[1]+1)*(gd[1]-rd[1]+1)/(pow(gd[1]+2,2)*(gd[1]+3))+rd_prime[1]/pow(gd[1],2)+pow(rd_prime[1],2)/pow(gd[1],3))) : 0));
              fAcceptance_yavg[i][k].tab[c][m][1][2] = sqrt(((rh[2]>0 && gh[2]>0 && gd[2]>0 && ((rd[2]+rd_prime[2])!=0)) ? double(pow(gd[2]/(rd_prime[2]+rd[2]),2)
                                                                                  *((rh[2]+1)*(gh[2]-rh[2]+1)/(pow(gh[2]+2,2)*(gh[2]+3))+rh_prime[2]/pow(gh[2],2)+pow(rh_prime[2],2)/pow(gh[2],3))
                                                                                  +pow(gd[2]/(rd[2]+rd_prime[2]),4)*pow((rh[2]+rh_prime[2])/gh[2],2)*((rd[2]+1)*(gd[2]-rd[2]+1)/(pow(gd[2]+2,2)*(gd[2]+3))+rd_prime[2]/pow(gd[2],2)+pow(rd_prime[2],2)/pow(gd[2],3))) : 0));
              fAcceptance_yavg[i][k].tab[c][m][1][3] = sqrt(((rh[3]>0 && gh[3]>0 && gd[0]>0 && ((rd[0]+rd_prime[0])!=0)) ? double(pow(gd[0]/(rd_prime[0]+rd[0]),2)
                                                                                  *((rh[3]+1)*(gh[3]-rh[3]+1)/(pow(gh[3]+2,2)*(gh[3]+3))+rh_prime[3]/pow(gh[3],2)+pow(rh_prime[3],2)/pow(gh[3],3))
                                                                                  +pow(gd[0]/(rd[0]+rd_prime[0]),4)*pow((rh[3]+rh_prime[3])/gh[3],2)*((rd[0]+1)*(gd[0]-rd[0]+1)/(pow(gd[0]+2,2)*(gd[0]+3))+rd_prime[0]/pow(gd[0],2)+pow(rd_prime[0],2)/pow(gd[0],3))) : 0));

              if(fAcceptance_yavg[i][k].tab[c][m][0][0]==0) fAcceptance_yavg[i][k].tab[c][m][1][0]=0;
              if(fAcceptance_yavg[i][k].tab[c][m][0][1]==0) fAcceptance_yavg[i][k].tab[c][m][1][1]=0;
              if(fAcceptance_yavg[i][k].tab[c][m][0][2]==0) fAcceptance_yavg[i][k].tab[c][m][1][2]=0;
              if(fAcceptance_yavg[i][k].tab[c][m][0][3]==0) fAcceptance_yavg[i][k].tab[c][m][1][3]=0;

              p_y.push_back(fAcceptance_yavg[i][k].tab[c][m][0][0]);
              k_y.push_back(fAcceptance_yavg[i][k].tab[c][m][0][1]);
              h_y.push_back(fAcceptance_yavg[i][k].tab[c][m][0][3]);

              p_y_err.push_back(fAcceptance_yavg[i][k].tab[c][m][1][0]);
              k_y_err.push_back(fAcceptance_yavg[i][k].tab[c][m][1][1]);
              h_y_err.push_back(fAcceptance_yavg[i][k].tab[c][m][1][3]);
            }

            for(int k=12; k>0; k--)
            {
              if(!p_y[k-1]) {p_y.erase(p_y.begin()+k-1); p_y_err.erase(p_y_err.begin()+k-1); z_range_p_y.erase(z_range_p_y.begin()+k-1);}
              if(!k_y[k-1]) {k_y.erase(k_y.begin()+k-1); k_y_err.erase(k_y_err.begin()+k-1); z_range_k_y.erase(z_range_k_y.begin()+k-1);}
              if(!h_y[k-1]) {h_y.erase(h_y.begin()+k-1); h_y_err.erase(h_y_err.begin()+k-1); z_range_h_y.erase(z_range_h_y.begin()+k-1);}
            }

            bool p_y_empty = 0;
            bool k_y_empty = 0;
            bool h_y_empty = 0;

            if(!(int(p_y.size()))) p_y_empty = 1;
            if(!(int(k_y.size()))) k_y_empty = 1;
            if(!(int(h_y.size()))) h_y_empty = 1;

            H_y[c][i][m] = new TGraphErrors(int(h_y.size()),&(z_range_h_y[0]),&(h_y[0]),0,&(h_y_err[0]));
            P_y[c][i][m] = new TGraphErrors(int(p_y.size()),&(z_range_p_y[0]),&(p_y[0]),0,&(p_y_err[0]));
            K_y[c][i][m] = new TGraphErrors(int(k_y.size()),&(z_range_k_y[0]),&(k_y[0]),0,&(k_y_err[0]));

            if(!c)
            {
              H_y[c][i][m]->SetMarkerColor(fMarkerColor[4]);
              P_y[c][i][m]->SetMarkerColor(fMarkerColor[4]);
              K_y[c][i][m]->SetMarkerColor(fMarkerColor[4]);
            }
            else
            {
              H_y[c][i][m]->SetMarkerColor(fMarkerColor[0]);
              P_y[c][i][m]->SetMarkerColor(fMarkerColor[0]);
              K_y[c][i][m]->SetMarkerColor(fMarkerColor[0]);
            }

            H_y[c][i][m]->SetMarkerSize(3);
            P_y[c][i][m]->SetMarkerSize(3);
            K_y[c][i][m]->SetMarkerSize(3);

            if(!m)
            {
              H_y[c][i][m]->SetMarkerStyle(fMarkerStyleb[1][c]);
              P_y[c][i][m]->SetMarkerStyle(fMarkerStyleb[1][c]);
              K_y[c][i][m]->SetMarkerStyle(fMarkerStyleb[1][c]);
            }
            else
            {
              H_y[c][i][m]->SetMarkerStyle(fMarkerStyleb[0][c]);
              P_y[c][i][m]->SetMarkerStyle(fMarkerStyleb[0][c]);
              K_y[c][i][m]->SetMarkerStyle(fMarkerStyleb[0][c]);
            }

            H_y[c][i][m]->GetYaxis()->SetTitle("");
            P_y[c][i][m]->GetYaxis()->SetTitle("");
            K_y[c][i][m]->GetYaxis()->SetTitle("");

            H_y[c][i][m]->GetXaxis()->SetTitle("");
            P_y[c][i][m]->GetXaxis()->SetTitle("");
            K_y[c][i][m]->GetXaxis()->SetTitle("");

            H_y[c][i][m]->SetTitle("");
            P_y[c][i][m]->SetTitle("");
            K_y[c][i][m]->SetTitle("");

            if(!h_y_empty)
            {
              c9.cd(i+1);
              gPad->SetFillStyle(4000);
              if(H_y[c][i][m])
              {
                if(!c && !m)
                {
                  H_y[c][i][m]->Draw("SAMEPA");
                  H_y[c][i][m]->GetXaxis()->SetLimits(0.1,0.9);
                  H_y[c][i][m]->SetMinimum(0.);
                  H_y[c][i][m]->SetMaximum(1.2);
                  H_y[c][i][m]->GetXaxis()->SetLabelSize(0.06);
                  H_y[c][i][m]->GetYaxis()->SetLabelSize(0.06);
                  H_y[c][i][m]->SetTitle("");
                  if(i>4) gPad->SetBottomMargin(.15);
                  if(i==0 || i==5) gPad->SetLeftMargin(.22);
                  if(i==8)
                  {
                    H_y[c][i][m]->GetXaxis()->SetTitle("#font[ 12]{z}");
                    H_y[c][i][m]->GetXaxis()->SetTitleSize(0.08);
                    H_y[c][i][m]->GetXaxis()->SetTitleOffset(.8);
                  }
                  H_y[c][i][m]->GetXaxis()->SetNdivisions(304,kTRUE);
                  H_y[c][i][m]->GetYaxis()->SetNdivisions(304,kTRUE);
                  if(i==0)
                  {
                    H_y[c][i][m]->GetYaxis()->SetTitle("#font[12]{acceptance}^{#font[ 12]{h}}");
                    H_y[c][i][m]->GetYaxis()->SetTitleSize(0.08);
                  }
                  H_y[c][i][m]->Draw("SAMEP");
                  H_y[c][i][m]->GetXaxis()->SetLimits(0.1,0.9);
                  H_y[c][i][m]->SetMinimum(0.);
                  H_y[c][i][m]->SetMaximum(1.2);
                  c9.Range(0.1,0.,0.9,1.2);
                }
                else
                {
                  H_y[c][i][m]->Draw("SAMEP");
                  H_y[c][i][m]->GetXaxis()->SetLimits(0.1,0.9);
                  H_y[c][i][m]->SetMinimum(0.);
                  H_y[c][i][m]->SetMaximum(1.2);
                }
              }
              c9.Update();
            }

            if(!p_y_empty)
            {
              c10.cd(i+1);
              if(P_y[c][i][m])
              {
                if(!c && !m)
                {
                  P_y[c][i][m]->Draw("SAMEPA");
                  P_y[c][i][m]->GetXaxis()->SetLimits(0.1,0.9);
                  P_y[c][i][m]->SetMinimum(0.);
                  P_y[c][i][m]->SetMaximum(1.2);
                  P_y[c][i][m]->GetXaxis()->SetLabelSize(0.06);
                  P_y[c][i][m]->GetYaxis()->SetLabelSize(0.06);
                  P_y[c][i][m]->SetTitle("");
                  if(i>4) gPad->SetBottomMargin(.15);
                  if(i==0 || i==5) gPad->SetLeftMargin(.22);
                  if(i==8)
                  {
                    P_y[c][i][m]->GetXaxis()->SetTitle("#font[ 12]{z}");
                    P_y[c][i][m]->GetXaxis()->SetTitleSize(0.08);
                    P_y[c][i][m]->GetXaxis()->SetTitleOffset(.8);
                  }
                  P_y[c][i][m]->GetXaxis()->SetNdivisions(304,kTRUE);
                  P_y[c][i][m]->GetYaxis()->SetNdivisions(304,kTRUE);
                  if(i==0)
                  {
                    P_y[c][i][m]->GetYaxis()->SetTitle("#font[12]{acceptance}^{#font[ 12]{#pi}}");
                    P_y[c][i][m]->GetYaxis()->SetTitleSize(0.08);
                  }
                  P_y[c][i][m]->Draw("SAMEP");
                  P_y[c][i][m]->GetXaxis()->SetLimits(0.1,0.9);
                  P_y[c][i][m]->SetMinimum(0.);
                  P_y[c][i][m]->SetMaximum(1.2);
                  c10.Range(0.,0.,1.,1.2);
                }
                else
                {
                  P_y[c][i][m]->Draw("SAMEP");
                  P_y[c][i][m]->GetXaxis()->SetLimits(0.1,0.9);
                  P_y[c][i][m]->SetMinimum(0.);
                  P_y[c][i][m]->SetMaximum(1.2);
                }
              }
              c10.Update();
            }

            if(!k_y_empty)
            {
              c11.cd(i+1);
              if(K_y[c][i][m])
              {
                if(!c && !m)
                {
                  K_y[c][i][m]->Draw("SAMEPA");
                  K_y[c][i][m]->GetXaxis()->SetLimits(0.1,0.9);
                  K_y[c][i][m]->SetMinimum(0.);
                  K_y[c][i][m]->SetMaximum(1.2);
                  K_y[c][i][m]->GetXaxis()->SetLabelSize(0.06);
                  K_y[c][i][m]->GetYaxis()->SetLabelSize(0.06);
                  K_y[c][i][m]->SetTitle("");
                  if(i>4) gPad->SetBottomMargin(.15);
                  if(i==0 || i==5) gPad->SetLeftMargin(.22);
                  if(i==8)
                  {
                    K_y[c][i][m]->GetXaxis()->SetTitle("#font[ 12]{z}");
                    K_y[c][i][m]->GetXaxis()->SetTitleSize(0.08);
                    K_y[c][i][m]->GetXaxis()->SetTitleOffset(.8);
                  }
                  K_y[c][i][m]->GetXaxis()->SetNdivisions(304,kTRUE);
                  K_y[c][i][m]->GetYaxis()->SetNdivisions(304,kTRUE);
                  if(i==0)
                  {
                    K_y[c][i][m]->GetYaxis()->SetTitle("#font[12]{acceptance}^{#font[ 12]{K}}");
                    K_y[c][i][m]->GetYaxis()->SetTitleSize(0.08);
                  }
                  K_y[c][i][m]->Draw("SAMEP");
                  K_y[c][i][m]->GetXaxis()->SetLimits(0.1,0.9);
                  K_y[c][i][m]->SetMinimum(0.);
                  K_y[c][i][m]->SetMaximum(1.2);
                  c11.Range(0.,0.,1.,1.2);
                }
                else
                {
                  K_y[c][i][m]->Draw("SAMEP");
                  K_y[c][i][m]->GetXaxis()->SetLimits(0.1,0.9);
                  K_y[c][i][m]->SetMinimum(0.);
                  K_y[c][i][m]->SetMaximum(1.2);
                }
              }
              c11.Update();
            }
          }
        }
      }

    TLatex fTitle;

    if(SPREAD)
    {

    }
    else
    {
      c5.cd(1);
      fTitle.SetTextSize(0.078);
      fTitle.SetTextAlign(21);
      fTitle.DrawLatex(0.5, 1.0,"0.004#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.01");

      c5.cd(2);
      fTitle.SetTextSize(0.078);
      fTitle.SetTextAlign(21);
      fTitle.DrawLatex(0.5, 1.0,"0.01#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.02");

      c5.cd(3);
      fTitle.SetTextSize(0.078);
      fTitle.SetTextAlign(21);
      fTitle.DrawLatex(0.5, 1.0,"0.02#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.03");

      c5.cd(4);
      fTitle.SetTextSize(0.078);
      fTitle.SetTextAlign(21);
      fTitle.DrawLatex(0.5, 1.0,"0.03#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.04");

      c5.cd(5);
      fTitle.SetTextSize(0.078);
      fTitle.SetTextAlign(21);
      fTitle.DrawLatex(0.5, 1.0,"0.04#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.06");

      c5.cd(6);
      fTitle.SetTextSize(0.078);
      fTitle.SetTextAlign(21);
      fTitle.DrawLatex(0.5, 1.0,"0.06#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.1");

      c5.cd(7);
      fTitle.SetTextSize(0.078);
      fTitle.SetTextAlign(21);
      fTitle.DrawLatex(0.5, 1.0,"0.1#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.14");

      c5.cd(8);
      fTitle.SetTextSize(0.078);
      fTitle.SetTextAlign(21);
      fTitle.DrawLatex(0.5, 1.0,"0.14#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.18");

      c5.cd(9);
      fTitle.SetTextSize(0.078);
      fTitle.SetTextAlign(21);
      fTitle.DrawLatex(0.5, 1.0,"0.18#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.4");

      c5.cd(10);
      fTitle.SetTextSize(0.095);
      fTitle.SetTextAlign(11);
      fTitle.DrawLatex(0.05, 0.72,"#color[221]{0.70#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.90}");
      fTitle.DrawLatex(0.05, 0.64,"#color[4]{0.50#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.70}");
      fTitle.DrawLatex(0.05, 0.56,"#color[226]{0.30#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.50}");
      fTitle.DrawLatex(0.05, 0.48,"#color[209]{0.20#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.30}");
      fTitle.DrawLatex(0.05, 0.40,"#color[95]{0.15#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.20}");
      fTitle.DrawLatex(0.05, 0.32,"#color[2]{0.10#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.15}");


      c6.cd(1);
      fTitle.SetTextSize(0.078);
      fTitle.SetTextAlign(21);
      fTitle.DrawLatex(0.5, 1.0,"0.004#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.01");

      c6.cd(2);
      fTitle.SetTextSize(0.078);
      fTitle.SetTextAlign(21);
      fTitle.DrawLatex(0.5, 1.0,"0.01#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.02");

      c6.cd(3);
      fTitle.SetTextSize(0.078);
      fTitle.SetTextAlign(21);
      fTitle.DrawLatex(0.5, 1.0,"0.02#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.03");

      c6.cd(4);
      fTitle.SetTextSize(0.078);
      fTitle.SetTextAlign(21);
      fTitle.DrawLatex(0.5, 1.0,"0.03#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.04");

      c6.cd(5);
      fTitle.SetTextSize(0.078);
      fTitle.SetTextAlign(21);
      fTitle.DrawLatex(0.5, 1.0,"0.04#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.06");

      c6.cd(6);
      fTitle.SetTextSize(0.078);
      fTitle.SetTextAlign(21);
      fTitle.DrawLatex(0.5, 1.0,"0.06#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.1");

      c6.cd(7);
      fTitle.SetTextSize(0.078);
      fTitle.SetTextAlign(21);
      fTitle.DrawLatex(0.5, 1.0,"0.1#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.14");

      c6.cd(8);
      fTitle.SetTextSize(0.078);
      fTitle.SetTextAlign(21);
      fTitle.DrawLatex(0.5, 1.0,"0.14#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.18");

      c6.cd(9);
      fTitle.SetTextSize(0.078);
      fTitle.SetTextAlign(21);
      fTitle.DrawLatex(0.5, 1.0,"0.18#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.4");

      c6.cd(10);
      fTitle.SetTextSize(0.095);
      fTitle.SetTextAlign(11);
      fTitle.DrawLatex(0.05, 0.72,"#color[221]{0.70#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.90}");
      fTitle.DrawLatex(0.05, 0.64,"#color[4]{0.50#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.70}");
      fTitle.DrawLatex(0.05, 0.56,"#color[226]{0.30#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.50}");
      fTitle.DrawLatex(0.05, 0.48,"#color[209]{0.20#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.30}");
      fTitle.DrawLatex(0.05, 0.40,"#color[95]{0.15#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.20}");
      fTitle.DrawLatex(0.05, 0.32,"#color[2]{0.10#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.15}");

      c7.cd(1);
      fTitle.SetTextSize(0.078);
      fTitle.SetTextAlign(21);
      fTitle.DrawLatex(0.5, 1.0,"0.004#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.01");

      c7.cd(2);
      fTitle.SetTextSize(0.078);
      fTitle.SetTextAlign(21);
      fTitle.DrawLatex(0.5, 1.0,"0.01#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.02");

      c7.cd(3);
      fTitle.SetTextSize(0.078);
      fTitle.SetTextAlign(21);
      fTitle.DrawLatex(0.5, 1.0,"0.02#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.03");

      c7.cd(4);
      fTitle.SetTextSize(0.078);
      fTitle.SetTextAlign(21);
      fTitle.DrawLatex(0.5, 1.0,"0.03#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.04");

      c7.cd(5);
      fTitle.SetTextSize(0.078);
      fTitle.SetTextAlign(21);
      fTitle.DrawLatex(0.5, 1.0,"0.04#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.06");

      c7.cd(6);
      fTitle.SetTextSize(0.078);
      fTitle.SetTextAlign(21);
      fTitle.DrawLatex(0.5, 1.0,"0.06#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.1");

      c7.cd(7);
      fTitle.SetTextSize(0.078);
      fTitle.SetTextAlign(21);
      fTitle.DrawLatex(0.5, 1.0,"0.1#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.14");

      c7.cd(8);
      fTitle.SetTextSize(0.078);
      fTitle.SetTextAlign(21);
      fTitle.DrawLatex(0.5, 1.0,"0.14#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.18");

      c7.cd(9);
      fTitle.SetTextSize(0.078);
      fTitle.SetTextAlign(21);
      fTitle.DrawLatex(0.5, 1.0,"0.18#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.4");

      c7.cd(10);
      fTitle.SetTextSize(0.095);
      fTitle.SetTextAlign(11);
      fTitle.DrawLatex(0.05, 0.72,"#color[221]{0.70#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.90}");
      fTitle.DrawLatex(0.05, 0.64,"#color[4]{0.50#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.70}");
      fTitle.DrawLatex(0.05, 0.56,"#color[226]{0.30#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.50}");
      fTitle.DrawLatex(0.05, 0.48,"#color[209]{0.20#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.30}");
      fTitle.DrawLatex(0.05, 0.40,"#color[95]{0.15#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.20}");
      fTitle.DrawLatex(0.05, 0.32,"#color[2]{0.10#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.15}");
    }

    for(int i=0; i<12; i++)
    {
      c8[i]->cd(1);
      fTitle.SetTextSize(0.078);
      fTitle.SetTextAlign(21);
      fTitle.DrawLatex(-200, 1.0,"0.004#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.01");

      c8[i]->cd(2);
      fTitle.SetTextSize(0.078);
      fTitle.SetTextAlign(21);
      fTitle.DrawLatex(-200, 1.0,"0.01#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.02");

      c8[i]->cd(3);
      fTitle.SetTextSize(0.078);
      fTitle.SetTextAlign(21);
      fTitle.DrawLatex(-200, 1.0,"0.02#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.03");

      c8[i]->cd(4);
      fTitle.SetTextSize(0.078);
      fTitle.SetTextAlign(21);
      fTitle.DrawLatex(-200, 1.0,"0.03#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.04");

      c8[i]->cd(5);
      fTitle.SetTextSize(0.078);
      fTitle.SetTextAlign(21);
      fTitle.DrawLatex(-200, 1.0,"0.04#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.06");

      c8[i]->cd(6);
      fTitle.SetTextSize(0.078);
      fTitle.SetTextAlign(21);
      fTitle.DrawLatex(-200, 1.0,"0.06#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.1");

      c8[i]->cd(7);
      fTitle.SetTextSize(0.078);
      fTitle.SetTextAlign(21);
      fTitle.DrawLatex(-200, 1.0,"0.1#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.14");

      c8[i]->cd(8);
      fTitle.SetTextSize(0.078);
      fTitle.SetTextAlign(21);
      fTitle.DrawLatex(-200, 1.0,"0.14#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.18");

      c8[i]->cd(9);
      fTitle.SetTextSize(0.078);
      fTitle.SetTextAlign(21);
      fTitle.DrawLatex(-200, 1.0,"0.18#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.4");

      c8[i]->cd(10);
      fTitle.SetTextSize(0.095);
      fTitle.SetTextAlign(11);
      fTitle.DrawLatex(0.05, 0.72,"#color[221]{0.70#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.90}");
      fTitle.DrawLatex(0.05, 0.64,"#color[4]{0.50#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.70}");
      fTitle.DrawLatex(0.05, 0.56,"#color[226]{0.30#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.50}");
      fTitle.DrawLatex(0.05, 0.48,"#color[209]{0.20#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.30}");
      fTitle.DrawLatex(0.05, 0.40,"#color[95]{0.15#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.20}");
      fTitle.DrawLatex(0.05, 0.32,"#color[2]{0.10#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.15}");

      c8[i]->Update();
    }

    c9.cd(1);
    fTitle.SetTextSize(0.078);
    fTitle.SetTextAlign(21);
    fTitle.DrawLatex(0.5, 1.0,"0.004#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.01");

    c9.cd(2);
    fTitle.SetTextSize(0.078);
    fTitle.SetTextAlign(21);
    fTitle.DrawLatex(0.5, 1.0,"0.01#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.02");

    c9.cd(3);
    fTitle.SetTextSize(0.078);
    fTitle.SetTextAlign(21);
    fTitle.DrawLatex(0.5, 1.0,"0.02#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.03");

    c9.cd(4);
    fTitle.SetTextSize(0.078);
    fTitle.SetTextAlign(21);
    fTitle.DrawLatex(0.5, 1.0,"0.03#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.04");

    c9.cd(5);
    fTitle.SetTextSize(0.078);
    fTitle.SetTextAlign(21);
    fTitle.DrawLatex(0.5, 1.0,"0.04#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.06");

    c9.cd(6);
    fTitle.SetTextSize(0.078);
    fTitle.SetTextAlign(21);
    fTitle.DrawLatex(0.5, 1.0,"0.06#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.1");

    c9.cd(7);
    fTitle.SetTextSize(0.078);
    fTitle.SetTextAlign(21);
    fTitle.DrawLatex(0.5, 1.0,"0.1#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.14");

    c9.cd(8);
    fTitle.SetTextSize(0.078);
    fTitle.SetTextAlign(21);
    fTitle.DrawLatex(0.5, 1.0,"0.14#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.18");

    c9.cd(9);
    fTitle.SetTextSize(0.078);
    fTitle.SetTextAlign(21);
    fTitle.DrawLatex(0.5, 1.0,"0.18#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.4");

    c10.cd(1);
    fTitle.SetTextSize(0.078);
    fTitle.SetTextAlign(21);
    fTitle.DrawLatex(0.5, 1.0,"0.004#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.01");

    c10.cd(2);
    fTitle.SetTextSize(0.078);
    fTitle.SetTextAlign(21);
    fTitle.DrawLatex(0.5, 1.0,"0.01#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.02");

    c10.cd(3);
    fTitle.SetTextSize(0.078);
    fTitle.SetTextAlign(21);
    fTitle.DrawLatex(0.5, 1.0,"0.02#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.03");

    c10.cd(4);
    fTitle.SetTextSize(0.078);
    fTitle.SetTextAlign(21);
    fTitle.DrawLatex(0.5, 1.0,"0.03#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.04");

    c10.cd(5);
    fTitle.SetTextSize(0.078);
    fTitle.SetTextAlign(21);
    fTitle.DrawLatex(0.5, 1.0,"0.04#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.06");

    c10.cd(6);
    fTitle.SetTextSize(0.078);
    fTitle.SetTextAlign(21);
    fTitle.DrawLatex(0.5, 1.0,"0.06#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.1");

    c10.cd(7);
    fTitle.SetTextSize(0.078);
    fTitle.SetTextAlign(21);
    fTitle.DrawLatex(0.5, 1.0,"0.1#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.14");

    c10.cd(8);
    fTitle.SetTextSize(0.078);
    fTitle.SetTextAlign(21);
    fTitle.DrawLatex(0.5, 1.0,"0.14#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.18");

    c10.cd(9);
    fTitle.SetTextSize(0.078);
    fTitle.SetTextAlign(21);
    fTitle.DrawLatex(0.5, 1.0,"0.18#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.4");

    c11.cd(1);
    fTitle.SetTextSize(0.078);
    fTitle.SetTextAlign(21);
    fTitle.DrawLatex(0.5, 1.0,"0.004#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.01");

    c11.cd(2);
    fTitle.SetTextSize(0.078);
    fTitle.SetTextAlign(21);
    fTitle.DrawLatex(0.5, 1.0,"0.01#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.02");

    c11.cd(3);
    fTitle.SetTextSize(0.078);
    fTitle.SetTextAlign(21);
    fTitle.DrawLatex(0.5, 1.0,"0.02#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.03");

    c11.cd(4);
    fTitle.SetTextSize(0.078);
    fTitle.SetTextAlign(21);
    fTitle.DrawLatex(0.5, 1.0,"0.03#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.04");

    c11.cd(5);
    fTitle.SetTextSize(0.078);
    fTitle.SetTextAlign(21);
    fTitle.DrawLatex(0.5, 1.0,"0.04#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.06");

    c11.cd(6);
    fTitle.SetTextSize(0.078);
    fTitle.SetTextAlign(21);
    fTitle.DrawLatex(0.5, 1.0,"0.06#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.1");

    c11.cd(7);
    fTitle.SetTextSize(0.078);
    fTitle.SetTextAlign(21);
    fTitle.DrawLatex(0.5, 1.0,"0.1#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.14");

    c11.cd(8);
    fTitle.SetTextSize(0.078);
    fTitle.SetTextAlign(21);
    fTitle.DrawLatex(0.5, 1.0,"0.14#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.18");

    c11.cd(9);
    fTitle.SetTextSize(0.078);
    fTitle.SetTextAlign(21);
    fTitle.DrawLatex(0.5, 1.0,"0.18#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.4");


    c5.Update();
    c6.Update();
    c7.Update();
    for(int i=0; i<12; i++) c8[i]->Update();
    c9.Update();
    c10.Update();
    c11.Update();

    //Output file
    //q_bin x_bin y_bin z_bin acc_pi acc_error_pi acc_k acc_error_k acc_p acc_error_p acc_h acc_error_h

    for(int c=0; c<2; c++)
    {
      for(int i=0; i<9; i++)
      {
        for(int j=0; j<6; j++)
        {
          for(int k=0; k<12; k++)
          {
            ofs << c << " " << fXrange[i] << " " << fYrange[j] << " " << fZrange[k] << " " <<
            fAcceptance[i][j][k].tab[c][1][0][0] << " " <<
            fAcceptance[i][j][k].tab[c][1][1][0] << " " <<
            fAcceptance[i][j][k].tab[c][1][0][1] << " " <<
            fAcceptance[i][j][k].tab[c][1][1][1] << " " <<
            fAcceptance[i][j][k].tab[c][1][0][2] << " " <<
            fAcceptance[i][j][k].tab[c][1][1][2] << " " <<
            fAcceptance[i][j][k].tab[c][1][0][3] << " " <<
            fAcceptance[i][j][k].tab[c][1][1][3] << " " <<
            fAcceptance[i][j][k].tab[c][0][0][0] << " " <<
            fAcceptance[i][j][k].tab[c][0][1][0] << " " <<
            fAcceptance[i][j][k].tab[c][0][0][1] << " " <<
            fAcceptance[i][j][k].tab[c][0][1][1] << " " <<
            fAcceptance[i][j][k].tab[c][0][0][2] << " " <<
            fAcceptance[i][j][k].tab[c][0][1][2] << " " <<
            fAcceptance[i][j][k].tab[c][0][0][3] << " " <<
            fAcceptance[i][j][k].tab[c][0][1][3] << endl;

            lepto << c << " " << fXrange[i] << " " << fYrange[j] << " " << fZrange[k] << " " <<
            (fAcceptance[i][j][k].tab[c][1][0][0] ? Double_t((fGnrt[i][j][k].tab[c][1][0][0]/pow(fNDIS_evt_MC[0][1][i][j][k],2)-pow(fGnrt[i][j][k].tab[c][1][0][0],2)/pow(fNDIS_evt_MC[0][1][i][j][k],3))*pow(fZrange[k],2))  : 0) << " " <<
            (fAcceptance[i][j][k].tab[c][1][0][1] ? Double_t((fGnrt[i][j][k].tab[c][1][0][1]/pow(fNDIS_evt_MC[0][1][i][j][k],2)-pow(fGnrt[i][j][k].tab[c][1][0][1],2)/pow(fNDIS_evt_MC[0][1][i][j][k],3))*pow(fZrange[k],2))  : 0) << " " <<
            (fAcceptance[i][j][k].tab[c][1][0][2] ? Double_t((fGnrt[i][j][k].tab[c][1][0][2]/pow(fNDIS_evt_MC[0][1][i][j][k],2)-pow(fGnrt[i][j][k].tab[c][1][0][2],2)/pow(fNDIS_evt_MC[0][1][i][j][k],3))*pow(fZrange[k],2))  : 0) << " " <<
            (fAcceptance[i][j][k].tab[c][1][0][3] ? Double_t((fGnrt[i][j][k].tab[c][1][0][3]/pow(fNDIS_evt_MC[0][1][i][j][k],2)-pow(fGnrt[i][j][k].tab[c][1][0][3],2)/pow(fNDIS_evt_MC[0][1][i][j][k],3))*pow(fZrange[k],2))  : 0) << " " <<
            (fAcceptance[i][j][k].tab[c][0][0][0] ? Double_t((fGnrt[i][j][k].tab[c][0][0][0]/pow(fNDIS_evt_MC[0][0][i][j][k],2)-pow(fGnrt[i][j][k].tab[c][0][0][0],2)/pow(fNDIS_evt_MC[0][0][i][j][k],3))*pow(fZrange[k],2))  : 0) << " " <<
            (fAcceptance[i][j][k].tab[c][0][0][1] ? Double_t((fGnrt[i][j][k].tab[c][0][0][1]/pow(fNDIS_evt_MC[0][0][i][j][k],2)-pow(fGnrt[i][j][k].tab[c][0][0][1],2)/pow(fNDIS_evt_MC[0][0][i][j][k],3))*pow(fZrange[k],2))  : 0) << " " <<
            (fAcceptance[i][j][k].tab[c][0][0][2] ? Double_t((fGnrt[i][j][k].tab[c][0][0][2]/pow(fNDIS_evt_MC[0][0][i][j][k],2)-pow(fGnrt[i][j][k].tab[c][0][0][2],2)/pow(fNDIS_evt_MC[0][0][i][j][k],3))*pow(fZrange[k],2))  : 0) << " " <<
            (fAcceptance[i][j][k].tab[c][0][0][3] ? Double_t((fGnrt[i][j][k].tab[c][0][0][3]/pow(fNDIS_evt_MC[0][0][i][j][k],2)-pow(fGnrt[i][j][k].tab[c][0][0][3],2)/pow(fNDIS_evt_MC[0][0][i][j][k],3))*pow(fZrange[k],2))  : 0) << " " << endl;

            ofs_zvtx << c << " " << fXrange[i] << " " << fYrange[j] << " " << fZrange[k];

            for(int l=0; l<4; l++)
            {
              ofs_zvtx << " " << fAcceptance_zvtx[i][j][k][l].tab[c][1][0][0] << " " <<
              fAcceptance_zvtx[i][j][k][l].tab[c][1][1][0] << " " <<
              fAcceptance_zvtx[i][j][k][l].tab[c][1][0][1] << " " <<
              fAcceptance_zvtx[i][j][k][l].tab[c][1][1][1] << " " <<
              fAcceptance_zvtx[i][j][k][l].tab[c][1][0][2] << " " <<
              fAcceptance_zvtx[i][j][k][l].tab[c][1][1][2] << " " <<
              fAcceptance_zvtx[i][j][k][l].tab[c][1][0][3] << " " <<
              fAcceptance_zvtx[i][j][k][l].tab[c][1][1][3] << " " <<
              fAcceptance_zvtx[i][j][k][l].tab[c][0][0][0] << " " <<
              fAcceptance_zvtx[i][j][k][l].tab[c][0][1][0] << " " <<
              fAcceptance_zvtx[i][j][k][l].tab[c][0][0][1] << " " <<
              fAcceptance_zvtx[i][j][k][l].tab[c][0][1][1] << " " <<
              fAcceptance_zvtx[i][j][k][l].tab[c][0][0][2] << " " <<
              fAcceptance_zvtx[i][j][k][l].tab[c][0][1][2] << " " <<
              fAcceptance_zvtx[i][j][k][l].tab[c][0][0][3] << " " <<
              fAcceptance_zvtx[i][j][k][l].tab[c][0][1][3];
            }
            ofs_zvtx << endl;

            for(int th=0; th<8; th++)
            {
              ofs_theta << c << " " << fXrange[i] << " " << fYrange[j] << " " << fZrange[k] << " " << fThrange[th];

              for(int l=0; l<4; l++)
              {
                ofs_theta << " " << fAcceptance_theta[i][j][k][th].tab[c][1][0][l] << " " <<
                fAcceptance_theta[i][j][k][th].tab[c][1][1][l];
              }
              for(int l=0; l<4; l++)
              {
                ofs_theta << " " << fAcceptance_theta[i][j][k][th].tab[c][0][0][l] << " " <<
                fAcceptance_theta[i][j][k][th].tab[c][0][1][l];
              }
              ofs_theta << endl;
            }

            for(int pt=0; pt<10; pt++)
            {
              ofs_pt << c << " " << fXrange[i] << " " << fYrange[j] << " " << fZrange[k] << " " << fpTrange[pt];

              for(int l=0; l<4; l++)
              {
                ofs_pt << " " << fAcceptance_pt[i][j][k][pt].tab[c][1][0][l] << " " <<
                fAcceptance_pt[i][j][k][pt].tab[c][1][1][l];
              }
              for(int l=0; l<4; l++)
              {
                ofs_pt << " " << fAcceptance_pt[i][j][k][pt].tab[c][0][0][l] << " " <<
                fAcceptance_pt[i][j][k][pt].tab[c][0][1][l];
              }
              ofs_pt << endl;
            }

            ofs_reld << c << " " << fXrange[i] << " " << fYrange[j] << " " << fZrange[k] << " " <<
            RelDiff(c,i,j,k,1,0) << " " << RelDiff(c,i,j,k,1,1) << " " << RelDiff(c,i,j,k,1,2) << " " << RelDiff(c,i,j,k,1,3) << " " <<
            RelDiff(c,i,j,k,0,0) << " " << RelDiff(c,i,j,k,0,1) << " " << RelDiff(c,i,j,k,0,2) << " " << RelDiff(c,i,j,k,0,3) <<   endl;
          }
        }
      }
    }

    for(int c=0; c<2; c++)
    {
      for(int i=0; i<9; i++)
      {
        for(int k=0; k<12; k++)
        {
          ofs_yavg << c << " " << fXrange[i] << " " << fZrange[k] << " " <<
          fAcceptance_yavg[i][k].tab[c][1][0][0] << " " <<
          fAcceptance_yavg[i][k].tab[c][1][1][0] << " " <<
          fAcceptance_yavg[i][k].tab[c][1][0][1] << " " <<
          fAcceptance_yavg[i][k].tab[c][1][1][1] << " " <<
          fAcceptance_yavg[i][k].tab[c][1][0][2] << " " <<
          fAcceptance_yavg[i][k].tab[c][1][1][2] << " " <<
          fAcceptance_yavg[i][k].tab[c][1][0][3] << " " <<
          fAcceptance_yavg[i][k].tab[c][1][1][3] << " " <<
          fAcceptance_yavg[i][k].tab[c][0][0][0] << " " <<
          fAcceptance_yavg[i][k].tab[c][0][1][0] << " " <<
          fAcceptance_yavg[i][k].tab[c][0][0][1] << " " <<
          fAcceptance_yavg[i][k].tab[c][0][1][1] << " " <<
          fAcceptance_yavg[i][k].tab[c][0][0][2] << " " <<
          fAcceptance_yavg[i][k].tab[c][0][1][2] << " " <<
          fAcceptance_yavg[i][k].tab[c][0][0][3] << " " <<
          fAcceptance_yavg[i][k].tab[c][0][1][3] <<  endl;
        }
      }
    }

    for(int c=1; c>=0; c--)
    {
      for(int i=0; i<9; i++)
      {
        for(int j=0; j<5; j++)
        {
          for(int k=0; k<12; k++)
          {
            ofs_test << fXrange[i] << " " << fYrange[j] << " " << fZrange[k] << " "
                     << fAcceptance[i][j][k].tab[c][1][0][3] << " "
                     << fGnrt[i][j][k].tab[c][1][0][3] << " " << fRcstr_c[i][j][k].tab[c][1][0][3] << " " << fRcstr[i][j][k].tab[c][1][0][3]-fRcstr_c[i][j][k].tab[c][1][0][3] << " "
                     << fNDIS_evt_MC[0][1][i][j][k] << " " << fNDIS_evt[0][1][i][j][k] << " "
                     << sqrt(1/fGnrt[i][j][k].tab[c][1][0][3]+1/fRcstr[i][j][k].tab[c][1][0][3]+1/fNDIS_evt_MC[0][1][i][j][k]+1/fNDIS_evt[0][1][i][j][k]) << " "
                     << sqrt(fAcceptance[i][j][k].tab[c][1][1][3]) << endl;
          }
        }
      }
    }

    c5.Print(Form("%s/%d/hadron_acceptance_%s.pdf",dirroot,year,periodName.c_str()));
    c6.Print(Form("%s/%d/pion_acceptance_%s.pdf",dirroot,year,periodName.c_str()));
    c7.Print(Form("%s/%d/kaon_acceptance_%s.pdf",dirroot,year,periodName.c_str()));
    c88.Print(Form("%s/%d/proton_acceptance_%s.pdf",dirroot,year,periodName.c_str()));
    c8[0]->Print(Form("%s/%d/hadron_acceptance_corr_%s.pdf(",dirroot,year,periodName.c_str()),"pdf");
    for(int i=1; i<11; i++)
      c8[i]->Print(Form("%s/%d/hadron_acceptance_corr_%s.pdf",dirroot,year,periodName.c_str()),"pdf");
    c8[11]->Print(Form("%s/%d/hadron_acceptance_corr_%s.pdf)",dirroot,year,periodName.c_str()),"pdf");
    c9.Print(Form("%s/%d/hadron_acceptance_yavg_%s.pdf",dirroot,year,periodName.c_str()));
    c10.Print(Form("%s/%d/pion_acceptance_yavg_%s.pdf",dirroot,year,periodName.c_str()));
    c11.Print(Form("%s/%d/kaon_acceptance_yavg_%s.pdf",dirroot,year,periodName.c_str()));

    ofs.close();
    ofs_zvtx.close();
    ofs_theta.close();
    ofs_pt.close();
    ofs_yavg.close();
    ofs_reld.close();
    ofs_test.close();

    for(int i=0; i<12; i++)
    {
      delete c8[i];
    }
  }

  return 0;
}

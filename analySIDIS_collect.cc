#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TMatrixD.h>
#include <TMatrixTUtils.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TLine.h>
#include <TGaxis.h>

#include "analySIDIS_collect.h"

#define data_path "/sps/compass/npierre/Multiplicities"
#define proton_sirc "data/proton_semi_inclusive_RC.txt"
#define DVM_2006 "data/DVM_2006.dat"

// Flags
#define Y2006 0
#define Y2012 0
#define Y2016 1
#define DVMC 1
#define SIRC 1
#define NO_ACC 0
#define YMULT 2 // 1: Mean, 2: Weighted Mean, 3: Integration (1 y-bin)
#define STAGGERED 0

using namespace std;

void fetch_acceptance(string pname, int np)
{
  ifstream acc_file(pname);
  double dummy;

  for(int c=0; c<2; c++)
  {
    for(int i=0; i<9; i++)
    {
      for(int j=0; j<6; j++)
      {
        for(int k=0; k<12; k++)
        {
          for(int l=0; l<4; l++)
          {
            acc_file >> dummy;
          }
          acc_file >> fAcceptance[np][i][j][k].tab[c][0][0];
          acc_file >> fAcceptance[np][i][j][k].tab[c][1][0];
          acc_file >> fAcceptance[np][i][j][k].tab[c][0][1];
          acc_file >> fAcceptance[np][i][j][k].tab[c][1][1];
          acc_file >> fAcceptance[np][i][j][k].tab[c][0][2];
          acc_file >> fAcceptance[np][i][j][k].tab[c][1][2];
          acc_file >> fAcceptance[np][i][j][k].tab[c][0][3];
          acc_file >> fAcceptance[np][i][j][k].tab[c][1][3];
#ifdef DEBUG
          cout << fAcceptance[np][i][j][k].tab[c][0][0] << " " <<
          fAcceptance[np][i][j][k].tab[c][1][0] << " " <<
          fAcceptance[np][i][j][k].tab[c][0][1] << " " <<
          fAcceptance[np][i][j][k].tab[c][1][1] << " " <<
          fAcceptance[np][i][j][k].tab[c][0][2] << " " <<
          fAcceptance[np][i][j][k].tab[c][1][2] << " " <<
          fAcceptance[np][i][j][k].tab[c][0][3] << " " <<
          fAcceptance[np][i][j][k].tab[c][1][3] << endl;
#endif
        }
      }
    }
  }

  acc_file.close();

}

void fetch_zvtx_acceptance(string pname, int np)
{
  ifstream acc_file(pname);
  double dummy;

  for(int c=0; c<2; c++)
  {
    for(int i=0; i<9; i++)
    {
      for(int j=0; j<6; j++)
      {
        for(int k=0; k<12; k++)
        {
          for(int l=0; l<4; l++)
          {
            acc_file >> dummy;
          }
          for(int zv=0; zv<4; zv++)
          {
            acc_file >> fAcceptance_zvtx[np][i][j][k][zv].tab[c][0][0];
            acc_file >> fAcceptance_zvtx[np][i][j][k][zv].tab[c][1][0];
            acc_file >> fAcceptance_zvtx[np][i][j][k][zv].tab[c][0][1];
            acc_file >> fAcceptance_zvtx[np][i][j][k][zv].tab[c][1][1];
            acc_file >> fAcceptance_zvtx[np][i][j][k][zv].tab[c][0][2];
            acc_file >> fAcceptance_zvtx[np][i][j][k][zv].tab[c][1][2];
            acc_file >> fAcceptance_zvtx[np][i][j][k][zv].tab[c][0][3];
            acc_file >> fAcceptance_zvtx[np][i][j][k][zv].tab[c][1][3];
#ifdef DEBUG
            cout << fAcceptance_zvtx[np][i][j][k][zv].tab[c][0][0] << " " <<
            fAcceptance_zvtx[np][i][j][k][zv].tab[c][1][0] << " " <<
            fAcceptance_zvtx[np][i][j][k][zv].tab[c][0][1] << " " <<
            fAcceptance_zvtx[np][i][j][k][zv].tab[c][1][1] << " " <<
            fAcceptance_zvtx[np][i][j][k][zv].tab[c][0][2] << " " <<
            fAcceptance_zvtx[np][i][j][k][zv].tab[c][1][2] << " " <<
            fAcceptance_zvtx[np][i][j][k][zv].tab[c][0][3] << " " <<
            fAcceptance_zvtx[np][i][j][k][zv].tab[c][1][3] << endl;
#endif
          }
        }
      }
    }
  }

  acc_file.close();

}

void fetch_yavg_acceptance(string pname, int np)
{
  ifstream acc_file(pname);
  double dummy;

  for(int c=0; c<2; c++)
  {
    for(int i=0; i<9; i++)
    {
      for(int k=0; k<12; k++)
      {
        for(int l=0; l<3; l++)
        {
          acc_file >> dummy;
        }
        acc_file >> fAcceptance_yavg[np][i][k].tab[c][0][0];
        acc_file >> fAcceptance_yavg[np][i][k].tab[c][1][0];
        acc_file >> fAcceptance_yavg[np][i][k].tab[c][0][1];
        acc_file >> fAcceptance_yavg[np][i][k].tab[c][1][1];
        acc_file >> fAcceptance_yavg[np][i][k].tab[c][0][2];
        acc_file >> fAcceptance_yavg[np][i][k].tab[c][1][2];
        acc_file >> fAcceptance_yavg[np][i][k].tab[c][0][3];
        acc_file >> fAcceptance_yavg[np][i][k].tab[c][1][3];
#ifdef DEBUG
        cout << fAcceptance_yavg[np][i][k].tab[c][0][0] << " " <<
        fAcceptance_yavg[np][i][k].tab[c][1][0] << " " <<
        fAcceptance_yavg[np][i][k].tab[c][0][1] << " " <<
        fAcceptance_yavg[np][i][k].tab[c][1][1] << " " <<
        fAcceptance_yavg[np][i][k].tab[c][0][2] << " " <<
        fAcceptance_yavg[np][i][k].tab[c][1][2] << " " <<
        fAcceptance_yavg[np][i][k].tab[c][0][3] << " " <<
        fAcceptance_yavg[np][i][k].tab[c][1][3] << endl;
#endif
      }
    }
  }

  acc_file.close();

}

void dummy_acceptance()
{
  for(int c=0; c<2; c++)
  {
    for(int i=0; i<9; i++)
    {
      for(int j=0; j<6; j++)
      {
        for(int k=0; k<12; k++)
        {
          for(int np=0; np<11; np++)
          {
            fAcceptance[np][i][j][k].tab[c][0][0]=1;
            fAcceptance[np][i][j][k].tab[c][1][0]=0;
            fAcceptance[np][i][j][k].tab[c][0][1]=1;
            fAcceptance[np][i][j][k].tab[c][1][1]=0;
            fAcceptance[np][i][j][k].tab[c][0][2]=1;
            fAcceptance[np][i][j][k].tab[c][1][2]=0;
            fAcceptance[np][i][j][k].tab[c][0][3]=1;
            fAcceptance[np][i][j][k].tab[c][1][3]=0;
            for(int zv=0; zv<4; zv++)
            {
              fAcceptance_zvtx[np][i][j][k][zv].tab[c][0][0]=1;
              fAcceptance_zvtx[np][i][j][k][zv].tab[c][1][0]=0;
              fAcceptance_zvtx[np][i][j][k][zv].tab[c][0][1]=1;
              fAcceptance_zvtx[np][i][j][k][zv].tab[c][1][1]=0;
              fAcceptance_zvtx[np][i][j][k][zv].tab[c][0][2]=1;
              fAcceptance_zvtx[np][i][j][k][zv].tab[c][1][2]=0;
              fAcceptance_zvtx[np][i][j][k][zv].tab[c][0][3]=1;
              fAcceptance_zvtx[np][i][j][k][zv].tab[c][1][3]=0;
            }
          }
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
        for(int np=0; np<11; np++)
        {
          fAcceptance_yavg[np][i][k].tab[c][0][0]=1;
          fAcceptance_yavg[np][i][k].tab[c][1][0]=0;
          fAcceptance_yavg[np][i][k].tab[c][0][1]=1;
          fAcceptance_yavg[np][i][k].tab[c][1][1]=0;
          fAcceptance_yavg[np][i][k].tab[c][0][2]=1;
          fAcceptance_yavg[np][i][k].tab[c][1][2]=0;
          fAcceptance_yavg[np][i][k].tab[c][0][3]=1;
          fAcceptance_yavg[np][i][k].tab[c][1][3]=0;
        }
      }
    }
  }
}

void LoadSemiInclusiveRadiativeCorrection()
{
  string sdum;

  ifstream proton(proton_sirc);

  for(int i=0; i<9; i++)
  {
    for(int j=0; j<6; j++)
    {
      for(int k=0; k<14; k++)
      {
        for(int l=0; l<7; l++)
        {
          proton >> sdum;
#ifdef DEBUG
          cout << sdum << "\t";
#endif
        }
        proton >> fSemiInclusiveRCproton[i][j][k];
#ifdef DEBUG
        cout << fSemiInclusiveRCproton[i][j][k] << "\t";
#endif
        proton >> sdum;
#ifdef DEBUG
        cout << sdum << endl;
#endif
      }
    }
  }
  proton.close();
}

Double_t GetSemiInclusiveRadiativeCorrection(int xb, int yb, int zb)
{
  if(Y2006 || !SIRC)
  {
    return 1;
  }
  else if(Y2012 || Y2016)
  {
    return fSemiInclusiveRCproton[xb][yb][zb];
  }
  else
  {
    cout << "ERROR in GetSemiInclusiveRadiativeCorrection : Year not recognized. No correction applied." << endl;
    return 1;
  }
}

void LoadDiffVectorMesonCorrection()
{
  int x,y,z;
  double dis,had;

  if(DVMC)
  {
    ifstream DVM(DVM_2006);

    while(DVM >> x)
    {
      DVM >> y >> z;
      DVM >> had >> dis;
      fDiffVectorMeson[1][x-1][y-1][z-1][0] = fDiffVectorMeson[1][x-1][y-1][z-1][3] = (DVMC ? had/dis : 1);
      DVM >> had >> dis;
      fDiffVectorMeson[0][x-1][y-1][z-1][0] = fDiffVectorMeson[0][x-1][y-1][z-1][3] = (DVMC ? had/dis : 1);
      DVM >> had >> dis;
      fDiffVectorMeson[1][x-1][y-1][z-1][1] = (DVMC ? had/dis : 1);
      DVM >> had >> dis;
      fDiffVectorMeson[0][x-1][y-1][z-1][1] = (DVMC ? had/dis : 1);
    }
    DVM.close();
  }
  else
  {
    for(int i=0; i<9;i++)
      for(int j=0; j<6;j++)
        for(int k=0; k<12;k++)
          for(int l=0; l<4;l++)
            fDiffVectorMeson[0][i][j][k][l] = fDiffVectorMeson[1][i][j][k][l] = 1;
  }
}

void yavg()
{
  int pMean;

  for(int c=0; c<2; c++)
  {
    for(int x=0; x<9; x++)
    {
      for(int z=0; z<12; z++)
      {
        for(int l=0; l<4; l++)
        {
          pMean = 0;
          for(int i=0; i<6; i++)
          {
            fMultiplicities_yavg[x][z].tab[c][0][l]+=fMultiplicities[x][i][z].tab[c][0][l];
            if(fMultiplicities[x][i][z].tab[c][0][l]) pMean++;
            fMultiplicities_yavg[x][z].tab[c][1][l]+=fMultiplicities[x][i][z].tab[c][1][l];
            fMultiplicities_yavg[x][z].tab[c][2][l]+=fMultiplicities[x][i][z].tab[c][2][l];
          }
          if(pMean)
          {
            fMultiplicities_yavg[x][z].tab[c][0][l]/=pMean;
            fMultiplicities_yavg[x][z].tab[c][1][l]/=pow(pMean,2);
            fMultiplicities_yavg[x][z].tab[c][2][l]/=pow(pMean,2);
          }
        }
      }
    }
  }
}

void yweightedavg()
{
  for(int c=0; c<2; c++)
  {
    for(int x=0; x<9; x++)
    {
      for(int z=0; z<12; z++)
      {
        for(int l=0; l<4; l++)
        {
          for(int i=0; i<6; i++)
          {
            if(fMultiplicities[x][i][z].tab[c][0][l])
            {
              fMultiplicities_yavg[x][z].tab[c][0][l]+=fMultiplicities[x][i][z].tab[c][0][l]/fMultiplicities[x][i][z].tab[c][1][l];
              fMultiplicities_yavg[x][z].tab[c][1][l]+=1/fMultiplicities[x][i][z].tab[c][1][l];
              fMultiplicities_yavg[x][z].tab[c][2][l]+=1/fMultiplicities[x][i][z].tab[c][2][l];
            }
            for(int ll=0; ll<4; ll++)
            {
              if(fMultiplicities_zvtx[x][i][z][ll].tab[c][0][l])
              {
                fMultiplicities_zvtx_yavg[x][z][ll].tab[c][0][l]+=fMultiplicities_zvtx[x][i][z][ll].tab[c][0][l]/fMultiplicities_zvtx[x][i][z][ll].tab[c][1][l];
                fMultiplicities_zvtx_yavg[x][z][ll].tab[c][1][l]+=1/fMultiplicities_zvtx[x][i][z][ll].tab[c][1][l];
                fMultiplicities_zvtx_yavg[x][z][ll].tab[c][2][l]+=1/fMultiplicities_zvtx[x][i][z][ll].tab[c][2][l];
              }
            }
          }
          if(fMultiplicities_yavg[x][z].tab[c][0][l])
          {
            fMultiplicities_yavg[x][z].tab[c][1][l]=1/fMultiplicities_yavg[x][z].tab[c][1][l];
            fMultiplicities_yavg[x][z].tab[c][2][l]=1/fMultiplicities_yavg[x][z].tab[c][2][l];
            fMultiplicities_yavg[x][z].tab[c][0][l]*=fMultiplicities_yavg[x][z].tab[c][1][l];
          }
          for(int ll=0; ll<4; ll++)
          {
            if(fMultiplicities_zvtx_yavg[x][z][ll].tab[c][0][l])
            {
              fMultiplicities_zvtx_yavg[x][z][ll].tab[c][1][l]=1/fMultiplicities_zvtx_yavg[x][z][ll].tab[c][1][l];
              fMultiplicities_zvtx_yavg[x][z][ll].tab[c][2][l]=1/fMultiplicities_zvtx_yavg[x][z][ll].tab[c][2][l];
              fMultiplicities_zvtx_yavg[x][z][ll].tab[c][0][l]*=fMultiplicities_zvtx_yavg[x][z][ll].tab[c][1][l];
            }
          }
        }
      }
    }
  }
}

void yintegrated()
{
  for(int c=0; c<2; c++)
  {
    for(int x=0; x<9; x++)
    {
      for(int z=0; z<12; z++)
      {
        for(int l=0; l<4; l++)
        {
          for(int i=0; i<6; i++)
          {
            for(auto period : fPeriods)
            {
              fBinning_yavg_period[period][x][z].tab[c][0][l]+=fBinning_period[period][x][i][z].tab[c][0][l];
              fBinning_yavg_period[period][x][z].tab[c][1][l]+=fBinning_period[period][x][i][z].tab[c][1][l];
            }
            fNDIS_evt_yavg[0][x][z]+=fNDIS_evt[0][x][i][z];
            fNDIS_evt_err_yavg[0][x][z]+=fNDIS_evt_err[0][x][i][z];
            fBinning_yavg[x][z].tab[c][0][l]+=fBinning[x][i][z].tab[c][0][l];
            fBinning_yavg[x][z].tab[c][1][l]+=fBinning[x][i][z].tab[c][1][l];
            fRich_sys_err_yavg[x][z].tab[c][1][l]+=pow(fRich_sys_err[x][i][z].tab[c][1][l],2);
          }
          fRich_sys_err_yavg[x][z].tab[c][1][l]=sqrt(fRich_sys_err_yavg[x][z].tab[c][1][l]);
        }
      }
    }
  }
}

void compMultiplicitiesIntegratedY()
{
  for(int c=0; c<2; c++)
  {
    for(int i=0; i<9; i++)
    {
      for(int k=0; k<12; k++)
      {
        for(int l=0; l<4; l++)
        {
          fMultiplicities_yavg[i][k].tab[c][0][l] = (fBinning_yavg[i][k].tab[c][0][l] && fNDIS_evt_yavg[0][i][k] && fAcceptance_yavg_weighted[i][k].tab[c][0][l] ?
                                                  fBinning_yavg[i][k].tab[c][0][l]/(fNDIS_evt_yavg[0][i][k]*fZ_bin_width[k]*fAcceptance_yavg_weighted[i][k].tab[c][0][l])
                                                  : 0);
          fMultiplicities_yavg[i][k].tab[c][1][l] = (fNDIS_evt_yavg[0][i][k] && fAcceptance_yavg_weighted[i][k].tab[c][0][l] ?
                                                  Double_t(((fBinning_yavg[i][k].tab[c][1][l]/pow(fNDIS_evt_yavg[0][i][k],2)-pow(fBinning_yavg[i][k].tab[c][0][l],2)*
                                                  fNDIS_evt_err_yavg[0][i][k]/pow(fNDIS_evt_yavg[0][i][k],4))/(pow(fZ_bin_width[k]*fAcceptance_yavg_weighted[i][k].tab[c][0][l],2)))
                                                  + fAcceptance_yavg_weighted[i][k].tab[c][1][l]*pow(fBinning_yavg[i][k].tab[c][0][l]/(fNDIS_evt_yavg[0][i][k]*fZ_bin_width[k]*pow(fAcceptance_yavg_weighted[i][k].tab[c][0][l],2)),2))
                                                  : 0);
          fMultiplicities_yavg[i][k].tab[c][2][l] = (fNDIS_evt_yavg[0][i][k] ?
                                                  Double_t(sqrt(pow(fRich_sys_err_yavg[i][k].tab[c][1][l],2)/pow(fNDIS_evt_yavg[0][i][k]*fZ_bin_width[k]*fAcceptance_yavg_weighted[i][k].tab[c][0][l],2)+
                                                  pow(0.05*sqrt(fAcceptance_yavg_weighted[i][k].tab[c][1][l])*fBinning_yavg[i][k].tab[c][0][l]/(fNDIS_evt_yavg[0][i][k]*fZ_bin_width[k]
                                                  *pow(fAcceptance_yavg_weighted[i][k].tab[c][0][l],2)),2)))
                                                  : 0);

          if(fMultiplicities_yavg[i][k].tab[c][0][l]<0 || fMultiplicities_yavg[i][k].tab[c][0][l]*0.9<fMultiplicities_yavg[i][k].tab[c][1][l])
          {
            fMultiplicities_yavg[i][k].tab[c][0][l] = 0 ;
            fMultiplicities_yavg[i][k].tab[c][1][l] = 0 ;
            fMultiplicities_yavg[i][k].tab[c][2][l] = 0 ;
          }
        }
      }
    }
  }
}

void weight_acceptance()
{
  for(int i=0; i<9; i++)
  {
    for(int k=0; k<12; k++)
    {
      for(int c=1; c>=0; c--)
      {
        for(int l=0; l<4; l++)
        {
          for(auto period : fPeriods)
          {
              for(int j=0; j<6; j++)
              {
                fAcceptance_weighted[i][j][k].tab[c][0][l] += fBinning[i][j][k].tab[c][0][l] ? fBinning_period[period][i][j][k].tab[c][0][l]*fAcceptance[period][i][j][k].tab[c][0][l]/fBinning[i][j][k].tab[c][0][l] : 0;
                fAcceptance_weighted[i][j][k].tab[c][1][l] += fBinning[i][j][k].tab[c][0][l] ? pow(fBinning_period[period][i][j][k].tab[c][0][l],2)*fAcceptance[period][i][j][k].tab[c][1][l]/pow(fBinning[i][j][k].tab[c][0][l],2) : 0;
                for(int zv=0; zv<4; zv++)
                {
                  fAcceptance_weighted_zvtx[i][j][k][zv].tab[c][0][l] += fBinning_zvtx[i][j][k][zv].tab[c][0][l] ? fBinning_period_zvtx[period][i][j][k][zv].tab[c][0][l]*fAcceptance_zvtx[period][i][j][k][zv].tab[c][0][l]/fBinning_zvtx[i][j][k][zv].tab[c][0][l] : 0;
                  // cout << fBinning_zvtx[i][j][k][zv].tab[c][0][l] << " " << fAcceptance_zvtx[period][i][j][k][zv].tab[c][0][l] << " " << fBinning_zvtx[i][j][k][zv].tab[c][0][l] << endl;
                  fAcceptance_weighted_zvtx[i][j][k][zv].tab[c][1][l] += fBinning_zvtx[i][j][k][zv].tab[c][0][l] ? pow(fBinning_period_zvtx[period][i][j][k][zv].tab[c][0][l],2)*fAcceptance_zvtx[period][i][j][k][zv].tab[c][1][l]/pow(fBinning_zvtx[i][j][k][zv].tab[c][0][l],2) : 0;
                }
              }
              fAcceptance_yavg_weighted[i][k].tab[c][0][l] += fBinning_yavg_period[period][i][k].tab[c][0][l]*fAcceptance_yavg[period][i][k].tab[c][0][l]/fBinning_yavg[i][k].tab[c][0][l];
              fAcceptance_yavg_weighted[i][k].tab[c][1][l] += fBinning_yavg_period[period][i][k].tab[c][0][l]*fAcceptance_yavg[period][i][k].tab[c][1][l]/fBinning_yavg[i][k].tab[c][0][l];
          }
        }
      }
    }
  }
}

void weight_meanvalues()
{
  for(int i=0; i<9; i++)
  {
    for(int j=0; j<6; j++)
    {
      for(int k=0; k<12; k++)
      {
        for(int c=1; c>=0; c--)
        {
          for(int ll=0; ll<4; ll++)
          {
            for(auto period : fPeriods)
            {
              fMeanvalues_size[i][j][k].tab[c][ll][0] += fMeanvalues_size_periods[period][i][j][k].tab[c][ll][0];
            }
            for(auto period : fPeriods)
            {
              fMeanvalues_data[i][j][k].tab[c][ll][0] += fMeanvalues_size_periods[period][i][j][k].tab[c][ll][0]*fMeanvalues_data_periods[fNumberPeriod-1][i][j][k].tab[c][ll][0]/fMeanvalues_size[i][j][k].tab[c][ll][0];
              fMeanvalues_data[i][j][k].tab[c][ll][1] += fMeanvalues_size_periods[period][i][j][k].tab[c][ll][0]*fMeanvalues_data_periods[fNumberPeriod-1][i][j][k].tab[c][ll][1]/fMeanvalues_size[i][j][k].tab[c][ll][0];
              fMeanvalues_data[i][j][k].tab[c][ll][2] += fMeanvalues_size_periods[period][i][j][k].tab[c][ll][0]*fMeanvalues_data_periods[fNumberPeriod-1][i][j][k].tab[c][ll][2]/fMeanvalues_size[i][j][k].tab[c][ll][0];
              fMeanvalues_data[i][j][k].tab[c][ll][3] += fMeanvalues_size_periods[period][i][j][k].tab[c][ll][0]*fMeanvalues_data_periods[fNumberPeriod-1][i][j][k].tab[c][ll][3]/fMeanvalues_size[i][j][k].tab[c][ll][0];
            }
          }
        }
      }
    }
  }
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
          for(int ll=0; ll<4; ll++)
          {
            fMeanvalues_temp[i][j][k].tab[c][ll][0] = 0;
            fMeanvalues_temp[i][j][k].tab[c][ll][1] = 0;
            fMeanvalues_temp[i][j][k].tab[c][ll][2] = 0;
            fMeanvalues_temp[i][j][k].tab[c][ll][3] = 0;
            fMeanvalues_size[i][j][k].tab[c][ll][0] = 0;
          }
        }
      }
    }
  }
}

Double_t RelDiff(int c, int x, int y, int z, int had)
{
  Double_t ups=fMultiplicities_zvtx[x][y][z][0].tab[c][0][had];
  Double_t dns=fMultiplicities_zvtx[x][y][z][3].tab[c][0][had];

  return (ups ? Double_t((dns-ups)/ups) : 0);
}

Double_t RelDiff_yavg(int c, int x,int z, int had)
{
  Double_t ups=fMultiplicities_zvtx_yavg[x][z][0].tab[c][0][had];
  Double_t dns=fMultiplicities_zvtx_yavg[x][z][3].tab[c][0][had];

  return (ups ? Double_t((dns-ups)/ups) : 0);
}

Double_t RelDiff_Err(int c, int x, int y, int z, int had)
{
  Double_t ups=fMultiplicities_zvtx[x][y][z][0].tab[c][1][had];
  Double_t dns=fMultiplicities_zvtx[x][y][z][3].tab[c][1][had];

  return (ups ? Double_t((dns-ups)/ups) : 0);
}

Double_t RelDiff_Err_yavg(int c, int x, int z, int had)
{
  Double_t ups=fMultiplicities_zvtx_yavg[x][z][0].tab[c][1][had];
  Double_t dns=fMultiplicities_zvtx_yavg[x][z][3].tab[c][1][had];

  return (ups ? Double_t((dns-ups)/ups) : 0);
}

int main(int argc, char **argv)
{

  if(argc < 2)
  {
    cout << "ERROR : Not enough arguments." << endl;
    cout << "Asked : at least 1 *** Received : " << argc-1 << endl;
    cout << "./analySIDIS_collect periodFile" << endl;

    return 1;
  }

  //q_bin x_bin y_bin z_bin acc_pi acc_error_pi acc_k acc_error_k acc_p acc_error_p acc_h acc_error_h

  int year=0;
  fNumberPeriod=0;

  if(Y2006) year=2006;
  else if(Y2012) year=2012;
  else if(Y2016) year=2016;

  double dummy;

  double rich_sys_err;

  for(int c=0; c<2; c++)
  {
    for(xbin=0; xbin<9; xbin++)
    {
      for(ybin=0; ybin<6; ybin++)
      {
        for(zbin=0; zbin<12; zbin++)
        {
          for(int ll=0; ll<4; ll++)
          {
            rich_sys_err = max(abs(sqrt(fBinning[xbin][ybin][zbin].tab[c][1][ll])),
                           max(abs(fBinning_loose[xbin][ybin][zbin].tab[c][0][ll]-fBinning[xbin][ybin][zbin].tab[c][0][ll]),
                               abs(fBinning_severe[xbin][ybin][zbin].tab[c][0][ll]-fBinning[xbin][ybin][zbin].tab[c][0][ll])));

            fRich_sys_err[xbin][ybin][zbin].tab[c][1][ll] = rich_sys_err;
          }
        }
      }
    }
  }

  LoadSemiInclusiveRadiativeCorrection();
  LoadDiffVectorMesonCorrection();

  ifstream periods(argv[1]);
  string filelist, periodName;
  int periodBit;
  while(periods >> periodName)
  {
    fNumberPeriod++;
    periods >> periodBit;

    if(!periodBit) continue;

    cout << periodName << " ";

    if(!NO_ACC)
    {
      fetch_acceptance(Form("acceptance/%d/acceptance_%s.txt",year,periodName.c_str()),fNumberPeriod-1);
      fetch_yavg_acceptance(Form("acceptance/%d/acceptance_yavg_%s.txt",year,periodName.c_str()),fNumberPeriod-1);
      fetch_zvtx_acceptance(Form("acceptance/%d/acceptance_vtx_%s.txt",year,periodName.c_str()),fNumberPeriod-1);
    }
    else
    {
      dummy_acceptance();
    }

    ifstream dis_file(Form("rawmult/%d/DIS_%s.txt",year,periodName.c_str()));
    ifstream had_file(Form("rawmult/%d/hadron_%s.txt",year,periodName.c_str()));
    ifstream dis_zvtx_file(Form("rawmult/%d/DIS_zvtx_%s.txt",year,periodName.c_str()));
    ifstream had_zvtx_file(Form("rawmult/%d/hadron_zvtx_%s.txt",year,periodName.c_str()));

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
              dis_file >> dummy;
              fNDIS_evt[0][i][j][k] += dummy;
              dis_file >> dummy;
              fNDIS_evt_err[0][i][j][k] += dummy;
              for(int zv=0; zv<4; zv++)
              {
                dis_zvtx_file >> dummy;
                fNDIS_evt_zvtx[0][i][j][k][zv] += dummy;
                dis_zvtx_file >> dummy;
                fNDIS_evt_err_zvtx[0][i][j][k][zv] += dummy;
              }
            }

            for(int ll=0; ll<4; ll++)
            {
              dis_file >> dummy;
              fMeanvalues_temp[i][j][k].tab[c][ll][0] += dummy;
              dis_file >> dummy;
              fMeanvalues_temp[i][j][k].tab[c][ll][1] += dummy;
              dis_file >> dummy;
              fMeanvalues_temp[i][j][k].tab[c][ll][2] += dummy;
              dis_file >> dummy;
              fMeanvalues_temp[i][j][k].tab[c][ll][3] += dummy;
              dis_file >> dummy;
              fMeanvalues_size[i][j][k].tab[c][ll][0] += dummy;
            }

            had_file >> dummy;
            fBinning_period[fNumberPeriod-1][i][j][k].tab[c][0][0] += dummy;
            had_file >> dummy;
            fBinning[i][j][k].tab[c][1][0] += dummy;
            had_file >> dummy;
            fBinning_loose[i][j][k].tab[c][0][0] += dummy;
            had_file >> dummy;
            fBinning_severe[i][j][k].tab[c][0][0] += dummy;
            had_file >> dummy;
            fBinning_period[fNumberPeriod-1][i][j][k].tab[c][0][1] += dummy;
            had_file >> dummy;
            fBinning[i][j][k].tab[c][1][1] += dummy;
            had_file >> dummy;
            fBinning_loose[i][j][k].tab[c][0][1] += dummy;
            had_file >> dummy;
            fBinning_severe[i][j][k].tab[c][0][1] += dummy;
            had_file >> dummy;
            fBinning_period[fNumberPeriod-1][i][j][k].tab[c][0][2] += dummy;
            had_file >> dummy;
            fBinning[i][j][k].tab[c][1][2] += dummy;
            had_file >> dummy;
            fBinning_loose[i][j][k].tab[c][0][2] += dummy;
            had_file >> dummy;
            fBinning_severe[i][j][k].tab[c][0][2] += dummy;
            had_file >> dummy;
            fBinning_period[fNumberPeriod-1][i][j][k].tab[c][0][3] += dummy;
            had_file >> dummy;
            fBinning[i][j][k].tab[c][1][3] += dummy;
            had_file >> dummy;
            fBinning_loose[i][j][k].tab[c][0][3] += dummy;
            had_file >> dummy;
            fBinning_severe[i][j][k].tab[c][0][3] += dummy;
            for(int zv=0; zv<4; zv++)
            {
              had_zvtx_file >> dummy;
              fBinning_period_zvtx[fNumberPeriod-1][i][j][k][zv].tab[c][0][0] += dummy;
              had_zvtx_file >> dummy;
              fBinning_zvtx[i][j][k][zv].tab[c][1][0] += dummy;
              had_zvtx_file >> dummy;
              fBinning_period_zvtx[fNumberPeriod-1][i][j][k][zv].tab[c][0][1] += dummy;
              had_zvtx_file >> dummy;
              fBinning_zvtx[i][j][k][zv].tab[c][1][1] += dummy;
              had_zvtx_file >> dummy;
              fBinning_period_zvtx[fNumberPeriod-1][i][j][k][zv].tab[c][0][2] += dummy;
              had_zvtx_file >> dummy;
              fBinning_zvtx[i][j][k][zv].tab[c][1][2] += dummy;
              had_zvtx_file >> dummy;
              fBinning_period_zvtx[fNumberPeriod-1][i][j][k][zv].tab[c][0][3] += dummy;
              had_zvtx_file >> dummy;
              fBinning_zvtx[i][j][k][zv].tab[c][1][3] += dummy;
            }
          }
        }
      }
    }


    for(int i=0; i<9; i++)
    {
      for(int j=0; j<6; j++)
      {
        for(int k=0; k<12; k++)
        {
          for(int c=0; c<2; c++)
          {
            for(int ll=0; ll<4; ll++)
            {
              fBinning[i][j][k].tab[c][0][ll] += fBinning_period[fNumberPeriod-1][i][j][k].tab[c][0][ll];
              for(int zv=0; zv<4; zv++)
                fBinning_zvtx[i][j][k][zv].tab[c][0][ll] += fBinning_period_zvtx[fNumberPeriod-1][i][j][k][zv].tab[c][0][ll];
              if(int(fMeanvalues_size[i][j][k].tab[c][ll][0]))
              {
                fMeanvalues_data_periods[fNumberPeriod-1][i][j][k].tab[c][ll][0] = fMeanvalues_temp[i][j][k].tab[c][ll][0]/int(fMeanvalues_size[i][j][k].tab[c][ll][0]);
                fMeanvalues_data_periods[fNumberPeriod-1][i][j][k].tab[c][ll][1] = fMeanvalues_temp[i][j][k].tab[c][ll][1]/int(fMeanvalues_size[i][j][k].tab[c][ll][0]);
                fMeanvalues_data_periods[fNumberPeriod-1][i][j][k].tab[c][ll][2] = fMeanvalues_temp[i][j][k].tab[c][ll][2]/int(fMeanvalues_size[i][j][k].tab[c][ll][0]);
                fMeanvalues_data_periods[fNumberPeriod-1][i][j][k].tab[c][ll][3] = fMeanvalues_temp[i][j][k].tab[c][ll][3]/int(fMeanvalues_size[i][j][k].tab[c][ll][0]);
              }
              fMeanvalues_size_periods[fNumberPeriod-1][i][j][k].tab[c][ll][0] = fMeanvalues_size[i][j][k].tab[c][ll][0];
            }
          }
        }
      }
    }

    fPeriods.push_back(fNumberPeriod-1);
    resetValues();
  }

  if(YMULT == 3) yintegrated();
  weight_acceptance();
  weight_meanvalues();
  if(YMULT == 3) compMultiplicitiesIntegratedY();

  TCanvas* c51;
  c51 = new TCanvas("Hadron_Multiplicities_plus","Hadron_Multiplicities_plus",3200,1600);
  TCanvas* c52;
  c52 = new TCanvas("Hadron_Multiplicities_minus","Hadron_Multiplicities_minus",3200,1600);

  TCanvas* c61;
  c61 = new TCanvas("Pion_Multiplicities_plus","Pion_Multiplicities_plus",3200,1600);
  TCanvas* c62;
  c62 = new TCanvas("Pion_Multiplicities_minus","Pion_Multiplicities_minus",3200,1600);

  TCanvas* c71;
  c71 = new TCanvas("Kaon_Multiplicities_plus","Kaon_Multiplicities_plus",3200,1600);
  TCanvas* c72;
  c72 = new TCanvas("Kaon_Multiplicities_minus","Kaon_Multiplicities_minus",3200,1600);

  TCanvas* c5;
  c5 = new TCanvas("Hadron_Multiplicities_zvtx","Hadron_Multiplicities_zvtx",3200,1600);
  TCanvas* c53;
  c53 = new TCanvas("Hadron_Multiplicities_reldiff","Hadron_Multiplicities_reldiff",3200,1600);
  TCanvas* c511;
  c511 = new TCanvas("Hadron_Multiplicities_y_zvtx","Hadron_Multiplicities_y_zvtx",3200,1600);
  TCanvas* c531;
  c531 = new TCanvas("Hadron_Multiplicities_y_reldiff","Hadron_Multiplicities_y_reldiff",3200,1600);

  TCanvas* c6;
  c6 = new TCanvas("Pion_Multiplicities_zvtx","Pion_Multiplicities_zvtx",3200,1600);
  TCanvas* c63;
  c63 = new TCanvas("Pion_Multiplicities_reldiff","Pion_Multiplicities_reldiff",3200,1600);
  TCanvas* c611;
  c611 = new TCanvas("Pion_Multiplicities_y_zvtx","Pion_Multiplicities_y_zvtx",3200,1600);
  TCanvas* c631;
  c631 = new TCanvas("Pion_Multiplicities_y_reldiff","Pion_Multiplicities_y_reldiff",3200,1600);

  TCanvas* c7;
  c7 = new TCanvas("Kaon_Multiplicities_zvtx","Kaon_Multiplicities_zvtx",3200,1600);
  TCanvas* c73;
  c73 = new TCanvas("Kaon_Multiplicities_reldiff","Kaon_Multiplicities_reldiff",3200,1600);
  TCanvas* c711;
  c711 = new TCanvas("Kaon_Multiplicities_y_zvtx","Kaon_Multiplicities_y_zvtx",3200,1600);
  TCanvas* c731;
  c731 = new TCanvas("Kaon_Multiplicities_y_reldiff","Kaon_Multiplicities_y_reldiff",3200,1600);

  TCanvas* c8;
  c8 = new TCanvas("Hadron_Multiplicities_yavg","Hadron_Multiplicities_yavg",3200,1600);

  TCanvas* c9;
  c9 = new TCanvas("Pion_Multiplicities_yavg","Pion_Multiplicities_yavg",3200,1600);

  TCanvas* c10;
  c10 = new TCanvas("Kaon_Multiplicities_yavg","Kaon_Multiplicities_yavg",3200,1600);

  TCanvas* c11;
  c11 = new TCanvas("Hadron_Multiplicities_sum","Hadron_Multiplicities_sum",3200,1600);

  TCanvas* c12;
  c12 = new TCanvas("Pion_Multiplicities_sum","Pion_Multiplicities_sum",3200,1600);

  TCanvas* c13;
  c13 = new TCanvas("Kaon_Multiplicities_sum","Kaon_Multiplicities_sum",3200,1600);

  TCanvas* c14;
  c14 = new TCanvas("Hadron_Multiplicities_ratio","Hadron_Multiplicities_ratio",3200,1600);

  TCanvas* c15;
  c15 = new TCanvas("Pion_Multiplicities_ratio","Pion_Multiplicities_ratio",3200,1600);

  TCanvas* c16;
  c16 = new TCanvas("Kaon_Multiplicities_ratio","Kaon_Multiplicities_ratio",3200,1600);

  c51->SetFillColor(0);
  c52->SetFillColor(0);
  c61->SetFillColor(0);
  c62->SetFillColor(0);
  c71->SetFillColor(0);
  c72->SetFillColor(0);
  c5->SetFillColor(0);
  c53->SetFillColor(0);
  c6->SetFillColor(0);
  c63->SetFillColor(0);
  c7->SetFillColor(0);
  c73->SetFillColor(0);
  c8->SetFillColor(0);
  c9->SetFillColor(0);
  c10->SetFillColor(0);
  c11->SetFillColor(0);
  c12->SetFillColor(0);
  c13->SetFillColor(0);
  c14->SetFillColor(0);
  c15->SetFillColor(0);
  c16->SetFillColor(0);

  c51->Divide(5,2,0,0);
  c52->Divide(5,2,0,0);
  c61->Divide(5,2,0,0);
  c62->Divide(5,2,0,0);
  c71->Divide(5,2,0,0);
  c72->Divide(5,2,0,0);
  c5->Divide(9,5,0,0);
  c53->Divide(9,5,0,0);
  c6->Divide(9,5,0,0);
  c63->Divide(9,5,0,0);
  c7->Divide(9,5,0,0);
  c73->Divide(9,5,0,0);
  c8->Divide(5,2,0,0);
  c9->Divide(5,2,0,0);
  c10->Divide(5,2,0,0);
  c11->Divide(1,1);
  c12->Divide(1,1);
  c13->Divide(1,1);
  c14->Divide(1,1);
  c15->Divide(1,1);
  c16->Divide(1,1);

  TGraphErrors* H_mult[2][9][6];
  TGraphErrors* P_mult[2][9][6];
  TGraphErrors* K_mult[2][9][6];
  TGraphErrors* H_zvtx[2][9][6][4];
  TGraphErrors* P_zvtx[2][9][6][4];
  TGraphErrors* K_zvtx[2][9][6][4];
  TGraphErrors* H_reldiff[2][9][6];
  TGraphErrors* P_reldiff[2][9][6];
  TGraphErrors* K_reldiff[2][9][6];
  TGraphErrors* H_y[2][9];
  TGraphErrors* P_y[2][9];
  TGraphErrors* K_y[2][9];
  TGraphErrors* H_y_zvtx[2][9][4];
  TGraphErrors* P_y_zvtx[2][9][4];
  TGraphErrors* K_y_zvtx[2][9][4];
  TGraphErrors* H_y_reldiff[2][9];
  TGraphErrors* P_y_reldiff[2][9];
  TGraphErrors* K_y_reldiff[2][9];
  TGraphErrors* sH_y;
  TGraphErrors* sP_y;
  TGraphErrors* sK_y;
  TGraphErrors* rH_y;
  TGraphErrors* rP_y;
  TGraphErrors* rK_y;

  TGraphAsymmErrors* H_sys[2][9][5];
  TGraphAsymmErrors* P_sys[2][9][5];
  TGraphAsymmErrors* K_sys[2][9][5];
  TGraphAsymmErrors* H_ysys[2][9];
  TGraphAsymmErrors* P_ysys[2][9];
  TGraphAsymmErrors* K_ysys[2][9];

  Double_t z_range[12] = {.225,.275,.325,.375,.425,.475,.525,.575,.625,.675,.725,.8};
  Double_t x_range[9] = {.008,.015,.025,.035,.05,.08,.12,.16,.29};
  Double_t errorx[12] = {0.05/2,0.05/2,0.05/2,0.05/2,0.05/2,0.05/2,0.05/2,0.05/2,0.05/2,0.05/2,0.05/2,0.1/2};
  Double_t h_yoffset[12] = {-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2};
  Double_t p_yoffset[12] = {-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2};
  Double_t k_yoffset[12] = {-0.037,-0.037,-0.037,-0.037,-0.037,-0.037,-0.037,-0.037,-0.037,-0.037,-0.037,-0.037};
  Double_t h_yoffset2[12] = {-0.3,-0.3,-0.3,-0.3,-0.3,-0.3,-0.3,-0.3,-0.3,-0.3,-0.3,-0.3};
  Double_t p_yoffset2[12] = {-0.3,-0.3,-0.3,-0.3,-0.3,-0.3,-0.3,-0.3,-0.3,-0.3,-0.3,-0.3};
  Double_t k_yoffset2[12] = {-0.065,-0.065,-0.065,-0.065,-0.065,-0.065,-0.065,-0.065,-0.065,-0.065,-0.065,-0.065};

  ofstream ofs_p(Form("%s/multiplicities_pion.txt",data_path), std::ofstream::out | std::ofstream::trunc);
  ofstream ofs_t(Form("%s/multiplicities_raw.txt",data_path), std::ofstream::out | std::ofstream::trunc);
  ofstream ofs_k(Form("%s/multiplicities_kaon.txt",data_path), std::ofstream::out | std::ofstream::trunc);
  ofstream ofs_h(Form("%s/multiplicities_hadron.txt",data_path), std::ofstream::out | std::ofstream::trunc);
  ofstream ofs_m(Form("%s/multiplicities_forMarcin.txt",data_path), std::ofstream::out | std::ofstream::trunc);
  ofstream ofs_mp(Form("%s/multiplicities_h+.txt",data_path), std::ofstream::out | std::ofstream::trunc);
  ofstream ofs_mm(Form("%s/multiplicities_h-.txt",data_path), std::ofstream::out | std::ofstream::trunc);

  std::vector<Double_t> p_m[2][9][6];
  std::vector<Double_t> k_m[2][9][6];
  std::vector<Double_t> h_m[2][9][6];
  std::vector<Double_t> p_err[2][9][6];
  std::vector<Double_t> k_err[2][9][6];
  std::vector<Double_t> h_err[2][9][6];
  std::vector<Double_t> p_sys[2][9][6];
  std::vector<Double_t> k_sys[2][9][6];
  std::vector<Double_t> h_sys[2][9][6];
  std::vector<Double_t> z_range_p[2][9][6];
  std::vector<Double_t> z_range_k[2][9][6];
  std::vector<Double_t> z_range_h[2][9][6];
  std::vector<Double_t> p_z[2][9][6][4];
  std::vector<Double_t> k_z[2][9][6][4];
  std::vector<Double_t> h_z[2][9][6][4];
  std::vector<Double_t> p_z_err[2][9][6][4];
  std::vector<Double_t> k_z_err[2][9][6][4];
  std::vector<Double_t> h_z_err[2][9][6][4];
  std::vector<Double_t> p_z_sys[2][9][6][4];
  std::vector<Double_t> k_z_sys[2][9][6][4];
  std::vector<Double_t> h_z_sys[2][9][6][4];
  std::vector<Double_t> z_range_p_z[2][9][6][4];
  std::vector<Double_t> z_range_k_z[2][9][6][4];
  std::vector<Double_t> z_range_h_z[2][9][6][4];
  std::vector<Double_t> p_reldiff[2][9][6];
  std::vector<Double_t> k_reldiff[2][9][6];
  std::vector<Double_t> h_reldiff[2][9][6];
  std::vector<Double_t> p_reldiff_err[2][9][6];
  std::vector<Double_t> k_reldiff_err[2][9][6];
  std::vector<Double_t> h_reldiff_err[2][9][6];
  std::vector<Double_t> z_range_p_reldiff[2][9][6];
  std::vector<Double_t> z_range_k_reldiff[2][9][6];
  std::vector<Double_t> z_range_h_reldiff[2][9][6];
  std::vector<Double_t> p_y[2][9];
  std::vector<Double_t> k_y[2][9];
  std::vector<Double_t> h_y[2][9];
  std::vector<Double_t> p_y_err[2][9];
  std::vector<Double_t> k_y_err[2][9];
  std::vector<Double_t> h_y_err[2][9];
  std::vector<Double_t> p_y_sys[2][9];
  std::vector<Double_t> k_y_sys[2][9];
  std::vector<Double_t> h_y_sys[2][9];
  std::vector<Double_t> z_range_p_y[2][9];
  std::vector<Double_t> z_range_k_y[2][9];
  std::vector<Double_t> z_range_h_y[2][9];
  std::vector<Double_t> p_y_z[2][9][4];
  std::vector<Double_t> k_y_z[2][9][4];
  std::vector<Double_t> h_y_z[2][9][4];
  std::vector<Double_t> p_y_z_err[2][9][4];
  std::vector<Double_t> k_y_z_err[2][9][4];
  std::vector<Double_t> h_y_z_err[2][9][4];
  std::vector<Double_t> p_y_z_sys[2][9][4];
  std::vector<Double_t> k_y_z_sys[2][9][4];
  std::vector<Double_t> h_y_z_sys[2][9][4];
  std::vector<Double_t> z_range_p_y_z[2][9][4];
  std::vector<Double_t> z_range_k_y_z[2][9][4];
  std::vector<Double_t> z_range_h_y_z[2][9][4];
  std::vector<Double_t> p_y_reldiff[2][9];
  std::vector<Double_t> k_y_reldiff[2][9];
  std::vector<Double_t> h_y_reldiff[2][9];
  std::vector<Double_t> p_y_reldiff_err[2][9];
  std::vector<Double_t> k_y_reldiff_err[2][9];
  std::vector<Double_t> h_y_reldiff_err[2][9];
  std::vector<Double_t> z_range_p_y_reldiff[2][9];
  std::vector<Double_t> z_range_k_y_reldiff[2][9];
  std::vector<Double_t> z_range_h_y_reldiff[2][9];
  std::vector<Double_t> sp_y;
  std::vector<Double_t> sk_y;
  std::vector<Double_t> sh_y;
  std::vector<Double_t> sp_y_err;
  std::vector<Double_t> sk_y_err;
  std::vector<Double_t> sh_y_err;
  std::vector<Double_t> sx_range_p_y;
  std::vector<Double_t> sx_range_k_y;
  std::vector<Double_t> sx_range_h_y;
  std::vector<Double_t> rp_y;
  std::vector<Double_t> rk_y;
  std::vector<Double_t> rh_y;
  std::vector<Double_t> rp_y_err;
  std::vector<Double_t> rk_y_err;
  std::vector<Double_t> rh_y_err;
  std::vector<Double_t> rx_range_p_y;
  std::vector<Double_t> rx_range_k_y;
  std::vector<Double_t> rx_range_h_y;

  TLine l1(0.1,0.9,0.9,0.9);
  TLine l2(0.1,1.1,0.9,1.1);
  TLine l3(0.1,0.95,0.9,0.95);
  TLine l4(0.1,1.05,0.9,1.05);
  l3.SetLineStyle(2); l4.SetLineStyle(2);

  TLine lsys(0.1,0,0.9,0);
  lsys.SetLineStyle(2);

  for(int i=0; i<9; i++)
  {
    int axisflagh1 = 0;
    int axisflagh2 = 0;
    int axisflagp1 = 0;
    int axisflagp2 = 0;
    int axisflagk1 = 0;
    int axisflagk2 = 0;

    for(int j=0; j<6; j++)
    {
      for(int k=0; k<12; k++)
      {
        for(int c=1; c>=0; c--)
        {
          for(int l=0; l<4; l++)
          {
            fMultiplicities[i][j][k].tab[c][0][l] = (fBinning[i][j][k].tab[c][0][l] && fNDIS_evt[0][i][j][k] && fAcceptance_weighted[i][j][k].tab[c][0][l] ?
                                                    Double_t(fDiffVectorMeson[c][i][j][k][l]*GetSemiInclusiveRadiativeCorrection(i,j,k+1)*fBinning[i][j][k].tab[c][0][l]/(fNDIS_evt[0][i][j][k]*fZ_bin_width[k]*fAcceptance_weighted[i][j][k].tab[c][0][l]))
                                                    : 0);
            fMultiplicities[i][j][k].tab[c][1][l] = (fNDIS_evt[0][i][j][k] && fAcceptance_weighted[i][j][k].tab[c][0][l] ?
                                                    Double_t(((fBinning[i][j][k].tab[c][1][l]/pow(fNDIS_evt[0][i][j][k],2)-pow(fBinning[i][j][k].tab[c][0][l],2)*
                                                    fNDIS_evt_err[0][i][j][k]/pow(fNDIS_evt[0][i][j][k],4))*pow(fDiffVectorMeson[c][i][j][k][l]*GetSemiInclusiveRadiativeCorrection(i,j,k+1)/(fZ_bin_width[k]*fAcceptance_weighted[i][j][k].tab[c][0][l]),2))
                                                    + fAcceptance_weighted[i][j][k].tab[c][1][l]*pow(fDiffVectorMeson[c][i][j][k][l]*GetSemiInclusiveRadiativeCorrection(i,j,k+1)*fBinning[i][j][k].tab[c][0][l]/(fNDIS_evt[0][i][j][k]*fZ_bin_width[k]*pow(fAcceptance_weighted[i][j][k].tab[c][0][l],2)),2))
                                                    : 0);
                                                    // if(l==3 && c==1) cout << i << " "  << j << " " << k << " " << ((fBinning[i][j][k].tab[c][1][l]/pow(fNDIS_evt[0][i][j][k],2)-pow(fBinning[i][j][k].tab[c][0][l],2)*
                                                    // fNDIS_evt_err[0][i][j][k]/pow(fNDIS_evt[0][i][j][k],4))*pow(fDiffVectorMeson[c][i][j][k][l]*GetSemiInclusiveRadiativeCorrection(i,j,k+1)/(fZ_bin_width[k]*fAcceptance_weighted[i][j][k].tab[c][0][l]),2)) << " "
                                                    // << fAcceptance_weighted[i][j][k].tab[c][1][l]*pow(fDiffVectorMeson[c][i][j][k][l]*GetSemiInclusiveRadiativeCorrection(i,j,k+1)*fBinning[i][j][k].tab[c][0][l]/(fNDIS_evt[0][i][j][k]*fZ_bin_width[k]*pow(fAcceptance_weighted[i][j][k].tab[c][0][l],2)),2)
                                                    // << " " << pow(fDiffVectorMeson[c][i][j][k][l]*GetSemiInclusiveRadiativeCorrection(i,j,k+1)*fBinning[i][j][k].tab[c][0][l]/(fNDIS_evt[0][i][j][k]*fZ_bin_width[k]*pow(fAcceptance_weighted[i][j][k].tab[c][0][l],2)),2) << " " << fAcceptance_weighted[i][j][k].tab[c][1][l] << endl;
            fMultiplicities[i][j][k].tab[c][2][l] = (fNDIS_evt[0][i][j][k] ?
                                                    Double_t(sqrt(pow(fDiffVectorMeson[c][i][j][k][l]*GetSemiInclusiveRadiativeCorrection(i,j,k+1)*fRich_sys_err[i][j][k].tab[c][1][l],2)/pow(fNDIS_evt[0][i][j][k]*fZ_bin_width[k]*fAcceptance_weighted[i][j][k].tab[c][0][l],2)+
                                                    pow(0.05*sqrt(fAcceptance_weighted[i][j][k].tab[c][1][l])*fDiffVectorMeson[c][i][j][k][l]*GetSemiInclusiveRadiativeCorrection(i,j,k+1)*fBinning[i][j][k].tab[c][0][l]/(fNDIS_evt[0][i][j][k]*fZ_bin_width[k]
                                                    *pow(fAcceptance_weighted[i][j][k].tab[c][0][l],2)),2)))
                                                    : 0);

            if(fMultiplicities[i][j][k].tab[c][0][l]<=0 || fMultiplicities[i][j][k].tab[c][0][l]*0.9<fMultiplicities[i][j][k].tab[c][1][l])
            {
              fMultiplicities[i][j][k].tab[c][0][l] = 0 ;
              fMultiplicities[i][j][k].tab[c][1][l] = 0 ;
              fMultiplicities[i][j][k].tab[c][2][l] = 0 ;
            }

            for(int zv=0; zv<4; zv++)
            {
              fMultiplicities_zvtx[i][j][k][zv].tab[c][0][l] = (fBinning_zvtx[i][j][k][zv].tab[c][0][l] && fNDIS_evt_zvtx[0][i][j][k][zv] && fAcceptance_weighted_zvtx[i][j][k][zv].tab[c][0][l] ?
                                                      Double_t(fDiffVectorMeson[c][i][j][k][l]*GetSemiInclusiveRadiativeCorrection(i,j,k+1)*fBinning_zvtx[i][j][k][zv].tab[c][0][l]/(fNDIS_evt_zvtx[0][i][j][k][zv]*fZ_bin_width[k]*fAcceptance_weighted_zvtx[i][j][k][zv].tab[c][0][l]))
                                                      : 0);
              fMultiplicities_zvtx[i][j][k][zv].tab[c][1][l] = (fNDIS_evt_zvtx[0][i][j][k][zv] && fAcceptance_weighted_zvtx[i][j][k][zv].tab[c][0][l] ?
                                                      Double_t(((fBinning_zvtx[i][j][k][zv].tab[c][1][l]/pow(fNDIS_evt_zvtx[0][i][j][k][zv],2)-pow(fBinning_zvtx[i][j][k][zv].tab[c][0][l],2)*
                                                      fNDIS_evt_err_zvtx[0][i][j][k][zv]/pow(fNDIS_evt_zvtx[0][i][j][k][zv],4))*pow(fDiffVectorMeson[c][i][j][k][l]*GetSemiInclusiveRadiativeCorrection(i,j,k+1)/(fZ_bin_width[k]*fAcceptance_weighted_zvtx[i][j][k][zv].tab[c][0][l]),2))
                                                      + fAcceptance_weighted_zvtx[i][j][k][zv].tab[c][1][l]*pow(fDiffVectorMeson[c][i][j][k][l]*GetSemiInclusiveRadiativeCorrection(i,j,k+1)*fBinning_zvtx[i][j][k][zv].tab[c][0][l]/(fNDIS_evt_zvtx[0][i][j][k][zv]*fZ_bin_width[k]*pow(fAcceptance_weighted_zvtx[i][j][k][zv].tab[c][0][l],2)),2))
                                                      : 0);

              if(fMultiplicities_zvtx[i][j][k][zv].tab[c][0][l]<=0 || fMultiplicities_zvtx[i][j][k][zv].tab[c][0][l]*0.9<fMultiplicities_zvtx[i][j][k][zv].tab[c][1][l])
              {
                fMultiplicities_zvtx[i][j][k][zv].tab[c][0][l] = 0 ;
                fMultiplicities_zvtx[i][j][k][zv].tab[c][1][l] = 0 ;
                fMultiplicities_zvtx[i][j][k][zv].tab[c][2][l] = 0 ;
              }

              // cout << c << " " << i << " " << j << " " << k << " " << l << " " << zv << " " << fNDIS_evt_zvtx[0][i][j][k][zv] << " " << fBinning_zvtx[i][j][k][zv].tab[c][0][l] << " " << fAcceptance_weighted_zvtx[i][j][k][zv].tab[c][0][l] << " " <<  fMultiplicities_zvtx[i][j][k][zv].tab[c][0][l] << endl;
            }
          }

          for(int zv=0; zv<4; zv++)
          {
            p_z[c][i][j][zv].push_back(fMultiplicities_zvtx[i][j][k][zv].tab[c][0][0]>0 ? fMultiplicities_zvtx[i][j][k][zv].tab[c][0][0] : 0);
            k_z[c][i][j][zv].push_back(fMultiplicities_zvtx[i][j][k][zv].tab[c][0][1]>0 ? fMultiplicities_zvtx[i][j][k][zv].tab[c][0][1] : 0);
            h_z[c][i][j][zv].push_back(fMultiplicities_zvtx[i][j][k][zv].tab[c][0][3]>0 ? fMultiplicities_zvtx[i][j][k][zv].tab[c][0][3] : 0);
            p_z_err[c][i][j][zv].push_back(fMultiplicities_zvtx[i][j][k][zv].tab[c][1][0] ? sqrt(fMultiplicities_zvtx[i][j][k][zv].tab[c][1][0]) : 0);
            k_z_err[c][i][j][zv].push_back(fMultiplicities_zvtx[i][j][k][zv].tab[c][1][1] ? sqrt(fMultiplicities_zvtx[i][j][k][zv].tab[c][1][1]) : 0);
            h_z_err[c][i][j][zv].push_back(fMultiplicities_zvtx[i][j][k][zv].tab[c][1][3] ? sqrt(fMultiplicities_zvtx[i][j][k][zv].tab[c][1][3]) : 0);
          }

          p_reldiff[c][i][j].push_back(RelDiff(c,i,j,k,0));
          k_reldiff[c][i][j].push_back(RelDiff(c,i,j,k,1));
          h_reldiff[c][i][j].push_back(RelDiff(c,i,j,k,3));
          p_reldiff_err[c][i][j].push_back(RelDiff_Err(c,i,j,k,0));
          k_reldiff_err[c][i][j].push_back(RelDiff_Err(c,i,j,k,1));
          h_reldiff_err[c][i][j].push_back(RelDiff_Err(c,i,j,k,3));

          if(c) ofs_p << fXrange[i] << " " << fYrange[j] << " " << fZrange[k] << " ";

          ofs_p <<
          fMeanvalues_data[i][j][k].tab[0][0][0] << " " << fMeanvalues_data[i][j][k].tab[0][0][1] << " " <<
          fMeanvalues_data[i][j][k].tab[0][0][2] << " " << fMeanvalues_data[i][j][k].tab[0][0][3] << " " <<
          fMultiplicities[i][j][k].tab[c][0][0] << " " <<
          fMultiplicities[i][j][k].tab[c][1][0] << " " <<
          fMultiplicities[i][j][k].tab[c][2][0] << " " <<
          (fMultiplicities[i][j][k].tab[c][0][0] ? 1 : 0) << " ";

          if(!c) ofs_p << endl;

          if(c) ofs_t << fXrange[i] << " " << fYrange[j] << " " << fZrange[k] << " ";

          ofs_t <<
          fMeanvalues_data[i][j][k].tab[0][0][0] << " " << fMeanvalues_data[i][j][k].tab[0][0][1] << " " <<
          fMeanvalues_data[i][j][k].tab[0][0][2] << " " << fMeanvalues_data[i][j][k].tab[0][0][3] << " " <<
          fMultiplicities[i][j][k].tab[c][0][0] << " ";

          if(!c) ofs_t << endl;

          if (c) ofs_k << fXrange[i] << " " << fYrange[j] << " " << fZrange[k] << " ";

          ofs_k <<
          fMeanvalues_data[i][j][k].tab[0][1][0] << " " << fMeanvalues_data[i][j][k].tab[0][1][1] << " " <<
          fMeanvalues_data[i][j][k].tab[0][1][2] << " " << fMeanvalues_data[i][j][k].tab[0][1][3] << " " <<
          fMultiplicities[i][j][k].tab[c][0][1] << " " <<
          fMultiplicities[i][j][k].tab[c][1][1] << " " <<
          fMultiplicities[i][j][k].tab[c][2][1] << " " <<
          (fMultiplicities[i][j][k].tab[c][0][1] ? 1 : 0) << " ";

          if(!c) ofs_k << endl;

          if (c) ofs_h << fXrange[i] << " " << fYrange[j] << " " << fZrange[k] << " ";

          ofs_h <<
          fMeanvalues_data[i][j][k].tab[0][3][0] << " " << fMeanvalues_data[i][j][k].tab[0][3][1] << " " <<
          fMeanvalues_data[i][j][k].tab[0][3][2] << " " << fMeanvalues_data[i][j][k].tab[0][3][3] << " " <<
          fMultiplicities[i][j][k].tab[c][0][3] << " " <<
          fMultiplicities[i][j][k].tab[c][1][3] << " " <<
          fMultiplicities[i][j][k].tab[c][2][3] << " " <<
          (fMultiplicities[i][j][k].tab[c][0][3] ? 1 : 0) << " ";

          if(!c) ofs_h << endl;

          if(c) ofs_m << fXrange[i] << " " << fYrange[j] << " " << fZrange[k] << " ";

          ofs_m <<
          fBinning[i][j][k].tab[c][0][3] << " " << fNDIS_evt[0][i][j][k] << " " <<
          fAcceptance_weighted[i][j][k].tab[c][0][3] << " " << GetSemiInclusiveRadiativeCorrection(i,j,k+1) << " " <<
          fDiffVectorMeson[c][i][j][k][3] << " ";

          if(!c) ofs_m << endl;

          if(c)
          {
            ofs_mp << fXrange[i] << " 0 0 " << fYrange[j] << " 0 0 0 " << fZrange[k] << " 0 0 " <<
            fMultiplicities[i][j][k].tab[c][0][3] << " " <<
            fMultiplicities[i][j][k].tab[c][1][3] << " " << fMultiplicities[i][j][k].tab[c][1][3] << " " <<
            fMultiplicities[i][j][k].tab[c][2][3] << " " << fMultiplicities[i][j][k].tab[c][2][3] << endl;
          }

          if(!c)
          {
            ofs_mm << fXrange[i] << " 0 0 " << fYrange[j] << " 0 0 0 " << fZrange[k] << " 0 0 " <<
            fMultiplicities[i][j][k].tab[c][0][3] << " " <<
            fMultiplicities[i][j][k].tab[c][1][3] << " " << fMultiplicities[i][j][k].tab[c][1][3] << " " <<
            fMultiplicities[i][j][k].tab[c][2][3] << " " << fMultiplicities[i][j][k].tab[c][2][3] << endl;
          }

          // if(c==1) cout << i << " " << j << " " << k << " " << fMultiplicities[i][j][k].tab[c][0][3] << " " << fMultiplicities[i][j][k].tab[c][1][3] << " " << sqrt(fMultiplicities[i][j][k].tab[c][1][3]) << endl;

          if(STAGGERED)
          {
            p_m[c][i][j].push_back(fMultiplicities[i][j][k].tab[c][0][0]>0 ? fMultiplicities[i][j][k].tab[c][0][0]+j*0.25 : 0);
            k_m[c][i][j].push_back(fMultiplicities[i][j][k].tab[c][0][1]>0 ? fMultiplicities[i][j][k].tab[c][0][1]+j*0.05 : 0);
            h_m[c][i][j].push_back(fMultiplicities[i][j][k].tab[c][0][3]>0 ? fMultiplicities[i][j][k].tab[c][0][3]+j*0.25 : 0);
          }
          else
          {
            p_m[c][i][j].push_back(fMultiplicities[i][j][k].tab[c][0][0]>0 ? fMultiplicities[i][j][k].tab[c][0][0] : 0);
            k_m[c][i][j].push_back(fMultiplicities[i][j][k].tab[c][0][1]>0 ? fMultiplicities[i][j][k].tab[c][0][1] : 0);
            h_m[c][i][j].push_back(fMultiplicities[i][j][k].tab[c][0][3]>0 ? fMultiplicities[i][j][k].tab[c][0][3] : 0);
          }
          p_err[c][i][j].push_back(fMultiplicities[i][j][k].tab[c][1][0] ? sqrt(fMultiplicities[i][j][k].tab[c][1][0]) : 0);
          k_err[c][i][j].push_back(fMultiplicities[i][j][k].tab[c][1][1] ? sqrt(fMultiplicities[i][j][k].tab[c][1][1]) : 0);
          h_err[c][i][j].push_back(fMultiplicities[i][j][k].tab[c][1][3] ? sqrt(fMultiplicities[i][j][k].tab[c][1][3]) : 0);
          p_sys[c][i][j].push_back(fMultiplicities[i][j][k].tab[c][2][0] ? sqrt(fMultiplicities[i][j][k].tab[c][2][0]) : 0);
          k_sys[c][i][j].push_back(fMultiplicities[i][j][k].tab[c][2][1] ? sqrt(fMultiplicities[i][j][k].tab[c][2][1]) : 0);
          h_sys[c][i][j].push_back(fMultiplicities[i][j][k].tab[c][2][3] ? sqrt(fMultiplicities[i][j][k].tab[c][2][3]) : 0);
        }
      }

      for(int c=0; c<2; c++)
      {
        for(int l=0; l<12; l++)
        {
          z_range_p[c][i][j].push_back(z_range[l]);
          z_range_k[c][i][j].push_back(z_range[l]);
          z_range_h[c][i][j].push_back(z_range[l]);
        }

        // cout << c << " " << i << " " << j << " ";
        //
        // for(int k=0; k<12; k++)
        // {
        //   cout << p_m[c][i][j][k] << " ";
        // }
        //
        // cout << endl;

        for(int k=12; k>0; k--)
        {
          if(!p_m[c][i][j][k-1]) {p_m[c][i][j].erase(p_m[c][i][j].begin()+k-1); p_err[c][i][j].erase(p_err[c][i][j].begin()+k-1); p_sys[c][i][j].erase(p_sys[c][i][j].begin()+k-1); z_range_p[c][i][j].erase(z_range_p[c][i][j].begin()+k-1);}
          if(!k_m[c][i][j][k-1]) {k_m[c][i][j].erase(k_m[c][i][j].begin()+k-1); k_err[c][i][j].erase(k_err[c][i][j].begin()+k-1); k_sys[c][i][j].erase(k_sys[c][i][j].begin()+k-1); z_range_k[c][i][j].erase(z_range_k[c][i][j].begin()+k-1);}
          if(!h_m[c][i][j][k-1]) {h_m[c][i][j].erase(h_m[c][i][j].begin()+k-1); h_err[c][i][j].erase(h_err[c][i][j].begin()+k-1); h_sys[c][i][j].erase(h_sys[c][i][j].begin()+k-1); z_range_h[c][i][j].erase(z_range_h[c][i][j].begin()+k-1);}
        }

        bool p_m_empty = 0;
        bool k_m_empty = 0;
        bool h_m_empty = 0;

        if(!(Int_t(p_m[c][i][j].size()))) p_m_empty = 1;
        if(!(Int_t(k_m[c][i][j].size()))) k_m_empty = 1;
        if(!(Int_t(h_m[c][i][j].size()))) h_m_empty = 1;

        H_mult[c][i][j] = new TGraphErrors(Int_t(h_m[c][i][j].size()),&(z_range_h[c][i][j][0]),&(h_m[c][i][j][0]),0,&(h_err[c][i][j][0]));
        P_mult[c][i][j] = new TGraphErrors(Int_t(p_m[c][i][j].size()),&(z_range_p[c][i][j][0]),&(p_m[c][i][j][0]),0,&(p_err[c][i][j][0]));
        K_mult[c][i][j] = new TGraphErrors(Int_t(k_m[c][i][j].size()),&(z_range_k[c][i][j][0]),&(k_m[c][i][j][0]),0,&(k_err[c][i][j][0]));
        H_sys[c][i][j] = new TGraphAsymmErrors(Int_t(h_m[c][i][j].size()),&(z_range_h[c][i][j][0]), &h_yoffset[0], &errorx[0], &errorx[0], 0, &(h_sys[c][i][j][0]));
        P_sys[c][i][j] = new TGraphAsymmErrors(Int_t(p_m[c][i][j].size()),&(z_range_p[c][i][j][0]), &p_yoffset[0], &errorx[0], &errorx[0], 0, &(p_sys[c][i][j][0]));
        K_sys[c][i][j] = new TGraphAsymmErrors(Int_t(k_m[c][i][j].size()),&(z_range_k[c][i][j][0]), &k_yoffset[0], &errorx[0], &errorx[0], 0, &(k_sys[c][i][j][0]));

        H_mult[c][i][j]->SetMarkerColor(fMarkerColor[j]);
        P_mult[c][i][j]->SetMarkerColor(fMarkerColor[j]);
        K_mult[c][i][j]->SetMarkerColor(fMarkerColor[j]);

        H_mult[c][i][j]->SetMarkerSize(2);
        P_mult[c][i][j]->SetMarkerSize(2);
        K_mult[c][i][j]->SetMarkerSize(2);

        H_mult[c][i][j]->SetMarkerStyle(fMarkerStyle[0][c]);
        P_mult[c][i][j]->SetMarkerStyle(fMarkerStyle[0][c]);
        K_mult[c][i][j]->SetMarkerStyle(fMarkerStyle[0][c]);

        H_mult[c][i][j]->SetTitle("");
        P_mult[c][i][j]->SetTitle("");
        K_mult[c][i][j]->SetTitle("");

        H_mult[c][i][j]->GetYaxis()->SetTitle("");
        P_mult[c][i][j]->GetYaxis()->SetTitle("");
        K_mult[c][i][j]->GetYaxis()->SetTitle("");

        H_sys[c][i][j]->SetFillColor(fMarkerColor[j]);
        P_sys[c][i][j]->SetFillColor(fMarkerColor[j]);
        K_sys[c][i][j]->SetFillColor(fMarkerColor[j]);

        if(!h_m_empty)
        {
          if(c) c51->cd(i+1);
          else c52->cd(i+1);
          if(H_mult[c][i][j])
          {
            if((!c && !axisflagh1) || (c && !axisflagh2))
            {
              H_mult[c][i][j]->Draw("SAMEPA");
              H_mult[c][i][j]->GetXaxis()->SetLimits(0.1,0.9);
              H_mult[c][i][j]->SetMinimum(-0.4);
              H_mult[c][i][j]->SetMaximum(5.);
              H_mult[c][i][j]->GetXaxis()->SetLabelSize(0.06);
              H_mult[c][i][j]->GetYaxis()->SetLabelSize(0.06);
              H_mult[c][i][j]->SetTitle("");
              if(i>4) gPad->SetBottomMargin(.15);
              if(i==0 || i==5) gPad->SetLeftMargin(.22);
              if(i==8)
              {
                H_mult[c][i][j]->GetXaxis()->SetTitle("#font[ 12]{z}");
                H_mult[c][i][j]->GetXaxis()->SetTitleSize(0.08);
                H_mult[c][i][j]->GetXaxis()->SetTitleOffset(.8);
              }
              H_mult[c][i][j]->GetXaxis()->SetNdivisions(304,kTRUE);
              H_mult[c][i][j]->GetYaxis()->SetNdivisions(304,kTRUE);
              if(i==0)
              {
                if(c) H_mult[c][i][j]->GetYaxis()->SetTitle("#font[12]{M}^{#font[ 12]{h^{+}}}+ #font[ 12]{#delta}");
                else H_mult[c][i][j]->GetYaxis()->SetTitle("#font[12]{M}^{#font[ 12]{h^{-}}}+ #font[ 12]{#delta}");
                H_mult[c][i][j]->GetYaxis()->SetTitleSize(0.08);
              }
              lsys.Draw();
              if(j==3) H_sys[c][i][j]->Draw("SAME3");
              if(!c) axisflagh1=1;
              else axisflagh2=1;
              if(c) c51->Range(0.1,-0.4,0.9,5.);
              else c52->Range(0.1,-0.4,0.9,5.);
            }
            else
            {
              H_mult[c][i][j]->Draw("SAMEP");
              if(j==3) H_sys[c][i][j]->Draw("SAME3");
              H_mult[c][i][j]->GetXaxis()->SetLimits(0.1,0.9);
              H_mult[c][i][j]->SetMinimum(-0.4);
              H_mult[c][i][j]->SetMaximum(5.);
            }
          }
          if(c) c51->Update();
          else c52->Update();
        }

        if(!p_m_empty)
        {
          if(c) c61->cd(i+1);
          else c62->cd(i+1);
          if(P_mult[c][i][j])
          {
            if((!c && !axisflagp1) || (c && !axisflagp2))
            {
              P_mult[c][i][j]->Draw("SAMEPA");
              P_mult[c][i][j]->GetXaxis()->SetLimits(0.1,0.9);
              P_mult[c][i][j]->SetMinimum(-0.4);
              P_mult[c][i][j]->SetMaximum(4.);
              P_mult[c][i][j]->GetXaxis()->SetLabelSize(0.06);
              P_mult[c][i][j]->GetYaxis()->SetLabelSize(0.06);
              P_mult[c][i][j]->SetTitle("");
              if(i>4) gPad->SetBottomMargin(.15);
              if(i==0 || i==5) gPad->SetLeftMargin(.22);
              if(i==8)
              {
                P_mult[c][i][j]->GetXaxis()->SetTitle("#font[ 12]{z}");
                P_mult[c][i][j]->GetXaxis()->SetTitleSize(0.08);
                P_mult[c][i][j]->GetXaxis()->SetTitleOffset(.8);
              }
              P_mult[c][i][j]->GetXaxis()->SetNdivisions(304,kTRUE);
              P_mult[c][i][j]->GetYaxis()->SetNdivisions(304,kTRUE);
              if(i==0)
              {
                if(c) P_mult[c][i][j]->GetYaxis()->SetTitle("#font[12]{M}^{#font[ 12]{#pi^{+}}}+ #font[ 12]{#delta}");
                else P_mult[c][i][j]->GetYaxis()->SetTitle("#font[12]{M}^{#font[ 12]{#pi^{-}}}+ #font[ 12]{#delta}");
                P_mult[c][i][j]->GetYaxis()->SetTitleSize(0.08);
              }
              lsys.Draw();
              if(j==3) P_sys[c][i][j]->Draw("SAME3");
              if(!c) axisflagp1=1;
              else axisflagp2=1;
              if(c) c61->Range(0.1,-0.4,0.9,4.);
              else c62->Range(0.1,-0.4,0.9,4.);
            }
            else
            {
              P_mult[c][i][j]->Draw("SAMEP");
              if(j==3) P_sys[c][i][j]->Draw("SAME3");
              P_mult[c][i][j]->GetXaxis()->SetLimits(0.1,0.9);
              P_mult[c][i][j]->SetMinimum(-0.4);
              P_mult[c][i][j]->SetMaximum(4.);
            }
          }
          if(c) c61->Update();
          else c62->Update();
        }

        if(!k_m_empty)
        {
          if(c) c71->cd(i+1);
          else c72->cd(i+1);
          if(K_mult[c][i][j])
          {
            if((!c && !axisflagk1) || (c && !axisflagk2))
            {
              K_mult[c][i][j]->Draw("SAMEPA");
              K_mult[c][i][j]->GetXaxis()->SetLimits(0.1,0.9);
              K_mult[c][i][j]->SetMinimum(-0.06);
              K_mult[c][i][j]->SetMaximum(0.8);
              K_mult[c][i][j]->GetXaxis()->SetLabelSize(0.06);
              K_mult[c][i][j]->GetYaxis()->SetLabelSize(0.06);
              K_mult[c][i][j]->SetTitle("");
              if(i>4) gPad->SetBottomMargin(.15);
              if(i==0 || i==5) gPad->SetLeftMargin(.22);
              if(i==8)
              {
                K_mult[c][i][j]->GetXaxis()->SetTitle("#font[ 12]{z}");
                K_mult[c][i][j]->GetXaxis()->SetTitleSize(0.08);
                K_mult[c][i][j]->GetXaxis()->SetTitleOffset(.8);
              }
              K_mult[c][i][j]->GetXaxis()->SetNdivisions(304,kTRUE);
              K_mult[c][i][j]->GetYaxis()->SetNdivisions(304,kTRUE);
              if(i==0)
              {
                if(c) K_mult[c][i][j]->GetYaxis()->SetTitle("#font[12]{M}^{#font[ 12]{K^{+}}}+ #font[ 12]{#delta}");
                else K_mult[c][i][j]->GetYaxis()->SetTitle("#font[12]{M}^{#font[ 12]{K^{-}}}+ #font[ 12]{#delta}");
                K_mult[c][i][j]->GetYaxis()->SetTitleSize(0.08);
              }
              lsys.Draw();
              if(j==3) K_sys[c][i][j]->Draw("SAME3");
              if(!c) axisflagk1=1;
              else axisflagk2=1;
              if(c) c71->Range(0.1,-0.06,0.9,0.8);
              else c72->Range(0.1,-0.06,0.9,0.8);
            }
            else
            {
              K_mult[c][i][j]->Draw("SAMEP");
              if(j==3) K_sys[c][i][j]->Draw("SAME3");
              K_mult[c][i][j]->GetXaxis()->SetLimits(0.1,0.9);
              K_mult[c][i][j]->SetMinimum(-0.06);
              K_mult[c][i][j]->SetMaximum(0.8);
            }
          }
          if(c) c71->Update();
          else c72->Update();
        }
        z_range_p[c][i][j].clear();
        z_range_k[c][i][j].clear();
        z_range_h[c][i][j].clear();

        for(int zv=0; zv<4; zv++)
        {
          for(int l=0; l<12; l++)
          {
            z_range_p_z[c][i][j][zv].push_back(z_range[l]);
            z_range_k_z[c][i][j][zv].push_back(z_range[l]);
            z_range_h_z[c][i][j][zv].push_back(z_range[l]);
          }

          for(int k=12; k>0; k--)
          {
            if(!p_z[c][i][j][zv][k-1]) {p_z[c][i][j][zv].erase(p_z[c][i][j][zv].begin()+k-1); p_z_err[c][i][j][zv].erase(p_z_err[c][i][j][zv].begin()+k-1); z_range_p_z[c][i][j][zv].erase(z_range_p_z[c][i][j][zv].begin()+k-1);}
            if(!k_z[c][i][j][zv][k-1]) {k_z[c][i][j][zv].erase(k_z[c][i][j][zv].begin()+k-1); k_z_err[c][i][j][zv].erase(k_z_err[c][i][j][zv].begin()+k-1); z_range_k_z[c][i][j][zv].erase(z_range_k_z[c][i][j][zv].begin()+k-1);}
            if(!h_z[c][i][j][zv][k-1]) {h_z[c][i][j][zv].erase(h_z[c][i][j][zv].begin()+k-1); h_z_err[c][i][j][zv].erase(h_z_err[c][i][j][zv].begin()+k-1); z_range_h_z[c][i][j][zv].erase(z_range_h_z[c][i][j][zv].begin()+k-1);}
          }

          bool p_z_empty = 0;
          bool k_z_empty = 0;
          bool h_z_empty = 0;

          if(!(Int_t(p_z[c][i][j][zv].size()))) p_z_empty = 1;
          if(!(Int_t(k_z[c][i][j][zv].size()))) k_z_empty = 1;
          if(!(Int_t(h_z[c][i][j][zv].size()))) h_z_empty = 1;

          H_zvtx[c][i][j][zv] = new TGraphErrors(Int_t(h_z[c][i][j][zv].size()),&(z_range_h_z[c][i][j][zv][0]),&(h_z[c][i][j][zv][0]),0,&(h_z_err[c][i][j][zv][0]));
          P_zvtx[c][i][j][zv] = new TGraphErrors(Int_t(p_z[c][i][j][zv].size()),&(z_range_p_z[c][i][j][zv][0]),&(p_z[c][i][j][zv][0]),0,&(p_z_err[c][i][j][zv][0]));
          K_zvtx[c][i][j][zv] = new TGraphErrors(Int_t(k_z[c][i][j][zv].size()),&(z_range_k_z[c][i][j][zv][0]),&(k_z[c][i][j][zv][0]),0,&(k_z_err[c][i][j][zv][0]));

          H_zvtx[c][i][j][zv]->SetMarkerColor(fMarkerColorZvtx[zv][c]);
          P_zvtx[c][i][j][zv]->SetMarkerColor(fMarkerColorZvtx[zv][c]);
          K_zvtx[c][i][j][zv]->SetMarkerColor(fMarkerColorZvtx[zv][c]);

          H_zvtx[c][i][j][zv]->SetMarkerSize(2);
          P_zvtx[c][i][j][zv]->SetMarkerSize(2);
          K_zvtx[c][i][j][zv]->SetMarkerSize(2);

          H_zvtx[c][i][j][zv]->SetMarkerStyle(fMarkerStyle[0][c]);
          P_zvtx[c][i][j][zv]->SetMarkerStyle(fMarkerStyle[0][c]);
          K_zvtx[c][i][j][zv]->SetMarkerStyle(fMarkerStyle[0][c]);

          H_zvtx[c][i][j][zv]->SetTitle("");
          P_zvtx[c][i][j][zv]->SetTitle("");
          K_zvtx[c][i][j][zv]->SetTitle("");

          H_zvtx[c][i][j][zv]->GetYaxis()->SetTitle("");
          P_zvtx[c][i][j][zv]->GetYaxis()->SetTitle("");
          K_zvtx[c][i][j][zv]->GetYaxis()->SetTitle("");

          if(!h_z_empty)
          {
            c5->cd(i+1+9*j);
            gPad->SetFillStyle(4000);
            if(H_zvtx[c][i][j][zv])
            {
              if(!c && !zv)
              {
                H_zvtx[c][i][j][zv]->Draw("SAMEPA");
                H_zvtx[c][i][j][zv]->GetXaxis()->SetLimits(0.1,0.9);
                H_zvtx[c][i][j][zv]->SetMinimum(0.);
                H_zvtx[c][i][j][zv]->SetMaximum(4.0);
                H_zvtx[c][i][j][zv]->GetXaxis()->SetLabelSize(0.06);
                H_zvtx[c][i][j][zv]->GetYaxis()->SetLabelSize(0.06);
                H_zvtx[c][i][j][zv]->SetTitle("");
                if(j==5) gPad->SetBottomMargin(.15);
                if(i==0) gPad->SetLeftMargin(.22);
                if(i==8 && j==5)
                {
                  H_zvtx[c][i][j][zv]->GetXaxis()->SetTitle("#font[ 12]{z}");
                  H_zvtx[c][i][j][zv]->GetXaxis()->SetTitleSize(0.08);
                  H_zvtx[c][i][j][zv]->GetXaxis()->SetTitleOffset(.8);
                }
                H_zvtx[c][i][j][zv]->GetXaxis()->SetNdivisions(304,kTRUE);
                H_zvtx[c][i][j][zv]->GetYaxis()->SetNdivisions(304,kTRUE);
                if(i==2 && j==0)
                {
                  H_zvtx[c][i][j][zv]->GetYaxis()->SetTitle("#font[12]{acceptance}^{#font[ 12]{h}}");
                  H_zvtx[c][i][j][zv]->GetYaxis()->SetTitleSize(0.08);
                }
                H_zvtx[c][i][j][zv]->Draw("SAMEP");
                H_zvtx[c][i][j][zv]->GetXaxis()->SetLimits(0.1,0.9);
                H_zvtx[c][i][j][zv]->SetMinimum(0.);
                H_zvtx[c][i][j][zv]->SetMaximum(4.0);
                c5->Range(0.1,0.,0.9,4.0);
              }
              else
              {
                H_zvtx[c][i][j][zv]->Draw("SAMEP");
                H_zvtx[c][i][j][zv]->GetXaxis()->SetLimits(0.1,0.9);
                H_zvtx[c][i][j][zv]->SetMinimum(0.);
                H_zvtx[c][i][j][zv]->SetMaximum(4.0);
              }
            }
            c5->Update();
          }


          if(!p_z_empty)
          {
            c6->cd(i+1+9*j);
            gPad->SetFillStyle(4000);
            if(P_zvtx[c][i][j][zv])
            {
              if(!c && !zv)
              {
                P_zvtx[c][i][j][zv]->Draw("SAMEPA");
                P_zvtx[c][i][j][zv]->GetXaxis()->SetLimits(0.1,0.9);
                P_zvtx[c][i][j][zv]->SetMinimum(0.);
                P_zvtx[c][i][j][zv]->SetMaximum(3.5);
                P_zvtx[c][i][j][zv]->GetXaxis()->SetLabelSize(0.06);
                P_zvtx[c][i][j][zv]->GetYaxis()->SetLabelSize(0.06);
                P_zvtx[c][i][j][zv]->SetTitle("");
                if(j==5) gPad->SetBottomMargin(.15);
                if(i==0) gPad->SetLeftMargin(.22);
                if(i==8 && j==5)
                {
                  P_zvtx[c][i][j][zv]->GetXaxis()->SetTitle("#font[ 12]{z}");
                  P_zvtx[c][i][j][zv]->GetXaxis()->SetTitleSize(0.08);
                  P_zvtx[c][i][j][zv]->GetXaxis()->SetTitleOffset(.8);
                }
                P_zvtx[c][i][j][zv]->GetXaxis()->SetNdivisions(304,kTRUE);
                P_zvtx[c][i][j][zv]->GetYaxis()->SetNdivisions(304,kTRUE);
                if(i==1 && j==0)
                {
                  P_zvtx[c][i][j][zv]->GetYaxis()->SetTitle("#font[12]{acceptance}^{#font[ 12]{h}}");
                  P_zvtx[c][i][j][zv]->GetYaxis()->SetTitleSize(0.08);
                }
                P_zvtx[c][i][j][zv]->Draw("SAMEP");
                P_zvtx[c][i][j][zv]->GetXaxis()->SetLimits(0.1,0.9);
                P_zvtx[c][i][j][zv]->SetMinimum(0.);
                P_zvtx[c][i][j][zv]->SetMaximum(3.5);
                c6->Range(0.1,0.,0.9,3.5);
              }
              else
              {
                P_zvtx[c][i][j][zv]->Draw("SAMEP");
                P_zvtx[c][i][j][zv]->GetXaxis()->SetLimits(0.1,0.9);
                P_zvtx[c][i][j][zv]->SetMinimum(0.);
                P_zvtx[c][i][j][zv]->SetMaximum(3.5);
              }
            }
            c6->Update();
          }

          if(!k_z_empty)
          {
            c7->cd(i+1+9*j);
            gPad->SetFillStyle(4000);
            if(P_zvtx[c][i][j][zv])
            {
              if(!c && !zv)
              {
                K_zvtx[c][i][j][zv]->Draw("SAMEPA");
                K_zvtx[c][i][j][zv]->GetXaxis()->SetLimits(0.1,0.9);
                K_zvtx[c][i][j][zv]->SetMinimum(0.);
                K_zvtx[c][i][j][zv]->SetMaximum(0.8);
                K_zvtx[c][i][j][zv]->GetXaxis()->SetLabelSize(0.06);
                K_zvtx[c][i][j][zv]->GetYaxis()->SetLabelSize(0.06);
                K_zvtx[c][i][j][zv]->SetTitle("");
                if(j==5) gPad->SetBottomMargin(.15);
                if(i==0) gPad->SetLeftMargin(.22);
                if(i==8 && j==5)
                {
                  K_zvtx[c][i][j][zv]->GetXaxis()->SetTitle("#font[ 12]{z}");
                  K_zvtx[c][i][j][zv]->GetXaxis()->SetTitleSize(0.08);
                  K_zvtx[c][i][j][zv]->GetXaxis()->SetTitleOffset(.8);
                }
                K_zvtx[c][i][j][zv]->GetXaxis()->SetNdivisions(304,kTRUE);
                K_zvtx[c][i][j][zv]->GetYaxis()->SetNdivisions(304,kTRUE);
                if(i==1 && j==0)
                {
                  K_zvtx[c][i][j][zv]->GetYaxis()->SetTitle("#font[12]{acceptance}^{#font[ 12]{h}}");
                  K_zvtx[c][i][j][zv]->GetYaxis()->SetTitleSize(0.08);
                }
                K_zvtx[c][i][j][zv]->Draw("SAMEP");
                K_zvtx[c][i][j][zv]->GetXaxis()->SetLimits(0.1,0.9);
                K_zvtx[c][i][j][zv]->SetMinimum(0.);
                K_zvtx[c][i][j][zv]->SetMaximum(0.8);
                c7->Range(0.1,0.,0.9,0.8);
              }
              else
              {
                K_zvtx[c][i][j][zv]->Draw("SAMEP");
                K_zvtx[c][i][j][zv]->GetXaxis()->SetLimits(0.1,0.9);
                K_zvtx[c][i][j][zv]->SetMinimum(0.);
                K_zvtx[c][i][j][zv]->SetMaximum(0.8);
              }
            }
            c7->Update();
          }
          z_range_p_z[c][i][j][zv].clear();
          z_range_k_z[c][i][j][zv].clear();
          z_range_h_z[c][i][j][zv].clear();
        }

        for(int l=0; l<12; l++)
        {
          z_range_p_reldiff[c][i][j].push_back(z_range[l]);
          z_range_k_reldiff[c][i][j].push_back(z_range[l]);
          z_range_h_reldiff[c][i][j].push_back(z_range[l]);
        }

        for(int k=12; k>0; k--)
        {
          if(!p_reldiff[c][i][j][k-1]) {p_reldiff[c][i][j].erase(p_reldiff[c][i][j].begin()+k-1); p_reldiff_err[c][i][j].erase(p_reldiff_err[c][i][j].begin()+k-1); z_range_p_reldiff[c][i][j].erase(z_range_p_reldiff[c][i][j].begin()+k-1);}
          if(!k_reldiff[c][i][j][k-1]) {k_reldiff[c][i][j].erase(k_reldiff[c][i][j].begin()+k-1); k_reldiff_err[c][i][j].erase(k_reldiff_err[c][i][j].begin()+k-1); z_range_k_reldiff[c][i][j].erase(z_range_k_reldiff[c][i][j].begin()+k-1);}
          if(!h_reldiff[c][i][j][k-1]) {h_reldiff[c][i][j].erase(h_reldiff[c][i][j].begin()+k-1); h_reldiff_err[c][i][j].erase(h_reldiff_err[c][i][j].begin()+k-1); z_range_h_reldiff[c][i][j].erase(z_range_h_reldiff[c][i][j].begin()+k-1);}
        }

        bool p_reldiff_empty = 0;
        bool k_reldiff_empty = 0;
        bool h_reldiff_empty = 0;

        if(!(int(p_reldiff[c][i][j].size()))) p_reldiff_empty = 1;
        if(!(int(k_reldiff[c][i][j].size()))) k_reldiff_empty = 1;
        if(!(int(h_reldiff[c][i][j].size()))) h_reldiff_empty = 1;

        P_reldiff[c][i][j] = new TGraphErrors(int(p_reldiff[c][i][j].size()),&(z_range_p_reldiff[c][i][j][0]),&(p_reldiff[c][i][j][0]),0,&(p_reldiff_err[c][i][j][0]));
        K_reldiff[c][i][j] = new TGraphErrors(int(k_reldiff[c][i][j].size()),&(z_range_k_reldiff[c][i][j][0]),&(k_reldiff[c][i][j][0]),0,&(k_reldiff_err[c][i][j][0]));
        H_reldiff[c][i][j] = new TGraphErrors(int(h_reldiff[c][i][j].size()),&(z_range_h_reldiff[c][i][j][0]),&(h_reldiff[c][i][j][0]),0,&(h_reldiff_err[c][i][j][0]));

        if(!c)
        {
          P_reldiff[c][i][j]->SetMarkerColor(fMarkerColor[4]);
          K_reldiff[c][i][j]->SetMarkerColor(fMarkerColor[4]);
          H_reldiff[c][i][j]->SetMarkerColor(fMarkerColor[4]);
        }
        else
        {
          P_reldiff[c][i][j]->SetMarkerColor(fMarkerColor[0]);
          K_reldiff[c][i][j]->SetMarkerColor(fMarkerColor[0]);
          H_reldiff[c][i][j]->SetMarkerColor(fMarkerColor[0]);
        }

        P_reldiff[c][i][j]->SetMarkerSize(2);
        K_reldiff[c][i][j]->SetMarkerSize(2);
        H_reldiff[c][i][j]->SetMarkerSize(2);

        P_reldiff[c][i][j]->SetMarkerStyle(fMarkerStyle[0][c]);
        K_reldiff[c][i][j]->SetMarkerStyle(fMarkerStyle[0][c]);
        H_reldiff[c][i][j]->SetMarkerStyle(fMarkerStyle[0][c]);

        P_reldiff[c][i][j]->GetYaxis()->SetTitle("");
        K_reldiff[c][i][j]->GetYaxis()->SetTitle("");
        H_reldiff[c][i][j]->GetYaxis()->SetTitle("");

        P_reldiff[c][i][j]->GetXaxis()->SetTitle("");
        K_reldiff[c][i][j]->GetXaxis()->SetTitle("");
        H_reldiff[c][i][j]->GetXaxis()->SetTitle("");

        P_reldiff[c][i][j]->SetTitle("");
        K_reldiff[c][i][j]->SetTitle("");
        H_reldiff[c][i][j]->SetTitle("");

        if(!h_reldiff_empty)
        {
          c53->cd(i+1+9*j);
          gPad->SetFillStyle(4000);
          if(H_reldiff[c][i][j])
          {
            if(!c)
            {
              H_reldiff[c][i][j]->Draw("SAMEPA");
              H_reldiff[c][i][j]->GetXaxis()->SetLimits(0.1,0.9);
              H_reldiff[c][i][j]->SetMinimum(-1);
              H_reldiff[c][i][j]->SetMaximum(1);
              H_reldiff[c][i][j]->GetXaxis()->SetLabelSize(0.06);
              H_reldiff[c][i][j]->GetYaxis()->SetLabelSize(0.06);
              H_reldiff[c][i][j]->SetTitle("");
              if(j==5) gPad->SetBottomMargin(.15);
              if(i==0) gPad->SetLeftMargin(.22);
              if(i==8 && j==5)
              {
                H_reldiff[c][i][j]->GetXaxis()->SetTitle("#font[ 12]{z}");
                H_reldiff[c][i][j]->GetXaxis()->SetTitleSize(0.08);
                H_reldiff[c][i][j]->GetXaxis()->SetTitleOffset(.8);
              }
              H_reldiff[c][i][j]->GetXaxis()->SetNdivisions(304,kTRUE);
              H_reldiff[c][i][j]->GetYaxis()->SetNdivisions(304,kTRUE);
              if(i==1 && j==0)
              {
                H_reldiff[c][i][j]->GetYaxis()->SetTitle("#font[12]{M}^{#font[ 12]{h}}_{#font[ 12]{ratio}}");
                H_reldiff[c][i][j]->GetYaxis()->SetTitleSize(0.08);
              }
              H_reldiff[c][i][j]->Draw("SAMEP");
              H_reldiff[c][i][j]->GetXaxis()->SetLimits(0.1,0.9);
              H_reldiff[c][i][j]->SetMinimum(-1);
              H_reldiff[c][i][j]->SetMaximum(1);
              // l1.Draw("SAME");
              // l2.Draw("SAME");
              // l3.Draw("SAME");
              // l4.Draw("SAME");
              c53->Range(0.1,-1,0.9,1);
            }
            else
            {
              H_reldiff[c][i][j]->Draw("SAMEP");
              H_reldiff[c][i][j]->GetXaxis()->SetLimits(0.1,0.9);
              H_reldiff[c][i][j]->SetMinimum(-1);
              H_reldiff[c][i][j]->SetMaximum(1);
            }
          }
          c53->Update();
        }
        if(!p_reldiff_empty)
        {
          c63->cd(i+1+9*j);
          gPad->SetFillStyle(4000);
          if(P_reldiff[c][i][j])
          {
            if(!c)
            {
              P_reldiff[c][i][j]->Draw("SAMEPA");
              P_reldiff[c][i][j]->GetXaxis()->SetLimits(0.1,0.9);
              P_reldiff[c][i][j]->SetMinimum(-1);
              P_reldiff[c][i][j]->SetMaximum(1);
              P_reldiff[c][i][j]->GetXaxis()->SetLabelSize(0.06);
              P_reldiff[c][i][j]->GetYaxis()->SetLabelSize(0.06);
              P_reldiff[c][i][j]->SetTitle("");
              if(j==5) gPad->SetBottomMargin(.15);
              if(i==0) gPad->SetLeftMargin(.22);
              if(i==8 && j==5)
              {
                P_reldiff[c][i][j]->GetXaxis()->SetTitle("#font[ 12]{z}");
                P_reldiff[c][i][j]->GetXaxis()->SetTitleSize(0.08);
                P_reldiff[c][i][j]->GetXaxis()->SetTitleOffset(.8);
              }
              P_reldiff[c][i][j]->GetXaxis()->SetNdivisions(304,kTRUE);
              P_reldiff[c][i][j]->GetYaxis()->SetNdivisions(304,kTRUE);
              if(i==1 && j==0)
              {
                P_reldiff[c][i][j]->GetYaxis()->SetTitle("#font[12]{M}^{#font[ 12]{h}}_{#font[ 12]{ratio}}");
                P_reldiff[c][i][j]->GetYaxis()->SetTitleSize(0.08);
              }
              P_reldiff[c][i][j]->Draw("SAMEP");
              P_reldiff[c][i][j]->GetXaxis()->SetLimits(0.1,0.9);
              P_reldiff[c][i][j]->SetMinimum(-1);
              P_reldiff[c][i][j]->SetMaximum(1);
              // l1.Draw("SAME");
              // l2.Draw("SAME");
              // l3.Draw("SAME");
              // l4.Draw("SAME");
              c63->Range(0.1,-1,0.9,1);
            }
            else
            {
              P_reldiff[c][i][j]->Draw("SAMEP");
              P_reldiff[c][i][j]->GetXaxis()->SetLimits(0.1,0.9);
              P_reldiff[c][i][j]->SetMinimum(-1);
              P_reldiff[c][i][j]->SetMaximum(1);
            }
          }
          c63->Update();
        }
        if(!k_reldiff_empty)
        {
          c73->cd(i+1+9*j);
          gPad->SetFillStyle(4000);
          if(K_reldiff[c][i][j])
          {
            if(!c)
            {
              K_reldiff[c][i][j]->Draw("SAMEPA");
              K_reldiff[c][i][j]->GetXaxis()->SetLimits(0.1,0.9);
              K_reldiff[c][i][j]->SetMinimum(-1);
              K_reldiff[c][i][j]->SetMaximum(1);
              K_reldiff[c][i][j]->GetXaxis()->SetLabelSize(0.06);
              K_reldiff[c][i][j]->GetYaxis()->SetLabelSize(0.06);
              K_reldiff[c][i][j]->SetTitle("");
              if(j==5) gPad->SetBottomMargin(.15);
              if(i==0) gPad->SetLeftMargin(.22);
              if(i==8 && j==5)
              {
                K_reldiff[c][i][j]->GetXaxis()->SetTitle("#font[ 12]{z}");
                K_reldiff[c][i][j]->GetXaxis()->SetTitleSize(0.08);
                K_reldiff[c][i][j]->GetXaxis()->SetTitleOffset(.8);
              }
              K_reldiff[c][i][j]->GetXaxis()->SetNdivisions(304,kTRUE);
              K_reldiff[c][i][j]->GetYaxis()->SetNdivisions(304,kTRUE);
              if(i==1 && j==0)
              {
                K_reldiff[c][i][j]->GetYaxis()->SetTitle("#font[12]{M}^{#font[ 12]{h}}_{#font[ 12]{ratio}}");
                K_reldiff[c][i][j]->GetYaxis()->SetTitleSize(0.08);
              }
              K_reldiff[c][i][j]->Draw("SAMEP");
              K_reldiff[c][i][j]->GetXaxis()->SetLimits(0.1,0.9);
              K_reldiff[c][i][j]->SetMinimum(-1);
              K_reldiff[c][i][j]->SetMaximum(1);
              // l1.Draw("SAME");
              // l2.Draw("SAME");
              // l3.Draw("SAME");
              // l4.Draw("SAME");
              c73->Range(0.1,-1,0.9,1);
            }
            else
            {
              K_reldiff[c][i][j]->Draw("SAMEP");
              K_reldiff[c][i][j]->GetXaxis()->SetLimits(0.1,0.9);
              K_reldiff[c][i][j]->SetMinimum(-1);
              K_reldiff[c][i][j]->SetMaximum(1);
            }
          }
          c73->Update();
        }
        z_range_p_reldiff[c][i][j].clear();
        z_range_k_reldiff[c][i][j].clear();
        z_range_h_reldiff[c][i][j].clear();
      }
    }

    Double_t MultiplicitiesSum[2][2][4];
    MultiplicitiesSum[0][0][0] = 0;
    MultiplicitiesSum[0][1][0] = 0;
    MultiplicitiesSum[0][0][1] = 0;
    MultiplicitiesSum[0][1][1] = 0;
    MultiplicitiesSum[0][0][3] = 0;
    MultiplicitiesSum[0][1][3] = 0;
    MultiplicitiesSum[1][0][0] = 0;
    MultiplicitiesSum[1][1][0] = 0;
    MultiplicitiesSum[1][0][1] = 0;
    MultiplicitiesSum[1][1][1] = 0;
    MultiplicitiesSum[1][0][3] = 0;
    MultiplicitiesSum[1][1][3] = 0;

    if(YMULT == 1) yavg();
    else if(YMULT == 2) yweightedavg();

    for(int c=0; c<2; c++)
    {
      for(int k=0; k<12; k++)
      {
        for(int l=0; l<4; l++)
        {
          if(fMultiplicities_yavg[i][k].tab[c][0][l]<0)
          {
            fMultiplicities_yavg[i][k].tab[c][0][l] = 0 ;
            fMultiplicities_yavg[i][k].tab[c][1][l] = 0 ;
            fMultiplicities_yavg[i][k].tab[c][2][l] = 0 ;
          }
        }

        // cout << c << " " << i << " " << k << " " << fMultiplicities_yavg[i][k].tab[c][0][0] << " " << fMultiplicities_yavg[i][k].tab[c][1][0] << " " << fMultiplicities_yavg[i][k].tab[c][2][0] << endl;

        p_y[c][i].push_back(fMultiplicities_yavg[i][k].tab[c][0][0]);
        k_y[c][i].push_back(fMultiplicities_yavg[i][k].tab[c][0][1]);
        h_y[c][i].push_back(fMultiplicities_yavg[i][k].tab[c][0][3]);
        p_y_err[c][i].push_back(sqrt(fMultiplicities_yavg[i][k].tab[c][1][0]));
        k_y_err[c][i].push_back(sqrt(fMultiplicities_yavg[i][k].tab[c][1][1]));
        h_y_err[c][i].push_back(sqrt(fMultiplicities_yavg[i][k].tab[c][1][3]));
        p_y_sys[c][i].push_back(sqrt(fMultiplicities_yavg[i][k].tab[c][2][0]));
        k_y_sys[c][i].push_back(sqrt(fMultiplicities_yavg[i][k].tab[c][2][1]));
        h_y_sys[c][i].push_back(sqrt(fMultiplicities_yavg[i][k].tab[c][2][3]));

        MultiplicitiesSum[0][c][0] += fMultiplicities_yavg[i][k].tab[c][0][0]*fZ_bin_width[k];
        MultiplicitiesSum[0][c][1] += fMultiplicities_yavg[i][k].tab[c][0][1]*fZ_bin_width[k];
        MultiplicitiesSum[0][c][3] += fMultiplicities_yavg[i][k].tab[c][0][3]*fZ_bin_width[k];
        MultiplicitiesSum[1][c][0] += fMultiplicities_yavg[i][k].tab[c][1][0]*pow(fZ_bin_width[k],2);
        MultiplicitiesSum[1][c][1] += fMultiplicities_yavg[i][k].tab[c][1][1]*pow(fZ_bin_width[k],2);
        MultiplicitiesSum[1][c][3] += fMultiplicities_yavg[i][k].tab[c][1][3]*pow(fZ_bin_width[k],2);

        for(int zv=0; zv<4; zv++)
        {
          p_z_y[c][i][zv].push_back(fMultiplicities_zvtx_yavg[i][k][zv].tab[c][0][0]>0 ? fMultiplicities_zvtx_yavg[i][k][zv].tab[c][0][0] : 0);
          k_z_y[c][i][zv].push_back(fMultiplicities_zvtx_yavg[i][k][zv].tab[c][0][1]>0 ? fMultiplicities_zvtx_yavg[i][k][zv].tab[c][0][1] : 0);
          h_z_y[c][i][zv].push_back(fMultiplicities_zvtx_yavg[i][k][zv].tab[c][0][3]>0 ? fMultiplicities_zvtx_yavg[i][k][zv].tab[c][0][3] : 0);
          p_z_y_err[c][i][zv].push_back(fMultiplicities_zvtx_yavg[i][k][zv].tab[c][1][0] ? sqrt(fMultiplicities_zvtx_yavg[i][k][zv].tab[c][1][0]) : 0);
          k_z_y_err[c][i][zv].push_back(fMultiplicities_zvtx_yavg[i][k][zv].tab[c][1][1] ? sqrt(fMultiplicities_zvtx_yavg[i][k][zv].tab[c][1][1]) : 0);
          h_z_y_err[c][i][zv].push_back(fMultiplicities_zvtx_yavg[i][k][zv].tab[c][1][3] ? sqrt(fMultiplicities_zvtx_yavg[i][k][zv].tab[c][1][3]) : 0);
        }

        p_y_reldiff[c][i].push_back(RelDiff_yavg(c,i,k,0));
        k_y_reldiff[c][i].push_back(RelDiff_yavg(c,i,k,1));
        h_y_reldiff[c][i].push_back(RelDiff_yavg(c,i,k,3));
        p_y_reldiff_err[c][i].push_back(RelDiff_Err_yavg(c,i,k,0));
        k_y_reldiff_err[c][i].push_back(RelDiff_Err_yavg(c,i,k,1));
        h_y_reldiff_err[c][i].push_back(RelDiff_Err_yavg(c,i,k,3));
      }

      for(int l=0; l<12; l++)
      {
        z_range_p_y[c][i].push_back(z_range[l]);
        z_range_k_y[c][i].push_back(z_range[l]);
        z_range_h_y[c][i].push_back(z_range[l]);
      }

      for(int k=12; k>0; k--)
      {
        if(!p_y[c][i][k-1]) {p_y[c][i].erase(p_y[c][i].begin()+k-1); p_y_err[c][i].erase(p_y_err[c][i].begin()+k-1); p_y_sys[c][i].erase(p_y_sys[c][i].begin()+k-1); z_range_p_y[c][i].erase(z_range_p_y[c][i].begin()+k-1);}
        if(!k_y[c][i][k-1]) {k_y[c][i].erase(k_y[c][i].begin()+k-1); k_y_err[c][i].erase(k_y_err[c][i].begin()+k-1); k_y_sys[c][i].erase(k_y_sys[c][i].begin()+k-1); z_range_k_y[c][i].erase(z_range_k_y[c][i].begin()+k-1);}
        if(!h_y[c][i][k-1]) {h_y[c][i].erase(h_y[c][i].begin()+k-1); h_y_err[c][i].erase(h_y_err[c][i].begin()+k-1); h_y_sys[c][i].erase(h_y_sys[c][i].begin()+k-1); z_range_h_y[c][i].erase(z_range_h_y[c][i].begin()+k-1);}
      }

      bool p_y_empty = 0;
      bool k_y_empty = 0;
      bool h_y_empty = 0;

      if(!(Int_t(p_y[c][i].size()))) p_y_empty = 1;
      if(!(Int_t(k_y[c][i].size()))) k_y_empty = 1;
      if(!(Int_t(h_y[c][i].size()))) h_y_empty = 1;

      H_y[c][i] = new TGraphErrors(Int_t(h_y[c][i].size()),&(z_range_h_y[c][i][0]),&(h_y[c][i][0]),0,&(h_y_err[c][i][0]));
      P_y[c][i] = new TGraphErrors(Int_t(p_y[c][i].size()),&(z_range_p_y[c][i][0]),&(p_y[c][i][0]),0,&(p_y_err[c][i][0]));
      K_y[c][i] = new TGraphErrors(Int_t(k_y[c][i].size()),&(z_range_k_y[c][i][0]),&(k_y[c][i][0]),0,&(k_y_err[c][i][0]));
      if(!c)
      {
        H_ysys[c][i] = new TGraphAsymmErrors(Int_t(h_y[c][i].size()),&(z_range_h_y[c][i][0]),&h_yoffset2[0], &errorx[0], &errorx[0], 0, &(h_y_sys[c][i][0]));
        P_ysys[c][i] = new TGraphAsymmErrors(Int_t(p_y[c][i].size()),&(z_range_p_y[c][i][0]),&p_yoffset2[0], &errorx[0], &errorx[0], 0, &(p_y_sys[c][i][0]));
        K_ysys[c][i] = new TGraphAsymmErrors(Int_t(k_y[c][i].size()),&(z_range_k_y[c][i][0]),&k_yoffset2[0], &errorx[0], &errorx[0], 0, &(k_y_sys[c][i][0]));
      }
      else
      {
        H_ysys[c][i] = new TGraphAsymmErrors(Int_t(h_y[c][i].size()),&(z_range_h_y[c][i][0]),&h_yoffset[0], &errorx[0], &errorx[0], 0, &(h_y_sys[c][i][0]));
        P_ysys[c][i] = new TGraphAsymmErrors(Int_t(p_y[c][i].size()),&(z_range_p_y[c][i][0]),&p_yoffset[0], &errorx[0], &errorx[0], 0, &(p_y_sys[c][i][0]));
        K_ysys[c][i] = new TGraphAsymmErrors(Int_t(k_y[c][i].size()),&(z_range_k_y[c][i][0]),&k_yoffset[0], &errorx[0], &errorx[0], 0, &(k_y_sys[c][i][0]));
      }

      if(!c)
      {
        H_y[c][i]->SetMarkerColor(fMarkerColor[4]);
        P_y[c][i]->SetMarkerColor(fMarkerColor[4]);
        K_y[c][i]->SetMarkerColor(fMarkerColor[4]);
        H_ysys[c][i]->SetFillColor(fMarkerColor[4]);
        P_ysys[c][i]->SetFillColor(fMarkerColor[4]);
        K_ysys[c][i]->SetFillColor(fMarkerColor[4]);
      }
      else
      {
        H_y[c][i]->SetMarkerColor(fMarkerColor[0]);
        P_y[c][i]->SetMarkerColor(fMarkerColor[0]);
        K_y[c][i]->SetMarkerColor(fMarkerColor[0]);
        H_ysys[c][i]->SetFillColor(fMarkerColor[0]);
        P_ysys[c][i]->SetFillColor(fMarkerColor[0]);
        K_ysys[c][i]->SetFillColor(fMarkerColor[0]);
      }

      H_y[c][i]->SetMarkerSize(2);
      P_y[c][i]->SetMarkerSize(2);
      K_y[c][i]->SetMarkerSize(2);

      H_y[c][i]->SetMarkerStyle(fMarkerStyle[0][c]);
      P_y[c][i]->SetMarkerStyle(fMarkerStyle[0][c]);
      K_y[c][i]->SetMarkerStyle(fMarkerStyle[0][c]);

      H_y[c][i]->SetTitle("");
      P_y[c][i]->SetTitle("");
      K_y[c][i]->SetTitle("");

      H_y[c][i]->GetYaxis()->SetTitle("");
      P_y[c][i]->GetYaxis()->SetTitle("");
      K_y[c][i]->GetYaxis()->SetTitle("");

      if(!h_y_empty)
      {
        c8->cd(i+1);
        gPad->SetFillStyle(4000);
        if(H_y[c][i])
        {
          if(!c)
          {
            H_y[c][i]->Draw("SAMEPA");
            H_y[c][i]->GetXaxis()->SetLimits(0.1,0.9);
            H_y[c][i]->SetMinimum(-0.5);
            H_y[c][i]->SetMaximum(4.);
            H_y[c][i]->GetXaxis()->SetLabelSize(0.06);
            H_y[c][i]->GetYaxis()->SetLabelSize(0.06);
            H_y[c][i]->SetTitle("");
            if(i>4) gPad->SetBottomMargin(.15);
            if(i==0 || i==5) gPad->SetLeftMargin(.22);
            if(i==8)
            {
              H_y[c][i]->GetXaxis()->SetTitle("#font[ 12]{z}");
              H_y[c][i]->GetXaxis()->SetTitleSize(0.08);
              H_y[c][i]->GetXaxis()->SetTitleOffset(.8);
            }
            H_y[c][i]->GetXaxis()->SetNdivisions(304,kTRUE);
            H_y[c][i]->GetYaxis()->SetNdivisions(304,kTRUE);
            if(i==0)
            {
              H_y[c][i]->GetYaxis()->SetTitle("#font[12]{M}^{#font[ 12]{h}}");
              H_y[c][i]->GetYaxis()->SetTitleSize(0.08);
            }
            lsys.Draw();
            H_y[c][i]->Draw("SAMEP");
            H_ysys[c][i]->Draw("SAME3");
            H_y[c][i]->GetXaxis()->SetLimits(0.1,0.9);
            H_y[c][i]->SetMinimum(-0.5);
            H_y[c][i]->SetMaximum(4.);
            c8->Range(0.1,-0.5,0.9,4.);
          }
          else
          {
            H_y[c][i]->Draw("SAMEP");
            H_ysys[c][i]->Draw("SAME3");
            H_y[c][i]->GetXaxis()->SetLimits(0.1,0.9);
            H_y[c][i]->SetMinimum(-0.5);
            H_y[c][i]->SetMaximum(4.);
          }
        }
        c8->Update();
      }

      if(!p_y_empty)
      {
        c9->cd(i+1);
        if(P_y[c][i])
        {
          if(!c)
          {
            P_y[c][i]->Draw("SAMEPA");
            P_y[c][i]->GetXaxis()->SetLimits(0.1,0.9);
            P_y[c][i]->SetMinimum(-0.5);
            P_y[c][i]->SetMaximum(3.5);
            P_y[c][i]->GetXaxis()->SetLabelSize(0.06);
            P_y[c][i]->GetYaxis()->SetLabelSize(0.06);
            P_y[c][i]->SetTitle("");
            if(i>4) gPad->SetBottomMargin(.15);
            if(i==0 || i==5) gPad->SetLeftMargin(.22);
            if(i==8)
            {
              P_y[c][i]->GetXaxis()->SetTitle("#font[ 12]{z}");
              P_y[c][i]->GetXaxis()->SetTitleSize(0.08);
              P_y[c][i]->GetXaxis()->SetTitleOffset(.8);
            }
            P_y[c][i]->GetXaxis()->SetNdivisions(304,kTRUE);
            P_y[c][i]->GetYaxis()->SetNdivisions(304,kTRUE);
            if(i==0)
            {
              P_y[c][i]->GetYaxis()->SetTitle("#font[12]{M}^{#font[ 12]{#pi}}");
              P_y[c][i]->GetYaxis()->SetTitleSize(0.08);
            }
            lsys.Draw();
            P_y[c][i]->Draw("SAMEP");
            P_ysys[c][i]->Draw("SAME3");
            P_y[c][i]->GetXaxis()->SetLimits(0.1,0.9);
            P_y[c][i]->SetMinimum(-0.5);
            P_y[c][i]->SetMaximum(3.5);
            c9->Range(0.1,-0.5,0.9,3.5);
          }
          else
          {
            P_y[c][i]->Draw("SAMEP");
            P_ysys[c][i]->Draw("SAME3");
            P_y[c][i]->GetXaxis()->SetLimits(0.1,0.9);
            P_y[c][i]->SetMinimum(-0.5);
            P_y[c][i]->SetMaximum(3.5);
          }
        }
        c9->Update();
      }

      if(!k_y_empty)
      {
        c10->cd(i+1);
        if(K_y[c][i])
        {
          if(!c)
          {
            K_y[c][i]->Draw("SAMEPA");
            K_y[c][i]->GetXaxis()->SetLimits(0.1,0.9);
            K_y[c][i]->SetMinimum(-0.1);
            K_y[c][i]->SetMaximum(0.8);
            K_y[c][i]->GetXaxis()->SetLabelSize(0.06);
            K_y[c][i]->GetYaxis()->SetLabelSize(0.06);
            K_y[c][i]->SetTitle("");
            if(i>4) gPad->SetBottomMargin(.15);
            if(i==0 || i==5) gPad->SetLeftMargin(.22);
            if(i==8)
            {
              K_y[c][i]->GetXaxis()->SetTitle("#font[ 12]{z}");
              K_y[c][i]->GetXaxis()->SetTitleSize(0.08);
              K_y[c][i]->GetXaxis()->SetTitleOffset(.8);
            }
            K_y[c][i]->GetXaxis()->SetNdivisions(304,kTRUE);
            K_y[c][i]->GetYaxis()->SetNdivisions(304,kTRUE);
            if(i==0)
            {
              K_y[c][i]->GetYaxis()->SetTitle("#font[12]{M}^{#font[ 12]{K}}");
              K_y[c][i]->GetYaxis()->SetTitleSize(0.08);
            }
            lsys.Draw();
            K_y[c][i]->Draw("SAMEP");
            K_ysys[c][i]->Draw("SAME3");
            K_y[c][i]->GetXaxis()->SetLimits(0.1,0.9);
            K_y[c][i]->SetMinimum(-0.1);
            K_y[c][i]->SetMaximum(0.8);
            c10->Range(0.1,-0.1,0.9,0.8);
          }
          else
          {
            K_y[c][i]->Draw("SAMEP");
            K_ysys[c][i]->Draw("SAME3");
            K_y[c][i]->GetXaxis()->SetLimits(0.1,0.9);
            K_y[c][i]->SetMinimum(-0.1);
            K_y[c][i]->SetMaximum(0.8);
          }
        }
        c10->Update();
      }
      z_range_p_y[c][i].clear();
      z_range_k_y[c][i].clear();
      z_range_h_y[c][i].clear();

      for(int zv=0; zv<4; zv++)
      {
        for(int l=0; l<12; l++)
        {
          z_range_p_y_z[c][i][zv].push_back(z_range[l]);
          z_range_k_y_z[c][i][zv].push_back(z_range[l]);
          z_range_h_y_z[c][i][zv].push_back(z_range[l]);
        }

        for(int k=12; k>0; k--)
        {
          if(!p_z[c][i][zv][k-1]) {p_y_z[c][i][zv].erase(p_y_z[c][i][zv].begin()+k-1); p_y_z_err[c][i][zv].erase(p_y_z_err[c][i][zv].begin()+k-1); z_range_p_y_z[c][i][zv].erase(z_range_p_y_z[c][i][zv].begin()+k-1);}
          if(!k_z[c][i][zv][k-1]) {k_y_z[c][i][zv].erase(k_y_z[c][i][zv].begin()+k-1); k_y_z_err[c][i][zv].erase(k_y_z_err[c][i][zv].begin()+k-1); z_range_k_y_z[c][i][zv].erase(z_range_k_y_z[c][i][zv].begin()+k-1);}
          if(!h_z[c][i][zv][k-1]) {h_y_z[c][i][zv].erase(h_y_z[c][i][zv].begin()+k-1); h_y_z_err[c][i][zv].erase(h_y_z_err[c][i][zv].begin()+k-1); z_range_h_y_z[c][i][zv].erase(z_range_h_y_z[c][i][zv].begin()+k-1);}
        }

        bool p_y_z_empty = 0;
        bool k_y_z_empty = 0;
        bool h_y_z_empty = 0;

        if(!(Int_t(p_y_z[c][i][zv].size()))) p_y_z_empty = 1;
        if(!(Int_t(k_y_z[c][i][zv].size()))) k_y_z_empty = 1;
        if(!(Int_t(h_y_z[c][i][zv].size()))) h_y_z_empty = 1;

        H_y_zvtx[c][i][zv] = new TGraphErrors(Int_t(h_y_z[c][i][zv].size()),&(z_range_h_y_z[c][i][zv][0]),&(h_y_z[c][i][zv][0]),0,&(h_y_z_err[c][i][zv][0]));
        P_y_zvtx[c][i][zv] = new TGraphErrors(Int_t(p_y_z[c][i][zv].size()),&(z_range_p_y_z[c][i][zv][0]),&(p_y_z[c][i][zv][0]),0,&(p_y_z_err[c][i][zv][0]));
        K_y_zvtx[c][i][zv] = new TGraphErrors(Int_t(k_y_z[c][i][zv].size()),&(z_range_k_y_z[c][i][zv][0]),&(k_y_z[c][i][zv][0]),0,&(k_y_z_err[c][i][zv][0]));

        H_y_zvtx[c][i][zv]->SetMarkerColor(fMarkerColorZvtx[zv][c]);
        P_y_zvtx[c][i][zv]->SetMarkerColor(fMarkerColorZvtx[zv][c]);
        K_y_zvtx[c][i][zv]->SetMarkerColor(fMarkerColorZvtx[zv][c]);

        H_y_zvtx[c][i][zv]->SetMarkerSize(2);
        P_y_zvtx[c][i][zv]->SetMarkerSize(2);
        K_y_zvtx[c][i][zv]->SetMarkerSize(2);

        H_y_zvtx[c][i][zv]->SetMarkerStyle(fMarkerStyle[0][c]);
        P_y_zvtx[c][i][zv]->SetMarkerStyle(fMarkerStyle[0][c]);
        K_y_zvtx[c][i][zv]->SetMarkerStyle(fMarkerStyle[0][c]);

        H_y_zvtx[c][i][zv]->SetTitle("");
        P_y_zvtx[c][i][zv]->SetTitle("");
        K_y_zvtx[c][i][zv]->SetTitle("");

        H_y_zvtx[c][i][zv]->GetYaxis()->SetTitle("");
        P_y_zvtx[c][i][zv]->GetYaxis()->SetTitle("");
        K_y_zvtx[c][i][zv]->GetYaxis()->SetTitle("");

        if(!h_y_z_empty)
        {
          c511->cd(i+1);
          gPad->SetFillStyle(4000);
          if(H_y_zvtx[c][i][zv])
          {
            if(!c && !zv)
            {
              H_y_zvtx[c][i][zv]->Draw("SAMEPA");
              H_y_zvtx[c][i][zv]->GetXaxis()->SetLimits(0.1,0.9);
              H_y_zvtx[c][i][zv]->SetMinimum(0.);
              H_y_zvtx[c][i][zv]->SetMaximum(4.0);
              H_y_zvtx[c][i][zv]->GetXaxis()->SetLabelSize(0.06);
              H_y_zvtx[c][i][zv]->GetYaxis()->SetLabelSize(0.06);
              H_y_zvtx[c][i][zv]->SetTitle("");
              if(i>4) gPad->SetBottomMargin(.15);
              if(i==0 || i==5) gPad->SetLeftMargin(.22);
              if(i==8)
              {
                H_y_zvtx[c][i][zv]->GetXaxis()->SetTitle("#font[ 12]{z}");
                H_y_zvtx[c][i][zv]->GetXaxis()->SetTitleSize(0.08);
                H_y_zvtx[c][i][zv]->GetXaxis()->SetTitleOffset(.8);
              }
              H_y_zvtx[c][i][zv]->GetXaxis()->SetNdivisions(304,kTRUE);
              H_y_zvtx[c][i][zv]->GetYaxis()->SetNdivisions(304,kTRUE);
              if(i==0)
              {
                H_y_zvtx[c][i][zv]->GetYaxis()->SetTitle("#font[12]{acceptance}^{#font[ 12]{h}}");
                H_y_zvtx[c][i][zv]->GetYaxis()->SetTitleSize(0.08);
              }
              H_y_zvtx[c][i][zv]->Draw("SAMEP");
              H_y_zvtx[c][i][zv]->GetXaxis()->SetLimits(0.1,0.9);
              H_y_zvtx[c][i][zv]->SetMinimum(0.);
              H_y_zvtx[c][i][zv]->SetMaximum(4.0);
              c511->Range(0.1,0.,0.9,4.0);
            }
            else
            {
              H_y_zvtx[c][i][zv]->Draw("SAMEP");
              H_y_zvtx[c][i][zv]->GetXaxis()->SetLimits(0.1,0.9);
              H_y_zvtx[c][i][zv]->SetMinimum(0.);
              H_y_zvtx[c][i][zv]->SetMaximum(4.0);
            }
          }
          c511->Update();
        }


        if(!p_y_z_empty)
        {
          c611->cd(i+1);
          gPad->SetFillStyle(4000);
          if(P_y_zvtx[c][i][zv])
          {
            if(!c && !zv)
            {
              P_y_zvtx[c][i][zv]->Draw("SAMEPA");
              P_y_zvtx[c][i][zv]->GetXaxis()->SetLimits(0.1,0.9);
              P_y_zvtx[c][i][zv]->SetMinimum(0.);
              P_y_zvtx[c][i][zv]->SetMaximum(3.5);
              P_y_zvtx[c][i][zv]->GetXaxis()->SetLabelSize(0.06);
              P_y_zvtx[c][i][zv]->GetYaxis()->SetLabelSize(0.06);
              P_y_zvtx[c][i][zv]->SetTitle("");
              if(i>4) gPad->SetBottomMargin(.15);
              if(i==0 || i==5) gPad->SetLeftMargin(.22);
              if(i==8)
              {
                P_y_zvtx[c][i][zv]->GetXaxis()->SetTitle("#font[ 12]{z}");
                P_y_zvtx[c][i][zv]->GetXaxis()->SetTitleSize(0.08);
                P_y_zvtx[c][i][zv]->GetXaxis()->SetTitleOffset(.8);
              }
              P_y_zvtx[c][i][zv]->GetXaxis()->SetNdivisions(304,kTRUE);
              P_y_zvtx[c][i][zv]->GetYaxis()->SetNdivisions(304,kTRUE);
              if(i==0)
              {
                P_y_zvtx[c][i][zv]->GetYaxis()->SetTitle("#font[12]{acceptance}^{#font[ 12]{h}}");
                P_y_zvtx[c][i][zv]->GetYaxis()->SetTitleSize(0.08);
              }
              P_y_zvtx[c][i][zv]->Draw("SAMEP");
              P_y_zvtx[c][i][zv]->GetXaxis()->SetLimits(0.1,0.9);
              P_y_zvtx[c][i][zv]->SetMinimum(0.);
              P_y_zvtx[c][i][zv]->SetMaximum(3.5);
              c611->Range(0.1,0.,0.9,3.5);
            }
            else
            {
              P_y_zvtx[c][i][zv]->Draw("SAMEP");
              P_y_zvtx[c][i][zv]->GetXaxis()->SetLimits(0.1,0.9);
              P_y_zvtx[c][i][zv]->SetMinimum(0.);
              P_y_zvtx[c][i][zv]->SetMaximum(3.5);
            }
          }
          c611->Update();
        }

        if(!k_y_z_empty)
        {
          c711->cd(i+1);
          gPad->SetFillStyle(4000);
          if(P_y_zvtx[c][i][zv])
          {
            if(!c && !zv)
            {
              K_y_zvtx[c][i][zv]->Draw("SAMEPA");
              K_y_zvtx[c][i][zv]->GetXaxis()->SetLimits(0.1,0.9);
              K_y_zvtx[c][i][zv]->SetMinimum(0.);
              K_y_zvtx[c][i][zv]->SetMaximum(0.8);
              K_y_zvtx[c][i][zv]->GetXaxis()->SetLabelSize(0.06);
              K_y_zvtx[c][i][zv]->GetYaxis()->SetLabelSize(0.06);
              K_y_zvtx[c][i][zv]->SetTitle("");
              if(i>4) gPad->SetBottomMargin(.15);
              if(i==0 || i==5) gPad->SetLeftMargin(.22);
              if(i==8)
              {
                K_y_zvtx[c][i][zv]->GetXaxis()->SetTitle("#font[ 12]{z}");
                K_y_zvtx[c][i][zv]->GetXaxis()->SetTitleSize(0.08);
                K_y_zvtx[c][i][zv]->GetXaxis()->SetTitleOffset(.8);
              }
              K_y_zvtx[c][i][zv]->GetXaxis()->SetNdivisions(304,kTRUE);
              K_y_zvtx[c][i][zv]->GetYaxis()->SetNdivisions(304,kTRUE);
              if(i==0)
              {
                K_y_zvtx[c][i][zv]->GetYaxis()->SetTitle("#font[12]{acceptance}^{#font[ 12]{h}}");
                K_y_zvtx[c][i][zv]->GetYaxis()->SetTitleSize(0.08);
              }
              K_y_zvtx[c][i][zv]->Draw("SAMEP");
              K_y_zvtx[c][i][zv]->GetXaxis()->SetLimits(0.1,0.9);
              K_y_zvtx[c][i][zv]->SetMinimum(0.);
              K_y_zvtx[c][i][zv]->SetMaximum(0.8);
              c711->Range(0.1,0.,0.9,0.8);
            }
            else
            {
              K_y_zvtx[c][i][zv]->Draw("SAMEP");
              K_y_zvtx[c][i][zv]->GetXaxis()->SetLimits(0.1,0.9);
              K_y_zvtx[c][i][zv]->SetMinimum(0.);
              K_y_zvtx[c][i][zv]->SetMaximum(0.8);
            }
          }
          c711->Update();
        }
        z_range_y_z[c][i][zv].clear();
        z_range_y_z[c][i][zv].clear();
        z_range_y_z[c][i][zv].clear();
      }

      for(int l=0; l<12; l++)
      {
        z_range_p_y_reldiff[c][i].push_back(z_range[l]);
        z_range_k_y_reldiff[c][i].push_back(z_range[l]);
        z_range_h_y_reldiff[c][i].push_back(z_range[l]);
      }

      for(int k=12; k>0; k--)
      {
        if(!p_y_reldiff[c][i][k-1]) {p_y_reldiff[c][i].erase(p_y_reldiff[c][i].begin()+k-1); p_y_reldiff_err[c][i].erase(p_y_reldiff_err[c][i].begin()+k-1); z_range_p_y_reldiff[c][i].erase(z_range_p_y_reldiff[c][i].begin()+k-1);}
        if(!k_y_reldiff[c][i][k-1]) {k_y_reldiff[c][i].erase(k_y_reldiff[c][i].begin()+k-1); k_y_reldiff_err[c][i].erase(k_y_reldiff_err[c][i].begin()+k-1); z_range_k_y_reldiff[c][i].erase(z_range_k_y_reldiff[c][i].begin()+k-1);}
        if(!h_y_reldiff[c][i][k-1]) {h_y_reldiff[c][i].erase(h_y_reldiff[c][i].begin()+k-1); h_y_reldiff_err[c][i].erase(h_y_reldiff_err[c][i].begin()+k-1); z_range_h_y_reldiff[c][i].erase(z_range_h_y_reldiff[c][i].begin()+k-1);}
      }

      bool p_y_reldiff_empty = 0;
      bool k_y_reldiff_empty = 0;
      bool h_y_reldiff_empty = 0;

      if(!(int(p_y_reldiff[c][i].size()))) p_y_reldiff_empty = 1;
      if(!(int(k_y_reldiff[c][i].size()))) k_y_reldiff_empty = 1;
      if(!(int(h_y_reldiff[c][i].size()))) h_y_reldiff_empty = 1;

      P_y_reldiff[c][i] = new TGraphErrors(int(p_y_reldiff[c][i].size()),&(z_range_p_y_reldiff[c][i][0]),&(p_y_reldiff[c][i][0]),0,&(p_y_reldiff_err[c][i][0]));
      K_y_reldiff[c][i] = new TGraphErrors(int(k_y_reldiff[c][i].size()),&(z_range_k_y_reldiff[c][i][0]),&(k_y_reldiff[c][i][0]),0,&(k_y_reldiff_err[c][i][0]));
      H_y_reldiff[c][i] = new TGraphErrors(int(h_y_reldiff[c][i].size()),&(z_range_h_y_reldiff[c][i][0]),&(h_y_reldiff[c][i][0]),0,&(h_y_reldiff_err[c][i][0]));

      if(!c)
      {
        P_y_reldiff[c][i]->SetMarkerColor(fMarkerColor[4]);
        K_y_reldiff[c][i]->SetMarkerColor(fMarkerColor[4]);
        H_y_reldiff[c][i]->SetMarkerColor(fMarkerColor[4]);
      }
      else
      {
        P_y_reldiff[c][i]->SetMarkerColor(fMarkerColor[0]);
        K_y_reldiff[c][i]->SetMarkerColor(fMarkerColor[0]);
        H_y_reldiff[c][i]->SetMarkerColor(fMarkerColor[0]);
      }

      P_y_reldiff[c][i]->SetMarkerSize(2);
      K_y_reldiff[c][i]->SetMarkerSize(2);
      H_y_reldiff[c][i]->SetMarkerSize(2);

      P_y_reldiff[c][i]->SetMarkerStyle(fMarkerStyle[0][c]);
      K_y_reldiff[c][i]->SetMarkerStyle(fMarkerStyle[0][c]);
      H_y_reldiff[c][i]->SetMarkerStyle(fMarkerStyle[0][c]);

      P_y_reldiff[c][i]->GetYaxis()->SetTitle("");
      K_y_reldiff[c][i]->GetYaxis()->SetTitle("");
      H_y_reldiff[c][i]->GetYaxis()->SetTitle("");

      P_y_reldiff[c][i]->GetXaxis()->SetTitle("");
      K_y_reldiff[c][i]->GetXaxis()->SetTitle("");
      H_y_reldiff[c][i]->GetXaxis()->SetTitle("");

      P_y_reldiff[c][i]->SetTitle("");
      K_y_reldiff[c][i]->SetTitle("");
      H_y_reldiff[c][i]->SetTitle("");

      if(!h_y_reldiff_empty)
      {
        c531->cd(i+1+9*j);
        gPad->SetFillStyle(4000);
        if(H_y_reldiff[c][i])
        {
          if(!c)
          {
            H_y_reldiff[c][i]->Draw("SAMEPA");
            H_y_reldiff[c][i]->GetXaxis()->SetLimits(0.1,0.9);
            H_y_reldiff[c][i]->SetMinimum(-1);
            H_y_reldiff[c][i]->SetMaximum(1);
            H_y_reldiff[c][i]->GetXaxis()->SetLabelSize(0.06);
            H_y_reldiff[c][i]->GetYaxis()->SetLabelSize(0.06);
            H_y_reldiff[c][i]->SetTitle("");
            if(i>4) gPad->SetBottomMargin(.15);
            if(i==0 || i==5) gPad->SetLeftMargin(.22);
            if(i==8)
            {
              H_y_reldiff[c][i]->GetXaxis()->SetTitle("#font[ 12]{z}");
              H_y_reldiff[c][i]->GetXaxis()->SetTitleSize(0.08);
              H_y_reldiff[c][i]->GetXaxis()->SetTitleOffset(.8);
            }
            H_y_reldiff[c][i]->GetXaxis()->SetNdivisions(304,kTRUE);
            H_y_reldiff[c][i]->GetYaxis()->SetNdivisions(304,kTRUE);
            if(i==1 && j==0)
            {
              H_y_reldiff[c][i]->GetYaxis()->SetTitle("#font[12]{M}^{#font[ 12]{h}}_{#font[ 12]{ratio}}");
              H_y_reldiff[c][i]->GetYaxis()->SetTitleSize(0.08);
            }
            H_y_reldiff[c][i]->Draw("SAMEP");
            H_y_reldiff[c][i]->GetXaxis()->SetLimits(0.1,0.9);
            H_y_reldiff[c][i]->SetMinimum(-1);
            H_y_reldiff[c][i]->SetMaximum(1);
            // l1.Draw("SAME");
            // l2.Draw("SAME");
            // l3.Draw("SAME");
            // l4.Draw("SAME");
            c531->Range(0.1,-1,0.9,1);
          }
          else
          {
            H_y_reldiff[c][i]->Draw("SAMEP");
            H_y_reldiff[c][i]->GetXaxis()->SetLimits(0.1,0.9);
            H_y_reldiff[c][i]->SetMinimum(-1);
            H_y_reldiff[c][i]->SetMaximum(1);
          }
        }
        c531->Update();
      }
      if(!p_y_reldiff_empty)
      {
        c631->cd(i+1);
        gPad->SetFillStyle(4000);
        if(P_y_reldiff[c][i])
        {
          if(!c)
          {
            P_y_reldiff[c][i]->Draw("SAMEPA");
            P_y_reldiff[c][i]->GetXaxis()->SetLimits(0.1,0.9);
            P_y_reldiff[c][i]->SetMinimum(-1);
            P_y_reldiff[c][i]->SetMaximum(1);
            P_y_reldiff[c][i]->GetXaxis()->SetLabelSize(0.06);
            P_y_reldiff[c][i]->GetYaxis()->SetLabelSize(0.06);
            P_y_reldiff[c][i]->SetTitle("");
            if(i>4) gPad->SetBottomMargin(.15);
            if(i==0 || i==5) gPad->SetLeftMargin(.22);
            if(i==8)
            {
              P_reldiff[c][i]->GetXaxis()->SetTitle("#font[ 12]{z}");
              P_reldiff[c][i]->GetXaxis()->SetTitleSize(0.08);
              P_reldiff[c][i]->GetXaxis()->SetTitleOffset(.8);
            }
            P_y_reldiff[c][i]->GetXaxis()->SetNdivisions(304,kTRUE);
            P_y_reldiff[c][i]->GetYaxis()->SetNdivisions(304,kTRUE);
            if(i==0)
            {
              P_reldiff[c][i]->GetYaxis()->SetTitle("#font[12]{M}^{#font[ 12]{h}}_{#font[ 12]{ratio}}");
              P_reldiff[c][i]->GetYaxis()->SetTitleSize(0.08);
            }
            P_y_reldiff[c][i]->Draw("SAMEP");
            P_y_reldiff[c][i]->GetXaxis()->SetLimits(0.1,0.9);
            P_y_reldiff[c][i]->SetMinimum(-1);
            P_y_reldiff[c][i]->SetMaximum(1);
            // l1.Draw("SAME");
            // l2.Draw("SAME");
            // l3.Draw("SAME");
            // l4.Draw("SAME");
            c631->Range(0.1,-1,0.9,1);
          }
          else
          {
            P_y_reldiff[c][i]->Draw("SAMEP");
            P_y_reldiff[c][i]->GetXaxis()->SetLimits(0.1,0.9);
            P_y_reldiff[c][i]->SetMinimum(-1);
            P_y_reldiff[c][i]->SetMaximum(1);
          }
        }
        c631->Update();
      }
      if(!k_y_reldiff_empty)
      {
        c731->cd(i+1);
        gPad->SetFillStyle(4000);
        if(K_y_reldiff[c][i])
        {
          if(!c)
          {
            K_y_reldiff[c][i]->Draw("SAMEPA");
            K_y_reldiff[c][i]->GetXaxis()->SetLimits(0.1,0.9);
            K_y_reldiff[c][i]->SetMinimum(-1);
            K_y_reldiff[c][i]->SetMaximum(1);
            K_y_reldiff[c][i]->GetXaxis()->SetLabelSize(0.06);
            K_y_reldiff[c][i]->GetYaxis()->SetLabelSize(0.06);
            K_y_reldiff[c][i]->SetTitle("");
            if(i>4) gPad->SetBottomMargin(.15);
            if(i==0 || i==5) gPad->SetLeftMargin(.22);
            if(i==8)
            {
              K_y_reldiff[c][i]->GetXaxis()->SetTitle("#font[ 12]{z}");
              K_y_reldiff[c][i]->GetXaxis()->SetTitleSize(0.08);
              K_y_reldiff[c][i]->GetXaxis()->SetTitleOffset(.8);
            }
            K_y_reldiff[c][i]->GetXaxis()->SetNdivisions(304,kTRUE);
            K_y_reldiff[c][i]->GetYaxis()->SetNdivisions(304,kTRUE);
            if(i==0)
            {
              K_y_reldiff[c][i]->GetYaxis()->SetTitle("#font[12]{M}^{#font[ 12]{h}}_{#font[ 12]{ratio}}");
              K_y_reldiff[c][i]->GetYaxis()->SetTitleSize(0.08);
            }
            K_y_reldiff[c][i]->Draw("SAMEP");
            K_y_reldiff[c][i]->GetXaxis()->SetLimits(0.1,0.9);
            K_y_reldiff[c][i]->SetMinimum(-1);
            K_y_reldiff[c][i]->SetMaximum(1);
            // l1.Draw("SAME");
            // l2.Draw("SAME");
            // l3.Draw("SAME");
            // l4.Draw("SAME");
            c731->Range(0.1,-1,0.9,1);
          }
          else
          {
            K_y_reldiff[c][i]->Draw("SAMEP");
            K_y_reldiff[c][i]->GetXaxis()->SetLimits(0.1,0.9);
            K_y_reldiff[c][i]->SetMinimum(-1);
            K_y_reldiff[c][i]->SetMaximum(1);
          }
        }
        c731->Update();
      }
      z_range_p_y_reldiff[c][i][j].clear();
      z_range_k_y_reldiff[c][i][j].clear();
      z_range_h_y_reldiff[c][i][j].clear();
    }
    sp_y.push_back(MultiplicitiesSum[0][0][0]+MultiplicitiesSum[0][1][0]);
    sk_y.push_back(MultiplicitiesSum[0][0][1]+MultiplicitiesSum[0][1][1]);
    sh_y.push_back(MultiplicitiesSum[0][0][3]+MultiplicitiesSum[0][1][3]);
    sp_y_err.push_back(sqrt(MultiplicitiesSum[1][0][0]+MultiplicitiesSum[1][1][0]));
    sk_y_err.push_back(sqrt(MultiplicitiesSum[1][0][1]+MultiplicitiesSum[1][1][1]));
    sh_y_err.push_back(sqrt(MultiplicitiesSum[1][0][3]+MultiplicitiesSum[1][1][3]));
    rp_y.push_back(MultiplicitiesSum[0][0][0] ? MultiplicitiesSum[0][1][0]/MultiplicitiesSum[0][0][0] : 0);
    rk_y.push_back(MultiplicitiesSum[0][0][1] ? MultiplicitiesSum[0][1][1]/MultiplicitiesSum[0][0][1] : 0);
    rh_y.push_back(MultiplicitiesSum[0][0][3] ? MultiplicitiesSum[0][1][3]/MultiplicitiesSum[0][0][3] : 0);
    rp_y_err.push_back(sqrt((MultiplicitiesSum[1][1][0]+pow(MultiplicitiesSum[0][1][0],2)*MultiplicitiesSum[1][0][0]/pow(MultiplicitiesSum[0][0][0],2))/pow(MultiplicitiesSum[0][0][0],2)));
    rk_y_err.push_back(sqrt((MultiplicitiesSum[1][1][1]+pow(MultiplicitiesSum[0][1][1],2)*MultiplicitiesSum[1][0][1]/pow(MultiplicitiesSum[0][0][1],2))/pow(MultiplicitiesSum[0][0][1],2)));
    rh_y_err.push_back(sqrt((MultiplicitiesSum[1][1][3]+pow(MultiplicitiesSum[0][1][3],2)*MultiplicitiesSum[1][0][3]/pow(MultiplicitiesSum[0][0][3],2))/pow(MultiplicitiesSum[0][0][3],2)));
  }

  for(int l=0; l<9; l++)
  {
    sx_range_p_y.push_back(x_range[l]);
    sx_range_k_y.push_back(x_range[l]);
    sx_range_h_y.push_back(x_range[l]);
    rx_range_p_y.push_back(x_range[l]);
    rx_range_k_y.push_back(x_range[l]);
    rx_range_h_y.push_back(x_range[l]);
  }

  for(int k=9; k>0; k--)
  {
    if(!sp_y[k-1]) {sp_y.erase(sp_y.begin()+k-1); sp_y_err.erase(sp_y_err.begin()+k-1); sx_range_p_y.erase(sx_range_p_y.begin()+k-1);}
    if(!sk_y[k-1]) {sk_y.erase(sk_y.begin()+k-1); sk_y_err.erase(sk_y_err.begin()+k-1); sx_range_k_y.erase(sx_range_k_y.begin()+k-1);}
    if(!sh_y[k-1]) {sh_y.erase(sh_y.begin()+k-1); sh_y_err.erase(sh_y_err.begin()+k-1); sx_range_h_y.erase(sx_range_h_y.begin()+k-1);}
    if(!rp_y[k-1]) {rp_y.erase(rp_y.begin()+k-1); rp_y_err.erase(rp_y_err.begin()+k-1); rx_range_p_y.erase(rx_range_p_y.begin()+k-1);}
    if(!rk_y[k-1]) {rk_y.erase(rk_y.begin()+k-1); rk_y_err.erase(rk_y_err.begin()+k-1); rx_range_k_y.erase(rx_range_k_y.begin()+k-1);}
    if(!rh_y[k-1]) {rh_y.erase(rh_y.begin()+k-1); rh_y_err.erase(rh_y_err.begin()+k-1); rx_range_h_y.erase(rx_range_h_y.begin()+k-1);}
  }

  sH_y = new TGraphErrors(Int_t(sh_y.size()),&(sx_range_h_y[0]),&(sh_y[0]),0,&(sh_y_err[0]));
  sP_y = new TGraphErrors(Int_t(sp_y.size()),&(sx_range_p_y[0]),&(sp_y[0]),0,&(sp_y_err[0]));
  sK_y = new TGraphErrors(Int_t(sk_y.size()),&(sx_range_k_y[0]),&(sk_y[0]),0,&(sk_y_err[0]));
  rH_y = new TGraphErrors(Int_t(rh_y.size()),&(rx_range_h_y[0]),&(rh_y[0]),0,&(rh_y_err[0]));
  rP_y = new TGraphErrors(Int_t(rp_y.size()),&(rx_range_p_y[0]),&(rp_y[0]),0,&(rp_y_err[0]));
  rK_y = new TGraphErrors(Int_t(rk_y.size()),&(rx_range_k_y[0]),&(rk_y[0]),0,&(rk_y_err[0]));

  sH_y->SetMarkerColor(fMarkerColor[0]);
  sP_y->SetMarkerColor(fMarkerColor[0]);
  sK_y->SetMarkerColor(fMarkerColor[0]);
  rH_y->SetMarkerColor(fMarkerColor[0]);
  rP_y->SetMarkerColor(fMarkerColor[0]);
  rK_y->SetMarkerColor(fMarkerColor[0]);

  sH_y->SetMarkerSize(3);
  sP_y->SetMarkerSize(3);
  sK_y->SetMarkerSize(3);
  rH_y->SetMarkerSize(3);
  rP_y->SetMarkerSize(3);
  rK_y->SetMarkerSize(3);

  sH_y->SetMarkerStyle(fMarkerStyle[0][1]);
  sP_y->SetMarkerStyle(fMarkerStyle[0][1]);
  sK_y->SetMarkerStyle(fMarkerStyle[0][1]);
  rH_y->SetMarkerStyle(fMarkerStyle[0][1]);
  rP_y->SetMarkerStyle(fMarkerStyle[0][1]);
  rK_y->SetMarkerStyle(fMarkerStyle[0][1]);

  sH_y->SetTitle("");
  sP_y->SetTitle("");
  sK_y->SetTitle("");
  rH_y->SetTitle("");
  rP_y->SetTitle("");
  rK_y->SetTitle("");

  sH_y->GetYaxis()->SetTitle("");
  sP_y->GetYaxis()->SetTitle("");
  sK_y->GetYaxis()->SetTitle("");
  rH_y->GetYaxis()->SetTitle("");
  rP_y->GetYaxis()->SetTitle("");
  rK_y->GetYaxis()->SetTitle("");

  sH_y->GetXaxis()->SetTitle("z");
  sP_y->GetXaxis()->SetTitle("z");
  sK_y->GetXaxis()->SetTitle("z");
  rH_y->GetXaxis()->SetTitle("z");
  rP_y->GetXaxis()->SetTitle("z");
  rK_y->GetXaxis()->SetTitle("z");

  c11->cd(1);
  gPad->SetFillStyle(4000);
  sH_y->Draw("PA");
  sH_y->GetXaxis()->SetLimits(0.01,1.);
  sH_y->SetMinimum(0.6);
  sH_y->SetMaximum(1.2);
  sH_y->SetTitle("");
  sH_y->GetXaxis()->SetTitle("#font[ 12]{x}");
  sH_y->GetXaxis()->SetNdivisions(304,kTRUE);
  sH_y->GetYaxis()->SetNdivisions(304,kTRUE);
  sH_y->GetYaxis()->SetTitle("#font[12]{M}^{#font[ 12]{h^{+}}}+#font[12]{M}^{#font[ 12]{h^{-}}}");
  c11->Range(0.1,0.6,0.9,1.2);
  gPad->SetLogx();
  c11->Update();

  c12->cd(1);
  gPad->SetFillStyle(4000);
  sP_y->Draw("PA");
  sP_y->GetXaxis()->SetLimits(0.01,1.);
  sP_y->SetMinimum(0.5);
  sP_y->SetMaximum(0.95);
  sP_y->SetTitle("");
  sP_y->GetXaxis()->SetTitle("#font[ 12]{x}");
  sP_y->GetXaxis()->SetNdivisions(304,kTRUE);
  sP_y->GetYaxis()->SetNdivisions(304,kTRUE);
  sP_y->GetYaxis()->SetTitle("#font[12]{M}^{#font[ 12]{#pi^{+}}}+#font[12]{M}^{#font[ 12]{#pi^{-}}}");
  c12->Range(0.1,0.5,0.9,0.95);
  gPad->SetLogx();
  c12->Update();

  c13->cd(1);
  gPad->SetFillStyle(4000);
  sK_y->Draw("PA");
  sK_y->GetXaxis()->SetLimits(0.01,1.);
  sK_y->SetMinimum(0.08);
  sK_y->SetMaximum(0.2);
  sK_y->SetTitle("");
  sK_y->GetXaxis()->SetTitle("#font[ 12]{x}");
  sK_y->GetXaxis()->SetNdivisions(304,kTRUE);
  sK_y->GetYaxis()->SetNdivisions(304,kTRUE);
  sK_y->GetYaxis()->SetTitle("#font[12]{M}^{#font[ 12]{K^{+}}}+#font[12]{M}^{#font[ 12]{K^{-}}}");
  c13->Range(0.1,0.08,0.9,0.2);
  gPad->SetLogx();
  c13->Update();

  c14->cd(1);
  gPad->SetFillStyle(4000);
  rH_y->Draw("PA");
  rH_y->GetXaxis()->SetLimits(0.01,1.);
  rH_y->SetMinimum(0.9);
  rH_y->SetMaximum(2.4);
  rH_y->SetTitle("");
  rH_y->GetXaxis()->SetTitle("#font[ 12]{x}");
  rH_y->GetXaxis()->SetNdivisions(304,kTRUE);
  rH_y->GetYaxis()->SetNdivisions(304,kTRUE);
  rH_y->GetYaxis()->SetTitle("#font[12]{M}^{#font[ 12]{h^{+}}}/#font[12]{M}^{#font[ 12]{h^{-}}}");
  c14->Range(0.1,0.9,0.9,2.4);
  gPad->SetLogx();
  c14->Update();

  c15->cd(1);
  gPad->SetFillStyle(4000);
  rP_y->Draw("PA");
  rP_y->GetXaxis()->SetLimits(0.01,1.);
  rP_y->SetMinimum(0.9);
  rP_y->SetMaximum(1.8);
  rP_y->SetTitle("");
  rP_y->GetXaxis()->SetTitle("#font[ 12]{x}");
  rP_y->GetXaxis()->SetNdivisions(304,kTRUE);
  rP_y->GetYaxis()->SetNdivisions(304,kTRUE);
  rP_y->GetYaxis()->SetTitle("#font[12]{M}^{#font[ 12]{#pi^{+}}}/#font[12]{M}^{#font[ 12]{#pi^{-}}}");
  c15->Range(0.1,0.9,0.9,1.8);
  gPad->SetLogx();
  c15->Update();

  c16->cd(1);
  gPad->SetFillStyle(4000);
  rK_y->Draw("PA");
  rK_y->GetXaxis()->SetLimits(0.01,1.);
  rK_y->SetMinimum(0.9);
  rK_y->SetMaximum(3.5);
  rK_y->SetTitle("");
  rK_y->GetXaxis()->SetTitle("#font[ 12]{x}");
  rK_y->GetXaxis()->SetNdivisions(304,kTRUE);
  rK_y->GetYaxis()->SetNdivisions(304,kTRUE);
  rK_y->GetYaxis()->SetTitle("#font[12]{M}^{#font[ 12]{K^{+}}}/#font[12]{M}^{#font[ 12]{K^{-}}}");
  c16->Range(0.1,0.9,0.9,3.5);
  gPad->SetLogx();
  c16->Update();

  TLatex fTitle;

  c51->cd(1);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 4.4,"0.004#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.01");
  c52->cd(1);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 4.4,"0.004#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.01");

  c51->cd(2);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 4.4,"0.01#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.02");
  c52->cd(2);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 4.4,"0.01#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.02");

  c51->cd(3);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 4.4,"0.02#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.03");
  c52->cd(3);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 4.4,"0.02#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.03");

  c51->cd(4);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 4.4,"0.03#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.04");
  c52->cd(4);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 4.4,"0.03#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.04");

  c51->cd(5);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 4.4,"0.04#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.06");
  c52->cd(5);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 4.4,"0.04#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.06");

  c51->cd(6);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 4.4,"0.06#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.1");
  c52->cd(6);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 4.4,"0.06#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.1");

  c51->cd(7);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 4.4,"0.1#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.14");
  c52->cd(7);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 4.4,"0.1#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.14");

  c51->cd(8);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 4.4,"0.14#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.18");
  c52->cd(8);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 4.4,"0.14#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.18");

  c51->cd(9);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 4.4,"0.18#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.4");
  c52->cd(9);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 4.4,"0.18#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.4");

  c51->cd(10);
  fTitle.SetTextSize(0.095);
  fTitle.SetTextAlign(11);
  fTitle.DrawLatex(0.05, 0.72,"#color[221]{0.70#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.90, #delta = 1.25}");
  fTitle.DrawLatex(0.05, 0.64,"#color[4]{0.50#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.70, #delta = 1.0}");
  fTitle.DrawLatex(0.05, 0.56,"#color[226]{0.30#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.50, #delta = 0.75}");
  fTitle.DrawLatex(0.05, 0.48,"#color[209]{0.20#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.30, #delta = 0.50}");
  fTitle.DrawLatex(0.05, 0.40,"#color[95]{0.15#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.20, #delta = 0.25}");
  fTitle.DrawLatex(0.05, 0.32,"#color[2]{0.10#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.15, #delta = 0}");
  c52->cd(10);
  fTitle.SetTextSize(0.095);
  fTitle.SetTextAlign(11);
  fTitle.DrawLatex(0.05, 0.72,"#color[221]{0.70#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.90, #delta = 1.25}");
  fTitle.DrawLatex(0.05, 0.64,"#color[4]{0.50#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.70, #delta = 1.0}");
  fTitle.DrawLatex(0.05, 0.56,"#color[226]{0.30#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.50, #delta = 0.75}");
  fTitle.DrawLatex(0.05, 0.48,"#color[209]{0.20#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.30, #delta = 0.50}");
  fTitle.DrawLatex(0.05, 0.40,"#color[95]{0.15#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.20, #delta = 0.25}");
  fTitle.DrawLatex(0.05, 0.32,"#color[2]{0.10#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.15, #delta = 0}");


  c61->cd(1);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 3.7,"0.004#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.01");
  c62->cd(1);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 3.7,"0.004#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.01");

  c61->cd(2);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 3.7,"0.01#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.02");
  c62->cd(2);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 3.7,"0.01#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.02");

  c61->cd(3);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 3.7,"0.02#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.03");
  c62->cd(3);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 3.7,"0.02#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.03");

  c61->cd(4);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 3.7,"0.03#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.04");
  c62->cd(4);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 3.7,"0.03#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.04");

  c61->cd(5);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 3.7,"0.04#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.06");
  c62->cd(5);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 3.7,"0.04#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.06");

  c61->cd(6);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 3.7,"0.06#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.1");
  c62->cd(6);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 3.7,"0.06#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.1");

  c61->cd(7);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 3.7,"0.1#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.14");
  c62->cd(7);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 3.7,"0.1#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.14");

  c61->cd(8);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 3.7,"0.14#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.18");
  c62->cd(8);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 3.7,"0.14#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.18");

  c61->cd(9);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 3.7,"0.18#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.4");
  c62->cd(9);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 3.7,"0.18#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.4");

  c61->cd(10);
  fTitle.SetTextSize(0.095);
  fTitle.SetTextAlign(11);
  fTitle.DrawLatex(0.05, 0.72,"#color[221]{0.70#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.90, #delta = 1.25}");
  fTitle.DrawLatex(0.05, 0.64,"#color[4]{0.50#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.70, #delta = 1.0}");
  fTitle.DrawLatex(0.05, 0.56,"#color[226]{0.30#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.50, #delta = 0.75}");
  fTitle.DrawLatex(0.05, 0.48,"#color[209]{0.20#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.30, #delta = 0.50}");
  fTitle.DrawLatex(0.05, 0.40,"#color[95]{0.15#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.20, #delta = 0.25}");
  fTitle.DrawLatex(0.05, 0.32,"#color[2]{0.10#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.15, #delta = 0}");
  c62->cd(10);
  fTitle.SetTextSize(0.095);
  fTitle.SetTextAlign(11);
  fTitle.DrawLatex(0.05, 0.72,"#color[221]{0.70#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.90, #delta = 1.25}");
  fTitle.DrawLatex(0.05, 0.64,"#color[4]{0.50#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.70, #delta = 1.0}");
  fTitle.DrawLatex(0.05, 0.56,"#color[226]{0.30#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.50, #delta = 0.75}");
  fTitle.DrawLatex(0.05, 0.48,"#color[209]{0.20#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.30, #delta = 0.50}");
  fTitle.DrawLatex(0.05, 0.40,"#color[95]{0.15#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.20, #delta = 0.25}");
  fTitle.DrawLatex(0.05, 0.32,"#color[2]{0.10#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.15, #delta = 0}");

  c71->cd(1);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 0.75,"0.004#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.01");
  c72->cd(1);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 0.75,"0.004#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.01");

  c71->cd(2);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 0.75,"0.01#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.02");
  c72->cd(2);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 0.75,"0.01#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.02");

  c71->cd(3);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 0.75,"0.02#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.03");
  c72->cd(3);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 0.75,"0.02#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.03");

  c71->cd(4);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 0.75,"0.03#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.04");
  c72->cd(4);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 0.75,"0.03#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.04");

  c71->cd(5);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 0.75,"0.04#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.06");
  c72->cd(5);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 0.75,"0.04#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.06");

  c71->cd(6);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 0.75,"0.06#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.1");
  c72->cd(6);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 0.75,"0.06#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.1");

  c71->cd(7);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 0.75,"0.1#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.14");
  c72->cd(7);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 0.75,"0.1#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.14");

  c71->cd(8);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 0.75,"0.14#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.18");
  c72->cd(8);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 0.75,"0.14#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.18");

  c71->cd(9);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 0.75,"0.18#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.4");
  c72->cd(9);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 0.75,"0.18#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.4");

  c71->cd(10);
  fTitle.SetTextSize(0.095);
  fTitle.SetTextAlign(11);
  fTitle.DrawLatex(0.05, 0.72,"#color[221]{0.70#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.90, #delta = 0.25}");
  fTitle.DrawLatex(0.05, 0.64,"#color[4]{0.50#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.70, #delta = 0.2}");
  fTitle.DrawLatex(0.05, 0.56,"#color[226]{0.30#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.50, #delta = 0.15}");
  fTitle.DrawLatex(0.05, 0.48,"#color[209]{0.20#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.30, #delta = 0.1}");
  fTitle.DrawLatex(0.05, 0.40,"#color[95]{0.15#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.20, #delta = 0.05}");
  fTitle.DrawLatex(0.05, 0.32,"#color[2]{0.10#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.15, #delta = 0}");
  c72->cd(10);
  fTitle.SetTextSize(0.095);
  fTitle.SetTextAlign(11);
  fTitle.DrawLatex(0.05, 0.72,"#color[221]{0.70#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.90, #delta = 0.25}");
  fTitle.DrawLatex(0.05, 0.64,"#color[4]{0.50#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.70, #delta = 0.2}");
  fTitle.DrawLatex(0.05, 0.56,"#color[226]{0.30#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.50, #delta = 0.15}");
  fTitle.DrawLatex(0.05, 0.48,"#color[209]{0.20#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.30, #delta = 0.1}");
  fTitle.DrawLatex(0.05, 0.40,"#color[95]{0.15#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.20, #delta = 0.05}");
  fTitle.DrawLatex(0.05, 0.32,"#color[2]{0.10#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.15, #delta = 0}");

  c8->cd(1);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 3.7,"0.004#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.01");

  c8->cd(2);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 3.7,"0.01#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.02");

  c8->cd(3);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 3.7,"0.02#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.03");

  c8->cd(4);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 3.7,"0.03#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.04");

  c8->cd(5);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 3.7,"0.04#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.06");

  c8->cd(6);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 3.7,"0.06#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.1");

  c8->cd(7);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 3.7,"0.1#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.14");

  c8->cd(8);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 3.7,"0.14#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.18");

  c8->cd(9);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 3.7,"0.18#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.4");

  c9->cd(1);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 3.2,"0.004#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.01");

  c9->cd(2);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 3.2,"0.01#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.02");

  c9->cd(3);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 3.2,"0.02#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.03");

  c9->cd(4);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 3.2,"0.03#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.04");

  c9->cd(5);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 3.2,"0.04#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.06");

  c9->cd(6);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 3.2,"0.06#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.1");

  c9->cd(7);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 3.2,"0.1#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.14");

  c9->cd(8);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 3.2,"0.14#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.18");

  c9->cd(9);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 3.2,"0.18#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.4");

  c10->cd(1);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 0.74,"0.004#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.01");

  c10->cd(2);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 0.74,"0.01#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.02");

  c10->cd(3);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 0.74,"0.02#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.03");

  c10->cd(4);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 0.74,"0.03#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.04");

  c10->cd(5);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 0.74,"0.04#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.06");

  c10->cd(6);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 0.74,"0.06#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.1");

  c10->cd(7);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 0.74,"0.1#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.14");

  c10->cd(8);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 0.74,"0.14#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.18");

  c10->cd(9);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 0.74,"0.18#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.4");


  c51->Update();
  c511->Update();
  c52->Update();
  c61->Update();
  c611->Update();
  c62->Update();
  c71->Update();
  c711->Update();
  c72->Update();
  c5->Update();
  c53->Update();
  c531->Update();
  c6->Update();
  c63->Update();
  c631->Update();
  c7->Update();
  c73->Update();
  c731->Update();
  c8->Update();
  c9->Update();
  c10->Update();
  c11->Update();
  c12->Update();
  c13->Update();
  c14->Update();
  c15->Update();
  c16->Update();

  c51->Print(Form("%s/hadron_multiplicity_file.pdf(",data_path),"pdf");
  c52->Print(Form("%s/hadron_multiplicity_file.pdf)",data_path),"pdf");
  c61->Print(Form("%s/pion_multiplicity_file.pdf(",data_path),"pdf");
  c62->Print(Form("%s/pion_multiplicity_file.pdf)",data_path),"pdf");
  c71->Print(Form("%s/kaon_multiplicity_file.pdf(",data_path),"pdf");
  c72->Print(Form("%s/kaon_multiplicity_file.pdf)",data_path),"pdf");
  c5->Print(Form("%s/hadron_multiplicity_zvtx_file.pdf(",data_path),"pdf");
  c53->Print(Form("%s/hadron_multiplicity_zvtx_file.pdf",data_path),"pdf");
  c511->Print(Form("%s/hadron_multiplicity_zvtx_file.pdf",data_path),"pdf");
  c531->Print(Form("%s/hadron_multiplicity_zvtx_file.pdf)",data_path),"pdf");
  c6->Print(Form("%s/pion_multiplicity_zvtx_file.pdf(",data_path),"pdf");
  c63->Print(Form("%s/pion_multiplicity_zvtx_file.pdf",data_path),"pdf");
  c611->Print(Form("%s/pion_multiplicity_zvtx_file.pdf",data_path),"pdf");
  c631->Print(Form("%s/pion_multiplicity_zvtx_file.pdf)",data_path),"pdf");
  c7->Print(Form("%s/kaon_multiplicity_zvtx_file.pdf(",data_path),"pdf");
  c73->Print(Form("%s/kaon_multiplicity_zvtx_file.pdf",data_path),"pdf");
  c711->Print(Form("%s/kaon_multiplicity_zvtx_file.pdf",data_path),"pdf");
  c731->Print(Form("%s/kaon_multiplicity_zvtx_file.pdf)",data_path),"pdf");
  c8->Print(Form("%s/hadron_multiplicity_yavg_file.pdf",data_path));
  c9->Print(Form("%s/pion_multiplicity_yavg_file.pdf",data_path));
  c10->Print(Form("%s/kaon_multiplicity_yavg_file.pdf",data_path));
  c11->Print(Form("%s/hadron_multiplicity_sum_file.pdf",data_path));
  c12->Print(Form("%s/pion_multiplicity_sum_file.pdf",data_path));
  c13->Print(Form("%s/kaon_multiplicity_sum_file.pdf",data_path));
  c14->Print(Form("%s/hadron_multiplicity_ratio_file.pdf",data_path));
  c15->Print(Form("%s/pion_multiplicity_ratio_file.pdf",data_path));
  c16->Print(Form("%s/kaon_multiplicity_ratio_file.pdf",data_path));

  ofs_p.close();
  ofs_k.close();
  ofs_h.close();
  ofs_m.close();
  ofs_mp.close();
  ofs_mm.close();

  return 0;
}

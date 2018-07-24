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
#include <TGaxis.h>

#include "analySIDIS_collect.h"

#define data_path "/sps/compass/npierre/Multiplicities"

// Flags
#define Y2006 0
#define Y2012 0
#define Y2016 1
#define DUMMY_ACC 0
#define RCUTSTUDY_ON 0

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

void yavg()
{
  for(int c=0; c<2; c++)
  {
    for(int x=0; x<9; x++)
    {
      for(int z=0; z<12; z++)
      {
        for(int i=0; i<6; i++)
        {
          for(auto period : fPeriods)
          {
            fBinning_yavg_period[period][x][z].tab[c][0][0]+=fBinning_period[period][x][i][z].tab[c][0][0];
            fBinning_yavg_period[period][x][z].tab[c][0][1]+=fBinning_period[period][x][i][z].tab[c][0][1];
            fBinning_yavg_period[period][x][z].tab[c][0][2]+=fBinning_period[period][x][i][z].tab[c][0][2];
            fBinning_yavg_period[period][x][z].tab[c][0][3]+=fBinning_period[period][x][i][z].tab[c][0][3];
            fBinning_yavg_period[period][x][z].tab[c][1][0]+=fBinning_period[period][x][i][z].tab[c][1][0];
            fBinning_yavg_period[period][x][z].tab[c][1][1]+=fBinning_period[period][x][i][z].tab[c][1][1];
            fBinning_yavg_period[period][x][z].tab[c][1][2]+=fBinning_period[period][x][i][z].tab[c][1][2];
            fBinning_yavg_period[period][x][z].tab[c][1][3]+=fBinning_period[period][x][i][z].tab[c][1][3];
          }
          fBinning_yavg[x][z].tab[c][0][0]+=fBinning[x][i][z].tab[c][0][0];
          fBinning_yavg[x][z].tab[c][0][1]+=fBinning[x][i][z].tab[c][0][1];
          fBinning_yavg[x][z].tab[c][0][2]+=fBinning[x][i][z].tab[c][0][2];
          fBinning_yavg[x][z].tab[c][0][3]+=fBinning[x][i][z].tab[c][0][3];
          fBinning_yavg[x][z].tab[c][1][0]+=fBinning[x][i][z].tab[c][1][0];
          fBinning_yavg[x][z].tab[c][1][1]+=fBinning[x][i][z].tab[c][1][1];
          fBinning_yavg[x][z].tab[c][1][2]+=fBinning[x][i][z].tab[c][1][2];
          fBinning_yavg[x][z].tab[c][1][3]+=fBinning[x][i][z].tab[c][1][3];
          if(!c && fBinning[x][i][z].tab[c][0][3])
          {
            fNDIS_evt_yavg[0][x][z]+=fNDIS_evt[0][x][i][z];
            fNDIS_evt_yavg[1][x][z]+=fNDIS_evt_err[0][x][i][z];
          }
          fRich_sys_err_yavg[x][z].tab[c][1][0]+=pow(fRich_sys_err[x][i][z].tab[c][1][0],2);
          fRich_sys_err_yavg[x][z].tab[c][1][1]+=pow(fRich_sys_err[x][i][z].tab[c][1][1],2);
          fRich_sys_err_yavg[x][z].tab[c][1][2]+=pow(fRich_sys_err[x][i][z].tab[c][1][2],2);
          fRich_sys_err_yavg[x][z].tab[c][1][3]+=pow(fRich_sys_err[x][i][z].tab[c][1][3],2);
        }
        fRich_sys_err_yavg[x][z].tab[c][1][0]=sqrt(fRich_sys_err_yavg[x][z].tab[c][1][0]);
        fRich_sys_err_yavg[x][z].tab[c][1][1]=sqrt(fRich_sys_err_yavg[x][z].tab[c][1][1]);
        fRich_sys_err_yavg[x][z].tab[c][1][2]=sqrt(fRich_sys_err_yavg[x][z].tab[c][1][2]);
        fRich_sys_err_yavg[x][z].tab[c][1][3]=sqrt(fRich_sys_err_yavg[x][z].tab[c][1][3]);
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
                fAcceptance_weighted[i][j][k].tab[c][0][l] += fBinning_period[period][i][j][k].tab[c][0][l]*fAcceptance[period][i][j][k].tab[c][0][l]/fBinning[i][j][k].tab[c][0][l];
                fAcceptance_weighted[i][j][k].tab[c][1][l] += fBinning_period[period][i][j][k].tab[c][0][l]*fAcceptance[period][i][j][k].tab[c][1][l]/fBinning[i][j][k].tab[c][0][l];
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

  ifstream periods(argv[1]);
  string filelist, periodName;
  int periodBit;
  while(periods >> periodName)
  {
    fNumberPeriod++;
    periods >> periodBit;

    if(!periodBit) continue;

    cout << periodName << " ";

    if(!DUMMY_ACC)
    {
      fetch_acceptance(Form("acceptance/%d/acceptance_%s.txt",year,periodName.c_str()),fNumberPeriod-1);
      fetch_yavg_acceptance(Form("acceptance/%d/acceptance_yavg_%s.txt",year,periodName.c_str()),fNumberPeriod-1);
    }
    else
    {
      dummy_acceptance();
    }

    ifstream dis_file(Form("rawmult/%d/DIS_%s.txt",year,periodName.c_str()));
    ifstream had_file(Form("rawmult/%d/hadron_%s.txt",year,periodName.c_str()));

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

  yavg();
  weight_acceptance();
  weight_meanvalues();

  TCanvas* c5;
  c5 = new TCanvas("Hadron_Multiplicities","Hadron_Multiplicities",3200,1600);

  TCanvas* c6;
  c6 = new TCanvas("Pion_Multiplicities","Pion_Multiplicities",3200,1600);

  TCanvas* c7;
  c7 = new TCanvas("Kaon_Multiplicities","Kaon_Multiplicities",3200,1600);

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

  c5->SetFillColor(0);
  //c5->SetFrameFillStyle(4000);
  c6->SetFillColor(0);
  //c6->SetFrameFillStyle(4000);
  c7->SetFillColor(0);
  //c7->SetFrameFillStyle(4000);
  c8->SetFillColor(0);
  //c8->SetFrameFillStyle(4000);
  c9->SetFillColor(0);
  //c9->SetFrameFillStyle(4000);
  c10->SetFillColor(0);
  //c10->SetFrameFillStyle(4000);
  c11->SetFillColor(0);
  //c11->SetFrameFillStyle(4000);
  c12->SetFillColor(0);
  //c12->SetFrameFillStyle(4000);
  c13->SetFillColor(0);
  //c13->SetFrameFillStyle(4000);
  c14->SetFillColor(0);
  //c14->SetFrameFillStyle(4000);
  c15->SetFillColor(0);
  //c15->SetFrameFillStyle(4000);
  c16->SetFillColor(0);
  //c16->SetFrameFillStyle(4000);

  c5->Divide(5,2,0,0);
  c6->Divide(5,2,0,0);
  c7->Divide(5,2,0,0);
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
  TGraphErrors* H_y[2][9];
  TGraphErrors* P_y[2][9];
  TGraphErrors* K_y[2][9];
  TGraphErrors* sH_y;
  TGraphErrors* sP_y;
  TGraphErrors* sK_y;
  TGraphErrors* rH_y;
  TGraphErrors* rP_y;
  TGraphErrors* rK_y;

  Double_t z_range[12] = {.225,.275,.325,.375,.425,.475,.525,.575,.625,.675,.725,.8};
  Double_t x_range[9] = {.008,.015,.025,.035,.05,.08,.12,.16,.29};

  ofstream ofs_p(Form("%s/multiplicities_pion.txt",data_path), std::ofstream::out | std::ofstream::trunc);
  ofstream ofs_t(Form("%s/multiplicities_raw.txt",data_path), std::ofstream::out | std::ofstream::trunc);
  ofstream ofs_k(Form("%s/multiplicities_kaon.txt",data_path), std::ofstream::out | std::ofstream::trunc);
  ofstream ofs_h(Form("%s/multiplicities_hadron.txt",data_path), std::ofstream::out | std::ofstream::trunc);

  std::vector<Double_t> p_m[2][9][6];
  std::vector<Double_t> k_m[2][9][6];
  std::vector<Double_t> h_m[2][9][6];
  std::vector<Double_t> p_err[2][9][6];
  std::vector<Double_t> k_err[2][9][6];
  std::vector<Double_t> h_err[2][9][6];
  std::vector<Double_t> z_range_p[2][9][6];
  std::vector<Double_t> z_range_k[2][9][6];
  std::vector<Double_t> z_range_h[2][9][6];
  std::vector<Double_t> p_y[2][9];
  std::vector<Double_t> k_y[2][9];
  std::vector<Double_t> h_y[2][9];
  std::vector<Double_t> p_y_err[2][9];
  std::vector<Double_t> k_y_err[2][9];
  std::vector<Double_t> h_y_err[2][9];
  std::vector<Double_t> z_range_p_y[2][9];
  std::vector<Double_t> z_range_k_y[2][9];
  std::vector<Double_t> z_range_h_y[2][9];
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

  for(int i=0; i<9; i++)
  {
    int axisflag = 0;
    for(int j=0; j<6; j++)
    {
      for(int k=0; k<12; k++)
      {
        for(int c=1; c>=0; c--)
        {
          for(int l=0; l<4; l++)
          {
            fMultiplicities[i][j][k].tab[c][0][l] = (fNDIS_evt[0][i][j][k] && fAcceptance_weighted[i][j][k].tab[c][0][l] ?
                                                    Double_t(fBinning[i][j][k].tab[c][0][l]/(fNDIS_evt[0][i][j][k]*fZ_bin_width[k]*fAcceptance_weighted[i][j][k].tab[c][0][l]))
                                                    : 0);
            fMultiplicities[i][j][k].tab[c][1][l] = (fNDIS_evt[0][i][j][k] && fAcceptance_weighted[i][j][k].tab[c][0][l] ?
                                                    Double_t(((fBinning[i][j][k].tab[c][1][l]/pow(fNDIS_evt[0][i][j][k],2)-pow(fBinning[i][j][k].tab[c][0][l],2)*
                                                    fNDIS_evt_err[0][i][j][k]/pow(fNDIS_evt[0][i][j][k],4))/(pow(fZ_bin_width[k]*fAcceptance_weighted[i][j][k].tab[c][0][l],2)))
                                                    + fAcceptance_weighted[i][j][k].tab[c][1][l]*pow(fBinning[i][j][k].tab[c][0][l]/(fNDIS_evt[0][i][j][k]*fZ_bin_width[k]*pow(fAcceptance_weighted[i][j][k].tab[c][0][l],2)),2))
                                                    : 0);
            fMultiplicities[i][j][k].tab[c][2][l] = (fNDIS_evt[0][i][j][k] ?
                                                    Double_t(sqrt(pow(fRich_sys_err[i][j][k].tab[c][1][l],2)/pow(fNDIS_evt[0][i][j][k]*fZ_bin_width[k]*fAcceptance_weighted[i][j][k].tab[c][0][l],2)+
                                                    pow(0.05*sqrt(fAcceptance_weighted[i][j][k].tab[c][1][l])*fBinning[i][j][k].tab[c][0][l]/(fNDIS_evt[0][i][j][k]*fZ_bin_width[k]
                                                    *pow(fAcceptance_weighted[i][j][k].tab[c][0][l],2)),2)))
                                                    : 0);

            if(fMultiplicities[i][j][k].tab[c][0][l]<0)
            {
              fMultiplicities[i][j][k].tab[c][0][l] = 0 ;
              fMultiplicities[i][j][k].tab[c][1][l] = 0 ;
              fMultiplicities[i][j][k].tab[c][2][l] = 0 ;
            }
          }


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

          // cout << c << " " << i << " " << j << " " << k << " " << fMultiplicities[i][j][k].tab[c][0][3] << endl;

          p_m[c][i][j].push_back(fMultiplicities[i][j][k].tab[c][0][0]>0 ? fMultiplicities[i][j][k].tab[c][0][0]+0.1*j : 0);
          k_m[c][i][j].push_back(fMultiplicities[i][j][k].tab[c][0][1]>0 ? fMultiplicities[i][j][k].tab[c][0][1]+0.1*j : 0);
          h_m[c][i][j].push_back(fMultiplicities[i][j][k].tab[c][0][3]>0 ? fMultiplicities[i][j][k].tab[c][0][3]+0.1*j : 0);
          p_err[c][i][j].push_back(fMultiplicities[i][j][k].tab[c][1][0] ? sqrt(fMultiplicities[i][j][k].tab[c][1][0]) : 0);
          k_err[c][i][j].push_back(fMultiplicities[i][j][k].tab[c][1][1] ? sqrt(fMultiplicities[i][j][k].tab[c][1][1]) : 0);
          h_err[c][i][j].push_back(fMultiplicities[i][j][k].tab[c][1][3] ? sqrt(fMultiplicities[i][j][k].tab[c][1][3]) : 0);
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
          if(!p_m[c][i][j][k-1]) {p_m[c][i][j].erase(p_m[c][i][j].begin()+k-1); p_err[c][i][j].erase(p_err[c][i][j].begin()+k-1); z_range_p[c][i][j].erase(z_range_p[c][i][j].begin()+k-1);}
          if(!k_m[c][i][j][k-1]) {k_m[c][i][j].erase(k_m[c][i][j].begin()+k-1); k_err[c][i][j].erase(k_err[c][i][j].begin()+k-1); z_range_k[c][i][j].erase(z_range_k[c][i][j].begin()+k-1);}
          if(!h_m[c][i][j][k-1]) {h_m[c][i][j].erase(h_m[c][i][j].begin()+k-1); h_err[c][i][j].erase(h_err[c][i][j].begin()+k-1); z_range_h[c][i][j].erase(z_range_h[c][i][j].begin()+k-1);}
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

        H_mult[c][i][j]->SetMarkerColor(fMarkerColor[j]);
        P_mult[c][i][j]->SetMarkerColor(fMarkerColor[j]);
        K_mult[c][i][j]->SetMarkerColor(fMarkerColor[j]);

        H_mult[c][i][j]->SetMarkerSize(3);
        P_mult[c][i][j]->SetMarkerSize(3);
        K_mult[c][i][j]->SetMarkerSize(3);

        H_mult[c][i][j]->SetMarkerStyle(fMarkerStyle[j][c]);
        P_mult[c][i][j]->SetMarkerStyle(fMarkerStyle[j][c]);
        K_mult[c][i][j]->SetMarkerStyle(fMarkerStyle[j][c]);

        H_mult[c][i][j]->SetTitle("");
        P_mult[c][i][j]->SetTitle("");
        K_mult[c][i][j]->SetTitle("");

        H_mult[c][i][j]->GetYaxis()->SetTitle("");
        P_mult[c][i][j]->GetYaxis()->SetTitle("");
        K_mult[c][i][j]->GetYaxis()->SetTitle("");

        H_mult[c][i][j]->GetXaxis()->SetTitle("z");
        P_mult[c][i][j]->GetXaxis()->SetTitle("z");
        K_mult[c][i][j]->GetXaxis()->SetTitle("z");


        if(!h_m_empty)
        {
          c5->cd(i+1);
          gPad->SetFillStyle(4000);
          if(H_mult[c][i][j])
          {
            if(!c && !axisflag)
            {
              H_mult[c][i][j]->Draw("SAMEPA");
              H_mult[c][i][j]->GetXaxis()->SetLimits(0.1,0.9);
              H_mult[c][i][j]->SetMinimum(0.);
              H_mult[c][i][j]->SetMaximum(4.);
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
                H_mult[c][i][j]->GetYaxis()->SetTitle("#font[12]{M}^{#font[ 12]{h}}+ #font[ 12]{#delta}");
                H_mult[c][i][j]->GetYaxis()->SetTitleSize(0.08);
              }
              // H_mult[c][i][0]->Draw("SAMEP");
              // H_mult[c][i][0]->GetXaxis()->SetLimits(0.1,0.9);
              // H_mult[c][i][0]->SetMinimum(0.);
              // H_mult[c][i][0]->SetMaximum(4.);
              // H_mult[c][i][1]->Draw("SAMEP");
              // H_mult[c][i][1]->GetXaxis()->SetLimits(0.1,0.9);
              // H_mult[c][i][1]->SetMinimum(0.);
              // H_mult[c][i][1]->SetMaximum(4.);
              // H_mult[c][i][2]->Draw("SAMEP");
              // H_mult[c][i][2]->GetXaxis()->SetLimits(0.1,0.9);
              // H_mult[c][i][2]->SetMinimum(0.);
              // H_mult[c][i][2]->SetMaximum(4.);
              axisflag=1;
              c5->Range(0.1,0.,0.9,4.);
            }
            else
            {
              H_mult[c][i][j]->Draw("SAMEP");
              H_mult[c][i][j]->GetXaxis()->SetLimits(0.1,0.9);
              H_mult[c][i][j]->SetMinimum(0.);
              H_mult[c][i][j]->SetMaximum(4.);
            }
          }
          c5->Update();
        }

        if(!p_m_empty)
        {
          c6->cd(i+1);
          if(P_mult[c][i][j])
          {
            if(!c && j==3)
            {
              P_mult[c][i][j]->Draw("SAMEPA");
              P_mult[c][i][j]->GetXaxis()->SetLimits(0.1,0.9);
              P_mult[c][i][j]->SetMinimum(0.);
              P_mult[c][i][j]->SetMaximum(2.0);
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
                P_mult[c][i][j]->GetYaxis()->SetTitle("#font[12]{M}^{#font[ 12]{#pi}}+ #font[ 12]{#delta}");
                P_mult[c][i][j]->GetYaxis()->SetTitleSize(0.08);
              }
              P_mult[c][i][0]->Draw("SAMEP");
              P_mult[c][i][0]->GetXaxis()->SetLimits(0.1,0.9);
              P_mult[c][i][0]->SetMinimum(0.);
              P_mult[c][i][0]->SetMaximum(2.);
              P_mult[c][i][1]->Draw("SAMEP");
              P_mult[c][i][1]->GetXaxis()->SetLimits(0.1,0.9);
              P_mult[c][i][1]->SetMinimum(0.);
              P_mult[c][i][1]->SetMaximum(2.);
              P_mult[c][i][2]->Draw("SAMEP");
              P_mult[c][i][2]->GetXaxis()->SetLimits(0.1,0.9);
              P_mult[c][i][2]->SetMinimum(0.);
              P_mult[c][i][2]->SetMaximum(2.);
              c6->Range(0.1,0.,0.9,2.);
            }
            else
            {
              P_mult[c][i][j]->Draw("SAMEP");
              P_mult[c][i][j]->GetXaxis()->SetLimits(0.1,0.9);
              P_mult[c][i][j]->SetMinimum(0.);
              P_mult[c][i][j]->SetMaximum(2.0);
            }
          }
          c6->Update();
        }

        if(!k_m_empty)
        {
          c7->cd(i+1);
          if(K_mult[c][i][j])
          {
            if(!c && j==3)
            {
              K_mult[c][i][j]->Draw("SAMEPA");
              K_mult[c][i][j]->GetXaxis()->SetLimits(0.1,0.9);
              K_mult[c][i][j]->SetMinimum(0.);
              K_mult[c][i][j]->SetMaximum(1.);
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
                K_mult[c][i][j]->GetYaxis()->SetTitle("#font[12]{M}^{#font[ 12]{K}}+ #font[ 12]{#delta}");
                K_mult[c][i][j]->GetYaxis()->SetTitleSize(0.08);
              }
              K_mult[c][i][0]->Draw("SAMEP");
              K_mult[c][i][0]->GetXaxis()->SetLimits(0.1,0.9);
              K_mult[c][i][0]->SetMinimum(0.);
              K_mult[c][i][0]->SetMaximum(1.);
              K_mult[c][i][1]->Draw("SAMEP");
              K_mult[c][i][1]->GetXaxis()->SetLimits(0.1,0.9);
              K_mult[c][i][1]->SetMinimum(0.);
              K_mult[c][i][1]->SetMaximum(1.);
              K_mult[c][i][2]->Draw("SAMEP");
              K_mult[c][i][2]->GetXaxis()->SetLimits(0.1,0.9);
              K_mult[c][i][2]->SetMinimum(0.);
              K_mult[c][i][2]->SetMaximum(1.);
              c7->Range(0.1,.0,0.9,1.);
            }
            else
            {
              K_mult[c][i][j]->Draw("SAMEP");
              K_mult[c][i][j]->GetXaxis()->SetLimits(0.1,0.9);
              K_mult[c][i][j]->SetMinimum(0.);
              K_mult[c][i][j]->SetMaximum(1.);
            }
          }
          c7->Update();
        }
        z_range_p[c][i][j].clear();
        z_range_k[c][i][j].clear();
        z_range_h[c][i][j].clear();
      }
    }

    Double_t MultiplicitiesSum[2][2][3];
    Double_t MultiplicitiesRatio[2][2][3];

    for(int c=0; c<2; c++)
    {
      for(int k=0; k<12; k++)
      {
        for(int l=0; l<4; l++)
        {
          cout << l << " " << fXrange[i] << " " << fZrange[k] << " " << fBinning_yavg[i][k].tab[c][0][l] << " " << fNDIS_evt_yavg[0][i][k] << " " << fAcceptance_yavg_weighted[i][k].tab[c][0][l] << endl;
          fMultiplicities_yavg[i][k].tab[c][0][l] = (fNDIS_evt_yavg[0][i][k] && fAcceptance_yavg_weighted[i][k].tab[c][0][l] ?
                                                  Double_t(fBinning_yavg[i][k].tab[c][0][l]/(fNDIS_evt_yavg[0][i][k]*fZ_bin_width[k]*fAcceptance_yavg_weighted[i][k].tab[c][0][l]))
                                                  : 0);
          fMultiplicities_yavg[i][k].tab[c][1][l] = (fNDIS_evt_yavg[0][i][k] && fAcceptance_yavg_weighted[i][k].tab[c][0][l] ?
                                                  Double_t(((fBinning_yavg[i][k].tab[c][1][l]/pow(fNDIS_evt_yavg[0][i][k],2)-pow(fBinning_yavg[i][k].tab[c][0][l],2)*
                                                  fNDIS_evt_yavg[1][i][k]/pow(fNDIS_evt_yavg[0][i][k],4))/(pow(fZ_bin_width[k]*fAcceptance_yavg_weighted[i][k].tab[c][0][l],2)))
                                                  + fAcceptance_yavg_weighted[i][k].tab[c][1][l]*pow(fBinning_yavg[i][k].tab[c][0][l]/(fNDIS_evt_yavg[0][i][k]*fZ_bin_width[k]*pow(fAcceptance_yavg_weighted[i][k].tab[c][0][l],2)),2))
                                                  : 0);
          fMultiplicities_yavg[i][k].tab[c][2][l] = (fNDIS_evt_yavg[0][i][k] ?
                                                  Double_t(sqrt(pow(fRich_sys_err_yavg[i][k].tab[c][1][l],2)/pow(fNDIS_evt_yavg[0][i][k]*fZ_bin_width[k]*fAcceptance_yavg_weighted[i][k].tab[c][0][l],2)+
                                                  pow(0.05*sqrt(fAcceptance_yavg_weighted[i][k].tab[c][1][l])*fBinning_yavg[i][k].tab[c][0][l]/(fNDIS_evt_yavg[0][i][k]*fZ_bin_width[k]
                                                  *pow(fAcceptance_yavg_weighted[i][k].tab[c][0][l],2)),2)))
                                                  : 0);
          if(fMultiplicities_yavg[i][k].tab[c][0][l]<0)
          {
            fMultiplicities_yavg[i][k].tab[c][0][l] = 0 ;
            fMultiplicities_yavg[i][k].tab[c][1][l] = 0 ;
            fMultiplicities_yavg[i][k].tab[c][2][l] = 0 ;
          }
        }

        // cout << c << " " << i << " " << k << " " << fMultiplicities_yavg[i][k].tab[c][0][3] << " " << fMultiplicities_yavg[i][k].tab[c][1][3] << " " << fMultiplicities_yavg[i][k].tab[c][2][3] << endl;

        p_y[c][i].push_back(fMultiplicities_yavg[i][k].tab[c][0][0]);
        k_y[c][i].push_back(fMultiplicities_yavg[i][k].tab[c][0][1]);
        h_y[c][i].push_back(fMultiplicities_yavg[i][k].tab[c][0][3]);
        p_y_err[c][i].push_back(sqrt(fMultiplicities_yavg[i][k].tab[c][1][0]));
        k_y_err[c][i].push_back(sqrt(fMultiplicities_yavg[i][k].tab[c][1][1]));
        h_y_err[c][i].push_back(sqrt(fMultiplicities_yavg[i][k].tab[c][1][3]));

        MultiplicitiesSum[0][c][0] += fMultiplicities_yavg[i][k].tab[c][0][0]*fZ_bin_width[k];
        MultiplicitiesSum[0][c][1] += fMultiplicities_yavg[i][k].tab[c][0][1]*fZ_bin_width[k];
        MultiplicitiesSum[0][c][3] += fMultiplicities_yavg[i][k].tab[c][0][3]*fZ_bin_width[k];
        MultiplicitiesRatio[0][c][0] += fMultiplicities_yavg[i][k].tab[c][0][0]*fZ_bin_width[k];
        MultiplicitiesRatio[0][c][1] += fMultiplicities_yavg[i][k].tab[c][0][1]*fZ_bin_width[k];
        MultiplicitiesRatio[0][c][3] += fMultiplicities_yavg[i][k].tab[c][0][3]*fZ_bin_width[k];
        MultiplicitiesSum[1][c][0] += fMultiplicities_yavg[i][k].tab[c][1][0]*fZ_bin_width[k];
        MultiplicitiesSum[1][c][1] += fMultiplicities_yavg[i][k].tab[c][1][1]*fZ_bin_width[k];
        MultiplicitiesSum[1][c][3] += fMultiplicities_yavg[i][k].tab[c][1][3]*fZ_bin_width[k];
        MultiplicitiesRatio[1][c][0] += fMultiplicities_yavg[i][k].tab[c][1][0]*fZ_bin_width[k];
        MultiplicitiesRatio[1][c][1] += fMultiplicities_yavg[i][k].tab[c][1][1]*fZ_bin_width[k];
        MultiplicitiesRatio[1][c][3] += fMultiplicities_yavg[i][k].tab[c][1][3]*fZ_bin_width[k];
      }

      for(int l=0; l<12; l++)
      {
        z_range_p_y[c][i].push_back(z_range[l]);
        z_range_k_y[c][i].push_back(z_range[l]);
        z_range_h_y[c][i].push_back(z_range[l]);
      }

      for(int k=12; k>0; k--)
      {
        if(!p_y[c][i][k-1]) {p_y[c][i].erase(p_y[c][i].begin()+k-1); p_y_err[c][i].erase(p_y_err[c][i].begin()+k-1); z_range_p_y[c][i].erase(z_range_p_y[c][i].begin()+k-1);}
        if(!k_y[c][i][k-1]) {k_y[c][i].erase(k_y[c][i].begin()+k-1); k_y_err[c][i].erase(k_y_err[c][i].begin()+k-1); z_range_k_y[c][i].erase(z_range_k_y[c][i].begin()+k-1);}
        if(!h_y[c][i][k-1]) {h_y[c][i].erase(h_y[c][i].begin()+k-1); h_y_err[c][i].erase(h_y_err[c][i].begin()+k-1); z_range_h_y[c][i].erase(z_range_h_y[c][i].begin()+k-1);}
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
        H_y[c][i]->SetMarkerColor(fMarkerColor[4]);
        P_y[c][i]->SetMarkerColor(fMarkerColor[4]);
        K_y[c][i]->SetMarkerColor(fMarkerColor[4]);
      }
      else
      {
        H_y[c][i]->SetMarkerColor(fMarkerColor[0]);
        P_y[c][i]->SetMarkerColor(fMarkerColor[0]);
        K_y[c][i]->SetMarkerColor(fMarkerColor[0]);
      }

      H_y[c][i]->SetMarkerSize(3);
      P_y[c][i]->SetMarkerSize(3);
      K_y[c][i]->SetMarkerSize(3);

      H_y[c][i]->SetMarkerStyle(fMarkerStyle[0][c]);
      P_y[c][i]->SetMarkerStyle(fMarkerStyle[0][c]);
      K_y[c][i]->SetMarkerStyle(fMarkerStyle[0][c]);

      H_y[c][i]->SetTitle("");
      P_y[c][i]->SetTitle("");
      K_y[c][i]->SetTitle("");

      H_y[c][i]->GetYaxis()->SetTitle("");
      P_y[c][i]->GetYaxis()->SetTitle("");
      K_y[c][i]->GetYaxis()->SetTitle("");

      H_y[c][i]->GetXaxis()->SetTitle("z");
      P_y[c][i]->GetXaxis()->SetTitle("z");
      K_y[c][i]->GetXaxis()->SetTitle("z");

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
            H_y[c][i]->SetMinimum(0.);
            H_y[c][i]->SetMaximum(5.);
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
            H_y[c][i]->Draw("SAMEP");
            H_y[c][i]->GetXaxis()->SetLimits(0.1,0.9);
            H_y[c][i]->SetMinimum(0.);
            H_y[c][i]->SetMaximum(5.);
            c8->Range(0.1,0.,0.9,5.);
          }
          else
          {
            H_y[c][i]->Draw("SAMEP");
            H_y[c][i]->GetXaxis()->SetLimits(0.1,0.9);
            H_y[c][i]->SetMinimum(0.);
            H_y[c][i]->SetMaximum(5.);
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
            P_y[c][i]->SetMinimum(0.);
            P_y[c][i]->SetMaximum(2.0);
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
            P_y[c][i]->Draw("SAMEP");
            P_y[c][i]->GetXaxis()->SetLimits(0.1,0.9);
            P_y[c][i]->SetMinimum(0.);
            P_y[c][i]->SetMaximum(2.);
            c9->Range(0.1,0.,0.9,2.);
          }
          else
          {
            P_y[c][i]->Draw("SAMEP");
            P_y[c][i]->GetXaxis()->SetLimits(0.1,0.9);
            P_y[c][i]->SetMinimum(0.);
            P_y[c][i]->SetMaximum(2.0);
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
            K_y[c][i]->SetMinimum(0.);
            K_y[c][i]->SetMaximum(1.);
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
            K_y[c][i]->Draw("SAMEP");
            K_y[c][i]->GetXaxis()->SetLimits(0.1,0.9);
            K_y[c][i]->SetMinimum(0.);
            K_y[c][i]->SetMaximum(1.);
            c10->Range(0.1,.0,0.9,1.);
          }
          else
          {
            K_y[c][i]->Draw("SAMEP");
            K_y[c][i]->GetXaxis()->SetLimits(0.1,0.9);
            K_y[c][i]->SetMinimum(0.);
            K_y[c][i]->SetMaximum(1.);
          }
        }
        c10->Update();
      }
      z_range_p_y[c][i].clear();
      z_range_k_y[c][i].clear();
      z_range_h_y[c][i].clear();
    }
    sp_y.push_back(MultiplicitiesSum[0][0][0]+MultiplicitiesSum[0][1][0]);
    sk_y.push_back(MultiplicitiesSum[0][0][1]+MultiplicitiesSum[0][1][1]);
    sh_y.push_back(MultiplicitiesSum[0][0][3]+MultiplicitiesSum[0][1][3]);
    sp_y_err.push_back(sqrt(MultiplicitiesSum[1][0][0]+MultiplicitiesSum[1][1][0]));
    sk_y_err.push_back(sqrt(MultiplicitiesSum[1][0][1]+MultiplicitiesSum[1][1][1]));
    sh_y_err.push_back(sqrt(MultiplicitiesSum[1][0][3]+MultiplicitiesSum[1][1][3]));
    rp_y.push_back(MultiplicitiesSum[0][1][0] ? MultiplicitiesSum[0][0][0]/MultiplicitiesSum[0][1][0] : 0);
    rk_y.push_back(MultiplicitiesSum[0][1][1] ? MultiplicitiesSum[0][0][1]/MultiplicitiesSum[0][1][1] : 0);
    rh_y.push_back(MultiplicitiesSum[0][1][3] ? MultiplicitiesSum[0][0][3]/MultiplicitiesSum[0][1][3] : 0);
    rp_y_err.push_back(sqrt(MultiplicitiesSum[1][0][0]+MultiplicitiesSum[1][1][0]));
    rk_y_err.push_back(sqrt(MultiplicitiesSum[1][0][1]+MultiplicitiesSum[1][1][1]));
    rh_y_err.push_back(sqrt(MultiplicitiesSum[1][0][3]+MultiplicitiesSum[1][1][3]));
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

  bool sp_y_empty = 0;
  bool rp_y_empty = 0;
  bool sk_y_empty = 0;
  bool rk_y_empty = 0;
  bool sh_y_empty = 0;
  bool rh_y_empty = 0;

  if(!(Int_t(sp_y[c][i].size()))) sp_y_empty = 1;
  if(!(Int_t(sk_y[c][i].size()))) sk_y_empty = 1;
  if(!(Int_t(sh_y[c][i].size()))) sh_y_empty = 1;
  if(!(Int_t(rp_y[c][i].size()))) rp_y_empty = 1;
  if(!(Int_t(rk_y[c][i].size()))) rk_y_empty = 1;
  if(!(Int_t(rh_y[c][i].size()))) rh_y_empty = 1;

  sH_y[c][i] = new TGraphErrors(Int_t(sh_y[c][i].size()),&(sx_range_h_y[c][i][0]),&(sh_y[c][i][0]),0,&(sh_y_err[c][i][0]));
  sP_y[c][i] = new TGraphErrors(Int_t(sp_y[c][i].size()),&(sx_range_p_y[c][i][0]),&(sp_y[c][i][0]),0,&(sp_y_err[c][i][0]));
  sK_y[c][i] = new TGraphErrors(Int_t(sk_y[c][i].size()),&(sx_range_k_y[c][i][0]),&(sk_y[c][i][0]),0,&(sk_y_err[c][i][0]));
  rH_y[c][i] = new TGraphErrors(Int_t(rh_y[c][i].size()),&(rx_range_h_y[c][i][0]),&(rh_y[c][i][0]),0,&(rh_y_err[c][i][0]));
  rP_y[c][i] = new TGraphErrors(Int_t(rp_y[c][i].size()),&(rx_range_p_y[c][i][0]),&(rp_y[c][i][0]),0,&(rp_y_err[c][i][0]));
  rK_y[c][i] = new TGraphErrors(Int_t(rk_y[c][i].size()),&(rx_range_k_y[c][i][0]),&(rk_y[c][i][0]),0,&(rk_y_err[c][i][0]));

  sH_y[c][i]->SetMarkerColor(fMarkerColor[0]);
  sP_y[c][i]->SetMarkerColor(fMarkerColor[0]);
  sK_y[c][i]->SetMarkerColor(fMarkerColor[0]);
  rH_y[c][i]->SetMarkerColor(fMarkerColor[0]);
  rP_y[c][i]->SetMarkerColor(fMarkerColor[0]);
  rK_y[c][i]->SetMarkerColor(fMarkerColor[0]);

  sH_y[c][i]->SetMarkerSize(3);
  sP_y[c][i]->SetMarkerSize(3);
  sK_y[c][i]->SetMarkerSize(3);
  rH_y[c][i]->SetMarkerSize(3);
  rP_y[c][i]->SetMarkerSize(3);
  rK_y[c][i]->SetMarkerSize(3);

  sH_y[c][i]->SetMarkerStyle(fMarkerStyle[0][1]);
  sP_y[c][i]->SetMarkerStyle(fMarkerStyle[0][1]);
  sK_y[c][i]->SetMarkerStyle(fMarkerStyle[0][1]);
  rH_y[c][i]->SetMarkerStyle(fMarkerStyle[0][1]);
  rP_y[c][i]->SetMarkerStyle(fMarkerStyle[0][1]);
  rK_y[c][i]->SetMarkerStyle(fMarkerStyle[0][1]);

  sH_y[c][i]->SetTitle("");
  sP_y[c][i]->SetTitle("");
  sK_y[c][i]->SetTitle("");
  rH_y[c][i]->SetTitle("");
  rP_y[c][i]->SetTitle("");
  rK_y[c][i]->SetTitle("");

  sH_y[c][i]->GetYaxis()->SetTitle("");
  sP_y[c][i]->GetYaxis()->SetTitle("");
  sK_y[c][i]->GetYaxis()->SetTitle("");
  rH_y[c][i]->GetYaxis()->SetTitle("");
  rP_y[c][i]->GetYaxis()->SetTitle("");
  rK_y[c][i]->GetYaxis()->SetTitle("");

  sH_y[c][i]->GetXaxis()->SetTitle("z");
  sP_y[c][i]->GetXaxis()->SetTitle("z");
  sK_y[c][i]->GetXaxis()->SetTitle("z");
  rH_y[c][i]->GetXaxis()->SetTitle("z");
  rP_y[c][i]->GetXaxis()->SetTitle("z");
  rK_y[c][i]->GetXaxis()->SetTitle("z");

  c11->cd(i+1);
  gPad->SetFillStyle(4000);
  sH_y->Draw("PA");
  sH_y->GetXaxis()->SetLimits(0.001,0.5);
  sH_y->SetMinimum(0.);
  sH_y->SetMaximum(5.);
  sH_y->GetXaxis()->SetLabelSize(0.06);
  sH_y->GetYaxis()->SetLabelSize(0.06);
  sH_y->SetTitle("");
  sH_y->GetXaxis()->SetTitle("#font[ 12]{x}");
  sH_y->GetXaxis()->SetTitleSize(0.08);
  sH_y->GetXaxis()->SetTitleOffset(.8);
  sH_y->GetXaxis()->SetNdivisions(304,kTRUE);
  sH_y->GetYaxis()->SetNdivisions(304,kTRUE);
  sH_y->GetYaxis()->SetTitle("#font[12]{M}^{#font[ 12]{h^+}}+#font[12]{M}^{#font[ 12]{h^-}}");
  sH_y->GetYaxis()->SetTitleSize(0.08);
  c11->Range(0.1,0.,0.9,5.);
  c11->Update();

  c12->cd(i+1);
  gPad->SetFillStyle(4000);
  sP_y->Draw("PA");
  sP_y->GetXaxis()->SetLimits(0.001,0.5);
  sP_y->SetMinimum(0.);
  sP_y->SetMaximum(5.);
  sP_y->GetXaxis()->SetLabelSize(0.06);
  sP_y->GetYaxis()->SetLabelSize(0.06);
  sP_y->SetTitle("");
  sP_y->GetXaxis()->SetTitle("#font[ 12]{x}");
  sP_y->GetXaxis()->SetTitleSize(0.08);
  sP_y->GetXaxis()->SetTitleOffset(.8);
  sP_y->GetXaxis()->SetNdivisions(304,kTRUE);
  sP_y->GetYaxis()->SetNdivisions(304,kTRUE);
  sP_y->GetYaxis()->SetTitle("#font[12]{M}^{#font[ 12]{#pi^+}}+#font[12]{M}^{#font[ 12]{#pi^-}}");
  sP_y->GetYaxis()->SetTitleSize(0.08);
  c12->Range(0.1,0.,0.9,5.);
  c12->Update();

  c13->cd(i+1);
  gPad->SetFillStyle(4000);
  sK_y->Draw("PA");
  sK_y->GetXaxis()->SetLimits(0.001,0.5);
  sK_y->SetMinimum(0.);
  sK_y->SetMaximum(5.);
  sK_y->GetXaxis()->SetLabelSize(0.06);
  sK_y->GetYaxis()->SetLabelSize(0.06);
  sK_y->SetTitle("");
  sK_y->GetXaxis()->SetTitle("#font[ 12]{x}");
  sK_y->GetXaxis()->SetTitleSize(0.08);
  sK_y->GetXaxis()->SetTitleOffset(.8);
  sK_y->GetXaxis()->SetNdivisions(304,kTRUE);
  sK_y->GetYaxis()->SetNdivisions(304,kTRUE);
  sK_y->GetYaxis()->SetTitle("#font[12]{M}^{#font[ 12]{#K^+}}+#font[12]{M}^{#font[ 12]{#K^-}}");
  sK_y->GetYaxis()->SetTitleSize(0.08);
  c13->Range(0.1,0.,0.9,5.);
  c13->Update();

  c14->cd(i+1);
  gPad->SetFillStyle(4000);
  rH_y->Draw("PA");
  rH_y->GetXaxis()->SetLimits(0.001,0.5);
  rH_y->SetMinimum(0.);
  rH_y->SetMaximum(5.);
  rH_y->GetXaxis()->SetLabelSize(0.06);
  rH_y->GetYaxis()->SetLabelSize(0.06);
  rH_y->SetTitle("");
  rH_y->GetXaxis()->SetTitle("#font[ 12]{x}");
  rH_y->GetXaxis()->SetTitleSize(0.08);
  rH_y->GetXaxis()->SetTitleOffset(.8);
  rH_y->GetXaxis()->SetNdivisions(304,kTRUE);
  rH_y->GetYaxis()->SetNdivisions(304,kTRUE);
  rH_y->GetYaxis()->SetTitle("#font[12]{M}^{#font[ 12]{h^+}}+#font[12]{M}^{#font[ 12]{h^-}}");
  rH_y->GetYaxis()->SetTitleSize(0.08);
  c14->Range(0.1,0.,0.9,5.);
  c14->Update();

  c15->cd(i+1);
  gPad->SetFillStyle(4000);
  rP_y->Draw("PA");
  rP_y->GetXaxis()->SetLimits(0.001,0.5);
  rP_y->SetMinimum(0.);
  rP_y->SetMaximum(5.);
  rP_y->GetXaxis()->SetLabelSize(0.06);
  rP_y->GetYaxis()->SetLabelSize(0.06);
  rP_y->SetTitle("");
  rP_y->GetXaxis()->SetTitle("#font[ 12]{x}");
  rP_y->GetXaxis()->SetTitleSize(0.08);
  rP_y->GetXaxis()->SetTitleOffset(.8);
  rP_y->GetXaxis()->SetNdivisions(304,kTRUE);
  rP_y->GetYaxis()->SetNdivisions(304,kTRUE);
  rP_y->GetYaxis()->SetTitle("#font[12]{M}^{#font[ 12]{#pi^+}}+#font[12]{M}^{#font[ 12]{#pi^-}}");
  rP_y->GetYaxis()->SetTitleSize(0.08);
  c15->Range(0.1,0.,0.9,5.);
  c15->Update();

  c16->cd(i+1);
  gPad->SetFillStyle(4000);
  rK_y->Draw("PA");
  rK_y->GetXaxis()->SetLimits(0.001,0.5);
  rK_y->SetMinimum(0.);
  rK_y->SetMaximum(5.);
  rK_y->GetXaxis()->SetLabelSize(0.06);
  rK_y->GetYaxis()->SetLabelSize(0.06);
  rK_y->SetTitle("");
  rK_y->GetXaxis()->SetTitle("#font[ 12]{x}");
  rK_y->GetXaxis()->SetTitleSize(0.08);
  rK_y->GetXaxis()->SetTitleOffset(.8);
  rK_y->GetXaxis()->SetNdivisions(304,kTRUE);
  rK_y->GetYaxis()->SetNdivisions(304,kTRUE);
  rK_y->GetYaxis()->SetTitle("#font[12]{M}^{#font[ 12]{#K^+}}+#font[12]{M}^{#font[ 12]{#K^-}}");
  rK_y->GetYaxis()->SetTitleSize(0.08);
  c16->Range(0.1,0.,0.9,5.);
  c16->Update();

  }

  TLatex fTitle;

  c5->cd(1);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 3.7,"0.004#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.01");

  c5->cd(2);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 3.7,"0.01#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.02");

  c5->cd(3);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 3.7,"0.02#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.03");

  c5->cd(4);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 3.7,"0.03#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.04");

  c5->cd(5);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 3.7,"0.04#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.06");

  c5->cd(6);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 3.7,"0.06#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.1");

  c5->cd(7);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 3.7,"0.1#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.14");

  c5->cd(8);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 3.7,"0.14#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.18");

  c5->cd(9);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 3.7,"0.18#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.4");

  c5->cd(10);
  fTitle.SetTextSize(0.095);
  fTitle.SetTextAlign(11);
  fTitle.DrawLatex(0.05, 0.72,"#color[221]{0.70#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.90, #delta = 0.5}");
  fTitle.DrawLatex(0.05, 0.64,"#color[4]{0.50#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.70, #delta = 0.4}");
  fTitle.DrawLatex(0.05, 0.56,"#color[226]{0.30#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.50, #delta = 0.3}");
  fTitle.DrawLatex(0.05, 0.48,"#color[209]{0.20#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.30, #delta = 0.2}");
  fTitle.DrawLatex(0.05, 0.40,"#color[95]{0.15#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.20, #delta = 0.1}");
  fTitle.DrawLatex(0.05, 0.32,"#color[2]{0.10#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.15, #delta = 0}");


  c6->cd(1);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.9,"0.004#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.01");

  c6->cd(2);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.9,"0.01#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.02");

  c6->cd(3);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.9,"0.02#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.03");

  c6->cd(4);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.9,"0.03#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.04");

  c6->cd(5);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.9,"0.04#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.06");

  c6->cd(6);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.9,"0.06#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.1");

  c6->cd(7);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.9,"0.1#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.14");

  c6->cd(8);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.9,"0.14#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.18");

  c6->cd(9);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.9,"0.18#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.4");

  c6->cd(10);
  fTitle.SetTextSize(0.095);
  fTitle.SetTextAlign(11);
  fTitle.DrawLatex(0.05, 0.72,"#color[221]{0.70#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.90, #delta = 0.5}");
  fTitle.DrawLatex(0.05, 0.64,"#color[4]{0.50#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.70, #delta = 0.4}");
  fTitle.DrawLatex(0.05, 0.56,"#color[226]{0.30#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.50, #delta = 0.3}");
  fTitle.DrawLatex(0.05, 0.48,"#color[209]{0.20#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.30, #delta = 0.2}");
  fTitle.DrawLatex(0.05, 0.40,"#color[95]{0.15#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.20, #delta = 0.1}");
  fTitle.DrawLatex(0.05, 0.32,"#color[2]{0.10#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.15, #delta = 0}");

  c7->cd(1);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.4,"0.004#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.01");

  c7->cd(2);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.4,"0.01#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.02");

  c7->cd(3);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.4,"0.02#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.03");

  c7->cd(4);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.4,"0.03#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.04");

  c7->cd(5);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.4,"0.04#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.06");

  c7->cd(6);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.4,"0.06#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.1");

  c7->cd(7);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.4,"0.1#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.14");

  c7->cd(8);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.4,"0.14#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.18");

  c7->cd(9);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.4,"0.18#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.4");

  c7->cd(10);
  fTitle.SetTextSize(0.095);
  fTitle.SetTextAlign(11);
  fTitle.DrawLatex(0.05, 0.72,"#color[221]{0.70#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.90, #delta = 0.5}");
  fTitle.DrawLatex(0.05, 0.64,"#color[4]{0.50#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.70, #delta = 0.4}");
  fTitle.DrawLatex(0.05, 0.56,"#color[226]{0.30#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.50, #delta = 0.3}");
  fTitle.DrawLatex(0.05, 0.48,"#color[209]{0.20#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.30, #delta = 0.2}");
  fTitle.DrawLatex(0.05, 0.40,"#color[95]{0.15#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.20, #delta = 0.1}");
  fTitle.DrawLatex(0.05, 0.32,"#color[2]{0.10#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.15, #delta = 0}");

  c8->cd(1);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 4.4,"0.004#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.01");

  c8->cd(2);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 4.4,"0.01#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.02");

  c8->cd(3);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 4.4,"0.02#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.03");

  c8->cd(4);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 4.4,"0.03#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.04");

  c8->cd(5);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 4.4,"0.04#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.06");

  c8->cd(6);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 4.4,"0.06#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.1");

  c8->cd(7);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 4.4,"0.1#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.14");

  c8->cd(8);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 4.4,"0.14#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.18");

  c8->cd(9);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 4.4,"0.18#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.4");

  c9->cd(1);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.4,"0.004#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.01");

  c9->cd(2);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.4,"0.01#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.02");

  c9->cd(3);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.4,"0.02#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.03");

  c9->cd(4);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.4,"0.03#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.04");

  c9->cd(5);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.4,"0.04#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.06");

  c9->cd(6);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.4,"0.06#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.1");

  c9->cd(7);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.4,"0.1#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.14");

  c9->cd(8);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.4,"0.14#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.18");

  c9->cd(9);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.4,"0.18#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.4");

  c10->cd(1);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.4,"0.004#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.01");

  c10->cd(2);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.4,"0.01#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.02");

  c10->cd(3);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.4,"0.02#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.03");

  c10->cd(4);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.4,"0.03#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.04");

  c10->cd(5);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.4,"0.04#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.06");

  c10->cd(6);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.4,"0.06#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.1");

  c10->cd(7);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.4,"0.1#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.14");

  c10->cd(8);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.4,"0.14#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.18");

  c10->cd(9);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.4,"0.18#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.4");


  c5->Update();
  c6->Update();
  c7->Update();
  c8->Update();
  c9->Update();
  c10->Update();
  c11->Update();
  c12->Update();
  c13->Update();
  c14->Update();
  c15->Update();
  c16->Update();

  c5->Print(Form("%s/hadron_multiplicity_file.pdf",data_path));
  c5->Print(Form("%s/hadron_multiplicity_file.root",data_path));
  c6->Print(Form("%s/pion_multiplicity_file.pdf",data_path));
  c6->Print(Form("%s/pion_multiplicity_file.root",data_path));
  c7->Print(Form("%s/kaon_multiplicity_file.pdf",data_path));
  c7->Print(Form("%s/kaon_multiplicity_file.root",data_path));
  c8->Print(Form("%s/hadron_multiplicity_yavg_file.pdf",data_path));
  c8->Print(Form("%s/hadron_multiplicity_yavg_file.root",data_path));
  c9->Print(Form("%s/pion_multiplicity_yavg_file.pdf",data_path));
  c9->Print(Form("%s/pion_multiplicity_yavg_file.root",data_path));
  c10->Print(Form("%s/kaon_multiplicity_yavg_file.pdf",data_path));
  c10->Print(Form("%s/kaon_multiplicity_yavg_file.root",data_path));
  c11->Print(Form("%s/hadron_multiplicity_sum_file.pdf",data_path));
  c11->Print(Form("%s/hadron_multiplicity_sum_file.root",data_path));
  c12->Print(Form("%s/pion_multiplicity_sum_file.pdf",data_path));
  c12->Print(Form("%s/pion_multiplicity_sum_file.root",data_path));
  c13->Print(Form("%s/kaon_multiplicity_sum_file.pdf",data_path));
  c13->Print(Form("%s/kaon_multiplicity_sum_file.root",data_path));
  c14->Print(Form("%s/hadron_multiplicity_ratio_file.pdf",data_path));
  c14->Print(Form("%s/hadron_multiplicity_ratio_file.root",data_path));
  c15->Print(Form("%s/pion_multiplicity_ratio_file.pdf",data_path));
  c15->Print(Form("%s/pion_multiplicity_ratio_file.root",data_path));
  c16->Print(Form("%s/kaon_multiplicity_ratio_file.pdf",data_path));
  c16->Print(Form("%s/kaon_multiplicity_ratio_file.root",data_path));

  ofs_p.close();
  ofs_k.close();
  ofs_h.close();

  return 0;
}

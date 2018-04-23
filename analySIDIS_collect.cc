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

#define dirroot "/sps/compass/npierre/Multiplicities"

// Flags
#define Y2006 0
#define Y2012 0
#define Y2016 1
#define RCUTSTUDY_ON 0

using namespace std;

void fetch_acceptance(string pname)
{
  ifstream acc_file(pname);
  double dummy;

  for(int c=0; c<2; c++)
  {
    for(int i=0; i<9; i++)
    {
      for(int j=0; j<5; j++)
      {
        for(int k=0; k<12; k++)
        {
          for(int l=0; l<4; l++)
          {
            acc_file >> dummy;
          }
          acc_file >> fAcceptance[i][j][k].tab[c][0][0];
          acc_file >> fAcceptance[i][j][k].tab[c][1][0];
          acc_file >> fAcceptance[i][j][k].tab[c][0][1];
          acc_file >> fAcceptance[i][j][k].tab[c][1][1];
          acc_file >> fAcceptance[i][j][k].tab[c][0][2];
          acc_file >> fAcceptance[i][j][k].tab[c][1][2];
          acc_file >> fAcceptance[i][j][k].tab[c][0][3];
          acc_file >> fAcceptance[i][j][k].tab[c][1][3];
#ifdef DEBUG
          cout << fAcceptance[i][j][k].tab[c][0][0] << " " <<
          fAcceptance[i][j][k].tab[c][1][0] << " " <<
          fAcceptance[i][j][k].tab[c][0][1] << " " <<
          fAcceptance[i][j][k].tab[c][1][1] << " " <<
          fAcceptance[i][j][k].tab[c][0][2] << " " <<
          fAcceptance[i][j][k].tab[c][1][2] << " " <<
          fAcceptance[i][j][k].tab[c][0][3] << " " <<
          fAcceptance[i][j][k].tab[c][1][3] << endl;
#endif
        }
      }
    }
  }

  acc_file.close();

}

int main()
{

  //q_bin x_bin y_bin z_bin acc_pi acc_error_pi acc_k acc_error_k acc_p acc_error_p acc_h acc_error_h

  int year=0;

  if(Y2006) year=2006;
  else if(Y2012) year=2012;
  else if(Y2016) year=2016;

  double dummy;

  fetch_acceptance(Form("acceptance/%d/acceptance.txt",year));

  for(int filen=0; filen<1/*5*/; filen++)
  {
    ifstream dis_file(Form("rawmult/%d/DIS_%d.txt",year,filen));
    ifstream had_file(Form("rawmult/%d/hadron_%d.txt",year,filen));

    for(int c=0; c<2; c++)
    {
      for(int i=0; i<9; i++)
      {
        for(int j=0; j<5; j++)
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
              fMeanvalues_data[i][j][k].tab[c][ll][0] += dummy;
              dis_file >> dummy;
              fMeanvalues_data[i][j][k].tab[c][ll][1] += dummy;
              dis_file >> dummy;
              fMeanvalues_data[i][j][k].tab[c][ll][2] += dummy;
              dis_file >> dummy;
              fMeanvalues_data[i][j][k].tab[c][ll][3] += dummy;
              dis_file >> dummy;
              fMeanvalues_size[i][j][k].tab[c][ll][0] += dummy;
            }

            had_file >> dummy;
            fBinning[i][j][k].tab[c][0][0] += dummy;
            had_file >> dummy;
            fBinning[i][j][k].tab[c][1][0] += dummy;
            had_file >> dummy;
            fBinning_loose[i][j][k].tab[c][0][0] += dummy;
            had_file >> dummy;
            fBinning_severe[i][j][k].tab[c][0][0] += dummy;
            had_file >> dummy;
            fBinning[i][j][k].tab[c][0][1] += dummy;
            had_file >> dummy;
            fBinning[i][j][k].tab[c][1][1] += dummy;
            had_file >> dummy;
            fBinning_loose[i][j][k].tab[c][0][1] += dummy;
            had_file >> dummy;
            fBinning_severe[i][j][k].tab[c][0][1] += dummy;
            had_file >> dummy;
            fBinning[i][j][k].tab[c][0][2] += dummy;
            had_file >> dummy;
            fBinning[i][j][k].tab[c][1][2] += dummy;
            had_file >> dummy;
            fBinning_loose[i][j][k].tab[c][0][2] += dummy;
            had_file >> dummy;
            fBinning_severe[i][j][k].tab[c][0][2] += dummy;
            had_file >> dummy;
            fBinning[i][j][k].tab[c][0][3] += dummy;
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
  }

  double rich_sys_err;

  for(int c=0; c<2; c++)
  {
    for(xbin=0; xbin<9; xbin++)
    {
      for(ybin=0; ybin<5; ybin++)
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

  for(int i=0; i<9; i++)
  {
    for(int j=0; j<5; j++)
    {
      for(int k=0; k<12; k++)
      {
        for(int c=0; c<2; c++)
        {
          for(int ll=0; ll<4; ll++)
          {
            if(int(fMeanvalues_size[i][j][k].tab[c][ll][0]))
            {
              fMeanvalues_data[i][j][k].tab[c][ll][0] /= int(fMeanvalues_size[i][j][k].tab[c][ll][0]);
              fMeanvalues_data[i][j][k].tab[c][ll][1] /= int(fMeanvalues_size[i][j][k].tab[c][ll][0]);
              fMeanvalues_data[i][j][k].tab[c][ll][2] /= int(fMeanvalues_size[i][j][k].tab[c][ll][0]);
              fMeanvalues_data[i][j][k].tab[c][ll][3] /= int(fMeanvalues_size[i][j][k].tab[c][ll][0]);
            }
          }
        }
      }
    }
  }

  TCanvas* c5;
  c5 = new TCanvas("Hadron_Multiplicities","Hadron_Multiplicities",3200,1600);

  TCanvas* c6;
  c6 = new TCanvas("Pion_Multiplicities","Pion_Multiplicities",3200,1600);

  TCanvas* c7;
  c7 = new TCanvas("Kaon_Multiplicities","Kaon_Multiplicities",3200,1600);

  c5->SetFillColor(0);
  //c5->SetFrameFillStyle(4000);
  c6->SetFillColor(0);
  //c6->SetFrameFillStyle(4000);
  c7->SetFillColor(0);
  //c7->SetFrameFillStyle(4000);

  c5->Divide(5,2,0,0);
  c6->Divide(5,2,0,0);
  c7->Divide(5,2,0,0);

  TGraphErrors* H_mult[2][9][5];
  TGraphErrors* P_mult[2][9][5];
  TGraphErrors* K_mult[2][9][5];

  Double_t z_range[12] = {.225,.275,.325,.375,.425,.475,.525,.575,.625,.675,.725,.8};

  ofstream ofs_p(Form("%s/multiplicities_pion.txt",dirroot), std::ofstream::out | std::ofstream::trunc);
  ofstream ofs_t(Form("%s/multiplicities_raw.txt",dirroot), std::ofstream::out | std::ofstream::trunc);
  ofstream ofs_k(Form("%s/multiplicities_kaon.txt",dirroot), std::ofstream::out | std::ofstream::trunc);
  ofstream ofs_h(Form("%s/multiplicities_hadron.txt",dirroot), std::ofstream::out | std::ofstream::trunc);

  std::vector<Double_t> p_m[2][9][5];
  std::vector<Double_t> k_m[2][9][5];
  std::vector<Double_t> h_m[2][9][5];
  std::vector<Double_t> p_err[2][9][5];
  std::vector<Double_t> k_err[2][9][5];
  std::vector<Double_t> h_err[2][9][5];
  std::vector<Double_t> z_range_p[2][9][5];
  std::vector<Double_t> z_range_k[2][9][5];
  std::vector<Double_t> z_range_h[2][9][5];

  for(int i=0; i<9; i++)
  {
    for(int j=0; j<5; j++)
    {
      for(int k=0; k<12; k++)
      {

        for(int c=1; c>=0; c--)
        {
          /*ofs << c << " " << i << " " << j << " " << k << "   " << (fNDIS_evt[0][i][j][k] ? Double_t(fBinning[i][j][k].tab[c][0][3]/(fNDIS_evt[0][i][j][k]*fZ_bin_width[k])) : 0) << " " <<
          fNDIS_evt[0][i][j][k] << " " << fNDIS_evt[0][i][j][k] << " " <<
          fBinning[i][j][k].tab[c][0][0] << " " << fBinning[i][j][k].tab[c][1][0] << " " <<
          fNDIS_evt[1][i][j][k] << " " << fNDIS_evt[1][i][j][k] << " " <<
          fBinning[i][j][k].tab[c][0][1] << " " << fBinning[i][j][k].tab[c][1][1] << " " <<
          fNDIS_evt[2][i][j][k] << " " << fNDIS_evt[2][i][j][k] << " " <<
          fBinning[i][j][k].tab[c][0][2] << " " << fBinning[i][j][k].tab[c][1][2] << " " <<
          fNDIS_evt[0][i][j][k] << " " << fNDIS_evt[0][i][j][k] << " " <<
          fBinning[i][j][k].tab[c][0][3] << " " << fBinning[i][j][k].tab[c][1][3] << endl;*/

          if((i==7 && j==4) || (i==8 && j==0) || (i==8 && j==4))
          {
            if(c) ofs_p << fXrange[i] << " " << fYrange[j] << " " << fZrange[k] << " ";

            ofs_p <<
            "0" << " " << "0" << " " <<
            "0" << " " << "0" << " " <<
            "0" << " " <<
            "0" << " " <<
            "0" << " " <<
            "0" << " ";

            if(!c) ofs_p << endl;

            if(c) ofs_k << fXrange[i] << " " << fYrange[j] << " " << fZrange[k] << " ";

            ofs_k <<
            "0" << " " << "0" << " " <<
            "0" << " " << "0" << " " <<
            "0" << " " <<
            "0" << " " <<
            "0" << " " <<
            "0" << " ";

            if(!c) ofs_k << endl;

            if(c) ofs_h << fXrange[i] << " " << fYrange[j] << " " << fZrange[k] << " ";

            ofs_h <<
            "0" << " " << "0" << " " <<
            "0" << " " << "0" << " " <<
            "0" << " " <<
            "0" << " " <<
            "0" << " " <<
            "0" << " ";

            if(!c) ofs_h << endl;

            p_m[c][i][j].push_back(0);
            k_m[c][i][j].push_back(0);
            h_m[c][i][j].push_back(0);
            p_err[c][i][j].push_back(0);
            k_err[c][i][j].push_back(0);
            h_err[c][i][j].push_back(0);
          }
          else if(!(fNDIS_evt[0][i][j][k] && fAcceptance[i][j][k].tab[c][0][0] ? Double_t(fBinning[i][j][k].tab[c][0][0]/(fNDIS_evt[0][i][j][k]*fZ_bin_width[k]*fAcceptance[i][j][k].tab[c][0][0])) : 0))
          {
            if(c) ofs_p << fXrange[i] << " " << fYrange[j] << " " << fZrange[k] << " ";

            ofs_p <<
            "0" << " " << "0" << " " <<
            "0" << " " << "0" << " " <<
            "0" << " " <<
            "0" << " " <<
            "0" << " " <<
            "0" << " ";

            if(!c) ofs_p << endl;

            if(c) ofs_k << fXrange[i] << " " << fYrange[j] << " " << fZrange[k] << " ";

            ofs_k <<
            "0" << " " << "0" << " " <<
            "0" << " " << "0" << " " <<
            "0" << " " <<
            "0" << " " <<
            "0" << " " <<
            "0" << " ";

            if(!c) ofs_k << endl;

            if(c) ofs_h << fXrange[i] << " " << fYrange[j] << " " << fZrange[k] << " ";

            ofs_h <<
            "0" << " " << "0" << " " <<
            "0" << " " << "0" << " " <<
            "0" << " " <<
            "0" << " " <<
            "0" << " " <<
            "0" << " ";

            if(!c) ofs_h << endl;

            p_m[c][i][j].push_back(0);
            k_m[c][i][j].push_back(0);
            h_m[c][i][j].push_back(0);
            p_err[c][i][j].push_back(0);
            k_err[c][i][j].push_back(0);
            h_err[c][i][j].push_back(0);
          }
          else
          {
            if(c) ofs_p << fXrange[i] << " " << fYrange[j] << " " << fZrange[k] << " ";

            ofs_p <<
            fMeanvalues_data[i][j][k].tab[0][0][0] << " " << fMeanvalues_data[i][j][k].tab[0][0][1] << " " <<
            fMeanvalues_data[i][j][k].tab[0][0][2] << " " << fMeanvalues_data[i][j][k].tab[0][0][3] << " " <<
            (fNDIS_evt[0][i][j][k] && fAcceptance[i][j][k].tab[c][0][0] ? Double_t(fBinning[i][j][k].tab[c][0][0]/(fNDIS_evt[0][i][j][k]*fZ_bin_width[k]*fAcceptance[i][j][k].tab[c][0][0])) : 0) << " " <<
            (fNDIS_evt[0][i][j][k] && fAcceptance[i][j][k].tab[c][0][0] ? Double_t(((fBinning[i][j][k].tab[c][1][0]/pow(fNDIS_evt[0][i][j][k],2)-pow(fBinning[i][j][k].tab[c][0][0],2)*fNDIS_evt_err[0][i][j][k]/pow(fNDIS_evt[0][i][j][k],4))/(pow(fZ_bin_width[k]*fAcceptance[i][j][k].tab[c][0][0],2)))
                                                                                    + fAcceptance[i][j][k].tab[c][1][0]*pow(fBinning[i][j][k].tab[c][0][0]/(fNDIS_evt[0][i][j][k]*fZ_bin_width[k]*pow(fAcceptance[i][j][k].tab[c][0][0],2)),2)) : 0) << " " <<
            (fNDIS_evt[0][i][j][k] ? Double_t(sqrt(pow(fRich_sys_err[i][j][k].tab[c][1][0],2)/pow(fNDIS_evt[0][i][j][k]*fZ_bin_width[k]*fAcceptance[i][j][k].tab[c][0][0],2)+pow(0.05*sqrt(fAcceptance[i][j][k].tab[c][1][0])*fBinning[i][j][k].tab[c][0][0]/(fNDIS_evt[0][i][j][k]*fZ_bin_width[k]*pow(fAcceptance[i][j][k].tab[c][0][0],2)),2))) : 0) << " " <<
            ((fNDIS_evt[0][i][j][k] && fAcceptance[i][j][k].tab[c][0][0] ? Double_t(fBinning[i][j][k].tab[c][0][0]/(fNDIS_evt[0][i][j][k]*fZ_bin_width[k]*fAcceptance[i][j][k].tab[c][0][0])) : 0) ? 1 : 0) << " ";

            if(!c) ofs_p << endl;

            if(c) ofs_t << fXrange[i] << " " << fYrange[j] << " " << fZrange[k] << " ";

            ofs_t <<
            fMeanvalues_data[i][j][k].tab[0][0][0] << " " << fMeanvalues_data[i][j][k].tab[0][0][1] << " " <<
            fMeanvalues_data[i][j][k].tab[0][0][2] << " " << fMeanvalues_data[i][j][k].tab[0][0][3] << " " <<
            (fNDIS_evt[0][i][j][k] && fAcceptance[i][j][k].tab[c][0][0] ? Double_t(fBinning[i][j][k].tab[c][0][0]/(fNDIS_evt[0][i][j][k]*fZ_bin_width[k])) : 0) << " ";

            if(!c) ofs_t << endl;

            if (c) ofs_k << fXrange[i] << " " << fYrange[j] << " " << fZrange[k] << " ";

            ofs_k <<
            fMeanvalues_data[i][j][k].tab[0][1][0] << " " << fMeanvalues_data[i][j][k].tab[0][1][1] << " " <<
            fMeanvalues_data[i][j][k].tab[0][1][2] << " " << fMeanvalues_data[i][j][k].tab[0][1][3] << " " <<
            (fNDIS_evt[0][i][j][k] && fAcceptance[i][j][k].tab[c][0][1] ? Double_t(fBinning[i][j][k].tab[c][0][1]/(fNDIS_evt[0][i][j][k]*fZ_bin_width[k]*fAcceptance[i][j][k].tab[c][0][1])) : 0) << " " <<
            (fNDIS_evt[0][i][j][k] && fAcceptance[i][j][k].tab[c][0][1] ? Double_t(((fBinning[i][j][k].tab[c][1][1]/pow(fNDIS_evt[0][i][j][k],2)-pow(fBinning[i][j][k].tab[c][0][0],2)*fNDIS_evt_err[0][i][j][k]/pow(fNDIS_evt[0][i][j][k],4))/(pow(fZ_bin_width[k]*fAcceptance[i][j][k].tab[c][0][1],2)))
                                                                                    + fAcceptance[i][j][k].tab[c][1][1]*pow(fBinning[i][j][k].tab[c][0][1]/(fNDIS_evt[0][i][j][k]*fZ_bin_width[k]*pow(fAcceptance[i][j][k].tab[c][0][1],2)),2)) : 0) << " " <<
            (fNDIS_evt[0][i][j][k] ? Double_t(sqrt(pow(fRich_sys_err[i][j][k].tab[c][1][1],2)/pow(fNDIS_evt[0][i][j][k]*fZ_bin_width[k]*fAcceptance[i][j][k].tab[c][0][1],2)+pow(0.05*sqrt(fAcceptance[i][j][k].tab[c][1][1])*fBinning[i][j][k].tab[c][0][1]/(fNDIS_evt[0][i][j][k]*fZ_bin_width[k]*pow(fAcceptance[i][j][k].tab[c][0][1],2)),2))) : 0) << " " <<
            ((fNDIS_evt[0][i][j][k] && fAcceptance[i][j][k].tab[c][0][1] ? Double_t(fBinning[i][j][k].tab[c][0][1]/(fNDIS_evt[0][i][j][k]*fZ_bin_width[k]*fAcceptance[i][j][k].tab[c][0][1])) : 0) ? 1 : 0) << " ";

            if(!c) ofs_k << endl;

            if (c) ofs_h << fXrange[i] << " " << fYrange[j] << " " << fZrange[k] << " ";

            ofs_h <<
            fMeanvalues_data[i][j][k].tab[0][3][0] << " " << fMeanvalues_data[i][j][k].tab[0][3][1] << " " <<
            fMeanvalues_data[i][j][k].tab[0][3][2] << " " << fMeanvalues_data[i][j][k].tab[0][3][3] << " " <<
            (fNDIS_evt[0][i][j][k] && fAcceptance[i][j][k].tab[c][0][3] ? Double_t(fBinning[i][j][k].tab[c][0][3]/(fNDIS_evt[0][i][j][k]*fZ_bin_width[k]*fAcceptance[i][j][k].tab[c][0][3])) : 0) << " " <<
            (fNDIS_evt[0][i][j][k] && fAcceptance[i][j][k].tab[c][0][3] ? Double_t(((fBinning[i][j][k].tab[c][1][3]/pow(fNDIS_evt[0][i][j][k],2)-pow(fBinning[i][j][k].tab[c][0][3],2)*fNDIS_evt_err[0][i][j][k]/pow(fNDIS_evt[0][i][j][k],4))/(pow(fZ_bin_width[k]*fAcceptance[i][j][k].tab[c][0][3],2)))
                                                                                    + fAcceptance[i][j][k].tab[c][1][0]*pow(fBinning[i][j][k].tab[c][0][0]/(fNDIS_evt[0][i][j][k]*fZ_bin_width[k]*pow(fAcceptance[i][j][k].tab[c][0][0],2)),2)) : 0) << " " <<
            (fNDIS_evt[0][i][j][k] ? Double_t(sqrt(pow(fRich_sys_err[i][j][k].tab[c][1][3],2)/pow(fNDIS_evt[0][i][j][k]*fZ_bin_width[k]*fAcceptance[i][j][k].tab[c][0][3],2)+pow(0.05*sqrt(fAcceptance[i][j][k].tab[c][1][3])*fBinning[i][j][k].tab[c][0][3]/(fNDIS_evt[0][i][j][k]*fZ_bin_width[k]*pow(fAcceptance[i][j][k].tab[c][0][3],2)),2))) : 0) << " " <<
            ((fNDIS_evt[0][i][j][k] && fAcceptance[i][j][k].tab[c][0][3] ? Double_t(fBinning[i][j][k].tab[c][0][3]/(fNDIS_evt[0][i][j][k]*fZ_bin_width[k]*fAcceptance[i][j][k].tab[c][0][3])) : 0) ? 1 : 0) << " ";

            if(!c) ofs_h << endl;

            p_m[c][i][j].push_back((fNDIS_evt[0][i][j][k] && fAcceptance[i][j][k].tab[c][0][0] ? Double_t(fBinning[i][j][k].tab[c][0][0]/(fNDIS_evt[0][i][j][k]*fZ_bin_width[k]*fAcceptance[i][j][k].tab[c][0][0])+j*0.1) : 0));
            k_m[c][i][j].push_back((fNDIS_evt[0][i][j][k] && fAcceptance[i][j][k].tab[c][0][1] ? Double_t(fBinning[i][j][k].tab[c][0][1]/(fNDIS_evt[0][i][j][k]*fZ_bin_width[k]*fAcceptance[i][j][k].tab[c][0][1])+j*0.1) : 0));
            h_m[c][i][j].push_back((fNDIS_evt[0][i][j][k] && fAcceptance[i][j][k].tab[c][0][3] ? Double_t(fBinning[i][j][k].tab[c][0][3]/(fNDIS_evt[0][i][j][k]*fZ_bin_width[k]*fAcceptance[i][j][k].tab[c][0][3])+j*0.1) : 0));
            p_err[c][i][j].push_back((fNDIS_evt[0][i][j][k] && fAcceptance[i][j][k].tab[c][0][0] ? Double_t(((fBinning[i][j][k].tab[c][1][0]/pow(fNDIS_evt[0][i][j][k],2)-pow(fBinning[i][j][k].tab[c][0][0],2)*fNDIS_evt_err[0][i][j][k]/pow(fNDIS_evt[0][i][j][k],4))/(pow(fZ_bin_width[k]*fAcceptance[i][j][k].tab[c][0][0],2)))
                                                                                    + fAcceptance[i][j][k].tab[c][1][0]*pow(fBinning[i][j][k].tab[c][0][0]/(fNDIS_evt[0][i][j][k]*fZ_bin_width[k]*pow(fAcceptance[i][j][k].tab[c][0][0],2)),2)) : 0));
            k_err[c][i][j].push_back((fNDIS_evt[0][i][j][k] && fAcceptance[i][j][k].tab[c][0][1] ? Double_t(((fBinning[i][j][k].tab[c][1][1]/pow(fNDIS_evt[0][i][j][k],2)-pow(fBinning[i][j][k].tab[c][0][1],2)*fNDIS_evt_err[0][i][j][k]/pow(fNDIS_evt[0][i][j][k],4))/(pow(fZ_bin_width[k]*fAcceptance[i][j][k].tab[c][0][1],2)))
                                                                                    + fAcceptance[i][j][k].tab[c][1][1]*pow(fBinning[i][j][k].tab[c][0][1]/(fNDIS_evt[0][i][j][k]*fZ_bin_width[k]*pow(fAcceptance[i][j][k].tab[c][0][1],2)),2)) : 0)););
            h_err[c][i][j].push_back((fNDIS_evt[0][i][j][k] && fAcceptance[i][j][k].tab[c][0][3] ? Double_t(((fBinning[i][j][k].tab[c][1][3]/pow(fNDIS_evt[0][i][j][k],2)-pow(fBinning[i][j][k].tab[c][0][3],2)*fNDIS_evt_err[0][i][j][k]/pow(fNDIS_evt[0][i][j][k],4))/(pow(fZ_bin_width[k]*fAcceptance[i][j][k].tab[c][0][3],2)))
                                                                                    + fAcceptance[i][j][k].tab[c][1][3]*pow(fBinning[i][j][k].tab[c][0][3]/(fNDIS_evt[0][i][j][k]*fZ_bin_width[k]*pow(fAcceptance[i][j][k].tab[c][0][3],2)),2)) : 0)););
          }
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
            if(!c && j==3)
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
              H_mult[c][i][0]->Draw("SAMEP");
              H_mult[c][i][0]->GetXaxis()->SetLimits(0.1,0.9);
              H_mult[c][i][0]->SetMinimum(0.);
              H_mult[c][i][0]->SetMaximum(4.);
              H_mult[c][i][1]->Draw("SAMEP");
              H_mult[c][i][1]->GetXaxis()->SetLimits(0.1,0.9);
              H_mult[c][i][1]->SetMinimum(0.);
              H_mult[c][i][1]->SetMaximum(4.);
              H_mult[c][i][2]->Draw("SAMEP");
              H_mult[c][i][2]->GetXaxis()->SetLimits(0.1,0.9);
              H_mult[c][i][2]->SetMinimum(0.);
              H_mult[c][i][2]->SetMaximum(4.);
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
      }
    }
  }

  TLatex fTitle;

  c5->cd(1);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 2.3,"0.004#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.01");

  c5->cd(2);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 2.3,"0.01#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.02");

  c5->cd(3);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 2.3,"0.02#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.03");

  c5->cd(4);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 2.3,"0.03#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.04");

  c5->cd(5);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 2.3,"0.04#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.06");

  c5->cd(6);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 2.3,"0.06#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.1");

  c5->cd(7);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 2.3,"0.1#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.14");

  c5->cd(8);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 2.3,"0.14#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.18");

  c5->cd(9);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 2.3,"0.18#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.4");

  c5->cd(10);
  fTitle.SetTextSize(0.095);
  fTitle.SetTextAlign(11);
  fTitle.DrawLatex(0.05, 0.64,"#color[221]{0.50#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.70, #delta = 0.4}");
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
  fTitle.DrawLatex(0.05, 0.64,"#color[221]{0.50#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.70, #delta = 0.4}");
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
  fTitle.DrawLatex(0.05, 0.64,"#color[221]{0.50#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.70, #delta = 0.4}");
  fTitle.DrawLatex(0.05, 0.56,"#color[226]{0.30#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.50, #delta = 0.3}");
  fTitle.DrawLatex(0.05, 0.48,"#color[209]{0.20#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.30, #delta = 0.2}");
  fTitle.DrawLatex(0.05, 0.40,"#color[95]{0.15#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.20, #delta = 0.1}");
  fTitle.DrawLatex(0.05, 0.32,"#color[2]{0.10#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.15, #delta = 0}");


  c5->Update();
  c6->Update();
  c7->Update();

  c5->Print(Form("%s/hadron_multiplicity_file.pdf",dirroot));
  c5->Print(Form("%s/hadron_multiplicity_file.root",dirroot));
  c6->Print(Form("%s/pion_multiplicity_file.pdf",dirroot));
  c6->Print(Form("%s/pion_multiplicity_file.root",dirroot));
  c7->Print(Form("%s/kaon_multiplicity_file.pdf",dirroot));
  c7->Print(Form("%s/kaon_multiplicity_file.root",dirroot));

  ofs_p.close();
  ofs_k.close();
  ofs_h.close();

  return 0;
}

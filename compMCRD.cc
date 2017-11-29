#include <iostream>
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

#include "compMCRD.h"

//Inputs
#define mat_RICH_name "rich_mat.txt"
#define err_RICH_name "rich_mat_error.txt"
#define target_file_2012 "target-107924-109081.dat"
#define target_file_2016 "target-274508-274901.dat"

// Flags
#define Y2006 0
#define Y2012 0
#define Y2016 1
#define MOMENTUM 3
#define RCUTSTUDY_ON 0
//#define SIZESPLIT 1
//#define OFFSET 0

using namespace std;

// Progress bar

# define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
# define PBWIDTH 60

void printProgress(int event, int total)
{
    string points[6] = {"   ",".  ",".. ","..."," ..","  ."};
    double percentage = double(event)/double(total);
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf ("\r Progress%s %3d%% [%.*s%*s] (%d/%d)", points[int(event/24)%6].c_str(), val, lpad, PBSTR, rpad, "", event, total);
    fflush (stdout);
}

// Target Management

void InitTargetFile(string pfile)
{
  char tstr[500];
  std::ifstream fin;
  sprintf(tstr,pfile.c_str());
  cout<<"INFO : Opening target cell description: "<<tstr<<"..."<<endl;
  fin.open(tstr);
  while(fin.is_open() && !fin.eof())
  {
    float z, x, y, dummy;
    fin >> z >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> x >> y;
    fZv.push_back(z);
    fXv.push_back(x);
    fYv.push_back(y);
  }
  cout<<"INFO : Target cell description loaded"<<endl;
}

void CellCenter(Double_t z, Double_t& xc, Double_t& yc)
{
  xc = 1000000;
  yc = 1000000;

  for(int i = 0; i < int(fZv.size()-1); i++)
  {
    Double_t z1 = fZv[i];
    Double_t z2 = fZv[i+1];

    if( z2 < z ) continue;
    if( z1 > z ) continue;

    Double_t xc1 = fXv[i];
    Double_t xc2 = fXv[i+1];

    Double_t yc1 = fYv[i];
    Double_t yc2 = fYv[i+1];

    Double_t dxcdz = (xc2-xc1)/(z2-z1);
    Double_t dycdz = (yc2-yc1)/(z2-z1);

    Double_t dz = z-z1;
    xc = xc1 + dxcdz*dz;
    yc = yc1 + dycdz*dz;

    break;
  }
}

bool InTarget(Double_t xvtx, Double_t yvtx, Double_t zvtx, Double_t R)
{
  Double_t xc, yc;
  CellCenter(zvtx, xc, yc);
  Double_t dx = xvtx-xc;
  Double_t dy = yvtx-yc;
  Double_t r = sqrt(dx*dx + dy*dy);

  return( r <= R );
}

void BinLogX(TH1*h)
{
   TAxis *axis = h->GetXaxis();
   int bins = axis->GetNbins();

   Axis_t from = axis->GetXmin();
   Axis_t to = axis->GetXmax();
   Axis_t width = (to - from) / bins;
   Axis_t *new_bins = new Axis_t[bins + 1];

   for (int i = 0; i <= bins; i++) {
     new_bins[i] = pow(10, from + i * width);
   }
   axis->Set(bins, new_bins);
   delete new_bins;
}

// Fusion sort

void fusion(Double_t* tab, Int_t beg1 , Int_t end1, Int_t end2)
{
   Double_t* tab2 = new Double_t[end1-beg1+1];
   Int_t beg2 = end1+1;
   Int_t loop1 = beg1;
   Int_t loop2 = beg2;

   for(int i=beg1; i<=end1; i++)
      tab2[i-beg1] = tab[i];

   for(int i=beg1; i<=end2; i++)
   {
      if(loop1==beg2)
         break;
      else if(loop2==(end2+1))
      {
         tab[i] = tab2[loop1-beg1];
         loop1++;
      }
      else if(tab2[loop1-beg1]<tab[loop2])
      {
         tab[i] = tab2[loop1-beg1];
         loop1++;
      }
      else
      {
         tab[i] = tab[loop2];
         loop2++;
      }
   }

   delete tab2;
}

void fusionSort2(Double_t* tab, Int_t begin, Int_t end)
{
   if(begin!=end)
   {
      Int_t mid = (begin+end)/2;
      fusionSort2(tab, begin, mid);
      fusionSort2(tab, mid+1, end);
      fusion(tab, begin, mid, end);
   }
}

void fusionSort(Double_t* tab, Int_t len)
{
   if(len>0)
      fusionSort2(tab, 0, len-1);
}

void load_rich_mat(string prich, string prich_err)
{

  pi_sigma_uni[0][0] = 1;
  k_sigma_uni[1][1] = 1;
  p_sigma_uni[2][2] = 1;
  pi_vect[0][0] = 1;
  k_vect[1][0] = 1;
  p_vect[2][0] = 1;

  for(int i=0; i<2; i++)
  {
    for(int j=0; j<10; j++)
    {
      rich_mat_p[i][j].ResizeTo(3,3);
      rich_mat_m[i][j].ResizeTo(3,3);
      inv_rich_p[i][j].ResizeTo(3,3);
      inv_rich_m[i][j].ResizeTo(3,3);
    }
  }

  for(int i=0; i<3; i++)
  {
    for(int j=0; j<3; j++)
    {
      err_rich_p[i][j].ResizeTo(3,3);
      err_rich_m[i][j].ResizeTo(3,3);
    }
  }

  ifstream matRICH(prich);

  for(int i=0; i<24; i++)
  {
    if(i < 4)
    {
      for(int j=0; j<20; j++)
        matRICH >> dummy;
    }
    else
    {
      matRICH >> mat_bin[0][i] >> mat_bin[1][i];
      for(int j=0; j<9; j++)
      {
        matRICH >> rich_mat_p[i%2][(i-4)/2][(int)j/3][j%3];
      }
      for(int j=0; j<9; j++)
      {
        matRICH >> rich_mat_m[i%2][(i-4)/2][(int)j/3][j%3];
      }
    }
  }

  matRICH.close();

  for(int i=0; i<10; i++)
  {
    inv_rich_p[0][i] = rich_mat_p[0][i].InvertFast();
    inv_rich_p[1][i] = rich_mat_p[1][i].InvertFast();
    inv_rich_m[0][i] = rich_mat_m[0][i].InvertFast();
    inv_rich_m[1][i] = rich_mat_m[1][i].InvertFast();
  }

  // Errors YODO

  ifstream errRICH(prich_err);


  for(int loop=0; loop<24; loop++)
  {
    if(loop < 4)
    {
      for(int j=0; j<38; j++)
        errRICH >> dummy;
    }
    else
    {
      errRICH >> err_bin[0][loop-4] >> err_bin[1][loop-4];

      errRICH >> err_rich_p[0][0][0][0];
      errRICH >> err_rich_p[0][1][0][1];
      errRICH >> err_rich_p[0][2][0][2];//3
      errRICH >> err_rich_p[0][0][1][0];
      errRICH >> err_rich_p[0][0][2][0];
      errRICH >> err_rich_p[1][0][2][0];//6
      errRICH >> err_rich_p[1][0][1][0];
      errRICH >> err_rich_p[1][1][1][1];
      errRICH >> err_rich_p[1][2][1][2];//9
      errRICH >> err_rich_p[0][1][1][1];
      errRICH >> err_rich_p[0][1][2][1];
      errRICH >> err_rich_p[1][1][2][1];//12
      errRICH >> err_rich_p[2][0][2][0];
      errRICH >> err_rich_p[2][1][2][1];
      errRICH >> err_rich_p[2][2][2][2];//15
      errRICH >> err_rich_p[0][2][1][2];
      errRICH >> err_rich_p[0][2][2][2];
      errRICH >> err_rich_p[1][2][2][2];//18

      err_rich_p[0][0][0][0] = pow(err_rich_p[0][0][0][0],2);
      err_rich_p[0][1][0][1] = pow(err_rich_p[0][1][0][1],2);
      err_rich_p[0][2][0][2] = pow(err_rich_p[0][2][0][2],2);
      err_rich_p[1][0][1][0] = pow(err_rich_p[1][0][1][0],2);
      err_rich_p[1][1][1][1] = pow(err_rich_p[1][1][1][1],2);
      err_rich_p[1][2][1][2] = pow(err_rich_p[1][2][1][2],2);
      err_rich_p[2][0][2][0] = pow(err_rich_p[2][0][2][0],2);
      err_rich_p[2][1][2][1] = pow(err_rich_p[2][1][2][1],2);
      err_rich_p[2][2][2][2] = pow(err_rich_p[2][2][2][2],2);

      err_rich_p[0][0][1][0] = err_rich_p[1][0][0][0];
      err_rich_p[0][0][2][0] = err_rich_p[2][0][0][0];
      err_rich_p[1][0][2][0] = err_rich_p[2][0][1][0];

      err_rich_p[0][1][1][1] = err_rich_p[1][1][0][1];
      err_rich_p[0][1][2][1] = err_rich_p[2][1][0][1];
      err_rich_p[1][1][2][1] = err_rich_p[2][1][1][1];

      err_rich_p[0][2][1][2] = err_rich_p[1][2][0][2];
      err_rich_p[0][2][2][2] = err_rich_p[2][2][0][2];
      err_rich_p[1][2][2][2] = err_rich_p[2][2][1][2];

      errRICH >> err_rich_m[0][0][0][0];
      errRICH >> err_rich_m[0][1][0][1];
      errRICH >> err_rich_m[0][2][0][2];//3
      errRICH >> err_rich_m[0][0][1][0];
      errRICH >> err_rich_m[0][0][2][0];
      errRICH >> err_rich_m[1][0][2][0];//6
      errRICH >> err_rich_m[1][0][1][0];
      errRICH >> err_rich_m[1][1][1][1];
      errRICH >> err_rich_m[1][2][1][2];//9
      errRICH >> err_rich_m[0][1][1][1];
      errRICH >> err_rich_m[0][1][2][1];
      errRICH >> err_rich_m[1][1][2][1];//12
      errRICH >> err_rich_m[2][0][2][0];
      errRICH >> err_rich_m[2][1][2][1];
      errRICH >> err_rich_m[2][2][2][2];//15
      errRICH >> err_rich_m[0][2][1][2];
      errRICH >> err_rich_m[0][2][2][2];
      errRICH >> err_rich_m[1][2][2][2];//18

      err_rich_m[0][0][0][0] = pow(err_rich_m[0][0][0][0],2);
      err_rich_m[0][1][0][1] = pow(err_rich_m[0][1][0][1],2);
      err_rich_m[0][2][0][2] = pow(err_rich_m[0][2][0][2],2);
      err_rich_m[1][0][1][0] = pow(err_rich_m[1][0][1][0],2);
      err_rich_m[1][1][1][1] = pow(err_rich_m[1][1][1][1],2);
      err_rich_m[1][2][1][2] = pow(err_rich_m[1][2][1][2],2);
      err_rich_m[2][0][2][0] = pow(err_rich_m[2][0][2][0],2);
      err_rich_m[2][1][2][1] = pow(err_rich_m[2][1][2][1],2);
      err_rich_m[2][2][2][2] = pow(err_rich_m[2][2][2][2],2);

      err_rich_m[0][0][1][0] = err_rich_m[1][0][0][0];
      err_rich_m[0][0][2][0] = err_rich_m[2][0][0][0];
      err_rich_m[1][0][2][0] = err_rich_m[2][0][1][0];

      err_rich_m[0][1][1][1] = err_rich_m[1][1][0][1];
      err_rich_m[0][1][2][1] = err_rich_m[2][1][0][1];
      err_rich_m[1][1][2][1] = err_rich_m[2][1][1][1];

      err_rich_m[0][2][1][2] = err_rich_m[1][2][0][2];
      err_rich_m[0][2][2][2] = err_rich_m[2][2][0][2];
      err_rich_m[1][2][2][2] = err_rich_m[2][2][1][2];

#ifdef DEBUG
      for(int i=0; i<9; i++)
      {
        for(int j=0; j<9; j++)
        {
          cout << err_rich_p[i][j] << " ";
        }
        cout << endl;
      }
#endif

      cout << "\n" << endl;

      for(int i=0; i<3; i++)
      {
        cov1_pi[0][i] = 0;
        cov1_pi[1][i] = 0;
        cov1_k[0][i] = 0;
        cov1_k[1][i] = 0;
        cov1_p[0][i] = 0;
        cov1_p[1][i] = 0;
        cov2[0][i] = 0;
        cov2[1][i] = 0;
      }

      for(int i=0; i<3; i++)
      {
        for(int j=0; j<3; j++)
        {
          cov1_pi[0][i] += pow(inv_rich_p[loop%2][(loop-4)/2][i][j]*pi_sigma_uni[j][j],2);
          cov1_pi[1][i] += pow(inv_rich_m[loop%2][(loop-4)/2][i][j]*pi_sigma_uni[j][j],2);
          cov1_k[0][i] += pow(inv_rich_p[loop%2][(loop-4)/2][i][j]*k_sigma_uni[j][j],2);
          cov1_k[1][i] += pow(inv_rich_m[loop%2][(loop-4)/2][i][j]*k_sigma_uni[j][j],2);
          cov1_p[0][i] += pow(inv_rich_p[loop%2][(loop-4)/2][i][j]*p_sigma_uni[j][j],2);
          cov1_p[1][i] += pow(inv_rich_m[loop%2][(loop-4)/2][i][j]*p_sigma_uni[j][j],2);

          for(int k=0; k<3; k++)
          {
            for(int l=0; l<3; l++)
            {
              for(int m=0; m<3; m++)
              {
                    cov2[0][i] += inv_rich_p[loop%2][(loop-4)/2][i][j]
                                   *inv_rich_p[loop%2][(loop-4)/2][i][l]
                                   *inv_rich_p[loop%2][(loop-4)/2][k][i]
                                   *inv_rich_p[loop%2][(loop-4)/2][m][i]
                                   *err_rich_p[j][k][l][m];
#ifdef DEBUG
                                  cout << err_rich_p[j*3+k][l*3+m] << endl;
                                  cout << j*3+k << " " << l*3+m << endl;
#endif
                    cov2[1][i] += inv_rich_m[loop%2][(loop-4)/2][i][j]
                                   *inv_rich_m[loop%2][(loop-4)/2][i][l]
                                   *inv_rich_m[loop%2][(loop-4)/2][k][i]
                                   *inv_rich_m[loop%2][(loop-4)/2][m][i]
                                   *err_rich_m[j][k][l][m];
              }
            }
          }
        }
        pi_unfolding_err_p[loop%2][(loop-4)/2][i] = cov1_pi[0][i] + cov2[0][i];
        pi_unfolding_err_m[loop%2][(loop-4)/2][i] = cov1_pi[1][i] + cov2[1][i];
        k_unfolding_err_p[loop%2][(loop-4)/2][i] = cov1_k[0][i] + cov2[0][i];
        k_unfolding_err_m[loop%2][(loop-4)/2][i] = cov1_k[1][i] + cov2[1][i];
        p_unfolding_err_p[loop%2][(loop-4)/2][i] = cov1_p[0][i] + cov2[0][i];
        p_unfolding_err_m[loop%2][(loop-4)/2][i] = cov1_p[1][i] + cov2[1][i];
      }
    }
  }

  errRICH.close();
}

void create_kin_plots()
{
  fKinematicsRD[0] = new TH1F("Q^{2}", "Q^{2}", 100, 0, 2);
  fKinematicsRD[1] = new TH1F("x_{Bj}", "x_{Bj}", 100, -3, 0);
  fKinematicsRD[2] = new TH1F("y", "y", 100, 0, 1);
  fKinematicsRD[3] = new TH1F("z", "z", 100, 0, 1);
  fKinematicsRD[4] = new TH1F("W", "W", 100, 2, 18);
  fKinematicsRD[5] = new TH1F("#nu", "#nu", 100, 0, 160);
  fKinematicsRD[6] = new TH1F("#Phi_h", "#Phi_h", 100, 0, 3);
  fKinematicsMC[0] = new TH1F("Q^{2} MC", "Q^{2} MC", 100, 0, 2);
  fKinematicsMC[1] = new TH1F("x_{Bj} MC", "x_{Bj} MC", 100, -3, 0);
  fKinematicsMC[2] = new TH1F("y MC", "y MC", 100, 0, 1);
  fKinematicsMC[3] = new TH1F("z MC", "z MC", 100, 0, 1);
  fKinematicsMC[4] = new TH1F("W MC", "W MC", 100, 2, 18);
  fKinematicsMC[5] = new TH1F("#nu MC", "#nu MC", 100, 0, 160);
  fKinematicsMC[6] = new TH1F("#Phi_h MC", "#Phi_h MC", 100, 0, 3);
  BinLogX(fKinematicsRD[0]);
  BinLogX(fKinematicsMC[0]);
  BinLogX(fKinematicsRD[1]);
  BinLogX(fKinematicsMC[1]);
}

void save_kin_plots()
{
  c1.Divide(1,1);
  c2.Divide(1,1);
  c3.Divide(1,1);
  c4.Divide(1,1);
  c5.Divide(1,1);
  c6.Divide(1,1);
  c7.Divide(1,1);
  c1.cd(1);
  fKinematicsRD[0]->Scale(1/fKinematicsRD[0]->GetEntries());
  fKinematicsRD[0]->SetLineColor(kRed);
  fKinematicsRD[0]->SetStats(0);
  fKinematicsRD[0]->SetMinimum(0.);
  fKinematicsRD[0]->SetMaximum(0.05);
  fKinematicsRD[0]->Draw();
  fKinematicsMC[0]->Scale(1/fKinematicsMC[0]->GetEntries());
  fKinematicsMC[0]->SetLineColor(kBlue);
  fKinematicsMC[0]->Draw("SAME");
  gPad->SetLogx();
  c1.Update();
  c2.cd(1);
  fKinematicsRD[1]->Scale(1/fKinematicsRD[1]->GetEntries());
  fKinematicsRD[1]->SetLineColor(kRed);
  fKinematicsRD[1]->SetStats(0);
  fKinematicsRD[1]->SetMinimum(0.);
  fKinematicsRD[1]->SetMaximum(0.035);
  fKinematicsRD[1]->Draw();
  fKinematicsMC[1]->Scale(1/fKinematicsMC[1]->GetEntries());
  fKinematicsMC[1]->SetLineColor(kBlue);
  fKinematicsMC[1]->Draw("SAME");
  gPad->SetLogx();
  c2.Update();
  c3.cd(1);
  fKinematicsRD[2]->Scale(1/fKinematicsRD[2]->Integral(), "width");
  fKinematicsRD[2]->SetLineColor(kRed);
  fKinematicsRD[2]->SetStats(0);
  fKinematicsRD[2]->SetMinimum(0.);
  fKinematicsRD[2]->SetMaximum(5.);
  fKinematicsRD[2]->Draw();
  fKinematicsMC[2]->Scale(1/fKinematicsMC[2]->Integral(), "width");
  fKinematicsMC[2]->SetLineColor(kBlue);
  fKinematicsMC[2]->Draw("SAME");
  c3.Update();
  c4.cd(1);
  fKinematicsRD[3]->Scale(1/fKinematicsRD[3]->Integral(), "width");
  fKinematicsRD[3]->SetLineColor(kRed);
  fKinematicsRD[3]->SetStats(0);
  fKinematicsRD[3]->Draw();
  fKinematicsMC[3]->Scale(1/fKinematicsMC[3]->Integral(), "width");
  fKinematicsMC[3]->SetLineColor(kBlue);
  fKinematicsMC[3]->Draw("SAME");
  c4.Update();
  c5.cd(1);
  fKinematicsRD[4]->Scale(1/fKinematicsRD[4]->Integral(), "width");
  fKinematicsRD[4]->SetLineColor(kRed);
  fKinematicsRD[4]->SetMinimum(0.);
  fKinematicsRD[4]->SetMaximum(0.22);
  fKinematicsRD[4]->Draw();
  fKinematicsMC[4]->Scale(1/fKinematicsMC[4]->Integral(), "width");
  fKinematicsMC[4]->SetLineColor(kBlue);
  fKinematicsMC[4]->Draw("SAME");
  c5.Update();
  c6.cd(1);
  fKinematicsRD[5]->Scale(1/fKinematicsRD[5]->Integral(), "width");
  fKinematicsRD[5]->SetLineColor(kRed);
  fKinematicsRD[5]->SetStats(0);
  fKinematicsRD[5]->SetMinimum(0.);
  fKinematicsRD[5]->SetMaximum(0.035);
  fKinematicsRD[5]->Draw();
  fKinematicsMC[5]->Scale(1/fKinematicsMC[5]->Integral(), "width");
  fKinematicsMC[5]->SetLineColor(kBlue);
  fKinematicsMC[5]->Draw("SAME");
  c6.Update();
  c7.cd(1);
  fKinematicsRD[6]->Scale(1/fKinematicsRD[6]->Integral(), "width");
  fKinematicsRD[6]->SetLineColor(kRed);
  fKinematicsRD[6]->SetStats(0);
  fKinematicsRD[6]->SetMinimum(0.);
  fKinematicsRD[6]->SetMaximum(0.4);
  fKinematicsRD[6]->Draw();
  fKinematicsMC[6]->Scale(1/fKinematicsMC[6]->Integral(), "width");
  fKinematicsMC[6]->SetLineColor(kBlue);
  fKinematicsMC[6]->Draw("SAME");
  c7.Update();

  c1.Print("kinMCRD.pdf(","pdf");
  c2.Print("kinMCRD.pdf","pdf");
  c3.Print("kinMCRD.pdf","pdf");
  c4.Print("kinMCRD.pdf","pdf");
  c5.Print("kinMCRD.pdf","pdf");
  c6.Print("kinMCRD.pdf","pdf");
  c7.Print("kinMCRD.pdf)","pdf");
}

void MCextraction(string pFilelist)
{

  //Kinematics
  Double_t Q2 = 0;
  Double_t xBj = 0;
  Double_t yBj = 0;
  Double_t zBj = 0;
  Double_t zBj_unid = 0;
  Double_t wBj = 0;
  Double_t nu = 0;

  Double_t MCE0 = 0;
  Double_t MCE1 = 0;
  Double_t Q2_MC = 0;
  Double_t xBj_MC = 0;
  Double_t yBj_MC = 0;
  Double_t zBj_MC = 0;
  Double_t zBj_MC_unid = 0;
  Double_t wBj_MC = 0;
  Double_t nu_MC = 0;

  // Target cells
  if(Y2012) InitTargetFile(target_file_2012);
  else if(Y2016) InitTargetFile(target_file_2016);

  //----------------------------------------------------------------------------
  //--------- nu cut prep ------------------------------------------------------
  //----------------------------------------------------------------------------

  for(int i=0; i<12; i++)
  {
    fNu_max[1][i] = sqrt(pow(40,2)+pow(fM_K,2))/fZrange[i+1];
    fNu_min[1][i] = sqrt(pow(MOMENTUM,2)+pow(fM_K,2))/fZrange[i];

    fNu_max[2][i] = sqrt(pow(40,2)+pow(fM_p,2))/fZrange[i+1];
    fNu_min[2][i] = sqrt(pow(MOMENTUM,2)+pow(fM_p,2))/fZrange[i];

    fNu_max[0][i] = sqrt(pow(40,2)+pow(fM_pi,2))/fZrange[i+1];
    fNu_min[0][i] = sqrt(pow(MOMENTUM,2)+pow(fM_pi,2))/fZrange[i];
  }

  // List of files

  ifstream list(pFilelist);
  string filename;

  while(list >> filename)
  {
    TFile *f;

    cout << ".. Processing file " << filename << " .." << endl;
    f = TFile::Open(filename.c_str());

    if(!f) continue;

    TTree* tree = (TTree*) f->Get("DISEvtTree");

    // ---------------------------------------------------------------------------
    // --------- Reading of the TTree --------------------------------------------
    // ---------------------------------------------------------------------------

    //DISEvt
    TBranch *runNo = (TBranch*) tree->FindBranch("runNo");
    TBranch *spillNo = (TBranch*) tree->FindBranch("spillNo");
    TBranch *evtInSpill = (TBranch*) tree->FindBranch("evtInSpill");
    TBranch *trigMask = (TBranch*) tree->FindBranch("trigMask");
    //TBranch *evNo = (TBranch*) tree->FindBranch("evNo");
    TBranch *x = (TBranch*) tree->FindBranch("x");
    TBranch *y = (TBranch*) tree->FindBranch("y");
    TBranch *z = (TBranch*) tree->FindBranch("z");
    TBranch *p0x = (TBranch*) tree->FindBranch("p0x");
    TBranch *p0y = (TBranch*) tree->FindBranch("p0y");
    TBranch *p0z = (TBranch*) tree->FindBranch("p0z");
    TBranch *p1x = (TBranch*) tree->FindBranch("p1x");
    TBranch *p1y = (TBranch*) tree->FindBranch("p1y");
    TBranch *p1z = (TBranch*) tree->FindBranch("p1z");
    TBranch *E_beam = (TBranch*) tree->FindBranch("E_beam");
    TBranch *E_mu_prim = (TBranch*) tree->FindBranch("E_mu_prim");
    TBranch *XX0 = (TBranch*) tree->FindBranch("XX0");
    TBranch *HM04x = (TBranch*) tree->FindBranch("HM04x");
    TBranch *HM04y = (TBranch*) tree->FindBranch("HM04y");
    TBranch *HM05x = (TBranch*) tree->FindBranch("HM05x");
    TBranch *HM05y = (TBranch*) tree->FindBranch("HM05y");
    //TBranch *HO03x = (TBranch*) tree->FindBranch("HO03x");
    //TBranch *HO03y = (TBranch*) tree->FindBranch("HO03y");
    //TBranch *HO04x = (TBranch*) tree->FindBranch("HO04x");
    //TBranch *HO04y = (TBranch*) tree->FindBranch("HO04y");
    //TBranch *saved = (TBranch*) tree->FindBranch("saved");
    TBranch *cellsCrossed = (TBranch*) tree->FindBranch("cellsCrossed");
    TBranch *backPropFlag = (TBranch*) tree->FindBranch("backPropFlag");

    //Hadrons
    TBranch *p = (TBranch*) tree->FindBranch("Hadrons.P");
    TBranch *th = (TBranch*) tree->FindBranch("Hadrons.th");
    TBranch *ph = (TBranch*) tree->FindBranch("Hadrons.ph");
    TBranch *hXX0 = (TBranch*) tree->FindBranch("Hadrons.XX0");
    TBranch *inHCALacc = (TBranch*) tree->FindBranch("Hadrons.inHCALacc");
    TBranch *HCAL = (TBranch*) tree->FindBranch("Hadrons.HCAL");
    TBranch *charge = (TBranch*) tree->FindBranch("Hadrons.charge");
    TBranch *thRICH = (TBranch*) tree->FindBranch("Hadrons.thRICH");
    //TBranch *LH = (TBranch*) tree->FindBranch("Hadrons.LH");
    TBranch *MCpid = (TBranch*) tree->FindBranch("Hadrons.MCpid");
    //TBranch *MM01x = (TBranch*) tree->FindBranch("Hadrons.MM01x");
    //TBranch *MM01y = (TBranch*) tree->FindBranch("Hadrons.MM01y");
    //TBranch *MM02x = (TBranch*) tree->FindBranch("Hadrons.MM02x");
    //TBranch *MM02y = (TBranch*) tree->FindBranch("Hadrons.MM02y");
    //TBranch *MM03x = (TBranch*) tree->FindBranch("Hadrons.MM03x");
    //TBranch *MM03y = (TBranch*) tree->FindBranch("Hadrons.MM03y");
    //TBranch *Z2Ax = (TBranch*) tree->FindBranch("Hadrons.Z2Ax");
    //TBranch *Z2Ay = (TBranch*) tree->FindBranch("Hadrons.Z2Ay");
    //TBranch *Z2Bx = (TBranch*) tree->FindBranch("Hadrons.Z2Bx");
    //TBranch *Z2By = (TBranch*) tree->FindBranch("Hadrons.Z2By");
    TBranch *RICHx = (TBranch*) tree->FindBranch("Hadrons.RICHx");
    TBranch *RICHy = (TBranch*) tree->FindBranch("Hadrons.RICHy");

    //DISMCEvt
    TBranch *MC_vx = (TBranch*) tree->FindBranch("MC_vx");
    TBranch *MC_vy = (TBranch*) tree->FindBranch("MC_vy");
    TBranch *MC_vz = (TBranch*) tree->FindBranch("MC_vz");
    TBranch *MC_p0x = (TBranch*) tree->FindBranch("MC_p0x");
    TBranch *MC_p0y = (TBranch*) tree->FindBranch("MC_p0y");
    TBranch *MC_p0z = (TBranch*) tree->FindBranch("MC_p0z");
    TBranch *MC_p1x = (TBranch*) tree->FindBranch("MC_p1x");
    TBranch *MC_p1y = (TBranch*) tree->FindBranch("MC_p1y");
    TBranch *MC_p1z = (TBranch*) tree->FindBranch("MC_p1z");
    TBranch *irad = (TBranch*) tree->FindBranch("irad");
    TBranch *MC_nuTr = (TBranch*) tree->FindBranch("MC_nuTr");
    TBranch *MC_Q2Tr = (TBranch*) tree->FindBranch("MC_Q2Tr");
    TBranch *MC_w = (TBranch*) tree->FindBranch("MC_w");
    TBranch *recons = (TBranch*) tree->FindBranch("recons");
    TBranch *MC_yTr = (TBranch*) tree->FindBranch("MC_yTr");
    TBranch *MC_xTr = (TBranch*) tree->FindBranch("MC_xTr");

    //MCHadrons
    TBranch *MC_p = (TBranch*) tree->FindBranch("MCHadrons.P");
    TBranch *MC_th = (TBranch*) tree->FindBranch("MCHadrons.th");
    TBranch *MC_ph = (TBranch*) tree->FindBranch("MCHadrons.ph");
    TBranch *MC_charge = (TBranch*) tree->FindBranch("MCHadrons.charge");
    TBranch *MC_pid = (TBranch*) tree->FindBranch("MCHadrons.pid");
    TBranch *MC_recons = (TBranch*) tree->FindBranch("MCHadrons.recons");
    TBranch *MC_recHadIdx = (TBranch*) tree->FindBranch("MCHadrons.recHadIdx");

    map<int,int> idMCrec[2][4];
    map<int,Double_t> pMCrec[2][4];
    map<int,Double_t> zMCrec[2][4];
    map<int,int> prevMCrec[2][4];
    map<int,int> hidMCrec[2][4];

    // Loopy loop over the events
    int N = (int) tree->GetEntries();

    for (int ip = 0; ip < N; ip++)
    {

      printProgress(ip,N);

      //DISEvt
      runNo->GetEntry(ip);
      spillNo->GetEntry(ip);
      evtInSpill->GetEntry(ip);
      trigMask->GetEntry(ip);
      //evNo->GetEntry(ip);
      x->GetEntry(ip);
      y->GetEntry(ip);
      z->GetEntry(ip);
      p0x->GetEntry(ip);
      p0y->GetEntry(ip);
      p0z->GetEntry(ip);
      p1x->GetEntry(ip);
      p1y->GetEntry(ip);
      p1z->GetEntry(ip);
      E_beam->GetEntry(ip);
      E_mu_prim->GetEntry(ip);
      XX0->GetEntry(ip);
      HM04x->GetEntry(ip);
      HM04y->GetEntry(ip);
      HM05x->GetEntry(ip);
      HM05y->GetEntry(ip);
      //HO03x->GetEntry(ip);
      //HO03y->GetEntry(ip);
      //HO04x->GetEntry(ip);
      //HO04y->GetEntry(ip);
      //saved->GetEntry(ip);
      cellsCrossed->GetEntry(ip);
      backPropFlag->GetEntry(ip);

      //Hadrons
      p->GetEntry(ip);
      th->GetEntry(ip);
      ph->GetEntry(ip);
      hXX0->GetEntry(ip);
      inHCALacc->GetEntry(ip);
      HCAL->GetEntry(ip);
      charge->GetEntry(ip);
      thRICH->GetEntry(ip);
      //LH->GetEntry(ip);
      MCpid->GetEntry(ip);
      //MM01x->GetEntry(ip);
      //MM01y->GetEntry(ip);
      //MM02x->GetEntry(ip);
      //MM02y->GetEntry(ip);
      //MM03x->GetEntry(ip);
      //MM03y->GetEntry(ip);
      //Z2Ax->GetEntry(ip);
      //Z2Ay->GetEntry(ip);
      //Z2Bx->GetEntry(ip);
      //Z2By->GetEntry(ip);
      RICHx->GetEntry(ip);
      RICHy->GetEntry(ip);

      //DISMCEvt
      MC_vx->GetEntry(ip);
      MC_vy->GetEntry(ip);
      MC_vz->GetEntry(ip);
      MC_p0x->GetEntry(ip);
      MC_p0y->GetEntry(ip);
      MC_p0z->GetEntry(ip);
      MC_p1x->GetEntry(ip);
      MC_p1y->GetEntry(ip);
      MC_p1z->GetEntry(ip);
      irad->GetEntry(ip);
      MC_nuTr->GetEntry(ip);
      MC_Q2Tr->GetEntry(ip);
      MC_w->GetEntry(ip);
      recons->GetEntry(ip);
      MC_yTr->GetEntry(ip);
      MC_xTr->GetEntry(ip);

      //MCHadrons
      MC_p->GetEntry(ip);
      MC_th->GetEntry(ip);
      MC_ph->GetEntry(ip);
      MC_charge->GetEntry(ip);
      MC_pid->GetEntry(ip);
      MC_recons->GetEntry(ip);
      MC_recHadIdx->GetEntry(ip);

      // -------------------------------------------------------------------------
      // --------- Calculation ---------------------------------------------------
      // -------------------------------------------------------------------------

      Double_t zlab = z->GetLeaf("z")->GetValue();
      int zlabbin;

      if(Y2012)
      {
        if(!(-311.2<zlab && zlab<-71.2)) continue;

        if(-311.2<=zlab && zlab<-301.2) zlabbin = 0;
        else if(-301.2<=zlab && zlab<-291.2) zlabbin = 1;
        else if(-291.2<=zlab && zlab<-281.2) zlabbin = 2;
        else if(-281.2<=zlab && zlab<-271.2) zlabbin = 3;
        else if(-271.2<=zlab && zlab<-261.2) zlabbin = 4;
        else if(-261.2<=zlab && zlab<-251.2) zlabbin = 5;
        else if(-251.2<=zlab && zlab<-241.2) zlabbin = 6;
        else if(-241.2<=zlab && zlab<-231.2) zlabbin = 7;
        else if(-231.2<=zlab && zlab<-221.2) zlabbin = 8;
        else if(-221.2<=zlab && zlab<-211.2) zlabbin = 9;
        else if(-211.2<=zlab && zlab<-201.2) zlabbin = 10;
        else if(-201.2<=zlab && zlab<-191.2) zlabbin = 11;
        else if(-191.2<=zlab && zlab<-181.2) zlabbin = 12;
        else if(-181.2<=zlab && zlab<-171.2) zlabbin = 13;
        else if(-171.2<=zlab && zlab<-161.2) zlabbin = 14;
        else if(-161.2<=zlab && zlab<-151.2) zlabbin = 15;
        else if(-151.2<=zlab && zlab<-141.2) zlabbin = 16;
        else if(-141.2<=zlab && zlab<-131.2) zlabbin = 17;
        else if(-131.2<=zlab && zlab<-121.2) zlabbin = 18;
        else if(-121.2<=zlab && zlab<-111.2) zlabbin = 19;
        else if(-111.2<=zlab && zlab<-101.2) zlabbin = 20;
        else if(-101.2<=zlab && zlab<-91.2) zlabbin = 21;
        else if(-91.2<=zlab && zlab<-81.2) zlabbin = 22;
        else zlabbin = 23;
      }

      //--------------------------------------------------------------------------
      //--------- Kinematics -----------------------------------------------------
      //--------------------------------------------------------------------------

      //Data
      Q2 = 2.*( E_beam->GetLeaf("E_beam")->GetValue()*E_mu_prim->GetLeaf("E_mu_prim")->GetValue()
           - p0x->GetLeaf("p0x")->GetValue()*p1x->GetLeaf("p1x")->GetValue()
           - p0y->GetLeaf("p0y")->GetValue()*p1y->GetLeaf("p1y")->GetValue()
           - p0z->GetLeaf("p0z")->GetValue()*p1z->GetLeaf("p1z")->GetValue()
           - pow(fM_mu,2));

      nu = E_beam->GetLeaf("E_beam")->GetValue() - E_mu_prim->GetLeaf("E_mu_prim")->GetValue();

      if(E_beam->GetLeaf("E_beam")->GetValue() != 0)
        yBj = nu/E_beam->GetLeaf("E_beam")->GetValue();
      else
        yBj = 0;

      if(nu != 0)
      {
        xBj = Q2/(2*fM_p*nu);
      }
      else
      {
        xBj = 0;
      }

      if(xBj != 0)
        wBj = pow(fM_p,2) + Q2*(1-xBj)/xBj;
      else
        wBj = 0;

      int trig= trigMask->GetLeaf("trigMask")->GetValue();


      //2006 ---

      if(Y2006)
      {
        if ((trig&256) && HM05x->GetLeaf("HM05x")->GetValue()<(HM05y->GetLeaf("HM05y")->GetValue()>0 ? 14.55-0.15 : 22.02864-0.12864) )
        {
          trig -= 256;
        }
      }

      //2006 ---

      //MC
      MCE0 = sqrt(pow(fM_mu,2)
                  +MC_p0x->GetLeaf("MC_p0x")->GetValue()*MC_p0x->GetLeaf("MC_p0x")->GetValue()
                  +MC_p0y->GetLeaf("MC_p0y")->GetValue()*MC_p0y->GetLeaf("MC_p0y")->GetValue()
                  +MC_p0z->GetLeaf("MC_p0z")->GetValue()*MC_p0z->GetLeaf("MC_p0z")->GetValue());

      MCE1 = sqrt(pow(fM_mu,2)
                  +MC_p1x->GetLeaf("MC_p1x")->GetValue()*MC_p1x->GetLeaf("MC_p1x")->GetValue()
                  +MC_p1y->GetLeaf("MC_p1y")->GetValue()*MC_p1y->GetLeaf("MC_p1y")->GetValue()
                  +MC_p1z->GetLeaf("MC_p1z")->GetValue()*MC_p1z->GetLeaf("MC_p1z")->GetValue());

      Q2_MC = 2.*( MCE0*MCE1
           - MC_p0x->GetLeaf("MC_p0x")->GetValue()*MC_p1x->GetLeaf("MC_p1x")->GetValue()
           - MC_p0y->GetLeaf("MC_p0y")->GetValue()*MC_p1y->GetLeaf("MC_p1y")->GetValue()
           - MC_p0z->GetLeaf("MC_p0z")->GetValue()*MC_p1z->GetLeaf("MC_p1z")->GetValue()
           - pow(fM_mu,2));

      nu_MC = MCE0 - MCE1;

      if(MCE0 != 0)
        yBj_MC = nu_MC/MCE0;
      else
        yBj_MC = 0;

      if(nu_MC != 0)
      {
        xBj_MC = Q2_MC/(2*fM_p*nu_MC);
      }
      else
      {
        xBj_MC = 0;
      }

      if(xBj_MC != 0)
        wBj_MC = pow(fM_p,2) + Q2_MC*(1-xBj_MC)/xBj_MC;
      else
        wBj_MC = 0;

      //--------------------------------------------------------------------------
      //--------- Target ---------------------------------------------------------
      //--------------------------------------------------------------------------

      //2006 ---

      //MC target position new
      static const Double_t dz = 2;

      static const Double_t mcxU = -0.085;
      static const Double_t mcyU = 0.33;
      static const Double_t mczU_1 = -65+dz+4;
      static const Double_t mcxD = -0.085;
      static const Double_t mcyD = 0.33;
      static const Double_t mczD_2 = 65+dz;

      Double_t mcR    = 1.4;

      //target position data 2006
      static const Double_t xU = -0.1;
      static const Double_t yU = 0.33;
      static const Double_t zU_1 = -65+dz+4;
      static const Double_t xD = -0.07;
      static const Double_t yD = 0.33;
      static const Double_t zD_2 =  65+dz;

      Double_t R    = 1.4;
      Double_t yCUT = 1.4;

      Double_t mcxC = (mcxD-mcxU) * (mczU_1-z->GetLeaf("z")->GetValue()) / (mczU_1-mczD_2) + mcxU;
      Double_t mcyC = (mcyD-mcyU) * (mczU_1-z->GetLeaf("z")->GetValue()) / (mczU_1-mczD_2) + mcyU;
      Double_t mcr = sqrt( (x->GetLeaf("x")->GetValue()-mcxC)*(x->GetLeaf("x")->GetValue()-mcxC)
                    + (y->GetLeaf("y")->GetValue()-mcyC)*(y->GetLeaf("y")->GetValue()-mcyC) );
      Double_t xC = (xD-xU) * (zU_1-z->GetLeaf("z")->GetValue()) / (zU_1-zD_2) + xU;
      Double_t yC = (yD-yU) * (zU_1-z->GetLeaf("z")->GetValue()) / (zU_1-zD_2) + yU;
      Double_t r = sqrt( (x->GetLeaf("x")->GetValue()-xC)*(x->GetLeaf("x")->GetValue()-xC)
                    + (y->GetLeaf("y")->GetValue()-yC)*(y->GetLeaf("y")->GetValue()-yC) );

      Double_t MC_mcxC = (mcxD-mcxU) * (mczU_1-MC_vz->GetLeaf("MC_vz")->GetValue()) / (mczU_1-mczD_2) + mcxU;
      Double_t MC_mcyC = (mcyD-mcyU) * (mczU_1-MC_vz->GetLeaf("MC_vz")->GetValue()) / (mczU_1-mczD_2) + mcyU;
      Double_t MC_mcr = sqrt( (MC_vx->GetLeaf("MC_vx")->GetValue()-MC_mcxC)*(MC_vx->GetLeaf("MC_vx")->GetValue()-MC_mcxC)
                    + (MC_vy->GetLeaf("MC_vy")->GetValue()-MC_mcyC)*(MC_vy->GetLeaf("MC_vy")->GetValue()-MC_mcyC) );
      Double_t MC_xC = (xD-xU) * (zU_1-MC_vz->GetLeaf("MC_vz")->GetValue()) / (zU_1-zD_2) + xU;
      Double_t MC_yC = (yD-yU) * (zU_1-MC_vz->GetLeaf("MC_vz")->GetValue()) / (zU_1-zD_2) + yU;
      Double_t MC_r = sqrt( (MC_vx->GetLeaf("MC_vx")->GetValue()-MC_xC)*(MC_vx->GetLeaf("MC_vx")->GetValue()-MC_xC)
                    + (MC_vy->GetLeaf("MC_vy")->GetValue()-MC_yC)*(MC_vy->GetLeaf("MC_vy")->GetValue()-MC_yC) );

      //2006 ---


      // -----------------------------------------------------------------------
      // --------- DIS Selection -----------------------------------------------
      // -----------------------------------------------------------------------

      Double_t mxc, myc;

      // -----------------------------------------------------------------------
      //  Data -----------------------------------------------------------------
      // -----------------------------------------------------------------------

      fAllDISflag = 0;
      int DIS_rec[3][12];

      // Best Primary Vertex
      fBP++;

      // Reconstructed muon
      if((0<E_beam->GetLeaf("E_beam")->GetValue()))
      {
        fRmu++;

        if(Y2012)
        {
          CellCenter(z->GetLeaf("z")->GetValue(), mxc, myc);
        }

        //BMS (reconstructed beam track)
        if(true/*(backPropFlag->GetLeaf("backPropFlag")->GetValue())*/) //not used in acceptance
        {
          fBMS++;

          // Energy of the muon beam
          if((140<E_beam->GetLeaf("E_beam")->GetValue() && E_beam->GetLeaf("E_beam")->GetValue()<180))
          {
            fBEC++;

            //2006 ---
            if(Y2006)
            {
              // Z coordinate within target regions
              if(((-56<z->GetLeaf("z")->GetValue() && z->GetLeaf("z")->GetValue()<-35)
                    ||(-20<z->GetLeaf("z")->GetValue() && z->GetLeaf("z")->GetValue()<31)
                    ||(43<z->GetLeaf("z")->GetValue() && z->GetLeaf("z")->GetValue()<66)))
              {
                if((mcr < mcR &&  (y->GetLeaf("y")->GetValue()-mcyC)<yCUT
                     && r < R
                     &&  (y->GetLeaf("y")->GetValue()-yC)<yCUT
                     && ((z->GetLeaf("z")->GetValue()>(-65+2+7) && z->GetLeaf("z")->GetValue()<(-35+2-2))
                          ||(z->GetLeaf("z")->GetValue() > (-30+2+8) && z->GetLeaf("z")->GetValue() < (30+2-1))
                          ||(z->GetLeaf("z")->GetValue() > (35+2+6) && z->GetLeaf("z")->GetValue() < (65+2-1)))))
                {
                  fTarg++;

                  // Cells crossing
                  if(true/*(cellsCrossed->GetLeaf("cellsCrossed")->GetValue())*/)
                  {
                    //fCell++;

                    // IM/O triggers
                    if((trig&8 || trig&256))
                    {
                      fTrig++;

                      // Q2 cut
                      if((Q2>1))
                      {
                        fQ2test++;

                        // y cut
                        if((0.1<yBj && yBj<0.7))
                        {
                          fYBjtest++;

                          // W cut
                          if((5<sqrt(wBj) && sqrt(wBj)<17))
                          {
                            fWBjtest++;

                            // x cut
                            if((0.004<xBj && xBj<0.4))
                            {
                              fXBjtest++;
                              fAllDISflag = 1;
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }

            }
            //2006 ---

            //2012 ---
            else if(Y2012)
            {
              if(InTarget(x->GetLeaf("x")->GetValue(),y->GetLeaf("y")->GetValue(),z->GetLeaf("z")->GetValue(),fRcutval[zlabbin]))
              {
                fTarg++;

                // Cells crossing
                if(true/*(cellsCrossed->GetLeaf("cellsCrossed")->GetValue())*/)
                {
                  //fCell++;

                  if((trig&2 || trig&4 || trig&8))
                  {
                    fTrig++;

                    // Q2 cut
                    if((Q2>1))
                    {
                      fQ2test++;

                      // y cut
                      if((0.1<yBj && yBj<0.7))
                      {
                        fYBjtest++;

                        // W cut
                        if((5<sqrt(wBj) && sqrt(wBj)<17))
                        {
                          fWBjtest++;
                          if((0.004<xBj && xBj<0.4))
                          {
                            fXBjtest++;
                            fAllDISflag = 1;
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
            //2012 ---

            //2016 ---
            else if(Y2016)
            {
              if(InTarget(x->GetLeaf("x")->GetValue(),y->GetLeaf("y")->GetValue(),z->GetLeaf("z")->GetValue(),1.5))
              {
                fTarg++;

                // Cells crossing
                if(true/*(cellsCrossed->GetLeaf("cellsCrossed")->GetValue())*/)
                {
                  fCell++;

                  if((trig&1 || trig&2 || trig&3 || trig&9))
                  {
                    fTrig++;

                    // Q2 cut
                    if((Q2>1))
                    {
                      fQ2test++;

                      // y cut
                      if((0<yBj && yBj<1))
                      {
                        fYBjtest++;

                        // W cut
                        if((5<sqrt(wBj) && sqrt(wBj)<17))
                        {
                          fWBjtest++;
                          if((0<xBj && xBj<1))
                          {
                            fXBjtest++;
                            fAllDISflag = 1;
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
            //2016 ---
          }
        }
      }


      // -----------------------------------------------------------------------
      // MC --------------------------------------------------------------------
      // -----------------------------------------------------------------------

      fAllDISflag_MC = 0;
      int DIS_MC[3][12];

      // Best Primary Vertex

      if((0<MCE0))
      {

        if(Y2012)
        {
          CellCenter(z->GetLeaf("z")->GetValue(), mxc, myc);
        }

        // Energy of the muon beam
        if((140<=MCE0 && MCE0<=180))
        {
          //2006 ---
          if(Y2006)
          {
            // Z coordinate within target regions
            if((((-56<MC_vz->GetLeaf("MC_vz")->GetValue() && MC_vz->GetLeaf("MC_vz")->GetValue()<-35)
                  ||(-20<MC_vz->GetLeaf("MC_vz")->GetValue() && MC_vz->GetLeaf("MC_vz")->GetValue()<31)
                  ||(43<MC_vz->GetLeaf("MC_vz")->GetValue() && MC_vz->GetLeaf("MC_vz")->GetValue()<66))))
            {
              if((MC_mcr < mcR &&  (MC_vy->GetLeaf("MC_vy")->GetValue()-MC_mcyC)<yCUT
                   && MC_r < R
                   &&  (MC_vy->GetLeaf("MC_vy")->GetValue()-yC)<yCUT
                   && ((MC_vz->GetLeaf("MC_vz")->GetValue()>(-65+2+7) && MC_vz->GetLeaf("MC_vz")->GetValue()<(-35+2-2))
                        ||(MC_vz->GetLeaf("MC_vz")->GetValue() > (-30+2+8) && MC_vz->GetLeaf("MC_vz")->GetValue() < (30+2-1))
                        ||(MC_vz->GetLeaf("MC_vz")->GetValue() > (35+2+6) && MC_vz->GetLeaf("MC_vz")->GetValue() < (65+2-1)))))
              {
                // Q2 cut
                if((Q2_MC>1))
                {
                  // y cut
                  if((0.1<yBj_MC && yBj_MC<0.7))
                  {
                    // W cut
                    if((5<sqrt(wBj_MC) && sqrt(wBj_MC)<17))
                    {
                      // x cut
                      if((0.004<xBj_MC && xBj_MC<0.4))
                      {
                        fAllDISflag_MC = 1;
                      }
                    }
                  }
                }
              }
            }
          }
          //2006 ---

          //2012 ---
          else if(Y2012)
          {
            if(InTarget(x->GetLeaf("x")->GetValue(),y->GetLeaf("y")->GetValue(),z->GetLeaf("z")->GetValue(),fRcutval[zlabbin]))
            {
              // Q2 cut
              if((Q2_MC>1))
              {
                // y cut
                if((0.1<yBj_MC && yBj_MC<0.7))
                {
                  // W cut
                  if((5<sqrt(wBj_MC) && sqrt(wBj_MC)<17))
                  {
                    // x cut
                    if((0.004<xBj_MC && xBj_MC<0.4))
                    {
                      fAllDISflag_MC = 1;
                    }
                  }
                }
              }
            }
          }
          //2012 ---

          //2016 ---
          else if(Y2016)
          {
            if(InTarget(MC_vx->GetLeaf("MC_vx")->GetValue(),MC_vy->GetLeaf("MC_vy")->GetValue(),MC_vz->GetLeaf("MC_vz")->GetValue(),1.5))
            {
              // Q2 cut
              if((Q2_MC>1))
              {
                // y cut
                if((0<yBj_MC && yBj_MC<1))
                {
                  // W cut
                  if((5<sqrt(wBj_MC) && sqrt(wBj_MC)<17))
                  {
                    // x cut
                    if((0<xBj_MC && xBj_MC<1))
                    {
                      fAllDISflag_MC = 1;
                    }
                  }
                }
              }
            }
          }
          //2016 ---
        }
      }

      if(fAllDISflag)
      {
        fQ2kinMC.push_back(Q2);
        fXBjkinMC.push_back(xBj);
        fYBjkinMC.push_back(yBj);
        fWBjkinMC.push_back(sqrt(wBj));
        fNukinMC.push_back(nu);
        fXMC.push_back(x->GetLeaf("x")->GetValue());
        fYMC.push_back(y->GetLeaf("y")->GetValue());
      }

      // -----------------------------------------------------------------------
      // -----------------------------------------------------------------------
      // --------- DIS event calculation ---------------------------------------
      // -----------------------------------------------------------------------
      // -----------------------------------------------------------------------

      for(int i=0; i<12; i++)
      {
        DIS_rec[0][i] = 0;
        DIS_rec[1][i] = 0;
        DIS_rec[2][i] = 0;
        DIS_MC[0][i] = 0;
        DIS_MC[1][i] = 0;
        DIS_MC[2][i] = 0;
      }

      // -----------------------------------------------------------------------
      // MC --------------------------------------------------------------------
      // -----------------------------------------------------------------------

      // x Binning

      if(0.004<xBj_MC && xBj_MC<0.01) xbin_MC = 0;
      else if(0.01<=xBj_MC && xBj_MC<0.02) xbin_MC = 1;
      else if(0.02<=xBj_MC && xBj_MC<0.03) xbin_MC = 2;
      else if(0.03<=xBj_MC && xBj_MC<0.04) xbin_MC = 3;
      else if(0.04<=xBj_MC && xBj_MC<0.06) xbin_MC = 4;
      else if(0.06<=xBj_MC && xBj_MC<0.1) xbin_MC = 5;
      else if(0.1<=xBj_MC && xBj_MC<0.14) xbin_MC = 6;
      else if(0.14<=xBj_MC && xBj_MC<0.18) xbin_MC = 7;
      else xbin_MC = 8;

      if(0.004<xBj && xBj<0.01) xbin = 0;
      else if(0.01<=xBj && xBj<0.02) xbin = 1;
      else if(0.02<=xBj && xBj<0.03) xbin = 2;
      else if(0.03<=xBj && xBj<0.04) xbin = 3;
      else if(0.04<=xBj && xBj<0.06) xbin = 4;
      else if(0.06<=xBj && xBj<0.1) xbin = 5;
      else if(0.1<=xBj && xBj<0.14) xbin = 6;
      else if(0.14<=xBj && xBj<0.18) xbin = 7;
      else xbin = 8;

      // y Binning

      if(0.1<yBj_MC && yBj_MC<0.15) ybin_MC = 0;
      else if(0.15<=yBj_MC && yBj_MC<0.2) ybin_MC = 1;
      else if(0.2<=yBj_MC && yBj_MC<0.3) ybin_MC = 2;
      else if(0.3<=yBj_MC && yBj_MC<0.5) ybin_MC = 3;
      else ybin_MC = 4;

      if(0.1<yBj && yBj<0.15) ybin = 0;
      else if(0.15<=yBj && yBj<0.2) ybin = 1;
      else if(0.2<=yBj && yBj<0.3) ybin = 2;
      else if(0.3<=yBj && yBj<0.5) ybin = 3;
      else ybin = 4;

      if(fAllDISflag_MC)
      {
        // z Binnig

        for(int i=0; i<12; i++)
        {
          fNDIS_evt_MC[0][xbin_MC][ybin_MC][i]++;
          fNDIS_evt_MC[1][xbin_MC][ybin_MC][i]++;
          fNDIS_evt_MC[2][xbin_MC][ybin_MC][i]++;

          fFlag_MC[0][xbin_MC][ybin_MC][i]=0;
          fFlag_MC[1][xbin_MC][ybin_MC][i]=0;
          fFlag_MC[2][xbin_MC][ybin_MC][i]=0;

          DIS_MC[0][i] = 1;
          DIS_MC[1][i] = 1;
          DIS_MC[2][i] = 1;

          // nu cut
          if(!(fNu_min[0][i]<nu_MC && nu_MC<fNu_max[0][i]))
          {
            fFlag_MC[0][xbin_MC][ybin_MC][i]=1;
          }
          if(!(fNu_min[1][i]<nu_MC && nu_MC<fNu_max[1][i]))
          {
            fFlag_MC[1][xbin_MC][ybin_MC][i]=1;
          }
          if(!(fNu_min[2][i]<nu_MC && nu_MC<fNu_max[2][i]))
          {
            fFlag_MC[2][xbin_MC][ybin_MC][i]=1;
          }
          if(fFlag_MC[0][xbin_MC][ybin_MC][i])
          {
            fNDIS_evt_MC[0][xbin_MC][ybin_MC][i]--; DIS_MC[0][i] = 0;
          }
          if(fFlag_MC[1][xbin_MC][ybin_MC][i])
          {
            fNDIS_evt_MC[1][xbin_MC][ybin_MC][i]--; DIS_MC[1][i] = 0;
          }
          if(fFlag_MC[2][xbin_MC][ybin_MC][i])
          {
            fNDIS_evt_MC[2][xbin_MC][ybin_MC][i]--; DIS_MC[2][i] = 0;
          }
        }
      }
      else
      {
        for(int i=0; i<12; i++)
        {
          fFlag_MC[0][xbin_MC][ybin_MC][i]=0;
          fFlag_MC[1][xbin_MC][ybin_MC][i]=0;
          fFlag_MC[2][xbin_MC][ybin_MC][i]=0;

          DIS_MC[0][i] = 1;
          DIS_MC[1][i] = 1;
          DIS_MC[2][i] = 1;

          // nu cut
          if(!(fNu_min[0][i]<nu_MC && nu_MC<fNu_max[0][i]))
          {
            fFlag_MC[0][xbin_MC][ybin_MC][i]=1;
          }
          if(!(fNu_min[1][i]<nu_MC && nu_MC<fNu_max[1][i]))
          {
            fFlag_MC[1][xbin_MC][ybin_MC][i]=1;
          }
          if(!(fNu_min[2][i]<nu_MC && nu_MC<fNu_max[2][i]))
          {
            fFlag_MC[2][xbin_MC][ybin_MC][i]=1;
          }
          if(fFlag_MC[0][xbin_MC][ybin_MC][i])
          {
            DIS_MC[0][i] = 0;
          }
          if(fFlag_MC[1][xbin_MC][ybin_MC][i])
          {
            DIS_MC[1][i] = 0;
          }
          if(fFlag_MC[2][xbin_MC][ybin_MC][i])
          {
            DIS_MC[2][i] = 0;
          }
        }
      }

      // -----------------------------------------------------------------------
      //  Data -----------------------------------------------------------------
      // -----------------------------------------------------------------------

      if(fAllDISflag)
      {
        // z Binning

        for(int i=0; i<12; i++)
        {
          fNDIS_evt[0][xbin][ybin][i]++;
          fNDIS_evt[1][xbin][ybin][i]++;
          fNDIS_evt[2][xbin][ybin][i]++;

          fFlag[0][xbin][ybin][i]=0;
          fFlag[1][xbin][ybin][i]=0;
          fFlag[2][xbin][ybin][i]=0;

          DIS_rec[0][i] = 1;
          DIS_rec[1][i] = 1;
          DIS_rec[2][i] = 1;

          // nu cut
          if(!(fNu_min[0][i]<nu && nu<fNu_max[0][i]))
          {
            fFlag[0][xbin][ybin][i]=1;
          }
          if(!(fNu_min[1][i]<nu && nu<fNu_max[1][i]))
          {
            fFlag[1][xbin][ybin][i]=1;
          }
          if(!(fNu_min[2][i]<nu && nu<fNu_max[2][i]))
          {
            fFlag[2][xbin][ybin][i]=1;
          }
          if(fFlag[0][xbin][ybin][i])
          {
            fNDIS_evt[0][xbin][ybin][i]--; DIS_rec[0][i] = 0;
          }
          if(fFlag[1][xbin][ybin][i])
          {
            fNDIS_evt[1][xbin][ybin][i]--; DIS_rec[1][i] = 0;
          }
          if(fFlag[2][xbin][ybin][i])
          {
            fNDIS_evt[2][xbin][ybin][i]--; DIS_rec[2][i] = 0;
          }
          if(xbin==xbin_MC && ybin==ybin_MC)
          {
            if(DIS_rec[0][i] && DIS_MC[0][i])
            {
              fNDIS_evt_c[0][xbin][ybin][i] += 1;
            }
            if(DIS_rec[1][i] && DIS_MC[1][i])
            {
              fNDIS_evt_c[1][xbin][ybin][i] += 1;
            }
            if(DIS_rec[2][i] && DIS_MC[2][i])
            {
              fNDIS_evt_c[2][xbin][ybin][i] += 1;
            }
          }
        }
      }

      // -----------------------------------------------------------------------
      // -----------------------------------------------------------------------
      // --------- Hadrons Selection -------------------------------------------
      // -----------------------------------------------------------------------
      // -----------------------------------------------------------------------


      // -----------------------------------------------------------------------
      //  MC -------------------------------------------------------------------
      // -----------------------------------------------------------------------

      if(fAllDISflag_MC)
      {

        for(int i=0; i<MC_p->GetLeaf("MCHadrons.P")->GetLen(); i++)
        {

          if(MC_pid->GetLeaf("MCHadrons.pid")->GetValue(i) == 8)//pi+
          {
            fId = 0;
          }
          else if(MC_pid->GetLeaf("MCHadrons.pid")->GetValue(i) == 9)//pi-
          {
            fId = 1;
          }
          else if(MC_pid->GetLeaf("MCHadrons.pid")->GetValue(i) == 11)//K+
          {
            fId = 2;
          }
          else if(MC_pid->GetLeaf("MCHadrons.pid")->GetValue(i) == 12)//K-
          {
            fId = 3;
          }
          else if(MC_pid->GetLeaf("MCHadrons.pid")->GetValue(i) == 14)//p
          {
            fId = 4;
          }
          else if(MC_pid->GetLeaf("MCHadrons.pid")->GetValue(i) == 15)//pb
          {
            fId = 5;
          }
          else//Hadron
          {
            if(MC_charge->GetLeaf("MCHadrons.charge")->GetValue(i)==1 && (MC_pid->GetLeaf("MCHadrons.pid")->GetValue(i)>7))
            {
              fId = 6;
            }
            else if(MC_charge->GetLeaf("MCHadrons.charge")->GetValue(i)==-1 && (MC_pid->GetLeaf("MCHadrons.pid")->GetValue(i)>7))
            {
              fId = 7;
            }
            else
            {
              continue;
            }
          }

          if(nu_MC)
          {
            if(fId == 2 || fId == 3)
              zBj_MC = sqrt(pow(MC_p->GetLeaf("MCHadrons.P")->GetValue(i),2)+pow(fM_K,2))/nu_MC;
            else if(fId == 4 || fId == 5)
              zBj_MC = sqrt(pow(MC_p->GetLeaf("MCHadrons.P")->GetValue(i),2)+pow(fM_p,2))/nu_MC;
            else
              zBj_MC = sqrt(pow(MC_p->GetLeaf("MCHadrons.P")->GetValue(i),2)+pow(fM_pi,2))/nu_MC;

            zBj_MC_unid = sqrt(pow(MC_p->GetLeaf("MCHadrons.P")->GetValue(i),2)+pow(fM_pi,2))/nu_MC;
          }
          else
          {
            zBj_MC = 0;
            zBj_MC_unid = 0;
          }

          if(!(0.2<zBj_MC && zBj_MC<0.85)) continue;

          if(0.2<zBj_MC && zBj_MC<0.25) zbin = 0;
          else if(0.25<zBj_MC && zBj_MC<0.30) zbin = 1;
          else if(0.30<zBj_MC && zBj_MC<0.35) zbin = 2;
          else if(0.35<zBj_MC && zBj_MC<0.40) zbin = 3;
          else if(0.40<zBj_MC && zBj_MC<0.45) zbin = 4;
          else if(0.45<zBj_MC && zBj_MC<0.50) zbin = 5;
          else if(0.50<zBj_MC && zBj_MC<0.55) zbin = 6;
          else if(0.55<zBj_MC && zBj_MC<0.60) zbin = 7;
          else if(0.60<zBj_MC && zBj_MC<0.65) zbin = 8;
          else if(0.65<zBj_MC && zBj_MC<0.70) zbin = 9;
          else if(0.70<zBj_MC && zBj_MC<0.75) zbin = 10;
          else zbin = 11;

          if(0.2<zBj_MC_unid && zBj_MC_unid<0.25) zbin_u = 0;
          else if(0.25<zBj_MC_unid && zBj_MC_unid<0.30) zbin_u = 1;
          else if(0.30<zBj_MC_unid && zBj_MC_unid<0.35) zbin_u = 2;
          else if(0.35<zBj_MC_unid && zBj_MC_unid<0.40) zbin_u = 3;
          else if(0.40<zBj_MC_unid && zBj_MC_unid<0.45) zbin_u = 4;
          else if(0.45<zBj_MC_unid && zBj_MC_unid<0.50) zbin_u = 5;
          else if(0.50<zBj_MC_unid && zBj_MC_unid<0.55) zbin_u = 6;
          else if(0.55<zBj_MC_unid && zBj_MC_unid<0.60) zbin_u = 7;
          else if(0.60<zBj_MC_unid && zBj_MC_unid<0.65) zbin_u = 8;
          else if(0.65<zBj_MC_unid && zBj_MC_unid<0.70) zbin_u = 9;
          else if(0.70<zBj_MC_unid && zBj_MC_unid<0.75) zbin_u = 10;
          else zbin_u = 11;

          // **********************************************************************

          // Save of hadrons

          if(fId==0)
          {
            if(!fFlag_MC[0][xbin_MC][ybin_MC][zbin_u])
            {
              fGnrt[xbin_MC][ybin_MC][zbin_u].tab[1][0][3] += 1;
              fMCHplus++;
              idMCrec[1][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zbin_u));
              prevMCrec[1][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),i));
              hidMCrec[1][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),fId));
              pMCrec[1][3].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),MC_p->GetLeaf("MCHadrons.P")->GetValue(i)));
              zMCrec[1][3].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zBj_MC_unid));
            }
            if(fFlag_MC[0][xbin_MC][ybin_MC][zbin]) continue;
            idMCrec[1][0].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zbin));
            prevMCrec[1][0].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),i));
            hidMCrec[1][0].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),fId));
            pMCrec[1][0].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),MC_p->GetLeaf("MCHadrons.P")->GetValue(i)));
            zMCrec[1][0].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zBj_MC));
            fMCPiplus++;
            fGnrt[xbin_MC][ybin_MC][zbin].tab[1][0][0] += 1;
          }
          else if(fId==1)
          {
            if(!fFlag_MC[0][xbin_MC][ybin_MC][zbin_u])
            {
              fGnrt[xbin_MC][ybin_MC][zbin_u].tab[0][0][3] += 1;
              fMCHminus++;
              idMCrec[0][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zbin_u));
              prevMCrec[0][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),i));
              hidMCrec[0][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),fId));
              pMCrec[0][3].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),MC_p->GetLeaf("MCHadrons.P")->GetValue(i)));
              zMCrec[0][3].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zBj_MC_unid));
            }
            if(fFlag_MC[0][xbin_MC][ybin_MC][zbin]) continue;
            idMCrec[0][0].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zbin));
            prevMCrec[0][0].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),i));
            hidMCrec[0][0].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),fId));
            pMCrec[0][0].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),MC_p->GetLeaf("MCHadrons.P")->GetValue(i)));
            zMCrec[0][0].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zBj_MC));
            fMCPiminus++;
            fGnrt[xbin_MC][ybin_MC][zbin].tab[0][0][0] += 1;
          }
          else if(fId==2)
          {
            if(!fFlag_MC[0][xbin_MC][ybin_MC][zbin_u])
            {
              fGnrt[xbin_MC][ybin_MC][zbin_u].tab[1][0][3] += 1;
              fMCHplus++;
              idMCrec[1][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zbin_u));
              prevMCrec[1][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),i));
              hidMCrec[1][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),fId));
              pMCrec[1][3].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),MC_p->GetLeaf("MCHadrons.P")->GetValue(i)));
              zMCrec[1][3].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zBj_MC_unid));
            }
            if(fFlag_MC[1][xbin_MC][ybin_MC][zbin]) continue;
            idMCrec[1][1].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zbin));
            prevMCrec[1][1].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),i));
            hidMCrec[1][1].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),fId));
            pMCrec[1][1].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),MC_p->GetLeaf("MCHadrons.P")->GetValue(i)));
            zMCrec[1][1].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zBj_MC));
            fMCKplus++;
            fGnrt[xbin_MC][ybin_MC][zbin].tab[1][0][1] += 1;
          }
          else if(fId==3)
          {
            if(!fFlag_MC[0][xbin_MC][ybin_MC][zbin_u])
            {
              fGnrt[xbin_MC][ybin_MC][zbin_u].tab[0][0][3] += 1;
              fMCHminus++;
              idMCrec[0][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zbin_u));
              prevMCrec[0][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),i));
              hidMCrec[0][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),fId));
              pMCrec[0][3].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),MC_p->GetLeaf("MCHadrons.P")->GetValue(i)));
              zMCrec[0][3].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zBj_MC_unid));
            }
            if(fFlag_MC[1][xbin_MC][ybin_MC][zbin]) continue;
            idMCrec[0][1].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zbin));
            prevMCrec[0][1].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),i));
            hidMCrec[0][1].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),fId));
            pMCrec[0][1].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),MC_p->GetLeaf("MCHadrons.P")->GetValue(i)));
            zMCrec[0][1].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zBj_MC));
            fMCKminus++;
            fGnrt[xbin_MC][ybin_MC][zbin].tab[0][0][1] += 1;
          }
          else if(fId==4)
          {
            if(!fFlag_MC[0][xbin_MC][ybin_MC][zbin_u])
            {
              fGnrt[xbin_MC][ybin_MC][zbin_u].tab[1][0][3] += 1;
              fMCHplus++;
              idMCrec[1][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zbin_u));
              prevMCrec[1][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),i));
              hidMCrec[1][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),fId));
              pMCrec[1][3].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),MC_p->GetLeaf("MCHadrons.P")->GetValue(i)));
              zMCrec[1][3].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zBj_MC_unid));
            }
            if(fFlag_MC[2][xbin_MC][ybin_MC][zbin]) continue;
            idMCrec[1][2].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zbin));
            prevMCrec[1][2].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),i));
            hidMCrec[1][2].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),fId));
            pMCrec[1][2].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),MC_p->GetLeaf("MCHadrons.P")->GetValue(i)));
            zMCrec[1][2].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zBj_MC));
            fMCPplus++;
            fGnrt[xbin_MC][ybin_MC][zbin].tab[1][0][2] += 1;
          }
          else if(fId==5)
          {
            if(!fFlag_MC[0][xbin_MC][ybin_MC][zbin_u])
            {
              fGnrt[xbin_MC][ybin_MC][zbin_u].tab[0][0][3] += 1;
              fMCHminus++;
              idMCrec[0][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zbin_u));
              prevMCrec[0][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),i));
              hidMCrec[0][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),fId));
              pMCrec[0][3].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),MC_p->GetLeaf("MCHadrons.P")->GetValue(i)));
              zMCrec[0][3].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zBj_MC_unid));
            }
            if(fFlag_MC[2][xbin_MC][ybin_MC][zbin]) continue;
            idMCrec[0][2].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zbin));
            prevMCrec[0][2].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),i));
            pMCrec[0][2].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),MC_p->GetLeaf("MCHadrons.P")->GetValue(i)));
            zMCrec[0][2].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zBj_MC));
            fMCPminus++;
            fGnrt[xbin_MC][ybin_MC][zbin].tab[0][0][2] += 1;
          }
          else if(fId==6)
          {
            if(fFlag_MC[0][xbin_MC][ybin_MC][zbin]) continue;
            idMCrec[1][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zbin));
            prevMCrec[1][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),i));
            hidMCrec[1][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),fId));
            pMCrec[1][3].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),MC_p->GetLeaf("MCHadrons.P")->GetValue(i)));
            zMCrec[1][3].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zBj_MC));
            fMCHplus++;
            fGnrt[xbin_MC][ybin_MC][zbin].tab[1][0][3] += 1;
          }
          else if(fId==7)
          {
            if(fFlag_MC[0][xbin_MC][ybin_MC][zbin]) continue;
            idMCrec[0][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zbin));
            prevMCrec[0][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),i));
            hidMCrec[0][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),fId));
            pMCrec[0][3].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),MC_p->GetLeaf("MCHadrons.P")->GetValue(i)));
            zMCrec[0][3].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zBj_MC));
            fMCHminus++;
            fGnrt[xbin_MC][ybin_MC][zbin].tab[0][0][3] += 1;
          }
          else
          {}

        }
      }
	    else
      {
	      for(int i=0; i<MC_p->GetLeaf("MCHadrons.P")->GetLen(); i++)
        {

          if(MC_pid->GetLeaf("MCHadrons.pid")->GetValue(i) == 8)//pi+
          {
            fId = 0;
          }
          else if(MC_pid->GetLeaf("MCHadrons.pid")->GetValue(i) == 9)//pi-
          {
            fId = 1;
          }
          else if(MC_pid->GetLeaf("MCHadrons.pid")->GetValue(i) == 11)//K+
          {
            fId = 2;
          }
          else if(MC_pid->GetLeaf("MCHadrons.pid")->GetValue(i) == 12)//K-
          {
            fId = 3;
          }
          else if(MC_pid->GetLeaf("MCHadrons.pid")->GetValue(i) == 14)//p
          {
            fId = 4;
          }
          else if(MC_pid->GetLeaf("MCHadrons.pid")->GetValue(i) == 15)//pb
          {
            fId = 5;
          }
          else//Hadron
          {
            if(MC_charge->GetLeaf("MCHadrons.charge")->GetValue(i)==1 && (MC_pid->GetLeaf("MCHadrons.pid")->GetValue(i)>7))
            {
              fId = 6;
            }
            else if(MC_charge->GetLeaf("MCHadrons.charge")->GetValue(i)==-1 && (MC_pid->GetLeaf("MCHadrons.pid")->GetValue(i)>7))
            {
              fId = 7;
            }
            else
            {
              continue;
            }
          }

          if(nu_MC)
          {
            if(fId == 2 || fId == 3)
              zBj_MC = sqrt(pow(MC_p->GetLeaf("MCHadrons.P")->GetValue(i),2)+pow(fM_K,2))/nu_MC;
            else if(fId == 4 || fId == 5)
              zBj_MC = sqrt(pow(MC_p->GetLeaf("MCHadrons.P")->GetValue(i),2)+pow(fM_p,2))/nu_MC;
            else
              zBj_MC = sqrt(pow(MC_p->GetLeaf("MCHadrons.P")->GetValue(i),2)+pow(fM_pi,2))/nu_MC;

              zBj_MC_unid = sqrt(pow(MC_p->GetLeaf("MCHadrons.P")->GetValue(i),2)+pow(fM_pi,2))/nu_MC;
          }
          else
          {
            zBj_MC = 0;
          }

          if(!(0.2<zBj_MC && zBj_MC<0.85)) continue;

          if(0.2<zBj_MC && zBj_MC<0.25) zbin = 0;
          else if(0.25<=zBj_MC && zBj_MC<0.30) zbin = 1;
          else if(0.30<=zBj_MC && zBj_MC<0.35) zbin = 2;
          else if(0.35<=zBj_MC && zBj_MC<0.40) zbin = 3;
          else if(0.40<=zBj_MC && zBj_MC<0.45) zbin = 4;
          else if(0.45<=zBj_MC && zBj_MC<0.50) zbin = 5;
          else if(0.50<=zBj_MC && zBj_MC<0.55) zbin = 6;
          else if(0.55<=zBj_MC && zBj_MC<0.60) zbin = 7;
          else if(0.60<=zBj_MC && zBj_MC<0.65) zbin = 8;
          else if(0.65<=zBj_MC && zBj_MC<0.70) zbin = 9;
          else if(0.70<=zBj_MC && zBj_MC<0.75) zbin = 10;
          else zbin = 11;

          if(0.2<zBj_MC_unid && zBj_MC_unid<0.25) zbin_u = 0;
          else if(0.25<zBj_MC_unid && zBj_MC_unid<0.30) zbin_u = 1;
          else if(0.30<zBj_MC_unid && zBj_MC_unid<0.35) zbin_u = 2;
          else if(0.35<zBj_MC_unid && zBj_MC_unid<0.40) zbin_u = 3;
          else if(0.40<zBj_MC_unid && zBj_MC_unid<0.45) zbin_u = 4;
          else if(0.45<zBj_MC_unid && zBj_MC_unid<0.50) zbin_u = 5;
          else if(0.50<zBj_MC_unid && zBj_MC_unid<0.55) zbin_u = 6;
          else if(0.55<zBj_MC_unid && zBj_MC_unid<0.60) zbin_u = 7;
          else if(0.60<zBj_MC_unid && zBj_MC_unid<0.65) zbin_u = 8;
          else if(0.65<zBj_MC_unid && zBj_MC_unid<0.70) zbin_u = 9;
          else if(0.70<zBj_MC_unid && zBj_MC_unid<0.75) zbin_u = 10;
          else zbin_u = 11;

          // **********************************************************************

          // Save of hadrons

          if(fId==0)
          {
            if(!fFlag_MC[0][xbin_MC][ybin_MC][zbin_u])
            {
              idMCrec[1][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zbin_u));
              prevMCrec[1][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),i));
              hidMCrec[1][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),fId));
              pMCrec[1][3].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),MC_p->GetLeaf("MCHadrons.P")->GetValue(i)));
              zMCrec[1][3].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zBj_MC_unid));
            }
            if(fFlag_MC[0][xbin_MC][ybin_MC][zbin]) continue;
            idMCrec[1][0].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zbin));
            prevMCrec[1][0].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),i));
            hidMCrec[1][0].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),fId));
            pMCrec[1][0].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),MC_p->GetLeaf("MCHadrons.P")->GetValue(i)));
            zMCrec[1][0].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zBj_MC));
          }
          else if(fId==1)
          {
            if(!fFlag_MC[0][xbin_MC][ybin_MC][zbin_u])
            {
              idMCrec[0][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zbin_u));
              prevMCrec[0][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),i));
              hidMCrec[0][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),fId));
              pMCrec[0][3].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),MC_p->GetLeaf("MCHadrons.P")->GetValue(i)));
              zMCrec[0][3].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zBj_MC_unid));
            }
            if(fFlag_MC[0][xbin_MC][ybin_MC][zbin]) continue;
            idMCrec[0][0].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zbin));
            prevMCrec[0][0].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),i));
            hidMCrec[0][0].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),fId));
            pMCrec[0][0].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),MC_p->GetLeaf("MCHadrons.P")->GetValue(i)));
            zMCrec[0][0].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zBj_MC));
          }
          else if(fId==2)
          {
            if(!fFlag_MC[0][xbin_MC][ybin_MC][zbin_u])
            {
              idMCrec[1][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zbin_u));
              prevMCrec[1][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),i));
              hidMCrec[1][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),fId));
              pMCrec[1][3].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),MC_p->GetLeaf("MCHadrons.P")->GetValue(i)));
              zMCrec[1][3].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zBj_MC_unid));
            }
            if(fFlag_MC[1][xbin_MC][ybin_MC][zbin]) continue;
            idMCrec[1][1].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zbin));
            prevMCrec[1][1].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),i));
            hidMCrec[1][1].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),fId));
            pMCrec[1][1].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),MC_p->GetLeaf("MCHadrons.P")->GetValue(i)));
            zMCrec[1][1].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zBj_MC));
          }
          else if(fId==3)
          {
            if(!fFlag_MC[0][xbin_MC][ybin_MC][zbin_u])
            {
              idMCrec[0][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zbin_u));
              prevMCrec[0][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),i));
              hidMCrec[0][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),fId));
              pMCrec[0][3].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),MC_p->GetLeaf("MCHadrons.P")->GetValue(i)));
              zMCrec[0][3].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zBj_MC_unid));
            }
            if(fFlag_MC[1][xbin_MC][ybin_MC][zbin]) continue;
            idMCrec[0][1].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zbin));
            prevMCrec[0][1].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),i));
            hidMCrec[0][1].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),fId));
            pMCrec[0][1].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),MC_p->GetLeaf("MCHadrons.P")->GetValue(i)));
            zMCrec[0][1].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zBj_MC));
          }
          else if(fId==4)
          {
            if(!fFlag_MC[0][xbin_MC][ybin_MC][zbin_u])
            {
              idMCrec[1][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zbin_u));
              prevMCrec[1][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),i));
              hidMCrec[1][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),fId));
              pMCrec[1][3].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),MC_p->GetLeaf("MCHadrons.P")->GetValue(i)));
              zMCrec[1][3].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zBj_MC_unid));
            }
            if(fFlag_MC[2][xbin_MC][ybin_MC][zbin]) continue;
            idMCrec[1][2].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zbin));
            prevMCrec[1][2].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),i));
            hidMCrec[1][2].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),fId));
            pMCrec[1][2].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),MC_p->GetLeaf("MCHadrons.P")->GetValue(i)));
            zMCrec[1][2].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zBj_MC));
          }
          else if(fId==5)
          {
            if(!fFlag_MC[0][xbin_MC][ybin_MC][zbin_u])
            {
              idMCrec[0][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zbin_u));
              prevMCrec[0][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),i));
              hidMCrec[0][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),fId));
              pMCrec[0][3].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),MC_p->GetLeaf("MCHadrons.P")->GetValue(i)));
              zMCrec[0][3].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zBj_MC_unid));
            }
            if(fFlag_MC[2][xbin_MC][ybin_MC][zbin]) continue;
            idMCrec[0][2].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zbin));
            prevMCrec[0][2].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),i));
            hidMCrec[0][2].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),fId));
            pMCrec[0][2].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),MC_p->GetLeaf("MCHadrons.P")->GetValue(i)));
            zMCrec[0][2].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zBj_MC));
          }
          else if(fId==6)
          {
            if(fFlag_MC[0][xbin_MC][ybin_MC][zbin]) continue;
            idMCrec[1][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zbin));
            prevMCrec[1][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),i));
            hidMCrec[1][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),fId));
            pMCrec[1][3].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),MC_p->GetLeaf("MCHadrons.P")->GetValue(i)));
            zMCrec[1][3].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zBj_MC));
          }
          else if(fId==7)
          {
            if(fFlag_MC[0][xbin_MC][ybin_MC][zbin]) continue;
            idMCrec[0][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zbin));
            prevMCrec[0][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),i));
            hidMCrec[0][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),fId));
            pMCrec[0][3].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),MC_p->GetLeaf("MCHadrons.P")->GetValue(i)));
            zMCrec[0][3].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zBj_MC));
          }
          else
          {}

        }
      }


      // -----------------------------------------------------------------------
      //  Data -----------------------------------------------------------------
      // -----------------------------------------------------------------------

      if(fAllDISflag)
      {
        for(int i=0; i<p->GetLeaf("Hadrons.P")->GetLen(); i++)
        {
          // **********************************************************************

          // Hadron identification cuts ------------------------------------------

          if(MCpid->GetLeaf("Hadrons.MCpid")->GetValue(i) == 8)//pi+
          {
            fId = 0;
          }
          else if(MCpid->GetLeaf("Hadrons.MCpid")->GetValue(i) == 9)//pi-
          {
            fId = 1;
          }
          else if(MCpid->GetLeaf("Hadrons.MCpid")->GetValue(i) == 11)//K+
          {
            fId = 2;
          }
          else if(MCpid->GetLeaf("Hadrons.MCpid")->GetValue(i) == 12)//K-
          {
            fId = 3;
          }
          else if(MCpid->GetLeaf("Hadrons.MCpid")->GetValue(i) == 14)//p
          {
            fId = 4;
          }
          else if(MCpid->GetLeaf("Hadrons.MCpid")->GetValue(i) == 15)//pb
          {
            fId = 5;
          }
          else if(MCpid->GetLeaf("Hadrons.MCpid")->GetValue(i) == 3)//e-
          {
            fId = 8;
          }
          else if(MCpid->GetLeaf("Hadrons.MCpid")->GetValue(i) == 2)//e+
          {
            fId = 9;
          }
          else//Hadron
          {
            if(charge->GetLeaf("Hadrons.charge")->GetValue(i)==1 && MCpid->GetLeaf("Hadrons.MCpid")->GetValue(i)>7)
            {
              fId = 6;
            }
            else if(charge->GetLeaf("Hadrons.charge")->GetValue(i)==-1 && MCpid->GetLeaf("Hadrons.MCpid")->GetValue(i)>7)
            {
              fId = 7;
            }
            else
            {
              continue;
            }
          }

          // **********************************************************************

          // z calculation
          if(nu)
          {
            if(fId == 2 || fId == 3)
              zBj = sqrt(pow(p->GetLeaf("Hadrons.P")->GetValue(i),2)+pow(fM_K,2))/nu;
            else if(fId == 4 || fId == 5)
              zBj = sqrt(pow(p->GetLeaf("Hadrons.P")->GetValue(i),2)+pow(fM_p,2))/nu;
            else
              zBj = sqrt(pow(p->GetLeaf("Hadrons.P")->GetValue(i),2)+pow(fM_pi,2))/nu;

            zBj_unid = sqrt(pow(p->GetLeaf("Hadrons.P")->GetValue(i),2)+pow(fM_pi,2))/nu;
          }
          else
          {
            zBj = 0;
            zBj_unid = 0;
          }

          if(0.1<zBj) fKinematicsMC[6]->Fill(abs(ph->GetLeaf("Hadrons.ph")->GetValue(i)));

          // Maximum radiation length cumulated
          if(!(hXX0->GetLeaf("Hadrons.XX0")->GetValue(i) < 15)) continue;
          fXX0test++;

          // Momentum cut (12 GeV to 40 GeV, increasing to 3 GeV to 40 GeV)
          if(!(MOMENTUM<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<40)) continue;
          fMom++;

          // Theta cut
          if(!(0.01<thRICH->GetLeaf("Hadrons.thRICH")->GetValue(i) && thRICH->GetLeaf("Hadrons.thRICH")->GetValue(i)<0.12)) continue;
          fTRICH++;

          // RICH position cut
          if(!(pow(RICHx->GetLeaf("Hadrons.RICHx")->GetValue(i),2)+pow(RICHy->GetLeaf("Hadrons.RICHy")->GetValue(i),2)>25)) continue;
          fPosRICH++;

          // z cut
          if(!(0.2<zBj && zBj<0.85)) continue;

          if(0.2<zBj && zBj<0.25) zbin = 0;
          else if(0.25<=zBj && zBj<0.30) zbin = 1;
          else if(0.30<=zBj && zBj<0.35) zbin = 2;
          else if(0.35<=zBj && zBj<0.40) zbin = 3;
          else if(0.40<=zBj && zBj<0.45) zbin = 4;
          else if(0.45<=zBj && zBj<0.50) zbin = 5;
          else if(0.50<=zBj && zBj<0.55) zbin = 6;
          else if(0.55<=zBj && zBj<0.60) zbin = 7;
          else if(0.60<=zBj && zBj<0.65) zbin = 8;
          else if(0.65<=zBj && zBj<0.70) zbin = 9;
          else if(0.70<=zBj && zBj<0.75) zbin = 10;
          else zbin = 11;

          if(0.2<zBj_unid && zBj_unid<0.25) zbin_u = 0;
          else if(0.25<=zBj_unid && zBj_unid<0.30) zbin_u = 1;
          else if(0.30<=zBj_unid && zBj_unid<0.35) zbin_u = 2;
          else if(0.35<=zBj_unid && zBj_unid<0.40) zbin_u = 3;
          else if(0.40<=zBj_unid && zBj_unid<0.45) zbin_u = 4;
          else if(0.45<=zBj_unid && zBj_unid<0.50) zbin_u = 5;
          else if(0.50<=zBj_unid && zBj_unid<0.55) zbin_u = 6;
          else if(0.55<=zBj_unid && zBj_unid<0.60) zbin_u = 7;
          else if(0.60<=zBj_unid && zBj_unid<0.65) zbin_u = 8;
          else if(0.65<=zBj_unid && zBj_unid<0.70) zbin_u = 9;
          else if(0.70<=zBj_unid && zBj_unid<0.75) zbin_u = 10;
          else zbin_u = 11;



          // **********************************************************************

          // Save of hadrons

          map<int,int>::iterator it;

          if(fId==0)
          {
            it = idMCrec[1][3].find(i);
            if(!fFlag[0][xbin][ybin][zbin_u])
            {
              fHplus++;
              fRcstr[xbin][ybin][zbin_u].tab[1][0][3] += 1;
              if(it!=idMCrec[1][3].end())
              {
                if(zbin_u == idMCrec[1][3][i] && xbin==xbin_MC && ybin==ybin_MC)
                {
                  fRcstr_c[xbin][ybin][zbin_u].tab[1][0][3]++;
                }
              }
            }
            it = idMCrec[1][0].find(i);
            if(fFlag[0][xbin][ybin][zbin]) continue;
            fPiplus++;
            fRcstr[xbin][ybin][zbin].tab[1][0][0] += 1;
            if(it!=idMCrec[1][0].end())
            {
              if(zbin == idMCrec[1][0][i] && xbin==xbin_MC && ybin==ybin_MC)
              {
                fRcstr_c[xbin][ybin][zbin].tab[1][0][0]++;
                idMCrec[1][0].erase(i);
              }
            }
          }
          else if(fId==1)
          {
            it = idMCrec[0][3].find(i);
            if(!fFlag[0][xbin][ybin][zbin_u])
            {
              fHminus++;
              fRcstr[xbin][ybin][zbin_u].tab[0][0][3] += 1;
              if(it!=idMCrec[0][3].end())
              {
                if(zbin_u == idMCrec[0][3][i] && xbin==xbin_MC && ybin==ybin_MC)
                {
                  fRcstr_c[xbin][ybin][zbin_u].tab[0][0][3]++;
                }
              }
            }
            it = idMCrec[0][0].find(i);
            if(fFlag[0][xbin][ybin][zbin]) continue;
            fPiminus++;
            fRcstr[xbin][ybin][zbin].tab[0][0][0] += 1;
            if(it!=idMCrec[0][0].end())
            {
              if(zbin == idMCrec[0][0][i] && xbin==xbin_MC && ybin==ybin_MC)
              {
                fRcstr_c[xbin][ybin][zbin].tab[0][0][0]++;
                idMCrec[0][0].erase(i);
              }
            }
          }
          else if(fId==2)
          {
            it = idMCrec[1][3].find(i);
            if(!fFlag[0][xbin][ybin][zbin_u])
            {
              fHplus++;
              fRcstr[xbin][ybin][zbin_u].tab[1][0][3] += 1;
              if(it!=idMCrec[1][3].end())
              {
                if(zbin_u == idMCrec[1][3][i] && xbin==xbin_MC && ybin==ybin_MC)
                {
                  fRcstr_c[xbin][ybin][zbin_u].tab[1][0][3]++;
                }
              }
            }
            it = idMCrec[1][1].find(i);
            if(fFlag[1][xbin][ybin][zbin]) continue;
            fKplus++;
            fRcstr[xbin][ybin][zbin].tab[1][0][1] += 1;
            if(it!=idMCrec[1][1].end())
            {
              if(zbin == idMCrec[1][1][i] && xbin==xbin_MC && ybin==ybin_MC)
              {
                fRcstr_c[xbin][ybin][zbin].tab[1][0][1]++;
                idMCrec[1][1].erase(i);
              }
            }
          }
          else if(fId==3)
          {
            it = idMCrec[0][3].find(i);
            if(!fFlag[0][xbin][ybin][zbin_u])
            {
              fHminus++;
              fRcstr[xbin][ybin][zbin_u].tab[0][0][3] += 1;
              if(it!=idMCrec[0][3].end())
              {
                if(zbin_u == idMCrec[0][3][i] && xbin==xbin_MC && ybin==ybin_MC)
                {
                  fRcstr_c[xbin][ybin][zbin_u].tab[0][0][3]++;
                }
              }
            }
            it = idMCrec[0][1].find(i);
            if(fFlag[1][xbin][ybin][zbin]) continue;
            fKminus++;
            fRcstr[xbin][ybin][zbin].tab[0][0][1] += 1;
            if(it!=idMCrec[0][1].end())
            {
              if(zbin == idMCrec[0][1][i] && xbin==xbin_MC && ybin==ybin_MC)
              {
                fRcstr_c[xbin][ybin][zbin].tab[0][0][1]++;
                idMCrec[0][1].erase(i);
              }
            }
          }
          else if(fId==4)
          {
            it = idMCrec[1][3].find(i);
            if(!fFlag[0][xbin][ybin][zbin_u])
            {
              fHplus++;
              fRcstr[xbin][ybin][zbin_u].tab[1][0][3] += 1;
              if(it!=idMCrec[1][3].end())
              {
                if(zbin_u == idMCrec[1][3][i] && xbin==xbin_MC && ybin==ybin_MC)
                {
                  fRcstr_c[xbin][ybin][zbin_u].tab[1][0][3]++;
                }
              }
            }
            it = idMCrec[1][2].find(i);
            if(fFlag[2][xbin][ybin][zbin]) continue;
            fPplus++;
            fRcstr[xbin][ybin][zbin].tab[1][0][2] += 1;
            if(it!=idMCrec[1][2].end())
            {
              if(zbin == idMCrec[1][2][i] && xbin==xbin_MC && ybin==ybin_MC)
              {
                fRcstr_c[xbin][ybin][zbin].tab[1][0][2]++;
                idMCrec[1][2].erase(i);
              }
            }
          }
          else if(fId==5)
          {
            it = idMCrec[0][3].find(i);
            if(!fFlag[0][xbin][ybin][zbin_u])
            {
              fHminus++;
              fRcstr[xbin][ybin][zbin_u].tab[0][0][3] += 1;
              if(it!=idMCrec[0][3].end())
              {
                if(zbin_u == idMCrec[0][3][i] && xbin==xbin_MC && ybin==ybin_MC)
                {
                  fRcstr_c[xbin][ybin][zbin_u].tab[0][0][3]++;
                }
              }
            }
            it = idMCrec[0][2].find(i);
            if(fFlag[2][xbin][ybin][zbin]) continue;
            fPminus++;
            fRcstr[xbin][ybin][zbin].tab[0][0][2] += 1;
            if(it!=idMCrec[0][2].end())
            {
              if(zbin == idMCrec[0][2][i] && xbin==xbin_MC && ybin==ybin_MC)
              {
                fRcstr_c[xbin][ybin][zbin].tab[0][0][2]++;
                idMCrec[0][2].erase(i);
              }
            }
          }
          else if(fId==6)
          {
            if(fFlag[0][xbin][ybin][zbin]) continue;
            fHplus++;
            fRcstr[xbin][ybin][zbin].tab[1][0][3] += 1;
            it = idMCrec[1][3].find(i);
            if(it!=idMCrec[1][3].end())
            {
              if(zbin == idMCrec[1][3][i] && xbin==xbin_MC && ybin==ybin_MC)
              {
                fRcstr_c[xbin][ybin][zbin].tab[1][0][3]++;
                idMCrec[1][3].erase(i);
              }
            }
          }
          else if(fId==7)
          {
            if(fFlag[0][xbin][ybin][zbin]) continue;
            fHminus++;
            fRcstr[xbin][ybin][zbin].tab[0][0][3] += 1;
            it = idMCrec[0][3].find(i);
            if(it!=idMCrec[0][3].end())
            {
              if(zbin == idMCrec[0][3][i] && xbin==xbin_MC && ybin==ybin_MC)
              {
                fRcstr_c[xbin][ybin][zbin].tab[0][0][3]++;
                idMCrec[0][3].erase(i);
              }
            }
          }
          else if(fId==8)
          {
            if(fFlag[0][xbin][ybin][zbin]) continue;
            fHminus++; fPiminus++;
            fRcstr[xbin][ybin][zbin].tab[0][0][0] += 1;
            fRcstr[xbin][ybin][zbin].tab[0][0][3] += 1;
          }
          else if(fId==9)
          {
            if(fFlag[0][xbin][ybin][zbin]) continue;
            fHplus++; fPiplus++;
            fRcstr[xbin][ybin][zbin].tab[1][0][0] += 1;
            fRcstr[xbin][ybin][zbin].tab[1][0][3] += 1;
          }
          else
          {
            continue;
          }

          fKinematicsMC[3]->Fill(zBj);
        }
      }

      for(int cc=0; cc<2; cc++)
      {
    	for(int hh=0; hh<4; hh++)
        {
          idMCrec[cc][hh].clear();
          hidMCrec[cc][hh].clear();
          pMCrec[cc][hh].clear();
          zMCrec[cc][hh].clear();
          prevMCrec[cc][hh].clear();
        }
      }


    }

    cout << "-> Finished processing file " << filename << " <-\n" << endl;

    delete f;
  }

  for(int i=0; i<int(fQ2kinMC.size()); i++)
  {
    fKinematicsMC[0]->Fill(fQ2kinMC[i]);
    fKinematicsMC[1]->Fill(fXBjkinMC[i]);
    fKinematicsMC[2]->Fill(fYBjkinMC[i]);
    fKinematicsMC[4]->Fill(fWBjkinMC[i]);
    fKinematicsMC[5]->Fill(fNukinMC[i]);
  }

}

void RDextraction(string pFilelist)
{

  //Kinematics
  Double_t Q2 = 0;
  Double_t xBj = 0;
  Double_t yBj = 0;
  Double_t zBj = 0;
  Double_t wBj = 0;
  Double_t nu = 0;

  load_rich_mat(mat_RICH_name, err_RICH_name);

  //cout << pi_sigma_uni[0][0] << " " << pi_sigma_uni[1][1] << " " << pi_sigma_uni[2][2] << endl;

  // Target cells
  if(Y2012) InitTargetFile(target_file_2012);
  else if(Y2016) InitTargetFile(target_file_2016);

  //----------------------------------------------------------------------------
  //--------- nu cut prep ------------------------------------------------------
  //----------------------------------------------------------------------------

  for(int i=0; i<12; i++)
  {
    fNu_max[1][i] = sqrt(pow(40,2)+pow(fM_K,2))/fZrange[i+1];
    fNu_min[1][i] = sqrt(pow(MOMENTUM,2)+pow(fM_K,2))/fZrange[i];

    fNu_max[2][i] = sqrt(pow(40,2)+pow(fM_p,2))/fZrange[i+1];
    fNu_min[2][i] = sqrt(pow(MOMENTUM,2)+pow(fM_p,2))/fZrange[i];

    fNu_max[0][i] = sqrt(pow(40,2)+pow(fM_pi,2))/fZrange[i+1];
    fNu_min[0][i] = sqrt(pow(MOMENTUM,2)+pow(fM_pi,2))/fZrange[i];
  }

  // List of files

  ifstream list(pFilelist);
  string filename;

  while(list >> filename)
  {
    TFile *f;

    cout << ".. Processing file " << filename << " .." << endl;
    f = TFile::Open(filename.c_str());

    if(!f) continue;

    TTree* tree = (TTree*) f->Get("DISEvtTree");

    // ---------------------------------------------------------------------------
    // --------- Reading of the TTree --------------------------------------------
    // ---------------------------------------------------------------------------

    //DISEvt
    TBranch *runNo = (TBranch*) tree->FindBranch("runNo");
    TBranch *spillNo = (TBranch*) tree->FindBranch("spillNo");
    TBranch *evtInSpill = (TBranch*) tree->FindBranch("evtInSpill");
    TBranch *trigMask = (TBranch*) tree->FindBranch("trigMask");
    TBranch *evNo = (TBranch*) tree->FindBranch("evNo");
    TBranch *x = (TBranch*) tree->FindBranch("x");
    TBranch *y = (TBranch*) tree->FindBranch("y");
    TBranch *z = (TBranch*) tree->FindBranch("z");
    TBranch *p0x = (TBranch*) tree->FindBranch("p0x");
    TBranch *p0y = (TBranch*) tree->FindBranch("p0y");
    TBranch *p0z = (TBranch*) tree->FindBranch("p0z");
    TBranch *p1x = (TBranch*) tree->FindBranch("p1x");
    TBranch *p1y = (TBranch*) tree->FindBranch("p1y");
    TBranch *p1z = (TBranch*) tree->FindBranch("p1z");
    TBranch *E_beam = (TBranch*) tree->FindBranch("E_beam");
    TBranch *E_mu_prim = (TBranch*) tree->FindBranch("E_mu_prim");
    TBranch *XX0 = (TBranch*) tree->FindBranch("XX0");
    TBranch *HM04x = (TBranch*) tree->FindBranch("HM04x");
    TBranch *HM04y = (TBranch*) tree->FindBranch("HM04y");
    TBranch *HM05x = (TBranch*) tree->FindBranch("HM05x");
    TBranch *HM05y = (TBranch*) tree->FindBranch("HM05y");
    TBranch *HO03x = (TBranch*) tree->FindBranch("HO03x");
    TBranch *HO03y = (TBranch*) tree->FindBranch("HO03y");
    TBranch *HO04x = (TBranch*) tree->FindBranch("HO04x");
    TBranch *HO04y = (TBranch*) tree->FindBranch("HO04y");
    TBranch *saved = (TBranch*) tree->FindBranch("saved");
    TBranch *cellsCrossed = (TBranch*) tree->FindBranch("cellsCrossed");
    TBranch *backPropFlag = (TBranch*) tree->FindBranch("backPropFlag");

    //Hadrons
    TBranch *p = (TBranch*) tree->FindBranch("Hadrons.P");
    TBranch *th = (TBranch*) tree->FindBranch("Hadrons.th");
    TBranch *ph = (TBranch*) tree->FindBranch("Hadrons.ph");
    TBranch *hXX0 = (TBranch*) tree->FindBranch("Hadrons.XX0");
    TBranch *inHCALacc = (TBranch*) tree->FindBranch("Hadrons.inHCALacc");
    TBranch *HCAL = (TBranch*) tree->FindBranch("Hadrons.HCAL");
    TBranch *charge = (TBranch*) tree->FindBranch("Hadrons.charge");
    TBranch *thRICH = (TBranch*) tree->FindBranch("Hadrons.thRICH");
    TBranch *LH = (TBranch*) tree->FindBranch("Hadrons.LH");
    TBranch *MCpid = (TBranch*) tree->FindBranch("Hadrons.MCpid");
    TBranch *MM01x = (TBranch*) tree->FindBranch("Hadrons.MM01x");
    TBranch *MM01y = (TBranch*) tree->FindBranch("Hadrons.MM01y");
    TBranch *MM02x = (TBranch*) tree->FindBranch("Hadrons.MM02x");
    TBranch *MM02y = (TBranch*) tree->FindBranch("Hadrons.MM02y");
    TBranch *MM03x = (TBranch*) tree->FindBranch("Hadrons.MM03x");
    TBranch *MM03y = (TBranch*) tree->FindBranch("Hadrons.MM03y");
    TBranch *Z2Ax = (TBranch*) tree->FindBranch("Hadrons.Z2Ax");
    TBranch *Z2Ay = (TBranch*) tree->FindBranch("Hadrons.Z2Ay");
    TBranch *Z2Bx = (TBranch*) tree->FindBranch("Hadrons.Z2Bx");
    TBranch *Z2By = (TBranch*) tree->FindBranch("Hadrons.Z2By");
    TBranch *RICHx = (TBranch*) tree->FindBranch("Hadrons.RICHx");
    TBranch *RICHy = (TBranch*) tree->FindBranch("Hadrons.RICHy");

    // Loopy loop over the events
    Int_t N = (Int_t) tree->GetEntries();

    vector<Pvsz> Pvszlocal;
    vector<Pvsz> Pvszloose;
    vector<Pvsz> Pvszsevere;
    vector<Pvsz> Pvsz_errlocal;
    vector<Double_t> XBjlocal;
    vector<Double_t> YBjlocal;
    vector<Double_t> Q2local;
    vector<Double_t> XBjloose;
    vector<Double_t> YBjloose;
    vector<Double_t> Q2loose;
    vector<Double_t> XBjsevere;
    vector<Double_t> YBjsevere;
    vector<Double_t> Q2severe;

    for (Int_t ip = 0; ip < N; ip++)
    {

      printProgress(ip,N);

      //DISEvt
      runNo->GetEntry(ip);
      spillNo->GetEntry(ip);
      evtInSpill->GetEntry(ip);
      trigMask->GetEntry(ip);
      evNo->GetEntry(ip);
      x->GetEntry(ip);
      y->GetEntry(ip);
      z->GetEntry(ip);
      p0x->GetEntry(ip);
      p0y->GetEntry(ip);
      p0z->GetEntry(ip);
      p1x->GetEntry(ip);
      p1y->GetEntry(ip);
      p1z->GetEntry(ip);
      E_beam->GetEntry(ip);
      E_mu_prim->GetEntry(ip);
      XX0->GetEntry(ip);
      HM04x->GetEntry(ip);
      HM04y->GetEntry(ip);
      HM05x->GetEntry(ip);
      HM05y->GetEntry(ip);
      HO03x->GetEntry(ip);
      HO03y->GetEntry(ip);
      HO04x->GetEntry(ip);
      HO04y->GetEntry(ip);
      saved->GetEntry(ip);
      cellsCrossed->GetEntry(ip);
      backPropFlag->GetEntry(ip);

      //Hadrons
      p->GetEntry(ip);
      th->GetEntry(ip);
      ph->GetEntry(ip);
      hXX0->GetEntry(ip);
      inHCALacc->GetEntry(ip);
      HCAL->GetEntry(ip);
      charge->GetEntry(ip);
      thRICH->GetEntry(ip);
      LH->GetEntry(ip);
      MCpid->GetEntry(ip);
      MM01x->GetEntry(ip);
      MM01y->GetEntry(ip);
      MM02x->GetEntry(ip);
      MM02y->GetEntry(ip);
      MM03x->GetEntry(ip);
      MM03y->GetEntry(ip);
      Z2Ax->GetEntry(ip);
      Z2Ay->GetEntry(ip);
      Z2Bx->GetEntry(ip);
      Z2By->GetEntry(ip);
      RICHx->GetEntry(ip);
      RICHy->GetEntry(ip);

      // -------------------------------------------------------------------------
      // --------- Calculation ---------------------------------------------------
      // -------------------------------------------------------------------------

      Double_t zlab = z->GetLeaf("z")->GetValue();
      Int_t zlabbin;

      if(Y2012 || RCUTSTUDY_ON)
      {
        if(!(-311.2<zlab && zlab<-71.2)) continue;

        if(-311.2<=zlab && zlab<-301.2) zlabbin = 0;
        else if(-301.2<=zlab && zlab<-291.2) zlabbin = 1;
        else if(-291.2<=zlab && zlab<-281.2) zlabbin = 2;
        else if(-281.2<=zlab && zlab<-271.2) zlabbin = 3;
        else if(-271.2<=zlab && zlab<-261.2) zlabbin = 4;
        else if(-261.2<=zlab && zlab<-251.2) zlabbin = 5;
        else if(-251.2<=zlab && zlab<-241.2) zlabbin = 6;
        else if(-241.2<=zlab && zlab<-231.2) zlabbin = 7;
        else if(-231.2<=zlab && zlab<-221.2) zlabbin = 8;
        else if(-221.2<=zlab && zlab<-211.2) zlabbin = 9;
        else if(-211.2<=zlab && zlab<-201.2) zlabbin = 10;
        else if(-201.2<=zlab && zlab<-191.2) zlabbin = 11;
        else if(-191.2<=zlab && zlab<-181.2) zlabbin = 12;
        else if(-181.2<=zlab && zlab<-171.2) zlabbin = 13;
        else if(-171.2<=zlab && zlab<-161.2) zlabbin = 14;
        else if(-161.2<=zlab && zlab<-151.2) zlabbin = 15;
        else if(-151.2<=zlab && zlab<-141.2) zlabbin = 16;
        else if(-141.2<=zlab && zlab<-131.2) zlabbin = 17;
        else if(-131.2<=zlab && zlab<-121.2) zlabbin = 18;
        else if(-121.2<=zlab && zlab<-111.2) zlabbin = 19;
        else if(-111.2<=zlab && zlab<-101.2) zlabbin = 20;
        else if(-101.2<=zlab && zlab<-91.2) zlabbin = 21;
        else if(-91.2<=zlab && zlab<-81.2) zlabbin = 22;
        else zlabbin = 23;
      }

      //--------------------------------------------------------------------------
      //--------- Kinematics -----------------------------------------------------
      //--------------------------------------------------------------------------

      Q2 = 2.*( E_beam->GetLeaf("E_beam")->GetValue()*E_mu_prim->GetLeaf("E_mu_prim")->GetValue()
           - p0x->GetLeaf("p0x")->GetValue()*p1x->GetLeaf("p1x")->GetValue()
           - p0y->GetLeaf("p0y")->GetValue()*p1y->GetLeaf("p1y")->GetValue()
           - p0z->GetLeaf("p0z")->GetValue()*p1z->GetLeaf("p1z")->GetValue()
           - pow(fM_mu,2));

      nu = E_beam->GetLeaf("E_beam")->GetValue() - E_mu_prim->GetLeaf("E_mu_prim")->GetValue();

      if(E_beam->GetLeaf("E_beam")->GetValue() != 0)
        yBj = nu/E_beam->GetLeaf("E_beam")->GetValue();
      else
        yBj = 0;

      if(nu != 0)
      {
        xBj = Q2/(2*fM_p*nu);
      }
      else
      {
        xBj = 0;
      }

      if(xBj != 0)
        wBj = pow(fM_p,2) + Q2*(1-xBj)/xBj;
      else
        wBj = 0;

      int trig= trigMask->GetLeaf("trigMask")->GetValue();


      //2006 ---

      if(Y2006)
      {
        if ((trig&256) && HM05x->GetLeaf("HM05x")->GetValue()<(HM05y->GetLeaf("HM05y")->GetValue()>0 ? 14.55-0.15 : 22.02864-0.12864) )
        {
          trig -= 256;
        }
      }

      //--------------------------------------------------------------------------
      //--------- Target ---------------------------------------------------------
      //--------------------------------------------------------------------------

      //MC target position new
      static const Double_t dz = 2;

      static const Double_t mcxU = -0.085;
      static const Double_t mcyU = 0.33;
      static const Double_t mczU_1 = -65+dz+4;
      //static const double mczU_2 = -35+dz;
      //static const double mczC_1 = -30+dz+8;
      //static const double mczC_2 = 30+dz;
      static const Double_t mcxD = -0.085;
      static const Double_t mcyD = 0.33;
      //static const double mczD_1 = 35+dz+2;
      static const Double_t mczD_2 = 65+dz;

      double mcR    = 1.4;
      //double mcyCUT = 1.4;

      //target position data 2006
      static const Double_t xU = -0.1;
      static const Double_t yU = 0.33;
      static const Double_t zU_1 = -65+dz+4;
      //static const double zU_2 = -35+dz;
      //static const double zC_1 = -30+dz+8;
      //static const double zC_2 =  30+dz;
      static const Double_t xD = -0.07;
      static const Double_t yD = 0.33;
      //static const double zD_1 =  35+dz+2;
      static const Double_t zD_2 =  65+dz;

      Double_t R    = 1.4;//1.4;
      Double_t yCUT = 1.4;

      Double_t mcxC = (mcxD-mcxU) * (mczU_1-z->GetLeaf("z")->GetValue()) / (mczU_1-mczD_2) + mcxU;
      Double_t mcyC = (mcyD-mcyU) * (mczU_1-z->GetLeaf("z")->GetValue()) / (mczU_1-mczD_2) + mcyU;
      Double_t mcr = sqrt( (x->GetLeaf("x")->GetValue()-mcxC)*(x->GetLeaf("x")->GetValue()-mcxC)
                    + (y->GetLeaf("y")->GetValue()-mcyC)*(y->GetLeaf("y")->GetValue()-mcyC) );
      Double_t xC = (xD-xU) * (zU_1-z->GetLeaf("z")->GetValue()) / (zU_1-zD_2) + xU;
      Double_t yC = (yD-yU) * (zU_1-z->GetLeaf("z")->GetValue()) / (zU_1-zD_2) + yU;
      Double_t r = sqrt( (x->GetLeaf("x")->GetValue()-xC)*(x->GetLeaf("x")->GetValue()-xC)
                    + (y->GetLeaf("y")->GetValue()-yC)*(y->GetLeaf("y")->GetValue()-yC) );

      //2006 ---


      // -------------------------------------------------------------------------
      // --------- DIS Selection -------------------------------------------------
      // -------------------------------------------------------------------------

      // Best Primary Vertex
      fBP++;

      // Reconstructed muon
      if(!(0<E_beam->GetLeaf("E_beam")->GetValue())) continue;
      fRmu++;

      Double_t mxc, myc;

      if(Y2012 || RCUTSTUDY_ON)
      {
        CellCenter(z->GetLeaf("z")->GetValue(), mxc, myc);
      }

      //Rcut study ---
      if(RCUTSTUDY_ON)
      {
        fRstudy[zlabbin].vec.push_back(sqrt((pow(x->GetLeaf("x")->GetValue()-mxc,2)+pow(y->GetLeaf("y")->GetValue()-myc,2))));
        fRstudy_xy[zlabbin].vec[0].push_back(x->GetLeaf("x")->GetValue());
        fRstudy_xy[zlabbin].vec[1].push_back(y->GetLeaf("y")->GetValue());

        if(!InTarget(x->GetLeaf("x")->GetValue(),y->GetLeaf("y")->GetValue(),z->GetLeaf("z")->GetValue(),fRcutval[zlabbin])) continue;
        fTarg++;

        fR_xy[zlabbin].vec[0].push_back(x->GetLeaf("x")->GetValue());
        fR_xy[zlabbin].vec[1].push_back(y->GetLeaf("y")->GetValue());
      }
      //Rcut study ---


      //BMS (reconstructed beam track)
      if((backPropFlag->GetLeaf("backPropFlag")->GetValue())) continue;
      fBMS++;

      // Energy of the muon beam
      if(!(140<E_beam->GetLeaf("E_beam")->GetValue() && E_beam->GetLeaf("E_beam")->GetValue()<180)) continue;
      fBEC++;

      //2006 ---
      if(Y2006)
      {
      // Z coordinate within target regions
      if(!((-56<z->GetLeaf("z")->GetValue() && z->GetLeaf("z")->GetValue()<-35)
            ||(-20<z->GetLeaf("z")->GetValue() && z->GetLeaf("z")->GetValue()<31)
            ||(43<z->GetLeaf("z")->GetValue() && z->GetLeaf("z")->GetValue()<66))) continue;
      if(!(mcr < mcR &&  (y->GetLeaf("y")->GetValue()-mcyC)<yCUT
           && r < R
           &&  (y->GetLeaf("y")->GetValue()-yC)<yCUT
           && ((z->GetLeaf("z")->GetValue()>(-65+2+7) && z->GetLeaf("z")->GetValue()<(-35+2-2))
                ||(z->GetLeaf("z")->GetValue() > (-30+2+8) && z->GetLeaf("z")->GetValue() < (30+2-1))
                ||(z->GetLeaf("z")->GetValue() > (35+2+6) && z->GetLeaf("z")->GetValue() < (65+2-1))))) continue;
      }
      //2006 ---
      //2012 ---
      else if(Y2012)
      {
        if(!InTarget(x->GetLeaf("x")->GetValue(),y->GetLeaf("y")->GetValue(),z->GetLeaf("z")->GetValue(),fRcutval[zlabbin])) continue;
      }
      //2012 ---
      //2016 ---
      else if(Y2016)
      {
        if(!InTarget(x->GetLeaf("x")->GetValue(),y->GetLeaf("y")->GetValue(),z->GetLeaf("z")->GetValue(),1.5)) continue;
      }
      //2016 ---
      fTarg++;

      // Cells crossing
      //if(!(cellsCrossed->GetLeaf("cellsCrossed")->GetValue())) continue;
      //fCell++;

      // IM/O triggers
      //2006 ---
      if(Y2006)
      {
        if(!(trig&8 || trig&256)) continue;
      }
      //2006 ---
      //2012 ---
      else if(Y2012)
      {
        if(!(trig&2 || trig&4 || trig&8)) continue;
      }
      //2012 ---
      //2016 ---
      else if(Y2016)
      {
        if(!(trig&1 || trig&2 || trig&3 || trig&9)) continue;
      }
      //2016 ---
      fTrig++;

      // Q2 cut
      if(!(Q2>1)) continue;
      fQ2test++;

      // y cut
      if(!(0<yBj && yBj<1)) continue;
      fYBjtest++;

      // W cut
      if(!(5<sqrt(wBj) && sqrt(wBj)<17)) continue;
      fWBjtest++;

      // x cut
      if(!(0<xBj && xBj<1)) continue;
      fXBjtest++;

      fQ2kin.push_back(Q2);
      fXBjkin.push_back(xBj);
      fYBjkin.push_back(yBj);
      fWBjkin.push_back(sqrt(wBj));
      fNukin.push_back(nu);
      fX.push_back(x->GetLeaf("x")->GetValue());
      fY.push_back(y->GetLeaf("y")->GetValue());


      // -------------------------------------------------------------------------
      // --------- DIS event calculation -----------------------------------------
      // -------------------------------------------------------------------------

      // x binning

      if(0.004<=xBj && xBj<0.01) xbin = 0;
      else if(0.01<=xBj && xBj<0.02) xbin = 1;
      else if(0.02<=xBj && xBj<0.03) xbin = 2;
      else if(0.03<=xBj && xBj<0.04) xbin = 3;
      else if(0.04<=xBj && xBj<0.06) xbin = 4;
      else if(0.06<=xBj && xBj<0.1) xbin = 5;
      else if(0.1<=xBj && xBj<0.14) xbin = 6;
      else if(0.14<=xBj && xBj<0.18) xbin = 7;
      else xbin = 8;

      // y binning

      if(0.1<yBj && yBj<0.15) ybin = 0;
      else if(0.15<yBj && yBj<0.2) ybin = 1;
      else if(0.2<yBj && yBj<0.3) ybin = 2;
      else if(0.3<yBj && yBj<0.5) ybin = 3;
      else ybin = 4;


      // z binning

      for(int i=0; i<12; i++)
      {
#ifdef NORC
        fNDIS_evt[0][xbin][ybin][i] += 1;
        fNDIS_evt[1][xbin][ybin][i] += 1;
        fNDIS_evt[2][xbin][ybin][i] += 1;

        fNDIS_evt_err[0][xbin][ybin][i] += 1;
        fNDIS_evt_err[1][xbin][ybin][i] += 1;
        fNDIS_evt_err[2][xbin][ybin][i] += 1;
#else
        fNDIS_evt[0][xbin][ybin][i] += 1*PaAlgo::GetRadiativeWeightLiD(xBj,yBj,1,z->GetLeaf("z")->GetValue());
        fNDIS_evt[1][xbin][ybin][i] += 1*PaAlgo::GetRadiativeWeightLiD(xBj,yBj,1,z->GetLeaf("z")->GetValue());
        fNDIS_evt[2][xbin][ybin][i] += 1*PaAlgo::GetRadiativeWeightLiD(xBj,yBj,1,z->GetLeaf("z")->GetValue());

        fNDIS_evt_err[0][xbin][ybin][i] += pow(PaAlgo::GetRadiativeWeightLiD(xBj,yBj,1,z->GetLeaf("z")->GetValue()),2);
        fNDIS_evt_err[1][xbin][ybin][i] += pow(PaAlgo::GetRadiativeWeightLiD(xBj,yBj,1,z->GetLeaf("z")->GetValue()),2);
        fNDIS_evt_err[2][xbin][ybin][i] += pow(PaAlgo::GetRadiativeWeightLiD(xBj,yBj,1,z->GetLeaf("z")->GetValue()),2);
#endif

        fFlag[0][xbin][ybin][i]=0;
        fFlag[1][xbin][ybin][i]=0;
        fFlag[2][xbin][ybin][i]=0;

        // nu cut
        if(!(fNu_min[0][i]<nu && nu<fNu_max[0][i]))
        {
          fFlag[0][xbin][ybin][i]=1;
        }
        if(!(fNu_min[1][i]<nu && nu<fNu_max[1][i]))
        {
          fFlag[1][xbin][ybin][i]=1;
        }
        if(!(fNu_min[2][i]<nu && nu<fNu_max[2][i]))
        {
          fFlag[2][xbin][ybin][i]=1;
        }
        if(fFlag[0][xbin][ybin][i] /*|| fFlag[1][xbin][ybin][i] || fFlag[2][xbin][ybin][i]*/)
        {
#ifdef NORC
          fNDIS_evt[0][xbin][ybin][i] -= 1;
        // fNDIS_evt[1][xbin][ybin][i] -= 1;
        // fNDIS_evt[2][xbin][ybin][i] -= 1;
          fNDIS_evt_err[0][xbin][ybin][i] -= 1;
        // fNDIS_evt_err[1][xbin][ybin][i] -= 1;
        // fNDIS_evt_err[2][xbin][ybin][i] -= 1;
#else
          fNDIS_evt[0][xbin][ybin][i] -= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,1,z->GetLeaf("z")->GetValue());
         // fNDIS_evt[1][xbin][ybin][i] -= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,1,z->GetLeaf("z")->GetValue());
         // fNDIS_evt[2][xbin][ybin][i] -= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,1,z->GetLeaf("z")->GetValue());
          fNDIS_evt_err[0][xbin][ybin][i] -= pow(PaAlgo::GetRadiativeWeightLiD(xBj,yBj,1,z->GetLeaf("z")->GetValue()),2);
         // fNDIS_evt_err[1][xbin][ybin][i] -= pow(PaAlgo::GetRadiativeWeightLiD(xBj,yBj,1,z->GetLeaf("z")->GetValue()),2);
         // fNDIS_evt_err[2][xbin][ybin][i] -= pow(PaAlgo::GetRadiativeWeightLiD(xBj,yBj,1,z->GetLeaf("z")->GetValue()),2);
#endif
        }
      }

      // -------------------------------------------------------------------------
      // --------- Hadrons Selection ---------------------------------------------
      // -------------------------------------------------------------------------

      Pvsz pzcontainer;
      Pvsz pzcontainer_loose;
      Pvsz pzcontainer_severe;
      Pvsz pzcontainer_err;
      hadiden hadcontainer;

      for(int i=0; i<p->GetLeaf("Hadrons.P")->GetLen(); i++)
      {

        fLHsec_set.clear();
        if(fLHsec_tab) delete fLHsec_tab;

        // Sorting of LH

        if(LH->GetLeaf("Hadrons.LH")->GetValue(3+6*i)>1.8*LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i))
        {
          fLHsec_set.insert(LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i));
          fLHsec_set.insert(LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i));
          fLHsec_set.insert(LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i));
          fLHsec_set.insert(LH->GetLeaf("Hadrons.LH")->GetValue(3+6*i));
          fLHsec_set.insert(LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i));

          fLHsec_tab = new Double_t[5];
          fLHsec_tab[0] = LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i);
          fLHsec_tab[1] = LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i);
          fLHsec_tab[2] = LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i);
          fLHsec_tab[3] = LH->GetLeaf("Hadrons.LH")->GetValue(3+6*i);
          fLHsec_tab[4] = LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i);
          fusionSort(fLHsec_tab,5);
          fLHsec = fLHsec_tab[3];
        }
        else
        {
          fLHsec_set.insert(LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i));
          fLHsec_set.insert(LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i));
          fLHsec_set.insert(LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i));
          fLHsec_set.insert(LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i));

          fLHsec_tab = new Double_t[4];
          fLHsec_tab[0] = LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i);
          fLHsec_tab[1] = LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i);
          fLHsec_tab[2] = LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i);
          fLHsec_tab[3] = LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i);
          fusionSort(fLHsec_tab,4);
          fLHsec = fLHsec_tab[2];
        }

        set<Double_t>::iterator it = fLHsec_set.begin();
        advance(it, fLHsec_set.size()-2);
#ifdef DEBUG
        if(*it != fLHsec) {cout << i << " : "
        <<  LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)  << " "
        << LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)  << " "
        << LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)  << " "
        << LH->GetLeaf("Hadrons.LH")->GetValue(3+6*i)  << " "
        << LH->GetLeaf("Hadrons.LH")->GetValue(4+6*i)  << " "
        << LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)  << " "
        << endl; cout << *it << " " << fLHsec << endl;}*/
#endif

        //**********************************************************************

        // Hadron identification cuts ------------------------------------------

        // Charge +

        if(charge->GetLeaf("Hadrons.charge")->GetValue(i) == 1)
        {
          // Normal cuts ---

          if((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)>0)
             && (LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i))
             && (LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i))
             && (LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)/fLHsec>1.02)
             && (LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)>2.02)) fId = 0;

          else if((LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)>0)
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)/fLHsec>1.08)
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)>2.08)) fId = 2;

          else if((8.9<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<=17.95-5)
                  && (((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)<2.2)
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)<2.9))
                    || ((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(3+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(4+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i) == 0)))) fId = 4;

          else if((p->GetLeaf("Hadrons.P")->GetValue(i)>(17.95+5))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>0)
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)/fLHsec>1)) fId = 4;

          else if(((17.95-5)<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<(17.95+5))
                  && (((LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>0)
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)/fLHsec>1))
                  || (((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)<2.2)
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)<2.9))
                  || ((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(3+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(4+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i) == 0))))) fId = 4;
          else fId = 6;

          // Normal cuts ---


          //**********************************************************************

          // Error associated to RICH unfolding ----------------------------------

          // Loose cuts ---

          if((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)>0)
             && (LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i))
             && (LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i))
             && (LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)/fLHsec>1.00)
             && (LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)>2.00)) fId_loose = 0;

          else if((LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)>0)
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)/fLHsec>1.06)
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)>2.00)) fId_loose = 2;

          else if((8.9<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<=17.95-5)
                  && (((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)<2.3)
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)<3.0))
                    || ((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(3+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(4+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i) == 0)))) fId_loose = 4;

          else if((p->GetLeaf("Hadrons.P")->GetValue(i)>(17.95+5))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>0)
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)/fLHsec>0.98)) fId_loose = 4;

          else if(((17.95-5)<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<(17.95+5))
                  && (((LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>0)
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)/fLHsec>0.98))
                  || (((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)<2.3)
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)<3.0))
                  || ((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(3+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(4+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i) == 0))))) fId_loose = 4;
          else fId = 6;

          // Loose cuts ---

          // Severe cuts ---

          if((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)>0)
             && (LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i))
             && (LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i))
             && (LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)/fLHsec>1.06)
             && (LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)>2.04)) fId_severe = 0;

          else if((LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)>0)
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)/fLHsec>1.10)
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)>2.16)) fId_severe = 2;

          else if((8.9<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<=17.95-5)
                  && (((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)<2.0)
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)<2.7))
                    || ((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(3+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(4+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i) == 0)))) fId_severe = 4;

          else if((p->GetLeaf("Hadrons.P")->GetValue(i)>(17.95+5))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>0)
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)/fLHsec>1.06)) fId_severe = 4;

          else if(((17.95-5)<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<(17.95+5))
                  && (((LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>0)
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)/fLHsec>1.06))
                  || (((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)<2.0)
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)<2.7))
                  || ((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(3+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(4+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i) == 0))))) fId_severe = 4;
          else fId_severe = 6;

          // Severe cuts ---
        }

        // Charge -

        else if(charge->GetLeaf("Hadrons.charge")->GetValue(i) == -1)
        {
          // Normal cuts ---

          if((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)>0)
              && (LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i))
              && (LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i))
              && (LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)/fLHsec>1.02)
              && (LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)>2.02)) fId = 1;

          else if((LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)>0)
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(2))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)/fLHsec>1.08)
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)>2.08)) fId = 3;

          else if((8.9<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<=17.95-5)
                  && (((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)<2.1)
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)<2.8))
                  || ((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(3+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(4+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i) == 0)))) fId = 5;

          else if((p->GetLeaf("Hadrons.P")->GetValue(i)>17.95+5)
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>0)
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)/fLHsec>1)) fId = 5;

          else if((17.95-5<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<17.95+5)
                  && (((LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>0)
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)/fLHsec>1))
                  || (((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)<2.1)
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)<2.8))
                  || ((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(3+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(4+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i) == 0))))) fId = 5;
          else fId = 7;

          // Normal cuts ---


          //**********************************************************************

          // Error associated to RICH unfolding ----------------------------------

          // Loose cuts ---

          if((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)>0)
             && (LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i))
             && (LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i))
             && (LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)/fLHsec>1.00)
             && (LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)>2.00)) fId_loose = 1;

          else if((LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)>0)
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)/fLHsec>1.06)
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)>2.00)) fId_loose = 3;

          else if((8.9<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<=17.95-5)
                  && (((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)<2.3)
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)<3.0))
                    || ((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(3+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(4+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i) == 0)))) fId_loose = 5;

          else if((p->GetLeaf("Hadrons.P")->GetValue(i)>(17.95+5))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>0)
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)/fLHsec>0.98)) fId_loose = 5;

          else if(((17.95-5)<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<(17.95+5))
                  && (((LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>0)
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)/fLHsec>0.98))
                  || (((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)<2.3)
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)<3.0))
                  || ((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(3+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(4+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i) == 0))))) fId_loose = 5;
          else fId = 7;

          // Loose cuts ---

          // Severe cuts ---

          if((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)>0)
             && (LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i))
             && (LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i))
             && (LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)/fLHsec>1.06)
             && (LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)>2.04)) fId_severe = 1;

          else if((LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)>0)
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)/fLHsec>1.10)
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)>2.16)) fId_severe = 3;

          else if((8.9<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<=17.95-5)
                  && (((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)<2.0)
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)<2.7))
                    || ((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(3+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(4+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i) == 0)))) fId_severe = 5;

          else if((p->GetLeaf("Hadrons.P")->GetValue(i)>(17.95+5))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>0)
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)/fLHsec>1.06)) fId_severe = 5;

          else if(((17.95-5)<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<(17.95+5))
                  && (((LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>0)
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i))
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)/fLHsec>1.06))
                  || (((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)<2.0)
                  && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)<2.7))
                  || ((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(3+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(4+6*i) == 0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i) == 0))))) fId_severe = 5;
          else fId_severe = 7;

          // Severe cuts ---
        }



        //**********************************************************************

        // z calculation
        if(nu)
        {
          if(fId == 2 || fId == 3)
            zBj = sqrt(pow(p->GetLeaf("Hadrons.P")->GetValue(i),2)+pow(fM_K,2))/nu;
          else if(fId == 4 || fId == 5)
            zBj = sqrt(pow(p->GetLeaf("Hadrons.P")->GetValue(i),2)+pow(fM_p,2))/nu;
          else
            zBj = sqrt(pow(p->GetLeaf("Hadrons.P")->GetValue(i),2)+pow(fM_pi,2))/nu;
        }
        else
        {
          zBj = 0;
        }

        if(0.1<zBj) fKinematicsRD[6]->Fill(abs(ph->GetLeaf("Hadrons.ph")->GetValue(i)));

        // Maximum radiation length cumulated
        if(!(hXX0->GetLeaf("Hadrons.XX0")->GetValue(i) < 15)) continue;
        fXX0test++;

        // Momentum cut (12 GeV to 40 GeV, increasing to 3 GeV to 40 GeV)
        if(!(MOMENTUM<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<40)) continue;
        fMom++;

        // Theta cut
        if(!(0.01<thRICH->GetLeaf("Hadrons.thRICH")->GetValue(i) && thRICH->GetLeaf("Hadrons.thRICH")->GetValue(i)<0.12)) continue;
        fTRICH++;

        // RICH position cut
        if(!(pow(RICHx->GetLeaf("Hadrons.RICHx")->GetValue(i),2)+pow(RICHy->GetLeaf("Hadrons.RICHy")->GetValue(i),2)>25)) continue;
        fPosRICH++;

        // Non null charge
        if(!charge->GetLeaf("Hadrons.charge")->GetValue(i)) continue;

        Int_t theta_bin, mom_bin;
        TMatrixD res_vect(3,1);
        Double_t res_vect_err[3];
        Double_t hadron_nb;

        // Theta and momentum binning

        if(0.01<thRICH->GetLeaf("Hadrons.thRICH")->GetValue(i) && thRICH->GetLeaf("Hadrons.thRICH")->GetValue(i)<0.04)
        {
          theta_bin = 0;
          if(MOMENTUM<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<13) mom_bin = 0; // Here from 3 to 13 GeV
          else if(13<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<15) mom_bin = 1;
          else if(15<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<17) mom_bin = 2;
          else if(17<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<19) mom_bin = 3;
          else if(19<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<22) mom_bin = 4;
          else if(22<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<25) mom_bin = 5;
          else if(25<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<27) mom_bin = 6;
          else if(27<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<30) mom_bin = 7;
          else if(30<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<35) mom_bin = 8;
          else mom_bin = 9;
        }
        else
        {
          theta_bin = 1;
          if(MOMENTUM<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<13) mom_bin = 0; // Here from 3 to 13 GeV
          else if(13<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<15) mom_bin = 1;
          else if(15<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<17) mom_bin = 2;
          else if(17<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<19) mom_bin = 3;
          else if(19<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<22) mom_bin = 4;
          else if(22<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<25) mom_bin = 5;
          else if(25<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<27) mom_bin = 6;
          else if(27<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<30) mom_bin = 7;
          else if(30<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<35) mom_bin = 8;
          else mom_bin = 9;
        }

        // z cut
        if(!(0.2<zBj && zBj<0.85)) continue;

        if(0.2<zBj && zBj<0.25) zbin = 0;
        else if(0.25<zBj && zBj<0.30) zbin = 1;
        else if(0.30<zBj && zBj<0.35) zbin = 2;
        else if(0.35<zBj && zBj<0.40) zbin = 3;
        else if(0.40<zBj && zBj<0.45) zbin = 4;
        else if(0.45<zBj && zBj<0.50) zbin = 5;
        else if(0.50<zBj && zBj<0.55) zbin = 6;
        else if(0.55<zBj && zBj<0.60) zbin = 7;
        else if(0.60<zBj && zBj<0.65) zbin = 8;
        else if(0.65<zBj && zBj<0.70) zbin = 9;
        else if(0.70<zBj && zBj<0.75) zbin = 10;
        else zbin = 11;

        //**********************************************************************

        // Save of hadrons

        int hadron_flag = 0;

        if(fId==0)
        {
          if(!fFlag[0][xbin][ybin][zbin])
          {
            fHplus++; fPiplus++;
            res_vect = inv_rich_p[theta_bin][mom_bin]*pi_vect;
            for(int rce=0; rce<3; rce++) res_vect_err[rce] = pi_unfolding_err_p[theta_bin][mom_bin][rce];
            hadron_nb = 1;
#ifdef NORC
#else
            res_vect[0][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect[1][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect[2][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect_err[0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect_err[1] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect_err[2] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            hadron_nb *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
#endif
            fPiplus_true += res_vect[0][0]; fKplus_true += res_vect[1][0]; fPplus_true += res_vect[2][0];
            fPiplus_err += pow(res_vect_err[0],2);
            fKplus_err += pow(res_vect_err[1],2);
            fPplus_err += pow(res_vect_err[2],2);

            pzcontainer.vec[1][0].push_back(zBj);
            pzcontainer.vec[1][1].push_back(res_vect[0][0]);
            pzcontainer.vec[1][2].push_back(res_vect[1][0]);
            pzcontainer.vec[1][3].push_back(res_vect[2][0]);

            pzcontainer_err.vec[1][0].push_back(zBj);
            pzcontainer_err.vec[1][1].push_back(pow(res_vect_err[0],2));
            pzcontainer_err.vec[1][2].push_back(pow(res_vect_err[1],2));
            pzcontainer_err.vec[1][3].push_back(pow(res_vect_err[2],2));

            pzcontainer.vec[1][4].push_back(hadron_nb);
            pzcontainer_err.vec[1][4].push_back(hadron_nb);

            hadcontainer.vec.push_back(0);
#ifdef DEBUG
            cout << res_vect[0][0] << " " << res_vect[1][0] << " " << res_vect[2][0] << endl;
#endif
          }
        }
        else if(fId==1)
        {
          if(!fFlag[0][xbin][ybin][zbin])
          {
            fHminus++; fPiminus++;
            res_vect = inv_rich_m[theta_bin][mom_bin]*pi_vect;
            for(int rce=0; rce<3; rce++) res_vect_err[rce] = pi_unfolding_err_m[theta_bin][mom_bin][rce];
            hadron_nb = 1;
#ifdef NORC
#else
            res_vect[0][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect[1][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect[2][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect_err[0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect_err[1] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect_err[2] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            hadron_nb *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
#endif
            fPiminus_true += res_vect[0][0]; fKminus_true += res_vect[1][0]; fPminus_true += res_vect[2][0];
            fPiminus_err += pow(res_vect_err[0],2);
            fKminus_err += pow(res_vect_err[1],2);
            fPminus_err += pow(res_vect_err[2],2);

            pzcontainer.vec[0][0].push_back(zBj);
            pzcontainer.vec[0][1].push_back(res_vect[0][0]);
            pzcontainer.vec[0][2].push_back(res_vect[1][0]);
            pzcontainer.vec[0][3].push_back(res_vect[2][0]);

            pzcontainer_err.vec[0][0].push_back(zBj);
            pzcontainer_err.vec[0][1].push_back(pow(res_vect_err[0],2));
            pzcontainer_err.vec[0][2].push_back(pow(res_vect_err[1],2));
            pzcontainer_err.vec[0][3].push_back(pow(res_vect_err[2],2));

            pzcontainer.vec[0][4].push_back(hadron_nb);
            pzcontainer_err.vec[0][4].push_back(hadron_nb);

            hadcontainer.vec.push_back(1);
#ifdef DEBUG
            cout << res_vect[0][0] << " " << res_vect[1][0] << " " << res_vect[2][0] << endl;
#endif
          }
        }
        else if(fId==2)
        {
          if(!fFlag[0][xbin][ybin][zbin])
          {
            fHplus++;
#ifdef NORC
            pzcontainer.vec[1][4].push_back(1);
            pzcontainer_err.vec[1][4].push_back(1);
#else
            pzcontainer.vec[1][4].push_back(1*PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue()));
            pzcontainer_err.vec[1][4].push_back(pow(PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue()),2));
#endif
            hadron_flag = 1;
          }
          if(!fFlag[1][xbin][ybin][zbin])
          {
            fKplus++;
            res_vect = inv_rich_p[theta_bin][mom_bin]*k_vect;
            for(int rce=0; rce<3; rce++) res_vect_err[rce] = k_unfolding_err_p[theta_bin][mom_bin][rce];
            hadron_nb = 1;
#ifdef NORC
#else
            res_vect[0][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect[1][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect[2][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect_err[0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect_err[1] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect_err[2] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            hadron_nb *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
#endif
            fPiplus_true += res_vect[0][0]; fKplus_true += res_vect[1][0]; fPplus_true += res_vect[2][0];
            fPiplus_err += pow(res_vect_err[0],2);
            fKplus_err += pow(res_vect_err[1],2);
            fPplus_err += pow(res_vect_err[2],2);

            pzcontainer.vec[1][0].push_back(zBj);
            pzcontainer.vec[1][1].push_back(res_vect[0][0]);
            pzcontainer.vec[1][2].push_back(res_vect[1][0]);
            pzcontainer.vec[1][3].push_back(res_vect[2][0]);

            pzcontainer_err.vec[1][0].push_back(zBj);
            pzcontainer_err.vec[1][1].push_back(pow(res_vect_err[0],2));
            pzcontainer_err.vec[1][2].push_back(pow(res_vect_err[1],2));
            pzcontainer_err.vec[1][3].push_back(pow(res_vect_err[2],2));

            if(!hadron_flag)
            {
              pzcontainer.vec[1][4].push_back(0);
              pzcontainer_err.vec[1][4].push_back(0);
            }

            hadcontainer.vec.push_back(2);
#ifdef DEBUG
            cout << res_vect[0][0] << " " << res_vect[1][0] << " " << res_vect[2][0] << endl;
#endif
          }
        }
        else if(fId==3)
        {
          if(!fFlag[0][xbin][ybin][zbin])
          {
            fHminus++;
#ifdef NORC
            pzcontainer.vec[0][4].push_back(1);
            pzcontainer_err.vec[0][4].push_back(1);
#else
            pzcontainer.vec[0][4].push_back(1*PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue()));
            pzcontainer_err.vec[0][4].push_back(pow(PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue()),2));
#endif
            hadron_flag = 1;
          }
          if(!fFlag[1][xbin][ybin][zbin])
          {
            fKminus++;
            res_vect = inv_rich_m[theta_bin][mom_bin]*k_vect;
            for(int rce=0; rce<3; rce++) res_vect_err[rce] = k_unfolding_err_m[theta_bin][mom_bin][rce];
            hadron_nb = 1;
#ifdef NORC
#else
            res_vect[0][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect[1][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect[2][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect_err[0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect_err[1] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect_err[2] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            hadron_nb *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
#endif
            fPiminus_true += res_vect[0][0]; fKminus_true += res_vect[1][0]; fPminus_true += res_vect[2][0];
            fPiminus_err += pow(res_vect_err[0],2);
            fKminus_err += pow(res_vect_err[1],2);
            fPminus_err += pow(res_vect_err[2],2);

            pzcontainer.vec[0][0].push_back(zBj);
            pzcontainer.vec[0][1].push_back(res_vect[0][0]);
            pzcontainer.vec[0][2].push_back(res_vect[1][0]);
            pzcontainer.vec[0][3].push_back(res_vect[2][0]);

            pzcontainer_err.vec[0][0].push_back(zBj);
            pzcontainer_err.vec[0][1].push_back(pow(res_vect_err[0],2));
            pzcontainer_err.vec[0][2].push_back(pow(res_vect_err[1],2));
            pzcontainer_err.vec[0][3].push_back(pow(res_vect_err[2],2));

            if(!hadron_flag)
            {
              pzcontainer.vec[0][4].push_back(0);
              pzcontainer_err.vec[0][4].push_back(0);
            }

            hadcontainer.vec.push_back(3);
#ifdef DEBUG
            cout << res_vect[0][0] << " " << res_vect[1][0] << " " << res_vect[2][0] << endl;
#endif
          }
        }
        else if(fId==4)
        {
          if(!fFlag[0][xbin][ybin][zbin])
          {
            fHplus++;
#ifdef NORC
            pzcontainer.vec[1][4].push_back(1);
            pzcontainer_err.vec[1][4].push_back(1);
#else
            pzcontainer.vec[1][4].push_back(1*PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue()));
            pzcontainer_err.vec[1][4].push_back(pow(PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue()),2));
#endif
            hadron_flag = 1;
          }
          if(!fFlag[2][xbin][ybin][zbin])
          {
            fPplus++;
            res_vect = inv_rich_p[theta_bin][mom_bin]*p_vect;
            for(int rce=0; rce<3; rce++) res_vect_err[rce] = p_unfolding_err_p[theta_bin][mom_bin][rce];
            hadron_nb = 1;
#ifdef NORC
#else
            res_vect[0][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect[1][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect[2][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect_err[0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect_err[1] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect_err[2] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            hadron_nb *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
#endif
            fPiplus_true += res_vect[0][0]; fKplus_true += res_vect[1][0]; fPplus_true += res_vect[2][0];
            fPiplus_err += pow(res_vect_err[0],2);
            fKplus_err += pow(res_vect_err[1],2);
            fPplus_err += pow(res_vect_err[2],2);

            pzcontainer.vec[1][0].push_back(zBj);
            pzcontainer.vec[1][1].push_back(res_vect[0][0]);
            pzcontainer.vec[1][2].push_back(res_vect[1][0]);
            pzcontainer.vec[1][3].push_back(res_vect[2][0]);

            pzcontainer_err.vec[1][0].push_back(zBj);
            pzcontainer_err.vec[1][1].push_back(pow(res_vect_err[0],2));
            pzcontainer_err.vec[1][2].push_back(pow(res_vect_err[1],2));
            pzcontainer_err.vec[1][3].push_back(pow(res_vect_err[2],2));

            if(!hadron_flag)
            {
              pzcontainer.vec[1][4].push_back(0);
              pzcontainer_err.vec[1][4].push_back(0);
            }

            hadcontainer.vec.push_back(0);
#ifdef DEBUG
            cout << res_vect[0][0] << " " << res_vect[1][0] << " " << res_vect[2][0] << endl;
#endif
          }
        }
        else if(fId==5)
        {
          if(!fFlag[0][xbin][ybin][zbin])
          {
            fHminus++;
#ifdef NORC
            pzcontainer.vec[0][4].push_back(1);
            pzcontainer_err.vec[0][4].push_back(1);
#else
            pzcontainer.vec[0][4].push_back(1*PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue()));
            pzcontainer_err.vec[0][4].push_back(pow(PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue()),2));
#endif
            hadron_flag = 1;
          }
          if(!fFlag[2][xbin][ybin][zbin])
          {
            fPminus++;
            res_vect = inv_rich_m[theta_bin][mom_bin]*p_vect;
            for(int rce=0; rce<3; rce++) res_vect_err[rce] = p_unfolding_err_m[theta_bin][mom_bin][rce];
            hadron_nb = 1;
#ifdef NORC
#else
            res_vect[0][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect[1][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect[2][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect_err[0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect_err[1] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect_err[2] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            hadron_nb *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
#endif
            fPiminus_true += res_vect[0][0]; fKminus_true += res_vect[1][0]; fPminus_true += res_vect[2][0];
            fPiminus_err += pow(res_vect_err[0],2);
            fKminus_err += pow(res_vect_err[1],2);
            fPminus_err += pow(res_vect_err[2],2);

            pzcontainer.vec[0][0].push_back(zBj);
            pzcontainer.vec[0][1].push_back(res_vect[0][0]);
            pzcontainer.vec[0][2].push_back(res_vect[1][0]);
            pzcontainer.vec[0][3].push_back(res_vect[2][0]);

            pzcontainer_err.vec[0][0].push_back(zBj);
            pzcontainer_err.vec[0][1].push_back(pow(res_vect_err[0],2));
            pzcontainer_err.vec[0][2].push_back(pow(res_vect_err[1],2));
            pzcontainer_err.vec[0][3].push_back(pow(res_vect_err[2],2));

            if(!hadron_flag)
            {
              pzcontainer.vec[0][4].push_back(0);
              pzcontainer_err.vec[0][4].push_back(0);
            }

            hadcontainer.vec.push_back(5);
#ifdef DEBUG
            cout << res_vect[0][0] << " " << res_vect[1][0] << " " << res_vect[2][0] << endl;
#endif
          }
        }
        else if(fId==6)
        {
          if(!fFlag[0][xbin][ybin][zbin])
          {
            fHplus++;
            pzcontainer.vec[1][0].push_back(zBj); pzcontainer_err.vec[1][0].push_back(zBj);
            pzcontainer.vec[1][1].push_back(0); pzcontainer_err.vec[1][1].push_back(0);
            pzcontainer.vec[1][2].push_back(0); pzcontainer_err.vec[1][2].push_back(0);
            pzcontainer.vec[1][3].push_back(0); pzcontainer_err.vec[1][3].push_back(0);
#ifdef NORC
            pzcontainer.vec[1][4].push_back(1);
            pzcontainer_err.vec[0][4].push_back(1);
#else
            pzcontainer.vec[1][4].push_back(1*PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue()));
            pzcontainer_err.vec[0][4].push_back(pow(PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue()),2));
#endif
            hadcontainer.vec.push_back(6);
          }
        }
        else if(fId==7)
        {
          if(!fFlag[0][xbin][ybin][zbin])
          {
            fHminus++;
            pzcontainer.vec[0][0].push_back(zBj); pzcontainer_err.vec[0][0].push_back(zBj);
            pzcontainer.vec[0][1].push_back(0); pzcontainer_err.vec[0][1].push_back(0);
            pzcontainer.vec[0][2].push_back(0); pzcontainer_err.vec[0][2].push_back(0);
            pzcontainer.vec[0][3].push_back(0); pzcontainer_err.vec[0][3].push_back(0);
#ifdef NORC
            pzcontainer.vec[0][4].push_back(1);
            pzcontainer_err.vec[0][4].push_back(1);
#else
            pzcontainer.vec[0][4].push_back(1*PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue()));
            pzcontainer_err.vec[0][4].push_back(pow(PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue()),2));
#endif
            hadcontainer.vec.push_back(7);
          }
        }
        else
        {

        }


        //**********************************************************************

        // Loose cut

        hadron_flag = 0;

        if(fId_loose==0)
        {
          if(!fFlag[0][xbin][ybin][zbin])
          {
            res_vect = inv_rich_p[theta_bin][mom_bin]*pi_vect;
            hadron_nb = 1;
#ifdef NORC
#else
            res_vect[0][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect[1][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect[2][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            hadron_nb *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
#endif
            pzcontainer_loose.vec[1][0].push_back(zBj);
            pzcontainer_loose.vec[1][1].push_back(res_vect[0][0]);
            pzcontainer_loose.vec[1][2].push_back(res_vect[1][0]);
            pzcontainer_loose.vec[1][3].push_back(res_vect[2][0]);
            pzcontainer_loose.vec[1][4].push_back(hadron_nb);
          }
        }
        else if(fId_loose==1)
        {
          if(!fFlag[0][xbin][ybin][zbin])
          {
            res_vect = inv_rich_m[theta_bin][mom_bin]*pi_vect;
            hadron_nb = 1;
#ifdef NORC
#else
            res_vect[0][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect[1][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect[2][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            hadron_nb *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
#endif
            pzcontainer_loose.vec[0][0].push_back(zBj);
            pzcontainer_loose.vec[0][1].push_back(res_vect[0][0]);
            pzcontainer_loose.vec[0][2].push_back(res_vect[1][0]);
            pzcontainer_loose.vec[0][3].push_back(res_vect[2][0]);
            pzcontainer_loose.vec[0][4].push_back(hadron_nb);
          }
        }
        else if(fId_loose==2)
        {
          if(!fFlag[0][xbin][ybin][zbin])
          {
#ifdef NORC
            pzcontainer_loose.vec[1][4].push_back(1);
#else
            pzcontainer_loose.vec[1][4].push_back(1*PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue()));
#endif
            hadron_flag = 1;
          }
          if(!fFlag[1][xbin][ybin][zbin])
          {
            res_vect = inv_rich_p[theta_bin][mom_bin]*k_vect;
#ifdef NORC
#else
            res_vect[0][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect[1][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect[2][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
#endif
            pzcontainer_loose.vec[1][0].push_back(zBj);
            pzcontainer_loose.vec[1][1].push_back(res_vect[0][0]);
            pzcontainer_loose.vec[1][2].push_back(res_vect[1][0]);
            pzcontainer_loose.vec[1][3].push_back(res_vect[2][0]);
            if(!hadron_flag) pzcontainer_loose.vec[1][4].push_back(0);
          }
        }
        else if(fId_loose==3)
        {
          if(!fFlag[0][xbin][ybin][zbin])
          {
#ifdef NORC
            pzcontainer_loose.vec[0][4].push_back(1);
#else
            pzcontainer_loose.vec[0][4].push_back(1*PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue()));
#endif
            hadron_flag = 1;
          }
          if(!fFlag[1][xbin][ybin][zbin])
          {
            res_vect = inv_rich_m[theta_bin][mom_bin]*k_vect;
#ifdef NORC
#else
            res_vect[0][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect[1][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect[2][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
#endif
            pzcontainer_loose.vec[0][0].push_back(zBj);
            pzcontainer_loose.vec[0][1].push_back(res_vect[0][0]);
            pzcontainer_loose.vec[0][2].push_back(res_vect[1][0]);
            pzcontainer_loose.vec[0][3].push_back(res_vect[2][0]);
            if(!hadron_flag) pzcontainer_loose.vec[0][4].push_back(0);
          }
        }
        else if(fId_loose==4)
        {
          if(!fFlag[0][xbin][ybin][zbin])
          {
#ifdef NORC
            pzcontainer_loose.vec[1][4].push_back(1);
#else
            pzcontainer_loose.vec[1][4].push_back(1*PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue()));
#endif
            hadron_flag = 1;
          }
          if(!fFlag[2][xbin][ybin][zbin])
          {
            res_vect = inv_rich_p[theta_bin][mom_bin]*p_vect;
#ifdef NORC
#else
            res_vect[0][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect[1][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect[2][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
#endif
            pzcontainer_loose.vec[1][0].push_back(zBj);
            pzcontainer_loose.vec[1][1].push_back(res_vect[0][0]);
            pzcontainer_loose.vec[1][2].push_back(res_vect[1][0]);
            pzcontainer_loose.vec[1][3].push_back(res_vect[2][0]);
            if(!hadron_flag) pzcontainer_loose.vec[1][4].push_back(0);
          }
        }
        else if(fId_loose==5)
        {
          if(!fFlag[0][xbin][ybin][zbin])
          {
#ifdef NORC
            pzcontainer_loose.vec[0][4].push_back(1);
#else
            pzcontainer_loose.vec[0][4].push_back(1*PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue()));
#endif
            hadron_flag = 1;
          }
          if(!fFlag[2][xbin][ybin][zbin])
          {
            res_vect = inv_rich_m[theta_bin][mom_bin]*p_vect;
#ifdef NORC
#else
            res_vect[0][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect[1][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect[2][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
#endif
            pzcontainer_loose.vec[0][0].push_back(zBj);
            pzcontainer_loose.vec[0][1].push_back(res_vect[0][0]);
            pzcontainer_loose.vec[0][2].push_back(res_vect[1][0]);
            pzcontainer_loose.vec[0][3].push_back(res_vect[2][0]);
            if(!hadron_flag) pzcontainer_loose.vec[0][4].push_back(0);
          }
        }
        else if(fId_loose==6)
        {
          if(!fFlag[0][xbin][ybin][zbin])
          {
            pzcontainer_loose.vec[1][0].push_back(zBj);
            pzcontainer_loose.vec[1][1].push_back(0);
            pzcontainer_loose.vec[1][2].push_back(0);
            pzcontainer_loose.vec[1][3].push_back(0);
            hadron_nb = 1;
#ifdef NORC
#else
            hadron_nb *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
#endif
            pzcontainer_loose.vec[1][4].push_back(hadron_nb);
          }
        }
        else if(fId_loose==7)
        {
          if(!fFlag[0][xbin][ybin][zbin])
          {
            pzcontainer_loose.vec[0][0].push_back(zBj);
            pzcontainer_loose.vec[0][1].push_back(0);
            pzcontainer_loose.vec[0][2].push_back(0);
            pzcontainer_loose.vec[0][3].push_back(0);
            hadron_nb = 1;
#ifdef NORC
#else
            hadron_nb *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
#endif
            pzcontainer_loose.vec[0][4].push_back(hadron_nb);
          }
        }
        else
        {

        }


        // Severe cut

        hadron_flag = 0;

        if(fId_severe==0)
        {
          if(!fFlag[0][xbin][ybin][zbin])
          {
            res_vect = inv_rich_p[theta_bin][mom_bin]*pi_vect;
            hadron_nb = 1;
#ifdef NORC
#else
            res_vect[0][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect[1][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect[2][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            hadron_nb *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
#endif
            pzcontainer_severe.vec[1][0].push_back(zBj);
            pzcontainer_severe.vec[1][1].push_back(res_vect[0][0]);
            pzcontainer_severe.vec[1][2].push_back(res_vect[1][0]);
            pzcontainer_severe.vec[1][3].push_back(res_vect[2][0]);
            pzcontainer_severe.vec[1][4].push_back(hadron_nb);
          }
        }
        else if(fId_severe==1)
        {
          if(!fFlag[0][xbin][ybin][zbin])
          {
            res_vect = inv_rich_m[theta_bin][mom_bin]*pi_vect;
            hadron_nb = 1;
#ifdef NORC
#else
            res_vect[0][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect[1][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect[2][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            hadron_nb *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
#endif
            pzcontainer_severe.vec[0][0].push_back(zBj);
            pzcontainer_severe.vec[0][1].push_back(res_vect[0][0]);
            pzcontainer_severe.vec[0][2].push_back(res_vect[1][0]);
            pzcontainer_severe.vec[0][3].push_back(res_vect[2][0]);
            pzcontainer_severe.vec[0][4].push_back(hadron_nb);
          }
        }
        else if(fId_severe==2)
        {
          if(!fFlag[0][xbin][ybin][zbin])
          {
#ifdef NORC
            pzcontainer_severe.vec[1][4].push_back(1);
#else
            pzcontainer_severe.vec[1][4].push_back(1*PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue()));
#endif
            hadron_flag = 1;
          }
          if(!fFlag[1][xbin][ybin][zbin])
          {
            res_vect = inv_rich_p[theta_bin][mom_bin]*k_vect;
#ifdef NORC
#else
            res_vect[0][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect[1][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect[2][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
#endif
            pzcontainer_severe.vec[1][0].push_back(zBj);
            pzcontainer_severe.vec[1][1].push_back(res_vect[0][0]);
            pzcontainer_severe.vec[1][2].push_back(res_vect[1][0]);
            pzcontainer_severe.vec[1][3].push_back(res_vect[2][0]);
            if(!hadron_flag) pzcontainer_severe.vec[1][4].push_back(0);
          }
        }
        else if(fId_severe==3)
        {
          if(!fFlag[0][xbin][ybin][zbin])
          {
#ifdef NORC
            pzcontainer_severe.vec[0][4].push_back(1);
#else
            pzcontainer_severe.vec[0][4].push_back(1*PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue()));
#endif
            hadron_flag = 1;
          }
          if(!fFlag[1][xbin][ybin][zbin])
          {
            res_vect = inv_rich_m[theta_bin][mom_bin]*k_vect;
#ifdef NORC
#else
            res_vect[0][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect[1][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect[2][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
#endif
            pzcontainer_severe.vec[0][0].push_back(zBj);
            pzcontainer_severe.vec[0][1].push_back(res_vect[0][0]);
            pzcontainer_severe.vec[0][2].push_back(res_vect[1][0]);
            pzcontainer_severe.vec[0][3].push_back(res_vect[2][0]);
            if(!hadron_flag) pzcontainer_severe.vec[0][4].push_back(0);
          }
        }
        else if(fId_severe==4)
        {
          if(!fFlag[0][xbin][ybin][zbin])
          {
#ifdef NORC
            pzcontainer_severe.vec[1][4].push_back(1);
#else
            pzcontainer_severe.vec[1][4].push_back(1*PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue()));
#endif
            hadron_flag = 1;
          }
          if(!fFlag[2][xbin][ybin][zbin])
          {
            res_vect = inv_rich_p[theta_bin][mom_bin]*p_vect;
#ifdef NORC
#else
            res_vect[0][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect[1][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect[2][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
#endif
            pzcontainer_severe.vec[1][0].push_back(zBj);
            pzcontainer_severe.vec[1][1].push_back(res_vect[0][0]);
            pzcontainer_severe.vec[1][2].push_back(res_vect[1][0]);
            pzcontainer_severe.vec[1][3].push_back(res_vect[2][0]);
            if(!hadron_flag) pzcontainer_severe.vec[1][4].push_back(0);
          }
        }
        else if(fId_severe==5)
        {
          if(!fFlag[0][xbin][ybin][zbin])
          {
#ifdef NORC
            pzcontainer_severe.vec[0][4].push_back(1);
#else
            pzcontainer_severe.vec[0][4].push_back(1*PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue()));
#endif
            hadron_flag = 1;
          }
          if(!fFlag[2][xbin][ybin][zbin])
          {
            res_vect = inv_rich_m[theta_bin][mom_bin]*p_vect;
#ifdef NORC
#else
            res_vect[0][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect[1][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
            res_vect[2][0] *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
#endif
            pzcontainer_severe.vec[0][0].push_back(zBj);
            pzcontainer_severe.vec[0][1].push_back(res_vect[0][0]);
            pzcontainer_severe.vec[0][2].push_back(res_vect[1][0]);
            pzcontainer_severe.vec[0][3].push_back(res_vect[2][0]);
            if(!hadron_flag) pzcontainer_severe.vec[0][4].push_back(0);
          }
        }
        else if(fId_severe==6)
        {
          if(!fFlag[0][xbin][ybin][zbin])
          {
            pzcontainer_severe.vec[1][0].push_back(zBj);
            pzcontainer_severe.vec[1][1].push_back(0);
            pzcontainer_severe.vec[1][2].push_back(0);
            pzcontainer_severe.vec[1][3].push_back(0);
            hadron_nb = 1;
#ifdef NORC
#else
            hadron_nb *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
#endif
            pzcontainer_severe.vec[1][4].push_back(hadron_nb);
          }
        }
        else if(fId_severe==7)
        {
          if(!fFlag[0][xbin][ybin][zbin])
          {
            pzcontainer_severe.vec[0][0].push_back(zBj);
            pzcontainer_severe.vec[0][1].push_back(0);
            pzcontainer_severe.vec[0][2].push_back(0);
            pzcontainer_severe.vec[0][3].push_back(0);
            hadron_nb = 1;
#ifdef NORC
#else
            hadron_nb *= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,2,z->GetLeaf("z")->GetValue());
#endif
            pzcontainer_severe.vec[0][4].push_back(hadron_nb);
          }
        }
        else
        {

        }

      }

      //Misc
      fQ2.push_back(Q2);
      fXBj.push_back(xBj);
      fYBj.push_back(yBj);
      fWBj.push_back(wBj);
      fNu.push_back(nu);

      fPvsz.push_back(pzcontainer);
      fPvsz_err.push_back(pzcontainer_err);
      fHadiden.push_back(hadcontainer);

      Q2local.push_back(Q2);
      Pvszlocal.push_back(pzcontainer);
      Pvsz_errlocal.push_back(pzcontainer_err);
      XBjlocal.push_back(xBj);
      YBjlocal.push_back(yBj);

      Q2loose.push_back(Q2);
      Pvszloose.push_back(pzcontainer_loose);
      XBjloose.push_back(xBj);
      YBjloose.push_back(yBj);

      Q2severe.push_back(Q2);
      Pvszsevere.push_back(pzcontainer_severe);
      XBjsevere.push_back(xBj);
      YBjsevere.push_back(yBj);

    }

    cout << "\n" << endl;

    // Loose cut

    for(int i=0; i<int(Q2loose.size()); i++)
    {
      if(0.004<=XBjloose[i] && XBjloose[i]<0.01) xbin = 0;
      else if(0.01<=XBjloose[i] && XBjloose[i]<0.02) xbin = 1;
      else if(0.02<=XBjloose[i] && XBjloose[i]<0.03) xbin = 2;
      else if(0.03<=XBjloose[i] && XBjloose[i]<0.04) xbin = 3;
      else if(0.04<=XBjloose[i] && XBjloose[i]<0.06) xbin = 4;
      else if(0.06<=XBjloose[i] && XBjloose[i]<0.1) xbin = 5;
      else if(0.1<=XBjloose[i] && XBjloose[i]<0.14) xbin = 6;
      else if(0.1<=XBjloose[i] && XBjloose[i]<0.18) xbin = 7;
      else xbin = 8;

      if(0.1<YBjloose[i] && YBjloose[i]<0.15) ybin = 0;
      else if(0.15<YBjloose[i] && YBjloose[i]<0.2) ybin = 1;
      else if(0.2<YBjloose[i] && YBjloose[i]<0.3) ybin = 2;
      else if(0.3<YBjloose[i] && YBjloose[i]<0.5) ybin = 3;
      else ybin = 4;

      for(int j=0; j<2; j++)
      {
        for(int l=0; l<int(Pvszloose[i].vec[j][0].size()); l++)
        {
          if(0.2<Pvszloose[i].vec[j][0][l] && Pvszloose[i].vec[j][0][l]<0.25) zbin = 0;
          else if(0.25<Pvszloose[i].vec[j][0][l] && Pvszloose[i].vec[j][0][l]<0.30) zbin = 1;
          else if(0.30<Pvszloose[i].vec[j][0][l] && Pvszloose[i].vec[j][0][l]<0.35) zbin = 2;
          else if(0.35<Pvszloose[i].vec[j][0][l] && Pvszloose[i].vec[j][0][l]<0.40) zbin = 3;
          else if(0.40<Pvszloose[i].vec[j][0][l] && Pvszloose[i].vec[j][0][l]<0.45) zbin = 4;
          else if(0.45<Pvszloose[i].vec[j][0][l] && Pvszloose[i].vec[j][0][l]<0.50) zbin = 5;
          else if(0.50<Pvszloose[i].vec[j][0][l] && Pvszloose[i].vec[j][0][l]<0.55) zbin = 6;
          else if(0.55<Pvszloose[i].vec[j][0][l] && Pvszloose[i].vec[j][0][l]<0.60) zbin = 7;
          else if(0.60<Pvszloose[i].vec[j][0][l] && Pvszloose[i].vec[j][0][l]<0.65) zbin = 8;
          else if(0.65<Pvszloose[i].vec[j][0][l] && Pvszloose[i].vec[j][0][l]<0.70) zbin = 9;
          else if(0.70<Pvszloose[i].vec[j][0][l] && Pvszloose[i].vec[j][0][l]<0.75) zbin = 10;
          else zbin = 11;

          for(int ll=0; ll<4; ll++)
          {
            fBinning_loose[xbin][ybin][zbin].tab[j][0][ll] += Pvszloose[i].vec[j][ll+1][l];
          }
        }
      }

    }

    // Severe cut

    for(int i=0; i<int(Q2severe.size()); i++)
    {
      if(0.004<=XBjsevere[i] && XBjsevere[i]<0.01) xbin = 0;
      else if(0.01<=XBjsevere[i] && XBjsevere[i]<0.02) xbin = 1;
      else if(0.02<=XBjsevere[i] && XBjsevere[i]<0.03) xbin = 2;
      else if(0.03<=XBjsevere[i] && XBjsevere[i]<0.04) xbin = 3;
      else if(0.04<=XBjsevere[i] && XBjsevere[i]<0.06) xbin = 4;
      else if(0.06<=XBjsevere[i] && XBjsevere[i]<0.1) xbin = 5;
      else if(0.1<=XBjsevere[i] && XBjsevere[i]<0.14) xbin = 6;
      else if(0.1<=XBjsevere[i] && XBjsevere[i]<0.18) xbin = 7;
      else xbin = 8;

      if(0.1<YBjsevere[i] && YBjsevere[i]<0.15) ybin = 0;
      else if(0.15<YBjsevere[i] && YBjsevere[i]<0.2) ybin = 1;
      else if(0.2<YBjsevere[i] && YBjsevere[i]<0.3) ybin = 2;
      else if(0.3<YBjsevere[i] && YBjsevere[i]<0.5) ybin = 3;
      else ybin = 4;

      for(int j=0; j<2; j++)
      {
        for(int l=0; l<int(Pvszsevere[i].vec[j][0].size()); l++)
        {
          if(0.2<Pvszsevere[i].vec[j][0][l] && Pvszsevere[i].vec[j][0][l]<0.25) zbin = 0;
          else if(0.25<Pvszsevere[i].vec[j][0][l] && Pvszsevere[i].vec[j][0][l]<0.30) zbin = 1;
          else if(0.30<Pvszsevere[i].vec[j][0][l] && Pvszsevere[i].vec[j][0][l]<0.35) zbin = 2;
          else if(0.35<Pvszsevere[i].vec[j][0][l] && Pvszsevere[i].vec[j][0][l]<0.40) zbin = 3;
          else if(0.40<Pvszsevere[i].vec[j][0][l] && Pvszsevere[i].vec[j][0][l]<0.45) zbin = 4;
          else if(0.45<Pvszsevere[i].vec[j][0][l] && Pvszsevere[i].vec[j][0][l]<0.50) zbin = 5;
          else if(0.50<Pvszsevere[i].vec[j][0][l] && Pvszsevere[i].vec[j][0][l]<0.55) zbin = 6;
          else if(0.55<Pvszsevere[i].vec[j][0][l] && Pvszsevere[i].vec[j][0][l]<0.60) zbin = 7;
          else if(0.60<Pvszsevere[i].vec[j][0][l] && Pvszsevere[i].vec[j][0][l]<0.65) zbin = 8;
          else if(0.65<Pvszsevere[i].vec[j][0][l] && Pvszsevere[i].vec[j][0][l]<0.70) zbin = 9;
          else if(0.70<Pvszsevere[i].vec[j][0][l] && Pvszsevere[i].vec[j][0][l]<0.75) zbin = 10;
          else zbin = 11;

          for(int ll=0; ll<4; ll++)
          {
            fBinning_severe[xbin][ybin][zbin].tab[j][0][ll] += Pvszsevere[i].vec[j][ll+1][l];
          }
        }
      }

    }

    for(int i=0; i<int(Q2local.size()); i++)
    {
      if(0.004<=XBjlocal[i] && XBjlocal[i]<0.01) xbin = 0;
      else if(0.01<=XBjlocal[i] && XBjlocal[i]<0.02) xbin = 1;
      else if(0.02<=XBjlocal[i] && XBjlocal[i]<0.03) xbin = 2;
      else if(0.03<=XBjlocal[i] && XBjlocal[i]<0.04) xbin = 3;
      else if(0.04<=XBjlocal[i] && XBjlocal[i]<0.06) xbin = 4;
      else if(0.06<=XBjlocal[i] && XBjlocal[i]<0.1) xbin = 5;
      else if(0.1<=XBjlocal[i] && XBjlocal[i]<0.14) xbin = 6;
      else if(0.1<=XBjlocal[i] && XBjlocal[i]<0.18) xbin = 7;
      else xbin = 8;

      if(0.1<YBjlocal[i] && YBjlocal[i]<0.15) ybin = 0;
      else if(0.15<YBjlocal[i] && YBjlocal[i]<0.2) ybin = 1;
      else if(0.2<YBjlocal[i] && YBjlocal[i]<0.3) ybin = 2;
      else if(0.3<YBjlocal[i] && YBjlocal[i]<0.5) ybin = 3;
      else ybin = 4;

      for(int j=0; j<2; j++)
      {
        for(int l=0; l<int(Pvszlocal[i].vec[j][0].size()); l++)
        {
          if(0.2<Pvszlocal[i].vec[j][0][l] && Pvszlocal[i].vec[j][0][l]<0.25) zbin = 0;
          else if(0.25<Pvszlocal[i].vec[j][0][l] && Pvszlocal[i].vec[j][0][l]<0.30) zbin = 1;
          else if(0.30<Pvszlocal[i].vec[j][0][l] && Pvszlocal[i].vec[j][0][l]<0.35) zbin = 2;
          else if(0.35<Pvszlocal[i].vec[j][0][l] && Pvszlocal[i].vec[j][0][l]<0.40) zbin = 3;
          else if(0.40<Pvszlocal[i].vec[j][0][l] && Pvszlocal[i].vec[j][0][l]<0.45) zbin = 4;
          else if(0.45<Pvszlocal[i].vec[j][0][l] && Pvszlocal[i].vec[j][0][l]<0.50) zbin = 5;
          else if(0.50<Pvszlocal[i].vec[j][0][l] && Pvszlocal[i].vec[j][0][l]<0.55) zbin = 6;
          else if(0.55<Pvszlocal[i].vec[j][0][l] && Pvszlocal[i].vec[j][0][l]<0.60) zbin = 7;
          else if(0.60<Pvszlocal[i].vec[j][0][l] && Pvszlocal[i].vec[j][0][l]<0.65) zbin = 8;
          else if(0.65<Pvszlocal[i].vec[j][0][l] && Pvszlocal[i].vec[j][0][l]<0.70) zbin = 9;
          else if(0.70<Pvszlocal[i].vec[j][0][l] && Pvszlocal[i].vec[j][0][l]<0.75) zbin = 10;
          else zbin = 11;

          fKinematicsRD[3]->Fill(Pvszlocal[i].vec[j][0][l]);

          for(int ll=0; ll<4; ll++)
          {
            fBinning[xbin][ybin][zbin].tab[j][0][ll] += Pvszlocal[i].vec[j][ll+1][l];
            //fBinning[xbin][ybin][zbin].tab[j][1][ll] += Pvsz_errlocal[i].vec[j][ll+1][l];
            fMeanvalues[xbin][ybin][zbin].vec[j][ll][2].push_back(Q2local[i]);
            fMeanvalues[xbin][ybin][zbin].vec[j][ll][0].push_back(XBjlocal[i]);
            fMeanvalues[xbin][ybin][zbin].vec[j][ll][1].push_back(YBjlocal[i]);
            fMeanvalues[xbin][ybin][zbin].vec[j][ll][3].push_back(Pvszlocal[i].vec[j][0][l]);
          }
        }
      }
    }

    Pvszlocal.clear();
    Pvsz_errlocal.clear();
    XBjlocal.clear();
    YBjlocal.clear();
    Q2local.clear();
    Pvszloose.clear();
    XBjloose.clear();
    YBjloose.clear();
    Q2loose.clear();
    Pvszsevere.clear();
    XBjsevere.clear();
    YBjsevere.clear();
    Q2severe.clear();
  }

  for(int i=0; i<int(fQ2.size()); i++)
  {
    fKinematicsRD[0]->Fill(fQ2kin[i]);
    fKinematicsRD[1]->Fill(fXBjkin[i]);
    fKinematicsRD[2]->Fill(fYBjkin[i]);
    fKinematicsRD[4]->Fill(fWBjkin[i]);
    fKinematicsRD[5]->Fill(fNukin[i]);
  }

}

int main(int argc, char **argv)
{

  create_kin_plots();
  RDextraction(argv[1]);
  MCextraction(argv[2]);
  save_kin_plots();

  return 0;
}

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

#include "compMCMC.h"

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
  for(int i=0; i<5; i++)
  {
    fKinematicsRD[i][0] = new TH1F(Form("Q^{2} %s",trigname[i].c_str()), Form("Q^{2} %s",trigname[i].c_str()), 100, -1, 2);
    fKinematicsRD[i][1] = new TH1F(Form("x_{Bj} %s",trigname[i].c_str()), Form("x_{Bj} %s",trigname[i].c_str()), 100, -3, 0);
    fKinematicsRD[i][2] = new TH1F(Form("y %s",trigname[i].c_str()), Form("y %s",trigname[i].c_str()), 100, 0, 1);
    fKinematicsRD[i][3] = new TH1F(Form("z %s",trigname[i].c_str()), Form("z %s",trigname[i].c_str()), 100, 0, 1);
    fKinematicsRD[i][4] = new TH1F(Form("W %s",trigname[i].c_str()), Form("W %s",trigname[i].c_str()), 100, 2, 18);
    fKinematicsRD[i][5] = new TH1F(Form("#nu %s",trigname[i].c_str()), Form("#nu %s",trigname[i].c_str()), 100, 0, 160);
    fKinematicsRD[i][6] = new TH1F(Form("E_{#mu} %s",trigname[i].c_str()), Form("E_{#mu} %s",trigname[i].c_str()), 100, 140, 180);
    fKinematicsRD[i][7] = new TH1F(Form("E_{#mu'} %s",trigname[i].c_str()), Form("E_{#mu'} %s",trigname[i].c_str()), 100, 0, 160);
    fKinematicsRD[i][8] = new TH1F(Form("#theta %s",trigname[i].c_str()), Form("#theta %s",trigname[i].c_str()), 100, 0, 0.1);
    fKinematicsRD[i][9] = new TH1F(Form("#phi %s",trigname[i].c_str()), Form("#phi %s",trigname[i].c_str()), 100, -1.7, 1.7);
    fKinematicsRD[i][10] = new TH1F(Form("Vertex %s",trigname[i].c_str()), Form("Vertex %s",trigname[i].c_str()), 100, -320, -70);
    fKinematicsRD[i][11] = new TH1F("#Phi_h","#Phi_h", 100, 0, 1);
    fKinematicsMC[i][0] = new TH1F(Form("Q^{2} MC %s",trigname[i].c_str()), Form("Q^{2} MC %s",trigname[i].c_str()), 100, -1, 2);
    fKinematicsMC[i][1] = new TH1F(Form("x_{Bj} MC %s",trigname[i].c_str()), Form("x_{Bj} MC %s",trigname[i].c_str()), 100, -3, 0);
    fKinematicsMC[i][2] = new TH1F(Form("y MC %s",trigname[i].c_str()), Form("y MC %s",trigname[i].c_str()), 100, 0, 1);
    fKinematicsMC[i][3] = new TH1F(Form("z MC %s",trigname[i].c_str()), Form("z MC %s",trigname[i].c_str()), 100, 0, 1);
    fKinematicsMC[i][4] = new TH1F(Form("W MC %s",trigname[i].c_str()), Form("W MC %s",trigname[i].c_str()), 100, 2, 18);
    fKinematicsMC[i][5] = new TH1F(Form("#nu MC %s",trigname[i].c_str()), Form("#nu MC %s",trigname[i].c_str()), 100, 0, 160);
    fKinematicsMC[i][6] = new TH1F(Form("E_{#mu} MC %s",trigname[i].c_str()), Form("E_{#mu} MC %s",trigname[i].c_str()), 100, 140, 180);
    fKinematicsMC[i][7] = new TH1F(Form("E_{#mu'} MC %s",trigname[i].c_str()), Form("E_{#mu'} MC %s",trigname[i].c_str()), 100, 0, 160);
    fKinematicsMC[i][8] = new TH1F(Form("#theta MC %s",trigname[i].c_str()), Form("#theta MC %s",trigname[i].c_str()), 100, 0, 0.1);
    fKinematicsMC[i][9] = new TH1F(Form("#phi MC %s",trigname[i].c_str()), Form("#phi MC %s",trigname[i].c_str()), 100, -1.7, 1.7);
    fKinematicsMC[i][10] = new TH1F(Form("Vertex MC %s",trigname[i].c_str()), Form("Vertex MC %s",trigname[i].c_str()), 100, -320, -70);
    fKinematicsMC[i][11] = new TH1F("#Phi_h MC","#Phi_h MC", 100, 0, 1);
    BinLogX(fKinematicsRD[i][0]);
    BinLogX(fKinematicsMC[i][0]);
    BinLogX(fKinematicsRD[i][1]);
    BinLogX(fKinematicsMC[i][1]);
  }
}

void save_kin_plots()
{
  c1.Divide(2,4);
  c2.Divide(2,4);
  c3.Divide(2,4);
  c4.Divide(2,4);
  c5.Divide(2,4);
  c6.Divide(2,4);
  c7.Divide(1,2);
  c8.Divide(1,2);
  c9.Divide(1,2);
  c10.Divide(1,2);
  c11.Divide(1,2);
  c12.Divide(1,2);
  c13.Divide(1,2);
  c14.Divide(2,4);
  c15.Divide(2,4);
  c16.Divide(2,4);
  c17.Divide(2,4);
  c18.Divide(2,4);

  for(int i=0; i<8; i++)
  {
    int idx=int(i/2);
    if(i%2)
    {
      c1.cd(idx+3+int(idx/2)*2);
      // TPad *pad1 = new TPad("pad1","pad1",0+i%2*0.5,0.6-i%2*0.5,0.5+i%2*0.5,0.7-i%2*0.5);
      fKinematicsRatio[idx][0] = (TH1F*)fKinematicsRD[idx][0]->Clone();
      fKinematicsRatio[idx][0]->SetStats(0);
      fKinematicsRatio[idx][0]->Divide(fKinematicsMC[idx][0]);
      fKinematicsRatio[idx][0]->SetMarkerStyle(21);
      fKinematicsRatio[idx][0]->SetMaximum(2.);
      fKinematicsRatio[idx][0]->Draw("P");
      gPad->SetLogx();
      c1.Update();

      c2.cd(idx+3+int(idx/2)*2);
      // TPad *pad2 = new TPad("pad2","pad2",0+i%2*0.5,0.6-i%2*0.5,0.5+i%2*0.5,0.7-i%2*0.5);
      fKinematicsRatio[idx][1] = (TH1F*)fKinematicsRD[idx][1]->Clone();
      fKinematicsRatio[idx][1]->SetStats(0);
      fKinematicsRatio[idx][1]->Divide(fKinematicsMC[idx][1]);
      fKinematicsRatio[idx][1]->SetMarkerStyle(21);
      fKinematicsRatio[idx][1]->SetMaximum(2.);
      fKinematicsRatio[idx][1]->Draw("P");
      gPad->SetLogx();
      c2.Update();

      c3.cd(idx+3+int(idx/2)*2);
      // TPad *pad3 = new TPad("pad3","pad3",0+i%2*0.5,0.6-i%2*0.5,0.5+i%2*0.5,0.7-i%2*0.5);
      fKinematicsRatio[idx][2] = (TH1F*)fKinematicsRD[idx][2]->Clone();
      fKinematicsRatio[idx][2]->SetStats(0);
      fKinematicsRatio[idx][2]->Divide(fKinematicsMC[idx][2]);
      fKinematicsRatio[idx][2]->SetMarkerStyle(21);
      fKinematicsRatio[idx][2]->SetMaximum(2.);
      fKinematicsRatio[idx][2]->Draw("P");
      c3.Update();

      c4.cd(idx+3+int(idx/2)*2);
      fKinematicsRatio[idx][3] = (TH1F*)fKinematicsRD[idx][3]->Clone();
      fKinematicsRatio[idx][3]->SetStats(0);
      fKinematicsRatio[idx][3]->Divide(fKinematicsMC[idx][3]);
      fKinematicsRatio[idx][3]->SetMarkerStyle(21);
      fKinematicsRatio[idx][3]->SetMaximum(2.);
      fKinematicsRatio[idx][3]->Draw("P");
      c4.Update();

      c5.cd(idx+3+int(idx/2)*2);
      // TPad *pad4 = new TPad("pad4","pad4",0+i%2*0.5,0.6-i%2*0.5,0.5+i%2*0.5,0.7-i%2*0.5);
      fKinematicsRatio[idx][4] = (TH1F*)fKinematicsRD[idx][4]->Clone();
      fKinematicsRatio[idx][4]->SetStats(0);
      fKinematicsRatio[idx][4]->Divide(fKinematicsMC[idx][4]);
      fKinematicsRatio[idx][4]->SetMarkerStyle(21);
      fKinematicsRatio[idx][4]->SetMaximum(2.);
      fKinematicsRatio[idx][4]->Draw("P");
      c5.Update();

      c6.cd(idx+3+int(idx/2)*2);
      // TPad *pad5 = new TPad("pad5","pad5",0+i%2*0.5,0.6-i%2*0.5,0.5+i%2*0.5,0.7-i%2*0.5);
      fKinematicsRatio[idx][5] = (TH1F*)fKinematicsRD[idx][5]->Clone();
      fKinematicsRatio[idx][5]->SetStats(0);
      fKinematicsRatio[idx][5]->Divide(fKinematicsMC[idx][5]);
      fKinematicsRatio[idx][5]->SetMarkerStyle(21);
      fKinematicsRatio[idx][5]->SetMaximum(2.);
      fKinematicsRatio[idx][5]->Draw("P");
      c6.Update();

      c14.cd(idx+3+int(idx/2)*2);
      // TPad *pad5 = new TPad("pad5","pad5",0+i%2*0.5,0.6-i%2*0.5,0.5+i%2*0.5,0.7-i%2*0.5);
      fKinematicsRatio[idx][6] = (TH1F*)fKinematicsRD[idx][6]->Clone();
      fKinematicsRatio[idx][6]->SetStats(0);
      fKinematicsRatio[idx][6]->Divide(fKinematicsMC[idx][6]);
      fKinematicsRatio[idx][6]->SetMarkerStyle(21);
      fKinematicsRatio[idx][6]->SetMaximum(2.);
      fKinematicsRatio[idx][6]->Draw("P");
      c14.Update();

      c15.cd(idx+3+int(idx/2)*2);
      // TPad *pad5 = new TPad("pad5","pad5",0+i%2*0.5,0.6-i%2*0.5,0.5+i%2*0.5,0.7-i%2*0.5);
      fKinematicsRatio[idx][7] = (TH1F*)fKinematicsRD[idx][7]->Clone();
      fKinematicsRatio[idx][7]->SetStats(0);
      fKinematicsRatio[idx][7]->Divide(fKinematicsMC[idx][7]);
      fKinematicsRatio[idx][7]->SetMarkerStyle(21);
      fKinematicsRatio[idx][7]->SetMaximum(2.);
      fKinematicsRatio[idx][7]->Draw("P");
      c15.Update();

      c16.cd(idx+3+int(idx/2)*2);
      // TPad *pad5 = new TPad("pad5","pad5",0+i%2*0.5,0.6-i%2*0.5,0.5+i%2*0.5,0.7-i%2*0.5);
      fKinematicsRatio[idx][8] = (TH1F*)fKinematicsRD[idx][8]->Clone();
      fKinematicsRatio[idx][8]->SetStats(0);
      fKinematicsRatio[idx][8]->Divide(fKinematicsMC[idx][8]);
      fKinematicsRatio[idx][8]->SetMarkerStyle(21);
      fKinematicsRatio[idx][8]->SetMaximum(2.);
      fKinematicsRatio[idx][8]->Draw("P");
      c16.Update();

      c17.cd(idx+3+int(idx/2)*2);
      // TPad *pad5 = new TPad("pad5","pad5",0+i%2*0.5,0.6-i%2*0.5,0.5+i%2*0.5,0.7-i%2*0.5);
      fKinematicsRatio[idx][9] = (TH1F*)fKinematicsRD[idx][9]->Clone();
      fKinematicsRatio[idx][9]->SetStats(0);
      fKinematicsRatio[idx][9]->Divide(fKinematicsMC[idx][9]);
      fKinematicsRatio[idx][9]->SetMarkerStyle(21);
      fKinematicsRatio[idx][9]->SetMaximum(2.);
      fKinematicsRatio[idx][9]->Draw("P");
      c17.Update();

      c18.cd(idx+3+int(idx/2)*2);
      // TPad *pad5 = new TPad("pad5","pad5",0+i%2*0.5,0.6-i%2*0.5,0.5+i%2*0.5,0.7-i%2*0.5);
      fKinematicsRatio[idx][10] = (TH1F*)fKinematicsRD[idx][10]->Clone();
      fKinematicsRatio[idx][10]->SetStats(0);
      fKinematicsRatio[idx][10]->Divide(fKinematicsMC[idx][10]);
      fKinematicsRatio[idx][10]->SetMarkerStyle(21);
      fKinematicsRatio[idx][10]->SetMaximum(2.);
      fKinematicsRatio[idx][10]->Draw("P");
      c18.Update();
    }
    else
    {
      c1.cd(idx+1+int(idx/2)*2);
      // TPad *pad1 = new TPad("pad1","pad1",0+i%2*0.5,0.7-i%2*0.5,0.5+i%2*0.5,1-i%2*0.5);
      fKinematicsRD[idx][0]->Scale(1/fKinematicsRD[2][0]->GetEntries());
      fKinematicsMC[idx][0]->Scale(1/fKinematicsMC[2][0]->GetEntries());
      fKinematicsRD[idx][0]->SetLineColor(kRed);
      fKinematicsRD[idx][0]->SetStats(0);
      fKinematicsRD[idx][0]->SetMinimum(0.);
      // fKinematicsRD[idx][0]->SetMaximum(0.06);
      fKinematicsRD[idx][0]->Draw();
      fKinematicsMC[idx][0]->SetLineColor(kBlue);
      fKinematicsMC[idx][0]->Draw("SAME");
      gPad->SetLogx();
      c1.Update();

      c2.cd(idx+1+int(idx/2)*2);
      // TPad *pad2 = new TPad("pad2","pad2",0+i%2*0.5,0.7-i%2*0.5,0.5+i%2*0.5,1-i%2*0.5);
      fKinematicsRD[idx][1]->Scale(1/fKinematicsRD[2][1]->GetEntries());
      fKinematicsMC[idx][1]->Scale(1/fKinematicsMC[2][1]->GetEntries());
      fKinematicsRD[idx][1]->SetLineColor(kRed);
      fKinematicsRD[idx][1]->SetStats(0);
      fKinematicsRD[idx][1]->SetMinimum(0.);
      // fKinematicsRD[idx][1]->SetMaximum(0.045);
      fKinematicsRD[idx][1]->Draw();
      fKinematicsMC[idx][1]->SetLineColor(kBlue);
      fKinematicsMC[idx][1]->Draw("SAME");
      gPad->SetLogx();
      c2.Update();

      c3.cd(idx+1+int(idx/2)*2);
      // TPad *pad3 = new TPad("pad3","pad3",0+i%2*0.5,0.7-i%2*0.5,0.5+i%2*0.5,1-i%2*0.5);
      fKinematicsRD[idx][2]->Scale(1/fKinematicsRD[2][2]->GetEntries());
      fKinematicsMC[idx][2]->Scale(1/fKinematicsMC[2][2]->GetEntries());
      fKinematicsRD[idx][2]->SetLineColor(kRed);
      fKinematicsRD[idx][2]->SetStats(0);
      fKinematicsRD[idx][2]->SetMinimum(0.);
      // fKinematicsRD[idx][2]->SetMaximum(0.06);
      fKinematicsRD[idx][2]->Draw();
      fKinematicsMC[idx][2]->SetLineColor(kBlue);
      fKinematicsMC[idx][2]->Draw("SAME");
      c3.Update();

      c4.cd(idx+1+int(idx/2)*2);
      fKinematicsRD[idx][3]->Scale(1/fKinematicsRD[idx][3]->GetEntries());
      fKinematicsMC[idx][3]->Scale(1/fKinematicsMC[idx][3]->GetEntries());
      fKinematicsRD[idx][3]->SetLineColor(kRed);
      fKinematicsRD[idx][3]->SetStats(0);
      fKinematicsRD[idx][2]->SetMinimum(0.);
      fKinematicsRD[idx][3]->Draw();
      fKinematicsMC[idx][3]->SetLineColor(kBlue);
      fKinematicsMC[idx][3]->Draw("SAME");
      c4.Update();

      c5.cd(idx+1+int(idx/2)*2);
      // TPad *pad4 = new TPad("pad4","pad4",0+i%2*0.5,0.7-i%2*0.5,0.5+i%2*0.5,1-i%2*0.5);
      fKinematicsRD[idx][4]->Scale(1/fKinematicsRD[2][4]->GetEntries());
      fKinematicsMC[idx][4]->Scale(1/fKinematicsMC[2][4]->GetEntries());
      fKinematicsRD[idx][4]->SetLineColor(kRed);
      fKinematicsRD[idx][4]->SetStats(0);
      fKinematicsRD[idx][4]->SetMinimum(0.);
      // fKinematicsRD[idx][4]->SetMaximum(0.05);
      fKinematicsRD[idx][4]->Draw();
      fKinematicsMC[idx][4]->SetLineColor(kBlue);
      fKinematicsMC[idx][4]->Draw("SAME");
      c5.Update();

      c6.cd(idx+1+int(idx/2)*2);
      // TPad *pad5 = new TPad("pad5","pad5",0+i%2*0.5,0.7-i%2*0.5,0.5+i%2*0.5,1-i%2*0.5);
      fKinematicsRD[idx][5]->Scale(1/fKinematicsRD[2][5]->GetEntries());
      fKinematicsMC[idx][5]->Scale(1/fKinematicsMC[2][5]->GetEntries());
      fKinematicsRD[idx][5]->SetLineColor(kRed);
      fKinematicsRD[idx][5]->SetStats(0);
      fKinematicsRD[idx][5]->SetMinimum(0.);
      // fKinematicsRD[idx][5]->SetMaximum(0.06);
      fKinematicsRD[idx][5]->Draw();
      fKinematicsMC[idx][5]->SetLineColor(kBlue);
      fKinematicsMC[idx][5]->Draw("SAME");
      c6.Update();

      c14.cd(idx+1+int(idx/2)*2);
      // TPad *pad5 = new TPad("pad5","pad5",0+i%2*0.5,0.7-i%2*0.5,0.5+i%2*0.5,1-i%2*0.5);
      fKinematicsRD[idx][6]->Scale(1/fKinematicsRD[2][6]->GetEntries());
      fKinematicsMC[idx][6]->Scale(1/fKinematicsMC[2][6]->GetEntries());
      fKinematicsRD[idx][6]->SetLineColor(kRed);
      fKinematicsRD[idx][6]->SetStats(0);
      fKinematicsRD[idx][6]->SetMinimum(0.);
      // fKinematicsRD[idx][5]->SetMaximum(0.06);
      fKinematicsRD[idx][6]->Draw();
      fKinematicsMC[idx][6]->SetLineColor(kBlue);
      fKinematicsMC[idx][6]->Draw("SAME");
      c14.Update();

      c15.cd(idx+1+int(idx/2)*2);
      // TPad *pad5 = new TPad("pad5","pad5",0+i%2*0.5,0.7-i%2*0.5,0.5+i%2*0.5,1-i%2*0.5);
      fKinematicsRD[idx][7]->Scale(1/fKinematicsRD[2][7]->GetEntries());
      fKinematicsMC[idx][7]->Scale(1/fKinematicsMC[2][7]->GetEntries());
      fKinematicsRD[idx][7]->SetLineColor(kRed);
      fKinematicsRD[idx][7]->SetStats(0);
      fKinematicsRD[idx][7]->SetMinimum(0.);
      // fKinematicsRD[idx][5]->SetMaximum(0.06);
      fKinematicsRD[idx][7]->Draw();
      fKinematicsMC[idx][7]->SetLineColor(kBlue);
      fKinematicsMC[idx][7]->Draw("SAME");
      c15.Update();

      c16.cd(idx+1+int(idx/2)*2);
      // TPad *pad5 = new TPad("pad5","pad5",0+i%2*0.5,0.7-i%2*0.5,0.5+i%2*0.5,1-i%2*0.5);
      fKinematicsRD[idx][8]->Scale(1/fKinematicsRD[2][8]->GetEntries());
      fKinematicsMC[idx][8]->Scale(1/fKinematicsMC[2][8]->GetEntries());
      fKinematicsRD[idx][8]->SetLineColor(kRed);
      fKinematicsRD[idx][8]->SetStats(0);
      fKinematicsRD[idx][8]->SetMinimum(0.);
      // fKinematicsRD[idx][5]->SetMaximum(0.06);
      fKinematicsRD[idx][8]->Draw();
      fKinematicsMC[idx][8]->SetLineColor(kBlue);
      fKinematicsMC[idx][8]->Draw("SAME");
      c16.Update();

      c17.cd(idx+1+int(idx/2)*2);
      // TPad *pad5 = new TPad("pad5","pad5",0+i%2*0.5,0.7-i%2*0.5,0.5+i%2*0.5,1-i%2*0.5);
      fKinematicsRD[idx][9]->Scale(1/fKinematicsRD[2][9]->GetEntries());
      fKinematicsMC[idx][9]->Scale(1/fKinematicsMC[2][9]->GetEntries());
      fKinematicsRD[idx][9]->SetLineColor(kRed);
      fKinematicsRD[idx][9]->SetStats(0);
      fKinematicsRD[idx][9]->SetMinimum(0.);
      // fKinematicsRD[idx][5]->SetMaximum(0.06);
      fKinematicsRD[idx][9]->Draw();
      fKinematicsMC[idx][9]->SetLineColor(kBlue);
      fKinematicsMC[idx][9]->Draw("SAME");
      c17.Update();

      c18.cd(idx+1+int(idx/2)*2);
      // TPad *pad5 = new TPad("pad5","pad5",0+i%2*0.5,0.7-i%2*0.5,0.5+i%2*0.5,1-i%2*0.5);
      fKinematicsRD[idx][10]->Scale(1/fKinematicsRD[2][10]->GetEntries());
      fKinematicsMC[idx][10]->Scale(1/fKinematicsMC[2][10]->GetEntries());
      fKinematicsRD[idx][10]->SetLineColor(kRed);
      fKinematicsRD[idx][10]->SetStats(0);
      fKinematicsRD[idx][10]->SetMinimum(0.);
      // fKinematicsRD[idx][5]->SetMaximum(0.06);
      fKinematicsRD[idx][10]->Draw();
      fKinematicsMC[idx][10]->SetLineColor(kBlue);
      fKinematicsMC[idx][10]->Draw("SAME");
      c18.Update();
    }
  }

  c7.cd(2);
  fKinematicsRD[0][11]->Scale(1/fKinematicsRD[0][11]->GetEntries());
  fKinematicsMC[0][11]->Scale(1/fKinematicsMC[0][11]->GetEntries());
  fKinematicsRatio[0][11] = (TH1F*)fKinematicsRD[0][11]->Clone();
  fKinematicsRatio[0][11]->SetStats(0);
  fKinematicsRatio[0][11]->Divide(fKinematicsMC[0][11]);
  fKinematicsRatio[0][11]->SetMarkerStyle(21);
  fKinematicsRatio[0][11]->Draw("P");
  c7.Update();
  c7.cd(1);
  fKinematicsRD[0][11]->SetLineColor(kRed);
  fKinematicsRD[0][11]->SetStats(0);
  fKinematicsRD[0][11]->SetMinimum(0.);
  fKinematicsRD[0][11]->SetMaximum(0.05);
  fKinematicsRD[0][11]->Draw();
  fKinematicsMC[0][11]->SetLineColor(kBlue);
  fKinematicsMC[0][11]->Draw("SAME");
  c7.Update();

  c8.cd(2);
  fKinematicsRD[4][0]->Scale(1/fKinematicsRD[4][0]->GetEntries());
  fKinematicsMC[4][0]->Scale(1/fKinematicsMC[4][0]->GetEntries());
  fKinematicsRatio[4][0] = (TH1F*)fKinematicsRD[4][0]->Clone();
  fKinematicsRatio[4][0]->SetStats(0);
  fKinematicsRatio[4][0]->Divide(fKinematicsMC[4][0]);
  fKinematicsRatio[4][0]->SetMarkerStyle(21);
  fKinematicsRatio[4][0]->SetMaximum(2.);
  fKinematicsRatio[4][0]->Draw("P");
  gPad->SetLogx();
  c8.Update();
  c8.cd(1);
  fKinematicsRD[4][0]->SetLineColor(kRed);
  fKinematicsRD[4][0]->SetStats(0);
  fKinematicsRD[4][0]->Draw();
  fKinematicsMC[4][0]->SetLineColor(kBlue);
  fKinematicsMC[4][0]->Draw("SAME");
  gPad->SetLogx();
  c8.Update();

  c9.cd(2);
  fKinematicsRD[4][1]->Scale(1/fKinematicsRD[4][1]->GetEntries());
  fKinematicsMC[4][1]->Scale(1/fKinematicsMC[4][1]->GetEntries());
  fKinematicsRatio[4][1] = (TH1F*)fKinematicsRD[4][1]->Clone();
  fKinematicsRatio[4][1]->SetStats(0);
  fKinematicsRatio[4][1]->Divide(fKinematicsMC[4][1]);
  fKinematicsRatio[4][1]->SetMarkerStyle(21);
  fKinematicsRatio[4][1]->SetMaximum(2.);
  fKinematicsRatio[4][1]->Draw("P");
  gPad->SetLogx();
  c9.Update();
  c9.cd(1);
  fKinematicsRD[4][1]->SetLineColor(kRed);
  fKinematicsRD[4][1]->SetStats(0);
  fKinematicsRD[4][1]->Draw();
  fKinematicsMC[4][1]->SetLineColor(kBlue);
  fKinematicsMC[4][1]->Draw("SAME");
  gPad->SetLogx();
  c9.Update();

  c10.cd(2);
  fKinematicsRD[4][2]->Scale(1/fKinematicsRD[4][2]->GetEntries());
  fKinematicsMC[4][2]->Scale(1/fKinematicsMC[4][2]->GetEntries());
  fKinematicsRatio[4][2] = (TH1F*)fKinematicsRD[4][2]->Clone();
  fKinematicsRatio[4][2]->SetStats(0);
  fKinematicsRatio[4][2]->Divide(fKinematicsMC[4][2]);
  fKinematicsRatio[4][2]->SetMarkerStyle(21);
  fKinematicsRatio[4][2]->SetMaximum(2.);
  fKinematicsRatio[4][2]->Draw("P");
  c10.Update();
  c10.cd(1);
  fKinematicsRD[4][2]->SetLineColor(kRed);
  fKinematicsRD[4][2]->SetStats(0);
  fKinematicsRD[4][2]->Draw();
  fKinematicsMC[4][2]->SetLineColor(kBlue);
  fKinematicsMC[4][2]->Draw("SAME");
  c10.Update();

  c11.cd(2);
  fKinematicsRD[4][3]->Scale(1/fKinematicsRD[4][3]->GetEntries());
  fKinematicsMC[4][3]->Scale(1/fKinematicsMC[4][3]->GetEntries());
  fKinematicsRatio[4][3] = (TH1F*)fKinematicsRD[4][3]->Clone();
  fKinematicsRatio[4][3]->SetStats(0);
  fKinematicsRatio[4][3]->Divide(fKinematicsMC[4][3]);
  fKinematicsRatio[4][3]->SetMarkerStyle(21);
  fKinematicsRatio[4][3]->SetMaximum(2.);
  fKinematicsRatio[4][3]->Draw("P");
  c11.Update();
  c11.cd(1);
  fKinematicsRD[4][3]->SetLineColor(kRed);
  fKinematicsRD[4][3]->SetStats(0);
  fKinematicsRD[4][3]->Draw();
  fKinematicsMC[4][3]->SetLineColor(kBlue);
  fKinematicsMC[4][3]->Draw("SAME");
  c11.Update();

  c12.cd(2);
  fKinematicsRD[4][4]->Scale(1/fKinematicsRD[4][4]->GetEntries());
  fKinematicsMC[4][4]->Scale(1/fKinematicsMC[4][4]->GetEntries());
  fKinematicsRatio[4][4] = (TH1F*)fKinematicsRD[4][4]->Clone();
  fKinematicsRatio[4][4]->SetStats(0);
  fKinematicsRatio[4][4]->Divide(fKinematicsMC[4][4]);
  fKinematicsRatio[4][4]->SetMarkerStyle(21);
  fKinematicsRatio[4][4]->SetMaximum(2.);
  fKinematicsRatio[4][4]->Draw("P");
  c12.Update();
  c12.cd(1);
  fKinematicsRD[4][4]->SetLineColor(kRed);
  fKinematicsRD[4][4]->SetStats(0);
  fKinematicsRD[4][4]->Draw();
  fKinematicsMC[4][4]->SetLineColor(kBlue);
  fKinematicsMC[4][4]->Draw("SAME");
  c12.Update();

  c13.cd(2);
  fKinematicsRD[4][5]->Scale(1/fKinematicsRD[4][5]->GetEntries());
  fKinematicsMC[4][5]->Scale(1/fKinematicsMC[4][5]->GetEntries());
  fKinematicsRatio[4][5] = (TH1F*)fKinematicsRD[4][5]->Clone();
  fKinematicsRatio[4][5]->SetStats(0);
  fKinematicsRatio[4][5]->Divide(fKinematicsMC[4][5]);
  fKinematicsRatio[4][5]->SetMarkerStyle(21);
  fKinematicsRatio[4][5]->SetMaximum(2.);
  fKinematicsRatio[4][5]->Draw("P");
  c13.Update();
  c13.cd(1);
  fKinematicsRD[4][5]->SetLineColor(kRed);
  fKinematicsRD[4][5]->SetStats(0);
  fKinematicsRD[4][5]->Draw();
  fKinematicsMC[4][5]->SetLineColor(kBlue);
  fKinematicsMC[4][5]->Draw("SAME");
  c13.Update();

  c1.Print("kinMCMC.pdf(","pdf");
  c2.Print("kinMCMC.pdf","pdf");
  c3.Print("kinMCMC.pdf","pdf");
  c4.Print("kinMCMC.pdf","pdf");
  c5.Print("kinMCMC.pdf","pdf");
  c6.Print("kinMCMC.pdf","pdf");
  c7.Print("kinMCMC.pdf","pdf");
  c8.Print("kinMCMC.pdf","pdf");
  c9.Print("kinMCMC.pdf","pdf");
  c10.Print("kinMCMC.pdf","pdf");
  c11.Print("kinMCMC.pdf","pdf");
  c12.Print("kinMCMC.pdf","pdf");
  c13.Print("kinMCMC.pdf","pdf");
  c14.Print("kinMCMC.pdf","pdf");
  c15.Print("kinMCMC.pdf","pdf");
  c16.Print("kinMCMC.pdf","pdf");
  c17.Print("kinMCMC.pdf","pdf");
  c18.Print("kinMCMC.pdf)","pdf");
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
      trig = (trig&2047);


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

                  if((trig&2 || trig&4 || trig&8))
                  // if((trig&2 || trig&4 || trig&8 || trig&512))
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
          //2016 ---
        }
      }

      if(fAllDISflag)
      {
        double theta_m = asin(sqrt(pow(p1x->GetLeaf("p1x")->GetValue()/sqrt(pow(E_mu_prim->GetLeaf("E_mu_prim")->GetValue(),2)-pow(fM_mu,2)),2)+pow(p1y->GetLeaf("p1y")->GetValue()/sqrt(pow(E_mu_prim->GetLeaf("E_mu_prim")->GetValue(),2)-pow(fM_mu,2)),2)));
        double phi_m = asin(p1x->GetLeaf("p1x")->GetValue()/sqrt(pow(p1x->GetLeaf("p1x")->GetValue(),2)+pow(p1y->GetLeaf("p1y")->GetValue(),2)));

        if(trig&2)
        {
          fQ2kinMC[0].push_back(Q2);
          fXBjkinMC[0].push_back(xBj);
          fYBjkinMC[0].push_back(yBj);
          fWBjkinMC[0].push_back(sqrt(wBj));
          fNukinMC[0].push_back(nu);
          fMuMC[0].push_back(E_beam->GetLeaf("E_beam")->GetValue());
          fMupMC[0].push_back(E_mu_prim->GetLeaf("E_mu_prim")->GetValue());
          fThetaMC[0].push_back(theta_m);
          fPhiMC[0].push_back(phi_m);
          fVertexMC[0].push_back(z->GetLeaf("z")->GetValue());
        }
        if(trig&4)
        {
          fQ2kinMC[1].push_back(Q2);
          fXBjkinMC[1].push_back(xBj);
          fYBjkinMC[1].push_back(yBj);
          fWBjkinMC[1].push_back(sqrt(wBj));
          fNukinMC[1].push_back(nu);
          fMuMC[1].push_back(E_beam->GetLeaf("E_beam")->GetValue());
          fMupMC[1].push_back(E_mu_prim->GetLeaf("E_mu_prim")->GetValue());
          fThetaMC[1].push_back(theta_m);
          fPhiMC[1].push_back(phi_m);
          fVertexMC[1].push_back(z->GetLeaf("z")->GetValue());
        }
        if(trig&8)
        {
          fQ2kinMC[2].push_back(Q2);
          fXBjkinMC[2].push_back(xBj);
          fYBjkinMC[2].push_back(yBj);
          fWBjkinMC[2].push_back(sqrt(wBj));
          fNukinMC[2].push_back(nu);
          fMuMC[2].push_back(E_beam->GetLeaf("E_beam")->GetValue());
          fMupMC[2].push_back(E_mu_prim->GetLeaf("E_mu_prim")->GetValue());
          fThetaMC[2].push_back(theta_m);
          fPhiMC[2].push_back(phi_m);
          fVertexMC[2].push_back(z->GetLeaf("z")->GetValue());
        }
        if(trig&512)
        {
          fQ2kinMC[3].push_back(Q2);
          fXBjkinMC[3].push_back(xBj);
          fYBjkinMC[3].push_back(yBj);
          fWBjkinMC[3].push_back(sqrt(wBj));
          fNukinMC[3].push_back(nu);
          fMuMC[3].push_back(E_beam->GetLeaf("E_beam")->GetValue());
          fMupMC[3].push_back(E_mu_prim->GetLeaf("E_mu_prim")->GetValue());
          fThetaMC[3].push_back(theta_m);
          fPhiMC[3].push_back(phi_m);
          fVertexMC[3].push_back(z->GetLeaf("z")->GetValue());
        }
        if(trig&2 || trig&4 || trig&8)
        // if(trig&2 || trig&4 || trig&8 || trig&512)
        {
          fQ2kinMC[4].push_back(Q2);
          fXBjkinMC[4].push_back(xBj);
          fYBjkinMC[4].push_back(yBj);
          fWBjkinMC[4].push_back(sqrt(wBj));
          fNukinMC[4].push_back(nu);
          fMuMC[4].push_back(E_beam->GetLeaf("E_beam")->GetValue());
          fMupMC[4].push_back(E_mu_prim->GetLeaf("E_mu_prim")->GetValue());
          fThetaMC[4].push_back(theta_m);
          fPhiMC[4].push_back(phi_m);
          fVertexMC[4].push_back(z->GetLeaf("z")->GetValue());
        }
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

      if(0.<xBj_MC && xBj_MC<0.01) xbin_MC = 0;
      else if(0.01<=xBj_MC && xBj_MC<0.02) xbin_MC = 1;
      else if(0.02<=xBj_MC && xBj_MC<0.03) xbin_MC = 2;
      else if(0.03<=xBj_MC && xBj_MC<0.04) xbin_MC = 3;
      else if(0.04<=xBj_MC && xBj_MC<0.06) xbin_MC = 4;
      else if(0.06<=xBj_MC && xBj_MC<0.1) xbin_MC = 5;
      else if(0.1<=xBj_MC && xBj_MC<0.14) xbin_MC = 6;
      else if(0.14<=xBj_MC && xBj_MC<0.18) xbin_MC = 7;
      else xbin_MC = 8;

      if(0.<xBj && xBj<0.01) xbin = 0;
      else if(0.01<=xBj && xBj<0.02) xbin = 1;
      else if(0.02<=xBj && xBj<0.03) xbin = 2;
      else if(0.03<=xBj && xBj<0.04) xbin = 3;
      else if(0.04<=xBj && xBj<0.06) xbin = 4;
      else if(0.06<=xBj && xBj<0.1) xbin = 5;
      else if(0.1<=xBj && xBj<0.14) xbin = 6;
      else if(0.14<=xBj && xBj<0.18) xbin = 7;
      else xbin = 8;

      // y Binning

      if(0.<yBj_MC && yBj_MC<0.15) ybin_MC = 0;
      else if(0.15<=yBj_MC && yBj_MC<0.2) ybin_MC = 1;
      else if(0.2<=yBj_MC && yBj_MC<0.3) ybin_MC = 2;
      else if(0.3<=yBj_MC && yBj_MC<0.5) ybin_MC = 3;
      else ybin_MC = 4;

      if(0.<yBj && yBj<0.15) ybin = 0;
      else if(0.15<=yBj && yBj<0.2) ybin = 1;
      else if(0.2<=yBj && yBj<0.3) ybin = 2;
      else if(0.3<=yBj && yBj<0.5) ybin = 3;
      else ybin = 4;

      // if(fAllDISflag_MC)
      // {
      //   // z Binnig
      //
      //   for(int i=0; i<12; i++)
      //   {
      //     fNDIS_evt_MC[0][xbin_MC][ybin_MC][i]++;
      //     fNDIS_evt_MC[1][xbin_MC][ybin_MC][i]++;
      //     fNDIS_evt_MC[2][xbin_MC][ybin_MC][i]++;
      //
      //     fFlag_MC[0][xbin_MC][ybin_MC][i]=0;
      //     fFlag_MC[1][xbin_MC][ybin_MC][i]=0;
      //     fFlag_MC[2][xbin_MC][ybin_MC][i]=0;
      //
      //     DIS_MC[0][i] = 1;
      //     DIS_MC[1][i] = 1;
      //     DIS_MC[2][i] = 1;
      //
      //     // nu cut
      //     if(!(fNu_min[0][i]<nu_MC && nu_MC<fNu_max[0][i]))
      //     {
      //       fFlag_MC[0][xbin_MC][ybin_MC][i]=1;
      //     }
      //     if(!(fNu_min[1][i]<nu_MC && nu_MC<fNu_max[1][i]))
      //     {
      //       fFlag_MC[1][xbin_MC][ybin_MC][i]=1;
      //     }
      //     if(!(fNu_min[2][i]<nu_MC && nu_MC<fNu_max[2][i]))
      //     {
      //       fFlag_MC[2][xbin_MC][ybin_MC][i]=1;
      //     }
      //     if(fFlag_MC[0][xbin_MC][ybin_MC][i])
      //     {
      //       fNDIS_evt_MC[0][xbin_MC][ybin_MC][i]--; DIS_MC[0][i] = 0;
      //     }
      //     if(fFlag_MC[1][xbin_MC][ybin_MC][i])
      //     {
      //       fNDIS_evt_MC[1][xbin_MC][ybin_MC][i]--; DIS_MC[1][i] = 0;
      //     }
      //     if(fFlag_MC[2][xbin_MC][ybin_MC][i])
      //     {
      //       fNDIS_evt_MC[2][xbin_MC][ybin_MC][i]--; DIS_MC[2][i] = 0;
      //     }
      //   }
      // }
      // else
      // {
      //   for(int i=0; i<12; i++)
      //   {
      //     fFlag_MC[0][xbin_MC][ybin_MC][i]=0;
      //     fFlag_MC[1][xbin_MC][ybin_MC][i]=0;
      //     fFlag_MC[2][xbin_MC][ybin_MC][i]=0;
      //
      //     DIS_MC[0][i] = 1;
      //     DIS_MC[1][i] = 1;
      //     DIS_MC[2][i] = 1;
      //
      //     // nu cut
      //     if(!(fNu_min[0][i]<nu_MC && nu_MC<fNu_max[0][i]))
      //     {
      //       fFlag_MC[0][xbin_MC][ybin_MC][i]=1;
      //     }
      //     if(!(fNu_min[1][i]<nu_MC && nu_MC<fNu_max[1][i]))
      //     {
      //       fFlag_MC[1][xbin_MC][ybin_MC][i]=1;
      //     }
      //     if(!(fNu_min[2][i]<nu_MC && nu_MC<fNu_max[2][i]))
      //     {
      //       fFlag_MC[2][xbin_MC][ybin_MC][i]=1;
      //     }
      //     if(fFlag_MC[0][xbin_MC][ybin_MC][i])
      //     {
      //       DIS_MC[0][i] = 0;
      //     }
      //     if(fFlag_MC[1][xbin_MC][ybin_MC][i])
      //     {
      //       DIS_MC[1][i] = 0;
      //     }
      //     if(fFlag_MC[2][xbin_MC][ybin_MC][i])
      //     {
      //       DIS_MC[2][i] = 0;
      //     }
      //   }
      // }

      // -----------------------------------------------------------------------
      //  Data -----------------------------------------------------------------
      // -----------------------------------------------------------------------

      // if(fAllDISflag)
      // {
      //   // z Binning
      //
      //   for(int i=0; i<12; i++)
      //   {
      //     fNDIS_evt[0][xbin][ybin][i]++;
      //     fNDIS_evt[1][xbin][ybin][i]++;
      //     fNDIS_evt[2][xbin][ybin][i]++;
      //
      //     fFlag[0][xbin][ybin][i]=0;
      //     fFlag[1][xbin][ybin][i]=0;
      //     fFlag[2][xbin][ybin][i]=0;
      //
      //     DIS_rec[0][i] = 1;
      //     DIS_rec[1][i] = 1;
      //     DIS_rec[2][i] = 1;
      //
      //     // nu cut
      //     if(!(fNu_min[0][i]<nu && nu<fNu_max[0][i]))
      //     {
      //       fFlag[0][xbin][ybin][i]=1;
      //     }
      //     if(!(fNu_min[1][i]<nu && nu<fNu_max[1][i]))
      //     {
      //       fFlag[1][xbin][ybin][i]=1;
      //     }
      //     if(!(fNu_min[2][i]<nu && nu<fNu_max[2][i]))
      //     {
      //       fFlag[2][xbin][ybin][i]=1;
      //     }
      //     if(fFlag[0][xbin][ybin][i])
      //     {
      //       fNDIS_evt[0][xbin][ybin][i]--; DIS_rec[0][i] = 0;
      //     }
      //     if(fFlag[1][xbin][ybin][i])
      //     {
      //       fNDIS_evt[1][xbin][ybin][i]--; DIS_rec[1][i] = 0;
      //     }
      //     if(fFlag[2][xbin][ybin][i])
      //     {
      //       fNDIS_evt[2][xbin][ybin][i]--; DIS_rec[2][i] = 0;
      //     }
      //     if(xbin==xbin_MC && ybin==ybin_MC)
      //     {
      //       if(DIS_rec[0][i] && DIS_MC[0][i])
      //       {
      //         fNDIS_evt_c[0][xbin][ybin][i] += 1;
      //       }
      //       if(DIS_rec[1][i] && DIS_MC[1][i])
      //       {
      //         fNDIS_evt_c[1][xbin][ybin][i] += 1;
      //       }
      //       if(DIS_rec[2][i] && DIS_MC[2][i])
      //       {
      //         fNDIS_evt_c[2][xbin][ybin][i] += 1;
      //       }
      //     }
      //   }
      // }

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

          if(0.1<zBj && (fId==8 || fId==9) && abs(ph->GetLeaf("Hadrons.ph")->GetValue(i))<1)
            fKinematicsMC[0][11]->Fill(abs(ph->GetLeaf("Hadrons.ph")->GetValue(i)));

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

          if(trig&2)
          {
            fKinematicsMC[0][3]->Fill(zBj);
          }
          if(trig&4)
          {
            fKinematicsMC[1][3]->Fill(zBj);
          }
          if(trig&8)
          {
            fKinematicsMC[2][3]->Fill(zBj);
          }
          if(trig&512)
          {
            fKinematicsMC[3][3]->Fill(zBj);
          }
          if(trig&2 || trig&4 || trig&8)
          // if(trig&2 || trig&4 || trig&8 || trig&512)
          {
            fKinematicsMC[4][3]->Fill(zBj);
          }
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

    cout << "\n-> Finished processing file " << filename << " <-\n" << endl;

    delete f;
  }

  for(int i=0; i<int(fQ2kinMC[0].size()); i++)
  {
    fKinematicsMC[0][0]->Fill(fQ2kinMC[0][i]);
    fKinematicsMC[0][1]->Fill(fXBjkinMC[0][i]);
    fKinematicsMC[0][2]->Fill(fYBjkinMC[0][i]);
    fKinematicsMC[0][4]->Fill(fWBjkinMC[0][i]);
    fKinematicsMC[0][5]->Fill(fNukinMC[0][i]);
    fKinematicsMC[0][6]->Fill(fMuMC[0][i]);
    fKinematicsMC[0][7]->Fill(fMupMC[0][i]);
    fKinematicsMC[0][8]->Fill(fThetaMC[0][i]);
    fKinematicsMC[0][9]->Fill(fPhiMC[0][i]);
    fKinematicsMC[0][10]->Fill(fVertexMC[0][i]);
  }
  for(int i=0; i<int(fQ2kinMC[1].size()); i++)
  {
    fKinematicsMC[1][0]->Fill(fQ2kinMC[1][i]);
    fKinematicsMC[1][1]->Fill(fXBjkinMC[1][i]);
    fKinematicsMC[1][2]->Fill(fYBjkinMC[1][i]);
    fKinematicsMC[1][4]->Fill(fWBjkinMC[1][i]);
    fKinematicsMC[1][5]->Fill(fNukinMC[1][i]);
    fKinematicsMC[1][6]->Fill(fMuMC[1][i]);
    fKinematicsMC[1][7]->Fill(fMupMC[1][i]);
    fKinematicsMC[1][8]->Fill(fThetaMC[1][i]);
    fKinematicsMC[1][9]->Fill(fPhiMC[1][i]);
    fKinematicsMC[1][10]->Fill(fVertexMC[1][i]);
  }
  for(int i=0; i<int(fQ2kinMC[2].size()); i++)
  {
    fKinematicsMC[2][0]->Fill(fQ2kinMC[2][i]);
    fKinematicsMC[2][1]->Fill(fXBjkinMC[2][i]);
    fKinematicsMC[2][2]->Fill(fYBjkinMC[2][i]);
    fKinematicsMC[2][4]->Fill(fWBjkinMC[2][i]);
    fKinematicsMC[2][5]->Fill(fNukinMC[2][i]);
    fKinematicsMC[2][6]->Fill(fMuMC[2][i]);
    fKinematicsMC[2][7]->Fill(fMupMC[2][i]);
    fKinematicsMC[2][8]->Fill(fThetaMC[2][i]);
    fKinematicsMC[2][9]->Fill(fPhiMC[2][i]);
    fKinematicsMC[2][10]->Fill(fVertexMC[2][i]);
  }
  for(int i=0; i<int(fQ2kinMC[3].size()); i++)
  {
    fKinematicsMC[3][0]->Fill(fQ2kinMC[3][i]);
    fKinematicsMC[3][1]->Fill(fXBjkinMC[3][i]);
    fKinematicsMC[3][2]->Fill(fYBjkinMC[3][i]);
    fKinematicsMC[3][4]->Fill(fWBjkinMC[3][i]);
    fKinematicsMC[3][5]->Fill(fNukinMC[3][i]);
    fKinematicsMC[3][6]->Fill(fMuMC[3][i]);
    fKinematicsMC[3][7]->Fill(fMupMC[3][i]);
    fKinematicsMC[3][8]->Fill(fThetaMC[3][i]);
    fKinematicsMC[3][9]->Fill(fPhiMC[3][i]);
    fKinematicsMC[3][10]->Fill(fVertexMC[3][i]);
  }
  for(int i=0; i<int(fQ2kinMC[4].size()); i++)
  {
    fKinematicsMC[4][0]->Fill(fQ2kinMC[4][i]);
    fKinematicsMC[4][1]->Fill(fXBjkinMC[4][i]);
    fKinematicsMC[4][2]->Fill(fYBjkinMC[4][i]);
    fKinematicsMC[4][4]->Fill(fWBjkinMC[4][i]);
    fKinematicsMC[4][5]->Fill(fNukinMC[4][i]);
    fKinematicsMC[4][6]->Fill(fMuMC[4][i]);
    fKinematicsMC[4][7]->Fill(fMupMC[4][i]);
    fKinematicsMC[4][8]->Fill(fThetaMC[4][i]);
    fKinematicsMC[4][9]->Fill(fPhiMC[4][i]);
    fKinematicsMC[4][10]->Fill(fVertexMC[4][i]);
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
      trig = (trig&2047);

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

      fAllDISflag = 0;
      int DIS_rec[3][12];

      // Best Primary Vertex
      fBP++;

      // Reconstructed muon
      if((0<E_beam->GetLeaf("E_beam")->GetValue()))
      {
        fRmu++;

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

                  if((trig&2 || trig&4 || trig&8))
                  // if((trig&2 || trig&4 || trig&8 || trig&512))
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
            //2016 ---
          }
        }
      }

      if(fAllDISflag)
      {
        double theta_m = asin(sqrt(pow(p1x->GetLeaf("p1x")->GetValue()/sqrt(pow(E_mu_prim->GetLeaf("E_mu_prim")->GetValue(),2)-pow(fM_mu,2)),2)+pow(p1y->GetLeaf("p1y")->GetValue()/sqrt(pow(E_mu_prim->GetLeaf("E_mu_prim")->GetValue(),2)-pow(fM_mu,2)),2)));
        double phi_m = asin(p1x->GetLeaf("p1x")->GetValue()/sqrt(pow(p1x->GetLeaf("p1x")->GetValue(),2)+pow(p1y->GetLeaf("p1y")->GetValue(),2)));

        if(trig&2)
        {
          fQ2kin[0].push_back(Q2);
          fXBjkin[0].push_back(xBj);
          fYBjkin[0].push_back(yBj);
          fWBjkin[0].push_back(sqrt(wBj));
          fNukin[0].push_back(nu);
          fMu[0].push_back(E_beam->GetLeaf("E_beam")->GetValue());
          fMup[0].push_back(E_mu_prim->GetLeaf("E_mu_prim")->GetValue());
          fTheta[0].push_back(theta_m);
          fPhi[0].push_back(phi_m);
          fVertex[0].push_back(z->GetLeaf("z")->GetValue());
        }
        if(trig&4)
        {
          fQ2kin[1].push_back(Q2);
          fXBjkin[1].push_back(xBj);
          fYBjkin[1].push_back(yBj);
          fWBjkin[1].push_back(sqrt(wBj));
          fNukin[1].push_back(nu);
          fMu[1].push_back(E_beam->GetLeaf("E_beam")->GetValue());
          fMup[1].push_back(E_mu_prim->GetLeaf("E_mu_prim")->GetValue());
          fTheta[1].push_back(theta_m);
          fPhi[1].push_back(phi_m);
          fVertex[1].push_back(z->GetLeaf("z")->GetValue());
        }
        if(trig&8)
        {
          fQ2kin[2].push_back(Q2);
          fXBjkin[2].push_back(xBj);
          fYBjkin[2].push_back(yBj);
          fWBjkin[2].push_back(sqrt(wBj));
          fNukin[2].push_back(nu);
          fMu[2].push_back(E_beam->GetLeaf("E_beam")->GetValue());
          fMup[2].push_back(E_mu_prim->GetLeaf("E_mu_prim")->GetValue());
          fTheta[2].push_back(theta_m);
          fPhi[2].push_back(phi_m);
          fVertex[2].push_back(z->GetLeaf("z")->GetValue());
        }
        if(trig&512)
        {
          fQ2kin[3].push_back(Q2);
          fXBjkin[3].push_back(xBj);
          fYBjkin[3].push_back(yBj);
          fWBjkin[3].push_back(sqrt(wBj));
          fNukin[3].push_back(nu);
          fMu[3].push_back(E_beam->GetLeaf("E_beam")->GetValue());
          fMup[3].push_back(E_mu_prim->GetLeaf("E_mu_prim")->GetValue());
          fTheta[3].push_back(theta_m);
          fPhi[3].push_back(phi_m);
          fVertex[3].push_back(z->GetLeaf("z")->GetValue());
        }
        if(trig&2 || trig&4 || trig&8)
        // if(trig&2 || trig&4 || trig&8 || trig&512)
        {
          fQ2kin[4].push_back(Q2);
          fXBjkin[4].push_back(xBj);
          fYBjkin[4].push_back(yBj);
          fWBjkin[4].push_back(sqrt(wBj));
          fNukin[4].push_back(nu);
          fMu[4].push_back(E_beam->GetLeaf("E_beam")->GetValue());
          fMup[4].push_back(E_mu_prim->GetLeaf("E_mu_prim")->GetValue());
          fTheta[4].push_back(theta_m);
          fPhi[4].push_back(phi_m);
          fVertex[4].push_back(z->GetLeaf("z")->GetValue());
        }
        fX.push_back(x->GetLeaf("x")->GetValue());
        fY.push_back(y->GetLeaf("y")->GetValue());
      }


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

//       for(int i=0; i<12; i++)
//       {
// #ifdef NORC
//         fNDIS_evt[0][xbin][ybin][i] += 1;
//         fNDIS_evt[1][xbin][ybin][i] += 1;
//         fNDIS_evt[2][xbin][ybin][i] += 1;
//
//         fNDIS_evt_err[0][xbin][ybin][i] += 1;
//         fNDIS_evt_err[1][xbin][ybin][i] += 1;
//         fNDIS_evt_err[2][xbin][ybin][i] += 1;
// #else
//         fNDIS_evt[0][xbin][ybin][i] += 1*PaAlgo::GetRadiativeWeightLiD(xBj,yBj,1,z->GetLeaf("z")->GetValue());
//         fNDIS_evt[1][xbin][ybin][i] += 1*PaAlgo::GetRadiativeWeightLiD(xBj,yBj,1,z->GetLeaf("z")->GetValue());
//         fNDIS_evt[2][xbin][ybin][i] += 1*PaAlgo::GetRadiativeWeightLiD(xBj,yBj,1,z->GetLeaf("z")->GetValue());
//
//         fNDIS_evt_err[0][xbin][ybin][i] += pow(PaAlgo::GetRadiativeWeightLiD(xBj,yBj,1,z->GetLeaf("z")->GetValue()),2);
//         fNDIS_evt_err[1][xbin][ybin][i] += pow(PaAlgo::GetRadiativeWeightLiD(xBj,yBj,1,z->GetLeaf("z")->GetValue()),2);
//         fNDIS_evt_err[2][xbin][ybin][i] += pow(PaAlgo::GetRadiativeWeightLiD(xBj,yBj,1,z->GetLeaf("z")->GetValue()),2);
// #endif
//
//         fFlag[0][xbin][ybin][i]=0;
//         fFlag[1][xbin][ybin][i]=0;
//         fFlag[2][xbin][ybin][i]=0;
//
//         // nu cut
//         if(!(fNu_min[0][i]<nu && nu<fNu_max[0][i]))
//         {
//           fFlag[0][xbin][ybin][i]=1;
//         }
//         if(!(fNu_min[1][i]<nu && nu<fNu_max[1][i]))
//         {
//           fFlag[1][xbin][ybin][i]=1;
//         }
//         if(!(fNu_min[2][i]<nu && nu<fNu_max[2][i]))
//         {
//           fFlag[2][xbin][ybin][i]=1;
//         }
//         if(fFlag[0][xbin][ybin][i] /*|| fFlag[1][xbin][ybin][i] || fFlag[2][xbin][ybin][i]*/)
//         {
// #ifdef NORC
//           fNDIS_evt[0][xbin][ybin][i] -= 1;
//         // fNDIS_evt[1][xbin][ybin][i] -= 1;
//         // fNDIS_evt[2][xbin][ybin][i] -= 1;
//           fNDIS_evt_err[0][xbin][ybin][i] -= 1;
//         // fNDIS_evt_err[1][xbin][ybin][i] -= 1;
//         // fNDIS_evt_err[2][xbin][ybin][i] -= 1;
// #else
//           fNDIS_evt[0][xbin][ybin][i] -= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,1,z->GetLeaf("z")->GetValue());
//          // fNDIS_evt[1][xbin][ybin][i] -= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,1,z->GetLeaf("z")->GetValue());
//          // fNDIS_evt[2][xbin][ybin][i] -= PaAlgo::GetRadiativeWeightLiD(xBj,yBj,1,z->GetLeaf("z")->GetValue());
//           fNDIS_evt_err[0][xbin][ybin][i] -= pow(PaAlgo::GetRadiativeWeightLiD(xBj,yBj,1,z->GetLeaf("z")->GetValue()),2);
//          // fNDIS_evt_err[1][xbin][ybin][i] -= pow(PaAlgo::GetRadiativeWeightLiD(xBj,yBj,1,z->GetLeaf("z")->GetValue()),2);
//          // fNDIS_evt_err[2][xbin][ybin][i] -= pow(PaAlgo::GetRadiativeWeightLiD(xBj,yBj,1,z->GetLeaf("z")->GetValue()),2);
// #endif
//         }
//       }

      // -------------------------------------------------------------------------
      // --------- Hadrons Selection ---------------------------------------------
      // -------------------------------------------------------------------------

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
          }
          else
          {
            zBj = 0;
          }

          if(0.1<zBj && (fId==8 || fId==9) && abs(ph->GetLeaf("Hadrons.ph")->GetValue(i))<1)
            fKinematicsMC[0][11]->Fill(abs(ph->GetLeaf("Hadrons.ph")->GetValue(i)));

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


          if(trig&2)
          {
            fKinematicsRD[0][3]->Fill(zBj);
          }
          if(trig&4)
          {
            fKinematicsRD[1][3]->Fill(zBj);
          }
          if(trig&8)
          {
            fKinematicsRD[2][3]->Fill(zBj);
          }
          if(trig&512)
          {
            fKinematicsRD[3][3]->Fill(zBj);
          }
          if(trig&2 || trig&4 || trig&8)
          // if(trig&2 || trig&4 || trig&8 || trig&512)
          {
            fKinematicsRD[4][3]->Fill(zBj);
          }
        }
      }
    }

    cout << "\n-> Finished processing file " << filename << " <-\n" << endl;

    delete f;
  }

  for(int i=0; i<int(fQ2kin[0].size()); i++)
  {
      fKinematicsRD[0][0]->Fill(fQ2kin[0][i]);
      fKinematicsRD[0][1]->Fill(fXBjkin[0][i]);
      fKinematicsRD[0][2]->Fill(fYBjkin[0][i]);
      fKinematicsRD[0][4]->Fill(fWBjkin[0][i]);
      fKinematicsRD[0][5]->Fill(fNukin[0][i]);
      fKinematicsRD[0][6]->Fill(fMu[0][i]);
      fKinematicsRD[0][7]->Fill(fMup[0][i]);
      fKinematicsRD[0][8]->Fill(fTheta[0][i]);
      fKinematicsRD[0][9]->Fill(fPhi[0][i]);
      fKinematicsRD[0][10]->Fill(fVertex[0][i]);
  }
  for(int i=0; i<int(fQ2kin[1].size()); i++)
  {
      fKinematicsRD[1][0]->Fill(fQ2kin[1][i]);
      fKinematicsRD[1][1]->Fill(fXBjkin[1][i]);
      fKinematicsRD[1][2]->Fill(fYBjkin[1][i]);
      fKinematicsRD[1][4]->Fill(fWBjkin[1][i]);
      fKinematicsRD[1][5]->Fill(fNukin[1][i]);
      fKinematicsRD[1][6]->Fill(fMu[1][i]);
      fKinematicsRD[1][7]->Fill(fMup[1][i]);
      fKinematicsRD[1][8]->Fill(fTheta[1][i]);
      fKinematicsRD[1][9]->Fill(fPhi[1][i]);
      fKinematicsRD[1][10]->Fill(fVertex[1][i]);
  }
  for(int i=0; i<int(fQ2kin[2].size()); i++)
  {
      fKinematicsRD[2][0]->Fill(fQ2kin[2][i]);
      fKinematicsRD[2][1]->Fill(fXBjkin[2][i]);
      fKinematicsRD[2][2]->Fill(fYBjkin[2][i]);
      fKinematicsRD[2][4]->Fill(fWBjkin[2][i]);
      fKinematicsRD[2][5]->Fill(fNukin[2][i]);
      fKinematicsRD[2][6]->Fill(fMu[2][i]);
      fKinematicsRD[2][7]->Fill(fMup[2][i]);
      fKinematicsRD[2][8]->Fill(fTheta[2][i]);
      fKinematicsRD[2][9]->Fill(fPhi[2][i]);
      fKinematicsRD[2][10]->Fill(fVertex[2][i]);
  }
  for(int i=0; i<int(fQ2kin[3].size()); i++)
  {
      fKinematicsRD[3][0]->Fill(fQ2kin[3][i]);
      fKinematicsRD[3][1]->Fill(fXBjkin[3][i]);
      fKinematicsRD[3][2]->Fill(fYBjkin[3][i]);
      fKinematicsRD[3][4]->Fill(fWBjkin[3][i]);
      fKinematicsRD[3][5]->Fill(fNukin[3][i]);
      fKinematicsRD[3][6]->Fill(fMu[3][i]);
      fKinematicsRD[3][7]->Fill(fMup[3][i]);
      fKinematicsRD[3][8]->Fill(fTheta[3][i]);
      fKinematicsRD[3][9]->Fill(fPhi[3][i]);
      fKinematicsRD[3][10]->Fill(fVertex[3][i]);
  }
  for(int i=0; i<int(fQ2kin[4].size()); i++)
  {
      fKinematicsRD[4][0]->Fill(fQ2kin[4][i]);
      fKinematicsRD[4][1]->Fill(fXBjkin[4][i]);
      fKinematicsRD[4][2]->Fill(fYBjkin[4][i]);
      fKinematicsRD[4][4]->Fill(fWBjkin[4][i]);
      fKinematicsRD[4][5]->Fill(fNukin[4][i]);
      fKinematicsRD[4][6]->Fill(fMu[4][i]);
      fKinematicsRD[4][7]->Fill(fMup[4][i]);
      fKinematicsRD[4][8]->Fill(fTheta[4][i]);
      fKinematicsRD[4][9]->Fill(fPhi[4][i]);
      fKinematicsRD[4][10]->Fill(fVertex[4][i]);
  }

}

int main(int argc, char **argv)
{

  if(argc < 2)
  {
    cout << "ERROR : Not enough arguments." << endl;
    cout << "Asked : 2 *** Received : " << argc-1 << endl;
    cout << "./compMCMC [RD filelist] [MC filelist]" << endl;

    return 1;
  }

  create_kin_plots();
  RDextraction(argv[1]);
  MCextraction(argv[2]);
  save_kin_plots();

  return 0;
}
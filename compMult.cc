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

#include "compMult.h"

using namespace std;

void LoadMultiplicityFiles(string pfile1, string pfile2)
{
  string sdum;
  Double_t x,y,z;

  ifstream mult1(pfile1);

  for(int i=0; i<9; i++)
  {
    for(int j=0; j<6; j++)
    {
      for(int k=0; k<12; k++)
      {
          mult1 >> x >> y >> z;
          // cout << x << "\t" << y << "\t" << z << "\t" << endl;
          for(int l=0; l<4; l++) mult1 >> sdum;
          mult1 >> fMultiplicities1[i][j][k].tab[1][0][0] >> fMultiplicities1[i][j][k].tab[1][1][0] >> fMultiplicities1[i][j][k].tab[1][2][0];
          for(int l=0; l<5; l++) mult1 >> sdum;
          mult1 >> fMultiplicities1[i][j][k].tab[0][0][0] >> fMultiplicities1[i][j][k].tab[0][1][0] >> fMultiplicities1[i][j][k].tab[0][2][0];
          mult1 >> sdum;
      }
    }
  }
  mult1.close();

  ifstream mult2(pfile2);

  for(int i=0; i<9; i++)
  {
    for(int j=0; j<6; j++)
    {
      for(int k=0; k<12; k++)
      {
          mult2 >> x >> y >> z;
          // cout << x << "\t" << y << "\t" << z << "\t" << endl;
          for(int l=0; l<4; l++) mult2 >> sdum;
          mult2 >> fMultiplicities2[i][j][k].tab[1][0][0] >> fMultiplicities2[i][j][k].tab[1][1][0] >> fMultiplicities2[i][j][k].tab[1][2][0];
          for(int l=0; l<5; l++) mult2 >> sdum;
          mult2 >> fMultiplicities2[i][j][k].tab[0][0][0] >> fMultiplicities2[i][j][k].tab[0][1][0] >> fMultiplicities2[i][j][k].tab[0][2][0];
          mult2 >> sdum;
      }
    }
  }
  mult2.close();
}

void yweightedavg()
{
  for(int c=0; c<2; c++)
  {
    for(int x=0; x<9; x++)
    {
      for(int z=0; z<12; z++)
      {
        for(int i=0; i<6; i++)
        {
          if(fMultiplicities1[x][i][z].tab[c][0][0])
          {
            fMultiplicities1_yavg[x][z].tab[c][0][0]+=fMultiplicities1[x][i][z].tab[c][0][0]/fMultiplicities1[x][i][z].tab[c][1][0];
            fMultiplicities1_yavg[x][z].tab[c][1][0]+=1/fMultiplicities1[x][i][z].tab[c][1][0];
            fMultiplicities1_yavg[x][z].tab[c][2][0]+=1/fMultiplicities1[x][i][z].tab[c][2][0];
          }
          if(fMultiplicities2[x][i][z].tab[c][0][0])
          {
            fMultiplicities2_yavg[x][z].tab[c][0][0]+=fMultiplicities2[x][i][z].tab[c][0][0]/fMultiplicities2[x][i][z].tab[c][1][0];
            fMultiplicities2_yavg[x][z].tab[c][1][0]+=1/fMultiplicities2[x][i][z].tab[c][1][0];
            fMultiplicities2_yavg[x][z].tab[c][2][0]+=1/fMultiplicities2[x][i][z].tab[c][2][0];
          }
        }
        if(fMultiplicities1_yavg[x][z].tab[c][0][0])
        {
          fMultiplicities1_yavg[x][z].tab[c][1][0]=1/fMultiplicities1_yavg[x][z].tab[c][1][0];
          fMultiplicities1_yavg[x][z].tab[c][2][0]=1/fMultiplicities1_yavg[x][z].tab[c][2][0];
          fMultiplicities1_yavg[x][z].tab[c][0][0]*=fMultiplicities1_yavg[x][z].tab[c][1][0];
        }
        if(fMultiplicities2_yavg[x][z].tab[c][0][0])
        {
          fMultiplicities2_yavg[x][z].tab[c][1][0]=1/fMultiplicities2_yavg[x][z].tab[c][1][0];
          fMultiplicities2_yavg[x][z].tab[c][2][0]=1/fMultiplicities2_yavg[x][z].tab[c][2][0];
          fMultiplicities2_yavg[x][z].tab[c][0][0]*=fMultiplicities2_yavg[x][z].tab[c][1][0];
        }
      }
    }
  }
}

int main(int argc, char **argv)
{

  if(argc < 3)
  {
    cout << "ERROR : Not enough arguments." << endl;
    cout << "Asked : 2 *** Received : " << argc-1 << endl;
    cout << "./compMult Mult1 Mult2" << endl;

    return 1;
  }

  LoadMultiplicityFiles(argv[1],argv[2]);

  TCanvas c1("Mult_comparison","Mult_comparison",3200,1600);
  TCanvas c2("Mult_comparison_yavg","Mult_comparison_yavg",3200,1600);
  TLine l1(0.1,0.9,0.9,0.9);
  TLine l2(0.1,1.1,0.9,1.1);
  TLine l3(0.1,0.95,0.9,0.95);
  TLine l4(0.1,1.05,0.9,1.05);
  l3.SetLineStyle(2); l4.SetLineStyle(2);
  c1.SetFillColor(0);
  c2.SetFillColor(0);
  c1.Divide(9,5,0,0);
  c2.Divide(5,2,0,0);

  ofstream ofs_ra(Form("%s/reldiff.txt",data_path), std::ofstream::out | std::ofstream::trunc);

  TGraphErrors* R[2][9][6];
  TGraphErrors* R_y[2][9];

  for(int c=0; c<2; c++)
  {
    for(int i=0; i<9; i++)
    {

      vector<double> r_y;
      vector<double> r_y_err;
      vector<double> z_range_r_y;

      for(int l=0; l<12; l++)
      {
        z_range_r_y.push_back(z_range[l]);
      }

      for(int j=0; j<6; j++)
      {
        std::vector<double> r;
        std::vector<double> r_err;
        std::vector<double> z_range_r;

        for(int l=0; l<12; l++)
        {
          z_range_r.push_back(z_range[l]);
        }

        for(int k=0; k<12; k++)
        {
          r.push_back(fMultiplicities2[i][j][k].tab[c][0][0] ? fMultiplicities1[i][j][k].tab[c][0][0]/fMultiplicities2[i][j][k].tab[c][0][0] : 0);
          r_err.push_back(sqrt((fMultiplicities1[i][j][k].tab[c][1][0]+pow(fMultiplicities2[i][j][k].tab[c][1][0],2)*fMultiplicities1[i][j][k].tab[c][0][0]
                                  /pow(fMultiplicities2[i][j][k].tab[c][0][0],2))/pow(fMultiplicities2[i][j][k].tab[c][0][0],2)));
        }

        for(int k=12; k>0; k--)
        {
          if(!r[k-1]) {r.erase(r.begin()+k-1); r_err.erase(r_err.begin()+k-1); z_range_r.erase(z_range_r.begin()+k-1);}
        }

        bool r_empty = 0;

        if(!(int(r.size()))) r_empty = 1;

        R[c][i][j] = new TGraphErrors(int(r.size()),&(z_range_r[0]),&(r[0]),0,&(r_err[0]));

        if(!c)
        {
          R[c][i][j]->SetMarkerColor(fMarkerColor[4]);
        }
        else
        {
          R[c][i][j]->SetMarkerColor(fMarkerColor[0]);
        }

        R[c][i][j]->SetMarkerSize(1);

        R[c][i][j]->SetMarkerStyle(fMarkerStyle[0][c]);
        R[c][i][j]->GetYaxis()->SetTitle("");

        R[c][i][j]->GetXaxis()->SetTitle("");

        R[c][i][j]->SetTitle("");

        if(!r_empty)
        {
          c1.cd(i+1+9*j);
          gPad->SetFillStyle(4000);
          if(R[c][i][j])
          {
            if(!c)
            {
              R[c][i][j]->Draw("SAMEPA");
              R[c][i][j]->GetXaxis()->SetLimits(0.1,0.9);
              R[c][i][j]->SetMinimum(0.);
              R[c][i][j]->SetMaximum(2.);
              R[c][i][j]->GetXaxis()->SetLabelSize(0.06);
              R[c][i][j]->GetYaxis()->SetLabelSize(0.06);
              R[c][i][j]->SetTitle("");
              if(j==5) gPad->SetBottomMargin(.15);
              if(i==0) gPad->SetLeftMargin(.22);
              if(i==8 && j==5)
              {
                R[c][i][j]->GetXaxis()->SetTitle("#font[ 12]{z}");
                R[c][i][j]->GetXaxis()->SetTitleSize(0.08);
                R[c][i][j]->GetXaxis()->SetTitleOffset(.8);
              }
              R[c][i][j]->GetXaxis()->SetNdivisions(304,kTRUE);
              R[c][i][j]->GetYaxis()->SetNdivisions(304,kTRUE);
              if(i==1 && j==0)
              {
                R[c][i][j]->GetYaxis()->SetTitle("#font[12]{acceptance}^{#font[ 12]{h}}");
                R[c][i][j]->GetYaxis()->SetTitleSize(0.08);
              }
              R[c][i][j]->Draw("SAMEP");
              R[c][i][j]->GetXaxis()->SetLimits(0.1,0.9);
              R[c][i][j]->SetMinimum(0.);
              R[c][i][j]->SetMaximum(2.);
              l1.Draw("SAME");
              l2.Draw("SAME");
              l3.Draw("SAME");
              l4.Draw("SAME");
              c1.Range(0.1,0.,0.9,2.);
            }
            else
            {
              R[c][i][j]->Draw("SAMEP");
              R[c][i][j]->GetXaxis()->SetLimits(0.1,0.9);
              R[c][i][j]->SetMinimum(0.);
              R[c][i][j]->SetMaximum(2.);
            }
          }
          c1.Update();
        }
      }

      yweightedavg();

      for(int k=0; k<12; k++)
      {
        double multy, multye;
        multy = fMultiplicities2_yavg[i][k].tab[c][0][0] ? fMultiplicities1_yavg[i][k].tab[c][0][0]/fMultiplicities2_yavg[i][k].tab[c][0][0] : 0;
        multye = fMultiplicities2_yavg[i][k].tab[c][0][0] ? sqrt((fMultiplicities1_yavg[i][k].tab[c][1][0]+pow(fMultiplicities2_yavg[i][k].tab[c][1][0],2)*fMultiplicities1_yavg[i][k].tab[c][0][0]
                                /pow(fMultiplicities2_yavg[i][k].tab[c][0][0],2))/pow(fMultiplicities2_yavg[i][k].tab[c][0][0],2)) : 0;
        r_y.push_back(multy);
        r_y_err.push_back(multye);
        int ratFlag;
        if(multy>=1)
        {
          if(multy+multye<=1) ratFlag = 1;
          else ratFlag = 0;
        }
        else
        {
          if(multy+multye>=1) ratFlag = 1;
          else ratFlag = 0;
        }
        ofs_ra << c << " " << i << " " << k << " " << multy << " " << multye << " " << ratFlag << endl;
      }
      for(int k=12; k>0; k--)
      {
        if(!r_y[k-1]) {r_y.erase(r_y.begin()+k-1); r_y_err.erase(r_y_err.begin()+k-1); z_range_r_y.erase(z_range_r_y.begin()+k-1);}
      }

      bool r_y_empty = 0;

      if(!(int(r_y.size()))) r_y_empty = 1;

      R_y[c][i] = new TGraphErrors(int(r_y.size()),&(z_range_r_y[0]),&(r_y[0]),0,&(r_y_err[0]));

      if(!c)
      {
        R_y[c][i]->SetMarkerColor(fMarkerColor[4]);
      }
      else
      {
        R_y[c][i]->SetMarkerColor(fMarkerColor[0]);
      }

      R_y[c][i]->SetMarkerSize(3);

      R_y[c][i]->SetMarkerStyle(fMarkerStyle[0][c]);
      R_y[c][i]->GetYaxis()->SetTitle("");

      R_y[c][i]->GetXaxis()->SetTitle("");

      R_y[c][i]->SetTitle("");

      if(!r_y_empty)
      {
        c2.cd(i+1);
        gPad->SetFillStyle(4000);
        if(R_y[c][i])
        {
          if(!c)
          {
            R_y[c][i]->Draw("SAMEPA");
            R_y[c][i]->GetXaxis()->SetLimits(0.1,0.9);
            R_y[c][i]->SetMinimum(0.);
            R_y[c][i]->SetMaximum(2.);
            R_y[c][i]->GetXaxis()->SetLabelSize(0.06);
            R_y[c][i]->GetYaxis()->SetLabelSize(0.06);
            R_y[c][i]->SetTitle("");
            if(i>4) gPad->SetBottomMargin(.15);
            if(i==0 || i==5) gPad->SetLeftMargin(.22);
            if(i==8)
            {
              R_y[c][i]->GetXaxis()->SetTitle("#font[ 12]{z}");
              R_y[c][i]->GetXaxis()->SetTitleSize(0.08);
              R_y[c][i]->GetXaxis()->SetTitleOffset(.8);
            }
            R_y[c][i]->GetXaxis()->SetNdivisions(304,kTRUE);
            R_y[c][i]->GetYaxis()->SetNdivisions(304,kTRUE);
            if(i==0)
            {
              R_y[c][i]->GetYaxis()->SetTitle("#font[12]{M}^{#font[ 12]{h}}_{#font[ 12]{ratio}}");
              R_y[c][i]->GetYaxis()->SetTitleSize(0.08);
            }
            R_y[c][i]->Draw("SAMEP");
            R_y[c][i]->GetXaxis()->SetLimits(0.1,0.9);
            R_y[c][i]->SetMinimum(0.);
            R_y[c][i]->SetMaximum(2.);
            l1.Draw("SAME");
            l2.Draw("SAME");
            l3.Draw("SAME");
            l4.Draw("SAME");
            c2.Range(0.1,0.,0.9,2.);
          }
          else
          {
            R_y[c][i]->Draw("SAMEP");
            R_y[c][i]->GetXaxis()->SetLimits(0.1,0.9);
            R_y[c][i]->SetMinimum(0.);
            R_y[c][i]->SetMaximum(2.0);
          }
        }
        c2.Update();
      }
    }
  }

  c1.Print("mult_ratio.pdf");
  c2.Print("mult_ratio_yavg.pdf");

  ofs_ra.close();

  return 0;
}

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

#include "plotMult.h"

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
          mult1 >> x;
          for(int l=0; l<2; l++) mult1 >> sdum;
          mult1 >> y;
          for(int l=0; l<3; l++) mult1 >> sdum;
          mult1 >> z;
          // cout << x << "\t" << y << "\t" << z << "\t" << endl;
          for(int l=0; l<2; l++) mult1 >> sdum;
          mult1 >> fMultiplicities[i][j][k].tab[1][0][0] >> fMultiplicities[i][j][k].tab[1][1][0] >> sdum >> fMultiplicities[i][j][k].tab[1][2][0] >> sdum;
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
        mult2 >> x;
        for(int l=0; l<2; l++) mult2 >> sdum;
        mult2 >> y;
        for(int l=0; l<3; l++) mult2 >> sdum;
        mult2 >> z;
        // cout << x << "\t" << y << "\t" << z << "\t" << endl;
        for(int l=0; l<2; l++) mult2 >> sdum;
        mult2 >> fMultiplicities[i][j][k].tab[0][0][0] >> fMultiplicities[i][j][k].tab[0][1][0] >> sdum >> fMultiplicities[i][j][k].tab[0][2][0] >> sdum;
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
          if(fMultiplicities[x][i][z].tab[c][0][0])
          {
            fMultiplicities_yavg[x][z].tab[c][0][0]+=fMultiplicities[x][i][z].tab[c][0][0]/fMultiplicities[x][i][z].tab[c][1][0];
            fMultiplicities_yavg[x][z].tab[c][1][0]+=1/fMultiplicities[x][i][z].tab[c][1][0];
            fMultiplicities_yavg[x][z].tab[c][2][0]+=1/fMultiplicities[x][i][z].tab[c][2][0];
          }
        }
        if(fMultiplicities_yavg[x][z].tab[c][0][0])
        {
          fMultiplicities_yavg[x][z].tab[c][1][0]=1/fMultiplicities_yavg[x][z].tab[c][1][0];
          fMultiplicities_yavg[x][z].tab[c][2][0]=1/fMultiplicities_yavg[x][z].tab[c][2][0];
          fMultiplicities_yavg[x][z].tab[c][0][0]*=fMultiplicities_yavg[x][z].tab[c][1][0];
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
    cout << "./compMult Mult+ Mult-" << endl;

    return 1;
  }

  LoadMultiplicityFiles(argv[1],argv[2]);

  TCanvas c1("Mult+","Mult+",3200,1600);
  TCanvas c2("Mult-","Mult-",3200,1600);
  TCanvas c3("Mult_yavg","Mult_yavg",3200,1600);
  c1.SetFillColor(0);
  c2.SetFillColor(0);
  c3.SetFillColor(0);
  c1.Divide(5,2,0,0);
  c2.Divide(5,2,0,0);
  c3.Divide(5,2,0,0);

  Double_t errorx[12] = {0.05/2,0.05/2,0.05/2,0.05/2,0.05/2,0.05/2,0.05/2,0.05/2,0.05/2,0.05/2,0.05/2,0.1/2};
  Double_t h_yoffset[12] = {-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2};
  Double_t h_yoffset2[12] = {-0.3,-0.3,-0.3,-0.3,-0.3,-0.3,-0.3,-0.3,-0.3,-0.3,-0.3,-0.3};

  TGraphErrors* mult[2][9][6];
  TGraphAsymmErrors* sys[2][9][6];
  TGraphErrors* mult_y[2][9];

  TLine lsys(0.1,0,0.9,0);
  lsys.SetLineStyle(2);

  for(int i=0; i<9; i++)
  {
    int axisflagh1 = 0;
    int axisflagh2 = 0;

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
      std::vector<double> r_sys;
      std::vector<double> z_range_r;

      for(int c=1; c>=0; c--)
      {
        for(int l=0; l<12; l++)
        {
          z_range_r.push_back(z_range[l]);
        }

        for(int k=0; k<12; k++)
        {
          r.push_back(fMultiplicities[i][j][k].tab[c][0][0]);
          r_err.push_back(fMultiplicities[i][j][k].tab[c][1][0]);
          r_sys.push_back(fMultiplicities[i][j][k].tab[c][2][0]);
        }

        for(int k=12; k>0; k--)
        {
          if(!r[k-1]) {r.erase(r.begin()+k-1); r_err.erase(r_err.begin()+k-1); r_sys.erase(r_sys.begin()+k-1); z_range_r.erase(z_range_r.begin()+k-1);}
        }

        bool r_empty = 0;

        if(!(int(r.size()))) r_empty = 1;

        mult[c][i][j] = new TGraphErrors(int(r.size()),&(z_range_r[0]),&(r[0]),0,&(r_err[0]));
        sys[c][i][j] = new TGraphAsymmErrors(Int_t(r.size()),&(z_range_r[0]), &h_yoffset[0], &errorx[0], &errorx[0], 0, &(r_sys[0]));

        if(!c)
        {
          mult[c][i][j]->SetMarkerColor(fMarkerColor[4]);
        }
        else
        {
          mult[c][i][j]->SetMarkerColor(fMarkerColor[0]);
        }

        mult[c][i][j]->SetMarkerSize(1);

        mult[c][i][j]->SetMarkerStyle(fMarkerStyle[0][c]);
        mult[c][i][j]->GetYaxis()->SetTitle("");

        mult[c][i][j]->GetXaxis()->SetTitle("");

        mult[c][i][j]->SetTitle("");

        sys[c][i][j]->SetFillColor(fMarkerColor[j]);

        if(!r_empty)
        {
          if(c) c1.cd(i+1);
          else c2.cd(i+1);
          if(mult[c][i][j])
          {
            if((!c && !axisflagh1) || (c && !axisflagh2))
            {
              mult[c][i][j]->Draw("SAMEPA");
              mult[c][i][j]->GetXaxis()->SetLimits(0.1,0.9);
              mult[c][i][j]->SetMinimum(-0.4);
              mult[c][i][j]->SetMaximum(5.);
              mult[c][i][j]->GetXaxis()->SetLabelSize(0.06);
              mult[c][i][j]->GetYaxis()->SetLabelSize(0.06);
              mult[c][i][j]->SetTitle("");
              if(i>4) gPad->SetBottomMargin(.15);
              if(i==0 || i==5) gPad->SetLeftMargin(.22);
              if(i==8)
              {
                mult[c][i][j]->GetXaxis()->SetTitle("#font[ 12]{z}");
                mult[c][i][j]->GetXaxis()->SetTitleSize(0.08);
                mult[c][i][j]->GetXaxis()->SetTitleOffset(.8);
              }
              mult[c][i][j]->GetXaxis()->SetNdivisions(304,kTRUE);
              mult[c][i][j]->GetYaxis()->SetNdivisions(304,kTRUE);
              if(i==0)
              {
                if(c) mult[c][i][j]->GetYaxis()->SetTitle("#font[12]{M}^{#font[ 12]{h^{+}}}+ #font[ 12]{#delta}");
                else mult[c][i][j]->GetYaxis()->SetTitle("#font[12]{M}^{#font[ 12]{h^{-}}}+ #font[ 12]{#delta}");
                mult[c][i][j]->GetYaxis()->SetTitleSize(0.08);
              }
              lsys.Draw();
              if(j==3) sys[c][i][j]->Draw("SAME3");
              if(!c) axisflagh1=1;
              else axisflagh2=1;
              if(c) c1.Range(0.1,-0.4,0.9,5.);
              else c2.Range(0.1,-0.4,0.9,5.);
            }
            else
            {
              mult[c][i][j]->Draw("SAMEP");
              if(j==3) sys[c][i][j]->Draw("SAME3");
              mult[c][i][j]->GetXaxis()->SetLimits(0.1,0.9);
              mult[c][i][j]->SetMinimum(-0.4);
              mult[c][i][j]->SetMaximum(5.);
            }
          }
          if(c) c1.Update();
          else c2.Update();
        }
      }
    }
  }

  return 0;
}

/*
LOG
16/03/16 : apparently working version for [3,40] range
22/03/16 : changing storage form Lyon to CASTOR
0X/04/16 : minor fixes
20/04/16 : added parallelized version
22/04/16 : lifted version
*/

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
#include <TGraphErrors.h>
#include <TArrow.h>

#include "acceptance_collect.h"

// Flags
#define Y2006 0
#define Y2012 0
#define Y2016 1

// Outputs
#define dirroot "/sps/compass/npierre/Multiplicities/acceptance"

using namespace std;

int main()
{

  int year=0;

  if(Y2006) year=2006;
  else if(Y2012) year=2012;
  else if(Y2016) year=2016;

  // Files

  for(int filen=0; filen<1; filen++)
  {

    double dummyd;

    ifstream DIS_file(Form("acceptance/%d/DIS_%d.txt",year,filen));
    ifstream had_file(Form("acceptance/%d/hadron_%d.txt",year,filen));

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
              DIS_file >> dummyd;
              fNDIS_evt[0][i][j][k] += dummyd;
              DIS_file >> dummyd;
              fNDIS_evt_c[0][i][j][k] += dummyd;
              DIS_file >> dummyd;
              fNDIS_evt_MC[0][i][j][k] += dummyd;
              DIS_file >> dummyd;
              fNDIS_evt[1][i][j][k] += dummyd;
              DIS_file >> dummyd;
              fNDIS_evt_c[1][i][j][k] += dummyd;
              DIS_file >> dummyd;
              fNDIS_evt_MC[1][i][j][k] += dummyd;
              DIS_file >> dummyd;
              fNDIS_evt[2][i][j][k] += dummyd;
              DIS_file >> dummyd;
              fNDIS_evt_c[2][i][j][k] += dummyd;
              DIS_file >> dummyd;
              fNDIS_evt_MC[2][i][j][k] += dummyd;
            }

            had_file >> dummyd;
            fRcstr[i][j][k].tab[c][0][0] += dummyd;
            had_file >> dummyd;
            fRcstr_c[i][j][k].tab[c][0][0] += dummyd;
            had_file >> dummyd;
            fGnrt[i][j][k].tab[c][0][0] += dummyd;

            had_file >> dummyd;
            fRcstr[i][j][k].tab[c][0][1] += dummyd;
            had_file >> dummyd;
            fRcstr_c[i][j][k].tab[c][0][1] += dummyd;
            had_file >> dummyd;
            fGnrt[i][j][k].tab[c][0][1] += dummyd;

            had_file >> dummyd;
            fRcstr[i][j][k].tab[c][0][2] += dummyd;
            had_file >> dummyd;
            fRcstr_c[i][j][k].tab[c][0][2] += dummyd;
            had_file >> dummyd;
            fGnrt[i][j][k].tab[c][0][2] += dummyd;

            had_file >> dummyd;
            fRcstr[i][j][k].tab[c][0][3] += dummyd;
            had_file >> dummyd;
            fRcstr_c[i][j][k].tab[c][0][3] += dummyd;
            had_file >> dummyd;
            fGnrt[i][j][k].tab[c][0][3] += dummyd;
          }
        }
      }
    }

    DIS_file.close();
    had_file.close();
  }

  //TArrow t1(0.1,0.9,0.1,0.9,0.1,"|>");
  //TArrow t2(0.1,0.1,0.1,0.9,0.1,"<|");

  TCanvas c5("Hadron_Acceptance","Hadron_Acceptance",3200,1600);
  //t1.Draw(); t2.Draw();

  TCanvas c6("Pion_Acceptance","Pion_Acceptance",3200,1600);
  //t1.Draw(); t2.Draw();

  TCanvas c7("Kaon_Acceptance","Kaon_Acceptance",3200,1600);
  //t1.Draw(); t2.Draw();

  c5.SetFillColor(0);
  //c5.SetFrameFillStyle(0);
  c6.SetFillColor(0);
  //c6.SetFrameFillStyle(0);
  c7.SetFillColor(0);
  //c7.SetFrameFillStyle(0);

  c5.Divide(9,5,0,0);
  c6.Divide(9,5,0,0);
  c7.Divide(9,5,0,0);

  TGraphErrors* H_acc[2][9][5];
  TGraphErrors* P_acc[2][9][5];
  TGraphErrors* K_acc[2][9][5];

  double z_range[12] = {.225,.275,.325,.375,.425,.475,.525,.575,.625,.675,.725,.8};

  ofstream ofs(Form("%s/%d/acceptance.txt",dirroot,year), std::ofstream::out | std::ofstream::trunc);
  ofstream lepto(Form("%s/%d/lepto.txt",dirroot,year), std::ofstream::out | std::ofstream::trunc);

  for(int c=0; c<2; c++)
  {
    for(int i=0; i<9; i++)
    {
      for(int j=0; j<5; j++)
      {

        std::vector<double> p_a;
        std::vector<double> k_a;
        std::vector<double> h_a;
        std::vector<double> p_err;
        std::vector<double> k_err;
        std::vector<double> h_err;

        std::vector<double> z_range_p;
        std::vector<double> z_range_k;
        std::vector<double> z_range_h;

        for(int l=0; l<12; l++)
        {
          z_range_p.push_back(z_range[l]);
          z_range_k.push_back(z_range[l]);
          z_range_h.push_back(z_range[l]);
        }

        for(int k=0; k<12; k++)
        {
          /*cout << "data count : " << fNDIS_evt_MC[0][i][j][k] << endl;
          cout << "MC count : " << fNDIS_evt_MC[0][i][j][k] << endl;
          cout << "Rcst count : " << fRcstr[i][j][k].tab[c][0][3] << endl;
          cout << "Gnrt count : " << fGnrt[i][j][k].tab[c][0][3] << endl;*/

          //cout << fRcstr[i][j][k].tab[c][0][2] << " " << fGnrt[i][j][k].tab[c][0][2] << endl;

          fAcceptance[i][j][k].tab[c][0][0] = ((fNDIS_evt[0][i][j][k] && fNDIS_evt_MC[0][i][j][k] && fGnrt[i][j][k].tab[c][0][0]) ? double((fRcstr[i][j][k].tab[c][0][0]/fNDIS_evt[0][i][j][k])/(fGnrt[i][j][k].tab[c][0][0]/fNDIS_evt_MC[0][i][j][k])) : 0);
          fAcceptance[i][j][k].tab[c][0][1] = ((fNDIS_evt[1][i][j][k] && fNDIS_evt_MC[1][i][j][k] && fGnrt[i][j][k].tab[c][0][1]) ? double((fRcstr[i][j][k].tab[c][0][1]/fNDIS_evt[1][i][j][k])/(fGnrt[i][j][k].tab[c][0][1]/fNDIS_evt_MC[1][i][j][k])) : 0);
          fAcceptance[i][j][k].tab[c][0][2] = ((fNDIS_evt[2][i][j][k] && fNDIS_evt_MC[2][i][j][k] && fGnrt[i][j][k].tab[c][0][2]) ? double((fRcstr[i][j][k].tab[c][0][2]/fNDIS_evt[2][i][j][k])/(fGnrt[i][j][k].tab[c][0][2]/fNDIS_evt_MC[2][i][j][k])) : 0);
          fAcceptance[i][j][k].tab[c][0][3] = ((fNDIS_evt[0][i][j][k] && fNDIS_evt_MC[0][i][j][k] && fGnrt[i][j][k].tab[c][0][3]) ? double((fRcstr[i][j][k].tab[c][0][3]/fNDIS_evt[0][i][j][k])/(fGnrt[i][j][k].tab[c][0][3]/fNDIS_evt_MC[0][i][j][k])) : 0);

          if(fAcceptance[i][j][k].tab[c][0][0]<0) fAcceptance[i][j][k].tab[c][0][0]=0;
          if(fAcceptance[i][j][k].tab[c][0][1]<0) fAcceptance[i][j][k].tab[c][0][1]=0;
          if(fAcceptance[i][j][k].tab[c][0][2]<0) fAcceptance[i][j][k].tab[c][0][2]=0;
          if(fAcceptance[i][j][k].tab[c][0][3]<0) fAcceptance[i][j][k].tab[c][0][3]=0;

          //cout << fAcceptance[i][j][k].tab[c][0][0] << " " << fAcceptance[i][j][k].tab[c][0][1] << " " << fAcceptance[i][j][k].tab[c][0][2] << " " << fAcceptance[i][j][k].tab[c][0][3] << endl;

          double rd_prime[3];
          double rd[3];
          double gd[3];
          double rh_prime[4];
          double rh[4];
          double gh[4];

          rd_prime[0] = fNDIS_evt[0][i][j][k] - fNDIS_evt_c[0][i][j][k];
          rd_prime[1] = fNDIS_evt[1][i][j][k] - fNDIS_evt_c[1][i][j][k];
          rd_prime[2] = fNDIS_evt[2][i][j][k] - fNDIS_evt_c[2][i][j][k];

          rd[0] = fNDIS_evt_c[0][i][j][k];
          rd[1] = fNDIS_evt_c[1][i][j][k];
          rd[2] = fNDIS_evt_c[2][i][j][k];

          gd[0] = fNDIS_evt_MC[0][i][j][k];
          gd[1] = fNDIS_evt_MC[1][i][j][k];
          gd[2] = fNDIS_evt_MC[2][i][j][k];

          rh_prime[0] = fRcstr[i][j][k].tab[c][0][0] - fRcstr_c[i][j][k].tab[c][0][0];
          rh_prime[1] = fRcstr[i][j][k].tab[c][0][1] - fRcstr_c[i][j][k].tab[c][0][1];
          rh_prime[2] = fRcstr[i][j][k].tab[c][0][2] - fRcstr_c[i][j][k].tab[c][0][2];
          rh_prime[3] = fRcstr[i][j][k].tab[c][0][3] - fRcstr_c[i][j][k].tab[c][0][3];

          rh[0] = fRcstr_c[i][j][k].tab[c][0][0];
          rh[1] = fRcstr_c[i][j][k].tab[c][0][1];
          rh[2] = fRcstr_c[i][j][k].tab[c][0][2];
          rh[3] = fRcstr_c[i][j][k].tab[c][0][3];

          gh[0] = fGnrt[i][j][k].tab[c][0][0];
          gh[1] = fGnrt[i][j][k].tab[c][0][1];
          gh[2] = fGnrt[i][j][k].tab[c][0][2];
          gh[3] = fGnrt[i][j][k].tab[c][0][3];

          fAcceptance[i][j][k].tab[c][1][0] = sqrt(((rh[0]>0 && gh[0]>0 && gd[0]>0 && ((rd[0]+rd_prime[0])!=0)) ? double(pow(gd[0]/(rd_prime[0]+rd[0]),2)
                                                                              *((rh[0]+1)*(gh[0]-rh[0]+1)/(pow(gh[0]+2,2)*(gh[0]+3))+rh_prime[0]/pow(gh[0],2)+pow(rh_prime[0],2)/pow(gh[0],3))
                                                                              +pow(gd[0]/(rd[0]+rd_prime[0]),4)*pow((rh[0]+rh_prime[0])/gh[0],2)*((rd[0]+1)*(gd[0]-rd[0]+1)/(pow(gd[0]+2,2)*(gd[0]+3))+rd_prime[0]/pow(gd[0],2)+pow(rd_prime[0],2)/pow(gd[0],3))) : 0));
          fAcceptance[i][j][k].tab[c][1][1] = sqrt(((rh[1]>0 && gh[1]>0 && gd[1]>0 && ((rd[1]+rd_prime[1])!=0)) ? double(pow(gd[1]/(rd_prime[1]+rd[1]),2)
                                                                              *((rh[1]+1)*(gh[1]-rh[1]+1)/(pow(gh[1]+2,2)*(gh[1]+3))+rh_prime[1]/pow(gh[1],2)+pow(rh_prime[1],2)/pow(gh[1],3))
                                                                              +pow(gd[1]/(rd[1]+rd_prime[1]),4)*pow((rh[1]+rh_prime[1])/gh[1],2)*((rd[1]+1)*(gd[1]-rd[1]+1)/(pow(gd[1]+2,2)*(gd[1]+3))+rd_prime[1]/pow(gd[1],2)+pow(rd_prime[1],2)/pow(gd[1],3))) : 0));
          fAcceptance[i][j][k].tab[c][1][2] = sqrt(((rh[2]>0 && gh[2]>0 && gd[2]>0 && ((rd[2]+rd_prime[2])!=0)) ? double(pow(gd[2]/(rd_prime[2]+rd[2]),2)
                                                                              *((rh[2]+1)*(gh[2]-rh[2]+1)/(pow(gh[2]+2,2)*(gh[2]+3))+rh_prime[2]/pow(gh[2],2)+pow(rh_prime[2],2)/pow(gh[2],3))
                                                                              +pow(gd[2]/(rd[2]+rd_prime[2]),4)*pow((rh[2]+rh_prime[2])/gh[2],2)*((rd[2]+1)*(gd[2]-rd[2]+1)/(pow(gd[2]+2,2)*(gd[2]+3))+rd_prime[2]/pow(gd[2],2)+pow(rd_prime[2],2)/pow(gd[2],3))) : 0));
          fAcceptance[i][j][k].tab[c][1][3] = sqrt(((rh[3]>0 && gh[3]>0 && gd[0]>0 && ((rd[0]+rd_prime[0])!=0)) ? double(pow(gd[0]/(rd_prime[0]+rd[0]),2)
                                                                              *((rh[3]+1)*(gh[3]-rh[3]+1)/(pow(gh[3]+2,2)*(gh[3]+3))+rh_prime[3]/pow(gh[3],2)+pow(rh_prime[3],2)/pow(gh[3],3))
                                                                              +pow(gd[0]/(rd[0]+rd_prime[0]),4)*pow((rh[3]+rh_prime[3])/gh[3],2)*((rd[0]+1)*(gd[0]-rd[0]+1)/(pow(gd[0]+2,2)*(gd[0]+3))+rd_prime[0]/pow(gd[0],2)+pow(rd_prime[0],2)/pow(gd[0],3))) : 0));

          if(fAcceptance[i][j][k].tab[c][0][0]==0) fAcceptance[i][j][k].tab[c][1][0]=0;
          if(fAcceptance[i][j][k].tab[c][0][1]==0) fAcceptance[i][j][k].tab[c][1][1]=0;
          if(fAcceptance[i][j][k].tab[c][0][2]==0) fAcceptance[i][j][k].tab[c][1][2]=0;
          if(fAcceptance[i][j][k].tab[c][0][3]==0) fAcceptance[i][j][k].tab[c][1][3]=0;

          /*if((fAcceptance[i][j][k].tab[c][0][3] != 0) && (fAcceptance[i][j][k].tab[c][1][3] > fAcceptance[i][j][k].tab[c][0][3]))
          {
            fAcceptance[i][j][k].tab[c][0][0] = 0;
            fAcceptance[i][j][k].tab[c][0][1] = 0;
            fAcceptance[i][j][k].tab[c][0][2] = 0;
            fAcceptance[i][j][k].tab[c][0][3] = 0;
            fAcceptance[i][j][k].tab[c][1][0] = 0;
            fAcceptance[i][j][k].tab[c][1][1] = 0;
            fAcceptance[i][j][k].tab[c][1][2] = 0;
            fAcceptance[i][j][k].tab[c][1][3] = 0;
          }

          if((fAcceptance[i][j][k].tab[c][0][0] != 0)
          && (fAcceptance[i][j][k].tab[c][1][0] != 0)
          && (fAcceptance[i][j][k].tab[c][1][0]/fAcceptance[i][j][k].tab[c][0][0] > 0.20))
          {
            fAcceptance[i][j][k].tab[c][0][0] = 0;
            fAcceptance[i][j][k].tab[c][0][1] = 0;
            fAcceptance[i][j][k].tab[c][0][2] = 0;
            fAcceptance[i][j][k].tab[c][0][3] = 0;
            fAcceptance[i][j][k].tab[c][1][0] = 0;
            fAcceptance[i][j][k].tab[c][1][1] = 0;
            fAcceptance[i][j][k].tab[c][1][2] = 0;
            fAcceptance[i][j][k].tab[c][1][3] = 0;
          }*/

          if((i==7 && j==4) || (i==8 && j==0) || (i==8 && j==4))
          {
            fAcceptance[i][j][k].tab[c][0][0] = 0;
            fAcceptance[i][j][k].tab[c][0][1] = 0;
            fAcceptance[i][j][k].tab[c][0][2] = 0;
            fAcceptance[i][j][k].tab[c][0][3] = 0;
            fAcceptance[i][j][k].tab[c][1][0] = 0;
            fAcceptance[i][j][k].tab[c][1][1] = 0;
            fAcceptance[i][j][k].tab[c][1][2] = 0;
            fAcceptance[i][j][k].tab[c][1][3] = 0;
          }

          p_a.push_back(fAcceptance[i][j][k].tab[c][0][0]);
          k_a.push_back(fAcceptance[i][j][k].tab[c][0][1]);
          h_a.push_back(fAcceptance[i][j][k].tab[c][0][3]);

          p_err.push_back(fAcceptance[i][j][k].tab[c][1][0]);
          k_err.push_back(fAcceptance[i][j][k].tab[c][1][1]);
          h_err.push_back(fAcceptance[i][j][k].tab[c][1][3]);

          //Output file
          //q_bin x_bin y_bin z_bin acc_pi acc_error_pi acc_k acc_error_k acc_p acc_error_p acc_h acc_error_h

          ofs << c << " " << fXrange[i] << " " << fYrange[j] << " " << fZrange[k] << " " <<
          fAcceptance[i][j][k].tab[c][0][0] << " " <<
          fAcceptance[i][j][k].tab[c][1][0] << " " <<
          fAcceptance[i][j][k].tab[c][0][1] << " " <<
          fAcceptance[i][j][k].tab[c][1][1] << " " <<
          fAcceptance[i][j][k].tab[c][0][2] << " " <<
          fAcceptance[i][j][k].tab[c][1][2] << " " <<
          fAcceptance[i][j][k].tab[c][0][3] << " " <<
          fAcceptance[i][j][k].tab[c][1][3] <<   endl;

          lepto << c << " " << fXrange[i] << " " << fYrange[j] << " " << fZrange[k] << " " <<
          (fAcceptance[i][j][k].tab[c][0][0] ? Double_t((fGnrt[i][j][k].tab[c][0][0]/pow(fNDIS_evt_MC[0][i][j][k],2)-pow(fGnrt[i][j][k].tab[c][0][0],2)/pow(fNDIS_evt_MC[0][i][j][k],3))*pow(fZrange[k],2))  : 0) << " " <<
          (fAcceptance[i][j][k].tab[c][0][1] ? Double_t((fGnrt[i][j][k].tab[c][0][1]/pow(fNDIS_evt_MC[0][i][j][k],2)-pow(fGnrt[i][j][k].tab[c][0][1],2)/pow(fNDIS_evt_MC[0][i][j][k],3))*pow(fZrange[k],2))  : 0) << " " <<
          (fAcceptance[i][j][k].tab[c][0][2] ? Double_t((fGnrt[i][j][k].tab[c][0][2]/pow(fNDIS_evt_MC[0][i][j][k],2)-pow(fGnrt[i][j][k].tab[c][0][2],2)/pow(fNDIS_evt_MC[0][i][j][k],3))*pow(fZrange[k],2))  : 0) << " " <<
          (fAcceptance[i][j][k].tab[c][0][3] ? Double_t((fGnrt[i][j][k].tab[c][0][3]/pow(fNDIS_evt_MC[0][i][j][k],2)-pow(fGnrt[i][j][k].tab[c][0][3],2)/pow(fNDIS_evt_MC[0][i][j][k],3))*pow(fZrange[k],2))  : 0) << " " << endl;

        }

        for(int k=12; k>0; k--)
        {
          if(!p_a[k-1]) {p_a.erase(p_a.begin()+k-1); p_err.erase(p_err.begin()+k-1); z_range_p.erase(z_range_p.begin()+k-1);}
          if(!k_a[k-1]) {k_a.erase(k_a.begin()+k-1); k_err.erase(k_err.begin()+k-1); z_range_k.erase(z_range_k.begin()+k-1);}
          if(!h_a[k-1]) {h_a.erase(h_a.begin()+k-1); h_err.erase(h_err.begin()+k-1); z_range_h.erase(z_range_h.begin()+k-1);}
        }

        bool p_a_empty = 0;
        bool k_a_empty = 0;
        bool h_a_empty = 0;

        if(!(int(p_a.size()))) p_a_empty = 1;
        if(!(int(k_a.size()))) k_a_empty = 1;
        if(!(int(h_a.size()))) h_a_empty = 1;

        H_acc[c][i][j] = new TGraphErrors(int(h_a.size()),&(z_range_h[0]),&(h_a[0]),0,&(h_err[0]));
        P_acc[c][i][j] = new TGraphErrors(int(p_a.size()),&(z_range_p[0]),&(p_a[0]),0,&(p_err[0]));
        K_acc[c][i][j] = new TGraphErrors(int(k_a.size()),&(z_range_k[0]),&(k_a[0]),0,&(k_err[0]));

        H_acc[c][i][j]->SetMarkerColor(fMarkerColor[c]);
        P_acc[c][i][j]->SetMarkerColor(fMarkerColor[c]);
        K_acc[c][i][j]->SetMarkerColor(fMarkerColor[c]);

        H_acc[c][i][j]->SetMarkerSize(2);
        P_acc[c][i][j]->SetMarkerSize(2);
        K_acc[c][i][j]->SetMarkerSize(2);

        H_acc[c][i][j]->SetMarkerStyle(fMarkerStyle[c]);
        P_acc[c][i][j]->SetMarkerStyle(fMarkerStyle[c]);
        K_acc[c][i][j]->SetMarkerStyle(fMarkerStyle[c]);

        H_acc[c][i][j]->GetYaxis()->SetTitle("");
        P_acc[c][i][j]->GetYaxis()->SetTitle("");
        K_acc[c][i][j]->GetYaxis()->SetTitle("");

        H_acc[c][i][j]->GetXaxis()->SetTitle("");
        P_acc[c][i][j]->GetXaxis()->SetTitle("");
        K_acc[c][i][j]->GetXaxis()->SetTitle("");

        H_acc[c][i][j]->SetTitle("");
        P_acc[c][i][j]->SetTitle("");
        K_acc[c][i][j]->SetTitle("");

        if(!h_a_empty)
        {
          c5.cd(i+1+j*9);
          if(H_acc[c][i][j])
          {
            if(!c)
            {
              H_acc[c][i][j]->Draw("SAMEPA");
              H_acc[c][i][j]->GetXaxis()->SetLimits(0.1,0.9);
              H_acc[c][i][j]->SetMinimum(0.2);
              H_acc[c][i][j]->SetMaximum(1.2);
            }
            else
            {
              H_acc[c][i][j]->Draw("SAMEP");
              H_acc[c][i][j]->GetXaxis()->SetLimits(0.1,0.9);
              H_acc[c][i][j]->SetMinimum(0.2);
              H_acc[c][i][j]->SetMaximum(1.2);
            }
          }
          c5.Update();
        }

        if(!p_a_empty)
        {
          c6.cd(i+1+j*9);
          if(P_acc[c][i][j])
          {
            if(!c)
            {
              P_acc[c][i][j]->Draw("SAMEPA");
              P_acc[c][i][j]->GetXaxis()->SetLimits(0.1,0.9);
              P_acc[c][i][j]->SetMinimum(0.3);
              P_acc[c][i][j]->SetMaximum(1.1);
            }
            else
            {
              P_acc[c][i][j]->Draw("SAMEP");
              P_acc[c][i][j]->GetXaxis()->SetLimits(0.1,0.9);
              P_acc[c][i][j]->SetMinimum(0.3);
              P_acc[c][i][j]->SetMaximum(1.1);
            }
          }
          c6.Update();
        }

        if(!k_a_empty)
        {
          c7.cd(i+1+j*9);
          if(K_acc[c][i][j])
          {
            if(!c)
            {
              K_acc[c][i][j]->Draw("SAMEPA");
              K_acc[c][i][j]->GetXaxis()->SetLimits(0.1,0.9);
              K_acc[c][i][j]->SetMinimum(0.2);
              K_acc[c][i][j]->SetMaximum(1.2);
            }
            else
            {
              K_acc[c][i][j]->Draw("SAMEP");
              K_acc[c][i][j]->GetXaxis()->SetLimits(0.1,0.9);
              K_acc[c][i][j]->SetMinimum(0.2);
              K_acc[c][i][j]->SetMaximum(1.2);
            }
          }
          c7.Update();
        }

      }
    }
  }

  c5.Update();
  c6.Update();
  c7.Update();

  c5.Print(Form("%s/%d/hadron_acceptance.pdf",dirroot,year));
  c6.Print(Form("%s/%d/pion_acceptance.pdf",dirroot,year));
  c7.Print(Form("%s/%d/kaon_acceptance.pdf",dirroot,year));

  ofs.close();

  return 0;
}

#include "compMCRD.h"

//Inputs
#define target_file_2012 "data/target-107924-109081.dat"
#define target_file_2016 "data/target-274508-274901.dat"

// Flags
#define Y2006 0
#define Y2012 0
#define Y2016 1
#define RCUTSTUDY_ON 0

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
    float z, x, y, r, dummy;
    fin >> z >> dummy >> dummy >> dummy >> dummy >> r >> dummy >> x >> y;
    fZv.push_back(z);
    fXv.push_back(x);
    fYv.push_back(y);
    fRv.push_back(r);
  }
  cout<<"INFO : Target cell description loaded"<<endl;
}

void CellCenter(Double_t z, Double_t& xc, Double_t& yc, Double_t& R)
{
  xc = 1000000;
  yc = 1000000;

  for(Int_t i = 0; i < int(fZv.size()-1); i++)
  {
    Double_t z1 = fZv[i];
    Double_t z2 = fZv[i+1];

    if( z2 < z ) continue;
    if( z1 > z ) continue;

    Double_t xc1 = fXv[i];
    Double_t xc2 = fXv[i+1];

    Double_t yc1 = fYv[i];
    Double_t yc2 = fYv[i+1];

    Double_t rc1 = fRv[i];
    Double_t rc2 = fRv[i+1];

    Double_t dxcdz = (xc2-xc1)/(z2-z1);
    Double_t dycdz = (yc2-yc1)/(z2-z1);
    Double_t drcdz = (rc2-rc1)/(z2-z1);

    Double_t dz = z-z1;
    xc = xc1 + dxcdz*dz;
    yc = yc1 + dycdz*dz;
    R = rc1 + drcdz*dz;

    break;
  }
}

bool InTarget(Double_t xvtx, Double_t yvtx, Double_t zvtx)
{
  Double_t xc, yc, R;
  CellCenter(zvtx, xc, yc, R);
  Double_t dx = xvtx-xc;
  Double_t dy = yvtx-yc;
  Double_t r = sqrt(dx*dx + dy*dy);

  return( r < 1.9 && yvtx < 1.2 );
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

void readKinCuts(string pFile)
{
  string dummy;
  ifstream list(pFile);
  list >> dummy >> fXmin;
  list >> dummy >> fXmax;
  list >> dummy >> fYmin;
  list >> dummy >> fYmax;
  list >> dummy >> fWmin;
  list >> dummy >> fWmax;
  list >> dummy >> fPmin;
  list >> dummy >> fPmax;
  list.close();
}

void create_kin_plots()
{
  for(int i=0; i<5; i++)
  {
    fKinematicsRD[i][0] = new TH1F(Form("Q^2 Ratio %s",trigname[i].c_str()), Form("Q2 %s",trigname[i].c_str()), 50, -1, 2);
    fKinematicsRD[i][1] = new TH1F(Form("x_{Bj} Ratio %s",trigname[i].c_str()), Form("x_{Bj} %s",trigname[i].c_str()), 50, -3, 0);
    fKinematicsRD[i][2] = new TH1F(Form("y Ratio %s",trigname[i].c_str()), Form("y %s",trigname[i].c_str()), 50, 0, 1);
    fKinematicsRD[i][3] = new TH1F(Form("z Ratio %s",trigname[i].c_str()), Form("z %s",trigname[i].c_str()), 50, 0, 1);
    fKinematicsRD[i][4] = new TH1F(Form("W Ratio %s",trigname[i].c_str()), Form("W %s",trigname[i].c_str()), 50, 2, 18);
    fKinematicsRD[i][5] = new TH1F(Form("#nu Ratio %s",trigname[i].c_str()), Form("nu %s",trigname[i].c_str()), 50, 0, 160);
    fKinematicsRD[i][6] = new TH1F(Form("E_{#mu} Ratio %s",trigname[i].c_str()), Form("E_{mu} %s",trigname[i].c_str()), 50, 140, 180);
    fKinematicsRD[i][7] = new TH1F(Form("E_{#mu'} Ratio %s",trigname[i].c_str()), Form("E_{mu'} %s",trigname[i].c_str()), 50, 0, 160);
    fKinematicsRD[i][8] = new TH1F(Form("#theta Ratio %s",trigname[i].c_str()), Form("theta %s",trigname[i].c_str()), 50, 0, 0.05);
    fKinematicsRD[i][9] = new TH1F(Form("#Phi Ratio %s",trigname[i].c_str()), Form("phi %s",trigname[i].c_str()), 50, -1.7, 1.7);
    fKinematicsRD[i][10] = new TH1F(Form("Vertex Ratio %s",trigname[i].c_str()), Form("Vertex %s",trigname[i].c_str()), 50, -340, -60);
    fKinematicsRD[i][12] = new TH1F(Form("p_{hadron+e} Ratio %s",trigname[i].c_str()), Form("p_{hadron+e} %s",trigname[i].c_str()), 50, 0, 40);
    fKinematicsRD[i][13] = new TH1F(Form("#theta_{hadron+e} Ratio %s",trigname[i].c_str()), Form("theta_{hadron+e} %s",trigname[i].c_str()), 50, 0, 0.25);
    fKinematicsRD[i][14] = new TH1F(Form("#Phi_{hadron+e,lab} Ratio %s",trigname[i].c_str()), Form("phi_{hadron+e,lab} %s",trigname[i].c_str()), 50, -3.5, 3.5);
    fKinematicsRD[i][15] = new TH1F(Form("#Phi_{hadron+e,prod.pl} Ratio %s",trigname[i].c_str()), Form("phi_{hadron+e,prod.pl} %s",trigname[i].c_str()), 50, 0, 3.5);
    fKinematicsRD[i][16] = new TH1F(Form("p_{T} Ratio %s",trigname[i].c_str()), Form("p_{T} %s",trigname[i].c_str()), 50, 0, 3);
    fKinematicsMC[i][0] = new TH1F(Form("Q^2 %s",trigname[i].c_str()), Form("Q2 Ratio %s",trigname[i].c_str()), 50, -1, 2);
    fKinematicsMC[i][1] = new TH1F(Form("x_{Bj} %s",trigname[i].c_str()), Form("x_{Bj} Ratio %s",trigname[i].c_str()), 50, -3, 0);
    fKinematicsMC[i][2] = new TH1F(Form("y %s",trigname[i].c_str()), Form("y Ratio %s",trigname[i].c_str()), 50, 0, 1);
    fKinematicsMC[i][3] = new TH1F(Form("z %s",trigname[i].c_str()), Form("z Ratio %s",trigname[i].c_str()), 50, 0, 1);
    fKinematicsMC[i][4] = new TH1F(Form("W %s",trigname[i].c_str()), Form("W Ratio %s",trigname[i].c_str()), 50, 2, 18);
    fKinematicsMC[i][5] = new TH1F(Form("#nu %s",trigname[i].c_str()), Form("nu Ratio %s",trigname[i].c_str()), 50, 0, 160);
    fKinematicsMC[i][6] = new TH1F(Form("E_{#mu} %s",trigname[i].c_str()), Form("E_{mu} Ratio %s",trigname[i].c_str()), 50, 140, 180);
    fKinematicsMC[i][7] = new TH1F(Form("E_{#mu'} %s",trigname[i].c_str()), Form("E_{mu'} Ratio %s",trigname[i].c_str()), 50, 0, 160);
    fKinematicsMC[i][8] = new TH1F(Form("#theta %s",trigname[i].c_str()), Form("theta Ratio %s",trigname[i].c_str()), 50, 0, 0.05);
    fKinematicsMC[i][9] = new TH1F(Form("#Phi %s",trigname[i].c_str()), Form("phi Ratio %s",trigname[i].c_str()), 50, -1.7, 1.7);
    fKinematicsMC[i][10] = new TH1F(Form("Vertex %s",trigname[i].c_str()), Form("Vertex Ratio %s",trigname[i].c_str()), 50, -340, -60);
    fKinematicsMC[i][12] = new TH1F(Form("p_{hadron+e} %s",trigname[i].c_str()), Form("p_{hadron+e} Ratio %s",trigname[i].c_str()), 50, 0, 40);
    fKinematicsMC[i][13] = new TH1F(Form("#theta_{hadron+e} %s",trigname[i].c_str()), Form("theta_{hadron+e} Ratio %s",trigname[i].c_str()), 50, 0, 0.25);
    fKinematicsMC[i][14] = new TH1F(Form("#Phi_{hadron+e,lab} %s",trigname[i].c_str()), Form("phi_{hadron+e,lab} Ratio %s",trigname[i].c_str()), 50, -3.5, 3.5);
    fKinematicsMC[i][15] = new TH1F(Form("#Phi_{hadron+e,prod.pl} %s",trigname[i].c_str()), Form("phi_{hadron+e,prod.pl} Ratio %s",trigname[i].c_str()), 50, 0, 3.5);
    fKinematicsMC[i][16] = new TH1F(Form("p_{T} %s",trigname[i].c_str()), Form("p_{T} Ratio %s",trigname[i].c_str()), 50, 0, 3);
    fKinematicsRD[i][17] = new TH1F(Form("Vertex (Hadron+e) Ratio %s",trigname[i].c_str()), Form("Vertex (Hadron+e) %s",trigname[i].c_str()), 50, -340, -60);
    fKinematicsRD[i][18] = new TH1F(Form("Vertex (e) Ratio %s",trigname[i].c_str()), Form("Vertex (e) %s",trigname[i].c_str()), 50, -340, -60);
    fKinematicsRD[i][19] = new TH1F(Form("Vertex (Hadron) Ratio %s",trigname[i].c_str()), Form("Vertex (Hadron) %s",trigname[i].c_str()), 50, -340, -60);
    fKinematicsMC[i][17] = new TH1F(Form("Vertex (Hadron+e) %s",trigname[i].c_str()), Form("Vertex (Hadron+e) Ratio %s",trigname[i].c_str()), 50, -340, -60);
    fKinematicsMC[i][18] = new TH1F(Form("Vertex (e) %s",trigname[i].c_str()), Form("Vertex (e) Ratio %s",trigname[i].c_str()), 50, -340, -60);
    fKinematicsMC[i][19] = new TH1F(Form("Vertex (Hadron) %s",trigname[i].c_str()), Form("Vertex (Hadron) Ratio %s",trigname[i].c_str()), 50, -340, -60);
    BinLogX(fKinematicsRD[i][0]);
    BinLogX(fKinematicsMC[i][0]);
    BinLogX(fKinematicsRD[i][1]);
    BinLogX(fKinematicsMC[i][1]);
  }

  fKinematicsRD[0][11] = new TH1F("phi_{e,prod.pl}","phi_{e,prod.pl}", 50, 0, 3.5);
  fKinematicsMC[0][11] = new TH1F("phi_{e,prod.pl} Ratio","phi_{e,prod.pl} Ratio", 50, 0, 3.5);
  fECAL0RD = new TH2F("ECAL0 Map","ECAL0 Map", 1000, -80, 80, 1000, -80, 80);
  fECAL0MC = new TH2F("ECAL0 Map MC","ECAL0 Map MC", 1000, -80, 80, 1000, -80, 80);
  fVertexRD[0] = new TH1F("Vertex Endpoint 0","Vertex Endpoint 0", 500, -1000, 2000);
  fVertexMCb[0] = new TH1F("Vertex Endpoint 0 MC","EVertex Endpoint 0 MC", 500, -1000, 2000);
  fVertexRD[1] = new TH1F("Vertex Endpoint 1","Vertex Endpoint 1", 500, -1000, 2000);
  fVertexMCb[1] = new TH1F("Vertex Endpoint 1 MC","EVertex Endpoint 1 MC", 500, -1000, 2000);
  fVertexRD[2] = new TH1F("Vertex Endpoint 2","Vertex Endpoint 2", 500, -1000, 2000);
  fVertexMCb[2] = new TH1F("Vertex Endpoint 2 MC","EVertex Endpoint 2 MC", 500, -1000, 2000);
  fVertexRD[3] = new TH1F("Vertex Endpoint 3","Vertex Endpoint 3", 500, -1000, 2000);
  fVertexMCb[3] = new TH1F("Vertex Endpoint 3 MC","EVertex Endpoint 3 MC", 500, -1000, 2000);
  fVertexRD[4] = new TH1F("Vertex Endpoint 4","Vertex Endpoint 4", 500, -1000, 2000);
  fVertexMCb[4] = new TH1F("Vertex Endpoint 4 MC","EVertex Endpoint 4 MC", 500, -1000, 2000);
  fThetaRDp[0] = new TH2F("theta_y RD", "theta_y RD", 100, -0.005, 0.005, 100, 140, 180);
  fThetaRDp[1] = new TH2F("theta_x RD", "theta_x RD", 100, -0.005, 0.005, 100, 140, 180);
  fThetaRDp[2] = new TH2F("theta_xy RD", "theta_xy RD", 100, -0.005, 0.005, 100, -0.005, 0.005);
  fThetaMCp[0] = new TH2F("theta_y MC", "theta_y MC", 100, -0.005, 0.005, 100, 140, 180);
  fThetaMCp[1] = new TH2F("theta_x MC", "theta_x MC", 100, -0.005, 0.005, 100, 140, 180);
  fThetaMCp[2] = new TH2F("theta_xy MC", "theta_xy MC", 100, -0.005, 0.005, 100, -0.005, 0.005);
  fTarget = new TH1F("Target yz", "Target yz", 1000, -500, 500);
  fTargetMC = new TH1F("TargetMC yz", "TargetMC yz", 1000, -500, 500);
  fTarget2D = new TH2F("Target2D yz", "Target2D yz", 1000, -400, 0, 1000, -5, 5);
  fTarget2DMC = new TH2F("Target2DMC yz", "Target2DMC yz", 1000, -400, 0, 1000, -5, 5);

  for(int i=0; i<7; i++)
  {
    l1[0][i] = new TLine(0.1,0.4+i*0.2,100,0.4+i*0.2);
    l1[1][i] = new TLine(0.001,0.4+i*0.2,1,0.4+i*0.2);
    l1[2][i] = new TLine(0,0.4+i*0.2,1,0.4+i*0.2);
    l1[3][i] = new TLine(0,0.4+i*0.2,1,0.4+i*0.2);
    l1[4][i] = new TLine(2,0.4+i*0.2,18,0.4+i*0.2);
    l1[5][i] = new TLine(0,0.4+i*0.2,160,0.4+i*0.2);
    l1[6][i] = new TLine(140,0.4+i*0.2,180,0.4+i*0.2);
    l1[7][i] = new TLine(0,0.4+i*0.2,160,0.4+i*0.2);
    l1[8][i] = new TLine(0,0.4+i*0.2,0.05,0.4+i*0.2);
    l1[9][i] = new TLine(-1.7,0.4+i*0.2,1.7,0.4+i*0.2);
    l1[10][i] = new TLine(-340,0.4+i*0.2,-60,0.4+i*0.2);
    l1[11][i] = new TLine(0,0.4+i*0.2,3.5,0.4+i*0.2);
    l1[12][i] = new TLine(0,0.4+i*0.2,40,0.4+i*0.2);
    l1[13][i] = new TLine(0,0.4+i*0.2,0.25,0.4+i*0.2);
    l1[14][i] = new TLine(-3.5,0.4+i*0.2,3.5,0.4+i*0.2);
    l1[15][i] = new TLine(0,0.4+i*0.2,3.5,0.4+i*0.2);
    l1[16][i] = new TLine(0,0.4+i*0.2,3,0.4+i*0.2);
    l1[17][i] = new TLine(-340,0.4+i*0.2,-60,0.4+i*0.2);
    l1[18][i] = new TLine(-340,0.4+i*0.2,-60,0.4+i*0.2);
    l1[19][i] = new TLine(-340,0.4+i*0.2,-60,0.4+i*0.2);
    for(int j=0; j<20; j++)
    {
      l1[j][i]->SetLineStyle(fLineStyle[i]);
      l1[j][i]->SetLineWidth(1);
    }
  }
}

void plotting_ratio(int i, int j)
{
  // for(int tt=0; tt<fKinematicsRD[idx][0]->GetNbinsX(); tt++)
  // {
  //   fError.push_back((fKinematicsRD[idx][0]->GetBinError(tt) && fKinematicsMC2[idx][0]->GetBinError(tt) ? sqrt(pow(1/fKinematicsMC1[idx][0]->GetBinError(tt),2)+pow(1/fKinematicsMC2[idx][0]->GetBinError(tt),2)):0));
  // }
  fKinematicsRD[i][j]->Sumw2();
  fKinematicsMC[i][j]->Sumw2();
  fCountingMC[i][j] = fKinematicsMC[i][j]->GetEntries();
  fCountingRD[i][j] = fKinematicsRD[i][j]->GetEntries();
  fKinematicsMC[i][j]->Scale(1/fKinematicsMC[i][j]->GetEntries());
  fKinematicsRD[i][j]->Scale(1/fKinematicsRD[i][j]->GetEntries());
  fKinematicsRatio[i][j] = (TH1F*)fKinematicsRD[i][j]->Clone();
  fKinematicsRatio[i][j]->SetStats(0);
  fKinematicsRatio[i][j]->Divide(fKinematicsMC[i][j]);
  // for(int tt=0; tt<fKinematicsRatio[idx][0]->GetNbinsX(); tt++)
  // {
  //   fKinematicsRatio[idx][0]->SetBinError(tt,fError[tt]);
  // }
  // fError.clear();
  fKinematicsRatio[i][j]->SetMarkerStyle(21);
  fKinematicsRatio[i][j]->SetFillColor(kYellow-7);
  fKinematicsRatio[i][j]->SetMaximum(2.);
  fKinematicsRatio[i][j]->SetMinimum(0.);
  fKinematicsRatio[i][j]->Draw("PE2");
  fKinematicsRatio[i][j]->GetXaxis()->SetLabelSize(0.08);
  fKinematicsRatio[i][j]->GetYaxis()->SetLabelSize(0.08);
  fKinematicsRatio[i][j]->GetYaxis()->SetNdivisions(2,kTRUE);
  for(int tt=0; tt<7; tt++)
  {
    l1[j][tt]->Draw();
  }
}

void plotting_device(int i, int j)
{
  fKinematicsMC[i][j]->SetLineColor(kBlue);
  fKinematicsMC[i][j]->SetFillColor(kBlue);
  fKinematicsMC[i][j]->SetStats(0);
  fKinematicsMC[i][j]->SetMinimum(0.);
  fKinematicsMC[i][j]->Draw();
  fKinematicsMC[i][j]->GetXaxis()->SetLabelSize(0.08);
  fKinematicsMC[i][j]->GetYaxis()->SetLabelSize(0.08);
  fKinematicsMC[i][j]->GetYaxis()->SetNdivisions(304,kTRUE);
  fKinematicsRD[i][j]->SetLineColor(kRed);
  fKinematicsRD[i][j]->Draw("SAME");
}

void plotting_ratio_vertex(int i, int j)
{
  // for(int tt=0; tt<fKinematicsRD[idx][0]->GetNbinsX(); tt++)
  // {
  //   fError.push_back((fKinematicsRD[idx][0]->GetBinError(tt) && fKinematicsMC2[idx][0]->GetBinError(tt) ? sqrt(pow(1/fKinematicsMC1[idx][0]->GetBinError(tt),2)+pow(1/fKinematicsMC2[idx][0]->GetBinError(tt),2)):0));
  // }
  fKinematicsRD[i][j]->Sumw2();
  fKinematicsMC[i][j]->Sumw2();
  fCountingMC[i][j] = fKinematicsMC[i][j]->GetEntries();
  fCountingRD[i][j] = fKinematicsRD[i][j]->GetEntries();
  fKinematicsMC[i][j]->Scale(1/fKinematicsMC[i][j]->GetEntries());
  fKinematicsRD[i][j]->Scale(1/fKinematicsRD[i][j]->GetEntries());
  fKinematicsRatio[i][j] = (TH1F*)fKinematicsRD[i][j]->Clone();
  fKinematicsRatio[i][j]->SetStats(0);
  fKinematicsRatio[i][j]->Divide(fKinematicsMC[i][j]);
  // for(int tt=0; tt<fKinematicsRatio[idx][0]->GetNbinsX(); tt++)
  // {
  //   fKinematicsRatio[idx][0]->SetBinError(tt,fError[tt]);
  // }
  // fError.clear();
  fKinematicsRatio[i][j]->SetMarkerStyle(21);
  fKinematicsRatio[i][j]->SetFillColor(kYellow-7);
  fKinematicsRatio[i][j]->SetMaximum(2.);
  fKinematicsRatio[i][j]->SetMinimum(0.);
  fKinematicsRatio[i][j]->Draw("PE2");
  fKinematicsRatio[i][j]->GetXaxis()->SetLabelSize(0.08);
  fKinematicsRatio[i][j]->GetYaxis()->SetLabelSize(0.08);
  fKinematicsRatio[i][j]->GetYaxis()->SetNdivisions(2,kTRUE);
  for(int tt=0; tt<7; tt++)
  {
    l1[j][tt]->Draw();
  }
}

void save_kin_plots()
{
  TFile *output = new TFile("kinMCRD.root","NEW");

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
  c19.Divide(2,4);
  c20.Divide(2,4);
  c21.Divide(2,4);
  c22.Divide(1,2);
  c23.Divide(1,2);
  c24.Divide(1,2);
  c25.Divide(2,4);
  c26.Divide(1,2);
  c27.Divide(2,4);
  c28.Divide(1,2);
  c29.Divide(1,2);
  c30.Divide(1,2);
  c31.Divide(1,2);
  c32.Divide(1,2);
  c33.Divide(1,2);
  c34.Divide(2,1);
  c35.Divide(2,1);
  c36.Divide(2,1);
  c37.Divide(2,4);
  c38.Divide(1,2);
  c39.Divide(2,4);
  c40.Divide(1,2);
  c41.Divide(2,4);
  c42.Divide(1,2);
  c43.Divide(2,1);
  c44.Divide(2,1);
  c45.Divide(2,2);

  int offset=0;

  for(int i=0; i<4; i++)
  {
    if(i<2) offset=0;
    else offset=2;

    c1.cd(i+offset+1+2);
    plotting_ratio(i,0);
    gPad->SetLogx();
    c1.Update();
    c1.cd(i+offset+1);
    plotting_device(i,0);
    gPad->SetLogx();
    fKinematicsRD[i][0]->GetXaxis()->SetTitle("Q^{2}");
    fKinematicsRD[i][0]->GetYaxis()->SetTitle("Entries");
    c1.Update();

    c2.cd(i+offset+1+2);
    plotting_ratio(i,1);
    gPad->SetLogx();
    c2.Update();
    c2.cd(i+offset+1);
    plotting_device(i,1);
    gPad->SetLogx();
    c2.Update();

    c3.cd(i+offset+1+2);
    plotting_ratio(i,2);
    c3.Update();
    c3.cd(i+offset+1);
    plotting_device(i,2);
    c3.Update();

    c4.cd(i+offset+1+2);
    plotting_ratio(i,3);
    c4.Update();
    c4.cd(i+offset+1);
    plotting_device(i,3);
    c4.Update();

    c5.cd(i+offset+1+2);
    plotting_ratio(i,4);
    c5.Update();
    c5.cd(i+offset+1);
    plotting_device(i,4);
    c5.Update();

    c6.cd(i+offset+1+2);
    plotting_ratio(i,5);
    c6.Update();
    c6.cd(i+offset+1);
    plotting_device(i,5);
    c6.Update();

    c14.cd(i+offset+1+2);
    plotting_ratio(i,6);
    c14.Update();
    c14.cd(i+offset+1);
    plotting_device(i,6);
    c14.Update();

    c15.cd(i+offset+1+2);
    plotting_ratio(i,7);
    c15.Update();
    c15.cd(i+offset+1);
    plotting_device(i,7);
    c15.Update();

    c16.cd(i+offset+1+2);
    plotting_ratio(i,8);
    c16.Update();
    c16.cd(i+offset+1);
    plotting_device(i,8);
    c16.Update();

    c17.cd(i+offset+1+2);
    plotting_ratio(i,9);
    c17.Update();
    c17.cd(i+offset+1);
    plotting_device(i,9);
    c17.Update();

    c18.cd(i+offset+1+2);
    plotting_ratio(i,10);
    c18.Update();
    c18.cd(i+offset+1);
    plotting_device(i,10);
    c18.Update();

    c19.cd(i+offset+1+2);
    plotting_ratio(i,12);
    c19.Update();
    c19.cd(i+offset+1);
    plotting_device(i,12);
    c19.Update();

    c20.cd(i+offset+1+2);
    plotting_ratio(i,13);
    c20.Update();
    c20.cd(i+offset+1);
    plotting_device(i,13);
    c20.Update();

    c21.cd(i+offset+1+2);
    plotting_ratio(i,14);
    c21.Update();
    c21.cd(i+offset+1);
    plotting_device(i,14);
    c21.Update();

    c25.cd(i+offset+1+2);
    plotting_ratio(i,15);
    c25.Update();
    c25.cd(i+offset+1);
    plotting_device(i,15);
    c25.Update();

    c27.cd(i+offset+1+2);
    plotting_ratio(i,16);
    c27.Update();
    c27.cd(i+offset+1);
    plotting_device(i,16);
    c27.Update();

    c37.cd(i+offset+1+2);
    plotting_ratio_vertex(i,17);
    c37.Update();
    c37.cd(i+offset+1);
    plotting_device(i,17);
    c37.Update();

    c39.cd(i+offset+1+2);
    plotting_ratio_vertex(i,18);
    c39.Update();
    c39.cd(i+offset+1);
    plotting_device(i,18);
    c39.Update();

    c41.cd(i+offset+1+2);
    plotting_ratio_vertex(i,19);
    c41.Update();
    c41.cd(i+offset+1);
    plotting_device(i,19);
    c41.Update();
  }

  c7.cd(2);
  plotting_ratio(0,11);
  c7.Update();
  c7.cd(1);
  plotting_device(0,11);
  c7.Update();

  c7.Write();

  c8.cd(2);
  plotting_ratio(4,0);
  gPad->SetLogx();
  c8.Update();
  c8.cd(1);
  plotting_device(4,0);
  gPad->SetLogx();
  c8.Update();

  c8.Write();

  c9.cd(2);
  plotting_ratio(4,1);
  gPad->SetLogx();
  c9.Update();
  c9.cd(1);
  plotting_device(4,1);
  gPad->SetLogx();
  c9.Update();

  c9.Write();

  c10.cd(2);
  plotting_ratio(4,2);
  c10.Update();
  c10.cd(1);
  plotting_device(4,2);
  c10.Update();

  c10.Write();

  c11.cd(2);
  plotting_ratio(4,3);
  c11.Update();
  c11.cd(1);
  plotting_device(4,3);
  c11.Update();

  c11.Write();

  c12.cd(2);
  plotting_ratio(4,4);
  c12.Update();
  c12.cd(1);
  plotting_device(4,4);
  c12.Update();

  c12.Write();

  c13.cd(2);
  plotting_ratio(4,5);
  c13.Update();
  c13.cd(1);
  plotting_device(4,5);
  c13.Update();

  c13.Write();

  c33.cd(2);
  plotting_ratio(4,6);
  c33.Update();
  c33.cd(1);
  plotting_device(4,6);
  c33.Update();

  c33.Write();

  c22.cd(2);
  plotting_ratio(4,12);
  c22.Update();
  c22.cd(1);
  plotting_device(4,12);
  c22.Update();

  c22.Write();

  c23.cd(2);
  plotting_ratio(4,13);
  c23.Update();
  c23.cd(1);
  plotting_device(4,13);
  c23.Update();

  c23.Write();

  c24.cd(2);
  plotting_ratio(4,14);
  c24.Update();
  c24.cd(1);
  plotting_device(4,14);
  c24.Update();

  c24.Write();

  c26.cd(2);
  plotting_ratio(4,15);
  c26.Update();
  c26.cd(1);
  plotting_device(4,15);
  c26.Update();

  c26.Write();

  c28.cd(2);
  plotting_ratio(4,16);
  c28.Update();
  c28.cd(1);
  plotting_device(4,16);
  c28.Update();

  c28.Write();

  c38.cd(2);
  plotting_ratio_vertex(4,17);
  c38.Update();
  c38.cd(1);
  plotting_device(4,17);
  c38.Update();

  c40.cd(2);
  plotting_ratio_vertex(4,18);
  c40.Update();
  c40.cd(1);
  plotting_device(4,18);
  c40.Update();

  c42.cd(2);
  plotting_ratio_vertex(4,19);
  c42.Update();
  c42.cd(1);
  plotting_device(4,19);
  c42.Update();

  c34.cd(1);
  fThetaRDp[0]->Draw("COLZ");
  fThetaRDp[0]->GetXaxis()->SetTitle("theta_y");
  fThetaRDp[0]->GetYaxis()->SetTitle("p");
  c34.Update();

  c34.cd(2);
  fThetaMCp[0]->Draw("COLZ");
  fThetaMCp[0]->GetXaxis()->SetTitle("theta_y");
  fThetaMCp[0]->GetYaxis()->SetTitle("p");
  c34.Update();

  c35.cd(1);
  fThetaRDp[1]->Draw("COLZ");
  fThetaRDp[1]->GetXaxis()->SetTitle("theta_x");
  fThetaRDp[1]->GetYaxis()->SetTitle("p");
  c35.Update();

  c35.cd(2);
  fThetaMCp[1]->Draw("COLZ");
  fThetaMCp[1]->GetXaxis()->SetTitle("theta_x");
  fThetaMCp[1]->GetYaxis()->SetTitle("p");
  c35.Update();

  c36.cd(1);
  fThetaRDp[2]->Draw("COLZ");
  fThetaRDp[2]->GetXaxis()->SetTitle("theta_y");
  fThetaRDp[2]->GetYaxis()->SetTitle("theta_x");
  c36.Update();

  c36.cd(2);
  fThetaMCp[2]->Draw("COLZ");
  fThetaMCp[2]->GetXaxis()->SetTitle("theta_y");
  fThetaMCp[2]->GetYaxis()->SetTitle("theta_x");
  c36.Update();

  c43.cd(1);
  fECAL0RD->Draw("COLZ");
  fECAL0RD->GetXaxis()->SetTitle("x");
  fECAL0RD->GetYaxis()->SetTitle("y");
  c43.Update();

  c43.cd(2);
  fECAL0MC->Draw("COLZ");
  fECAL0MC->GetXaxis()->SetTitle("x");
  fECAL0MC->GetYaxis()->SetTitle("y");
  c43.Update();

  c43.Write();

  c44.cd(1);
  fVertexRD[4]->SetLineColor(kRed);
  // fVertexRD[4]->Scale(1/fHadronRD);
  fVertexRD[4]->Draw("");
  fVertexRD[3]->SetLineColor(kGreen);
  // fVertexRD[3]->Scale(1/fHadronRD);
  fVertexRD[3]->Draw("SAMES");
  fVertexRD[0]->SetLineColor(kMagenta);
  // fVertexRD[0]->Scale(1/fHadronRD);
  fVertexRD[0]->Draw("SAMES");
  fVertexRD[1]->SetLineColor(kBlue);
  // fVertexRD[1]->Scale(1/fHadronRD);
  fVertexRD[1]->Draw("SAMES");
  fVertexRD[2]->SetLineColor(kCyan);
  // fVertexRD[2]->Scale(1/fHadronRD);
  fVertexRD[2]->Draw("SAMES");
  c44.Update();

  c44.cd(2);
  fVertexMCb[4]->SetLineColor(kRed);
  fVertexMCb[4]->Scale(fHadronRD/fHadronMC);
  fVertexMCb[4]->Draw("");
  fVertexMCb[3]->SetLineColor(kGreen);
  fVertexMCb[3]->Scale(fHadronRD/fHadronMC);
  fVertexMCb[3]->Draw("SAMES");
  fVertexMCb[0]->SetLineColor(kMagenta);
  fVertexMCb[0]->Scale(fHadronRD/fHadronMC);
  fVertexMCb[0]->Draw("SAMES");
  fVertexMCb[1]->SetLineColor(kBlue);
  fVertexMCb[1]->Scale(fHadronRD/fHadronMC);
  fVertexMCb[1]->Draw("SAMES");
  fVertexMCb[2]->SetLineColor(kCyan);
  fVertexMCb[2]->Scale(fHadronRD/fHadronMC);
  fVertexMCb[2]->Draw("SAMES");
  c44.Update();

  c44.Write();

  c45.cd(1);
  fTarget->SetLineColor(kRed);
  fTarget->Draw("");
  c45.Update();

  c45.cd(2);
  fTargetMC->SetLineColor(kBlue);
  fTargetMC->Draw("");
  c45.Update();

  c45.cd(3);
  fTarget2D->Draw("COLZ");
  c45.Update();

  c45.cd(4);
  fTarget2DMC->Draw("COLZ");
  c45.Update();

  c45.Write();

  c1.Print("kinMCRD.pdf(","pdf");
  c2.Print("kinMCRD.pdf","pdf");
  c3.Print("kinMCRD.pdf","pdf");
  c4.Print("kinMCRD.pdf","pdf");
  c5.Print("kinMCRD.pdf","pdf");
  c6.Print("kinMCRD.pdf","pdf");
  c19.Print("kinMCRD.pdf","pdf");
  c20.Print("kinMCRD.pdf","pdf");
  c21.Print("kinMCRD.pdf","pdf");
  c25.Print("kinMCRD.pdf","pdf");
  c27.Print("kinMCRD.pdf","pdf");
  c7.Print("kinMCRD.pdf","pdf");
  c8.Print("kinMCRD.pdf","pdf");
  c9.Print("kinMCRD.pdf","pdf");
  c10.Print("kinMCRD.pdf","pdf");
  c11.Print("kinMCRD.pdf","pdf");
  c12.Print("kinMCRD.pdf","pdf");
  c13.Print("kinMCRD.pdf","pdf");
  c33.Print("kinMCRD.pdf","pdf");
  c34.Print("kinMCRD.pdf","pdf");
  c35.Print("kinMCRD.pdf","pdf");
  c36.Print("kinMCRD.pdf","pdf");
  c22.Print("kinMCRD.pdf","pdf");
  c23.Print("kinMCRD.pdf","pdf");
  c24.Print("kinMCRD.pdf","pdf");
  c26.Print("kinMCRD.pdf","pdf");
  c28.Print("kinMCRD.pdf","pdf");
  c14.Print("kinMCRD.pdf","pdf");
  c15.Print("kinMCRD.pdf","pdf");
  c16.Print("kinMCRD.pdf","pdf");
  c17.Print("kinMCRD.pdf","pdf");
  c18.Print("kinMCRD.pdf","pdf");
  c37.Print("kinMCRD.pdf","pdf");
  c38.Print("kinMCRD.pdf","pdf");
  c39.Print("kinMCRD.pdf","pdf");
  c40.Print("kinMCRD.pdf","pdf");
  c41.Print("kinMCRD.pdf","pdf");
  c42.Print("kinMCRD.pdf","pdf");
  c43.Print("kinMCRD.pdf","pdf");
  c44.Print("kinMCRD.pdf","pdf");
  c45.Print("kinMCRD.pdf)","pdf");
}

void MCextraction(string pFilelist)
{

  //Kinematics
  Double_t Q2 = 0;
  Double_t xBj = 0;
  Double_t yBj = 0;
  Double_t zBj = 0;
  Double_t wBj = 0;
  Double_t nu = 0;


  // Target cells
  if(Y2012) InitTargetFile(target_file_2012);
  else if(Y2016) InitTargetFile(target_file_2016);

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
    TBranch *Charge = (TBranch*) tree->FindBranch("Charge");
    TBranch *XX0 = (TBranch*) tree->FindBranch("XX0");
    TBranch *HM04x = (TBranch*) tree->FindBranch("HM04x");
    TBranch *HM04y = (TBranch*) tree->FindBranch("HM04y");
    TBranch *HM05x = (TBranch*) tree->FindBranch("HM05x");
    TBranch *HM05y = (TBranch*) tree->FindBranch("HM05y");
    TBranch *HL04x = (TBranch*) tree->FindBranch("HL04x");
    TBranch *HL04y = (TBranch*) tree->FindBranch("HL04y");
    TBranch *HL05x = (TBranch*) tree->FindBranch("HL05x");
    TBranch *HL05y = (TBranch*) tree->FindBranch("HL05y");
    TBranch *HO03x = (TBranch*) tree->FindBranch("HO03x");
    TBranch *HO03y = (TBranch*) tree->FindBranch("HO03y");
    TBranch *HO04x = (TBranch*) tree->FindBranch("HO04x");
    TBranch *HO04y = (TBranch*) tree->FindBranch("HO04y");
    TBranch *HG01x = (TBranch*) tree->FindBranch("HG01x");
    TBranch *HG01y = (TBranch*) tree->FindBranch("HG01y");
    TBranch *HG021x = (TBranch*) tree->FindBranch("HG021x");
    TBranch *HG021y = (TBranch*) tree->FindBranch("HG021y");
    TBranch *HG022x = (TBranch*) tree->FindBranch("HG022x");
    TBranch *HG022y = (TBranch*) tree->FindBranch("HG022y");
    TBranch *saved = (TBranch*) tree->FindBranch("saved");
    TBranch *BPV = (TBranch*) tree->FindBranch("BPV");
    TBranch *isMuPrim = (TBranch*) tree->FindBranch("isMuPrim");
    TBranch *MZfirst = (TBranch*) tree->FindBranch("MZfirst");
    TBranch *beam_chi2 = (TBranch*) tree->FindBranch("beam_chi2");
    TBranch *mu_prim_chi2 = (TBranch*) tree->FindBranch("mu_prim_chi2");
    TBranch *cellsCrossed = (TBranch*) tree->FindBranch("cellsCrossed");
    TBranch *BMS = (TBranch*) tree->FindBranch("BMS");

    //Hadrons
    TBranch *p = (TBranch*) tree->FindBranch("Hadrons.P");
    TBranch *pt = (TBranch*) tree->FindBranch("Hadrons.pt");
    TBranch *th = (TBranch*) tree->FindBranch("Hadrons.th");
    TBranch *ph = (TBranch*) tree->FindBranch("Hadrons.ph");
    TBranch *ph_pl = (TBranch*) tree->FindBranch("Hadrons.ph_pl");
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
    TBranch *chi2_hadron = (TBranch*) tree->FindBranch("Hadrons.chi2_hadron");
    TBranch *HZfirst = (TBranch*) tree->FindBranch("Hadrons.HZfirst");
    TBranch *HZlast = (TBranch*) tree->FindBranch("Hadrons.HZlast");

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
    TBranch *MC_HM04x = (TBranch*) tree->FindBranch("MC_HM04x");
    TBranch *MC_HM04y = (TBranch*) tree->FindBranch("MC_HM04y");
    TBranch *MC_HM05x = (TBranch*) tree->FindBranch("MC_HM05x");
    TBranch *MC_HM05y = (TBranch*) tree->FindBranch("MC_HM05y");
    TBranch *MC_HL04x = (TBranch*) tree->FindBranch("MC_HL04x");
    TBranch *MC_HL04y = (TBranch*) tree->FindBranch("MC_HL04y");
    TBranch *MC_HL05x = (TBranch*) tree->FindBranch("MC_HL05x");
    TBranch *MC_HL05y = (TBranch*) tree->FindBranch("MC_HL05y");
    TBranch *MC_HO03x = (TBranch*) tree->FindBranch("MC_HO03x");
    TBranch *MC_HO03y = (TBranch*) tree->FindBranch("MC_HO03y");
    TBranch *MC_HO04x = (TBranch*) tree->FindBranch("MC_HO04x");
    TBranch *MC_HO04y = (TBranch*) tree->FindBranch("MC_HO04y");
    TBranch *MC_HG01x = (TBranch*) tree->FindBranch("MC_HG01x");
    TBranch *MC_HG01y = (TBranch*) tree->FindBranch("MC_HG01y");
    TBranch *MC_HG021x = (TBranch*) tree->FindBranch("MC_HG021x");
    TBranch *MC_HG021y = (TBranch*) tree->FindBranch("MC_HG021y");
    TBranch *MC_HG022x = (TBranch*) tree->FindBranch("MC_HG022x");
    TBranch *MC_HG022y = (TBranch*) tree->FindBranch("MC_HG022y");
    TBranch *recons = (TBranch*) tree->FindBranch("recons");
    TBranch *MC_yTr = (TBranch*) tree->FindBranch("MC_yTr");
    TBranch *MC_xTr = (TBranch*) tree->FindBranch("MC_xTr");
    TBranch *MC_TCx = (TBranch*) tree->FindBranch("MC_TCx");
    TBranch *MC_TCy = (TBranch*) tree->FindBranch("MC_TCy");

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
      Charge->GetEntry(ip);
      XX0->GetEntry(ip);
      HM04x->GetEntry(ip);
      HM04y->GetEntry(ip);
      HM05x->GetEntry(ip);
      HM05y->GetEntry(ip);
      HL04x->GetEntry(ip);
      HL04y->GetEntry(ip);
      HL05x->GetEntry(ip);
      HL05y->GetEntry(ip);
      HO03x->GetEntry(ip);
      HO03y->GetEntry(ip);
      HO04x->GetEntry(ip);
      HO04y->GetEntry(ip);
      HG01x->GetEntry(ip);
      HG01y->GetEntry(ip);
      HG021x->GetEntry(ip);
      HG021y->GetEntry(ip);
      HG022x->GetEntry(ip);
      HG022y->GetEntry(ip);
      saved->GetEntry(ip);
      BPV->GetEntry(ip);
      isMuPrim->GetEntry(ip);
      MZfirst->GetEntry(ip);
      beam_chi2->GetEntry(ip);
      mu_prim_chi2->GetEntry(ip);
      cellsCrossed->GetEntry(ip);
      BMS->GetEntry(ip);

      //Hadrons
      p->GetEntry(ip);
      pt->GetEntry(ip);
      th->GetEntry(ip);
      ph->GetEntry(ip);
      ph_pl->GetEntry(ip);
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
      chi2_hadron->GetEntry(ip);
      HZfirst->GetEntry(ip);
      HZlast->GetEntry(ip);

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
      MC_HM04x->GetEntry(ip);
      MC_HM04y->GetEntry(ip);
      MC_HM05x->GetEntry(ip);
      MC_HM05y->GetEntry(ip);
      MC_HL04x->GetEntry(ip);
      MC_HL04y->GetEntry(ip);
      MC_HL05x->GetEntry(ip);
      MC_HL05y->GetEntry(ip);
      MC_HO03x->GetEntry(ip);
      MC_HO03y->GetEntry(ip);
      MC_HO04x->GetEntry(ip);
      MC_HO04y->GetEntry(ip);
      MC_HG01x->GetEntry(ip);
      MC_HG01y->GetEntry(ip);
      MC_HG021x->GetEntry(ip);
      MC_HG021y->GetEntry(ip);
      MC_HG022x->GetEntry(ip);
      MC_HG022y->GetEntry(ip);
      MC_yTr->GetEntry(ip);
      MC_xTr->GetEntry(ip);
      MC_TCx->GetEntry(ip);
      MC_TCy->GetEntry(ip);

      //MCHadrons
      MC_p->GetEntry(ip);
      MC_th->GetEntry(ip);
      MC_ph->GetEntry(ip);
      MC_charge->GetEntry(ip);
      MC_pid->GetEntry(ip);
      MC_recons->GetEntry(ip);
      MC_recHadIdx->GetEntry(ip);

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


      //2006 ---


      // -----------------------------------------------------------------------
      // --------- DIS Selection -----------------------------------------------
      // -----------------------------------------------------------------------

      double theta_m = asin(sqrt(pow(p1x->GetLeaf("p1x")->GetValue()/sqrt(pow(E_mu_prim->GetLeaf("E_mu_prim")->GetValue(),2)-pow(fM_mu,2)),2)+pow(p1y->GetLeaf("p1y")->GetValue()/sqrt(pow(E_mu_prim->GetLeaf("E_mu_prim")->GetValue(),2)-pow(fM_mu,2)),2)));
      double phi_m = asin(p1x->GetLeaf("p1x")->GetValue()/sqrt(pow(p1x->GetLeaf("p1x")->GetValue(),2)+pow(p1y->GetLeaf("p1y")->GetValue(),2)));
      double thetay_b = atan(p0y->GetLeaf("p0y")->GetValue()/sqrt(pow(E_beam->GetLeaf("E_beam")->GetValue(),2)-pow(fM_mu,2)));
      double thetax_b = atan(p0x->GetLeaf("p0x")->GetValue()/sqrt(pow(E_beam->GetLeaf("E_beam")->GetValue(),2)-pow(fM_mu,2)));


      // -----------------------------------------------------------------------
      //  Data -----------------------------------------------------------------
      // -----------------------------------------------------------------------

      fAllDISflag = 0;

      // Best Primary Vertex

      // Reconstructed muon
      if((0<E_beam->GetLeaf("E_beam")->GetValue()))
      {

        //BMS (reconstructed beam track)
        if(true) //not used in acceptance
        {

          // Energy of the muon beam
          if((140<E_beam->GetLeaf("E_beam")->GetValue() && E_beam->GetLeaf("E_beam")->GetValue()<180 && !(E_beam->GetLeaf("E_beam")->GetValue()==160)))
          {

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

                  // Cells crossing
                  if(true)
                  {

                    // IM/O triggers
                    if((trig&8 || trig&256))
                    {

                      // Q2 cut
                      if((Q2>1))
                      {

                        // y cut
                        if((fYmin<yBj && yBj<fYmax))
                        {

                          // W cut
                          if((fWmin<sqrt(wBj) && sqrt(wBj)<fWmax))
                          {

                            // x cut
                            if((fXmin<xBj && xBj<fXmax))
                            {
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
              if(InTarget(x->GetLeaf("x")->GetValue(),y->GetLeaf("y")->GetValue(),z->GetLeaf("z")->GetValue()))
              {

                // Cells crossing
                if(true)
                {

                  if((trig&2 || trig&4 || trig&8))
                  {

                    // Q2 cut
                    if((Q2>1))
                    {

                      // y cut
                      if((fYmin<yBj && yBj<fYmax))
                      {

                        // W cut
                        if((fWmin<sqrt(wBj) && sqrt(wBj)<fWmax))
                        {
                          if((fXmin<xBj && xBj<fXmax))
                          {
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
              if( (inTarget->GetLeaf("inTarget")->GetValue())
                  && (-325<z->GetLeaf("z")->GetValue() && z->GetLeaf("z")->GetValue()<-71))
              {

                if((beam_chi2->GetLeaf("beam_chi2")->GetValue()<10))
                {

                  // Cells crossing
                  if((cellsCrossed->GetLeaf("cellsCrossed")->GetValue()))
                  {

                    if((mu_prim_chi2->GetLeaf("mu_prim_chi2")->GetValue()<10))
                    {

                      if((MZfirst->GetLeaf("MZfirst")->GetValue()<350))
                      {

                        if((trig&2 || trig&4 || trig&8 || trig&512))
                        {

                          // Q2 cut
                          if((Q2>1))
                          {
                            fMuMC[4].push_back(E_beam->GetLeaf("E_beam")->GetValue());
                            fThetaMCMu[2].push_back(sqrt(pow(p0x->GetLeaf("p0x")->GetValue(),2)
                                                        +pow(p0y->GetLeaf("p0y")->GetValue(),2)
                                                        +pow(p0z->GetLeaf("p0z")->GetValue(),2)));
                            fThetaMCMu[1].push_back(thetax_b);
                            fThetaMCMu[0].push_back(thetay_b);

                            // y cut
                            if((fYmin<yBj && yBj<fYmax))
                            {

                              // W cut
                              if((fWmin<sqrt(wBj) && sqrt(wBj)<fWmax))
                              {
                                if((fXmin<xBj && xBj<fXmax))
                                {
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
              }
            }
            //2016 ---
          }
        }
      }

      // If all DIS tests are good, then event is saved
      if(fAllDISflag)
      {
        fNEventsMC++;

        fTargetMC->Fill(z->GetLeaf("z")->GetValue());
        fTarget2DMC->Fill(z->GetLeaf("z")->GetValue(),y->GetLeaf("y")->GetValue());
        // MT
        if(int(trig&2) && !int(trig&4) && !int(trig&8) && !int(trig&512))
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
        // LT
        if(int(trig&4) && !int(trig&2) && !int(trig&8)&& !int(trig&512))
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
        // OT
        if(int(trig&8) && !int(trig&2) && !int(trig&4) && !int(trig&512))
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
        // LAST
        if(int(trig&512) && !int(trig&4) && !int(trig&8) && !int(trig&2))
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

        // ALL TRIGGERS
        if(int(trig&2) || int(trig&4) || int(trig&8) || int(trig&512))
        {
          fQ2kinMC[4].push_back(Q2);
          fXBjkinMC[4].push_back(xBj);
          fYBjkinMC[4].push_back(yBj);
          fWBjkinMC[4].push_back(sqrt(wBj));
          fNukinMC[4].push_back(nu);
          // fMuMC[4].push_back(E_beam->GetLeaf("E_beam")->GetValue());
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

      // x Binning

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

      if(0.<yBj && yBj<0.15) ybin = 0;
      else if(0.15<=yBj && yBj<0.2) ybin = 1;
      else if(0.2<=yBj && yBj<0.3) ybin = 2;
      else if(0.3<=yBj && yBj<0.5) ybin = 3;
      else ybin = 4;

      // -----------------------------------------------------------------------
      // -----------------------------------------------------------------------
      // --------- Hadrons Selection -------------------------------------------
      // -----------------------------------------------------------------------
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
          }
          else
          {
            zBj = 0;
          }

          // /phi_plane for electron (Radiative correction test for electro-production from real photons)
          // Has to be done before Hadron cuts
          if(0.1<zBj && (fId==8 || fId==9))
            fKinematicsMC[0][11]->Fill(abs(ph_pl->GetLeaf("Hadrons.ph_pl")->GetValue(i)));

          // Maximum radiation length cumulated
          if(!(hXX0->GetLeaf("Hadrons.XX0")->GetValue(i) < 15)) continue;

          // Chi2/ndf
          if(!(chi2_hadron->GetLeaf("Hadrons.chi2_hadron")->GetValue(i) < 10)) continue;

          // Zfirst
          if(!(HZfirst->GetLeaf("Hadrons.HZfirst")->GetValue(i)<350)) continue;

          // Zlast
          //if(!(350<HZlast->GetLeaf("Hadrons.HZlast")->GetValue(i))) continue;

          // Momentum cut (12 GeV to 40 GeV, increasing to 3 GeV to 40 GeV)
          if(!(fPmin<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<fPmax)) continue;

          // Theta cut
          if(!(0.01<thRICH->GetLeaf("Hadrons.thRICH")->GetValue(i) && thRICH->GetLeaf("Hadrons.thRICH")->GetValue(i)<0.12)) continue;

          // RICH position cut
          if(!(pow(RICHx->GetLeaf("Hadrons.RICHx")->GetValue(i),2)+pow(RICHy->GetLeaf("Hadrons.RICHy")->GetValue(i),2)>25)) continue;

          // z cut
          if(!(0.2<zBj && zBj<0.85)) continue;

          int dz = abs(z->GetLeaf("z")->GetValue()-70);
          int ydy = y->GetLeaf("y")->GetValue()+dz*tan(th->GetLeaf("Hadrons.th")->GetValue(i))*sin(ph->GetLeaf("Hadrons.ph")->GetValue(i));
          int xdx = x->GetLeaf("x")->GetValue()+dz*tan(th->GetLeaf("Hadrons.th")->GetValue(i))*cos(ph->GetLeaf("Hadrons.ph")->GetValue(i));
          // if(!( ( -35 < xdx && xdx < 35 ) && ( -25 < ydy && ydy < 25 ) )) continue;
          fECAL0MC->Fill(xdx,ydy);

          if(-325<=z->GetLeaf("z")->GetValue() && z->GetLeaf("z")->GetValue()<-261.5) fVertexMCb[0]->Fill(HZlast->GetLeaf("Hadrons.HZlast")->GetValue(i));
          else if(-261.5<=z->GetLeaf("z")->GetValue() && z->GetLeaf("z")->GetValue()<-198) fVertexMCb[1]->Fill(HZlast->GetLeaf("Hadrons.HZlast")->GetValue(i));
          else if(-198<=z->GetLeaf("z")->GetValue() && z->GetLeaf("z")->GetValue()<-134.5) fVertexMCb[2]->Fill(HZlast->GetLeaf("Hadrons.HZlast")->GetValue(i));
          else if(-134.5<=z->GetLeaf("z")->GetValue() && z->GetLeaf("z")->GetValue()<=-71) fVertexMCb[3]->Fill(HZlast->GetLeaf("Hadrons.HZlast")->GetValue(i));
          fVertexMCb[4]->Fill(HZlast->GetLeaf("Hadrons.HZlast")->GetValue(i));
          fHadronMC++;

          // fXBjkinMC[4].push_back(xBj);

          // MT
          if(int(trig&2) && !int(trig&4) && !int(trig&8) && !int(trig&512))
          {
            fKinematicsMC[0][3]->Fill(zBj);
            fKinematicsMC[0][12]->Fill(p->GetLeaf("Hadrons.P")->GetValue(i));
            fKinematicsMC[0][13]->Fill(th->GetLeaf("Hadrons.th")->GetValue(i));
            fKinematicsMC[0][14]->Fill(ph->GetLeaf("Hadrons.ph")->GetValue(i));
            fKinematicsMC[0][15]->Fill(ph_pl->GetLeaf("Hadrons.ph_pl")->GetValue(i));
            fKinematicsMC[0][16]->Fill(pt->GetLeaf("Hadrons.pt")->GetValue(i));
            fKinematicsMC[0][17]->Fill(z->GetLeaf("z")->GetValue());
            if(fId == 8 || fId == 9) fKinematicsMC[0][18]->Fill(z->GetLeaf("z")->GetValue());
            else fKinematicsMC[0][19]->Fill(z->GetLeaf("z")->GetValue());
          }
          // LT
          if(int(trig&4) && !int(trig&2) && !int(trig&8) && !int(trig&512))
          {
            fKinematicsMC[1][3]->Fill(zBj);
            fKinematicsMC[1][12]->Fill(p->GetLeaf("Hadrons.P")->GetValue(i));
            fKinematicsMC[1][13]->Fill(th->GetLeaf("Hadrons.th")->GetValue(i));
            fKinematicsMC[1][14]->Fill(ph->GetLeaf("Hadrons.ph")->GetValue(i));
            fKinematicsMC[1][15]->Fill(ph_pl->GetLeaf("Hadrons.ph_pl")->GetValue(i));
            fKinematicsMC[1][16]->Fill(pt->GetLeaf("Hadrons.pt")->GetValue(i));
            fKinematicsMC[1][17]->Fill(z->GetLeaf("z")->GetValue());
            if(fId == 8 || fId == 9) fKinematicsMC[1][18]->Fill(z->GetLeaf("z")->GetValue());
            else fKinematicsMC[1][19]->Fill(z->GetLeaf("z")->GetValue());
          }
          // OT
          if(int(trig&8) && !int(trig&2) && !int(trig&4) && !int(trig&512))
          {
            fKinematicsMC[2][3]->Fill(zBj);
            fKinematicsMC[2][12]->Fill(p->GetLeaf("Hadrons.P")->GetValue(i));
            fKinematicsMC[2][13]->Fill(th->GetLeaf("Hadrons.th")->GetValue(i));
            fKinematicsMC[2][14]->Fill(ph->GetLeaf("Hadrons.ph")->GetValue(i));
            fKinematicsMC[2][15]->Fill(ph_pl->GetLeaf("Hadrons.ph_pl")->GetValue(i));
            fKinematicsMC[2][16]->Fill(pt->GetLeaf("Hadrons.pt")->GetValue(i));
            fKinematicsMC[2][17]->Fill(z->GetLeaf("z")->GetValue());
            if(fId == 8 || fId == 9) fKinematicsMC[2][18]->Fill(z->GetLeaf("z")->GetValue());
            else fKinematicsMC[2][19]->Fill(z->GetLeaf("z")->GetValue());
          }
          // LAST
          if(int(trig&512) && !int(trig&4) && !int(trig&8) && !int(trig&2))
          {
            fKinematicsMC[3][3]->Fill(zBj);
            fKinematicsMC[3][12]->Fill(p->GetLeaf("Hadrons.P")->GetValue(i));
            fKinematicsMC[3][13]->Fill(th->GetLeaf("Hadrons.th")->GetValue(i));
            fKinematicsMC[3][14]->Fill(ph->GetLeaf("Hadrons.ph")->GetValue(i));
            fKinematicsMC[3][15]->Fill(ph_pl->GetLeaf("Hadrons.ph_pl")->GetValue(i));
            fKinematicsMC[3][16]->Fill(pt->GetLeaf("Hadrons.pt")->GetValue(i));
            fKinematicsMC[3][17]->Fill(z->GetLeaf("z")->GetValue());
            if(fId == 8 || fId == 9) fKinematicsMC[3][18]->Fill(z->GetLeaf("z")->GetValue());
            else fKinematicsMC[3][19]->Fill(z->GetLeaf("z")->GetValue());
          }

          // ALL TRIGGERS
          if(int(trig&2) || int(trig&4) || int(trig&8) || int(trig&512))
          {
            fKinematicsMC[4][3]->Fill(zBj);
            fKinematicsMC[4][12]->Fill(p->GetLeaf("Hadrons.P")->GetValue(i));
            fKinematicsMC[4][13]->Fill(th->GetLeaf("Hadrons.th")->GetValue(i));
            fKinematicsMC[4][14]->Fill(ph->GetLeaf("Hadrons.ph")->GetValue(i));
            fKinematicsMC[4][15]->Fill(ph_pl->GetLeaf("Hadrons.ph_pl")->GetValue(i));
            fKinematicsMC[4][16]->Fill(pt->GetLeaf("Hadrons.pt")->GetValue(i));
            fKinematicsMC[4][17]->Fill(z->GetLeaf("z")->GetValue());
            if(fId == 8 || fId == 9) fKinematicsMC[4][18]->Fill(z->GetLeaf("z")->GetValue());
            else fKinematicsMC[4][19]->Fill(z->GetLeaf("z")->GetValue());
          }
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
    // fKinematicsMC[4][1]->Fill(fXBjkinMC[4][i]);
    fKinematicsMC[4][2]->Fill(fYBjkinMC[4][i]);
    fKinematicsMC[4][4]->Fill(fWBjkinMC[4][i]);
    fKinematicsMC[4][5]->Fill(fNukinMC[4][i]);
    // fKinematicsMC[4][6]->Fill(fMuMC[4][i]);
    fKinematicsMC[4][7]->Fill(fMupMC[4][i]);
    fKinematicsMC[4][8]->Fill(fThetaMC[4][i]);
    fKinematicsMC[4][9]->Fill(fPhiMC[4][i]);
    fKinematicsMC[4][10]->Fill(fVertexMC[4][i]);
  }
  for(int i=0; i<int(fMuMC[4].size()); i++)
  {
    fKinematicsMC[4][6]->Fill(fMuMC[4][i]);
  }
  for(int i=0; i<int(fXBjkinMC[4].size()); i++)
  {
    fKinematicsMC[4][1]->Fill(fXBjkinMC[4][i]);
  }
  for(int i=0; i<int(fThetaMCMu[0].size()); i++)
  {
    fThetaMCp[0]->Fill(fThetaMCMu[0][i],fThetaMCMu[2][i]);
    fThetaMCp[1]->Fill(fThetaMCMu[1][i],fThetaMCMu[2][i]);
    fThetaMCp[2]->Fill(fThetaMCMu[0][i],fThetaMCMu[1][i]);
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

  // Target cells
  if(Y2012) InitTargetFile(target_file_2012);
  else if(Y2016) InitTargetFile(target_file_2016);

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
    TBranch *Charge = (TBranch*) tree->FindBranch("Charge");
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
    TBranch *BPV = (TBranch*) tree->FindBranch("BPV");
    TBranch *isMuPrim = (TBranch*) tree->FindBranch("isMuPrim");
    TBranch *MZfirst = (TBranch*) tree->FindBranch("MZfirst");
    TBranch *beam_chi2 = (TBranch*) tree->FindBranch("beam_chi2");
    TBranch *mu_prim_chi2 = (TBranch*) tree->FindBranch("mu_prim_chi2");
    TBranch *cellsCrossed = (TBranch*) tree->FindBranch("cellsCrossed");
    TBranch *BMS = (TBranch*) tree->FindBranch("BMS");

    //Hadrons
    TBranch *p = (TBranch*) tree->FindBranch("Hadrons.P");
    TBranch *pt = (TBranch*) tree->FindBranch("Hadrons.pt");
    TBranch *th = (TBranch*) tree->FindBranch("Hadrons.th");
    TBranch *ph = (TBranch*) tree->FindBranch("Hadrons.ph");
    TBranch *ph_pl = (TBranch*) tree->FindBranch("Hadrons.ph_pl");
    TBranch *hXX0 = (TBranch*) tree->FindBranch("Hadrons.XX0");
    TBranch *inHCALacc = (TBranch*) tree->FindBranch("Hadrons.inHCALacc");
    TBranch *HCAL = (TBranch*) tree->FindBranch("Hadrons.HCAL");
    TBranch *charge = (TBranch*) tree->FindBranch("Hadrons.charge");
    TBranch *thRICH = (TBranch*) tree->FindBranch("Hadrons.thRICH");
    TBranch *thC = (TBranch*) tree->FindBranch("Hadrons.thC");
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
    TBranch *chi2_hadron = (TBranch*) tree->FindBranch("Hadrons.chi2_hadron");
    TBranch *HZfirst = (TBranch*) tree->FindBranch("Hadrons.HZfirst");
    TBranch *HZlast = (TBranch*) tree->FindBranch("Hadrons.HZlast");

    // Loopy loop over the events
    Int_t N = (Int_t) tree->GetEntries();

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
      Charge->GetEntry(ip);
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
      BPV->GetEntry(ip);
      isMuPrim->GetEntry(ip);
      MZfirst->GetEntry(ip);
      beam_chi2->GetEntry(ip);
      mu_prim_chi2->GetEntry(ip);
      cellsCrossed->GetEntry(ip);
      BMS->GetEntry(ip);

      //Hadrons
      p->GetEntry(ip);
      pt->GetEntry(ip);
      th->GetEntry(ip);
      ph->GetEntry(ip);
      ph_pl->GetEntry(ip);
      hXX0->GetEntry(ip);
      inHCALacc->GetEntry(ip);
      HCAL->GetEntry(ip);
      charge->GetEntry(ip);
      thRICH->GetEntry(ip);
      thC->GetEntry(ip);
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
      chi2_hadron->GetEntry(ip);
      HZfirst->GetEntry(ip);
      HZlast->GetEntry(ip);

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

      double theta_m = asin(sqrt(pow(p1x->GetLeaf("p1x")->GetValue()/sqrt(pow(E_mu_prim->GetLeaf("E_mu_prim")->GetValue(),2)-pow(fM_mu,2)),2)+pow(p1y->GetLeaf("p1y")->GetValue()/sqrt(pow(E_mu_prim->GetLeaf("E_mu_prim")->GetValue(),2)-pow(fM_mu,2)),2)));
      double phi_m = asin(p1x->GetLeaf("p1x")->GetValue()/sqrt(pow(p1x->GetLeaf("p1x")->GetValue(),2)+pow(p1y->GetLeaf("p1y")->GetValue(),2)));
      double thetay_b = asin(p0y->GetLeaf("p0y")->GetValue()/sqrt(pow(E_beam->GetLeaf("E_beam")->GetValue(),2)-pow(fM_mu,2)));
      double thetax_b = atan(p0x->GetLeaf("p0x")->GetValue()/sqrt(pow(E_beam->GetLeaf("E_beam")->GetValue(),2)-pow(fM_mu,2)));

      // Best Primary Vertex

      // IsMuPrim
      if(!(0<isMuPrim->GetLeaf("isMuPrim")->GetValue())) continue;

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
        if(!InTarget(x->GetLeaf("x")->GetValue(),y->GetLeaf("y")->GetValue(),z->GetLeaf("z")->GetValue())) continue;
      }
      //2012 ---
      //2016 ---
      else if(Y2016)
      {
        if(!inTarget->GetLeaf("inTarget")->GetValue()) continue;
        if(!(-325<z->GetLeaf("z")->GetValue() && z->GetLeaf("z")->GetValue()<-71)) continue;
      }
      //2016 ---

      // Energy of the muon beam
      if(!(140<E_beam->GetLeaf("E_beam")->GetValue() && E_beam->GetLeaf("E_beam")->GetValue()<180)) continue;

      // BMS
      if(!(BMS->GetLeaf("BMS")->GetValue()>3)) continue;

      // Chi2 beam
      if(!(beam_chi2->GetLeaf("beam_chi2")->GetValue()<10)) continue;

      // Cells crossing
      if(!(cellsCrossed->GetLeaf("cellsCrossed")->GetValue())) continue;

      if(!(mu_prim_chi2->GetLeaf("mu_prim_chi2")->GetValue()<10)) continue;

      if(!(MZfirst->GetLeaf("MZfirst")->GetValue()<350)) continue;

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
        if(!(int(trig&2) || int(trig&4) || int(trig&8) || int(trig&512))) continue;
      }
      //2016 ---

      // Q2 cut
      if(!(Q2>1)) continue;
      fMu[4].push_back(E_beam->GetLeaf("E_beam")->GetValue());

      fThetaMu[2].push_back(sqrt(pow(p0x->GetLeaf("p0x")->GetValue(),2)
                                +pow(p0y->GetLeaf("p0y")->GetValue(),2)
                                +pow(p0z->GetLeaf("p0z")->GetValue(),2)));
      fThetaMu[1].push_back(thetax_b);
      fThetaMu[0].push_back(thetay_b);

      // y cut
      if(!(fYmin<yBj && yBj<fYmax)) continue;

      // W cut
      if(!(fWmin<sqrt(wBj) && sqrt(wBj)<fWmax)) continue;

      // x cut
      if(!(fXmin<xBj && xBj<fXmax)) continue;

      fTarget->Fill(z->GetLeaf("z")->GetValue());
      fTarget2D->Fill(z->GetLeaf("z")->GetValue(),y->GetLeaf("y")->GetValue());

      fNEventsRD++;
      // MT
      if(int(trig&2) && !int(trig&4) && !int(trig&8) && !int(trig&512))
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
      // LT
      if(int(trig&4) && !int(trig&2) && !int(trig&8)&& !int(trig&512))
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
      // OT
      if(int(trig&8) && !int(trig&2) && !int(trig&4) && !int(trig&512))
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
      // LAST
      if(int(trig&512) && !int(trig&4) && !int(trig&8) && !int(trig&2))
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

      // ALL TRIGGERS
      // if(trig&2 || trig&4 || trig&8)
      if(int(trig&2) || int(trig&4) || int(trig&8) || int(trig&512))
      {
        fQ2kin[4].push_back(Q2);
        fXBjkin[4].push_back(xBj);
        fYBjkin[4].push_back(yBj);
        fWBjkin[4].push_back(sqrt(wBj));
        fNukin[4].push_back(nu);
        // fMu[4].push_back(E_beam->GetLeaf("E_beam")->GetValue());
        fMup[4].push_back(E_mu_prim->GetLeaf("E_mu_prim")->GetValue());
        fTheta[4].push_back(theta_m);
        fPhi[4].push_back(phi_m);
        fVertex[4].push_back(z->GetLeaf("z")->GetValue());
      }
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

      // -------------------------------------------------------------------------
      // --------- Hadrons Selection ---------------------------------------------
      // -------------------------------------------------------------------------

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
        }
        // Charge -

        else if(charge->GetLeaf("Hadrons.charge")->GetValue(i) == -1)
        {
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

        // /phi_plane for electron (Radiative correction test for electro-production from real photons)
        // Has to be done before Hadron cuts
        if(0.1<zBj && (LH->GetLeaf("Hadrons.LH")->GetValue(3+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(4+6*i)) && (LH->GetLeaf("Hadrons.LH")->GetValue(3+6*i)>fLHsec_tab[3]))
          fKinematicsRD[0][11]->Fill(abs(ph_pl->GetLeaf("Hadrons.ph_pl")->GetValue(i)));

        // Maximum radiation length cumulated
        if(!(hXX0->GetLeaf("Hadrons.XX0")->GetValue(i) < 15)) continue;

        // Chi2/ndf
        if(!(chi2_hadron->GetLeaf("Hadrons.chi2_hadron")->GetValue(i) < 10)) continue;

        // Zfirst
        if(!(HZfirst->GetLeaf("Hadrons.HZfirst")->GetValue(i)<350)) continue;

        // Zlast
        //if(!(350<HZlast->GetLeaf("Hadrons.HZlast")->GetValue(i))) continue;

        // Theta cut
        if(!(0.01<thRICH->GetLeaf("Hadrons.thRICH")->GetValue(i) && thRICH->GetLeaf("Hadrons.thRICH")->GetValue(i)<0.12)) continue;

        // RICH position cut
        if(!(pow(RICHx->GetLeaf("Hadrons.RICHx")->GetValue(i),2)+pow(RICHy->GetLeaf("Hadrons.RICHy")->GetValue(i),2)>25)) continue;

        // Momentum cut (12 GeV to 40 GeV, increasing to 3 GeV to 40 GeV)
        if(!(fPmin<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<fPmax)) continue;

        // z cut
        if(!(0.2<zBj && zBj<0.85)) continue;

        int dz = abs(z->GetLeaf("z")->GetValue()-70);
        int ydy = y->GetLeaf("y")->GetValue()+dz*tan(th->GetLeaf("Hadrons.th")->GetValue(i))*sin(ph->GetLeaf("Hadrons.ph")->GetValue(i));
        int xdx = x->GetLeaf("x")->GetValue()+dz*tan(th->GetLeaf("Hadrons.th")->GetValue(i))*cos(ph->GetLeaf("Hadrons.ph")->GetValue(i));
        // if(!( ( -35 < xdx && xdx < 35 ) && ( -25 < ydy && ydy < 25 ) )) continue;
        fECAL0RD->Fill(xdx,ydy);

        if(-325<=z->GetLeaf("z")->GetValue() && z->GetLeaf("z")->GetValue()<-261.5) fVertexRD[0]->Fill(HZlast->GetLeaf("Hadrons.HZlast")->GetValue(i));
        else if(-261.5<=z->GetLeaf("z")->GetValue() && z->GetLeaf("z")->GetValue()<-198) fVertexRD[1]->Fill(HZlast->GetLeaf("Hadrons.HZlast")->GetValue(i));
        else if(-198<=z->GetLeaf("z")->GetValue() && z->GetLeaf("z")->GetValue()<-134.5) fVertexRD[2]->Fill(HZlast->GetLeaf("Hadrons.HZlast")->GetValue(i));
        else if(-134.5<=z->GetLeaf("z")->GetValue() && z->GetLeaf("z")->GetValue()<=-71) fVertexRD[3]->Fill(HZlast->GetLeaf("Hadrons.HZlast")->GetValue(i));
        fVertexRD[4]->Fill(HZlast->GetLeaf("Hadrons.HZlast")->GetValue(i));
        fHadronRD++;

        // fXBjkin[4].push_back(xBj);

        // Non null charge
        if(!charge->GetLeaf("Hadrons.charge")->GetValue(i)) continue;

        if(int(trig&2) && !int(trig&4) && !int(trig&8) && !int(trig&512))
        {
          fKinematicsRD[0][3]->Fill(zBj);
          fKinematicsRD[0][12]->Fill(p->GetLeaf("Hadrons.P")->GetValue(i));
          fKinematicsRD[0][13]->Fill(th->GetLeaf("Hadrons.th")->GetValue(i));
          fKinematicsRD[0][14]->Fill(ph->GetLeaf("Hadrons.ph")->GetValue(i));
          fKinematicsRD[0][15]->Fill(ph_pl->GetLeaf("Hadrons.ph_pl")->GetValue(i));
          fKinematicsRD[0][16]->Fill(pt->GetLeaf("Hadrons.pt")->GetValue(i));
          fKinematicsRD[0][17]->Fill(z->GetLeaf("z")->GetValue());
          fKinematicsRD[0][18]->Fill(z->GetLeaf("z")->GetValue());
          fKinematicsRD[0][19]->Fill(z->GetLeaf("z")->GetValue());

        }
        if(int(trig&4) && !int(trig&2) && !int(trig&8)&& !int(trig&512))
        {
          fKinematicsRD[1][3]->Fill(zBj);
          fKinematicsRD[1][12]->Fill(p->GetLeaf("Hadrons.P")->GetValue(i));
          fKinematicsRD[1][13]->Fill(th->GetLeaf("Hadrons.th")->GetValue(i));
          fKinematicsRD[1][14]->Fill(ph->GetLeaf("Hadrons.ph")->GetValue(i));
          fKinematicsRD[1][15]->Fill(ph_pl->GetLeaf("Hadrons.ph_pl")->GetValue(i));
          fKinematicsRD[1][16]->Fill(pt->GetLeaf("Hadrons.pt")->GetValue(i));
          fKinematicsRD[1][17]->Fill(z->GetLeaf("z")->GetValue());
          fKinematicsRD[1][18]->Fill(z->GetLeaf("z")->GetValue());
          fKinematicsRD[1][19]->Fill(z->GetLeaf("z")->GetValue());
        }
        if(int(trig&8) && !int(trig&2) && !int(trig&4) && !int(trig&512))
        {
          fKinematicsRD[2][3]->Fill(zBj);
          fKinematicsRD[2][12]->Fill(p->GetLeaf("Hadrons.P")->GetValue(i));
          fKinematicsRD[2][13]->Fill(th->GetLeaf("Hadrons.th")->GetValue(i));
          fKinematicsRD[2][14]->Fill(ph->GetLeaf("Hadrons.ph")->GetValue(i));
          fKinematicsRD[2][15]->Fill(ph_pl->GetLeaf("Hadrons.ph_pl")->GetValue(i));
          fKinematicsRD[2][16]->Fill(pt->GetLeaf("Hadrons.pt")->GetValue(i));
          fKinematicsRD[2][17]->Fill(z->GetLeaf("z")->GetValue());
          fKinematicsRD[2][18]->Fill(z->GetLeaf("z")->GetValue());
          fKinematicsRD[2][19]->Fill(z->GetLeaf("z")->GetValue());
        }
        if(int(trig&512) && !int(trig&4) && !int(trig&8) && !int(trig&2))
        {
          fKinematicsRD[3][3]->Fill(zBj);
          fKinematicsRD[3][12]->Fill(p->GetLeaf("Hadrons.P")->GetValue(i));
          fKinematicsRD[3][13]->Fill(th->GetLeaf("Hadrons.th")->GetValue(i));
          fKinematicsRD[3][14]->Fill(ph->GetLeaf("Hadrons.ph")->GetValue(i));
          fKinematicsRD[3][15]->Fill(ph_pl->GetLeaf("Hadrons.ph_pl")->GetValue(i));
          fKinematicsRD[3][16]->Fill(pt->GetLeaf("Hadrons.pt")->GetValue(i));
          fKinematicsRD[3][17]->Fill(z->GetLeaf("z")->GetValue());
          fKinematicsRD[3][18]->Fill(z->GetLeaf("z")->GetValue());
          fKinematicsRD[3][19]->Fill(z->GetLeaf("z")->GetValue());
        }

        // if(trig&2 || trig&4 || trig&8)
        if(int(trig&2) || int(trig&4) || int(trig&8) || int(trig&512))
        {
          fKinematicsRD[4][3]->Fill(zBj);
          fKinematicsRD[4][12]->Fill(p->GetLeaf("Hadrons.P")->GetValue(i));
          fKinematicsRD[4][13]->Fill(th->GetLeaf("Hadrons.th")->GetValue(i));
          fKinematicsRD[4][14]->Fill(ph->GetLeaf("Hadrons.ph")->GetValue(i));
          fKinematicsRD[4][15]->Fill(ph_pl->GetLeaf("Hadrons.ph_pl")->GetValue(i));
          fKinematicsRD[4][16]->Fill(pt->GetLeaf("Hadrons.pt")->GetValue(i));
          fKinematicsRD[4][17]->Fill(z->GetLeaf("z")->GetValue());
          fKinematicsRD[4][18]->Fill(z->GetLeaf("z")->GetValue());
          fKinematicsRD[4][19]->Fill(z->GetLeaf("z")->GetValue());
        }
      }

      //Misc
      fQ2.push_back(Q2);
      fXBj.push_back(xBj);
      fYBj.push_back(yBj);
      fWBj.push_back(wBj);
      fNu.push_back(nu);
    }

    cout << "\n" << endl;
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
      // fKinematicsRD[4][1]->Fill(fXBjkin[4][i]);
      fKinematicsRD[4][2]->Fill(fYBjkin[4][i]);
      fKinematicsRD[4][4]->Fill(fWBjkin[4][i]);
      fKinematicsRD[4][5]->Fill(fNukin[4][i]);
      // fKinematicsRD[4][6]->Fill(fMu[4][i]);
      fKinematicsRD[4][7]->Fill(fMup[4][i]);
      fKinematicsRD[4][8]->Fill(fTheta[4][i]);
      fKinematicsRD[4][9]->Fill(fPhi[4][i]);
      fKinematicsRD[4][10]->Fill(fVertex[4][i]);
  }
  for(int i=0; i<int(fMu[4].size()); i++)
  {
    fKinematicsRD[4][6]->Fill(fMu[4][i]);
  }
  for(int i=0; i<int(fXBjkin[4].size()); i++)
  {
    fKinematicsRD[4][1]->Fill(fXBjkin[4][i]);
  }
  for(int i=0; i<int(fThetaMu[0].size()); i++)
  {
    fThetaRDp[0]->Fill(fThetaMu[0][i],fThetaMu[2][i]);
    fThetaRDp[1]->Fill(fThetaMu[1][i],fThetaMu[2][i]);
    fThetaRDp[2]->Fill(fThetaMu[0][i],fThetaMu[1][i]);
  }
}

int main(int argc, char **argv)
{

  if(argc < 3)
  {
    cout << "ERROR : Not enough arguments." << endl;
    cout << "Asked : 3 *** Received : " << argc-1 << endl;
    cout << "./compMCRD [RD filelist] [MC filelist] [Cutfile]" << endl;

    return 1;
  }

  create_kin_plots();
  readKinCuts(argv[3]);
  cout << "... Real Data treatment ..." << endl;
  RDextraction(argv[1]);
  cout << "... Monte-Carlo treatment ..." << endl;
  MCextraction(argv[2]);
  cout << "... Saving plots ..." << endl;
  save_kin_plots();

  cout << "\n\n";
  cout << "             ********* Event distribution within MT/LT/OT/LAST in percentage of total ********* " << endl;
  cout << "             ---------------------------------------------------------------------------------- " << endl;

  cout << "\n ==> Real Data <==" << endl;

  cout << '|' << setw(15) << "Type" << '|' << setw(15) << "All" << '|' << setw(15) << "MT" << '|' << setw(15) << "LT" << '|' << setw(15) << "OT" << '|' << setw(15) << "LAST" << '|' << endl;
  cout << '|' << setw(15) << "DIS" << '|' << setw(15) << fCountingRD[4][0] << '|' << setw(15) << float(fCountingRD[0][0])/float(fCountingRD[4][0])*100
                                                                    << '|' << setw(15) << float(fCountingRD[1][0])/float(fCountingRD[4][0])*100
                                                                    << '|' << setw(15) << float(fCountingRD[2][0])/float(fCountingRD[4][0])*100
                                                                    << '|' << setw(15) << float(fCountingRD[3][0])/float(fCountingRD[4][0])*100 << '|' << endl;
  cout << '|' << setw(15) << "Hadron" << '|' << setw(15) << fCountingRD[4][3] << '|' << setw(15) << float(fCountingRD[0][3])/float(fCountingRD[4][3])*100
                                                                    << '|' << setw(15) << float(fCountingRD[1][3])/float(fCountingRD[4][3])*100
                                                                    << '|' << setw(15) << float(fCountingRD[2][3])/float(fCountingRD[4][3])*100
                                                                    << '|' << setw(15) << float(fCountingRD[3][3])/float(fCountingRD[4][3])*100 << '|' << endl;

  cout << "\n ==> Monte Carlo <==" << endl;

  cout << '|' << setw(15) << "Type" << '|' << setw(15) << "All" << '|' << setw(15) << "MT" << '|' << setw(15) << "LT" << '|' << setw(15) << "OT" << '|' << setw(15) << "LAST" << '|' << endl;
  cout << '|' << setw(15) << "DIS" << '|' << setw(15) << fCountingMC[4][0] << '|' << setw(15) << float(fCountingMC[0][0])/float(fCountingMC[4][0])*100
                                                                    << '|' << setw(15) << float(fCountingMC[1][0])/float(fCountingMC[4][0])*100
                                                                    << '|' << setw(15) << float(fCountingMC[2][0])/float(fCountingMC[4][0])*100
                                                                    << '|' << setw(15) << float(fCountingMC[3][0])/float(fCountingMC[4][0])*100 << '|' << endl;
  cout << '|' << setw(15) << "Hadron" << '|' << setw(15) << fCountingMC[4][3] << '|' << setw(15) << float(fCountingMC[0][3])/float(fCountingMC[4][3])*100
                                                                    << '|' << setw(15) << float(fCountingMC[1][3])/float(fCountingMC[4][3])*100
                                                                    << '|' << setw(15) << float(fCountingMC[2][3])/float(fCountingMC[4][3])*100
                                                                    << '|' << setw(15) << float(fCountingMC[3][3])/float(fCountingMC[4][3])*100 << '|' << endl;

  return 0;
}

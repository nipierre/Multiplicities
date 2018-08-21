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
    fKinematicsRD[i][0] = new TH1F(Form("Q^{2} %s",trigname[i].c_str()), Form("Q^{2} %s",trigname[i].c_str()), 50, -1, 2);
    fKinematicsRD[i][1] = new TH1F(Form("x_{Bj} %s",trigname[i].c_str()), Form("x_{Bj} %s",trigname[i].c_str()), 50, -3, 0);
    fKinematicsRD[i][2] = new TH1F(Form("y %s",trigname[i].c_str()), Form("y %s",trigname[i].c_str()), 50, 0, 1);
    fKinematicsRD[i][3] = new TH1F(Form("z %s",trigname[i].c_str()), Form("z %s",trigname[i].c_str()), 50, 0, 1);
    fKinematicsRD[i][4] = new TH1F(Form("W %s",trigname[i].c_str()), Form("W %s",trigname[i].c_str()), 50, 2, 18);
    fKinematicsRD[i][5] = new TH1F(Form("#nu %s",trigname[i].c_str()), Form("#nu %s",trigname[i].c_str()), 50, 0, 160);
    fKinematicsRD[i][6] = new TH1F(Form("E_{#mu} %s",trigname[i].c_str()), Form("E_{#mu} %s",trigname[i].c_str()), 50, 140, 180);
    fKinematicsRD[i][7] = new TH1F(Form("E_{#mu'} %s",trigname[i].c_str()), Form("E_{#mu'} %s",trigname[i].c_str()), 50, 0, 160);
    fKinematicsRD[i][8] = new TH1F(Form("#theta %s",trigname[i].c_str()), Form("#theta %s",trigname[i].c_str()), 50, 0, 0.05);
    fKinematicsRD[i][9] = new TH1F(Form("#phi %s",trigname[i].c_str()), Form("#phi %s",trigname[i].c_str()), 50, -1.7, 1.7);
    fKinematicsRD[i][10] = new TH1F(Form("Vertex %s",trigname[i].c_str()), Form("Vertex %s",trigname[i].c_str()), 50, -320, -70);
    fKinematicsRD[i][12] = new TH1F(Form("p_{hadron+e} %s",trigname[i].c_str()), Form("p_{hadron+e} %s",trigname[i].c_str()), 50, 0, 40);
    fKinematicsRD[i][13] = new TH1F(Form("#theta_{hadron+e} %s",trigname[i].c_str()), Form("#theta_{hadron+e} %s",trigname[i].c_str()), 50, 0, 0.25);
    fKinematicsRD[i][14] = new TH1F(Form("#phi_{hadron+e,lab} %s",trigname[i].c_str()), Form("#phi_{hadron+e,lab} %s",trigname[i].c_str()), 50, -3.5, 3.5);
    fKinematicsRD[i][15] = new TH1F(Form("#phi_{hadron+e,prod.pl} %s",trigname[i].c_str()), Form("#phi_{hadron+e,prod.pl} %s",trigname[i].c_str()), 50, 0, 3.5);
    fKinematicsRD[i][16] = new TH1F(Form("p_{T} %s",trigname[i].c_str()), Form("p_{T} %s",trigname[i].c_str()), 50, 0, 3);
    fKinematicsMC[i][0] = new TH1F(Form("Q^{2} Ratio %s",trigname[i].c_str()), Form("Q^{2} Ratio %s",trigname[i].c_str()), 50, -1, 2);
    fKinematicsMC[i][1] = new TH1F(Form("x_{Bj} Ratio %s",trigname[i].c_str()), Form("x_{Bj} Ratio %s",trigname[i].c_str()), 50, -3, 0);
    fKinematicsMC[i][2] = new TH1F(Form("y Ratio %s",trigname[i].c_str()), Form("y Ratio %s",trigname[i].c_str()), 50, 0, 1);
    fKinematicsMC[i][3] = new TH1F(Form("z Ratio %s",trigname[i].c_str()), Form("z Ratio %s",trigname[i].c_str()), 50, 0, 1);
    fKinematicsMC[i][4] = new TH1F(Form("W Ratio %s",trigname[i].c_str()), Form("W Ratio %s",trigname[i].c_str()), 50, 2, 18);
    fKinematicsMC[i][5] = new TH1F(Form("#nu Ratio %s",trigname[i].c_str()), Form("#nu Ratio %s",trigname[i].c_str()), 50, 0, 160);
    fKinematicsMC[i][6] = new TH1F(Form("E_{#mu} Ratio %s",trigname[i].c_str()), Form("E_{#mu} Ratio %s",trigname[i].c_str()), 50, 140, 180);
    fKinematicsMC[i][7] = new TH1F(Form("E_{#mu'} Ratio %s",trigname[i].c_str()), Form("E_{#mu'} Ratio %s",trigname[i].c_str()), 50, 0, 160);
    fKinematicsMC[i][8] = new TH1F(Form("#theta Ratio %s",trigname[i].c_str()), Form("#theta Ratio %s",trigname[i].c_str()), 50, 0, 0.05);
    fKinematicsMC[i][9] = new TH1F(Form("#phi Ratio %s",trigname[i].c_str()), Form("#phi Ratio %s",trigname[i].c_str()), 50, -1.7, 1.7);
    fKinematicsMC[i][10] = new TH1F(Form("Vertex Ratio %s",trigname[i].c_str()), Form("Vertex Ratio %s",trigname[i].c_str()), 50, -320, -70);
    fKinematicsMC[i][12] = new TH1F(Form("p_{hadron+e} Ratio %s",trigname[i].c_str()), Form("p_{hadron+e} Ratio %s",trigname[i].c_str()), 50, 0, 40);
    fKinematicsMC[i][13] = new TH1F(Form("#theta_{hadron+e} Ratio %s",trigname[i].c_str()), Form("#theta_{hadron+e} Ratio %s",trigname[i].c_str()), 50, 0, 0.25);
    fKinematicsMC[i][14] = new TH1F(Form("#phi_{hadron+e,lab} Ratio %s",trigname[i].c_str()), Form("#phi_{hadron+e,lab} Ratio %s",trigname[i].c_str()), 50, -3.5, 3.5);
    fKinematicsMC[i][15] = new TH1F(Form("#phi_{hadron+e,prod.pl} Ratio %s",trigname[i].c_str()), Form("#phi_{hadron+e,prod.pl} Ratio %s",trigname[i].c_str()), 50, 0, 3.5);
    fKinematicsMC[i][16] = new TH1F(Form("p_{T} Ratio %s",trigname[i].c_str()), Form("p_{T} Ratio %s",trigname[i].c_str()), 50, 0, 3);
    BinLogX(fKinematicsRD[i][0]);
    BinLogX(fKinematicsMC[i][0]);
    BinLogX(fKinematicsRD[i][1]);
    BinLogX(fKinematicsMC[i][1]);
  }
  fKinematicsRD[0][11] = new TH1F("#phi_{e,prod.pl}","#phi_{e,prod.pl}", 50, 0, 3.5);
  fKinematicsMC[0][11] = new TH1F("#phi_{e,prod.pl} Ratio","#phi_{e,prod.pl} Ratio", 50, 0, 3.5);
  fThetaRDp[0] = new TH2F("#theta_y RD", "#theta_y RD", 100, -0.005, 0.005, 100, 140, 180);
  fThetaRDp[1] = new TH2F("#theta_x RD", "#theta_x RD", 100, -0.005, 0.005, 100, 140, 180);
  fThetaRDp[2] = new TH2F("#theta_xy RD", "#theta_xy RD", 100, -0.005, 0.005, 100, -0.005, 0.005);
  fThetaMCp[0] = new TH2F("#theta_y MC", "#theta_y MC", 100, -0.005, 0.005, 100, 140, 180);
  fThetaMCp[1] = new TH2F("#theta_x MC", "#theta_x MC", 100, -0.005, 0.005, 100, 140, 180);
  fThetaMCp[2] = new TH2F("#theta_xy MC", "#theta_xy MC", 100, -0.005, 0.005, 100, -0.005, 0.005);
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
    l1[10][i] = new TLine(-320,0.4+i*0.2,-70,0.4+i*0.2);
    l1[11][i] = new TLine(0,0.4+i*0.2,3.5,0.4+i*0.2);
    l1[12][i] = new TLine(0,0.4+i*0.2,40,0.4+i*0.2);
    l1[13][i] = new TLine(0,0.4+i*0.2,0.25,0.4+i*0.2);
    l1[14][i] = new TLine(-3.5,0.4+i*0.2,3.5,0.4+i*0.2);
    l1[15][i] = new TLine(0,0.4+i*0.2,3.5,0.4+i*0.2);
    l1[16][i] = new TLine(0,0.4+i*0.2,3,0.4+i*0.2);
    for(int j=0; j<17; j++)
    {
      l1[j][i]->SetLineStyle(fLineStyle[i]);
      l1[j][i]->SetLineWidth(1);
    }
  }
}

void save_kin_plots()
{
  c1.Divide(2,2);
  c2.Divide(2,2);
  c3.Divide(2,2);
  c4.Divide(2,2);
  c5.Divide(2,2);
  c6.Divide(2,2);
  c7.Divide(1,1);
  c8.Divide(1,1);
  c9.Divide(1,1);
  c10.Divide(1,1);
  c11.Divide(1,1);
  c12.Divide(1,1);
  c13.Divide(1,1);
  c14.Divide(2,2);
  c15.Divide(2,2);
  c16.Divide(2,2);
  c17.Divide(2,2);
  c18.Divide(2,2);
  c19.Divide(2,2);
  c20.Divide(2,2);
  c21.Divide(2,2);
  c22.Divide(1,1);
  c23.Divide(1,1);
  c24.Divide(1,1);
  c25.Divide(2,2);
  c26.Divide(1,1);
  c27.Divide(2,2);
  c28.Divide(1,1);
  c29.Divide(1,1);
  c30.Divide(1,1);
  c31.Divide(1,1);
  c32.Divide(1,1);
  c33.Divide(1,1);
  c34.Divide(2,1);
  c35.Divide(2,1);
  c36.Divide(2,1);

  for(int i=0; i<4; i++)
  {
    cout << endl;
    // i=int(i/2);
    // if(i%2)
    // {
    //   c1.cd(i+3+int(i/2)*2);
    //   for(int tt=0; tt<fKinematicsRD[i][0]->GetNbinsX(); tt++)
    //   {
    //     fError.push_back((fKinematicsRD[i][0]->GetBinError(tt) && fKinematicsMC[i][0]->GetBinError(tt) ? sqrt(pow(1/fKinematicsRD[i][0]->GetBinError(tt),2)+pow(1/fKinematicsMC[i][0]->GetBinError(tt),2)):0));
    //   }
    //   fKinematicsRD[i][0]->Scale(1/fKinematicsRD[2][0]->GetEntries());
    //   fKinematicsMC[i][0]->Scale(1/fKinematicsMC[2][0]->GetEntries());
    //   fKinematicsRatio[i][0] = (TH1F*)fKinematicsRD[i][0]->Clone();
    //   fKinematicsRatio[i][0]->SetStats(0);
    //   fKinematicsRatio[i][0]->Divide(fKinematicsMC[i][0]);
    //   for(int tt=0; tt<fKinematicsRatio[i][0]->GetNbinsX(); tt++)
    //   {
    //     fKinematicsRatio[i][0]->SetBinError(tt,fError[tt]);
    //   }
    //   fError.clear();
    //   fKinematicsRatio[i][0]->SetMarkerStyle(21);
    //   fKinematicsRatio[i][0]->SetFillColor(kYellow-7);
    //   fKinematicsRatio[i][0]->SetMaximum(2.);
    //   fKinematicsRatio[i][0]->SetMinimum(0.);
    //   fKinematicsRatio[i][0]->Draw("PE2");
    //   fKinematicsRatio[i][0]->GetXaxis()->SetLabelSize(0.08);
    //   fKinematicsRatio[i][0]->GetYaxis()->SetLabelSize(0.08);
    //   fKinematicsRatio[i][0]->GetYaxis()->SetNdivisions(2,kFALSE);
    //   for(int tt=0; tt<7; tt++)
    //   {
    //     l1[0][tt]->Draw();
    //   }
    //   gPad->SetLogx();
    //   c1.Update();
    //   cout << ".";
    //   c2.cd(i+3+int(i/2)*2);
    //   for(int tt=0; tt<fKinematicsRD[i][1]->GetNbinsX(); tt++)
    //   {
    //     fError.push_back((fKinematicsRD[i][1]->GetBinError(tt) && fKinematicsMC[i][1]->GetBinError(tt) ? sqrt(pow(1/fKinematicsRD[i][1]->GetBinError(tt),2)+pow(1/fKinematicsMC[i][1]->GetBinError(tt),2)):0));
    //   }
    //   fKinematicsRD[i][1]->Scale(1/fKinematicsRD[2][1]->GetEntries());
    //   fKinematicsMC[i][1]->Scale(1/fKinematicsMC[2][1]->GetEntries());
    //   fKinematicsRatio[i][1] = (TH1F*)fKinematicsRD[i][1]->Clone();
    //   fKinematicsRatio[i][1]->SetStats(0);
    //   fKinematicsRatio[i][1]->Divide(fKinematicsMC[i][1]);
    //   for(int tt=0; tt<fKinematicsRatio[i][1]->GetNbinsX(); tt++)
    //   {
    //     fKinematicsRatio[i][1]->SetBinError(tt,fError[tt]);
    //   }
    //   fError.clear();
    //   fKinematicsRatio[i][1]->SetMarkerStyle(21);
    //   fKinematicsRatio[i][1]->SetFillColor(kYellow-7);
    //   fKinematicsRatio[i][1]->SetMaximum(2.);
    //   fKinematicsRatio[i][1]->SetMinimum(0.);
    //   fKinematicsRatio[i][1]->Draw("PE2");
    //   fKinematicsRatio[i][1]->GetXaxis()->SetLabelSize(0.08);
    //   fKinematicsRatio[i][1]->GetYaxis()->SetLabelSize(0.08);
    //   fKinematicsRatio[i][1]->GetYaxis()->SetNdivisions(2,kFALSE);
    //   for(int tt=0; tt<7; tt++)
    //   {
    //     l1[1][tt]->Draw();
    //   }
    //   gPad->SetLogx();
    //   c2.Update();
    //   cout << ".";
    //   c3.cd(i+3+int(i/2)*2);
    //   for(int tt=0; tt<fKinematicsRD[i][2]->GetNbinsX(); tt++)
    //   {
    //     fError.push_back((fKinematicsRD[i][2]->GetBinError(tt) && fKinematicsMC[i][2]->GetBinError(tt) ? sqrt(pow(1/fKinematicsRD[i][2]->GetBinError(tt),2)+pow(1/fKinematicsMC[i][2]->GetBinError(tt),2)):0));
    //   }
    //   fKinematicsRD[i][2]->Scale(1/fKinematicsRD[2][2]->GetEntries());
    //   fKinematicsMC[i][2]->Scale(1/fKinematicsMC[2][2]->GetEntries());
    //   fKinematicsRatio[i][2] = (TH1F*)fKinematicsRD[i][2]->Clone();
    //   fKinematicsRatio[i][2]->SetStats(0);
    //   fKinematicsRatio[i][2]->Divide(fKinematicsMC[i][2]);
    //   for(int tt=0; tt<fKinematicsRatio[i][2]->GetNbinsX(); tt++)
    //   {
    //     fKinematicsRatio[i][2]->SetBinError(tt,fError[tt]);
    //   }
    //   fError.clear();
    //   fKinematicsRatio[i][2]->SetMarkerStyle(21);
    //   fKinematicsRatio[i][2]->SetFillColor(kYellow-7);
    //   fKinematicsRatio[i][2]->SetMaximum(2.);
    //   fKinematicsRatio[i][2]->SetMinimum(0.);
    //   fKinematicsRatio[i][2]->Draw("PE2");
    //   fKinematicsRatio[i][2]->GetXaxis()->SetLabelSize(0.08);
    //   fKinematicsRatio[i][2]->GetYaxis()->SetLabelSize(0.08);
    //   fKinematicsRatio[i][2]->GetYaxis()->SetNdivisions(2,kFALSE);
    //   for(int tt=0; tt<7; tt++)
    //   {
    //     l1[2][tt]->Draw();
    //   }
    //   c3.Update();
    //   cout << ".";
    //   c4.cd(i+3+int(i/2)*2);
    //   for(int tt=0; tt<fKinematicsRD[i][3]->GetNbinsX(); tt++)
    //   {
    //     fError.push_back((fKinematicsRD[i][3]->GetBinError(tt) && fKinematicsMC[i][3]->GetBinError(tt) ? sqrt(pow(1/fKinematicsRD[i][3]->GetBinError(tt),2)+pow(1/fKinematicsMC[i][3]->GetBinError(tt),2)):0));
    //   }
    //   fKinematicsRD[i][3]->Scale(1/fKinematicsRD[2][3]->GetEntries());
    //   fKinematicsMC[i][3]->Scale(1/fKinematicsMC[2][3]->GetEntries());
    //   fKinematicsRatio[i][3] = (TH1F*)fKinematicsRD[i][3]->Clone();
    //   fKinematicsRatio[i][3]->SetStats(0);
    //   fKinematicsRatio[i][3]->Divide(fKinematicsMC[i][3]);
    //   for(int tt=0; tt<fKinematicsRatio[i][3]->GetNbinsX(); tt++)
    //   {
    //     fKinematicsRatio[i][3]->SetBinError(tt,fError[tt]);
    //   }
    //   fError.clear();
    //   fKinematicsRatio[i][3]->SetMarkerStyle(21);
    //   fKinematicsRatio[i][3]->SetFillColor(kYellow-7);
    //   fKinematicsRatio[i][3]->SetMaximum(2.);
    //   fKinematicsRatio[i][3]->SetMinimum(0.);
    //   fKinematicsRatio[i][3]->Draw("PE2");
    //   fKinematicsRatio[i][3]->GetXaxis()->SetLabelSize(0.08);
    //   fKinematicsRatio[i][3]->GetYaxis()->SetLabelSize(0.08);
    //   fKinematicsRatio[i][3]->GetYaxis()->SetNdivisions(2,kFALSE);
    //   for(int tt=0; tt<7; tt++)
    //   {
    //     l1[3][tt]->Draw();
    //   }
    //   c4.Update();
    //   cout << ".";
    //   c5.cd(i+3+int(i/2)*2);
    //   for(int tt=0; tt<fKinematicsRD[i][4]->GetNbinsX(); tt++)
    //   {
    //     fError.push_back((fKinematicsRD[i][4]->GetBinError(tt) && fKinematicsMC[i][4]->GetBinError(tt) ? sqrt(pow(1/fKinematicsRD[i][4]->GetBinError(tt),2)+pow(1/fKinematicsMC[i][4]->GetBinError(tt),2)):0));
    //   }
    //   fKinematicsRD[i][4]->Scale(1/fKinematicsRD[2][4]->GetEntries());
    //   fKinematicsMC[i][4]->Scale(1/fKinematicsMC[2][4]->GetEntries());
    //   fKinematicsRatio[i][4] = (TH1F*)fKinematicsRD[i][4]->Clone();
    //   fKinematicsRatio[i][4]->SetStats(0);
    //   fKinematicsRatio[i][4]->Divide(fKinematicsMC[i][4]);
    //   for(int tt=0; tt<fKinematicsRatio[i][4]->GetNbinsX(); tt++)
    //   {
    //     fKinematicsRatio[i][4]->SetBinError(tt,fError[tt]);
    //   }
    //   fError.clear();
    //   fKinematicsRatio[i][4]->SetMarkerStyle(21);
    //   fKinematicsRatio[i][4]->SetFillColor(kYellow-7);
    //   fKinematicsRatio[i][4]->SetMaximum(2.);
    //   fKinematicsRatio[i][4]->SetMinimum(0.);
    //   fKinematicsRatio[i][4]->Draw("PE2");
    //   fKinematicsRatio[i][4]->GetXaxis()->SetLabelSize(0.08);
    //   fKinematicsRatio[i][4]->GetYaxis()->SetLabelSize(0.08);
    //   fKinematicsRatio[i][4]->GetYaxis()->SetNdivisions(2,kFALSE);
    //   for(int tt=0; tt<7; tt++)
    //   {
    //     l1[4][tt]->Draw();
    //   }
    //   c5.Update();
    //   cout << ".";
    //   c6.cd(i+3+int(i/2)*2);
    //   for(int tt=0; tt<fKinematicsRD[i][5]->GetNbinsX(); tt++)
    //   {
    //     fError.push_back((fKinematicsRD[i][5]->GetBinError(tt) && fKinematicsMC[i][5]->GetBinError(tt) ? sqrt(pow(1/fKinematicsRD[i][5]->GetBinError(tt),2)+pow(1/fKinematicsMC[i][5]->GetBinError(tt),2)):0));
    //   }
    //   fKinematicsRD[i][5]->Scale(1/fKinematicsRD[2][5]->GetEntries());
    //   fKinematicsMC[i][5]->Scale(1/fKinematicsMC[2][5]->GetEntries());
    //   fKinematicsRatio[i][5] = (TH1F*)fKinematicsRD[i][5]->Clone();
    //   fKinematicsRatio[i][5]->SetStats(0);
    //   fKinematicsRatio[i][5]->Divide(fKinematicsMC[i][5]);
    //   for(int tt=0; tt<fKinematicsRatio[i][5]->GetNbinsX(); tt++)
    //   {
    //     fKinematicsRatio[i][5]->SetBinError(tt,fError[tt]);
    //   }
    //   fError.clear();
    //   fKinematicsRatio[i][5]->SetMarkerStyle(21);
    //   fKinematicsRatio[i][5]->SetFillColor(kYellow-7);
    //   fKinematicsRatio[i][5]->SetMaximum(2.);
    //   fKinematicsRatio[i][5]->SetMinimum(0.);
    //   fKinematicsRatio[i][5]->Draw("PE2");
    //   fKinematicsRatio[i][5]->GetXaxis()->SetLabelSize(0.08);
    //   fKinematicsRatio[i][5]->GetYaxis()->SetLabelSize(0.08);
    //   fKinematicsRatio[i][5]->GetYaxis()->SetNdivisions(2,kFALSE);
    //   for(int tt=0; tt<7; tt++)
    //   {
    //     l1[5][tt]->Draw();
    //   }
    //   c6.Update();
    //   cout << ".";
    //   c14.cd(i+3+int(i/2)*2);
    //   for(int tt=0; tt<fKinematicsRD[i][6]->GetNbinsX(); tt++)
    //   {
    //     fError.push_back((fKinematicsRD[i][6]->GetBinError(tt) && fKinematicsMC[i][6]->GetBinError(tt) ? sqrt(pow(1/fKinematicsRD[i][6]->GetBinError(tt),2)+pow(1/fKinematicsMC[i][6]->GetBinError(tt),2)):0));
    //   }
    //   fKinematicsRD[i][6]->Scale(1/fKinematicsRD[2][6]->GetEntries());
    //   fKinematicsMC[i][6]->Scale(1/fKinematicsMC[2][6]->GetEntries());
    //   fKinematicsRatio[i][6] = (TH1F*)fKinematicsRD[i][6]->Clone();
    //   fKinematicsRatio[i][6]->SetStats(0);
    //   fKinematicsRatio[i][6]->Divide(fKinematicsMC[i][6]);
    //   for(int tt=0; tt<fKinematicsRatio[i][6]->GetNbinsX(); tt++)
    //   {
    //     fKinematicsRatio[i][6]->SetBinError(tt,fError[tt]);
    //   }
    //   fError.clear();
    //   fKinematicsRatio[i][6]->SetMarkerStyle(21);
    //   fKinematicsRatio[i][6]->SetFillColor(kYellow-7);
    //   fKinematicsRatio[i][6]->SetMaximum(2.);
    //   fKinematicsRatio[i][6]->SetMinimum(0.);
    //   fKinematicsRatio[i][6]->Draw("PE2");
    //   fKinematicsRatio[i][6]->GetXaxis()->SetLabelSize(0.08);
    //   fKinematicsRatio[i][6]->GetYaxis()->SetLabelSize(0.08);
    //   fKinematicsRatio[i][6]->GetYaxis()->SetNdivisions(2,kFALSE);
    //   for(int tt=0; tt<7; tt++)
    //   {
    //     l1[6][tt]->Draw();
    //   }
    //   c14.Update();
    //   cout << ".";
    //   c15.cd(i+3+int(i/2)*2);
    //   for(int tt=0; tt<fKinematicsRD[i][7]->GetNbinsX(); tt++)
    //   {
    //     fError.push_back((fKinematicsRD[i][7]->GetBinError(tt) && fKinematicsMC[i][7]->GetBinError(tt) ? sqrt(pow(1/fKinematicsRD[i][7]->GetBinError(tt),2)+pow(1/fKinematicsMC[i][7]->GetBinError(tt),2)):0));
    //   }
    //   fKinematicsRD[i][7]->Scale(1/fKinematicsRD[2][7]->GetEntries());
    //   fKinematicsMC[i][7]->Scale(1/fKinematicsMC[2][7]->GetEntries());
    //   fKinematicsRatio[i][7] = (TH1F*)fKinematicsRD[i][7]->Clone();
    //   fKinematicsRatio[i][7]->SetStats(0);
    //   fKinematicsRatio[i][7]->Divide(fKinematicsMC[i][7]);
    //   for(int tt=0; tt<fKinematicsRatio[i][7]->GetNbinsX(); tt++)
    //   {
    //     fKinematicsRatio[i][7]->SetBinError(tt,fError[tt]);
    //   }
    //   fError.clear();
    //   fKinematicsRatio[i][7]->SetMarkerStyle(21);
    //   fKinematicsRatio[i][7]->SetFillColor(kYellow-7);
    //   fKinematicsRatio[i][7]->SetMaximum(2.);
    //   fKinematicsRatio[i][7]->SetMinimum(0.);
    //   fKinematicsRatio[i][7]->Draw("PE2");
    //   fKinematicsRatio[i][7]->GetXaxis()->SetLabelSize(0.08);
    //   fKinematicsRatio[i][7]->GetYaxis()->SetLabelSize(0.08);
    //   fKinematicsRatio[i][7]->GetYaxis()->SetNdivisions(2,kFALSE);
    //   for(int tt=0; tt<7; tt++)
    //   {
    //     l1[7][tt]->Draw();
    //   }
    //   c15.Update();
    //   cout << ".";
    //   c16.cd(i+3+int(i/2)*2);
    //   for(int tt=0; tt<fKinematicsRD[i][8]->GetNbinsX(); tt++)
    //   {
    //     fError.push_back((fKinematicsRD[i][8]->GetBinError(tt) && fKinematicsMC[i][8]->GetBinError(tt) ? sqrt(pow(1/fKinematicsRD[i][8]->GetBinError(tt),2)+pow(1/fKinematicsMC[i][8]->GetBinError(tt),2)):0));
    //   }
    //   fKinematicsRD[i][8]->Scale(1/fKinematicsRD[2][8]->GetEntries());
    //   fKinematicsMC[i][8]->Scale(1/fKinematicsMC[2][8]->GetEntries());
    //   fKinematicsRatio[i][8] = (TH1F*)fKinematicsRD[i][8]->Clone();
    //   fKinematicsRatio[i][8]->SetStats(0);
    //   fKinematicsRatio[i][8]->Divide(fKinematicsMC[i][8]);
    //   for(int tt=0; tt<fKinematicsRatio[i][8]->GetNbinsX(); tt++)
    //   {
    //     fKinematicsRatio[i][8]->SetBinError(tt,fError[tt]);
    //   }
    //   fError.clear();
    //   fKinematicsRatio[i][8]->SetMarkerStyle(21);
    //   fKinematicsRatio[i][8]->SetFillColor(kYellow-7);
    //   fKinematicsRatio[i][8]->SetMaximum(2.);
    //   fKinematicsRatio[i][8]->SetMinimum(0.);
    //   fKinematicsRatio[i][8]->Draw("PE2");
    //   fKinematicsRatio[i][8]->GetXaxis()->SetLabelSize(0.08);
    //   fKinematicsRatio[i][8]->GetYaxis()->SetLabelSize(0.08);
    //   fKinematicsRatio[i][8]->GetYaxis()->SetNdivisions(2,kFALSE);
    //   for(int tt=0; tt<7; tt++)
    //   {
    //     l1[8][tt]->Draw();
    //   }
    //   c16.Update();
    //   cout << ".";
    //   c17.cd(i+3+int(i/2)*2);
    //   for(int tt=0; tt<fKinematicsRD[i][9]->GetNbinsX(); tt++)
    //   {
    //     fError.push_back((fKinematicsRD[i][9]->GetBinError(tt) && fKinematicsMC[i][9]->GetBinError(tt) ? sqrt(pow(1/fKinematicsRD[i][9]->GetBinError(tt),2)+pow(1/fKinematicsMC[i][9]->GetBinError(tt),2)):0));
    //   }
    //   fKinematicsRD[i][9]->Scale(1/fKinematicsRD[2][9]->GetEntries());
    //   fKinematicsMC[i][9]->Scale(1/fKinematicsMC[2][9]->GetEntries());
    //   fKinematicsRatio[i][9] = (TH1F*)fKinematicsRD[i][9]->Clone();
    //   fKinematicsRatio[i][9]->SetStats(0);
    //   fKinematicsRatio[i][9]->Divide(fKinematicsMC[i][9]);
    //   for(int tt=0; tt<fKinematicsRatio[i][9]->GetNbinsX(); tt++)
    //   {
    //     fKinematicsRatio[i][9]->SetBinError(tt,fError[tt]);
    //   }
    //   fError.clear();
    //   fKinematicsRatio[i][9]->SetMarkerStyle(21);
    //   fKinematicsRatio[i][9]->SetFillColor(kYellow-7);
    //   fKinematicsRatio[i][9]->SetMaximum(2.);
    //   fKinematicsRatio[i][9]->SetMinimum(0.);
    //   fKinematicsRatio[i][9]->Draw("PE2");
    //   fKinematicsRatio[i][9]->GetXaxis()->SetLabelSize(0.08);
    //   fKinematicsRatio[i][9]->GetYaxis()->SetLabelSize(0.08);
    //   fKinematicsRatio[i][9]->GetYaxis()->SetNdivisions(2,kFALSE);
    //   for(int tt=0; tt<7; tt++)
    //   {
    //     l1[9][tt]->Draw();
    //   }
    //   c17.Update();
    //   cout << ".";
    //   c18.cd(i+3+int(i/2)*2);
    //   for(int tt=0; tt<fKinematicsRD[i][10]->GetNbinsX(); tt++)
    //   {
    //     fError.push_back((fKinematicsRD[i][10]->GetBinError(tt) && fKinematicsMC[i][10]->GetBinError(tt) ? sqrt(pow(1/fKinematicsRD[i][10]->GetBinError(tt),2)+pow(1/fKinematicsMC[i][10]->GetBinError(tt),2)):0));
    //   }
    //   fKinematicsRD[i][10]->Scale(1/fKinematicsRD[2][10]->GetEntries());
    //   fKinematicsMC[i][10]->Scale(1/fKinematicsMC[2][10]->GetEntries());
    //   fKinematicsRatio[i][10] = (TH1F*)fKinematicsRD[i][10]->Clone();
    //   fKinematicsRatio[i][10]->SetStats(0);
    //   fKinematicsRatio[i][10]->Divide(fKinematicsMC[i][10]);
    //   for(int tt=0; tt<fKinematicsRatio[i][10]->GetNbinsX(); tt++)
    //   {
    //     fKinematicsRatio[i][10]->SetBinError(tt,fError[tt]);
    //   }
    //   fError.clear();
    //   fKinematicsRatio[i][10]->SetMarkerStyle(21);
    //   fKinematicsRatio[i][10]->SetFillColor(kYellow-7);
    //   fKinematicsRatio[i][10]->SetMaximum(2.);
    //   fKinematicsRatio[i][10]->SetMinimum(0.);
    //   fKinematicsRatio[i][10]->Draw("PE2");
    //   fKinematicsRatio[i][10]->GetXaxis()->SetLabelSize(0.08);
    //   fKinematicsRatio[i][10]->GetYaxis()->SetLabelSize(0.08);
    //   fKinematicsRatio[i][10]->GetYaxis()->SetNdivisions(2,kFALSE);
    //   for(int tt=0; tt<7; tt++)
    //   {
    //     l1[10][tt]->Draw();
    //   }
    //   c18.Update();
    //   cout << ".";
    //   c19.cd(i+3+int(i/2)*2);
    //   for(int tt=0; tt<fKinematicsRD[i][12]->GetNbinsX(); tt++)
    //   {
    //     fError.push_back((fKinematicsRD[i][12]->GetBinError(tt) && fKinematicsMC[i][12]->GetBinError(tt) ? sqrt(pow(1/fKinematicsRD[i][12]->GetBinError(tt),2)+pow(1/fKinematicsMC[i][12]->GetBinError(tt),2)):0));
    //   }
    //   fKinematicsRD[i][12]->Scale(1/fKinematicsRD[2][12]->GetEntries());
    //   fKinematicsMC[i][12]->Scale(1/fKinematicsMC[2][12]->GetEntries());
    //   fKinematicsRatio[i][12] = (TH1F*)fKinematicsRD[i][12]->Clone();
    //   fKinematicsRatio[i][12]->SetStats(0);
    //   fKinematicsRatio[i][12]->Divide(fKinematicsMC[i][12]);
    //   for(int tt=0; tt<fKinematicsRatio[i][12]->GetNbinsX(); tt++)
    //   {
    //     fKinematicsRatio[i][12]->SetBinError(tt,fError[tt]);
    //   }
    //   fError.clear();
    //   fKinematicsRatio[i][12]->SetMarkerStyle(21);
    //   fKinematicsRatio[i][12]->SetFillColor(kYellow-7);
    //   fKinematicsRatio[i][12]->SetMaximum(2.);
    //   fKinematicsRatio[i][12]->SetMinimum(0.);
    //   fKinematicsRatio[i][12]->Draw("PE2");
    //   fKinematicsRatio[i][12]->GetXaxis()->SetLabelSize(0.08);
    //   fKinematicsRatio[i][12]->GetYaxis()->SetLabelSize(0.08);
    //   fKinematicsRatio[i][12]->GetYaxis()->SetNdivisions(2,kFALSE);
    //   for(int tt=0; tt<7; tt++)
    //   {
    //     l1[12][tt]->Draw();
    //   }
    //   c19.Update();
    //   cout << ".";
    //   c20.cd(i+3+int(i/2)*2);
    //   for(int tt=0; tt<fKinematicsRD[i][13]->GetNbinsX(); tt++)
    //   {
    //     fError.push_back((fKinematicsRD[i][13]->GetBinError(tt) && fKinematicsMC[i][13]->GetBinError(tt) ? sqrt(pow(1/fKinematicsRD[i][13]->GetBinError(tt),2)+pow(1/fKinematicsMC[i][13]->GetBinError(tt),2)):0));
    //   }
    //   fKinematicsRD[i][13]->Scale(1/fKinematicsRD[2][13]->GetEntries());
    //   fKinematicsMC[i][13]->Scale(1/fKinematicsMC[2][13]->GetEntries());
    //   fKinematicsRatio[i][13] = (TH1F*)fKinematicsRD[i][13]->Clone();
    //   fKinematicsRatio[i][13]->SetStats(0);
    //   fKinematicsRatio[i][13]->Divide(fKinematicsMC[i][13]);
    //   for(int tt=0; tt<fKinematicsRatio[i][13]->GetNbinsX(); tt++)
    //   {
    //     fKinematicsRatio[i][13]->SetBinError(tt,fError[tt]);
    //   }
    //   fError.clear();
    //   fKinematicsRatio[i][13]->SetMarkerStyle(21);
    //   fKinematicsRatio[i][13]->SetFillColor(kYellow-7);
    //   fKinematicsRatio[i][13]->SetMaximum(2.);
    //   fKinematicsRatio[i][13]->SetMinimum(0.);
    //   fKinematicsRatio[i][13]->Draw("PE2");
    //   fKinematicsRatio[i][13]->GetXaxis()->SetLabelSize(0.08);
    //   fKinematicsRatio[i][13]->GetYaxis()->SetLabelSize(0.08);
    //   fKinematicsRatio[i][13]->GetYaxis()->SetNdivisions(2,kFALSE);
    //   for(int tt=0; tt<7; tt++)
    //   {
    //     l1[13][tt]->Draw();
    //   }
    //   c20.Update();
    //   cout << ".";
    //   c21.cd(i+3+int(i/2)*2);
    //   for(int tt=0; tt<fKinematicsRD[i][14]->GetNbinsX(); tt++)
    //   {
    //     fError.push_back((fKinematicsRD[i][14]->GetBinError(tt) && fKinematicsMC[i][14]->GetBinError(tt) ? sqrt(pow(1/fKinematicsRD[i][14]->GetBinError(tt),2)+pow(1/fKinematicsMC[i][14]->GetBinError(tt),2)):0));
    //   }
    //   fKinematicsRD[i][14]->Scale(1/fKinematicsRD[2][14]->GetEntries());
    //   fKinematicsMC[i][14]->Scale(1/fKinematicsMC[2][14]->GetEntries());
    //   fKinematicsRatio[i][14] = (TH1F*)fKinematicsRD[i][14]->Clone();
    //   fKinematicsRatio[i][14]->SetStats(0);
    //   fKinematicsRatio[i][14]->Divide(fKinematicsMC[i][14]);
    //   for(int tt=0; tt<fKinematicsRatio[i][14]->GetNbinsX(); tt++)
    //   {
    //     fKinematicsRatio[i][14]->SetBinError(tt,fError[tt]);
    //   }
    //   fError.clear();
    //   fKinematicsRatio[i][14]->SetMarkerStyle(21);
    //   fKinematicsRatio[i][14]->SetFillColor(kYellow-7);
    //   fKinematicsRatio[i][14]->SetMaximum(2.);
    //   fKinematicsRatio[i][14]->SetMinimum(0.);
    //   fKinematicsRatio[i][14]->Draw("PE2");
    //   fKinematicsRatio[i][14]->GetXaxis()->SetLabelSize(0.08);
    //   fKinematicsRatio[i][14]->GetYaxis()->SetLabelSize(0.08);
    //   fKinematicsRatio[i][14]->GetYaxis()->SetNdivisions(2,kFALSE);
    //   for(int tt=0; tt<7; tt++)
    //   {
    //     l1[14][tt]->Draw();
    //   }
    //   c21.Update();
    //   cout << ".";
    //   c25.cd(i+3+int(i/2)*2);
    //   for(int tt=0; tt<fKinematicsRD[i][15]->GetNbinsX(); tt++)
    //   {
    //     fError.push_back((fKinematicsRD[i][15]->GetBinError(tt) && fKinematicsMC[i][15]->GetBinError(tt) ? sqrt(pow(1/fKinematicsRD[i][15]->GetBinError(tt),2)+pow(1/fKinematicsMC[i][15]->GetBinError(tt),2)):0));
    //   }
    //   fKinematicsRD[i][15]->Scale(1/fKinematicsRD[2][15]->GetEntries());
    //   fKinematicsMC[i][15]->Scale(1/fKinematicsMC[2][15]->GetEntries());
    //   fKinematicsRatio[i][15] = (TH1F*)fKinematicsRD[i][15]->Clone();
    //   fKinematicsRatio[i][15]->SetStats(0);
    //   fKinematicsRatio[i][15]->Divide(fKinematicsMC[i][15]);
    //   for(int tt=0; tt<fKinematicsRatio[i][15]->GetNbinsX(); tt++)
    //   {
    //     fKinematicsRatio[i][15]->SetBinError(tt,fError[tt]);
    //   }
    //   fError.clear();
    //   fKinematicsRatio[i][15]->SetMarkerStyle(21);
    //   fKinematicsRatio[i][15]->SetFillColor(kYellow-7);
    //   fKinematicsRatio[i][15]->SetMaximum(2.);
    //   fKinematicsRatio[i][15]->SetMinimum(0.);
    //   fKinematicsRatio[i][15]->Draw("PE2");
    //   fKinematicsRatio[i][15]->GetXaxis()->SetLabelSize(0.08);
    //   fKinematicsRatio[i][15]->GetYaxis()->SetLabelSize(0.08);
    //   fKinematicsRatio[i][15]->GetYaxis()->SetNdivisions(2,kFALSE);
    //   for(int tt=0; tt<7; tt++)
    //   {
    //     l1[15][tt]->Draw();
    //   }
    //   c25.Update();
    //   cout << ".";
    //   c27.cd(i+3+int(i/2)*2);
    //   for(int tt=0; tt<fKinematicsRD[i][16]->GetNbinsX(); tt++)
    //   {
    //     fError.push_back((fKinematicsRD[i][16]->GetBinError(tt) && fKinematicsMC[i][16]->GetBinError(tt) ? sqrt(pow(1/fKinematicsRD[i][16]->GetBinError(tt),2)+pow(1/fKinematicsMC[i][16]->GetBinError(tt),2)):0));
    //   }
    //   fKinematicsRD[i][16]->Scale(1/fKinematicsRD[2][16]->GetEntries());
    //   fKinematicsMC[i][16]->Scale(1/fKinematicsMC[2][16]->GetEntries());
    //   fKinematicsRatio[i][16] = (TH1F*)fKinematicsRD[i][16]->Clone();
    //   fKinematicsRatio[i][16]->SetStats(0);
    //   fKinematicsRatio[i][16]->Divide(fKinematicsMC[i][16]);
    //   for(int tt=0; tt<fKinematicsRatio[i][16]->GetNbinsX(); tt++)
    //   {
    //     fKinematicsRatio[i][16]->SetBinError(tt,fError[tt]);
    //   }
    //   fError.clear();
    //   fKinematicsRatio[i][16]->SetMarkerStyle(21);
    //   fKinematicsRatio[i][16]->SetFillColor(kYellow-7);
    //   fKinematicsRatio[i][16]->SetMaximum(2.);
    //   fKinematicsRatio[i][16]->SetMinimum(0.);
    //   fKinematicsRatio[i][16]->Draw("PE2");
    //   fKinematicsRatio[i][16]->GetXaxis()->SetLabelSize(0.08);
    //   fKinematicsRatio[i][16]->GetYaxis()->SetLabelSize(0.08);
    //   fKinematicsRatio[i][16]->GetYaxis()->SetNdivisions(2,kFALSE);
    //   for(int tt=0; tt<7; tt++)
    //   {
    //     l1[16][tt]->Draw();
    //   }
    //   c27.Update();
    // }
    // else
    // {
      c1.cd(i+1);
      fKinematicsMC[i][0]->Scale(fKinematicsRD[i][0]->GetEntries()/fKinematicsMC[i][0]->GetEntries());
      fKinematicsRD[i][0]->SetLineColor(kRed);
      fKinematicsMC[i][0]->SetLineColor(kBlue);
      fKinematicsRD[i][0]->SetMinimum(0.);
      fKinematicsRD[i][0]->SetMaximum(max(fKinematicsRD[i][0]->GetMaximum()*1.2,fKinematicsMC[i][0]->GetMaximum()*1.2));
      fKinematicsRD[i][0]->GetYaxis()->SetNdivisions(304,kTRUE);
      for(int tt=0; tt<fKinematicsRD[i][0]->GetNbinsX(); tt++)
      {
        fKinematicsRD[i][0]->SetBinError(tt,sqrt(fKinematicsRD[i][0]->GetBinContent(tt)));
      }
      fKinematicsRD[i][0]->Draw("E2");
      fKinematicsRD[i][0]->SetMarkerStyle(22);
      fKinematicsRD[i][0]->Draw("PSAME");
      fKinematicsRD[i][0]->GetXaxis()->SetLabelSize(0.03);
      fKinematicsRD[i][0]->GetYaxis()->SetLabelSize(0.03);
      for(int tt=0; tt<fKinematicsMC[i][0]->GetNbinsX(); tt++)
      {
        fKinematicsMC[i][0]->SetBinError(tt,sqrt(fKinematicsMC[i][0]->GetBinContent(tt)));
      }
      fKinematicsMC[i][0]->Draw("E2SAME");
      fKinematicsMC[i][0]->Draw("SAME");
      gPad->SetLogx();
      c1.Update();


      c2.cd(i+1);
      fKinematicsMC[i][1]->Scale(fKinematicsRD[i][1]->GetEntries()/fKinematicsMC[i][1]->GetEntries());
      fKinematicsRD[i][1]->SetLineColor(kRed);
      fKinematicsMC[i][1]->SetLineColor(kBlue);
      fKinematicsRD[i][1]->SetMinimum(0.);
      fKinematicsRD[i][1]->SetMaximum(max(fKinematicsRD[i][1]->GetMaximum()*1.2,fKinematicsMC[i][1]->GetMaximum()*1.2));
      fKinematicsRD[i][1]->GetYaxis()->SetNdivisions(304,kTRUE);
      for(int tt=0; tt<fKinematicsRD[i][1]->GetNbinsX(); tt++)
      {
        fKinematicsRD[i][1]->SetBinError(tt,sqrt(fKinematicsRD[i][1]->GetBinContent(tt)));
      }
      fKinematicsRD[i][1]->Draw("E2");
      fKinematicsRD[i][1]->SetMarkerStyle(22);
      fKinematicsRD[i][1]->Draw("PSAME");
      fKinematicsRD[i][1]->GetXaxis()->SetLabelSize(0.03);
      fKinematicsRD[i][1]->GetYaxis()->SetLabelSize(0.03);
      for(int tt=0; tt<fKinematicsMC[i][1]->GetNbinsX(); tt++)
      {
        fKinematicsMC[i][1]->SetBinError(tt,sqrt(fKinematicsMC[i][1]->GetBinContent(tt)));
      }
      fKinematicsMC[i][1]->Draw("E2SAME");
      fKinematicsMC[i][1]->Draw("SAME");
      gPad->SetLogx();
      c2.Update();

      c3.cd(i+1);
      fKinematicsMC[i][2]->Scale(fKinematicsRD[i][2]->GetEntries()/fKinematicsMC[i][2]->GetEntries());
      fKinematicsRD[i][2]->SetLineColor(kRed);
      fKinematicsMC[i][2]->SetLineColor(kBlue);
      fKinematicsRD[i][2]->SetMinimum(0.);
      fKinematicsRD[i][2]->SetMaximum(max(fKinematicsRD[i][2]->GetMaximum()*1.2,fKinematicsMC[i][2]->GetMaximum()*1.2));
      fKinematicsRD[i][2]->GetYaxis()->SetNdivisions(304,kTRUE);
      for(int tt=0; tt<fKinematicsRD[i][2]->GetNbinsX(); tt++)
      {
        fKinematicsRD[i][2]->SetBinError(tt,sqrt(fKinematicsRD[i][2]->GetBinContent(tt)));
      }
      fKinematicsRD[i][2]->Draw("E2");
      fKinematicsRD[i][2]->SetMarkerStyle(22);
      fKinematicsRD[i][2]->Draw("PSAME");
      fKinematicsRD[i][2]->GetXaxis()->SetLabelSize(0.03);
      fKinematicsRD[i][2]->GetYaxis()->SetLabelSize(0.03);
      for(int tt=0; tt<fKinematicsMC[i][2]->GetNbinsX(); tt++)
      {
        fKinematicsMC[i][2]->SetBinError(tt,sqrt(fKinematicsMC[i][2]->GetBinContent(tt)));
      }
      fKinematicsMC[i][2]->Draw("E2SAME");
      fKinematicsMC[i][2]->Draw("SAME");
      c3.Update();

      c4.cd(i+1);
      fKinematicsMC[i][3]->Scale(fKinematicsRD[i][3]->GetEntries()/fKinematicsMC[i][3]->GetEntries());
      fKinematicsRD[i][3]->SetLineColor(kRed);
      fKinematicsMC[i][3]->SetLineColor(kBlue);
      fKinematicsRD[i][3]->SetMinimum(0.);
      fKinematicsRD[i][3]->SetMaximum(max(fKinematicsRD[i][3]->GetMaximum()*1.2,fKinematicsMC[i][3]->GetMaximum()*1.2));
      fKinematicsRD[i][3]->GetYaxis()->SetNdivisions(304,kTRUE);
      for(int tt=0; tt<fKinematicsRD[i][3]->GetNbinsX(); tt++)
      {
        fKinematicsRD[i][3]->SetBinError(tt,sqrt(fKinematicsRD[i][3]->GetBinContent(tt)));
      }
      fKinematicsRD[i][3]->Draw("E2");
      fKinematicsRD[i][3]->SetMarkerStyle(22);
      fKinematicsRD[i][3]->Draw("PSAME");
      fKinematicsRD[i][3]->GetXaxis()->SetLabelSize(0.03);
      fKinematicsRD[i][3]->GetYaxis()->SetLabelSize(0.03);
      for(int tt=0; tt<fKinematicsMC[i][3]->GetNbinsX(); tt++)
      {
        fKinematicsMC[i][3]->SetBinError(tt,sqrt(fKinematicsMC[i][3]->GetBinContent(tt)));
      }
      fKinematicsMC[i][3]->Draw("E2SAME");
      fKinematicsMC[i][3]->Draw("SAME");
      c4.Update();

      c5.cd(i+1);
      fKinematicsMC[i][4]->Scale(fKinematicsRD[i][4]->GetEntries()/fKinematicsMC[i][4]->GetEntries());
      fKinematicsRD[i][4]->SetLineColor(kRed);
      fKinematicsMC[i][4]->SetLineColor(kBlue);
      fKinematicsRD[i][4]->SetMinimum(0.);
      fKinematicsRD[i][4]->SetMaximum(max(fKinematicsRD[i][4]->GetMaximum()*1.2,fKinematicsMC[i][4]->GetMaximum()*1.2));
      fKinematicsRD[i][4]->GetYaxis()->SetNdivisions(304,kTRUE);
      for(int tt=0; tt<fKinematicsRD[i][4]->GetNbinsX(); tt++)
      {
        fKinematicsRD[i][4]->SetBinError(tt,sqrt(fKinematicsRD[i][4]->GetBinContent(tt)));
      }
      fKinematicsRD[i][4]->Draw("E2");
      fKinematicsRD[i][4]->SetMarkerStyle(22);
      fKinematicsRD[i][4]->Draw("PSAME");
      fKinematicsRD[i][4]->GetXaxis()->SetLabelSize(0.03);
      fKinematicsRD[i][4]->GetYaxis()->SetLabelSize(0.03);
      for(int tt=0; tt<fKinematicsMC[i][4]->GetNbinsX(); tt++)
      {
        fKinematicsMC[i][4]->SetBinError(tt,sqrt(fKinematicsMC[i][4]->GetBinContent(tt)));
      }
      fKinematicsMC[i][4]->Draw("E2SAME");
      fKinematicsMC[i][4]->Draw("SAME");
      c5.Update();

      c6.cd(i+1);
      fKinematicsMC[i][5]->Scale(fKinematicsRD[i][5]->GetEntries()/fKinematicsMC[i][5]->GetEntries());
      fKinematicsRD[i][5]->SetLineColor(kRed);
      fKinematicsMC[i][5]->SetLineColor(kBlue);
      fKinematicsRD[i][5]->SetMinimum(0.);
      fKinematicsRD[i][5]->SetMaximum(max(fKinematicsRD[i][5]->GetMaximum()*1.2,fKinematicsMC[i][5]->GetMaximum()*1.2));
      fKinematicsRD[i][5]->GetYaxis()->SetNdivisions(304,kTRUE);
      for(int tt=0; tt<fKinematicsRD[i][5]->GetNbinsX(); tt++)
      {
        fKinematicsRD[i][5]->SetBinError(tt,sqrt(fKinematicsRD[i][5]->GetBinContent(tt)));
      }
      fKinematicsRD[i][5]->Draw("E2");
      fKinematicsRD[i][5]->SetMarkerStyle(22);
      fKinematicsRD[i][5]->Draw("PSAME");
      fKinematicsRD[i][5]->GetXaxis()->SetLabelSize(0.03);
      fKinematicsRD[i][5]->GetYaxis()->SetLabelSize(0.03);
      for(int tt=0; tt<fKinematicsMC[i][5]->GetNbinsX(); tt++)
      {
        fKinematicsMC[i][5]->SetBinError(tt,sqrt(fKinematicsMC[i][5]->GetBinContent(tt)));
      }
      fKinematicsMC[i][5]->Draw("E2SAME");
      fKinematicsMC[i][5]->Draw("SAME");
      c6.Update();

      c14.cd(i+1);
      fKinematicsMC[i][6]->Scale(fKinematicsRD[i][6]->GetEntries()/fKinematicsMC[i][6]->GetEntries());
      fKinematicsRD[i][6]->SetLineColor(kRed);
      fKinematicsMC[i][6]->SetLineColor(kBlue);
      fKinematicsRD[i][6]->SetMinimum(0.);
      fKinematicsRD[i][6]->SetMaximum(max(fKinematicsRD[i][6]->GetMaximum()*1.2,fKinematicsMC[i][6]->GetMaximum()*1.2));
      fKinematicsRD[i][6]->GetYaxis()->SetNdivisions(304,kTRUE);
      for(int tt=0; tt<fKinematicsRD[i][6]->GetNbinsX(); tt++)
      {
        fKinematicsRD[i][6]->SetBinError(tt,sqrt(fKinematicsRD[i][6]->GetBinContent(tt)));
      }
      fKinematicsRD[i][6]->Draw("E2");
      fKinematicsRD[i][6]->SetMarkerStyle(22);
      fKinematicsRD[i][6]->Draw("PSAME");
      fKinematicsRD[i][6]->GetXaxis()->SetLabelSize(0.03);
      fKinematicsRD[i][6]->GetYaxis()->SetLabelSize(0.03);
      for(int tt=0; tt<fKinematicsMC[i][6]->GetNbinsX(); tt++)
      {
        fKinematicsMC[i][6]->SetBinError(tt,sqrt(fKinematicsMC[i][6]->GetBinContent(tt)));
      }
      fKinematicsMC[i][6]->Draw("E2SAME");
      fKinematicsMC[i][6]->Draw("SAME");
      c14.Update();

      c15.cd(i+1);
      fKinematicsMC[i][7]->Scale(fKinematicsRD[i][7]->GetEntries()/fKinematicsMC[i][7]->GetEntries());
      fKinematicsRD[i][7]->SetLineColor(kRed);
      fKinematicsMC[i][7]->SetLineColor(kBlue);
      fKinematicsRD[i][7]->SetMinimum(0.);
      fKinematicsRD[i][7]->SetMaximum(max(fKinematicsRD[i][7]->GetMaximum()*1.2,fKinematicsMC[i][7]->GetMaximum()*1.2));
      fKinematicsRD[i][7]->GetYaxis()->SetNdivisions(304,kTRUE);
      for(int tt=0; tt<fKinematicsRD[i][7]->GetNbinsX(); tt++)
      {
        fKinematicsRD[i][7]->SetBinError(tt,sqrt(fKinematicsRD[i][7]->GetBinContent(tt)));
      }
      fKinematicsRD[i][7]->Draw("E2");
      fKinematicsRD[i][7]->SetMarkerStyle(22);
      fKinematicsRD[i][7]->Draw("PSAME");
      fKinematicsRD[i][7]->GetXaxis()->SetLabelSize(0.03);
      fKinematicsRD[i][7]->GetYaxis()->SetLabelSize(0.03);
      for(int tt=0; tt<fKinematicsMC[i][7]->GetNbinsX(); tt++)
      {
        fKinematicsMC[i][7]->SetBinError(tt,sqrt(fKinematicsMC[i][7]->GetBinContent(tt)));
      }
      fKinematicsMC[i][7]->Draw("E2SAME");
      fKinematicsMC[i][7]->Draw("SAME");
      c15.Update();

      c16.cd(i+1);
      fKinematicsMC[i][8]->Scale(fKinematicsRD[i][8]->GetEntries()/fKinematicsMC[i][8]->GetEntries());
      fKinematicsRD[i][8]->SetLineColor(kRed);
      fKinematicsMC[i][8]->SetLineColor(kBlue);
      fKinematicsRD[i][8]->SetMinimum(0.);
      fKinematicsRD[i][8]->SetMaximum(max(fKinematicsRD[i][8]->GetMaximum()*1.2,fKinematicsMC[i][8]->GetMaximum()*1.2));
      fKinematicsRD[i][8]->GetYaxis()->SetNdivisions(304,kTRUE);
      for(int tt=0; tt<fKinematicsRD[i][8]->GetNbinsX(); tt++)
      {
        fKinematicsRD[i][8]->SetBinError(tt,sqrt(fKinematicsRD[i][8]->GetBinContent(tt)));
      }
      fKinematicsRD[i][8]->Draw("E2");
      fKinematicsRD[i][8]->SetMarkerStyle(22);
      fKinematicsRD[i][8]->Draw("PSAME");
      fKinematicsRD[i][8]->GetXaxis()->SetLabelSize(0.03);
      fKinematicsRD[i][8]->GetYaxis()->SetLabelSize(0.03);
      for(int tt=0; tt<fKinematicsMC[i][8]->GetNbinsX(); tt++)
      {
        fKinematicsMC[i][8]->SetBinError(tt,sqrt(fKinematicsMC[i][8]->GetBinContent(tt)));
      }
      fKinematicsMC[i][8]->Draw("E2SAME");
      fKinematicsMC[i][8]->Draw("SAME");
      c16.Update();

      c17.cd(i+1);
      fKinematicsMC[i][9]->Scale(fKinematicsRD[i][9]->GetEntries()/fKinematicsMC[i][9]->GetEntries());
      fKinematicsRD[i][9]->SetLineColor(kRed);
      fKinematicsMC[i][9]->SetLineColor(kBlue);
      fKinematicsRD[i][9]->SetMinimum(0.);
      fKinematicsRD[i][9]->SetMaximum(max(fKinematicsRD[i][9]->GetMaximum()*1.2,fKinematicsMC[i][9]->GetMaximum()*1.2));
      fKinematicsRD[i][9]->GetYaxis()->SetNdivisions(304,kTRUE);
      for(int tt=0; tt<fKinematicsRD[i][9]->GetNbinsX(); tt++)
      {
        fKinematicsRD[i][9]->SetBinError(tt,sqrt(fKinematicsRD[i][9]->GetBinContent(tt)));
      }
      fKinematicsRD[i][9]->Draw("E2");
      fKinematicsRD[i][9]->SetMarkerStyle(22);
      fKinematicsRD[i][9]->Draw("PSAME");
      fKinematicsRD[i][9]->GetXaxis()->SetLabelSize(0.03);
      fKinematicsRD[i][9]->GetYaxis()->SetLabelSize(0.03);
      for(int tt=0; tt<fKinematicsMC[i][9]->GetNbinsX(); tt++)
      {
        fKinematicsMC[i][9]->SetBinError(tt,sqrt(fKinematicsMC[i][9]->GetBinContent(tt)));
      }
      fKinematicsMC[i][9]->Draw("E2SAME");
      fKinematicsMC[i][9]->Draw("SAME");
      c17.Update();

      c18.cd(i+1);
      fKinematicsMC[i][10]->Scale(fKinematicsRD[i][10]->GetEntries()/fKinematicsMC[i][10]->GetEntries());
      fKinematicsRD[i][10]->SetLineColor(kRed);
      fKinematicsMC[i][10]->SetLineColor(kBlue);
      fKinematicsRD[i][10]->SetMinimum(0.);
      fKinematicsRD[i][10]->SetMaximum(max(fKinematicsRD[i][10]->GetMaximum()*1.2,fKinematicsMC[i][10]->GetMaximum()*1.2));
      fKinematicsRD[i][10]->GetYaxis()->SetNdivisions(304,kTRUE);
      for(int tt=0; tt<fKinematicsRD[i][10]->GetNbinsX(); tt++)
      {
        fKinematicsRD[i][10]->SetBinError(tt,sqrt(fKinematicsRD[i][10]->GetBinContent(tt)));
      }
      fKinematicsRD[i][10]->Draw("E2");
      fKinematicsRD[i][10]->SetMarkerStyle(22);
      fKinematicsRD[i][10]->Draw("PSAME");
      fKinematicsRD[i][10]->GetXaxis()->SetLabelSize(0.03);
      fKinematicsRD[i][10]->GetYaxis()->SetLabelSize(0.03);
      for(int tt=0; tt<fKinematicsMC[i][10]->GetNbinsX(); tt++)
      {
        fKinematicsMC[i][10]->SetBinError(tt,sqrt(fKinematicsMC[i][10]->GetBinContent(tt)));
      }
      fKinematicsMC[i][10]->Draw("E2SAME");
      fKinematicsMC[i][10]->Draw("SAME");
      c18.Update();

      c19.cd(i+1);
      fKinematicsMC[i][12]->Scale(fKinematicsRD[i][12]->GetEntries()/fKinematicsMC[i][12]->GetEntries());
      fKinematicsRD[i][12]->SetLineColor(kRed);
      fKinematicsMC[i][12]->SetLineColor(kBlue);
      fKinematicsRD[i][12]->SetMinimum(0.);
      fKinematicsRD[i][12]->SetMaximum(max(fKinematicsRD[i][12]->GetMaximum()*1.2,fKinematicsMC[i][12]->GetMaximum()*1.2));
      fKinematicsRD[i][12]->GetYaxis()->SetNdivisions(304,kTRUE);
      for(int tt=0; tt<fKinematicsRD[i][12]->GetNbinsX(); tt++)
      {
        fKinematicsRD[i][12]->SetBinError(tt,sqrt(fKinematicsRD[i][12]->GetBinContent(tt)));
      }
      fKinematicsRD[i][12]->Draw("E2");
      fKinematicsRD[i][12]->SetMarkerStyle(22);
      fKinematicsRD[i][12]->Draw("PSAME");
      fKinematicsRD[i][12]->GetXaxis()->SetLabelSize(0.03);
      fKinematicsRD[i][12]->GetYaxis()->SetLabelSize(0.03);
      for(int tt=0; tt<fKinematicsMC[i][12]->GetNbinsX(); tt++)
      {
        fKinematicsMC[i][12]->SetBinError(tt,sqrt(fKinematicsMC[i][12]->GetBinContent(tt)));
      }
      fKinematicsMC[i][12]->Draw("E2SAME");
      fKinematicsMC[i][12]->Draw("SAME");
      c19.Update();

      c20.cd(i+1);
      fKinematicsMC[i][13]->Scale(fKinematicsRD[i][13]->GetEntries()/fKinematicsMC[i][13]->GetEntries());
      fKinematicsRD[i][13]->SetLineColor(kRed);
      fKinematicsMC[i][13]->SetLineColor(kBlue);
      fKinematicsRD[i][13]->SetMinimum(0.);
      fKinematicsRD[i][13]->SetMaximum(max(fKinematicsRD[i][13]->GetMaximum()*1.2,fKinematicsMC[i][13]->GetMaximum()*1.2));
      fKinematicsRD[i][13]->GetYaxis()->SetNdivisions(304,kTRUE);
      for(int tt=0; tt<fKinematicsRD[i][13]->GetNbinsX(); tt++)
      {
        fKinematicsRD[i][13]->SetBinError(tt,sqrt(fKinematicsRD[i][13]->GetBinContent(tt)));
      }
      fKinematicsRD[i][13]->Draw("E2");
      fKinematicsRD[i][13]->SetMarkerStyle(22);
      fKinematicsRD[i][13]->Draw("PSAME");
      fKinematicsRD[i][13]->GetXaxis()->SetLabelSize(0.03);
      fKinematicsRD[i][13]->GetYaxis()->SetLabelSize(0.03);
      for(int tt=0; tt<fKinematicsMC[i][13]->GetNbinsX(); tt++)
      {
        fKinematicsMC[i][13]->SetBinError(tt,sqrt(fKinematicsMC[i][13]->GetBinContent(tt)));
      }
      fKinematicsMC[i][13]->Draw("E2SAME");
      fKinematicsMC[i][13]->Draw("SAME");
      c20.Update();

      c21.cd(i+1);
      fKinematicsMC[i][14]->Scale(fKinematicsRD[i][14]->GetEntries()/fKinematicsMC[i][14]->GetEntries());
      fKinematicsRD[i][14]->SetLineColor(kRed);
      fKinematicsMC[i][14]->SetLineColor(kBlue);
      fKinematicsRD[i][14]->SetMinimum(0.);
      fKinematicsRD[i][14]->SetMaximum(max(fKinematicsRD[i][14]->GetMaximum()*1.2,fKinematicsMC[i][14]->GetMaximum()*1.2));
      fKinematicsRD[i][14]->GetYaxis()->SetNdivisions(304,kTRUE);
      for(int tt=0; tt<fKinematicsRD[i][14]->GetNbinsX(); tt++)
      {
        fKinematicsRD[i][14]->SetBinError(tt,sqrt(fKinematicsRD[i][14]->GetBinContent(tt)));
      }
      fKinematicsRD[i][14]->Draw("E2");
      fKinematicsRD[i][14]->SetMarkerStyle(22);
      fKinematicsRD[i][14]->Draw("PSAME");
      fKinematicsRD[i][14]->GetXaxis()->SetLabelSize(0.03);
      fKinematicsRD[i][14]->GetYaxis()->SetLabelSize(0.03);
      for(int tt=0; tt<fKinematicsMC[i][14]->GetNbinsX(); tt++)
      {
        fKinematicsMC[i][14]->SetBinError(tt,sqrt(fKinematicsMC[i][14]->GetBinContent(tt)));
      }
      fKinematicsMC[i][14]->Draw("E2SAME");
      fKinematicsMC[i][14]->Draw("SAME");
      c21.Update();

      c25.cd(i+1);
      fKinematicsMC[i][15]->Scale(fKinematicsRD[i][15]->GetEntries()/fKinematicsMC[i][15]->GetEntries());
      fKinematicsRD[i][15]->SetLineColor(kRed);
      fKinematicsMC[i][15]->SetLineColor(kBlue);
      fKinematicsRD[i][15]->SetMinimum(0.);
      fKinematicsRD[i][15]->SetMaximum(max(fKinematicsRD[i][15]->GetMaximum()*1.2,fKinematicsMC[i][15]->GetMaximum()*1.2));
      fKinematicsRD[i][15]->GetYaxis()->SetNdivisions(304,kTRUE);
      for(int tt=0; tt<fKinematicsRD[i][15]->GetNbinsX(); tt++)
      {
        fKinematicsRD[i][15]->SetBinError(tt,sqrt(fKinematicsRD[i][15]->GetBinContent(tt)));
      }
      fKinematicsRD[i][15]->Draw("E2");
      fKinematicsRD[i][15]->SetMarkerStyle(22);
      fKinematicsRD[i][15]->Draw("PSAME");
      fKinematicsRD[i][15]->GetXaxis()->SetLabelSize(0.03);
      fKinematicsRD[i][15]->GetYaxis()->SetLabelSize(0.03);
      for(int tt=0; tt<fKinematicsMC[i][15]->GetNbinsX(); tt++)
      {
        fKinematicsMC[i][15]->SetBinError(tt,sqrt(fKinematicsMC[i][15]->GetBinContent(tt)));
      }
      fKinematicsMC[i][15]->Draw("E2SAME");
      fKinematicsMC[i][15]->Draw("SAME");
      c25.Update();

      c27.cd(i+1);
      fKinematicsMC[i][16]->Scale(fKinematicsRD[i][16]->GetEntries()/fKinematicsMC[i][16]->GetEntries());
      fKinematicsRD[i][16]->SetLineColor(kRed);
      fKinematicsMC[i][16]->SetLineColor(kBlue);
      fKinematicsRD[i][16]->SetMinimum(0.);
      fKinematicsRD[i][16]->SetMaximum(max(fKinematicsRD[i][16]->GetMaximum()*1.2,fKinematicsMC[i][16]->GetMaximum()*1.2));
      fKinematicsRD[i][16]->GetYaxis()->SetNdivisions(304,kTRUE);
      for(int tt=0; tt<fKinematicsRD[i][16]->GetNbinsX(); tt++)
      {
        fKinematicsRD[i][16]->SetBinError(tt,sqrt(fKinematicsRD[i][16]->GetBinContent(tt)));
      }
      fKinematicsRD[i][16]->Draw("E2");
      fKinematicsRD[i][16]->SetMarkerStyle(22);
      fKinematicsRD[i][16]->Draw("PSAME");
      fKinematicsRD[i][16]->GetXaxis()->SetLabelSize(0.03);
      fKinematicsRD[i][16]->GetYaxis()->SetLabelSize(0.03);
      for(int tt=0; tt<fKinematicsMC[i][16]->GetNbinsX(); tt++)
      {
        fKinematicsMC[i][16]->SetBinError(tt,sqrt(fKinematicsMC[i][16]->GetBinContent(tt)));
      }
      fKinematicsMC[i][16]->Draw("E2SAME");
      fKinematicsMC[i][16]->Draw("SAME");
      c27.Update();
    // }
  }

  // c7.cd(2);
  // for(int tt=0; tt<fKinematicsRD[0][11]->GetNbinsX(); tt++)
  // {
  //   fError.push_back((fKinematicsRD[0][11]->GetBinError(tt) && fKinematicsMC[0][11]->GetBinError(tt) ? sqrt(pow(1/fKinematicsRD[0][11]->GetBinError(tt),2)+pow(1/fKinematicsMC[0][11]->GetBinError(tt),2)):0));
  // }
  // fKinematicsRD[0][11]->Scale(1/fKinematicsRD[0][11]->GetEntries());
  // fKinematicsMC[0][11]->Scale(1/fKinematicsMC[0][11]->GetEntries());
  // fKinematicsRatio[0][11] = (TH1F*)fKinematicsRD[0][11]->Clone();
  // fKinematicsRatio[0][11]->SetStats(0);
  // fKinematicsRatio[0][11]->Divide(fKinematicsMC[0][11]);
  // for(int tt=0; tt<fKinematicsRatio[0][11]->GetNbinsX(); tt++)
  // {
  //   fKinematicsRatio[0][11]->SetBinError(tt,fError[tt]);
  // }
  // fError.clear();
  // fKinematicsRatio[0][11]->SetMarkerStyle(21);
  // fKinematicsRatio[0][11]->SetFillColor(kYellow-7);
  // fKinematicsRatio[0][11]->SetMaximum(2.);
  // fKinematicsRatio[0][11]->SetMinimum(0.);
  // fKinematicsRatio[0][11]->Draw("PE2");
  // fKinematicsRatio[0][11]->GetXaxis()->SetLabelSize(0.08);
  // fKinematicsRatio[0][11]->GetYaxis()->SetLabelSize(0.08);
  // fKinematicsRatio[0][11]->GetYaxis()->SetNdivisions(2,kTRUE);
  // for(int tt=0; tt<7; tt++)
  // {
  //   l1[11][tt]->Draw();
  // }
  // c7.Update();
  c7.cd(1);
  fKinematicsMC[0][11]->Scale(fKinematicsRD[0][11]->GetEntries()/fKinematicsMC[0][11]->GetEntries());
  fKinematicsRD[0][11]->SetLineColor(kRed);
  fKinematicsMC[0][11]->SetLineColor(kBlue);
  fKinematicsRD[0][11]->SetMinimum(0.);
  fKinematicsRD[0][11]->SetMaximum(max(fKinematicsRD[0][11]->GetMaximum()*1.2,fKinematicsMC[0][11]->GetMaximum()*1.2));
  fKinematicsRD[0][11]->GetYaxis()->SetNdivisions(304,kTRUE);
  for(int tt=0; tt<fKinematicsRD[0][11]->GetNbinsX(); tt++)
  {
    fKinematicsRD[0][11]->SetBinError(tt,sqrt(fKinematicsRD[0][11]->GetBinContent(tt)));
  }
  fKinematicsRD[0][11]->Draw("E2");
  fKinematicsRD[0][11]->SetMarkerStyle(22);
  fKinematicsRD[0][11]->Draw("PSAME");
  fKinematicsRD[0][11]->GetXaxis()->SetLabelSize(0.03);
  fKinematicsRD[0][11]->GetYaxis()->SetLabelSize(0.03);
  for(int tt=0; tt<fKinematicsMC[0][11]->GetNbinsX(); tt++)
  {
    fKinematicsMC[0][11]->SetBinError(tt,sqrt(fKinematicsMC[0][11]->GetBinContent(tt)));
  }
  fKinematicsMC[0][11]->Draw("E2SAME");
  fKinematicsMC[0][11]->Draw("SAME");
  c7.Update();

  // c8.cd(2);
  // for(int tt=0; tt<fKinematicsRD[4][0]->GetNbinsX(); tt++)
  // {
  //   fError.push_back((fKinematicsRD[4][0]->GetBinError(tt) && fKinematicsMC[4][0]->GetBinError(tt) ? sqrt(pow(1/fKinematicsRD[4][0]->GetBinError(tt),2)+pow(1/fKinematicsMC[4][0]->GetBinError(tt),2)):0));
  // }
  // fKinematicsRD[4][0]->Scale(1/fKinematicsRD[4][0]->GetEntries());
  // fKinematicsMC[4][0]->Scale(1/fKinematicsMC[4][0]->GetEntries());
  // fKinematicsRatio[4][0] = (TH1F*)fKinematicsRD[4][0]->Clone();
  // fKinematicsRatio[4][0]->SetStats(0);
  // fKinematicsRatio[4][0]->Divide(fKinematicsMC[4][0]);
  // for(int tt=0; tt<fKinematicsRatio[4][0]->GetNbinsX(); tt++)
  // {
  //   fKinematicsRatio[4][0]->SetBinError(tt,fError[tt]);
  // }
  // fError.clear();
  // fKinematicsRatio[4][0]->SetMarkerStyle(21);
  // fKinematicsRatio[4][0]->SetFillColor(kYellow-7);
  // fKinematicsRatio[4][0]->SetMaximum(2.);
  // fKinematicsRatio[4][0]->SetMinimum(0.);
  // fKinematicsRatio[4][0]->Draw("PE2");
  // fKinematicsRatio[4][0]->GetXaxis()->SetLabelSize(0.08);
  // fKinematicsRatio[4][0]->GetYaxis()->SetLabelSize(0.08);
  // fKinematicsRatio[4][0]->GetYaxis()->SetNdivisions(2,kFALSE);
  // for(int tt=0; tt<7; tt++)
  // {
  //   l1[0][tt]->Draw();
  // }
  // gPad->SetLogx();
  // c8.Update();
  c8.cd(1);
  fKinematicsMC[4][0]->Scale(fKinematicsRD[4][0]->GetEntries()/fKinematicsMC[4][0]->GetEntries());
  fKinematicsRD[4][0]->SetLineColor(kRed);
  fKinematicsMC[4][0]->SetLineColor(kBlue);
  fKinematicsRD[4][0]->SetMinimum(0.);
  fKinematicsRD[4][0]->SetMaximum(max(fKinematicsRD[4][0]->GetMaximum()*1.2,fKinematicsMC[4][0]->GetMaximum()*1.2));
  fKinematicsRD[4][0]->GetYaxis()->SetNdivisions(304,kTRUE);
  for(int tt=0; tt<fKinematicsRD[4][0]->GetNbinsX(); tt++)
  {
    fKinematicsRD[4][0]->SetBinError(tt,sqrt(fKinematicsRD[4][0]->GetBinContent(tt)));
  }
  fKinematicsRD[4][0]->Draw("E2");
  fKinematicsRD[4][0]->SetMarkerStyle(22);
  fKinematicsRD[4][0]->Draw("PSAME");
  fKinematicsRD[4][0]->GetXaxis()->SetLabelSize(0.03);
  fKinematicsRD[4][0]->GetYaxis()->SetLabelSize(0.03);
  for(int tt=0; tt<fKinematicsMC[4][0]->GetNbinsX(); tt++)
  {
    fKinematicsMC[4][0]->SetBinError(tt,sqrt(fKinematicsMC[4][0]->GetBinContent(tt)));
  }
  fKinematicsMC[4][0]->Draw("E2SAME");
  fKinematicsMC[4][0]->Draw("SAME");
  gPad->SetLogx();
  c8.Update();

  // c30.cd(1);
  // for(int tt=0; tt<fKinematicsRD[4][0]->GetNbinsX(); tt++)
  // {
  //   fKinematicsRD[4][0]->SetBinError(tt,sqrt(fKinematicsRD[4][0]->GetBinContent(tt)));
  // }
  // fKinematicsRD[4][0]->Draw("E2");
  // fKinematicsRD[4][0]->Draw("SAME");
  // fKinematicsRD[4][0]->GetXaxis()->SetLabelSize(0.03);
  // fKinematicsRD[4][0]->GetYaxis()->SetLabelSize(0.03);
  // for(int tt=0; tt<fKinematicsMC[4][0]->GetNbinsX(); tt++)
  // {
  //   fKinematicsMC[4][0]->SetBinError(tt,sqrt(fKinematicsMC[4][0]->GetBinContent(tt)));
  // }
  // fKinematicsMC[4][0]->Draw("E2SAME");
  // fKinematicsMC[4][0]->Draw("SAME");
  // gPad->SetLogx();
  // c30.Update();

  // c9.cd(2);
  // for(int tt=0; tt<fKinematicsRD[4][1]->GetNbinsX(); tt++)
  // {
  //   fError.push_back((fKinematicsRD[4][1]->GetBinError(tt) && fKinematicsMC[4][1]->GetBinError(tt) ? sqrt(pow(1/fKinematicsRD[4][1]->GetBinError(tt),2)+pow(1/fKinematicsMC[4][1]->GetBinError(tt),2)):0));
  // }
  // fKinematicsRD[4][1]->Scale(1/fKinematicsRD[4][1]->GetEntries());
  // fKinematicsMC[4][1]->Scale(1/fKinematicsMC[4][1]->GetEntries());
  // fKinematicsRatio[4][1] = (TH1F*)fKinematicsRD[4][1]->Clone();
  // fKinematicsRatio[4][1]->SetStats(0);
  // fKinematicsRatio[4][1]->Divide(fKinematicsMC[4][1]);
  // for(int tt=0; tt<fKinematicsRatio[4][1]->GetNbinsX(); tt++)
  // {
  //   fKinematicsRatio[4][1]->SetBinError(tt,fError[tt]);
  // }
  // fError.clear();
  // fKinematicsRatio[4][1]->SetMarkerStyle(21);
  // fKinematicsRatio[4][1]->SetFillColor(kYellow-7);
  // fKinematicsRatio[4][1]->SetMaximum(2.);
  // fKinematicsRatio[4][1]->SetMinimum(0.);
  // fKinematicsRatio[4][1]->Draw("PE2");
  // fKinematicsRatio[4][1]->GetXaxis()->SetLabelSize(0.08);
  // fKinematicsRatio[4][1]->GetYaxis()->SetLabelSize(0.08);
  // fKinematicsRatio[4][1]->GetYaxis()->SetNdivisions(2,kFALSE);
  // for(int tt=0; tt<7; tt++)
  // {
  //   l1[1][tt]->Draw();
  // }
  // gPad->SetLogx();
  // c9.Update();
  c9.cd(1);
  fKinematicsMC[4][1]->Scale(fKinematicsRD[4][1]->GetEntries()/fKinematicsMC[4][1]->GetEntries());
  fKinematicsRD[4][1]->SetLineColor(kRed);
  fKinematicsMC[4][1]->SetLineColor(kBlue);
  fKinematicsRD[4][1]->SetMinimum(0.);
  fKinematicsRD[4][1]->SetMaximum(max(fKinematicsRD[4][1]->GetMaximum()*1.2,fKinematicsMC[4][1]->GetMaximum()*1.2));
  fKinematicsRD[4][1]->GetYaxis()->SetNdivisions(304,kTRUE);
  for(int tt=0; tt<fKinematicsRD[4][1]->GetNbinsX(); tt++)
  {
    fKinematicsRD[4][1]->SetBinError(tt,sqrt(fKinematicsRD[4][1]->GetBinContent(tt)));
  }
  fKinematicsRD[4][1]->Draw("E2");
  fKinematicsRD[4][1]->SetMarkerStyle(22);
  fKinematicsRD[4][1]->Draw("PSAME");
  fKinematicsRD[4][1]->GetXaxis()->SetLabelSize(0.03);
  fKinematicsRD[4][1]->GetYaxis()->SetLabelSize(0.03);
  for(int tt=0; tt<fKinematicsMC[4][1]->GetNbinsX(); tt++)
  {
    fKinematicsMC[4][1]->SetBinError(tt,sqrt(fKinematicsMC[4][1]->GetBinContent(tt)));
  }
  fKinematicsMC[4][1]->Draw("E2SAME");
  fKinematicsMC[4][1]->Draw("SAME");
  gPad->SetLogx();
  c9.Update();

  // c10.cd(2);
  // for(int tt=0; tt<fKinematicsRD[4][2]->GetNbinsX(); tt++)
  // {
  //   fError.push_back((fKinematicsRD[4][2]->GetBinError(tt) && fKinematicsMC[4][2]->GetBinError(tt) ? sqrt(pow(1/fKinematicsRD[4][2]->GetBinError(tt),2)+pow(1/fKinematicsMC[4][2]->GetBinError(tt),2)):0));
  // }
  // fKinematicsRD[4][2]->Scale(1/fKinematicsRD[4][2]->GetEntries());
  // fKinematicsMC[4][2]->Scale(1/fKinematicsMC[4][2]->GetEntries());
  // fKinematicsRatio[4][2] = (TH1F*)fKinematicsRD[4][2]->Clone();
  // fKinematicsRatio[4][2]->SetStats(0);
  // fKinematicsRatio[4][2]->Divide(fKinematicsMC[4][2]);
  // for(int tt=0; tt<fKinematicsRatio[4][2]->GetNbinsX(); tt++)
  // {
  //   fKinematicsRatio[4][2]->SetBinError(tt,fError[tt]);
  // }
  // fError.clear();
  // fKinematicsRatio[4][2]->SetMarkerStyle(21);
  // fKinematicsRatio[4][2]->SetFillColor(kYellow-7);
  // fKinematicsRatio[4][2]->SetMaximum(2.);
  // fKinematicsRatio[4][2]->SetMinimum(0.);
  // fKinematicsRatio[4][2]->Draw("PE2");
  // fKinematicsRatio[4][2]->GetXaxis()->SetLabelSize(0.08);
  // fKinematicsRatio[4][2]->GetYaxis()->SetLabelSize(0.08);
  // fKinematicsRatio[4][2]->GetYaxis()->SetNdivisions(2,kFALSE);
  // for(int tt=0; tt<7; tt++)
  // {
  //   l1[2][tt]->Draw();
  // }
  // c10.Update();
  c10.cd(1);
  fKinematicsMC[4][2]->Scale(fKinematicsRD[4][2]->GetEntries()/fKinematicsMC[4][2]->GetEntries());
  fKinematicsRD[4][2]->SetLineColor(kRed);
  fKinematicsMC[4][2]->SetLineColor(kBlue);
  fKinematicsRD[4][2]->SetMinimum(0.);
  fKinematicsRD[4][2]->SetMaximum(max(fKinematicsRD[4][2]->GetMaximum()*1.2,fKinematicsMC[4][2]->GetMaximum()*1.2));
  fKinematicsRD[4][2]->GetYaxis()->SetNdivisions(304,kTRUE);
  for(int tt=0; tt<fKinematicsRD[4][2]->GetNbinsX(); tt++)
  {
    fKinematicsRD[4][2]->SetBinError(tt,sqrt(fKinematicsRD[4][2]->GetBinContent(tt)));
  }
  fKinematicsRD[4][2]->Draw("E2");
  fKinematicsRD[4][2]->SetMarkerStyle(22);
  fKinematicsRD[4][2]->Draw("PSAME");
  fKinematicsRD[4][2]->GetXaxis()->SetLabelSize(0.03);
  fKinematicsRD[4][2]->GetYaxis()->SetLabelSize(0.03);
  for(int tt=0; tt<fKinematicsMC[4][2]->GetNbinsX(); tt++)
  {
    fKinematicsMC[4][2]->SetBinError(tt,sqrt(fKinematicsMC[4][2]->GetBinContent(tt)));
  }
  fKinematicsMC[4][2]->Draw("E2SAME");
  fKinematicsMC[4][2]->Draw("SAME");
  c10.Update();

  // c11.cd(2);
  // for(int tt=0; tt<fKinematicsRD[4][3]->GetNbinsX(); tt++)
  // {
  //   fError.push_back((fKinematicsRD[4][3]->GetBinError(tt) && fKinematicsMC[4][3]->GetBinError(tt) ? sqrt(pow(1/fKinematicsRD[4][3]->GetBinError(tt),2)+pow(1/fKinematicsMC[4][3]->GetBinError(tt),2)):0));
  // }
  // fKinematicsRD[4][3]->Scale(1/fKinematicsRD[4][3]->GetEntries());
  // fKinematicsMC[4][3]->Scale(1/fKinematicsMC[4][3]->GetEntries());
  // fKinematicsRatio[4][3] = (TH1F*)fKinematicsRD[4][3]->Clone();
  // fKinematicsRatio[4][3]->SetStats(0);
  // fKinematicsRatio[4][3]->Divide(fKinematicsMC[4][3]);
  // for(int tt=0; tt<fKinematicsRatio[4][3]->GetNbinsX(); tt++)
  // {
  //   fKinematicsRatio[4][3]->SetBinError(tt,fError[tt]);
  // }
  // fError.clear();
  // fKinematicsRatio[4][3]->SetMarkerStyle(21);
  // fKinematicsRatio[4][3]->SetFillColor(kYellow-7);
  // fKinematicsRatio[4][3]->SetMaximum(2.);
  // fKinematicsRatio[4][3]->SetMinimum(0.);
  // fKinematicsRatio[4][3]->Draw("PE2");
  // fKinematicsRatio[4][3]->GetXaxis()->SetLabelSize(0.08);
  // fKinematicsRatio[4][3]->GetYaxis()->SetLabelSize(0.08);
  // fKinematicsRatio[4][3]->GetYaxis()->SetNdivisions(2,kFALSE);
  // for(int tt=0; tt<7; tt++)
  // {
  //   l1[3][tt]->Draw();
  // }
  // c11.Update();
  c11.cd(1);
  fKinematicsMC[4][3]->Scale(fKinematicsRD[4][3]->GetEntries()/fKinematicsMC[4][3]->GetEntries());
  fKinematicsRD[4][3]->SetLineColor(kRed);
  fKinematicsMC[4][3]->SetLineColor(kBlue);
  fKinematicsRD[4][3]->SetMinimum(0.);
  fKinematicsRD[4][3]->SetMaximum(max(fKinematicsRD[4][3]->GetMaximum()*1.2,fKinematicsMC[4][3]->GetMaximum()*1.2));
  fKinematicsRD[4][3]->GetYaxis()->SetNdivisions(304,kTRUE);
  for(int tt=0; tt<fKinematicsRD[4][3]->GetNbinsX(); tt++)
  {
    fKinematicsRD[4][3]->SetBinError(tt,sqrt(fKinematicsRD[4][3]->GetBinContent(tt)));
  }
  fKinematicsRD[4][3]->Draw("E2");
  fKinematicsRD[4][3]->SetMarkerStyle(22);
  fKinematicsRD[4][3]->Draw("PSAME");
  fKinematicsRD[4][3]->GetXaxis()->SetLabelSize(0.03);
  fKinematicsRD[4][3]->GetYaxis()->SetLabelSize(0.03);
  for(int tt=0; tt<fKinematicsMC[4][3]->GetNbinsX(); tt++)
  {
    fKinematicsMC[4][3]->SetBinError(tt,sqrt(fKinematicsMC[4][3]->GetBinContent(tt)));
  }
  fKinematicsMC[4][3]->Draw("E2SAME");
  fKinematicsMC[4][3]->Draw("SAME");
  c11.Update();

  // c31.cd(1);
  // for(int tt=0; tt<fKinematicsRD[4][3]->GetNbinsX(); tt++)
  // {
  //   fKinematicsRD[4][3]->SetBinError(tt,sqrt(fKinematicsRD[4][3]->GetBinContent(tt)));
  // }
  // fKinematicsRD[4][3]->Draw("E2");
  // fKinematicsRD[4][3]->Draw("SAME");
  // fKinematicsRD[4][3]->GetXaxis()->SetLabelSize(0.03);
  // fKinematicsRD[4][3]->GetYaxis()->SetLabelSize(0.03);
  // for(int tt=0; tt<fKinematicsMC[4][3]->GetNbinsX(); tt++)
  // {
  //   fKinematicsMC[4][3]->SetBinError(tt,sqrt(fKinematicsMC[4][3]->GetBinContent(tt)));
  // }
  // fKinematicsMC[4][3]->Draw("E2SAME");
  // fKinematicsMC[4][3]->Draw("SAME");
  // c31.Update();

  // c12.cd(2);
  // for(int tt=0; tt<fKinematicsRD[4][4]->GetNbinsX(); tt++)
  // {
  //   fError.push_back((fKinematicsRD[4][4]->GetBinError(tt) && fKinematicsMC[4][4]->GetBinError(tt) ? sqrt(pow(1/fKinematicsRD[4][4]->GetBinError(tt),2)+pow(1/fKinematicsMC[4][4]->GetBinError(tt),2)):0));
  // }
  // fKinematicsRD[4][4]->Scale(1/fKinematicsRD[4][4]->GetEntries());
  // fKinematicsMC[4][4]->Scale(1/fKinematicsMC[4][4]->GetEntries());
  // fKinematicsRatio[4][4] = (TH1F*)fKinematicsRD[4][4]->Clone();
  // fKinematicsRatio[4][4]->SetStats(0);
  // fKinematicsRatio[4][4]->Divide(fKinematicsMC[4][4]);
  // for(int tt=0; tt<fKinematicsRatio[4][4]->GetNbinsX(); tt++)
  // {
  //   fKinematicsRatio[4][4]->SetBinError(tt,fError[tt]);
  // }
  // fError.clear();
  // fKinematicsRatio[4][4]->SetMarkerStyle(21);
  // fKinematicsRatio[4][4]->SetFillColor(kYellow-7);
  // fKinematicsRatio[4][4]->SetMaximum(2.);
  // fKinematicsRatio[4][4]->SetMinimum(0.);
  // fKinematicsRatio[4][4]->Draw("PE2");
  // fKinematicsRatio[4][4]->GetXaxis()->SetLabelSize(0.08);
  // fKinematicsRatio[4][4]->GetYaxis()->SetLabelSize(0.08);
  // fKinematicsRatio[4][4]->GetYaxis()->SetNdivisions(2,kFALSE);
  // for(int tt=0; tt<7; tt++)
  // {
  //   l1[4][tt]->Draw();
  // }
  // c12.Update();
  c12.cd(1);
  fKinematicsMC[4][4]->Scale(fKinematicsRD[4][4]->GetEntries()/fKinematicsMC[4][4]->GetEntries());
  fKinematicsRD[4][4]->SetLineColor(kRed);
  fKinematicsMC[4][4]->SetLineColor(kBlue);
  fKinematicsRD[4][4]->SetMinimum(0.);
  fKinematicsRD[4][4]->SetMaximum(max(fKinematicsRD[4][4]->GetMaximum()*1.2,fKinematicsMC[4][4]->GetMaximum()*1.2));
  fKinematicsRD[4][4]->GetYaxis()->SetNdivisions(304,kTRUE);
  for(int tt=0; tt<fKinematicsRD[4][4]->GetNbinsX(); tt++)
  {
    fKinematicsRD[4][4]->SetBinError(tt,sqrt(fKinematicsRD[4][4]->GetBinContent(tt)));
  }
  fKinematicsRD[4][4]->Draw("E2");
  fKinematicsRD[4][4]->SetMarkerStyle(22);
  fKinematicsRD[4][4]->Draw("PSAME");
  fKinematicsRD[4][4]->GetXaxis()->SetLabelSize(0.03);
  fKinematicsRD[4][4]->GetYaxis()->SetLabelSize(0.03);
  for(int tt=0; tt<fKinematicsMC[4][4]->GetNbinsX(); tt++)
  {
    fKinematicsMC[4][4]->SetBinError(tt,sqrt(fKinematicsMC[4][4]->GetBinContent(tt)));
  }
  fKinematicsMC[4][4]->Draw("E2SAME");
  fKinematicsMC[4][4]->Draw("SAME");
  c12.Update();

  // c13.cd(2);
  // for(int tt=0; tt<fKinematicsRD[4][5]->GetNbinsX(); tt++)
  // {
  //   fError.push_back((fKinematicsRD[4][5]->GetBinError(tt) && fKinematicsMC[4][5]->GetBinError(tt) ? sqrt(pow(1/fKinematicsRD[4][5]->GetBinError(tt),2)+pow(1/fKinematicsMC[4][5]->GetBinError(tt),2)):0));
  // }
  // fKinematicsRD[4][5]->Scale(1/fKinematicsRD[4][5]->GetEntries());
  // fKinematicsMC[4][5]->Scale(1/fKinematicsMC[4][5]->GetEntries());
  // fKinematicsRatio[4][5] = (TH1F*)fKinematicsRD[4][5]->Clone();
  // fKinematicsRatio[4][5]->SetStats(0);
  // fKinematicsRatio[4][5]->Divide(fKinematicsMC[4][5]);
  // for(int tt=0; tt<fKinematicsRatio[4][5]->GetNbinsX(); tt++)
  // {
  //   fKinematicsRatio[4][5]->SetBinError(tt,fError[tt]);
  // }
  // fError.clear();
  // fKinematicsRatio[4][5]->SetMarkerStyle(21);
  // fKinematicsRatio[4][5]->SetFillColor(kYellow-7);
  // fKinematicsRatio[4][5]->SetMaximum(2.);
  // fKinematicsRatio[4][5]->SetMinimum(0.);
  // fKinematicsRatio[4][5]->Draw("PE2");
  // fKinematicsRatio[4][5]->GetXaxis()->SetLabelSize(0.08);
  // fKinematicsRatio[4][5]->GetYaxis()->SetLabelSize(0.08);
  // fKinematicsRatio[4][5]->GetYaxis()->SetNdivisions(2,kFALSE);
  // for(int tt=0; tt<7; tt++)
  // {
  //   l1[5][tt]->Draw();
  // }
  // c13.Update();
  c13.cd(1);
  fKinematicsMC[4][5]->Scale(fKinematicsRD[4][5]->GetEntries()/fKinematicsMC[4][5]->GetEntries());
  fKinematicsRD[4][5]->SetLineColor(kRed);
  fKinematicsMC[4][5]->SetLineColor(kBlue);
  fKinematicsRD[4][5]->SetMinimum(0.);
  fKinematicsRD[4][5]->SetMaximum(max(fKinematicsRD[4][5]->GetMaximum()*1.2,fKinematicsMC[4][5]->GetMaximum()*1.2));
  fKinematicsRD[4][5]->GetYaxis()->SetNdivisions(304,kTRUE);
  for(int tt=0; tt<fKinematicsRD[4][5]->GetNbinsX(); tt++)
  {
    fKinematicsRD[4][5]->SetBinError(tt,sqrt(fKinematicsRD[4][5]->GetBinContent(tt)));
  }
  fKinematicsRD[4][5]->Draw("E2");
  fKinematicsRD[4][5]->SetMarkerStyle(22);
  fKinematicsRD[4][5]->Draw("PSAME");
  fKinematicsRD[4][5]->GetXaxis()->SetLabelSize(0.03);
  fKinematicsRD[4][5]->GetYaxis()->SetLabelSize(0.03);
  for(int tt=0; tt<fKinematicsMC[4][5]->GetNbinsX(); tt++)
  {
    fKinematicsMC[4][5]->SetBinError(tt,sqrt(fKinematicsMC[4][5]->GetBinContent(tt)));
  }
  fKinematicsMC[4][5]->Draw("E2SAME");
  fKinematicsMC[4][5]->Draw("SAME");
  c13.Update();

  c33.cd(1);
  fKinematicsMC[4][6]->Scale(fKinematicsRD[4][6]->GetEntries()/fKinematicsMC[4][6]->GetEntries());
  fKinematicsRD[4][6]->SetLineColor(kRed);
  fKinematicsMC[4][6]->SetLineColor(kBlue);
  fKinematicsRD[4][6]->SetMinimum(0.);
  fKinematicsRD[4][6]->SetMaximum(max(fKinematicsRD[4][6]->GetMaximum()*1.2,fKinematicsMC[4][6]->GetMaximum()*1.2));
  fKinematicsRD[4][6]->GetYaxis()->SetNdivisions(304,kTRUE);
  for(int tt=0; tt<fKinematicsRD[4][6]->GetNbinsX(); tt++)
  {
    fKinematicsRD[4][6]->SetBinError(tt,sqrt(fKinematicsRD[4][6]->GetBinContent(tt)));
  }
  fKinematicsRD[4][6]->Draw("E2");
  fKinematicsRD[4][6]->SetMarkerStyle(22);
  fKinematicsRD[4][6]->Draw("PSAME");
  fKinematicsRD[4][6]->GetXaxis()->SetLabelSize(0.03);
  fKinematicsRD[4][6]->GetYaxis()->SetLabelSize(0.03);
  for(int tt=0; tt<fKinematicsMC[4][6]->GetNbinsX(); tt++)
  {
    fKinematicsMC[4][6]->SetBinError(tt,sqrt(fKinematicsMC[4][6]->GetBinContent(tt)));
  }
  fKinematicsMC[4][6]->Draw("E2SAME");
  fKinematicsMC[4][6]->Draw("SAME");
  c33.Update();

  // c32.cd(1);
  // for(int tt=0; tt<fKinematicsRD[4][5]->GetNbinsX(); tt++)
  // {
  //   fKinematicsRD[4][5]->SetBinError(tt,sqrt(fKinematicsRD[4][5]->GetBinContent(tt)));
  // }
  // fKinematicsRD[4][5]->Draw("E2");
  // fKinematicsRD[4][5]->Draw("SAME");
  // fKinematicsRD[4][5]->GetXaxis()->SetLabelSize(0.03);
  // fKinematicsRD[4][5]->GetYaxis()->SetLabelSize(0.03);
  // for(int tt=0; tt<fKinematicsMC[4][5]->GetNbinsX(); tt++)
  // {
  //   fKinematicsMC[4][5]->SetBinError(tt,sqrt(fKinematicsMC[4][5]->GetBinContent(tt)));
  // }
  // fKinematicsMC[4][5]->Draw("E2SAME");
  // fKinematicsMC[4][5]->Draw("SAME");
  // c32.Update();

  // c22.cd(2);
  // for(int tt=0; tt<fKinematicsRD[4][12]->GetNbinsX(); tt++)
  // {
  //   fError.push_back((fKinematicsRD[4][12]->GetBinError(tt) && fKinematicsMC[4][12]->GetBinError(tt) ? sqrt(pow(1/fKinematicsRD[4][12]->GetBinError(tt),2)+pow(1/fKinematicsMC[4][12]->GetBinError(tt),2)):0));
  // }
  // fKinematicsRD[4][12]->Scale(1/fKinematicsRD[4][12]->GetEntries());
  // fKinematicsMC[4][12]->Scale(1/fKinematicsMC[4][12]->GetEntries());
  // fKinematicsRatio[4][12] = (TH1F*)fKinematicsRD[4][12]->Clone();
  // fKinematicsRatio[4][12]->SetStats(0);
  // fKinematicsRatio[4][12]->Divide(fKinematicsMC[4][12]);
  // for(int tt=0; tt<fKinematicsRatio[4][12]->GetNbinsX(); tt++)
  // {
  //   fKinematicsRatio[4][12]->SetBinError(tt,fError[tt]);
  // }
  // fError.clear();
  // fKinematicsRatio[4][12]->SetMarkerStyle(21);
  // fKinematicsRatio[4][12]->SetFillColor(kYellow-7);
  // fKinematicsRatio[4][12]->SetMaximum(2.);
  // fKinematicsRatio[4][12]->SetMinimum(0.);
  // fKinematicsRatio[4][12]->Draw("PE2");
  // fKinematicsRatio[4][12]->GetXaxis()->SetLabelSize(0.08);
  // fKinematicsRatio[4][12]->GetYaxis()->SetLabelSize(0.08);
  // fKinematicsRatio[4][12]->GetYaxis()->SetNdivisions(2,kFALSE);
  // for(int tt=0; tt<7; tt++)
  // {
  //   l1[12][tt]->Draw();
  // }
  // c22.Update();
  c22.cd(1);
  fKinematicsMC[4][12]->Scale(fKinematicsRD[4][12]->GetEntries()/fKinematicsMC[4][12]->GetEntries());
  fKinematicsRD[4][12]->SetLineColor(kRed);
  fKinematicsMC[4][12]->SetLineColor(kBlue);
  fKinematicsRD[4][12]->SetMinimum(0.);
  fKinematicsRD[4][12]->SetMaximum(max(fKinematicsRD[4][12]->GetMaximum()*1.2,fKinematicsMC[4][12]->GetMaximum()*1.2));
  fKinematicsRD[4][12]->GetYaxis()->SetNdivisions(304,kTRUE);
  for(int tt=0; tt<fKinematicsRD[4][12]->GetNbinsX(); tt++)
  {
    fKinematicsRD[4][12]->SetBinError(tt,sqrt(fKinematicsRD[4][12]->GetBinContent(tt)));
  }
  fKinematicsRD[4][12]->Draw("E2");
  fKinematicsRD[4][12]->SetMarkerStyle(22);
  fKinematicsRD[4][12]->Draw("PSAME");
  fKinematicsRD[4][12]->GetXaxis()->SetLabelSize(0.03);
  fKinematicsRD[4][12]->GetYaxis()->SetLabelSize(0.03);
  for(int tt=0; tt<fKinematicsMC[4][12]->GetNbinsX(); tt++)
  {
    fKinematicsMC[4][12]->SetBinError(tt,sqrt(fKinematicsMC[4][12]->GetBinContent(tt)));
  }
  fKinematicsMC[4][12]->Draw("E2SAME");
  fKinematicsMC[4][12]->Draw("SAME");
  c22.Update();

  // c23.cd(2);
  // for(int tt=0; tt<fKinematicsRD[4][13]->GetNbinsX(); tt++)
  // {
  //   fError.push_back((fKinematicsRD[4][13]->GetBinError(tt) && fKinematicsMC[4][13]->GetBinError(tt) ? sqrt(pow(1/fKinematicsRD[4][13]->GetBinError(tt),2)+pow(1/fKinematicsMC[4][13]->GetBinError(tt),2)):0));
  // }
  // fKinematicsRD[4][13]->Scale(1/fKinematicsRD[4][13]->GetEntries());
  // fKinematicsMC[4][13]->Scale(1/fKinematicsMC[4][13]->GetEntries());
  // fKinematicsRatio[4][13] = (TH1F*)fKinematicsRD[4][13]->Clone();
  // fKinematicsRatio[4][13]->SetStats(0);
  // fKinematicsRatio[4][13]->Divide(fKinematicsMC[4][13]);
  // for(int tt=0; tt<fKinematicsRatio[4][13]->GetNbinsX(); tt++)
  // {
  //   fKinematicsRatio[4][13]->SetBinError(tt,fError[tt]);
  // }
  // fError.clear();
  // fKinematicsRatio[4][13]->SetMarkerStyle(21);
  // fKinematicsRatio[4][13]->SetFillColor(kYellow-7);
  // fKinematicsRatio[4][13]->SetMaximum(2.);
  // fKinematicsRatio[4][13]->SetMinimum(0.);
  // fKinematicsRatio[4][13]->Draw("PE2");
  // fKinematicsRatio[4][13]->GetXaxis()->SetLabelSize(0.08);
  // fKinematicsRatio[4][13]->GetYaxis()->SetLabelSize(0.08);
  // fKinematicsRatio[4][13]->GetYaxis()->SetNdivisions(2,kFALSE);
  // for(int tt=0; tt<7; tt++)
  // {
  //   l1[13][tt]->Draw();
  // }
  // c23.Update();
  c23.cd(1);
  fKinematicsMC[4][13]->Scale(fKinematicsRD[4][13]->GetEntries()/fKinematicsMC[4][13]->GetEntries());
  fKinematicsRD[4][13]->SetLineColor(kRed);
  fKinematicsMC[4][13]->SetLineColor(kBlue);
  fKinematicsRD[4][13]->SetMinimum(0.);
  fKinematicsRD[4][13]->SetMaximum(max(fKinematicsRD[4][13]->GetMaximum()*1.2,fKinematicsMC[4][13]->GetMaximum()*1.2));
  fKinematicsRD[4][13]->GetYaxis()->SetNdivisions(304,kTRUE);
  for(int tt=0; tt<fKinematicsRD[4][13]->GetNbinsX(); tt++)
  {
    fKinematicsRD[4][13]->SetBinError(tt,sqrt(fKinematicsRD[4][13]->GetBinContent(tt)));
  }
  fKinematicsRD[4][13]->Draw("E2");
  fKinematicsRD[4][13]->SetMarkerStyle(22);
  fKinematicsRD[4][13]->Draw("PSAME");
  fKinematicsRD[4][13]->GetXaxis()->SetLabelSize(0.03);
  fKinematicsRD[4][13]->GetYaxis()->SetLabelSize(0.03);
  for(int tt=0; tt<fKinematicsMC[4][13]->GetNbinsX(); tt++)
  {
    fKinematicsMC[4][13]->SetBinError(tt,sqrt(fKinematicsMC[4][13]->GetBinContent(tt)));
  }
  fKinematicsMC[4][13]->Draw("E2SAME");
  fKinematicsMC[4][13]->Draw("SAME");
  c23.Update();

  // c24.cd(2);
  // for(int tt=0; tt<fKinematicsRD[4][14]->GetNbinsX(); tt++)
  // {
  //   fError.push_back((fKinematicsRD[4][14]->GetBinError(tt) && fKinematicsMC[4][14]->GetBinError(tt) ? sqrt(pow(1/fKinematicsRD[4][14]->GetBinError(tt),2)+pow(1/fKinematicsMC[4][14]->GetBinError(tt),2)):0));
  // }
  // fKinematicsRD[4][14]->Scale(1/fKinematicsRD[4][14]->GetEntries());
  // fKinematicsMC[4][14]->Scale(1/fKinematicsMC[4][14]->GetEntries());
  // fKinematicsRatio[4][14] = (TH1F*)fKinematicsRD[4][14]->Clone();
  // fKinematicsRatio[4][14]->SetStats(0);
  // fKinematicsRatio[4][14]->Divide(fKinematicsMC[4][14]);
  // for(int tt=0; tt<fKinematicsRatio[4][14]->GetNbinsX(); tt++)
  // {
  //   fKinematicsRatio[4][14]->SetBinError(tt,fError[tt]);
  // }
  // fError.clear();
  // fKinematicsRatio[4][14]->SetMarkerStyle(21);
  // fKinematicsRatio[4][14]->SetFillColor(kYellow-7);
  // fKinematicsRatio[4][14]->SetMaximum(2.);
  // fKinematicsRatio[4][14]->SetMinimum(0.);
  // fKinematicsRatio[4][14]->Draw("PE2");
  // fKinematicsRatio[4][14]->GetXaxis()->SetLabelSize(0.08);
  // fKinematicsRatio[4][14]->GetYaxis()->SetLabelSize(0.08);
  // fKinematicsRatio[4][14]->GetYaxis()->SetNdivisions(2,kFALSE);
  // for(int tt=0; tt<7; tt++)
  // {
  //   l1[14][tt]->Draw();
  // }
  // c24.Update();
  c24.cd(1);
  fKinematicsMC[4][14]->Scale(fKinematicsRD[4][14]->GetEntries()/fKinematicsMC[4][14]->GetEntries());
  fKinematicsRD[4][14]->SetLineColor(kRed);
  fKinematicsMC[4][14]->SetLineColor(kBlue);
  fKinematicsRD[4][14]->SetMinimum(0.);
  fKinematicsRD[4][14]->SetMaximum(max(fKinematicsRD[4][14]->GetMaximum()*1.2,fKinematicsMC[4][14]->GetMaximum()*1.2));
  fKinematicsRD[4][14]->GetYaxis()->SetNdivisions(304,kTRUE);
  for(int tt=0; tt<fKinematicsRD[4][14]->GetNbinsX(); tt++)
  {
    fKinematicsRD[4][14]->SetBinError(tt,sqrt(fKinematicsRD[4][14]->GetBinContent(tt)));
  }
  fKinematicsRD[4][14]->Draw("E2");
  fKinematicsRD[4][14]->SetMarkerStyle(22);
  fKinematicsRD[4][14]->Draw("PSAME");
  fKinematicsRD[4][14]->GetXaxis()->SetLabelSize(0.03);
  fKinematicsRD[4][14]->GetYaxis()->SetLabelSize(0.03);
  for(int tt=0; tt<fKinematicsMC[4][14]->GetNbinsX(); tt++)
  {
    fKinematicsMC[4][14]->SetBinError(tt,sqrt(fKinematicsMC[4][14]->GetBinContent(tt)));
  }
  fKinematicsMC[4][14]->Draw("E2SAME");
  fKinematicsMC[4][14]->Draw("SAME");
  c24.Update();

  // c26.cd(2);
  // for(int tt=0; tt<fKinematicsRD[4][15]->GetNbinsX(); tt++)
  // {
  //   fError.push_back((fKinematicsRD[4][15]->GetBinError(tt) && fKinematicsMC[4][15]->GetBinError(tt) ? sqrt(pow(1/fKinematicsRD[4][15]->GetBinError(tt),2)+pow(1/fKinematicsMC[4][15]->GetBinError(tt),2)):0));
  // }
  // fKinematicsRD[4][15]->Scale(1/fKinematicsRD[4][15]->GetEntries());
  // fKinematicsMC[4][15]->Scale(1/fKinematicsMC[4][15]->GetEntries());
  // fKinematicsRatio[4][15] = (TH1F*)fKinematicsRD[4][15]->Clone();
  // fKinematicsRatio[4][15]->SetStats(0);
  // fKinematicsRatio[4][15]->Divide(fKinematicsMC[4][15]);
  // for(int tt=0; tt<fKinematicsRatio[4][15]->GetNbinsX(); tt++)
  // {
  //   fKinematicsRatio[4][15]->SetBinError(tt,fError[tt]);
  // }
  // fError.clear();
  // fKinematicsRatio[4][15]->SetMarkerStyle(21);
  // fKinematicsRatio[4][15]->SetFillColor(kYellow-7);
  // fKinematicsRatio[4][15]->SetMaximum(2.);
  // fKinematicsRatio[4][15]->SetMinimum(0.);
  // fKinematicsRatio[4][15]->Draw("PE2");
  // fKinematicsRatio[4][15]->GetXaxis()->SetLabelSize(0.08);
  // fKinematicsRatio[4][15]->GetYaxis()->SetLabelSize(0.08);
  // fKinematicsRatio[4][15]->GetYaxis()->SetNdivisions(2,kFALSE);
  // for(int tt=0; tt<7; tt++)
  // {
  //   l1[15][tt]->Draw();
  // }
  // c26.Update();
  c26.cd(1);
  fKinematicsMC[4][15]->Scale(fKinematicsRD[4][15]->GetEntries()/fKinematicsMC[4][15]->GetEntries());
  fKinematicsRD[4][15]->SetLineColor(kRed);
  fKinematicsMC[4][15]->SetLineColor(kBlue);
  fKinematicsRD[4][15]->SetMinimum(0.);
  fKinematicsRD[4][15]->SetMaximum(max(fKinematicsRD[4][15]->GetMaximum()*1.2,fKinematicsMC[4][15]->GetMaximum()*1.2));
  fKinematicsRD[4][15]->GetYaxis()->SetNdivisions(304,kTRUE);
  for(int tt=0; tt<fKinematicsRD[4][15]->GetNbinsX(); tt++)
  {
    fKinematicsRD[4][15]->SetBinError(tt,sqrt(fKinematicsRD[4][15]->GetBinContent(tt)));
  }
  fKinematicsRD[4][15]->Draw("E2");
  fKinematicsRD[4][15]->SetMarkerStyle(22);
  fKinematicsRD[4][15]->Draw("PSAME");
  fKinematicsRD[4][15]->GetXaxis()->SetLabelSize(0.03);
  fKinematicsRD[4][15]->GetYaxis()->SetLabelSize(0.03);
  for(int tt=0; tt<fKinematicsMC[4][15]->GetNbinsX(); tt++)
  {
    fKinematicsMC[4][15]->SetBinError(tt,sqrt(fKinematicsMC[4][15]->GetBinContent(tt)));
  }
  fKinematicsMC[4][15]->Draw("E2SAME");
  fKinematicsMC[4][15]->Draw("SAME");
  c26.Update();

  // c28.cd(2);
  // for(int tt=0; tt<fKinematicsRD[4][16]->GetNbinsX(); tt++)
  // {
  //   fError.push_back((fKinematicsRD[4][16]->GetBinError(tt) && fKinematicsMC[4][16]->GetBinError(tt) ? sqrt(pow(1/fKinematicsRD[4][16]->GetBinError(tt),2)+pow(1/fKinematicsMC[4][16]->GetBinError(tt),2)):0));
  // }
  // fKinematicsRD[4][16]->Scale(1/fKinematicsRD[4][16]->GetEntries());
  // fKinematicsMC[4][16]->Scale(1/fKinematicsMC[4][16]->GetEntries());
  // fKinematicsRatio[4][16] = (TH1F*)fKinematicsRD[4][16]->Clone();
  // fKinematicsRatio[4][16]->SetStats(0);
  // fKinematicsRatio[4][16]->Divide(fKinematicsMC[4][16]);
  // for(int tt=0; tt<fKinematicsRatio[4][16]->GetNbinsX(); tt++)
  // {
  //   fKinematicsRatio[4][16]->SetBinError(tt,fError[tt]);
  // }
  // fError.clear();
  // fKinematicsRatio[4][16]->SetMarkerStyle(21);
  // fKinematicsRatio[4][16]->SetFillColor(kYellow-7);
  // fKinematicsRatio[4][16]->SetMaximum(2.);
  // fKinematicsRatio[4][16]->SetMinimum(0.);
  // fKinematicsRatio[4][16]->Draw("PE2");
  // fKinematicsRatio[4][16]->GetXaxis()->SetLabelSize(0.08);
  // fKinematicsRatio[4][16]->GetYaxis()->SetLabelSize(0.08);
  // fKinematicsRatio[4][16]->GetYaxis()->SetNdivisions(2,kFALSE);
  // for(int tt=0; tt<7; tt++)
  // {
  //   l1[16][tt]->Draw();
  // }
  // c28.Update();
  c28.cd(1);
  fKinematicsMC[4][16]->Scale(fKinematicsRD[4][16]->GetEntries()/fKinematicsMC[4][16]->GetEntries());
  fKinematicsRD[4][16]->SetLineColor(kRed);
  fKinematicsMC[4][16]->SetLineColor(kBlue);
  fKinematicsRD[4][16]->SetMinimum(0.);
  fKinematicsRD[4][16]->SetMaximum(max(fKinematicsRD[4][16]->GetMaximum()*1.2,fKinematicsMC[4][16]->GetMaximum()*1.2));
  fKinematicsRD[4][16]->GetYaxis()->SetNdivisions(304,kTRUE);
  for(int tt=0; tt<fKinematicsRD[4][16]->GetNbinsX(); tt++)
  {
    fKinematicsRD[4][16]->SetBinError(tt,(fKinematicsRD[4][16]->GetBinContent(tt) ? 1/sqrt(fKinematicsRD[4][16]->GetBinContent(tt)):0));
  }
  fKinematicsRD[4][16]->Draw("E2");
  fKinematicsRD[4][16]->SetMarkerStyle(22);
  fKinematicsRD[4][16]->Draw("PSAME");
  fKinematicsRD[4][16]->GetXaxis()->SetLabelSize(0.03);
  fKinematicsRD[4][16]->GetYaxis()->SetLabelSize(0.03);
  for(int tt=0; tt<fKinematicsMC[4][16]->GetNbinsX(); tt++)
  {
    fKinematicsMC[4][16]->SetBinError(tt,(fKinematicsMC[4][16]->GetBinContent(tt) ? 1/sqrt(fKinematicsMC[4][16]->GetBinContent(tt)):0));
  }
  fKinematicsMC[4][16]->Draw("E2SAME");
  fKinematicsMC[4][16]->Draw("SAME");
  c28.Update();

  c34.cd(1);
  fThetaRDp[0]->Draw("COLZ");
  fThetaRDp[0]->GetXaxis()->SetTitle("#theta_y");
  fThetaRDp[0]->GetYaxis()->SetTitle("p");
  c34.Update();

  c34.cd(2);
  fThetaMCp[0]->Draw("COLZ");
  fThetaMCp[0]->GetXaxis()->SetTitle("#theta_y");
  fThetaMCp[0]->GetYaxis()->SetTitle("p");
  c34.Update();

  c35.cd(1);
  fThetaRDp[1]->Draw("COLZ");
  fThetaRDp[1]->GetXaxis()->SetTitle("#theta_x");
  fThetaRDp[1]->GetYaxis()->SetTitle("p");
  c35.Update();

  c35.cd(2);
  fThetaMCp[1]->Draw("COLZ");
  fThetaMCp[1]->GetXaxis()->SetTitle("#theta_x");
  fThetaMCp[1]->GetYaxis()->SetTitle("p");
  c35.Update();

  c36.cd(1);
  fThetaRDp[2]->Draw("COLZ");
  fThetaRDp[2]->GetXaxis()->SetTitle("#theta_y");
  fThetaRDp[2]->GetYaxis()->SetTitle("#theta_x");
  c36.Update();

  c36.cd(2);
  fThetaMCp[2]->Draw("COLZ");
  fThetaMCp[2]->GetXaxis()->SetTitle("#theta_y");
  fThetaMCp[2]->GetYaxis()->SetTitle("#theta_x");
  c36.Update();

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
  c18.Print("kinMCRD.pdf)","pdf");
  c29.Print("photoelec.pdf","pdf");
  c30.Print("Q2.pdf","pdf");
  c31.Print("z.pdf","pdf");
  c32.Print("nu.pdf","pdf");
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
      fBP++;

      // Reconstructed muon
      if((0<E_beam->GetLeaf("E_beam")->GetValue()))
      {
        fRmu++;

        //BMS (reconstructed beam track)
        if(true) //not used in acceptance
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
                  if(true)
                  {
                    fCell++;

                    // IM/O triggers
                    if((trig&8 || trig&256))
                    {
                      fTrig++;

                      // Q2 cut
                      if((Q2>1))
                      {
                        fQ2test++;

                        // y cut
                        if((fYmin<yBj && yBj<fYmax))
                        {
                          fYBjtest++;

                          // W cut
                          if((fWmin<sqrt(wBj) && sqrt(wBj)<fWmax))
                          {
                            fWBjtest++;

                            // x cut
                            if((fXmin<xBj && xBj<fXmax))
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
              if(InTarget(x->GetLeaf("x")->GetValue(),y->GetLeaf("y")->GetValue(),z->GetLeaf("z")->GetValue()))
              {
                fTarg++;

                // Cells crossing
                if(true)
                {
                  fCell++;

                  if((trig&2 || trig&4 || trig&8))
                  {
                    fTrig++;

                    // Q2 cut
                    if((Q2>1))
                    {
                      fQ2test++;

                      // y cut
                      if((fYmin<yBj && yBj<fYmax))
                      {
                        fYBjtest++;

                        // W cut
                        if((fWmin<sqrt(wBj) && sqrt(wBj)<fWmax))
                        {
                          fWBjtest++;
                          if((fXmin<xBj && xBj<fXmax))
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
              if(InTarget(x->GetLeaf("x")->GetValue(),y->GetLeaf("y")->GetValue(),z->GetLeaf("z")->GetValue()))
              {
                fTarg++;

                // Cells crossing
                if(true)
                {
                  fCell++;

                  if((int(trig&2) || int(trig&4) || int(trig&8) || int(trig&512)))
                  {
                    fTrig++;

                    // Q2 cut
                    if((Q2>1))
                    {
                      fQ2test++;
                      fMuMC[4].push_back(E_beam->GetLeaf("E_beam")->GetValue());
                      fThetaMCMu[2].push_back(sqrt(pow(p0x->GetLeaf("p0x")->GetValue(),2)
                                                  +pow(p0y->GetLeaf("p0y")->GetValue(),2)
                                                  +pow(p0z->GetLeaf("p0z")->GetValue(),2)));
                      fThetaMCMu[1].push_back(thetax_b);
                      fThetaMCMu[0].push_back(thetay_b);

                      // y cut
                      if((fYmin<yBj && yBj<fYmax))
                      {
                        fYBjtest++;

                        // W cut
                        if((fWmin<sqrt(wBj) && sqrt(wBj)<fWmax))
                        {
                          fWBjtest++;
                          if((fXmin<xBj && xBj<fXmax))
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

      // If all DIS tests are good, then event is saved
      if(fAllDISflag)
      {
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
          fXX0test++;

          // Momentum cut (12 GeV to 40 GeV, increasing to 3 GeV to 40 GeV)
          if(!(fPmin<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<fPmax)) continue;
          fMom++;

          // Theta cut
          if(!(0.01<thRICH->GetLeaf("Hadrons.thRICH")->GetValue(i) && thRICH->GetLeaf("Hadrons.thRICH")->GetValue(i)<0.12)) continue;
          fTRICH++;

          // RICH position cut
          if(!(pow(RICHx->GetLeaf("Hadrons.RICHx")->GetValue(i),2)+pow(RICHy->GetLeaf("Hadrons.RICHy")->GetValue(i),2)>25)) continue;
          fPosRICH++;

          // MT
          if(int(trig&2) && !int(trig&4) && !int(trig&8) && !int(trig&512))
          {
            fKinematicsMC[0][3]->Fill(zBj);
            fKinematicsMC[0][12]->Fill(p->GetLeaf("Hadrons.P")->GetValue(i));
            fKinematicsMC[0][13]->Fill(th->GetLeaf("Hadrons.th")->GetValue(i));
            fKinematicsMC[0][14]->Fill(ph->GetLeaf("Hadrons.ph")->GetValue(i));
            fKinematicsMC[0][15]->Fill(ph_pl->GetLeaf("Hadrons.ph_pl")->GetValue(i));
            fKinematicsMC[0][16]->Fill(pt->GetLeaf("Hadrons.pt")->GetValue(i));
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
    fKinematicsMC[4][1]->Fill(fXBjkinMC[4][i]);
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
    TBranch *pt = (TBranch*) tree->FindBranch("Hadrons.pt");
    TBranch *th = (TBranch*) tree->FindBranch("Hadrons.th");
    TBranch *ph = (TBranch*) tree->FindBranch("Hadrons.ph");
    TBranch *ph_pl = (TBranch*) tree->FindBranch("Hadrons.ph_pl");
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
      pt->GetEntry(ip);
      th->GetEntry(ip);
      ph->GetEntry(ip);
      ph_pl->GetEntry(ip);
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
      fBP++;

      // Reconstructed muon
      if(!(0<E_beam->GetLeaf("E_beam")->GetValue())) continue;
      fRmu++;

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
        if(!InTarget(x->GetLeaf("x")->GetValue(),y->GetLeaf("y")->GetValue(),z->GetLeaf("z")->GetValue())) continue;
      }
      //2012 ---
      //2016 ---
      else if(Y2016)
      {
        if(!InTarget(x->GetLeaf("x")->GetValue(),y->GetLeaf("y")->GetValue(),z->GetLeaf("z")->GetValue())) continue;
      }
      //2016 ---
      fTarg++;

      // Cells crossing
      if(/*!(cellsCrossed->GetLeaf("cellsCrossed")->GetValue())*/false) continue;
      fCell++;

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
      fTrig++;

      // Q2 cut
      if(!(Q2>1)) continue;
      fQ2test++;
      fMu[4].push_back(E_beam->GetLeaf("E_beam")->GetValue());

      fThetaMu[2].push_back(sqrt(pow(p0x->GetLeaf("p0x")->GetValue(),2)
                                +pow(p0y->GetLeaf("p0y")->GetValue(),2)
                                +pow(p0z->GetLeaf("p0z")->GetValue(),2)));
      fThetaMu[1].push_back(thetax_b);
      fThetaMu[0].push_back(thetay_b);

      // y cut
      if(!(fYmin<yBj && yBj<fYmax)) continue;
      fYBjtest++;

      // W cut
      if(!(fWmin<sqrt(wBj) && sqrt(wBj)<fWmax)) continue;
      fWBjtest++;

      // x cut
      if(!(fXmin<xBj && xBj<fXmax)) continue;
      fXBjtest++;

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
        fXX0test++;

        // Momentum cut (12 GeV to 40 GeV, increasing to 3 GeV to 40 GeV)
        if(!(fPmin<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<fPmax)) continue;
        fMom++;

        // Theta cut
        if(!(0.01<thRICH->GetLeaf("Hadrons.thRICH")->GetValue(i) && thRICH->GetLeaf("Hadrons.thRICH")->GetValue(i)<0.12)) continue;
        fTRICH++;

        // RICH position cut
        if(!(pow(RICHx->GetLeaf("Hadrons.RICHx")->GetValue(i),2)+pow(RICHy->GetLeaf("Hadrons.RICHy")->GetValue(i),2)>25)) continue;
        fPosRICH++;

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
        }
        if(int(trig&4) && !int(trig&2) && !int(trig&8)&& !int(trig&512))
        {
          fKinematicsRD[1][3]->Fill(zBj);
          fKinematicsRD[1][12]->Fill(p->GetLeaf("Hadrons.P")->GetValue(i));
          fKinematicsRD[1][13]->Fill(th->GetLeaf("Hadrons.th")->GetValue(i));
          fKinematicsRD[1][14]->Fill(ph->GetLeaf("Hadrons.ph")->GetValue(i));
          fKinematicsRD[1][15]->Fill(ph_pl->GetLeaf("Hadrons.ph_pl")->GetValue(i));
          fKinematicsRD[1][16]->Fill(pt->GetLeaf("Hadrons.pt")->GetValue(i));
        }
        if(int(trig&8) && !int(trig&2) && !int(trig&4) && !int(trig&512))
        {
          fKinematicsRD[2][3]->Fill(zBj);
          fKinematicsRD[2][12]->Fill(p->GetLeaf("Hadrons.P")->GetValue(i));
          fKinematicsRD[2][13]->Fill(th->GetLeaf("Hadrons.th")->GetValue(i));
          fKinematicsRD[2][14]->Fill(ph->GetLeaf("Hadrons.ph")->GetValue(i));
          fKinematicsRD[2][15]->Fill(ph_pl->GetLeaf("Hadrons.ph_pl")->GetValue(i));
          fKinematicsRD[2][16]->Fill(pt->GetLeaf("Hadrons.pt")->GetValue(i));
        }
        if(int(trig&512) && !int(trig&4) && !int(trig&8) && !int(trig&2))
        {
          fKinematicsRD[3][3]->Fill(zBj);
          fKinematicsRD[3][12]->Fill(p->GetLeaf("Hadrons.P")->GetValue(i));
          fKinematicsRD[3][13]->Fill(th->GetLeaf("Hadrons.th")->GetValue(i));
          fKinematicsRD[3][14]->Fill(ph->GetLeaf("Hadrons.ph")->GetValue(i));
          fKinematicsRD[3][15]->Fill(ph_pl->GetLeaf("Hadrons.ph_pl")->GetValue(i));
          fKinematicsRD[3][16]->Fill(pt->GetLeaf("Hadrons.pt")->GetValue(i));
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
      fKinematicsRD[4][1]->Fill(fXBjkin[4][i]);
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

  return 0;
}

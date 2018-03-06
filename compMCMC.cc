#include "compMCMC.h"

//Inputs
#define target_file_2012 "target-107924-109081.dat"
#define target_file_2016 "target-274508-274901.dat"

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
    fKinematicsMC2[i][0] = new TH1F(Form("Q^{2} %s",trigname[i].c_str()), Form("Q^{2} %s",trigname[i].c_str()), 100, -1, 2);
    fKinematicsMC2[i][1] = new TH1F(Form("x_{Bj} %s",trigname[i].c_str()), Form("x_{Bj} %s",trigname[i].c_str()), 100, -3, 0);
    fKinematicsMC2[i][2] = new TH1F(Form("y %s",trigname[i].c_str()), Form("y %s",trigname[i].c_str()), 100, 0, 1);
    fKinematicsMC2[i][3] = new TH1F(Form("z %s",trigname[i].c_str()), Form("z %s",trigname[i].c_str()), 100, 0, 1);
    fKinematicsMC2[i][4] = new TH1F(Form("W %s",trigname[i].c_str()), Form("W %s",trigname[i].c_str()), 100, 2, 18);
    fKinematicsMC2[i][5] = new TH1F(Form("#nu %s",trigname[i].c_str()), Form("#nu %s",trigname[i].c_str()), 100, 0, 160);
    fKinematicsMC2[i][6] = new TH1F(Form("E_{#mu} %s",trigname[i].c_str()), Form("E_{#mu} %s",trigname[i].c_str()), 100, 140, 180);
    fKinematicsMC2[i][7] = new TH1F(Form("E_{#mu'} %s",trigname[i].c_str()), Form("E_{#mu'} %s",trigname[i].c_str()), 100, 0, 160);
    fKinematicsMC2[i][8] = new TH1F(Form("#theta %s",trigname[i].c_str()), Form("#theta %s",trigname[i].c_str()), 100, 0, 0.05);
    fKinematicsMC2[i][9] = new TH1F(Form("#phi %s",trigname[i].c_str()), Form("#phi %s",trigname[i].c_str()), 100, -1.7, 1.7);
    fKinematicsMC2[i][10] = new TH1F(Form("Vertex %s",trigname[i].c_str()), Form("Vertex %s",trigname[i].c_str()), 100, -320, -70);
    fKinematicsMC2[i][12] = new TH1F(Form("p_{hadron+e} %s",trigname[i].c_str()), Form("p_{hadron+e} %s",trigname[i].c_str()), 100, 0, 40);
    fKinematicsMC2[i][13] = new TH1F(Form("#theta_{hadron+e} %s",trigname[i].c_str()), Form("#theta_{hadron+e} %s",trigname[i].c_str()), 100, 0, 0.25);
    fKinematicsMC2[i][14] = new TH1F(Form("#phi_{hadron+e,lab} %s",trigname[i].c_str()), Form("#phi_{hadron+e,lab} %s",trigname[i].c_str()), 100, -3.5, 3.5);
    fKinematicsMC2[i][15] = new TH1F(Form("#phi_{hadron+e,prod.pl} %s",trigname[i].c_str()), Form("#phi_{hadron+e,prod.pl} %s",trigname[i].c_str()), 100, 0, 3.5);
    fKinematicsMC2[i][16] = new TH1F(Form("p_{T} %s",trigname[i].c_str()), Form("p_{T} %s",trigname[i].c_str()), 100, 0, 3);
    fKinematicsMC1[i][0] = new TH1F(Form("Q^{2} Ratio %s",trigname[i].c_str()), Form("Q^{2} Ratio %s",trigname[i].c_str()), 100, -1, 2);
    fKinematicsMC1[i][1] = new TH1F(Form("x_{Bj} Ratio %s",trigname[i].c_str()), Form("x_{Bj} Ratio %s",trigname[i].c_str()), 100, -3, 0);
    fKinematicsMC1[i][2] = new TH1F(Form("y Ratio %s",trigname[i].c_str()), Form("y Ratio %s",trigname[i].c_str()), 100, 0, 1);
    fKinematicsMC1[i][3] = new TH1F(Form("z Ratio %s",trigname[i].c_str()), Form("z Ratio %s",trigname[i].c_str()), 100, 0, 1);
    fKinematicsMC1[i][4] = new TH1F(Form("W Ratio %s",trigname[i].c_str()), Form("W Ratio %s",trigname[i].c_str()), 100, 2, 18);
    fKinematicsMC1[i][5] = new TH1F(Form("#nu Ratio %s",trigname[i].c_str()), Form("#nu Ratio %s",trigname[i].c_str()), 100, 0, 160);
    fKinematicsMC1[i][6] = new TH1F(Form("E_{#mu} Ratio %s",trigname[i].c_str()), Form("E_{#mu} Ratio %s",trigname[i].c_str()), 100, 140, 180);
    fKinematicsMC1[i][7] = new TH1F(Form("E_{#mu'} Ratio %s",trigname[i].c_str()), Form("E_{#mu'} Ratio %s",trigname[i].c_str()), 100, 0, 160);
    fKinematicsMC1[i][8] = new TH1F(Form("#theta Ratio %s",trigname[i].c_str()), Form("#theta Ratio %s",trigname[i].c_str()), 100, 0, 0.05);
    fKinematicsMC1[i][9] = new TH1F(Form("#phi Ratio %s",trigname[i].c_str()), Form("#phi Ratio %s",trigname[i].c_str()), 100, -1.7, 1.7);
    fKinematicsMC1[i][10] = new TH1F(Form("Vertex Ratio %s",trigname[i].c_str()), Form("Vertex Ratio %s",trigname[i].c_str()), 100, -320, -70);
    fKinematicsMC1[i][12] = new TH1F(Form("p_{hadron+e} Ratio %s",trigname[i].c_str()), Form("p_{hadron+e} Ratio %s",trigname[i].c_str()), 100, 0, 40);
    fKinematicsMC1[i][13] = new TH1F(Form("#theta_{hadron+e} Ratio %s",trigname[i].c_str()), Form("#theta_{hadron+e} Ratio %s",trigname[i].c_str()), 100, 0, 0.25);
    fKinematicsMC1[i][14] = new TH1F(Form("#phi_{hadron+e,lab} Ratio %s",trigname[i].c_str()), Form("#phi_{hadron+e,lab} Ratio %s",trigname[i].c_str()), 100, -3.5, 3.5);
    fKinematicsMC1[i][15] = new TH1F(Form("#phi_{hadron+e,prod.pl} Ratio %s",trigname[i].c_str()), Form("#phi_{hadron+e,prod.pl} Ratio %s",trigname[i].c_str()), 100, 0, 3.5);
    fKinematicsMC1[i][16] = new TH1F(Form("p_{T} Ratio %s",trigname[i].c_str()), Form("p_{T} Ratio %s",trigname[i].c_str()), 100, 0, 3);
    BinLogX(fKinematicsMC2[i][0]);
    BinLogX(fKinematicsMC1[i][0]);
    BinLogX(fKinematicsMC2[i][1]);
    BinLogX(fKinematicsMC1[i][1]);
  }
  fKinematicsMC2[0][11] = new TH1F("#phi_{e,prod.pl}","#phi_{e,prod.pl}", 100, 0, 3.5);
  fKinematicsMC1[0][11] = new TH1F("#phi_{e,prod.pl} Ratio","#phi_{e,prod.pl} Ratio", 100, 0, 3.5);
  for(int i=0; i<7; i++)
  {
    l2[0][i] = new TLine(0.1,0.4+i*0.2,100,0.4+i*0.2);
    l2[1][i] = new TLine(0.001,0.4+i*0.2,1,0.4+i*0.2);
    l2[2][i] = new TLine(0,0.4+i*0.2,1,0.4+i*0.2);
    l2[3][i] = new TLine(0,0.4+i*0.2,1,0.4+i*0.2);
    l2[4][i] = new TLine(2,0.4+i*0.2,18,0.4+i*0.2);
    l2[5][i] = new TLine(0,0.4+i*0.2,160,0.4+i*0.2);
    l2[6][i] = new TLine(140,0.4+i*0.2,180,0.4+i*0.2);
    l2[7][i] = new TLine(0,0.4+i*0.2,160,0.4+i*0.2);
    l2[8][i] = new TLine(0,0.4+i*0.2,0.05,0.4+i*0.2);
    l2[9][i] = new TLine(-1.7,0.4+i*0.2,1.7,0.4+i*0.2);
    l2[10][i] = new TLine(-320,0.4+i*0.2,-70,0.4+i*0.2);
    l2[11][i] = new TLine(0,0.4+i*0.2,3.5,0.4+i*0.2);
    l2[12][i] = new TLine(0,0.4+i*0.2,40,0.4+i*0.2);
    l2[13][i] = new TLine(0,0.4+i*0.2,0.25,0.4+i*0.2);
    l2[14][i] = new TLine(-3.5,0.4+i*0.2,3.5,0.4+i*0.2);
    l2[15][i] = new TLine(0,0.4+i*0.2,3.5,0.4+i*0.2);
    l2[16][i] = new TLine(0,0.4+i*0.2,3,0.4+i*0.2);
    for(int j=0; j<17; j++)
    {
      l2[j][i]->SetLineStyle(fLineStyle[i]);
      l2[j][i]->SetLineWidth(1);
    }
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
  c22.Divide(2,4);
  c23.Divide(2,4);
  c24.Divide(2,4);
  c25.Divide(1,2);
  c26.Divide(1,2);
  c27.Divide(1,2);
  c28.Divide(2,4);
  c29.Divide(1,2);
  c30.Divide(1,2);

  for(int i=0; i<8; i++)
  {
    int idx=int(i/2);
    if(i%2)
    {
      c1.cd(idx+3+int(idx/2)*2);
      for(int tt=0; tt<fKinematicsMC1[idx][0]->GetNbinsX(); tt++)
      {
        fError.push_back((fKinematicsMC1[idx][0]->GetBinError(tt) && fKinematicsMC2[idx][0]->GetBinError(tt) ? sqrt(pow(1/fKinematicsMC1[idx][0]->GetBinError(tt),2)+pow(1/fKinematicsMC2[idx][0]->GetBinError(tt),2)):0));
      }
      fKinematicsMC2[idx][0]->Scale(1/fKinematicsMC2[2][0]->GetEntries());
      fKinematicsMC1[idx][0]->Scale(1/fKinematicsMC1[2][0]->GetEntries());
      fKinematicsRatio[idx][0] = (TH1F*)fKinematicsMC1[idx][0]->Clone();
      fKinematicsRatio[idx][0]->SetStats(0);
      fKinematicsRatio[idx][0]->Divide(fKinematicsMC2[idx][0]);
      for(int tt=0; tt<fKinematicsRatio[idx][0]->GetNbinsX(); tt++)
      {
        fKinematicsRatio[idx][0]->SetBinError(tt,fError[tt]);
      }
      fError.clear();
      fKinematicsRatio[idx][0]->SetMarkerStyle(21);
      fKinematicsRatio[idx][0]->SetFillColor(kYellow-7);
      fKinematicsRatio[idx][0]->SetMaximum(2.);
      fKinematicsRatio[idx][0]->SetMinimum(0.);
      fKinematicsRatio[idx][0]->Draw("PE2");
      fKinematicsRatio[idx][0]->GetXaxis()->SetLabelSize(0.08);
      fKinematicsRatio[idx][0]->GetYaxis()->SetLabelSize(0.08);
      fKinematicsRatio[idx][0]->GetYaxis()->SetNdivisions(2,kTRUE);
      for(int tt=0; tt<7; tt++)
      {
        l2[0][tt]->Draw();
      }
      gPad->SetLogx();
      c1.Update();

      c2.cd(idx+3+int(idx/2)*2);
      for(int tt=0; tt<fKinematicsMC1[idx][1]->GetNbinsX(); tt++)
      {
        fError.push_back((fKinematicsMC1[idx][1]->GetBinError(tt) && fKinematicsMC2[idx][1]->GetBinError(tt) ? sqrt(pow(1/fKinematicsMC1[idx][1]->GetBinError(tt),2)+pow(1/fKinematicsMC2[idx][1]->GetBinError(tt),2)):0));
      }
      fKinematicsMC2[idx][1]->Scale(1/fKinematicsMC2[2][1]->GetEntries());
      fKinematicsMC1[idx][1]->Scale(1/fKinematicsMC1[2][1]->GetEntries());
      fKinematicsRatio[idx][1] = (TH1F*)fKinematicsMC1[idx][1]->Clone();
      fKinematicsRatio[idx][1]->SetStats(0);
      fKinematicsRatio[idx][1]->Divide(fKinematicsMC2[idx][1]);
      for(int tt=0; tt<fKinematicsRatio[idx][1]->GetNbinsX(); tt++)
      {
        fKinematicsRatio[idx][1]->SetBinError(tt,fError[tt]);
      }
      fError.clear();
      fKinematicsRatio[idx][1]->SetMarkerStyle(21);
      fKinematicsRatio[idx][1]->SetFillColor(kYellow-7);
      fKinematicsRatio[idx][1]->SetMaximum(2.);
      fKinematicsRatio[idx][1]->SetMinimum(0.);
      fKinematicsRatio[idx][1]->Draw("PE2");
      fKinematicsRatio[idx][1]->GetXaxis()->SetLabelSize(0.08);
      fKinematicsRatio[idx][1]->GetYaxis()->SetLabelSize(0.08);
      fKinematicsRatio[idx][1]->GetYaxis()->SetNdivisions(2,kTRUE);
      for(int tt=0; tt<7; tt++)
      {
        l2[1][tt]->Draw();
      }
      gPad->SetLogx();
      c2.Update();

      c3.cd(idx+3+int(idx/2)*2);
      for(int tt=0; tt<fKinematicsMC1[idx][2]->GetNbinsX(); tt++)
      {
        fError.push_back((fKinematicsMC1[idx][2]->GetBinError(tt) && fKinematicsMC2[idx][2]->GetBinError(tt) ? sqrt(pow(1/fKinematicsMC1[idx][2]->GetBinError(tt),2)+pow(1/fKinematicsMC2[idx][2]->GetBinError(tt),2)):0));
      }
      fKinematicsMC2[idx][2]->Scale(1/fKinematicsMC2[2][2]->GetEntries());
      fKinematicsMC1[idx][2]->Scale(1/fKinematicsMC1[2][2]->GetEntries());
      fKinematicsRatio[idx][2] = (TH1F*)fKinematicsMC1[idx][2]->Clone();
      fKinematicsRatio[idx][2]->SetStats(0);
      fKinematicsRatio[idx][2]->Divide(fKinematicsMC2[idx][2]);
      for(int tt=0; tt<fKinematicsRatio[idx][2]->GetNbinsX(); tt++)
      {
        fKinematicsRatio[idx][2]->SetBinError(tt,fError[tt]);
      }
      fError.clear();
      fKinematicsRatio[idx][2]->SetMarkerStyle(21);
      fKinematicsRatio[idx][2]->SetFillColor(kYellow-7);
      fKinematicsRatio[idx][2]->SetMaximum(2.);
      fKinematicsRatio[idx][2]->SetMinimum(0.);
      fKinematicsRatio[idx][2]->Draw("PE2");
      fKinematicsRatio[idx][2]->GetXaxis()->SetLabelSize(0.08);
      fKinematicsRatio[idx][2]->GetYaxis()->SetLabelSize(0.08);
      fKinematicsRatio[idx][2]->GetYaxis()->SetNdivisions(2,kTRUE);
      for(int tt=0; tt<7; tt++)
      {
        l2[2][tt]->Draw();
      }
      c3.Update();

      c4.cd(idx+3+int(idx/2)*2);
      for(int tt=0; tt<fKinematicsMC1[idx][3]->GetNbinsX(); tt++)
      {
        fError.push_back((fKinematicsMC1[idx][3]->GetBinError(tt) && fKinematicsMC2[idx][3]->GetBinError(tt) ? sqrt(pow(1/fKinematicsMC1[idx][3]->GetBinError(tt),2)+pow(1/fKinematicsMC2[idx][3]->GetBinError(tt),2)):0));
      }
      fKinematicsMC2[idx][3]->Scale(1/fKinematicsMC2[2][3]->GetEntries());
      fKinematicsMC1[idx][3]->Scale(1/fKinematicsMC1[2][3]->GetEntries());
      fKinematicsRatio[idx][3] = (TH1F*)fKinematicsMC1[idx][3]->Clone();
      fKinematicsRatio[idx][3]->SetStats(0);
      fKinematicsRatio[idx][3]->Divide(fKinematicsMC2[idx][3]);
      for(int tt=0; tt<fKinematicsRatio[idx][3]->GetNbinsX(); tt++)
      {
        fKinematicsRatio[idx][3]->SetBinError(tt,fError[tt]);
      }
      fError.clear();
      fKinematicsRatio[idx][3]->SetMarkerStyle(21);
      fKinematicsRatio[idx][3]->SetFillColor(kYellow-7);
      fKinematicsRatio[idx][3]->SetMaximum(2.);
      fKinematicsRatio[idx][3]->SetMinimum(0.);
      fKinematicsRatio[idx][3]->Draw("PE2");
      fKinematicsRatio[idx][3]->GetXaxis()->SetLabelSize(0.08);
      fKinematicsRatio[idx][3]->GetYaxis()->SetLabelSize(0.08);
      fKinematicsRatio[idx][3]->GetYaxis()->SetNdivisions(2,kTRUE);
      for(int tt=0; tt<7; tt++)
      {
        l2[3][tt]->Draw();
      }
      c4.Update();

      c5.cd(idx+3+int(idx/2)*2);
      for(int tt=0; tt<fKinematicsMC1[idx][4]->GetNbinsX(); tt++)
      {
        fError.push_back((fKinematicsMC1[idx][4]->GetBinError(tt) && fKinematicsMC2[idx][4]->GetBinError(tt) ? sqrt(pow(1/fKinematicsMC1[idx][4]->GetBinError(tt),2)+pow(1/fKinematicsMC2[idx][4]->GetBinError(tt),2)):0));
      }
      fKinematicsMC2[idx][4]->Scale(1/fKinematicsMC2[2][4]->GetEntries());
      fKinematicsMC1[idx][4]->Scale(1/fKinematicsMC1[2][4]->GetEntries());
      fKinematicsRatio[idx][4] = (TH1F*)fKinematicsMC1[idx][4]->Clone();
      fKinematicsRatio[idx][4]->SetStats(0);
      fKinematicsRatio[idx][4]->Divide(fKinematicsMC2[idx][4]);
      for(int tt=0; tt<fKinematicsRatio[idx][4]->GetNbinsX(); tt++)
      {
        fKinematicsRatio[idx][4]->SetBinError(tt,fError[tt]);
      }
      fError.clear();
      fKinematicsRatio[idx][4]->SetMarkerStyle(21);
      fKinematicsRatio[idx][4]->SetFillColor(kYellow-7);
      fKinematicsRatio[idx][4]->SetMaximum(2.);
      fKinematicsRatio[idx][4]->SetMinimum(0.);
      fKinematicsRatio[idx][4]->Draw("PE2");
      fKinematicsRatio[idx][4]->GetXaxis()->SetLabelSize(0.08);
      fKinematicsRatio[idx][4]->GetYaxis()->SetLabelSize(0.08);
      fKinematicsRatio[idx][4]->GetYaxis()->SetNdivisions(2,kTRUE);
      for(int tt=0; tt<7; tt++)
      {
        l2[4][tt]->Draw();
      }
      c5.Update();

      c6.cd(idx+3+int(idx/2)*2);
      for(int tt=0; tt<fKinematicsMC1[idx][5]->GetNbinsX(); tt++)
      {
        fError.push_back((fKinematicsMC1[idx][5]->GetBinError(tt) && fKinematicsMC2[idx][5]->GetBinError(tt) ? sqrt(pow(1/fKinematicsMC1[idx][5]->GetBinError(tt),2)+pow(1/fKinematicsMC2[idx][5]->GetBinError(tt),2)):0));
      }
      fKinematicsMC2[idx][5]->Scale(1/fKinematicsMC2[2][5]->GetEntries());
      fKinematicsMC1[idx][5]->Scale(1/fKinematicsMC1[2][5]->GetEntries());
      fKinematicsRatio[idx][5] = (TH1F*)fKinematicsMC1[idx][5]->Clone();
      fKinematicsRatio[idx][5]->SetStats(0);
      fKinematicsRatio[idx][5]->Divide(fKinematicsMC2[idx][5]);
      for(int tt=0; tt<fKinematicsRatio[idx][5]->GetNbinsX(); tt++)
      {
        fKinematicsRatio[idx][5]->SetBinError(tt,fError[tt]);
      }
      fError.clear();
      fKinematicsRatio[idx][5]->SetMarkerStyle(21);
      fKinematicsRatio[idx][5]->SetFillColor(kYellow-7);
      fKinematicsRatio[idx][5]->SetMaximum(2.);
      fKinematicsRatio[idx][5]->SetMinimum(0.);
      fKinematicsRatio[idx][5]->Draw("PE2");
      fKinematicsRatio[idx][5]->GetXaxis()->SetLabelSize(0.08);
      fKinematicsRatio[idx][5]->GetYaxis()->SetLabelSize(0.08);
      fKinematicsRatio[idx][5]->GetYaxis()->SetNdivisions(2,kTRUE);
      for(int tt=0; tt<7; tt++)
      {
        l2[5][tt]->Draw();
      }
      c6.Update();

      c14.cd(idx+3+int(idx/2)*2);
      for(int tt=0; tt<fKinematicsMC1[idx][6]->GetNbinsX(); tt++)
      {
        fError.push_back((fKinematicsMC1[idx][6]->GetBinError(tt) && fKinematicsMC2[idx][6]->GetBinError(tt) ? sqrt(pow(1/fKinematicsMC1[idx][6]->GetBinError(tt),2)+pow(1/fKinematicsMC2[idx][6]->GetBinError(tt),2)):0));
      }
      fKinematicsMC2[idx][6]->Scale(1/fKinematicsMC2[2][6]->GetEntries());
      fKinematicsMC1[idx][6]->Scale(1/fKinematicsMC1[2][6]->GetEntries());
      fKinematicsRatio[idx][6] = (TH1F*)fKinematicsMC1[idx][6]->Clone();
      fKinematicsRatio[idx][6]->SetStats(0);
      fKinematicsRatio[idx][6]->Divide(fKinematicsMC2[idx][6]);
      for(int tt=0; tt<fKinematicsRatio[idx][6]->GetNbinsX(); tt++)
      {
        fKinematicsRatio[idx][6]->SetBinError(tt,fError[tt]);
      }
      fError.clear();
      fKinematicsRatio[idx][6]->SetMarkerStyle(21);
      fKinematicsRatio[idx][6]->SetFillColor(kYellow-7);
      fKinematicsRatio[idx][6]->SetMaximum(2.);
      fKinematicsRatio[idx][6]->SetMinimum(0.);
      fKinematicsRatio[idx][6]->Draw("PE2");
      fKinematicsRatio[idx][6]->GetXaxis()->SetLabelSize(0.08);
      fKinematicsRatio[idx][6]->GetYaxis()->SetLabelSize(0.08);
      fKinematicsRatio[idx][6]->GetYaxis()->SetNdivisions(2,kTRUE);
      for(int tt=0; tt<7; tt++)
      {
        l2[6][tt]->Draw();
      }
      c14.Update();

      c15.cd(idx+3+int(idx/2)*2);
      for(int tt=0; tt<fKinematicsMC1[idx][7]->GetNbinsX(); tt++)
      {
        fError.push_back((fKinematicsMC1[idx][7]->GetBinError(tt) && fKinematicsMC2[idx][7]->GetBinError(tt) ? sqrt(pow(1/fKinematicsMC1[idx][7]->GetBinError(tt),2)+pow(1/fKinematicsMC2[idx][7]->GetBinError(tt),2)):0));
      }
      fKinematicsMC2[idx][7]->Scale(1/fKinematicsMC2[2][7]->GetEntries());
      fKinematicsMC1[idx][7]->Scale(1/fKinematicsMC1[2][7]->GetEntries());
      fKinematicsRatio[idx][7] = (TH1F*)fKinematicsMC1[idx][7]->Clone();
      fKinematicsRatio[idx][7]->SetStats(0);
      fKinematicsRatio[idx][7]->Divide(fKinematicsMC2[idx][7]);
      for(int tt=0; tt<fKinematicsRatio[idx][7]->GetNbinsX(); tt++)
      {
        fKinematicsRatio[idx][7]->SetBinError(tt,fError[tt]);
      }
      fError.clear();
      fKinematicsRatio[idx][7]->SetMarkerStyle(21);
      fKinematicsRatio[idx][7]->SetFillColor(kYellow-7);
      fKinematicsRatio[idx][7]->SetMaximum(2.);
      fKinematicsRatio[idx][7]->SetMinimum(0.);
      fKinematicsRatio[idx][7]->Draw("PE2");
      fKinematicsRatio[idx][7]->GetXaxis()->SetLabelSize(0.08);
      fKinematicsRatio[idx][7]->GetYaxis()->SetLabelSize(0.08);
      fKinematicsRatio[idx][7]->GetYaxis()->SetNdivisions(2,kTRUE);
      for(int tt=0; tt<7; tt++)
      {
        l2[7][tt]->Draw();
      }
      c15.Update();

      c16.cd(idx+3+int(idx/2)*2);
      for(int tt=0; tt<fKinematicsMC1[idx][8]->GetNbinsX(); tt++)
      {
        fError.push_back((fKinematicsMC1[idx][8]->GetBinError(tt) && fKinematicsMC2[idx][8]->GetBinError(tt) ? sqrt(pow(1/fKinematicsMC1[idx][8]->GetBinError(tt),2)+pow(1/fKinematicsMC2[idx][8]->GetBinError(tt),2)):0));
      }
      fKinematicsMC2[idx][8]->Scale(1/fKinematicsMC2[2][8]->GetEntries());
      fKinematicsMC1[idx][8]->Scale(1/fKinematicsMC1[2][8]->GetEntries());
      fKinematicsRatio[idx][8] = (TH1F*)fKinematicsMC1[idx][8]->Clone();
      fKinematicsRatio[idx][8]->SetStats(0);
      fKinematicsRatio[idx][8]->Divide(fKinematicsMC2[idx][8]);
      for(int tt=0; tt<fKinematicsRatio[idx][8]->GetNbinsX(); tt++)
      {
        fKinematicsRatio[idx][8]->SetBinError(tt,fError[tt]);
      }
      fError.clear();
      fKinematicsRatio[idx][8]->SetMarkerStyle(21);
      fKinematicsRatio[idx][8]->SetFillColor(kYellow-7);
      fKinematicsRatio[idx][8]->SetMaximum(2.);
      fKinematicsRatio[idx][8]->SetMinimum(0.);
      fKinematicsRatio[idx][8]->Draw("PE2");
      fKinematicsRatio[idx][8]->GetXaxis()->SetLabelSize(0.08);
      fKinematicsRatio[idx][8]->GetYaxis()->SetLabelSize(0.08);
      fKinematicsRatio[idx][8]->GetYaxis()->SetNdivisions(2,kTRUE);
      for(int tt=0; tt<7; tt++)
      {
        l2[8][tt]->Draw();
      }
      c16.Update();

      c17.cd(idx+3+int(idx/2)*2);
      for(int tt=0; tt<fKinematicsMC1[idx][9]->GetNbinsX(); tt++)
      {
        fError.push_back((fKinematicsMC1[idx][9]->GetBinError(tt) && fKinematicsMC2[idx][9]->GetBinError(tt) ? sqrt(pow(1/fKinematicsMC1[idx][9]->GetBinError(tt),2)+pow(1/fKinematicsMC2[idx][9]->GetBinError(tt),2)):0));
      }
      fKinematicsMC2[idx][9]->Scale(1/fKinematicsMC2[2][9]->GetEntries());
      fKinematicsMC1[idx][9]->Scale(1/fKinematicsMC1[2][9]->GetEntries());
      fKinematicsRatio[idx][9] = (TH1F*)fKinematicsMC1[idx][9]->Clone();
      fKinematicsRatio[idx][9]->SetStats(0);
      fKinematicsRatio[idx][9]->Divide(fKinematicsMC2[idx][9]);
      for(int tt=0; tt<fKinematicsRatio[idx][9]->GetNbinsX(); tt++)
      {
        fKinematicsRatio[idx][9]->SetBinError(tt,fError[tt]);
      }
      fError.clear();
      fKinematicsRatio[idx][9]->SetMarkerStyle(21);
      fKinematicsRatio[idx][9]->SetFillColor(kYellow-7);
      fKinematicsRatio[idx][9]->SetMaximum(2.);
      fKinematicsRatio[idx][9]->SetMinimum(0.);
      fKinematicsRatio[idx][9]->Draw("PE2");
      fKinematicsRatio[idx][9]->GetXaxis()->SetLabelSize(0.08);
      fKinematicsRatio[idx][9]->GetYaxis()->SetLabelSize(0.08);
      fKinematicsRatio[idx][9]->GetYaxis()->SetNdivisions(2,kTRUE);
      for(int tt=0; tt<7; tt++)
      {
        l2[9][tt]->Draw();
      }
      c17.Update();

      c18.cd(idx+3+int(idx/2)*2);
      for(int tt=0; tt<fKinematicsMC1[idx][10]->GetNbinsX(); tt++)
      {
        fError.push_back((fKinematicsMC1[idx][10]->GetBinError(tt) && fKinematicsMC2[idx][10]->GetBinError(tt) ? sqrt(pow(1/fKinematicsMC1[idx][10]->GetBinError(tt),2)+pow(1/fKinematicsMC2[idx][10]->GetBinError(tt),2)):0));
      }
      fKinematicsMC2[idx][10]->Scale(1/fKinematicsMC2[2][10]->GetEntries());
      fKinematicsMC1[idx][10]->Scale(1/fKinematicsMC1[2][10]->GetEntries());
      fKinematicsRatio[idx][10] = (TH1F*)fKinematicsMC1[idx][10]->Clone();
      fKinematicsRatio[idx][10]->SetStats(0);
      fKinematicsRatio[idx][10]->Divide(fKinematicsMC2[idx][10]);
      for(int tt=0; tt<fKinematicsRatio[idx][10]->GetNbinsX(); tt++)
      {
        fKinematicsRatio[idx][10]->SetBinError(tt,fError[tt]);
      }
      fError.clear();
      fKinematicsRatio[idx][10]->SetMarkerStyle(21);
      fKinematicsRatio[idx][10]->SetFillColor(kYellow-7);
      fKinematicsRatio[idx][10]->SetMaximum(2.);
      fKinematicsRatio[idx][10]->SetMinimum(0.);
      fKinematicsRatio[idx][10]->Draw("PE2");
      fKinematicsRatio[idx][10]->GetXaxis()->SetLabelSize(0.08);
      fKinematicsRatio[idx][10]->GetYaxis()->SetLabelSize(0.08);
      fKinematicsRatio[idx][10]->GetYaxis()->SetNdivisions(2,kTRUE);
      for(int tt=0; tt<7; tt++)
      {
        l2[10][tt]->Draw();
      }
      c18.Update();

      c22.cd(idx+3+int(idx/2)*2);
      for(int tt=0; tt<fKinematicsMC1[idx][12]->GetNbinsX(); tt++)
      {
        fError.push_back((fKinematicsMC1[idx][12]->GetBinError(tt) && fKinematicsMC2[idx][12]->GetBinError(tt) ? sqrt(pow(1/fKinematicsMC1[idx][12]->GetBinError(tt),2)+pow(1/fKinematicsMC2[idx][12]->GetBinError(tt),2)):0));
      }
      fKinematicsMC2[idx][12]->Scale(1/fKinematicsMC2[2][12]->GetEntries());
      fKinematicsMC1[idx][12]->Scale(1/fKinematicsMC1[2][12]->GetEntries());
      fKinematicsRatio[idx][12] = (TH1F*)fKinematicsMC1[idx][12]->Clone();
      fKinematicsRatio[idx][12]->SetStats(0);
      fKinematicsRatio[idx][12]->Divide(fKinematicsMC2[idx][12]);
      for(int tt=0; tt<fKinematicsRatio[idx][12]->GetNbinsX(); tt++)
      {
        fKinematicsRatio[idx][12]->SetBinError(tt,fError[tt]);
      }
      fError.clear();
      fKinematicsRatio[idx][12]->SetMarkerStyle(21);
      fKinematicsRatio[idx][12]->SetFillColor(kYellow-7);
      fKinematicsRatio[idx][12]->SetMaximum(2.);
      fKinematicsRatio[idx][12]->SetMinimum(0.);
      fKinematicsRatio[idx][12]->Draw("PE2");
      fKinematicsRatio[idx][12]->GetXaxis()->SetLabelSize(0.08);
      fKinematicsRatio[idx][12]->GetYaxis()->SetLabelSize(0.08);
      fKinematicsRatio[idx][12]->GetYaxis()->SetNdivisions(2,kTRUE);
      for(int tt=0; tt<7; tt++)
      {
        l2[12][tt]->Draw();
      }
      c22.Update();

      c23.cd(idx+3+int(idx/2)*2);
      for(int tt=0; tt<fKinematicsMC1[idx][13]->GetNbinsX(); tt++)
      {
        fError.push_back((fKinematicsMC1[idx][13]->GetBinError(tt) && fKinematicsMC2[idx][13]->GetBinError(tt) ? sqrt(pow(1/fKinematicsMC1[idx][13]->GetBinError(tt),2)+pow(1/fKinematicsMC2[idx][13]->GetBinError(tt),2)):0));
      }
      fKinematicsMC2[idx][13]->Scale(1/fKinematicsMC2[2][13]->GetEntries());
      fKinematicsMC1[idx][13]->Scale(1/fKinematicsMC1[2][13]->GetEntries());
      fKinematicsRatio[idx][13] = (TH1F*)fKinematicsMC1[idx][13]->Clone();
      fKinematicsRatio[idx][13]->SetStats(0);
      fKinematicsRatio[idx][13]->Divide(fKinematicsMC2[idx][13]);
      for(int tt=0; tt<fKinematicsRatio[idx][13]->GetNbinsX(); tt++)
      {
        fKinematicsRatio[idx][13]->SetBinError(tt,fError[tt]);
      }
      fError.clear();
      fKinematicsRatio[idx][13]->SetMarkerStyle(21);
      fKinematicsRatio[idx][13]->SetFillColor(kYellow-7);
      fKinematicsRatio[idx][13]->SetMaximum(2.);
      fKinematicsRatio[idx][13]->SetMinimum(0.);
      fKinematicsRatio[idx][13]->Draw("PE2");
      fKinematicsRatio[idx][13]->GetXaxis()->SetLabelSize(0.08);
      fKinematicsRatio[idx][13]->GetYaxis()->SetLabelSize(0.08);
      fKinematicsRatio[idx][13]->GetYaxis()->SetNdivisions(2,kTRUE);
      for(int tt=0; tt<7; tt++)
      {
        l2[13][tt]->Draw();
      }
      c23.Update();

      c24.cd(idx+3+int(idx/2)*2);
      for(int tt=0; tt<fKinematicsMC1[idx][14]->GetNbinsX(); tt++)
      {
        fError.push_back((fKinematicsMC1[idx][14]->GetBinError(tt) && fKinematicsMC2[idx][14]->GetBinError(tt) ? sqrt(pow(1/fKinematicsMC1[idx][14]->GetBinError(tt),2)+pow(1/fKinematicsMC2[idx][14]->GetBinError(tt),2)):0));
      }
      fKinematicsMC2[idx][14]->Scale(1/fKinematicsMC2[2][14]->GetEntries());
      fKinematicsMC1[idx][14]->Scale(1/fKinematicsMC1[2][14]->GetEntries());
      fKinematicsRatio[idx][14] = (TH1F*)fKinematicsMC1[idx][14]->Clone();
      fKinematicsRatio[idx][14]->SetStats(0);
      fKinematicsRatio[idx][14]->Divide(fKinematicsMC2[idx][14]);
      for(int tt=0; tt<fKinematicsRatio[idx][14]->GetNbinsX(); tt++)
      {
        fKinematicsRatio[idx][14]->SetBinError(tt,fError[tt]);
      }
      fError.clear();
      fKinematicsRatio[idx][14]->SetMarkerStyle(21);
      fKinematicsRatio[idx][14]->SetFillColor(kYellow-7);
      fKinematicsRatio[idx][14]->SetMaximum(2.);
      fKinematicsRatio[idx][14]->SetMinimum(0.);
      fKinematicsRatio[idx][14]->Draw("PE2");
      fKinematicsRatio[idx][14]->GetXaxis()->SetLabelSize(0.08);
      fKinematicsRatio[idx][14]->GetYaxis()->SetLabelSize(0.08);
      fKinematicsRatio[idx][14]->GetYaxis()->SetNdivisions(2,kTRUE);
      for(int tt=0; tt<7; tt++)
      {
        l2[14][tt]->Draw();
      }
      c24.Update();

      c28.cd(idx+3+int(idx/2)*2);
      for(int tt=0; tt<fKinematicsMC1[idx][15]->GetNbinsX(); tt++)
      {
        fError.push_back((fKinematicsMC1[idx][15]->GetBinError(tt) && fKinematicsMC2[idx][15]->GetBinError(tt) ? sqrt(pow(1/fKinematicsMC1[idx][15]->GetBinError(tt),2)+pow(1/fKinematicsMC2[idx][15]->GetBinError(tt),2)):0));
      }
      fKinematicsMC2[idx][15]->Scale(1/fKinematicsMC2[2][15]->GetEntries());
      fKinematicsMC1[idx][15]->Scale(1/fKinematicsMC1[2][15]->GetEntries());
      fKinematicsRatio[idx][15] = (TH1F*)fKinematicsMC1[idx][15]->Clone();
      fKinematicsRatio[idx][15]->SetStats(0);
      fKinematicsRatio[idx][15]->Divide(fKinematicsMC2[idx][15]);
      for(int tt=0; tt<fKinematicsRatio[idx][15]->GetNbinsX(); tt++)
      {
        fKinematicsRatio[idx][15]->SetBinError(tt,fError[tt]);
      }
      fError.clear();
      fKinematicsRatio[idx][15]->SetMarkerStyle(21);
      fKinematicsRatio[idx][15]->SetFillColor(kYellow-7);
      fKinematicsRatio[idx][15]->SetMaximum(2.);
      fKinematicsRatio[idx][15]->SetMinimum(0.);
      fKinematicsRatio[idx][15]->Draw("PE2");
      fKinematicsRatio[idx][15]->GetXaxis()->SetLabelSize(0.08);
      fKinematicsRatio[idx][15]->GetYaxis()->SetLabelSize(0.08);
      fKinematicsRatio[idx][15]->GetYaxis()->SetNdivisions(2,kTRUE);
      for(int tt=0; tt<7; tt++)
      {
        l2[15][tt]->Draw();
      }
      c28.Update();
    }
    else
    {
      c1.cd(idx+1+int(idx/2)*2);
      fKinematicsMC2[idx][0]->SetLineColor(kGreen);
      fKinematicsMC2[idx][0]->SetStats(0);
      fKinematicsMC2[idx][0]->SetMinimum(0.);
      // fKinematicsMC2[idx][0]->SetMaximum(0.08);
      fKinematicsMC2[idx][0]->Draw();
      fKinematicsMC2[idx][0]->GetXaxis()->SetLabelSize(0.08);
      fKinematicsMC2[idx][0]->GetYaxis()->SetLabelSize(0.08);
      fKinematicsMC2[idx][0]->GetYaxis()->SetNdivisions(304,kTRUE);
      fKinematicsMC1[idx][0]->SetLineColor(kBlue);
      fKinematicsMC1[idx][0]->Draw("SAME");
      gPad->SetLogx();
      c1.Update();

      c2.cd(idx+1+int(idx/2)*2);
      fKinematicsMC2[idx][1]->SetLineColor(kGreen);
      fKinematicsMC2[idx][1]->SetStats(0);
      fKinematicsMC2[idx][1]->SetMinimum(0.);
      // fKinematicsMC2[idx][1]->SetMaximum(0.045);
      fKinematicsMC2[idx][1]->Draw();
      fKinematicsMC2[idx][1]->GetXaxis()->SetLabelSize(0.08);
      fKinematicsMC2[idx][1]->GetYaxis()->SetLabelSize(0.08);
      fKinematicsMC2[idx][1]->GetYaxis()->SetNdivisions(304,kTRUE);
      fKinematicsMC1[idx][1]->SetLineColor(kBlue);
      fKinematicsMC1[idx][1]->Draw("SAME");
      gPad->SetLogx();
      c2.Update();

      c3.cd(idx+1+int(idx/2)*2);
      fKinematicsMC2[idx][2]->SetLineColor(kGreen);
      fKinematicsMC2[idx][2]->SetStats(0);
      fKinematicsMC2[idx][2]->SetMinimum(0.);
      // fKinematicsMC2[idx][2]->SetMaximum(0.08);
      fKinematicsMC2[idx][2]->Draw();
      fKinematicsMC2[idx][2]->GetXaxis()->SetLabelSize(0.08);
      fKinematicsMC2[idx][2]->GetYaxis()->SetLabelSize(0.08);
      fKinematicsMC2[idx][2]->GetYaxis()->SetNdivisions(304,kTRUE);
      fKinematicsMC1[idx][2]->SetLineColor(kBlue);
      fKinematicsMC1[idx][2]->Draw("SAME");
      c3.Update();

      c4.cd(idx+1+int(idx/2)*2);
      fKinematicsMC2[idx][3]->SetLineColor(kGreen);
      fKinematicsMC2[idx][3]->SetStats(0);
      fKinematicsMC2[idx][3]->SetMinimum(0.);
      fKinematicsMC2[idx][3]->Draw();
      fKinematicsMC2[idx][3]->GetXaxis()->SetLabelSize(0.08);
      fKinematicsMC2[idx][3]->GetYaxis()->SetLabelSize(0.08);
      fKinematicsMC2[idx][3]->GetYaxis()->SetNdivisions(304,kTRUE);
      fKinematicsMC1[idx][3]->SetLineColor(kBlue);
      fKinematicsMC1[idx][3]->Draw("SAME");
      c4.Update();

      c5.cd(idx+1+int(idx/2)*2);
      fKinematicsMC2[idx][4]->SetLineColor(kGreen);
      fKinematicsMC2[idx][4]->SetStats(0);
      fKinematicsMC2[idx][4]->SetMinimum(0.);
      // fKinematicsMC2[idx][4]->SetMaximum(0.05);
      fKinematicsMC2[idx][4]->Draw();
      fKinematicsMC2[idx][4]->GetXaxis()->SetLabelSize(0.08);
      fKinematicsMC2[idx][4]->GetYaxis()->SetLabelSize(0.08);
      fKinematicsMC2[idx][4]->GetYaxis()->SetNdivisions(304,kTRUE);
      fKinematicsMC1[idx][4]->SetLineColor(kBlue);
      fKinematicsMC1[idx][4]->Draw("SAME");
      c5.Update();

      c6.cd(idx+1+int(idx/2)*2);
      fKinematicsMC2[idx][5]->SetLineColor(kGreen);
      fKinematicsMC2[idx][5]->SetStats(0);
      fKinematicsMC2[idx][5]->SetMinimum(0.);
      // fKinematicsMC2[idx][5]->SetMaximum(0.08);
      fKinematicsMC2[idx][5]->Draw();
      fKinematicsMC2[idx][5]->GetXaxis()->SetLabelSize(0.08);
      fKinematicsMC2[idx][5]->GetYaxis()->SetLabelSize(0.08);
      fKinematicsMC2[idx][5]->GetYaxis()->SetNdivisions(304,kTRUE);
      fKinematicsMC1[idx][5]->SetLineColor(kBlue);
      fKinematicsMC1[idx][5]->Draw("SAME");
      c6.Update();

      c14.cd(idx+1+int(idx/2)*2);
      fKinematicsMC2[idx][6]->SetLineColor(kGreen);
      fKinematicsMC2[idx][6]->SetStats(0);
      fKinematicsMC2[idx][6]->SetMinimum(0.);
      // fKinematicsMC2[idx][5]->SetMaximum(0.08);
      fKinematicsMC2[idx][6]->Draw();
      fKinematicsMC2[idx][6]->GetXaxis()->SetLabelSize(0.08);
      fKinematicsMC2[idx][6]->GetYaxis()->SetLabelSize(0.08);
      fKinematicsMC2[idx][6]->GetYaxis()->SetNdivisions(304,kTRUE);
      fKinematicsMC1[idx][6]->SetLineColor(kBlue);
      fKinematicsMC1[idx][6]->Draw("SAME");
      c14.Update();

      c15.cd(idx+1+int(idx/2)*2);
      fKinematicsMC2[idx][7]->SetLineColor(kGreen);
      fKinematicsMC2[idx][7]->SetStats(0);
      fKinematicsMC2[idx][7]->SetMinimum(0.);
      // fKinematicsMC2[idx][5]->SetMaximum(0.08);
      fKinematicsMC2[idx][7]->Draw();
      fKinematicsMC2[idx][7]->GetXaxis()->SetLabelSize(0.08);
      fKinematicsMC2[idx][7]->GetYaxis()->SetLabelSize(0.08);
      fKinematicsMC2[idx][7]->GetYaxis()->SetNdivisions(304,kTRUE);
      fKinematicsMC1[idx][7]->SetLineColor(kBlue);
      fKinematicsMC1[idx][7]->Draw("SAME");
      c15.Update();

      c16.cd(idx+1+int(idx/2)*2);
      fKinematicsMC2[idx][8]->SetLineColor(kGreen);
      fKinematicsMC2[idx][8]->SetStats(0);
      fKinematicsMC2[idx][8]->SetMinimum(0.);
      // fKinematicsMC2[idx][5]->SetMaximum(0.08);
      fKinematicsMC2[idx][8]->Draw();
      fKinematicsMC2[idx][8]->GetXaxis()->SetLabelSize(0.08);
      fKinematicsMC2[idx][8]->GetYaxis()->SetLabelSize(0.08);
      fKinematicsMC2[idx][8]->GetYaxis()->SetNdivisions(304,kTRUE);
      fKinematicsMC1[idx][8]->SetLineColor(kBlue);
      fKinematicsMC1[idx][8]->Draw("SAME");
      c16.Update();

      c17.cd(idx+1+int(idx/2)*2);
      fKinematicsMC2[idx][9]->SetLineColor(kGreen);
      fKinematicsMC2[idx][9]->SetStats(0);
      fKinematicsMC2[idx][9]->SetMinimum(0.);
      // fKinematicsMC2[idx][5]->SetMaximum(0.08);
      fKinematicsMC2[idx][9]->Draw();
      fKinematicsMC2[idx][9]->GetXaxis()->SetLabelSize(0.08);
      fKinematicsMC2[idx][9]->GetYaxis()->SetLabelSize(0.08);
      fKinematicsMC2[idx][9]->GetYaxis()->SetNdivisions(304,kTRUE);
      fKinematicsMC1[idx][9]->SetLineColor(kBlue);
      fKinematicsMC1[idx][9]->Draw("SAME");
      c17.Update();

      c18.cd(idx+1+int(idx/2)*2);
      fKinematicsMC2[idx][10]->SetLineColor(kGreen);
      fKinematicsMC2[idx][10]->SetStats(0);
      fKinematicsMC2[idx][10]->SetMinimum(0.);
      // fKinematicsMC2[idx][5]->SetMaximum(0.08);
      fKinematicsMC2[idx][10]->Draw();
      fKinematicsMC2[idx][10]->GetXaxis()->SetLabelSize(0.08);
      fKinematicsMC2[idx][10]->GetYaxis()->SetLabelSize(0.08);
      fKinematicsMC2[idx][10]->GetYaxis()->SetNdivisions(304,kTRUE);
      fKinematicsMC1[idx][10]->SetLineColor(kBlue);
      fKinematicsMC1[idx][10]->Draw("SAME");
      c18.Update();

      c22.cd(idx+1+int(idx/2)*2);
      fKinematicsMC2[idx][12]->SetLineColor(kGreen);
      fKinematicsMC2[idx][12]->SetStats(0);
      fKinematicsMC2[idx][12]->SetMinimum(0.);
      // fKinematicsMC2[idx][5]->SetMaximum(0.08);
      fKinematicsMC2[idx][12]->Draw();
      fKinematicsMC2[idx][12]->GetXaxis()->SetLabelSize(0.08);
      fKinematicsMC2[idx][12]->GetYaxis()->SetLabelSize(0.08);
      fKinematicsMC2[idx][12]->GetYaxis()->SetNdivisions(304,kTRUE);
      fKinematicsMC1[idx][12]->SetLineColor(kBlue);
      fKinematicsMC1[idx][12]->Draw("SAME");
      c22.Update();

      c23.cd(idx+1+int(idx/2)*2);
      fKinematicsMC2[idx][13]->SetLineColor(kGreen);
      fKinematicsMC2[idx][13]->SetStats(0);
      fKinematicsMC2[idx][13]->SetMinimum(0.);
      // fKinematicsMC2[idx][5]->SetMaximum(0.08);
      fKinematicsMC2[idx][13]->Draw();
      fKinematicsMC2[idx][13]->GetXaxis()->SetLabelSize(0.08);
      fKinematicsMC2[idx][13]->GetYaxis()->SetLabelSize(0.08);
      fKinematicsMC2[idx][13]->GetYaxis()->SetNdivisions(304,kTRUE);
      fKinematicsMC1[idx][13]->SetLineColor(kBlue);
      fKinematicsMC1[idx][13]->Draw("SAME");
      c23.Update();

      c24.cd(idx+1+int(idx/2)*2);
      fKinematicsMC2[idx][14]->SetLineColor(kGreen);
      fKinematicsMC2[idx][14]->SetStats(0);
      fKinematicsMC2[idx][14]->SetMinimum(0.);
      // fKinematicsMC2[idx][5]->SetMaximum(0.08);
      fKinematicsMC2[idx][14]->Draw();
      fKinematicsMC2[idx][14]->GetXaxis()->SetLabelSize(0.08);
      fKinematicsMC2[idx][14]->GetYaxis()->SetLabelSize(0.08);
      fKinematicsMC2[idx][14]->GetYaxis()->SetNdivisions(304,kTRUE);
      fKinematicsMC1[idx][14]->SetLineColor(kBlue);
      fKinematicsMC1[idx][14]->Draw("SAME");
      c24.Update();

      c28.cd(idx+1+int(idx/2)*2);
      fKinematicsMC2[idx][15]->SetLineColor(kGreen);
      fKinematicsMC2[idx][15]->SetStats(0);
      fKinematicsMC2[idx][15]->SetMinimum(0.);
      // fKinematicsMC2[idx][5]->SetMaximum(0.08);
      fKinematicsMC2[idx][15]->Draw();
      fKinematicsMC2[idx][15]->GetXaxis()->SetLabelSize(0.08);
      fKinematicsMC2[idx][15]->GetYaxis()->SetLabelSize(0.08);
      fKinematicsMC2[idx][15]->GetYaxis()->SetNdivisions(304,kTRUE);
      fKinematicsMC1[idx][15]->SetLineColor(kBlue);
      fKinematicsMC1[idx][15]->Draw("SAME");
      c28.Update();
    }
  }

  c7.cd(2);
  for(int tt=0; tt<fKinematicsMC1[0][11]->GetNbinsX(); tt++)
  {
    fError.push_back((fKinematicsMC1[0][11]->GetBinError(tt) && fKinematicsMC2[0][11]->GetBinError(tt) ? sqrt(pow(1/fKinematicsMC1[0][11]->GetBinError(tt),2)+pow(1/fKinematicsMC2[0][11]->GetBinError(tt),2)):0));
  }
  fKinematicsMC2[0][11]->Scale(1/fKinematicsMC2[0][11]->GetEntries());
  fKinematicsMC1[0][11]->Scale(1/fKinematicsMC1[0][11]->GetEntries());
  fKinematicsRatio[0][11] = (TH1F*)fKinematicsMC1[0][11]->Clone();
  fKinematicsRatio[0][11]->SetStats(0);
  fKinematicsRatio[0][11]->Divide(fKinematicsMC2[0][11]);
  for(int tt=0; tt<fKinematicsRatio[0][11]->GetNbinsX(); tt++)
  {
    fKinematicsRatio[0][11]->SetBinError(tt,fError[tt]);
  }
  fError.clear();
  fKinematicsRatio[0][11]->SetMarkerStyle(21);
  fKinematicsRatio[0][11]->SetMaximum(2.);
  fKinematicsRatio[0][11]->SetFillColor(kYellow-7);
  fKinematicsRatio[0][11]->Draw("PE2");
  fKinematicsRatio[0][11]->GetXaxis()->SetLabelSize(0.08);
  fKinematicsRatio[0][11]->GetYaxis()->SetLabelSize(0.08);
  fKinematicsRatio[0][11]->GetYaxis()->SetNdivisions(2,kTRUE);
  for(int tt=0; tt<7; tt++)
  {
    l2[11][tt]->Draw();
  }
  c7.Update();
  c7.cd(1);
  fKinematicsMC2[0][11]->SetLineColor(kGreen);
  fKinematicsMC2[0][11]->SetStats(0);
  fKinematicsMC2[0][11]->SetMinimum(0.);
  fKinematicsMC2[0][11]->SetMaximum(max(fKinematicsMC2[0][11]->GetMaximum()*1.2,fKinematicsMC1[0][11]->GetMaximum()*1.2));
  fKinematicsMC2[0][11]->Draw();
  fKinematicsMC2[0][11]->GetXaxis()->SetLabelSize(0.08);
  fKinematicsMC2[0][11]->GetYaxis()->SetLabelSize(0.08);
  fKinematicsMC2[0][11]->GetYaxis()->SetNdivisions(304,kTRUE);
  fKinematicsMC1[0][11]->SetLineColor(kBlue);
  fKinematicsMC1[0][11]->Draw("SAME");
  c7.Update();

  c8.cd(2);
  for(int tt=0; tt<fKinematicsMC1[4][0]->GetNbinsX(); tt++)
  {
    fError.push_back((fKinematicsMC1[4][0]->GetBinError(tt) && fKinematicsMC2[4][0]->GetBinError(tt) ? sqrt(pow(1/fKinematicsMC1[4][0]->GetBinError(tt),2)+pow(1/fKinematicsMC2[4][0]->GetBinError(tt),2)):0));
  }
  fKinematicsMC2[4][0]->Scale(1/fKinematicsMC2[4][0]->GetEntries());
  fKinematicsMC1[4][0]->Scale(1/fKinematicsMC1[4][0]->GetEntries());
  fKinematicsRatio[4][0] = (TH1F*)fKinematicsMC1[4][0]->Clone();
  fKinematicsRatio[4][0]->SetStats(0);
  fKinematicsRatio[4][0]->Divide(fKinematicsMC2[4][0]);
  for(int tt=0; tt<fKinematicsRatio[4][0]->GetNbinsX(); tt++)
  {
    fKinematicsRatio[4][0]->SetBinError(tt,fError[tt]);
  }
  fError.clear();
  fKinematicsRatio[4][0]->SetMarkerStyle(21);
  fKinematicsRatio[4][0]->SetMaximum(2.);
  fKinematicsRatio[4][0]->SetMinimum(0.);
  fKinematicsRatio[4][0]->SetFillColor(kYellow-7);
  fKinematicsRatio[4][0]->Draw("PE2");
  fKinematicsRatio[4][0]->GetXaxis()->SetLabelSize(0.08);
  fKinematicsRatio[4][0]->GetYaxis()->SetLabelSize(0.08);
  fKinematicsRatio[4][0]->GetYaxis()->SetNdivisions(2,kTRUE);
  for(int tt=0; tt<7; tt++)
  {
    l2[0][tt]->Draw();
  }
  gPad->SetLogx();
  c8.Update();
  c8.cd(1);
  fKinematicsMC2[4][0]->SetLineColor(kGreen);
  fKinematicsMC2[4][0]->SetStats(0);
  fKinematicsMC2[4][0]->Draw();
  fKinematicsMC2[4][0]->GetXaxis()->SetLabelSize(0.08);
  fKinematicsMC2[4][0]->GetYaxis()->SetLabelSize(0.08);
  fKinematicsMC2[4][0]->GetYaxis()->SetNdivisions(304,kTRUE);
  fKinematicsMC1[4][0]->SetLineColor(kBlue);
  fKinematicsMC1[4][0]->Draw("SAME");
  gPad->SetLogx();
  c8.Update();

  c9.cd(2);
  for(int tt=0; tt<fKinematicsMC1[4][1]->GetNbinsX(); tt++)
  {
    fError.push_back((fKinematicsMC1[4][1]->GetBinError(tt) && fKinematicsMC2[4][1]->GetBinError(tt) ? sqrt(pow(1/fKinematicsMC1[4][1]->GetBinError(tt),2)+pow(1/fKinematicsMC2[4][1]->GetBinError(tt),2)):0));
  }
  fKinematicsMC2[4][1]->Scale(1/fKinematicsMC2[4][1]->GetEntries());
  fKinematicsMC1[4][1]->Scale(1/fKinematicsMC1[4][1]->GetEntries());
  fKinematicsRatio[4][1] = (TH1F*)fKinematicsMC1[4][1]->Clone();
  fKinematicsRatio[4][1]->SetStats(0);
  fKinematicsRatio[4][1]->Divide(fKinematicsMC2[4][1]);
  for(int tt=0; tt<fKinematicsRatio[4][1]->GetNbinsX(); tt++)
  {
    fKinematicsRatio[4][1]->SetBinError(tt,fError[tt]);
  }
  fError.clear();
  fKinematicsRatio[4][1]->SetMarkerStyle(21);
  fKinematicsRatio[4][1]->SetMaximum(2.);
  fKinematicsRatio[4][1]->SetMinimum(0.);
  fKinematicsRatio[4][1]->SetFillColor(kYellow-7);
  fKinematicsRatio[4][1]->Draw("PE2");
  fKinematicsRatio[4][1]->GetXaxis()->SetLabelSize(0.08);
  fKinematicsRatio[4][1]->GetYaxis()->SetLabelSize(0.08);
  fKinematicsRatio[4][1]->GetYaxis()->SetNdivisions(2,kTRUE);
  for(int tt=0; tt<7; tt++)
  {
    l2[1][tt]->Draw();
  }
  gPad->SetLogx();
  c9.Update();
  c9.cd(1);
  fKinematicsMC2[4][1]->SetLineColor(kGreen);
  fKinematicsMC2[4][1]->SetStats(0);
  fKinematicsMC2[4][1]->Draw();
  fKinematicsMC2[4][1]->GetXaxis()->SetLabelSize(0.08);
  fKinematicsMC2[4][1]->GetYaxis()->SetLabelSize(0.08);
  fKinematicsMC2[4][1]->GetYaxis()->SetNdivisions(304,kTRUE);
  fKinematicsMC1[4][1]->SetLineColor(kBlue);
  fKinematicsMC1[4][1]->Draw("SAME");
  gPad->SetLogx();
  c9.Update();

  c10.cd(2);
  for(int tt=0; tt<fKinematicsMC1[4][2]->GetNbinsX(); tt++)
  {
    fError.push_back((fKinematicsMC1[4][2]->GetBinError(tt) && fKinematicsMC2[4][2]->GetBinError(tt) ? sqrt(pow(1/fKinematicsMC1[4][2]->GetBinError(tt),2)+pow(1/fKinematicsMC2[4][2]->GetBinError(tt),2)):0));
  }
  fKinematicsMC2[4][2]->Scale(1/fKinematicsMC2[4][2]->GetEntries());
  fKinematicsMC1[4][2]->Scale(1/fKinematicsMC1[4][2]->GetEntries());
  fKinematicsRatio[4][2] = (TH1F*)fKinematicsMC1[4][2]->Clone();
  fKinematicsRatio[4][2]->SetStats(0);
  fKinematicsRatio[4][2]->Divide(fKinematicsMC2[4][2]);
  for(int tt=0; tt<fKinematicsRatio[4][2]->GetNbinsX(); tt++)
  {
    fKinematicsRatio[4][2]->SetBinError(tt,fError[tt]);
  }
  fError.clear();
  fKinematicsRatio[4][2]->SetMarkerStyle(21);
  fKinematicsRatio[4][2]->SetFillColor(kYellow-7);
  fKinematicsRatio[4][2]->SetMaximum(2.);
  fKinematicsRatio[4][2]->SetMinimum(0.);
  fKinematicsRatio[4][2]->Draw("PE2");
  fKinematicsRatio[4][2]->GetXaxis()->SetLabelSize(0.08);
  fKinematicsRatio[4][2]->GetYaxis()->SetLabelSize(0.08);
  fKinematicsRatio[4][2]->GetYaxis()->SetNdivisions(2,kTRUE);
  for(int tt=0; tt<7; tt++)
  {
    l2[2][tt]->Draw();
  }
  c10.Update();
  c10.cd(1);
  fKinematicsMC2[4][2]->SetLineColor(kGreen);
  fKinematicsMC2[4][2]->SetStats(0);
  fKinematicsMC2[4][2]->Draw();
  fKinematicsMC2[4][2]->GetXaxis()->SetLabelSize(0.08);
  fKinematicsMC2[4][2]->GetYaxis()->SetLabelSize(0.08);
  fKinematicsMC2[4][2]->GetYaxis()->SetNdivisions(304,kTRUE);
  fKinematicsMC1[4][2]->SetLineColor(kBlue);
  fKinematicsMC1[4][2]->Draw("SAME");
  c10.Update();

  c11.cd(2);
  for(int tt=0; tt<fKinematicsMC1[4][3]->GetNbinsX(); tt++)
  {
    fError.push_back((fKinematicsMC1[4][3]->GetBinError(tt) && fKinematicsMC2[4][3]->GetBinError(tt) ? sqrt(pow(1/fKinematicsMC1[4][3]->GetBinError(tt),2)+pow(1/fKinematicsMC2[4][3]->GetBinError(tt),2)):0));
  }
  fKinematicsMC2[4][3]->Scale(1/fKinematicsMC2[4][3]->GetEntries());
  fKinematicsMC1[4][3]->Scale(1/fKinematicsMC1[4][3]->GetEntries());
  fKinematicsRatio[4][3] = (TH1F*)fKinematicsMC1[4][3]->Clone();
  fKinematicsRatio[4][3]->SetStats(0);
  fKinematicsRatio[4][3]->Divide(fKinematicsMC2[4][3]);
  for(int tt=0; tt<fKinematicsRatio[4][3]->GetNbinsX(); tt++)
  {
    fKinematicsRatio[4][3]->SetBinError(tt,fError[tt]);
  }
  fError.clear();
  fKinematicsRatio[4][3]->SetMarkerStyle(21);
  fKinematicsRatio[4][3]->SetFillColor(kYellow-7);
  fKinematicsRatio[4][3]->SetMaximum(2.);
  fKinematicsRatio[4][3]->SetMinimum(0.);
  fKinematicsRatio[4][3]->Draw("PE2");
  fKinematicsRatio[4][3]->GetXaxis()->SetLabelSize(0.08);
  fKinematicsRatio[4][3]->GetYaxis()->SetLabelSize(0.08);
  fKinematicsRatio[4][3]->GetYaxis()->SetNdivisions(2,kTRUE);
  for(int tt=0; tt<7; tt++)
  {
    l2[3][tt]->Draw();
  }
  c11.Update();
  c11.cd(1);
  fKinematicsMC2[4][3]->SetLineColor(kGreen);
  fKinematicsMC2[4][3]->SetStats(0);
  fKinematicsMC2[4][3]->Draw();
  fKinematicsMC2[4][3]->GetXaxis()->SetLabelSize(0.08);
  fKinematicsMC2[4][3]->GetYaxis()->SetLabelSize(0.08);
  fKinematicsMC2[4][3]->GetYaxis()->SetNdivisions(304,kTRUE);
  fKinematicsMC1[4][3]->SetLineColor(kBlue);
  fKinematicsMC1[4][3]->Draw("SAME");
  c11.Update();

  c12.cd(2);
  for(int tt=0; tt<fKinematicsMC1[4][4]->GetNbinsX(); tt++)
  {
    fError.push_back((fKinematicsMC1[4][4]->GetBinError(tt) && fKinematicsMC2[4][4]->GetBinError(tt) ? sqrt(pow(1/fKinematicsMC1[4][4]->GetBinError(tt),2)+pow(1/fKinematicsMC2[4][4]->GetBinError(tt),2)):0));
  }
  fKinematicsMC2[4][4]->Scale(1/fKinematicsMC2[4][4]->GetEntries());
  fKinematicsMC1[4][4]->Scale(1/fKinematicsMC1[4][4]->GetEntries());
  fKinematicsRatio[4][4] = (TH1F*)fKinematicsMC1[4][4]->Clone();
  fKinematicsRatio[4][4]->SetStats(0);
  fKinematicsRatio[4][4]->Divide(fKinematicsMC2[4][4]);
  for(int tt=0; tt<fKinematicsRatio[4][4]->GetNbinsX(); tt++)
  {
    fKinematicsRatio[4][4]->SetBinError(tt,fError[tt]);
  }
  fError.clear();
  fKinematicsRatio[4][4]->SetMarkerStyle(21);
  fKinematicsRatio[4][4]->SetFillColor(kYellow-7);
  fKinematicsRatio[4][4]->SetMaximum(2.);
  fKinematicsRatio[4][4]->SetMinimum(0.);
  fKinematicsRatio[4][4]->Draw("PE2");
  fKinematicsRatio[4][4]->GetXaxis()->SetLabelSize(0.08);
  fKinematicsRatio[4][4]->GetYaxis()->SetLabelSize(0.08);
  fKinematicsRatio[4][4]->GetYaxis()->SetNdivisions(2,kTRUE);
  for(int tt=0; tt<7; tt++)
  {
    l2[4][tt]->Draw();
  }
  c12.Update();
  c12.cd(1);
  fKinematicsMC2[4][4]->SetLineColor(kGreen);
  fKinematicsMC2[4][4]->SetStats(0);
  fKinematicsMC2[4][4]->Draw();
  fKinematicsMC2[4][4]->GetXaxis()->SetLabelSize(0.08);
  fKinematicsMC2[4][4]->GetYaxis()->SetLabelSize(0.08);
  fKinematicsMC2[4][4]->GetYaxis()->SetNdivisions(304,kTRUE);
  fKinematicsMC1[4][4]->SetLineColor(kBlue);
  fKinematicsMC1[4][4]->Draw("SAME");
  c12.Update();

  c13.cd(2);
  for(int tt=0; tt<fKinematicsMC1[4][5]->GetNbinsX(); tt++)
  {
    fError.push_back((fKinematicsMC1[4][5]->GetBinError(tt) && fKinematicsMC2[4][5]->GetBinError(tt) ? sqrt(pow(1/fKinematicsMC1[4][5]->GetBinError(tt),2)+pow(1/fKinematicsMC2[4][5]->GetBinError(tt),2)):0));
  }
  fKinematicsMC2[4][5]->Scale(1/fKinematicsMC2[4][5]->GetEntries());
  fKinematicsMC1[4][5]->Scale(1/fKinematicsMC1[4][5]->GetEntries());
  fKinematicsRatio[4][5] = (TH1F*)fKinematicsMC1[4][5]->Clone();
  fKinematicsRatio[4][5]->SetStats(0);
  fKinematicsRatio[4][5]->Divide(fKinematicsMC2[4][5]);
  for(int tt=0; tt<fKinematicsRatio[4][5]->GetNbinsX(); tt++)
  {
    fKinematicsRatio[4][5]->SetBinError(tt,fError[tt]);
  }
  fError.clear();
  fKinematicsRatio[4][5]->SetMarkerStyle(21);
  fKinematicsRatio[4][5]->SetFillColor(kYellow-7);
  fKinematicsRatio[4][5]->SetMaximum(2.);
  fKinematicsRatio[4][5]->SetMinimum(0.);
  fKinematicsRatio[4][5]->Draw("PE2");
  fKinematicsRatio[4][5]->GetXaxis()->SetLabelSize(0.08);
  fKinematicsRatio[4][5]->GetYaxis()->SetLabelSize(0.08);
  fKinematicsRatio[4][5]->GetYaxis()->SetNdivisions(2,kTRUE);
  for(int tt=0; tt<7; tt++)
  {
    l2[5][tt]->Draw();
  }
  c13.Update();
  c13.cd(1);
  fKinematicsMC2[4][5]->SetLineColor(kGreen);
  fKinematicsMC2[4][5]->SetStats(0);
  fKinematicsMC2[4][5]->Draw();
  fKinematicsMC2[4][5]->GetXaxis()->SetLabelSize(0.08);
  fKinematicsMC2[4][5]->GetYaxis()->SetLabelSize(0.08);
  fKinematicsMC2[4][5]->GetYaxis()->SetNdivisions(304,kTRUE);
  fKinematicsMC1[4][5]->SetLineColor(kBlue);
  fKinematicsMC1[4][5]->Draw("SAME");
  c13.Update();

  c25.cd(2);
  for(int tt=0; tt<fKinematicsMC1[4][12]->GetNbinsX(); tt++)
  {
    fError.push_back((fKinematicsMC1[4][12]->GetBinError(tt) && fKinematicsMC2[4][12]->GetBinError(tt) ? sqrt(pow(1/fKinematicsMC1[4][12]->GetBinError(tt),2)+pow(1/fKinematicsMC2[4][12]->GetBinError(tt),2)):0));
  }
  fKinematicsMC2[4][12]->Scale(1/fKinematicsMC2[4][12]->GetEntries());
  fKinematicsMC1[4][12]->Scale(1/fKinematicsMC1[4][12]->GetEntries());
  fKinematicsRatio[4][12] = (TH1F*)fKinematicsMC1[4][12]->Clone();
  fKinematicsRatio[4][12]->SetStats(0);
  fKinematicsRatio[4][12]->Divide(fKinematicsMC2[4][12]);
  for(int tt=0; tt<fKinematicsRatio[4][12]->GetNbinsX(); tt++)
  {
    fKinematicsRatio[4][12]->SetBinError(tt,fError[tt]);
  }
  fError.clear();
  fKinematicsRatio[4][12]->SetMarkerStyle(21);
  fKinematicsRatio[4][12]->SetFillColor(kYellow-7);
  fKinematicsRatio[4][12]->SetMaximum(2.);
  fKinematicsRatio[4][12]->SetMinimum(0.);
  fKinematicsRatio[4][12]->Draw("PE2");
  fKinematicsRatio[4][12]->GetXaxis()->SetLabelSize(0.08);
  fKinematicsRatio[4][12]->GetYaxis()->SetLabelSize(0.08);
  fKinematicsRatio[4][12]->GetYaxis()->SetNdivisions(2,kTRUE);
  for(int tt=0; tt<7; tt++)
  {
    l2[12][tt]->Draw();
  }
  c25.Update();
  c25.cd(1);
  fKinematicsMC2[4][12]->SetLineColor(kGreen);
  fKinematicsMC2[4][12]->SetStats(0);
  fKinematicsMC2[4][12]->Draw();
  fKinematicsMC2[4][12]->GetXaxis()->SetLabelSize(0.08);
  fKinematicsMC2[4][12]->GetYaxis()->SetLabelSize(0.08);
  fKinematicsMC2[4][12]->GetYaxis()->SetNdivisions(304,kTRUE);
  fKinematicsMC1[4][12]->SetLineColor(kBlue);
  fKinematicsMC1[4][12]->Draw("SAME");
  c25.Update();

  c26.cd(2);
  for(int tt=0; tt<fKinematicsMC1[4][13]->GetNbinsX(); tt++)
  {
    fError.push_back((fKinematicsMC1[4][13]->GetBinError(tt) && fKinematicsMC2[4][13]->GetBinError(tt) ? sqrt(pow(1/fKinematicsMC1[4][13]->GetBinError(tt),2)+pow(1/fKinematicsMC2[4][13]->GetBinError(tt),2)):0));
  }
  fKinematicsMC2[4][13]->Scale(1/fKinematicsMC2[4][13]->GetEntries());
  fKinematicsMC1[4][13]->Scale(1/fKinematicsMC1[4][13]->GetEntries());
  fKinematicsRatio[4][13] = (TH1F*)fKinematicsMC1[4][13]->Clone();
  fKinematicsRatio[4][13]->SetStats(0);
  fKinematicsRatio[4][13]->Divide(fKinematicsMC2[4][13]);
  for(int tt=0; tt<fKinematicsRatio[4][13]->GetNbinsX(); tt++)
  {
    fKinematicsRatio[4][13]->SetBinError(tt,fError[tt]);
  }
  fError.clear();
  fKinematicsRatio[4][13]->SetMarkerStyle(21);
  fKinematicsRatio[4][13]->SetFillColor(kYellow-7);
  fKinematicsRatio[4][13]->SetMaximum(2.);
  fKinematicsRatio[4][13]->SetMinimum(0.);
  fKinematicsRatio[4][13]->Draw("PE2");
  fKinematicsRatio[4][13]->GetXaxis()->SetLabelSize(0.08);
  fKinematicsRatio[4][13]->GetYaxis()->SetLabelSize(0.08);
  fKinematicsRatio[4][13]->GetYaxis()->SetNdivisions(2,kTRUE);
  for(int tt=0; tt<7; tt++)
  {
    l2[13][tt]->Draw();
  }
  c26.Update();
  c26.cd(1);
  fKinematicsMC2[4][13]->SetLineColor(kGreen);
  fKinematicsMC2[4][13]->SetStats(0);
  fKinematicsMC2[4][13]->Draw();
  fKinematicsMC2[4][13]->GetXaxis()->SetLabelSize(0.08);
  fKinematicsMC2[4][13]->GetYaxis()->SetLabelSize(0.08);
  fKinematicsMC2[4][13]->GetYaxis()->SetNdivisions(304,kTRUE);
  fKinematicsMC1[4][13]->SetLineColor(kBlue);
  fKinematicsMC1[4][13]->Draw("SAME");
  c26.Update();

  c27.cd(2);
  for(int tt=0; tt<fKinematicsMC1[4][14]->GetNbinsX(); tt++)
  {
    fError.push_back((fKinematicsMC1[4][14]->GetBinError(tt) && fKinematicsMC2[4][14]->GetBinError(tt) ? sqrt(pow(1/fKinematicsMC1[4][14]->GetBinError(tt),2)+pow(1/fKinematicsMC2[4][14]->GetBinError(tt),2)):0));
  }
  fKinematicsMC2[4][14]->Scale(1/fKinematicsMC2[4][14]->GetEntries());
  fKinematicsMC1[4][14]->Scale(1/fKinematicsMC1[4][14]->GetEntries());
  fKinematicsRatio[4][14] = (TH1F*)fKinematicsMC1[4][14]->Clone();
  fKinematicsRatio[4][14]->SetStats(0);
  fKinematicsRatio[4][14]->Divide(fKinematicsMC2[4][14]);
  for(int tt=0; tt<fKinematicsRatio[4][14]->GetNbinsX(); tt++)
  {
    fKinematicsRatio[4][14]->SetBinError(tt,fError[tt]);
  }
  fError.clear();
  fKinematicsRatio[4][14]->SetMarkerStyle(21);
  fKinematicsRatio[4][14]->SetFillColor(kYellow-7);
  fKinematicsRatio[4][14]->SetMaximum(2.);
  fKinematicsRatio[4][14]->SetMinimum(0.);
  fKinematicsRatio[4][14]->Draw("PE2");
  fKinematicsRatio[4][14]->GetXaxis()->SetLabelSize(0.08);
  fKinematicsRatio[4][14]->GetYaxis()->SetLabelSize(0.08);
  fKinematicsRatio[4][14]->GetYaxis()->SetNdivisions(2,kTRUE);
  for(int tt=0; tt<7; tt++)
  {
    l2[14][tt]->Draw();
  }
  c27.Update();
  c27.cd(1);
  fKinematicsMC2[4][14]->SetLineColor(kGreen);
  fKinematicsMC2[4][14]->SetStats(0);
  fKinematicsMC2[4][14]->Draw();
  fKinematicsMC2[4][14]->GetXaxis()->SetLabelSize(0.08);
  fKinematicsMC2[4][14]->GetYaxis()->SetLabelSize(0.08);
  fKinematicsMC2[4][14]->GetYaxis()->SetNdivisions(304,kTRUE);
  fKinematicsMC1[4][14]->SetLineColor(kBlue);
  fKinematicsMC1[4][14]->Draw("SAME");
  c27.Update();

  c29.cd(2);
  for(int tt=0; tt<fKinematicsMC1[4][15]->GetNbinsX(); tt++)
  {
    fError.push_back((fKinematicsMC1[4][15]->GetBinError(tt) && fKinematicsMC2[4][15]->GetBinError(tt) ? sqrt(pow(1/fKinematicsMC1[4][15]->GetBinError(tt),2)+pow(1/fKinematicsMC2[4][15]->GetBinError(tt),2)):0));
  }
  fKinematicsMC2[4][15]->Scale(1/fKinematicsMC2[4][15]->GetEntries());
  fKinematicsMC1[4][15]->Scale(1/fKinematicsMC1[4][15]->GetEntries());
  fKinematicsRatio[4][15] = (TH1F*)fKinematicsMC1[4][15]->Clone();
  fKinematicsRatio[4][15]->SetStats(0);
  fKinematicsRatio[4][15]->Divide(fKinematicsMC2[4][15]);
  for(int tt=0; tt<fKinematicsRatio[4][15]->GetNbinsX(); tt++)
  {
    fKinematicsRatio[4][15]->SetBinError(tt,fError[tt]);
  }
  fError.clear();
  fKinematicsRatio[4][15]->SetMarkerStyle(21);
  fKinematicsRatio[4][15]->SetFillColor(kYellow-7);
  fKinematicsRatio[4][15]->SetMaximum(2.);
  fKinematicsRatio[4][15]->SetMinimum(0.);
  fKinematicsRatio[4][15]->Draw("PE2");
  fKinematicsRatio[4][15]->GetXaxis()->SetLabelSize(0.08);
  fKinematicsRatio[4][15]->GetYaxis()->SetLabelSize(0.08);
  fKinematicsRatio[4][15]->GetYaxis()->SetNdivisions(2,kTRUE);
  for(int tt=0; tt<7; tt++)
  {
    l2[15][tt]->Draw();
  }
  c29.Update();
  c29.cd(1);
  fKinematicsMC2[4][15]->SetLineColor(kGreen);
  fKinematicsMC2[4][15]->SetStats(0);
  fKinematicsMC2[4][15]->Draw();
  fKinematicsMC2[4][15]->GetXaxis()->SetLabelSize(0.08);
  fKinematicsMC2[4][15]->GetYaxis()->SetLabelSize(0.08);
  fKinematicsMC2[4][15]->GetYaxis()->SetNdivisions(304,kTRUE);
  fKinematicsMC1[4][15]->SetLineColor(kBlue);
  fKinematicsMC1[4][15]->Draw("SAME");
  c29.Update();

  c30.cd(2);
  for(int tt=0; tt<fKinematicsMC1[4][16]->GetNbinsX(); tt++)
  {
    fError.push_back((fKinematicsMC1[4][16]->GetBinError(tt) && fKinematicsMC2[4][16]->GetBinError(tt) ? sqrt(pow(1/fKinematicsMC1[4][16]->GetBinError(tt),2)+pow(1/fKinematicsMC2[4][16]->GetBinError(tt),2)):0));
  }
  fKinematicsMC2[4][16]->Scale(1/fKinematicsMC2[4][16]->GetEntries());
  fKinematicsMC1[4][16]->Scale(1/fKinematicsMC1[4][16]->GetEntries());
  fKinematicsRatio[4][16] = (TH1F*)fKinematicsMC1[4][16]->Clone();
  fKinematicsRatio[4][16]->SetStats(0);
  fKinematicsRatio[4][16]->Divide(fKinematicsMC2[4][16]);
  for(int tt=0; tt<fKinematicsRatio[4][16]->GetNbinsX(); tt++)
  {
    fKinematicsRatio[4][16]->SetBinError(tt,fError[tt]);
  }
  fError.clear();
  fKinematicsRatio[4][16]->SetMarkerStyle(21);
  fKinematicsRatio[4][16]->SetFillColor(kYellow-7);
  fKinematicsRatio[4][16]->SetMaximum(2.);
  fKinematicsRatio[4][16]->SetMinimum(0.);
  fKinematicsRatio[4][16]->Draw("PE2");
  fKinematicsRatio[4][16]->GetXaxis()->SetLabelSize(0.08);
  fKinematicsRatio[4][16]->GetYaxis()->SetLabelSize(0.08);
  fKinematicsRatio[4][16]->GetYaxis()->SetNdivisions(2,kTRUE);
  for(int tt=0; tt<7; tt++)
  {
    l2[16][tt]->Draw();
  }
  c30.Update();
  c30.cd(1);
  fKinematicsMC2[4][16]->SetLineColor(kGreen);
  fKinematicsMC2[4][16]->SetStats(0);
  fKinematicsMC2[4][16]->Draw();
  fKinematicsMC2[4][16]->GetXaxis()->SetLabelSize(0.08);
  fKinematicsMC2[4][16]->GetYaxis()->SetLabelSize(0.08);
  fKinematicsMC2[4][16]->GetYaxis()->SetNdivisions(304,kTRUE);
  fKinematicsMC1[4][16]->SetLineColor(kBlue);
  fKinematicsMC1[4][16]->Draw("SAME");
  c30.Update();


  c1.Print("kinMCMC.pdf(","pdf");
  c2.Print("kinMCMC.pdf","pdf");
  c3.Print("kinMCMC.pdf","pdf");
  c4.Print("kinMCMC.pdf","pdf");
  c5.Print("kinMCMC.pdf","pdf");
  c6.Print("kinMCMC.pdf","pdf");
  c22.Print("kinMCMC.pdf","pdf");
  c23.Print("kinMCMC.pdf","pdf");
  c24.Print("kinMCMC.pdf","pdf");
  c28.Print("kinMCMC.pdf","pdf");
  c7.Print("kinMCMC.pdf","pdf");
  c8.Print("kinMCMC.pdf","pdf");
  c9.Print("kinMCMC.pdf","pdf");
  c10.Print("kinMCMC.pdf","pdf");
  c11.Print("kinMCMC.pdf","pdf");
  c12.Print("kinMCMC.pdf","pdf");
  c13.Print("kinMCMC.pdf","pdf");
  c25.Print("kinMCMC.pdf","pdf");
  c26.Print("kinMCMC.pdf","pdf");
  c27.Print("kinMCMC.pdf","pdf");
  c29.Print("kinMCMC.pdf","pdf");
  c30.Print("kinMCMC.pdf","pdf");
  c14.Print("kinMCMC.pdf","pdf");
  c15.Print("kinMCMC.pdf","pdf");
  c16.Print("kinMCMC.pdf","pdf");
  c17.Print("kinMCMC.pdf","pdf");
  c18.Print("kinMCMC.pdf)","pdf");
}

void MC1extraction(string pFilelist)
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

      Double_t mxc, myc;

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

        if(Y2012)
        {
          CellCenter(z->GetLeaf("z")->GetValue(), mxc, myc);
        }

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
              if(InTarget(x->GetLeaf("x")->GetValue(),y->GetLeaf("y")->GetValue(),z->GetLeaf("z")->GetValue(),fRcutval[zlabbin]))
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
              if(InTarget(x->GetLeaf("x")->GetValue(),y->GetLeaf("y")->GetValue(),z->GetLeaf("z")->GetValue(),1.5))
              {
                fTarg++;

                // Cells crossing
                if(true)
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

      if(fAllDISflag)
      {
        double theta_m = asin(sqrt(pow(p1x->GetLeaf("p1x")->GetValue()/sqrt(pow(E_mu_prim->GetLeaf("E_mu_prim")->GetValue(),2)-pow(fM_mu,2)),2)+pow(p1y->GetLeaf("p1y")->GetValue()/sqrt(pow(E_mu_prim->GetLeaf("E_mu_prim")->GetValue(),2)-pow(fM_mu,2)),2)));
        double phi_m = asin(p1x->GetLeaf("p1x")->GetValue()/sqrt(pow(p1x->GetLeaf("p1x")->GetValue(),2)+pow(p1y->GetLeaf("p1y")->GetValue(),2)));

        // MT
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
        // LT
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
        // OT
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
        // LAST
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
        // ALL TRIGGERS
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

      // -----------------------------------------------------------------------
      // MC --------------------------------------------------------------------
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

      // -----------------------------------------------------------------------
      //  Data -----------------------------------------------------------------
      // -----------------------------------------------------------------------

      if(fAllDISflag)
      {
        for(int i=0; i<12; i++)
          fNDIS_evt_MC1[xbin][ybin][i]++;

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
            fKinematicsMC1[0][11]->Fill(abs(ph_pl->GetLeaf("Hadrons.ph_pl")->GetValue(i)));

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
          if(trig&2)
          {
            fKinematicsMC1[0][3]->Fill(zBj);
            fKinematicsMC1[0][12]->Fill(p->GetLeaf("Hadrons.P")->GetValue(i));
            fKinematicsMC1[0][13]->Fill(th->GetLeaf("Hadrons.th")->GetValue(i));
            fKinematicsMC1[0][14]->Fill(ph->GetLeaf("Hadrons.ph")->GetValue(i));
            fKinematicsMC1[0][15]->Fill(ph_pl->GetLeaf("Hadrons.ph_pl")->GetValue(i));
            fKinematicsMC1[0][16]->Fill(pt->GetLeaf("Hadrons.pt")->GetValue(i));
          }
          // LT
          if(trig&4)
          {
            fKinematicsMC1[1][3]->Fill(zBj);
            fKinematicsMC1[1][12]->Fill(p->GetLeaf("Hadrons.P")->GetValue(i));
            fKinematicsMC1[1][13]->Fill(th->GetLeaf("Hadrons.th")->GetValue(i));
            fKinematicsMC1[1][14]->Fill(ph->GetLeaf("Hadrons.ph")->GetValue(i));
            fKinematicsMC1[1][15]->Fill(ph_pl->GetLeaf("Hadrons.ph_pl")->GetValue(i));
            fKinematicsMC1[1][16]->Fill(pt->GetLeaf("Hadrons.pt")->GetValue(i));
          }
          // OT
          if(trig&8)
          {
            fKinematicsMC1[2][3]->Fill(zBj);
            fKinematicsMC1[2][12]->Fill(p->GetLeaf("Hadrons.P")->GetValue(i));
            fKinematicsMC1[2][13]->Fill(th->GetLeaf("Hadrons.th")->GetValue(i));
            fKinematicsMC1[2][14]->Fill(ph->GetLeaf("Hadrons.ph")->GetValue(i));
            fKinematicsMC1[2][15]->Fill(ph_pl->GetLeaf("Hadrons.ph_pl")->GetValue(i));
            fKinematicsMC1[2][16]->Fill(pt->GetLeaf("Hadrons.pt")->GetValue(i));
          }
          // LAST
          if(trig&512)
          {
            fKinematicsMC1[3][3]->Fill(zBj);
            fKinematicsMC1[3][12]->Fill(p->GetLeaf("Hadrons.P")->GetValue(i));
            fKinematicsMC1[3][13]->Fill(th->GetLeaf("Hadrons.th")->GetValue(i));
            fKinematicsMC1[3][14]->Fill(ph->GetLeaf("Hadrons.ph")->GetValue(i));
            fKinematicsMC1[3][15]->Fill(ph_pl->GetLeaf("Hadrons.ph_pl")->GetValue(i));
            fKinematicsMC1[3][16]->Fill(pt->GetLeaf("Hadrons.pt")->GetValue(i));
          }
          // ALL TRIGGERS
          if(trig&2 || trig&4 || trig&8)
          // if(trig&2 || trig&4 || trig&8 || trig&512)
          {
            fKinematicsMC1[4][3]->Fill(zBj);
            fKinematicsMC1[4][12]->Fill(p->GetLeaf("Hadrons.P")->GetValue(i));
            fKinematicsMC1[4][13]->Fill(th->GetLeaf("Hadrons.th")->GetValue(i));
            fKinematicsMC1[4][14]->Fill(ph->GetLeaf("Hadrons.ph")->GetValue(i));
            fKinematicsMC1[4][15]->Fill(ph_pl->GetLeaf("Hadrons.ph_pl")->GetValue(i));
            fKinematicsMC1[4][16]->Fill(pt->GetLeaf("Hadrons.pt")->GetValue(i));
          }

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


          // **********************************************************************

          // Save of hadrons

          if(fId==0)
          {
            fMC1[xbin][ybin][zbin].tab[1][0][0] += 1;
            fMC1[xbin][ybin][zbin].tab[1][0][3] += 1;
          }
          else if(fId==1)
          {
            fMC1[xbin][ybin][zbin].tab[0][0][0] += 1;
            fMC1[xbin][ybin][zbin].tab[0][0][3] += 1;
          }
          else if(fId==2)
          {
            fMC1[xbin][ybin][zbin].tab[1][0][1] += 1;
            fMC1[xbin][ybin][zbin].tab[1][0][3] += 1;
          }
          else if(fId==3)
          {
            fMC1[xbin][ybin][zbin].tab[0][0][1] += 1;
            fMC1[xbin][ybin][zbin].tab[0][0][3] += 1;
          }
          else if(fId==4)
          {
            fMC1[xbin][ybin][zbin].tab[1][0][2] += 1;
            fMC1[xbin][ybin][zbin].tab[1][0][3] += 1;
          }
          else if(fId==5)
          {
            fMC1[xbin][ybin][zbin].tab[0][0][2] += 1;
            fMC1[xbin][ybin][zbin].tab[0][0][3] += 1;
          }
          else if(fId==6)
          {
            fMC1[xbin][ybin][zbin].tab[1][0][3] += 1;
          }
          else if(fId==7)
          {
            fMC1[xbin][ybin][zbin].tab[0][0][3] += 1;
          }
          else
          {
            continue;
          }
        }
      }
    }

    cout << "\n-> Finished processing file " << filename << " <-\n" << endl;

    delete f;
  }

  for(int i=0; i<int(fQ2kinMC[0].size()); i++)
  {
    fKinematicsMC1[0][0]->Fill(fQ2kinMC[0][i]);
    fKinematicsMC1[0][1]->Fill(fXBjkinMC[0][i]);
    fKinematicsMC1[0][2]->Fill(fYBjkinMC[0][i]);
    fKinematicsMC1[0][4]->Fill(fWBjkinMC[0][i]);
    fKinematicsMC1[0][5]->Fill(fNukinMC[0][i]);
    fKinematicsMC1[0][6]->Fill(fMuMC[0][i]);
    fKinematicsMC1[0][7]->Fill(fMupMC[0][i]);
    fKinematicsMC1[0][8]->Fill(fThetaMC[0][i]);
    fKinematicsMC1[0][9]->Fill(fPhiMC[0][i]);
    fKinematicsMC1[0][10]->Fill(fVertexMC[0][i]);
  }
  for(int i=0; i<int(fQ2kinMC[1].size()); i++)
  {
    fKinematicsMC1[1][0]->Fill(fQ2kinMC[1][i]);
    fKinematicsMC1[1][1]->Fill(fXBjkinMC[1][i]);
    fKinematicsMC1[1][2]->Fill(fYBjkinMC[1][i]);
    fKinematicsMC1[1][4]->Fill(fWBjkinMC[1][i]);
    fKinematicsMC1[1][5]->Fill(fNukinMC[1][i]);
    fKinematicsMC1[1][6]->Fill(fMuMC[1][i]);
    fKinematicsMC1[1][7]->Fill(fMupMC[1][i]);
    fKinematicsMC1[1][8]->Fill(fThetaMC[1][i]);
    fKinematicsMC1[1][9]->Fill(fPhiMC[1][i]);
    fKinematicsMC1[1][10]->Fill(fVertexMC[1][i]);
  }
  for(int i=0; i<int(fQ2kinMC[2].size()); i++)
  {
    fKinematicsMC1[2][0]->Fill(fQ2kinMC[2][i]);
    fKinematicsMC1[2][1]->Fill(fXBjkinMC[2][i]);
    fKinematicsMC1[2][2]->Fill(fYBjkinMC[2][i]);
    fKinematicsMC1[2][4]->Fill(fWBjkinMC[2][i]);
    fKinematicsMC1[2][5]->Fill(fNukinMC[2][i]);
    fKinematicsMC1[2][6]->Fill(fMuMC[2][i]);
    fKinematicsMC1[2][7]->Fill(fMupMC[2][i]);
    fKinematicsMC1[2][8]->Fill(fThetaMC[2][i]);
    fKinematicsMC1[2][9]->Fill(fPhiMC[2][i]);
    fKinematicsMC1[2][10]->Fill(fVertexMC[2][i]);
  }
  for(int i=0; i<int(fQ2kinMC[3].size()); i++)
  {
    fKinematicsMC1[3][0]->Fill(fQ2kinMC[3][i]);
    fKinematicsMC1[3][1]->Fill(fXBjkinMC[3][i]);
    fKinematicsMC1[3][2]->Fill(fYBjkinMC[3][i]);
    fKinematicsMC1[3][4]->Fill(fWBjkinMC[3][i]);
    fKinematicsMC1[3][5]->Fill(fNukinMC[3][i]);
    fKinematicsMC1[3][6]->Fill(fMuMC[3][i]);
    fKinematicsMC1[3][7]->Fill(fMupMC[3][i]);
    fKinematicsMC1[3][8]->Fill(fThetaMC[3][i]);
    fKinematicsMC1[3][9]->Fill(fPhiMC[3][i]);
    fKinematicsMC1[3][10]->Fill(fVertexMC[3][i]);
  }
  for(int i=0; i<int(fQ2kinMC[4].size()); i++)
  {
    fKinematicsMC1[4][0]->Fill(fQ2kinMC[4][i]);
    fKinematicsMC1[4][1]->Fill(fXBjkinMC[4][i]);
    fKinematicsMC1[4][2]->Fill(fYBjkinMC[4][i]);
    fKinematicsMC1[4][4]->Fill(fWBjkinMC[4][i]);
    fKinematicsMC1[4][5]->Fill(fNukinMC[4][i]);
    fKinematicsMC1[4][6]->Fill(fMuMC[4][i]);
    fKinematicsMC1[4][7]->Fill(fMupMC[4][i]);
    fKinematicsMC1[4][8]->Fill(fThetaMC[4][i]);
    fKinematicsMC1[4][9]->Fill(fPhiMC[4][i]);
    fKinematicsMC1[4][10]->Fill(fVertexMC[4][i]);
  }
}

void MC2extraction(string pFilelist)
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
              if(InTarget(x->GetLeaf("x")->GetValue(),y->GetLeaf("y")->GetValue(),z->GetLeaf("z")->GetValue(),fRcutval[zlabbin]))
              {
                fTarg++;

                // Cells crossing
                if(true)
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
              if(InTarget(x->GetLeaf("x")->GetValue(),y->GetLeaf("y")->GetValue(),z->GetLeaf("z")->GetValue(),1.5))
              {
                fTarg++;

                // Cells crossing
                if(true)
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

      if(fAllDISflag)
      {
        double theta_m = asin(sqrt(pow(p1x->GetLeaf("p1x")->GetValue()/sqrt(pow(E_mu_prim->GetLeaf("E_mu_prim")->GetValue(),2)-pow(fM_mu,2)),2)+pow(p1y->GetLeaf("p1y")->GetValue()/sqrt(pow(E_mu_prim->GetLeaf("E_mu_prim")->GetValue(),2)-pow(fM_mu,2)),2)));
        double phi_m = asin(p1x->GetLeaf("p1x")->GetValue()/sqrt(pow(p1x->GetLeaf("p1x")->GetValue(),2)+pow(p1y->GetLeaf("p1y")->GetValue(),2)));

        // MT
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
        // LT
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
        // OT
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
        // LAST
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
        // ALL TRIGGERS
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

      // -------------------------------------------------------------------------
      // --------- Hadrons Selection ---------------------------------------------
      // -------------------------------------------------------------------------

      if(fAllDISflag)
      {
        for(int i=0; i<12; i++)
          fNDIS_evt_MC2[xbin][ybin][i]++;

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
            fKinematicsMC2[0][11]->Fill(abs(ph_pl->GetLeaf("Hadrons.ph_pl")->GetValue(i)));

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

          if(trig&2)
          {
            fKinematicsMC2[0][3]->Fill(zBj);
            fKinematicsMC2[0][12]->Fill(p->GetLeaf("Hadrons.P")->GetValue(i));
            fKinematicsMC2[0][13]->Fill(th->GetLeaf("Hadrons.th")->GetValue(i));
            fKinematicsMC2[0][14]->Fill(ph->GetLeaf("Hadrons.ph")->GetValue(i));
            fKinematicsMC2[0][15]->Fill(ph_pl->GetLeaf("Hadrons.ph_pl")->GetValue(i));
            fKinematicsMC2[0][16]->Fill(pt->GetLeaf("Hadrons.pt")->GetValue(i));
          }
          if(trig&4)
          {
            fKinematicsMC2[1][3]->Fill(zBj);
            fKinematicsMC2[1][12]->Fill(p->GetLeaf("Hadrons.P")->GetValue(i));
            fKinematicsMC2[1][13]->Fill(th->GetLeaf("Hadrons.th")->GetValue(i));
            fKinematicsMC2[1][14]->Fill(ph->GetLeaf("Hadrons.ph")->GetValue(i));
            fKinematicsMC2[1][15]->Fill(ph_pl->GetLeaf("Hadrons.ph_pl")->GetValue(i));
            fKinematicsMC2[1][16]->Fill(pt->GetLeaf("Hadrons.pt")->GetValue(i));
          }
          if(trig&8)
          {
            fKinematicsMC2[2][3]->Fill(zBj);
            fKinematicsMC2[2][12]->Fill(p->GetLeaf("Hadrons.P")->GetValue(i));
            fKinematicsMC2[2][13]->Fill(th->GetLeaf("Hadrons.th")->GetValue(i));
            fKinematicsMC2[2][14]->Fill(ph->GetLeaf("Hadrons.ph")->GetValue(i));
            fKinematicsMC2[2][15]->Fill(ph_pl->GetLeaf("Hadrons.ph_pl")->GetValue(i));
            fKinematicsMC2[2][16]->Fill(pt->GetLeaf("Hadrons.pt")->GetValue(i));
          }
          if(trig&512)
          {
            fKinematicsMC2[3][3]->Fill(zBj);
            fKinematicsMC2[3][12]->Fill(p->GetLeaf("Hadrons.P")->GetValue(i));
            fKinematicsMC2[3][13]->Fill(th->GetLeaf("Hadrons.th")->GetValue(i));
            fKinematicsMC2[3][14]->Fill(ph->GetLeaf("Hadrons.ph")->GetValue(i));
            fKinematicsMC2[3][15]->Fill(ph_pl->GetLeaf("Hadrons.ph_pl")->GetValue(i));
            fKinematicsMC2[3][16]->Fill(pt->GetLeaf("Hadrons.pt")->GetValue(i));
          }
          if(trig&2 || trig&4 || trig&8)
          // if(trig&2 || trig&4 || trig&8 || trig&512)
          {
            fKinematicsMC2[4][3]->Fill(zBj);
            fKinematicsMC2[4][12]->Fill(p->GetLeaf("Hadrons.P")->GetValue(i));
            fKinematicsMC2[4][13]->Fill(th->GetLeaf("Hadrons.th")->GetValue(i));
            fKinematicsMC2[4][14]->Fill(ph->GetLeaf("Hadrons.ph")->GetValue(i));
            fKinematicsMC2[4][15]->Fill(ph_pl->GetLeaf("Hadrons.ph_pl")->GetValue(i));
            fKinematicsMC2[4][16]->Fill(pt->GetLeaf("Hadrons.pt")->GetValue(i));
          }

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

          // **********************************************************************

          // Save of hadrons

          if(fId==0)
          {
            fMC2[xbin][ybin][zbin].tab[1][0][0] += 1;
            fMC2[xbin][ybin][zbin].tab[1][0][3] += 1;
          }
          else if(fId==1)
          {
            fMC2[xbin][ybin][zbin].tab[0][0][0] += 1;
            fMC2[xbin][ybin][zbin].tab[0][0][3] += 1;
          }
          else if(fId==2)
          {
            fMC2[xbin][ybin][zbin].tab[1][0][1] += 1;
            fMC2[xbin][ybin][zbin].tab[1][0][3] += 1;
          }
          else if(fId==3)
          {
            fMC2[xbin][ybin][zbin].tab[0][0][1] += 1;
            fMC2[xbin][ybin][zbin].tab[0][0][3] += 1;
          }
          else if(fId==4)
          {
            fMC2[xbin][ybin][zbin].tab[1][0][2] += 1;
            fMC2[xbin][ybin][zbin].tab[1][0][3] += 1;
          }
          else if(fId==5)
          {
            fMC2[xbin][ybin][zbin].tab[0][0][2] += 1;
            fMC2[xbin][ybin][zbin].tab[0][0][3] += 1;
          }
          else if(fId==6)
          {
            fMC2[xbin][ybin][zbin].tab[1][0][3] += 1;
          }
          else if(fId==7)
          {
            fMC2[xbin][ybin][zbin].tab[0][0][3] += 1;
          }
          else
          {
            continue;
          }
        }
      }
    }

    cout << "\n-> Finished processing file " << filename << " <-\n" << endl;

    delete f;
  }

  for(int i=0; i<int(fQ2kin[0].size()); i++)
  {
      fKinematicsMC2[0][0]->Fill(fQ2kin[0][i]);
      fKinematicsMC2[0][1]->Fill(fXBjkin[0][i]);
      fKinematicsMC2[0][2]->Fill(fYBjkin[0][i]);
      fKinematicsMC2[0][4]->Fill(fWBjkin[0][i]);
      fKinematicsMC2[0][5]->Fill(fNukin[0][i]);
      fKinematicsMC2[0][6]->Fill(fMu[0][i]);
      fKinematicsMC2[0][7]->Fill(fMup[0][i]);
      fKinematicsMC2[0][8]->Fill(fTheta[0][i]);
      fKinematicsMC2[0][9]->Fill(fPhi[0][i]);
      fKinematicsMC2[0][10]->Fill(fVertex[0][i]);
  }
  for(int i=0; i<int(fQ2kin[1].size()); i++)
  {
      fKinematicsMC2[1][0]->Fill(fQ2kin[1][i]);
      fKinematicsMC2[1][1]->Fill(fXBjkin[1][i]);
      fKinematicsMC2[1][2]->Fill(fYBjkin[1][i]);
      fKinematicsMC2[1][4]->Fill(fWBjkin[1][i]);
      fKinematicsMC2[1][5]->Fill(fNukin[1][i]);
      fKinematicsMC2[1][6]->Fill(fMu[1][i]);
      fKinematicsMC2[1][7]->Fill(fMup[1][i]);
      fKinematicsMC2[1][8]->Fill(fTheta[1][i]);
      fKinematicsMC2[1][9]->Fill(fPhi[1][i]);
      fKinematicsMC2[1][10]->Fill(fVertex[1][i]);
  }
  for(int i=0; i<int(fQ2kin[2].size()); i++)
  {
      fKinematicsMC2[2][0]->Fill(fQ2kin[2][i]);
      fKinematicsMC2[2][1]->Fill(fXBjkin[2][i]);
      fKinematicsMC2[2][2]->Fill(fYBjkin[2][i]);
      fKinematicsMC2[2][4]->Fill(fWBjkin[2][i]);
      fKinematicsMC2[2][5]->Fill(fNukin[2][i]);
      fKinematicsMC2[2][6]->Fill(fMu[2][i]);
      fKinematicsMC2[2][7]->Fill(fMup[2][i]);
      fKinematicsMC2[2][8]->Fill(fTheta[2][i]);
      fKinematicsMC2[2][9]->Fill(fPhi[2][i]);
      fKinematicsMC2[2][10]->Fill(fVertex[2][i]);
  }
  for(int i=0; i<int(fQ2kin[3].size()); i++)
  {
      fKinematicsMC2[3][0]->Fill(fQ2kin[3][i]);
      fKinematicsMC2[3][1]->Fill(fXBjkin[3][i]);
      fKinematicsMC2[3][2]->Fill(fYBjkin[3][i]);
      fKinematicsMC2[3][4]->Fill(fWBjkin[3][i]);
      fKinematicsMC2[3][5]->Fill(fNukin[3][i]);
      fKinematicsMC2[3][6]->Fill(fMu[3][i]);
      fKinematicsMC2[3][7]->Fill(fMup[3][i]);
      fKinematicsMC2[3][8]->Fill(fTheta[3][i]);
      fKinematicsMC2[3][9]->Fill(fPhi[3][i]);
      fKinematicsMC2[3][10]->Fill(fVertex[3][i]);
  }
  for(int i=0; i<int(fQ2kin[4].size()); i++)
  {
      fKinematicsMC2[4][0]->Fill(fQ2kin[4][i]);
      fKinematicsMC2[4][1]->Fill(fXBjkin[4][i]);
      fKinematicsMC2[4][2]->Fill(fYBjkin[4][i]);
      fKinematicsMC2[4][4]->Fill(fWBjkin[4][i]);
      fKinematicsMC2[4][5]->Fill(fNukin[4][i]);
      fKinematicsMC2[4][6]->Fill(fMu[4][i]);
      fKinematicsMC2[4][7]->Fill(fMup[4][i]);
      fKinematicsMC2[4][8]->Fill(fTheta[4][i]);
      fKinematicsMC2[4][9]->Fill(fPhi[4][i]);
      fKinematicsMC2[4][10]->Fill(fVertex[4][i]);
  }

}

void MultRatio()
{
  c19.Divide(5,2,0,0);
  c20.Divide(5,2,0,0);
  c21.Divide(5,2,0,0);

  TGraphErrors* H[2][9][5];
  TGraphErrors* P[2][9][5];
  TGraphErrors* K[2][9][5];

  for(int i=0; i<5; i++)
  {
    l1[i] = new TLine(0.1,1+i*0.1,0.9,1+i*0.1);
    l1[i]->SetLineStyle(4);
    l1[i]->SetLineColor(fMarkerColorAlt[i]);
  }

  double z_range[12] = {.225,.275,.325,.375,.425,.475,.525,.575,.625,.675,.725,.8};

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
          fMultRatio[i][j][k].tab[c][0][0] = ((fNDIS_evt_MC1[i][j][k] && fNDIS_evt_MC2[i][j][k] && fMC2[i][j][k].tab[c][0][0]) ?
                                                double((fMC1[i][j][k].tab[c][0][0]/fNDIS_evt_MC1[i][j][k])/(fMC2[i][j][k].tab[c][0][0]/fNDIS_evt_MC2[i][j][k])) : 0);
          fMultRatio[i][j][k].tab[c][0][1] = ((fNDIS_evt_MC1[i][j][k] && fNDIS_evt_MC2[i][j][k] && fMC2[i][j][k].tab[c][0][1]) ?
                                                double((fMC1[i][j][k].tab[c][0][1]/fNDIS_evt_MC1[i][j][k])/(fMC2[i][j][k].tab[c][0][1]/fNDIS_evt_MC2[i][j][k])) : 0);
          fMultRatio[i][j][k].tab[c][0][2] = ((fNDIS_evt_MC1[i][j][k] && fNDIS_evt_MC2[i][j][k] && fMC2[i][j][k].tab[c][0][2]) ?
                                                double((fMC1[i][j][k].tab[c][0][2]/fNDIS_evt_MC1[i][j][k])/(fMC2[i][j][k].tab[c][0][2]/fNDIS_evt_MC2[i][j][k])) : 0);
          fMultRatio[i][j][k].tab[c][0][3] = ((fNDIS_evt_MC1[i][j][k] && fNDIS_evt_MC2[i][j][k] && fMC2[i][j][k].tab[c][0][3]) ?
                                                double((fMC1[i][j][k].tab[c][0][3]/fNDIS_evt_MC1[i][j][k])/(fMC2[i][j][k].tab[c][0][3]/fNDIS_evt_MC2[i][j][k])) : 0);

          fMultRatio[i][j][k].tab[c][1][0] = ((fNDIS_evt_MC1[i][j][k] && fMC1[i][j][k].tab[c][0][0] && fNDIS_evt_MC2[i][j][k]
                                                && fMC2[i][j][k].tab[c][0][0]) ?
                                                double(sqrt(1/fMC1[i][j][k].tab[c][0][0]+1/fNDIS_evt_MC1[i][j][k]+1/fMC2[i][j][k].tab[c][0][0]+1/fNDIS_evt_MC2[i][j][k])) : 0);
          fMultRatio[i][j][k].tab[c][1][1] = ((fNDIS_evt_MC1[i][j][k] && fMC1[i][j][k].tab[c][0][1] && fNDIS_evt_MC2[i][j][k]
                                                && fMC2[i][j][k].tab[c][0][1]) ?
                                                double(sqrt(1/fMC1[i][j][k].tab[c][0][1]+1/fNDIS_evt_MC1[i][j][k]+1/fMC2[i][j][k].tab[c][0][1]+1/fNDIS_evt_MC2[i][j][k])) : 0);
          fMultRatio[i][j][k].tab[c][1][2] = ((fNDIS_evt_MC1[i][j][k] && fMC1[i][j][k].tab[c][0][2] && fNDIS_evt_MC2[i][j][k]
                                                && fMC2[i][j][k].tab[c][0][2]) ?
                                                double(sqrt(1/fMC1[i][j][k].tab[c][0][2]+1/fNDIS_evt_MC1[i][j][k]+1/fMC2[i][j][k].tab[c][0][2]+1/fNDIS_evt_MC2[i][j][k])) : 0);
          fMultRatio[i][j][k].tab[c][1][3] = ((fNDIS_evt_MC1[i][j][k] && fMC1[i][j][k].tab[c][0][3] && fNDIS_evt_MC2[i][j][k]
                                                && fMC2[i][j][k].tab[c][0][3]) ?
                                                double(sqrt(1/fMC1[i][j][k].tab[c][0][3]+1/fNDIS_evt_MC1[i][j][k]+1/fMC2[i][j][k].tab[c][0][3]+1/fNDIS_evt_MC2[i][j][k])) : 0);

          if(fMultRatio[i][j][k].tab[c][0][0]<0) fMultRatio[i][j][k].tab[c][0][0]=0;
          if(fMultRatio[i][j][k].tab[c][0][1]<0) fMultRatio[i][j][k].tab[c][0][1]=0;
          if(fMultRatio[i][j][k].tab[c][0][2]<0) fMultRatio[i][j][k].tab[c][0][2]=0;
          if(fMultRatio[i][j][k].tab[c][0][3]<0) fMultRatio[i][j][k].tab[c][0][3]=0;

          p_a.push_back(fMultRatio[i][j][k].tab[c][0][0]);
          k_a.push_back(fMultRatio[i][j][k].tab[c][0][1]);
          h_a.push_back(fMultRatio[i][j][k].tab[c][0][3]);

          p_err.push_back(fMultRatio[i][j][k].tab[c][1][0]);
          k_err.push_back(fMultRatio[i][j][k].tab[c][1][1]);
          h_err.push_back(fMultRatio[i][j][k].tab[c][1][3]);
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

        H[c][i][j] = new TGraphErrors(int(h_a.size()),&(z_range_h[0]),&(h_a[0]),0,&(h_err[0]));
        P[c][i][j] = new TGraphErrors(int(p_a.size()),&(z_range_p[0]),&(p_a[0]),0,&(p_err[0]));
        K[c][i][j] = new TGraphErrors(int(k_a.size()),&(z_range_k[0]),&(k_a[0]),0,&(k_err[0]));

        H[c][i][j]->SetMarkerColor(fMarkerColorAlt[j]);
        P[c][i][j]->SetMarkerColor(fMarkerColorAlt[j]);
        K[c][i][j]->SetMarkerColor(fMarkerColorAlt[j]);

        H[c][i][j]->SetMarkerSize(3);
        P[c][i][j]->SetMarkerSize(3);
        K[c][i][j]->SetMarkerSize(3);

        H[c][i][j]->SetMarkerStyle(fMarkerStyleAlt[j][c]);
        P[c][i][j]->SetMarkerStyle(fMarkerStyleAlt[j][c]);
        K[c][i][j]->SetMarkerStyle(fMarkerStyleAlt[j][c]);

        H[c][i][j]->GetYaxis()->SetTitle("");
        P[c][i][j]->GetYaxis()->SetTitle("");
        K[c][i][j]->GetYaxis()->SetTitle("");

        H[c][i][j]->GetXaxis()->SetTitle("");
        P[c][i][j]->GetXaxis()->SetTitle("");
        K[c][i][j]->GetXaxis()->SetTitle("");

        H[c][i][j]->SetTitle("");
        P[c][i][j]->SetTitle("");
        K[c][i][j]->SetTitle("");

        if(!h_a_empty)
        {
          c19.cd(i+1);
          gPad->SetFillStyle(4000);
          if(H[c][i][j])
          {
            if(!c && j==3)
            {
              H[c][i][j]->Draw("SAMEPA");
              H[c][i][j]->GetXaxis()->SetLimits(-0.05,1.05);
              H[c][i][j]->SetMinimum(0.85);
              H[c][i][j]->SetMaximum(1.55);
              H[c][i][j]->GetXaxis()->SetLabelSize(0.06);
              H[c][i][j]->GetYaxis()->SetLabelSize(0.06);
              H[c][i][j]->SetTitle("");
              if(i>4) gPad->SetBottomMargin(.15);
              if(i==0 || i==5) gPad->SetLeftMargin(.22);
              if(i==8)
              {
                H[c][i][j]->GetXaxis()->SetTitle("#font[ 12]{z}");
                H[c][i][j]->GetXaxis()->SetTitleSize(0.08);
                H[c][i][j]->GetXaxis()->SetTitleOffset(.8);
              }
              H[c][i][j]->GetXaxis()->SetNdivisions(304,kTRUE);
              H[c][i][j]->GetYaxis()->SetNdivisions(304,kTRUE);
              if(i==0)
              {
                H[c][i][j]->GetYaxis()->SetTitle("#font[ 12]{#eta}^{#font[ 12]{h}}+ #font[ 12]{#delta}");
                H[c][i][j]->GetYaxis()->SetTitleSize(0.08);
              }
              H[c][i][0]->Draw("SAMEP");
              H[c][i][0]->GetXaxis()->SetLimits(-0.05,1.05);
              H[c][i][0]->SetMinimum(0.85);
              H[c][i][0]->SetMaximum(1.55);
              H[c][i][1]->Draw("SAMEP");
              H[c][i][1]->GetXaxis()->SetLimits(-0.05,1.05);
              H[c][i][1]->SetMinimum(0.85);
              H[c][i][1]->SetMaximum(1.55);
              H[c][i][2]->Draw("SAMEP");
              H[c][i][2]->GetXaxis()->SetLimits(-0.05,1.05);
              H[c][i][2]->SetMinimum(0.85);
              H[c][i][2]->SetMaximum(1.55);
              c7.Range(.0,.85,1.,1.55);
              l1[0]->Draw();
              l1[1]->Draw();
              l1[2]->Draw();
              l1[3]->Draw();
              l1[4]->Draw();
            }
            else
            {
              H[c][i][j]->Draw("SAMEP");
              H[c][i][j]->GetXaxis()->SetLimits(-0.05,1.05);
              H[c][i][j]->SetMinimum(0.85);
              H[c][i][j]->SetMaximum(1.55);
            }
          }
          c19.Update();
        }

        if(!p_a_empty)
        {
          c20.cd(i+1);
          gPad->SetFillStyle(4000);
          if(P[c][i][j])
          {
            if(!c && j==3)
            {
              P[c][i][j]->Draw("SAMEPA");
              P[c][i][j]->GetXaxis()->SetLimits(-0.05,1.05);
              P[c][i][j]->SetMinimum(0.85);
              P[c][i][j]->SetMaximum(1.55);
              P[c][i][j]->GetXaxis()->SetLabelSize(0.06);
              P[c][i][j]->GetYaxis()->SetLabelSize(0.06);
              P[c][i][j]->SetTitle("");
              if(i>4) gPad->SetBottomMargin(.15);
              if(i==0 || i==5) gPad->SetLeftMargin(.22);
              if(i==8)
              {
                P[c][i][j]->GetXaxis()->SetTitle("#font[ 12]{z}");
                P[c][i][j]->GetXaxis()->SetTitleSize(0.08);
                P[c][i][j]->GetXaxis()->SetTitleOffset(.8);
              }
              P[c][i][j]->GetXaxis()->SetNdivisions(304,kTRUE);
              P[c][i][j]->GetYaxis()->SetNdivisions(304,kTRUE);
              if(i==0)
              {
                P[c][i][j]->GetYaxis()->SetTitle("#font[ 12]{#eta}^{#font[ 12]{#pi}}+ #font[ 12]{#delta}");
                P[c][i][j]->GetYaxis()->SetTitleSize(0.08);
              }
              P[c][i][0]->Draw("SAMEP");
              P[c][i][0]->GetXaxis()->SetLimits(-0.05,1.05);
              P[c][i][0]->SetMinimum(0.85);
              P[c][i][0]->SetMaximum(1.55);
              P[c][i][1]->Draw("SAMEP");
              P[c][i][1]->GetXaxis()->SetLimits(-0.05,1.05);
              P[c][i][1]->SetMinimum(0.85);
              P[c][i][1]->SetMaximum(1.55);
              P[c][i][2]->Draw("SAMEP");
              P[c][i][2]->GetXaxis()->SetLimits(-0.05,1.05);
              P[c][i][2]->SetMinimum(0.85);
              P[c][i][2]->SetMaximum(1.55);
              c8.Range(0.,.85,1.,1.55);
              l1[0]->Draw();
              l1[1]->Draw();
              l1[2]->Draw();
              l1[3]->Draw();
              l1[4]->Draw();
            }
            else
            {
              P[c][i][j]->Draw("SAMEP");
              P[c][i][j]->GetXaxis()->SetLimits(-0.05,1.05);
              P[c][i][j]->SetMinimum(0.85);
              P[c][i][j]->SetMaximum(1.55);
            }
          }
          c20.Update();
        }

        if(!k_a_empty)
        {
          c21.cd(i+1);
          gPad->SetFillStyle(4000);
          if(K[c][i][j])
          {
            if(!c && j==3)
            {
              K[c][i][j]->Draw("SAMEPA");
              K[c][i][j]->GetXaxis()->SetLimits(-0.05,1.05);
              K[c][i][j]->SetMinimum(0.85);
              K[c][i][j]->SetMaximum(1.55);
              K[c][i][j]->GetXaxis()->SetLabelSize(0.06);
              K[c][i][j]->GetYaxis()->SetLabelSize(0.06);
              K[c][i][j]->SetTitle("");
              if(i>4) gPad->SetBottomMargin(.15);
              if(i==0 || i==5) gPad->SetLeftMargin(.22);
              if(i==8)
              {
                K[c][i][j]->GetXaxis()->SetTitle("#font[ 12]{z}");
                K[c][i][j]->GetXaxis()->SetTitleSize(0.08);
                K[c][i][j]->GetXaxis()->SetTitleOffset(.8);
              }
              K[c][i][j]->GetXaxis()->SetNdivisions(304,kTRUE);
              K[c][i][j]->GetYaxis()->SetNdivisions(304,kTRUE);
              if(i==0)
              {
                K[c][i][j]->GetYaxis()->SetTitle("#font[ 12]{#eta}^{#font[ 12]{K}}+ #font[ 12]{#delta}");
                K[c][i][j]->GetYaxis()->SetTitleSize(0.08);
              }
              K[c][i][0]->Draw("SAMEP");
              K[c][i][0]->GetXaxis()->SetLimits(-0.05,1.05);
              K[c][i][0]->SetMinimum(0.85);
              K[c][i][0]->SetMaximum(1.55);
              K[c][i][1]->Draw("SAMEP");
              K[c][i][1]->GetXaxis()->SetLimits(-0.05,1.05);
              K[c][i][1]->SetMinimum(0.85);
              K[c][i][1]->SetMaximum(1.55);
              K[c][i][2]->Draw("SAMEP");
              K[c][i][2]->GetXaxis()->SetLimits(-0.05,1.05);
              K[c][i][2]->SetMinimum(0.85);
              K[c][i][2]->SetMaximum(1.55);
              c9.Range(0.,.85,1.,1.55);
              l1[0]->Draw();
              l1[1]->Draw();
              l1[2]->Draw();
              l1[3]->Draw();
              l1[4]->Draw();
            }
            else
            {
              K[c][i][j]->Draw("SAMEP");
              K[c][i][j]->GetXaxis()->SetLimits(-0.05,1.05);
              K[c][i][j]->SetMinimum(0.85);
              K[c][i][j]->SetMaximum(1.55);
            }
          }
          c21.Update();
        }
      }
    }
  }
  c19.Update();
  c20.Update();
  c21.Update();

  TLatex fTitle;

  c19.cd(1);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.5,"0.004#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.01");

  c19.cd(2);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.5,"0.01#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.02");

  c19.cd(3);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.5,"0.02#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.03");

  c19.cd(4);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.5,"0.03#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.04");

  c19.cd(5);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.5,"0.04#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.06");

  c19.cd(6);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.5,"0.06#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.1");

  c19.cd(7);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.5,"0.1#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.14");

  c19.cd(8);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.5,"0.14#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.18");

  c19.cd(9);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.5,"0.18#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.4");

  c19.cd(10);
  fTitle.SetTextSize(0.095);
  fTitle.SetTextAlign(11);
  fTitle.DrawLatex(0.05, 0.64,"#color[221]{0.50#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.70, #delta = 0.4}");
  fTitle.DrawLatex(0.05, 0.56,"#color[226]{0.30#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.50, #delta = 0.3}");
  fTitle.DrawLatex(0.05, 0.48,"#color[209]{0.20#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.30, #delta = 0.2}");
  fTitle.DrawLatex(0.05, 0.40,"#color[95]{0.15#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.20, #delta = 0.1}");
  fTitle.DrawLatex(0.05, 0.32,"#color[2]{0.10#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.15, #delta = 0}");

  c20.cd(1);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.5,"0.004#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.01");

  c20.cd(2);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.5,"0.01#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.02");

  c20.cd(3);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.5,"0.02#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.03");

  c20.cd(4);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.5,"0.03#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.04");

  c20.cd(5);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.5,"0.04#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.06");

  c20.cd(6);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.5,"0.06#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.1");

  c20.cd(7);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.5,"0.1#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.14");

  c20.cd(8);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.5,"0.14#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.18");

  c20.cd(9);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.5,"0.18#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.4");

  c20.cd(10);
  fTitle.SetTextSize(0.095);
  fTitle.SetTextAlign(11);
  fTitle.DrawLatex(0.05, 0.64,"#color[221]{0.50#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.70, #delta = 0.4}");
  fTitle.DrawLatex(0.05, 0.56,"#color[226]{0.30#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.50, #delta = 0.3}");
  fTitle.DrawLatex(0.05, 0.48,"#color[209]{0.20#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.30, #delta = 0.2}");
  fTitle.DrawLatex(0.05, 0.40,"#color[95]{0.15#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.20, #delta = 0.1}");
  fTitle.DrawLatex(0.05, 0.32,"#color[2]{0.10#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.15, #delta = 0}");

  c21.cd(1);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.5,"0.004#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.01");

  c21.cd(2);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.5,"0.01#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.02");

  c21.cd(3);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.5,"0.02#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.03");

  c21.cd(4);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.5,"0.03#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.04");

  c21.cd(5);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.5,"0.04#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.06");

  c21.cd(6);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.5,"0.06#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.1");

  c19.cd(7);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.5,"0.1#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.14");

  c21.cd(8);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.5,"0.14#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.18");

  c21.cd(9);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.5,"0.18#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.4");

  c21.cd(10);
  fTitle.SetTextSize(0.095);
  fTitle.SetTextAlign(11);
  fTitle.DrawLatex(0.05, 0.64,"#color[221]{0.50#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.70, #delta = 0.4}");
  fTitle.DrawLatex(0.05, 0.56,"#color[226]{0.30#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.50, #delta = 0.3}");
  fTitle.DrawLatex(0.05, 0.48,"#color[209]{0.20#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.30, #delta = 0.2}");
  fTitle.DrawLatex(0.05, 0.40,"#color[95]{0.15#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.20, #delta = 0.1}");
  fTitle.DrawLatex(0.05, 0.32,"#color[2]{0.10#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.15, #delta = 0}");

  c19.Update();
  c20.Update();
  c21.Update();

  c19.Print("MC_ratio_hadron.pdf");
  c20.Print("MC_ratio_pion.pdf");
  c21.Print("MC_ratio_kaon.pdf");
}

int main(int argc, char **argv)
{

  if(argc < 3)
  {
    cout << "ERROR : Not enough arguments." << endl;
    cout << "Asked : 3 *** Received : " << argc-1 << endl;
    cout << "./compMCMC [MC1 filelist] [MC2 filelist] [Cutfile]" << endl;

    return 1;
  }

  create_kin_plots();
  readKinCuts(argv[3]);
  MC1extraction(argv[1]);
  MC2extraction(argv[2]);
  MultRatio();
  save_kin_plots();

  return 0;
}

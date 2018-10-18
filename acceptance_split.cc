#include "acceptance_split.h"

//Inputs
#define target_file_2012 "data/target-107924-109081.dat"
#define target_file_2016 "data/target-274508-274901.dat"

//User dependant input
#define data_path "/sps/compass/npierre/MC"

// Flags
#define Y2006 0
#define Y2012 0
#define Y2016 1
#define MOMENTUM 12
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
          for(int n=0; n<3; n++)
          {
            fNDIS_evt[n][i][j][k] = 0;
            fNDIS_evt_c[n][i][j][k] = 0;
            fNDIS_evt_MC[n][i][j][k] = 0;
            for(int l=0; l<4; l++)
            {
              fNDIS_evt_zvtx[n][i][j][k][l] = 0;
              fNDIS_evt_MC_zvtx[n][i][j][k][l] = 0;
            }
          }
          for(int n=0; n<4; n++)
          {
            fRcstr[i][j][k].tab[c][0][n] = 0;
            fRcstr_c[i][j][k].tab[c][0][n] = 0;
            fGnrt[i][j][k].tab[c][0][n] = 0;
            for(int l=0; l<4; l++)
            {
              fRcstr_zvtx[i][j][k][l].tab[c][0][n] = 0;
              fGnrt_zvtx[i][j][k][l].tab[c][0][n] = 0;
            }
          }
        }
      }
    }
  }
}

void create_kin_plots()
{
  fKinematics[0] = new TH1F("Q^{2}", "Q^{2}", 100, 0, 2);
  fKinematics[1] = new TH1F("x_{Bj}", "x_{Bj}", 100, -3, 0);
  fKinematics[2] = new TH1F("y", "y", 100, 0, 1);
  fKinematics[3] = new TH1F("z", "z", 100, 0, 1);
  fKinematics[4] = new TH1F("W", "W", 100, 2, 18);
  fKinematics[5] = new TH1F("#nu", "#nu", 100, 0, 160);
  fKinematics[6] = new TH1F("E_{#mu}", "E_{#mu}", 100, 140, 180);
  fKinematics2D = new TH2F("DIS kin space", "DIS kin space", 100, -3, 0, 100, 0.1, 0.7);
  fTarget2D = new TH2F("Target xy", "Target xy", 100, -3, 3, 100, -3, 3);
  fHM04 = new TH2F("HM04Y1", "HM04Y1", 100, 0, 120, 100, -60, 60);
  fHM05 = new TH2F("HM05Y1", "HM05Y1", 100, 0, 120, 100, -60, 60);
  fHL04 = new TH2F("HL04X1", "HL04X1", 100, 50, 190, 100, -25, 25);
  fHL05 = new TH2F("HL05X1", "HL05X1", 100, 60, 240, 100, -25, 25);
  fHL04x1D = new TH1F("HL04X1x", "HL04X1x", 1000, 50, 190);
  fHL05x1D = new TH1F("HL05X1x", "HL05X1x", 1000, 60, 240);
  fHL04MCx1D = new TH1F("HL04MCX1x", "HL04MCX1x", 1000, 50, 190);
  fHL05MCx1D = new TH1F("HL05MCX1x", "HL05MCX1x", 1000, 60, 240);
  fHLx2D = new TH2F("HLx2D", "HLx2D", 1000, 50, 190, 1000, 60, 240);
  fHO03 = new TH2F("HO03Y1", "HO03Y1", 100, -60, 90, 100, -60, 60);
  fHO04 = new TH2F("HO04Y1", "HO04Y1", 100, -100, 250, 100, -90, 90);
  fHG01 = new TH2F("HG01Y1", "HG01Y1", 100, -100, 100, 100, -50, 50);
  fHG021 = new TH2F("HG02Y1", "HG02Y1", 100, -20, 250, 100, -230, 230);
  fHG022 = new TH2F("HG02Y2", "HG02Y2", 100, -250, 20, 100, -230, 230);
  fHM04MC = new TH2F("HM04MCY1", "HM04MCY1", 100, 0, 120, 100, -60, 60);
  fHM05MC = new TH2F("HM05MCY1", "HM05MCY1", 100, 0, 120, 100, -60, 60);
  fHL04MC = new TH2F("HL04MCX1", "HL04MCX1", 100, 50, 190, 100, -25, 25);
  fHL05MC = new TH2F("HL05MCX1", "HL05MCX1", 100, 60, 240, 100, -25, 25);
  fHL04MCx2D = new TH2F("HL04MCX2x", "HL04MCX2x", 100, -10, 10, 100, -10, 10);
  fHL05MCx2D = new TH2F("HL05MCX2x", "HL05MCX2x", 100, -10, 10, 100, -10, 10);
  fHO03MC = new TH2F("HO03MCY1", "HO03MCY1", 100, -60, 90, 100, -60, 60);
  fHO04MC = new TH2F("HO04MCY1", "HO04MCY1", 100, -100, 250, 100, -90, 90);
  fHG01MC = new TH2F("HG01MCY1", "HG01MCY1", 100, -100, 100, 100, -50, 50);
  fHG021MC = new TH2F("HG02MCY1", "HG02MCY1", 100, -20, 250, 100, -230, 230);
  fHG022MC = new TH2F("HG02MCY2", "HG02MCY2", 100, -250, 20, 100, -230, 230);
  fTrigCov = new TH2F("Trigger_Coverage", "Trigger_Coverage", 100, -227, 258, 100, -112, 112);
  fVertexHadron[0] = new TH1F("Vertex Hadron+e", "Vertex Hadron+e", 50, -325, -71);
  fVertexHadron[1] = new TH1F("Vertex e", "Vertex e", 50, -325, -71);
  fVertexHadron[2] = new TH1F("Vertex Hadron", "Vertex Hadron", 50, -325, -71);
  BinLogX(fKinematics[0]);
  BinLogX(fKinematics[1]);
  BinLogX(fKinematics2D);
  fKinematicsMC[0] = new TH1F("Q^{2} MC", "Q^{2} MC", 100, 0, 2);
  fKinematicsMC[1] = new TH1F("x_{Bj} MC", "x_{Bj} MC", 100, -3, 0);
  fKinematicsMC[2] = new TH1F("y MC", "y MC", 100, 0, 1);
  fKinematicsMC[3] = new TH1F("z MC", "z MC", 100, 0, 1);
  fKinematicsMC[4] = new TH1F("W MC", "W MC", 100, 2, 18);
  fKinematicsMC[5] = new TH1F("#nu MC", "#nu MC", 100, 0, 160);
  fKinematicsMC[6] = new TH1F("E_{#mu,MC}", "E_{#mu,MC}", 100, 140, 180);
  fKinematics2DMC = new TH2F("DIS kin space MC", "DIS kin space MC", 100, -3, 0, 100, 0.1, 0.7);
  fTarget2DMC = new TH2F("Target xy MC", "Target xy MC", 100, -3, 3, 100, -3, 3);
  fVertexHadronMC[0] = new TH1F("Vertex Hadron+e MC", "Vertex Hadron+e MC", 50, -325, -71);
  fVertexHadronMC[1] = new TH1F("Vertex e MC", "Vertex e MC", 50, -325, -71);
  fVertexHadronMC[2] = new TH1F("Vertex Hadron MC", "Vertex Hadron MC", 50, -325, -71);
  BinLogX(fKinematicsMC[0]);
  BinLogX(fKinematicsMC[1]);
  BinLogX(fKinematics2DMC);
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
  c8.Divide(1,1);
  c9.Divide(1,1);
  c10.Divide(1,1);
  c11.Divide(1,1);
  c12.Divide(1,1);
  c13.Divide(1,1);
  c14.Divide(1,1);
  c15.Divide(1,1);
  c16.Divide(1,1);
  c17.Divide(1,1);
  c26.Divide(1,1);
  c27.Divide(1,1);
  c28.Divide(1,1);
  c29.Divide(1,1);
  c30.Divide(1,1);
  c31.Divide(1,1);
  c32.Divide(1,1);
  c33.Divide(1,1);
  c34.Divide(1,1);
  c35.Divide(1,1);
  c36.Divide(1,1);
  c37.Divide(1,1);
  c38.Divide(1,1);
  c39.Divide(1,1);
  c40.Divide(1,1);
  c41.Divide(1,1);
  c42.Divide(1,1);
  c43.Divide(1,1);
  c44.Divide(1,1);
  c1.cd(1);
  fKinematics[0]->Draw();
  gPad->SetLogx();
  c1.Update();
  c2.cd(1);
  fKinematics[1]->Draw();
  gPad->SetLogx();
  c2.Update();
  c3.cd(1);
  fKinematics[2]->Draw();
  c3.Update();
  c4.cd(1);
  fKinematics[3]->Draw();
  c4.Update();
  c5.cd(1);
  fKinematics[4]->Draw();
  c5.Update();
  c6.cd(1);
  fKinematics[5]->Draw();
  c6.Update();
  c7.cd(1);
  fKinematics2D->Draw("COLZ");
  gPad->SetLogx();
  c7.Update();
  c8.cd(1);
  fTarget2D->Draw("COLZ");
  c8.Update();
  c9.cd(1);
  fHM04->Draw("COLZ");
  c9.Update();
  c10.cd(1);
  fHM05->Draw("COLZ");
  c10.Update();
  c11.cd(1);
  fHL04->Draw("COLZ");
  c11.Update();
  c12.cd(1);
  fHL05->Draw("COLZ");
  c12.Update();
  c13.cd(1);
  fHO03->Draw("COLZ");
  c13.Update();
  c14.cd(1);
  fHO04->Draw("COLZ");
  c14.Update();
  c15.cd(1);
  fHG01->Draw("COLZ");
  c15.Update();
  c16.cd(1);
  fHG021->Draw("COLZ");
  c16.Update();
  c17.cd(1);
  fHG022->Draw("COLZ");
  c17.Update();
  c26.cd(1);
  fHL04MCx1D->SetLineColor(kRed);
  fHL04MCx1D->Draw();
  fHL04x1D->Draw("SAME");
  c26.Update();
  c27.cd(1);
  fHL05MCx1D->SetLineColor(kRed);
  fHL05MCx1D->Draw();
  fHL05x1D->Draw("SAME");
  c27.Update();
  c28.cd(1);
  fHLx2D->Draw("COLZ");
  c28.Update();
  c29.cd(1);
  fHM04MC->Draw("COLZ");
  c29.Update();
  c30.cd(1);
  fHM05MC->Draw("COLZ");
  c30.Update();
  c31.cd(1);
  fHL04MC->Draw("COLZ");
  c31.Update();
  c32.cd(1);
  fHL05MC->Draw("COLZ");
  c32.Update();
  c33.cd(1);
  fHO03MC->Draw("COLZ");
  c33.Update();
  c34.cd(1);
  fHO04MC->Draw("COLZ");
  c34.Update();
  c35.cd(1);
  fHG01MC->Draw("COLZ");
  c35.Update();
  c36.cd(1);
  fHG021MC->Draw("COLZ");
  c36.Update();
  c37.cd(1);
  fHG022MC->Draw("COLZ");
  c37.Update();
  c38.cd(1);
  fHL04MCx2D->Draw("COLZ");
  c38.Update();
  c39.cd(1);
  fHL05MCx2D->Draw("COLZ");
  c39.Update();
  c40.cd(1);
  fKinematics[6]->Draw();
  fKinematicsMC[6]->SetLineColor(kRed);
  fKinematicsMC[6]->Draw("SAME");
  c40.Update();
  c41.cd(1);
  fTrigCov->SetStats(0);
  fTrigCov->Draw("COLZ");

  // MIDDLE
  TLine* M04[8];
  M04[0] = new TLine(0,-24.1719-4*(2.1-.1)+.1+7*(2.1-.1)+2.1,120,-24.1719-4*(2.1-.1)+.1+7*(2.1-.1)+2.1);
  M04[1] = new TLine(0,-41.5719-4*(2.5-.1)+.1,0,-24.1719-4*(2.1-.1)+.1+7*(2.1-.1)+2.1);
  M04[2] = new TLine(0,-41.5719-4*(2.5-.1)+.1,120,-41.5719-4*(2.5-.1)+.1);
  M04[3] = new TLine(120,-41.5719-4*(2.5-.1)+.1,120,-24.1719-4*(2.1-.1)+.1+7*(2.1-.1)+2.1);

  M04[4] = new TLine(0,23.7216-4*(2.1-.1)+.1,120,23.7216-4*(2.1-.1)+.1);
  M04[5] = new TLine(0,23.7216-4*(2.1-.1)+.1,0,41.1216-4*(2.5-.1)+.1+7*(2.5-.1)+2.5);
  M04[6] = new TLine(0,41.1216-4*(2.5-.1)+.1+7*(2.5-.1)+2.5,120,41.1216-4*(2.5-.1)+.1+7*(2.5-.1)+2.5);
  M04[7] = new TLine(120,23.7216-4*(2.1-.1)+.1,120,41.1216-4*(2.5-.1)+.1+7*(2.5-.1)+2.5);

  for(int i=0; i<8; i++) {M04[i]->SetLineColor(kCyan-4); M04[i]->Draw("SAME");}

  TLine* m1[8]; TLine* m2[8]; TLine* m3[8]; TLine* m4[8];
  TLine* m5[8]; TLine* m6[8]; TLine* m7[8]; TLine* m8[8];

  double subdetsize, detectpos;
  for(int i=0; i<8; i++)
  {
    detectpos = -24.1719;
    subdetsize = 4*(2.1-.1)+.1;
    m1[i] = new TLine(0,detectpos-subdetsize+i*(2.1-.1),120,detectpos-subdetsize+i*(2.1-.1));
    m2[i] = new TLine(0,detectpos-subdetsize+i*(2.1-.1)+2.1,120,detectpos-subdetsize+i*(2.1-.1)+2.1);
    m1[i]->SetLineColor(kCyan-4);
    m2[i]->SetLineColor(kCyan-4);
    m1[i]->Draw("SAME");
    m2[i]->Draw("SAME");
  }
  for(int i=0; i<8; i++)
  {
    detectpos = -41.5719;
    subdetsize = 4*(2.5-.1)+.1;
    m3[i] = new TLine(0,detectpos-subdetsize+i*(2.5-.1),120,detectpos-subdetsize+i*(2.5-.1));
    m4[i] = new TLine(0,detectpos-subdetsize+i*(2.5-.1)+2.5,120,detectpos-subdetsize+i*(2.5-.1)+2.5);
    m3[i]->SetLineColor(kCyan-4);
    m4[i]->SetLineColor(kCyan-4);
    m3[i]->Draw("SAME");
    m4[i]->Draw("SAME");
  }
  for(int i=0; i<8; i++)
  {
    detectpos = 23.7216;
    subdetsize = 4*(2.1-.1)+.1;
    m5[i] = new TLine(0,detectpos-subdetsize+i*(2.1-.1),120,detectpos-subdetsize+i*(2.1-.1));
    m6[i] = new TLine(0,detectpos-subdetsize+i*(2.1-.1)+2.1,120,detectpos-subdetsize+i*(2.1-.1)+2.1);
    m5[i]->SetLineColor(kCyan-4);
    m6[i]->SetLineColor(kCyan-4);
    m5[i]->Draw("SAME");
    m6[i]->Draw("SAME");
  }
  for(int i=0; i<8; i++)
  {
    detectpos = 41.1216;
    subdetsize = 4*(2.5-.1)+.1;
    m7[i] = new TLine(0,detectpos-subdetsize+i*(2.5-.1),120,detectpos-subdetsize+i*(2.5-.1));
    m8[i] = new TLine(0,detectpos-subdetsize+i*(2.5-.1)+2.5,120,detectpos-subdetsize+i*(2.5-.1)+2.5);
    m7[i]->SetLineColor(kCyan-4);
    m8[i]->SetLineColor(kCyan-4);
    m7[i]->Draw("SAME");
    m8[i]->Draw("SAME");
  }

  // LADDER
  TLine* L04[2];
  L04[0] = new TLine(62.0055-4*(2.2-.1)+.1,20,156.0055-4*(6.7-.3)+.3+7*(6.7-.3)+6.7,20);
  L04[1] = new TLine(62.0055-4*(2.2-.1)+.1,-20,156.0055-4*(6.7-.3)+.3+7*(6.7-.3)+6.7,-20);
  L04[0]->SetLineColor(kMagenta-4); L04[0]->Draw("SAME"); L04[1]->SetLineColor(kMagenta-4); L04[1]->Draw("SAME");

  TLine* l1[8]; TLine* l2[8]; TLine* l3[8]; TLine* l4[8];
  TLine* l5[8]; TLine* l6[8]; TLine* l7[8]; TLine* l8[8];

  for(int i=0; i<8; i++)
  {
    detectpos = 62.0055;
    subdetsize = 4*(2.2-.1)+.1;
    l1[i] = new TLine(detectpos-subdetsize+i*(2.2-.1),-20,detectpos-subdetsize+i*(2.2-.1),20);
    l2[i] = new TLine(detectpos-subdetsize+i*(2.2-.1)+2.2,-20,detectpos-subdetsize+i*(2.2-.1)+2.2,20);
    l1[i]->SetLineColor(kMagenta-4);
    l2[i]->SetLineColor(kMagenta-4);
    l1[i]->Draw("SAME");
    l2[i]->Draw("SAME");
  }
  for(int i=0; i<8; i++)
  {
    detectpos = 82.0055;
    subdetsize = 4*(3.2-.3)+.3;
    l3[i] = new TLine(detectpos-subdetsize+i*(3.2-.3),-20,detectpos-subdetsize+i*(3.2-.3),20);
    l4[i] = new TLine(detectpos-subdetsize+i*(3.2-.3)+3.2,-20,detectpos-subdetsize+i*(3.2-.3)+3.2,20);
    l3[i]->SetLineColor(kMagenta-4);
    l4[i]->SetLineColor(kMagenta-4);
    l3[i]->Draw("SAME");
    l4[i]->Draw("SAME");
  }
  for(int i=0; i<8; i++)
  {
    detectpos = 112.0055;
    subdetsize = 4*(4.7-.1)+.1;
    l5[i] = new TLine(detectpos-subdetsize+i*(4.7-.1),-20,detectpos-subdetsize+i*(4.7-.1),20);
    l6[i] = new TLine(detectpos-subdetsize+i*(4.7-.1)+4.7,-20,detectpos-subdetsize+i*(4.7-.1)+4.7,20);
    l5[i]->SetLineColor(kMagenta-4);
    l6[i]->SetLineColor(kMagenta-4);
    l5[i]->Draw("SAME");
    l6[i]->Draw("SAME");
  }
  for(int i=0; i<8; i++)
  {
    detectpos = 156.0055;
    subdetsize = 4*(6.7-.3)+.3;
    l7[i] = new TLine(detectpos-subdetsize+i*(6.7-.3),-20,detectpos-subdetsize+i*(6.7-.3),20);
    l8[i] = new TLine(detectpos-subdetsize+i*(6.7-.3)+6.7,-20,detectpos-subdetsize+i*(6.7-.3)+6.7,20);
    l7[i]->SetLineColor(kMagenta-4);
    l8[i]->SetLineColor(kMagenta-4);
    l7[i]->Draw("SAME");
    l8[i]->Draw("SAME");
  }

  // OUTER
  TLine* O04[6];
  O04[0] = new TLine(113,13.5,113,41.5);
  O04[1] = new TLine(113,13.5,168,13.5);
  O04[2] = new TLine(168,-13.5,168,13.5);
  O04[3] = new TLine(113,-13.5,168,-13.5);
  O04[4] = new TLine(113,-41.5,113,-13.5);
  O04[5] = new TLine(4,-41.5,4,41.5);

  for(int i=0; i<6; i++) {O04[i]->SetLineColor(kRed-7); O04[i]->Draw("SAME");}

  TLine* o1[4]; TLine* o2[8]; TLine* o3[20];
  TLine* o4[12]; TLine* o5[20];

  o1[0] = new TLine(168,-14.5,258,-14.5);
  o1[1] = new TLine(168,-.5,258,-.5);
  o1[2] = new TLine(168,.5,258,.5);
  o1[3] = new TLine(168,14.5,258,14.5);

  for(int i=0; i<4; i++) {o1[i]->SetLineColor(kRed-7); o1[i]->Draw("SAME");}

  o2[0] = new TLine(113,-13.5,258,-14.5);
  o2[1] = new TLine(113,-28.5,258,-28.5);
  o2[2] = new TLine(113,-27.5,258,-27.5);
  o2[3] = new TLine(113,-42.5,258,-42.5);
  o2[4] = new TLine(113,13.5,258,14.5);
  o2[5] = new TLine(113,28.5,258,28.5);
  o2[6] = new TLine(113,27.5,258,27.5);
  o2[7] = new TLine(113,42.5,258,42.5);

  for(int i=0; i<8; i++) {o2[i]->SetLineColor(kRed-7); o2[i]->Draw("SAME");}

  o3[0] = new TLine(8,-41.5,258,-41.5);
  o3[1] = new TLine(8,-56.5,258,-56.5);
  o3[2] = new TLine(8,-55.5,258,-55.5);
  o3[3] = new TLine(8,-70.5,258,-70.5);
  o3[4] = new TLine(8,-69.5,258,-69.5);
  o3[5] = new TLine(8,-84.5,258,-84.5);
  o3[6] = new TLine(8,-83.5,258,-83.5);
  o3[7] = new TLine(8,-98.5,258,-98.5);
  o3[8] = new TLine(8,-97.5,258,-97.5);
  o3[9] = new TLine(8,-112.5,258,-112.5);
  o3[10] = new TLine(8,41.5,258,41.5);
  o3[11] = new TLine(8,56.5,258,56.5);
  o3[12] = new TLine(8,55.5,258,55.5);
  o3[13] = new TLine(8,70.5,258,70.5);
  o3[14] = new TLine(8,69.5,258,69.5);
  o3[15] = new TLine(8,84.5,258,84.5);
  o3[16] = new TLine(8,83.5,258,83.5);
  o3[17] = new TLine(8,98.5,258,98.5);
  o3[18] = new TLine(8,97.5,258,97.5);
  o3[19] = new TLine(8,112.5,258,112.5);

  for(int i=0; i<20; i++) {o3[i]->SetLineColor(kRed-7); o3[i]->Draw("SAME");}

  o4[0] = new TLine(-227,-.5,4,-.5);
  o4[1] = new TLine(-227,.5,4,.5);
  o4[2] = new TLine(-227,-14.5,4,-14.5);
  o4[3] = new TLine(-227,-13.5,4,-13.5);
  o4[4] = new TLine(-227,14.5,4,14.5);
  o4[5] = new TLine(-227,13.5,4,13.5);
  o4[6] = new TLine(-227,-28.5,4,-28.5);
  o4[7] = new TLine(-227,-27.5,4,-27.5);
  o4[8] = new TLine(-227,28.5,4,28.5);
  o4[9] = new TLine(-227,27.5,4,27.5);
  o4[10] = new TLine(-227,-42.5,4,-42.5);
  o4[11] = new TLine(-227,42.5,4,42.5);

  for(int i=0; i<12; i++) {o4[i]->SetLineColor(kRed-7); o4[i]->Draw("SAME");}

  o5[0] = new TLine(-227,-41.5,23,-41.5);
  o5[1] = new TLine(-227,-56.5,23,-56.5);
  o5[2] = new TLine(-227,-55.5,23,-55.5);
  o5[3] = new TLine(-227,-70.5,23,-70.5);
  o5[4] = new TLine(-227,-69.5,23,-69.5);
  o5[5] = new TLine(-227,-84.5,23,-84.5);
  o5[6] = new TLine(-227,-83.5,23,-83.5);
  o5[7] = new TLine(-227,-98.5,23,-98.5);
  o5[8] = new TLine(-227,-97.5,23,-97.5);
  o5[9] = new TLine(-227,-112.5,23,-112.5);
  o5[10] = new TLine(-227,41.5,23,41.5);
  o5[11] = new TLine(-227,56.5,23,56.5);
  o5[12] = new TLine(-227,55.5,23,55.5);
  o5[13] = new TLine(-227,70.5,23,70.5);
  o5[14] = new TLine(-227,69.5,23,69.5);
  o5[15] = new TLine(-227,84.5,23,84.5);
  o5[16] = new TLine(-227,83.5,23,83.5);
  o5[17] = new TLine(-227,98.5,23,98.5);
  o5[18] = new TLine(-227,97.5,23,97.5);
  o5[19] = new TLine(-227,112.5,23,112.5);

  for(int i=0; i<20; i++) {o5[i]->SetLineColor(kRed-7); o5[i]->Draw("SAME");}

  // LAST
  // TLine* G02[4];
  // G02[0] = new TLine(-70,-40,70,-40);
  // G02[1] = new TLine(-70,-40,-70,40);
  // G02[2] = new TLine(-70,40,70,40);
  // G02[3] = new TLine(70,-40,70,40);
  //
  // for(int i=0; i<4; i++) G02[i]->Draw("SAME");

  c41.Update();

  c42.cd(1);
  fVertexHadronMC[0]->Scale(1/fMCDIS);
  fVertexHadron[0]->Scale(1/fXBjtest);
  fVertexHadronMC[0]->SetLineColor(kBlue);
  fVertexHadron[0]->SetLineColor(kRed);
  fVertexHadron[0]->Draw();
  fVertexHadronMC[0]->Draw("SAME");
  c42.Update();

  c43.cd(1);
  fVertexHadronMC[1]->Scale(1/fMCDIS);
  fVertexHadron[1]->Scale(1/fXBjtest);
  fVertexHadronMC[1]->SetLineColor(kBlue);
  fVertexHadron[1]->SetLineColor(kRed);
  fVertexHadron[1]->Draw();
  fVertexHadronMC[1]->Draw("SAME");
  c43.Update();

  c44.cd(1);
  fVertexHadronMC[2]->Scale(1/fMCDIS);
  fVertexHadron[2]->Scale(1/fXBjtest);
  fVertexHadronMC[2]->SetLineColor(kBlue);
  fVertexHadron[2]->SetLineColor(kRed);
  fVertexHadron[2]->Draw();
  fVertexHadronMC[2]->Draw("SAME");
  c44.Update();

  c1.Print("kinMC.pdf(","pdf");
  c2.Print("kinMC.pdf","pdf");
  c3.Print("kinMC.pdf","pdf");
  c4.Print("kinMC.pdf","pdf");
  c5.Print("kinMC.pdf","pdf");
  c6.Print("kinMC.pdf","pdf");
  c7.Print("kinMC.pdf","pdf");
  c8.Print("kinMC.pdf","pdf");
  c9.Print("kinMC.pdf","pdf");
  c10.Print("kinMC.pdf","pdf");
  c11.Print("kinMC.pdf","pdf");
  c12.Print("kinMC.pdf","pdf");
  c13.Print("kinMC.pdf","pdf");
  c14.Print("kinMC.pdf","pdf");
  c15.Print("kinMC.pdf","pdf");
  c16.Print("kinMC.pdf","pdf");
  c17.Print("kinMC.pdf","pdf");
  c26.Print("kinMC.pdf","pdf");
  c27.Print("kinMC.pdf","pdf");
  c28.Print("kinMC.pdf","pdf");
  c29.Print("kinMC.pdf","pdf");
  c30.Print("kinMC.pdf","pdf");
  c31.Print("kinMC.pdf","pdf");
  c32.Print("kinMC.pdf","pdf");
  c33.Print("kinMC.pdf","pdf");
  c34.Print("kinMC.pdf","pdf");
  c35.Print("kinMC.pdf","pdf");
  c36.Print("kinMC.pdf","pdf");
  c37.Print("kinMC.pdf","pdf");
  c38.Print("kinMC.pdf","pdf");
  c39.Print("kinMC.pdf","pdf");
  c40.Print("kinMC.pdf","pdf");

  c18.Divide(1,1);
  c19.Divide(1,1);
  c20.Divide(1,1);
  c21.Divide(1,1);
  c22.Divide(1,1);
  c23.Divide(1,1);
  c24.Divide(1,1);
  c25.Divide(1,1);
  c18.cd(1);
  fKinematicsMC[0]->Draw();
  gPad->SetLogx();
  c18.Update();
  c19.cd(1);
  fKinematicsMC[1]->Draw();
  gPad->SetLogx();
  c19.Update();
  c20.cd(1);
  fKinematicsMC[2]->Draw();
  c20.Update();
  c21.cd(1);
  fKinematicsMC[3]->Draw();
  c21.Update();
  c22.cd(1);
  fKinematicsMC[4]->Draw();
  c22.Update();
  c23.cd(1);
  fKinematicsMC[5]->Draw();
  c23.Update();
  c24.cd(1);
  fKinematics2DMC->Draw("COLZ");
  gPad->SetLogx();
  c24.Update();
  c25.cd(1);
  fTarget2DMC->Draw("COLZ");
  c25.Update();

  c18.Print("kinMC.pdf","pdf");
  c19.Print("kinMC.pdf","pdf");
  c20.Print("kinMC.pdf","pdf");
  c21.Print("kinMC.pdf","pdf");
  c22.Print("kinMC.pdf","pdf");
  c23.Print("kinMC.pdf","pdf");
  c24.Print("kinMC.pdf","pdf");
  c25.Print("kinMC.pdf","pdf");
  c42.Print("kinMC.pdf","pdf");
  c43.Print("kinMC.pdf","pdf");
  c44.Print("kinMC.pdf)","pdf");
  c41.Print("Trigger_Coverage.pdf");
}

int main(int argc, char **argv)
{

  if(argc < 2)
  {
    cout << "ERROR : Not enough arguments." << endl;
    cout << "Asked : at least 1 *** Received : " << argc-1 << endl;
    cout << "./acceptance_split filelist" << endl;

    return 1;
  }

  int kin_flag=0;

  for (int i = 1; i < argc; i++)
  {
    if (string(argv[i]) == "-h")
    {
      cout << FCYN("HELP : available flags :") << endl;
      cout << FCYN("-k") << endl;
      return 0;
    }

    if (string(argv[i]) == "-k")
    {
      kin_flag=1;
    }
  }

  int year=0;

  if(Y2006) year=2006;
  else if(Y2012) year=2012;
  else if(Y2016) year=2016;

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

  if(kin_flag) create_kin_plots();

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

  ofstream trigmaskout("trigmask.dat", std::ofstream::out | std::ofstream::trunc);

  ifstream periods(argv[1]);
  string filelist, periodName;
  int periodBit;
  while(periods >> periodName)
  {
    periods >> periodBit;
    if(!periodBit) continue;

    ifstream list(Form("%s/%s/filelist.txt",data_path,periodName.c_str()));
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

        trigmaskout << int(trigMask->GetLeaf("trigMask")->GetValue()) << endl;

        // -------------------------------------------------------------------------
        // --------- Calculation ---------------------------------------------------
        // -------------------------------------------------------------------------

        Double_t zlab = z->GetLeaf("z")->GetValue();
        Double_t zlab_MC = MC_vz->GetLeaf("MC_vz")->GetValue();
        int zlabbin=-1;
        int zlabbin_MC=-1;

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

        if(Y2016)
        {
          if(!(-311.19<zlab && zlab<-71.19)) continue;

          if(-311.19<=zlab && zlab<-251.19) zlabbin = 0;
          else if(-251.19<=zlab && zlab<-191.19) zlabbin = 1;
          else if(-191.19<=zlab && zlab<-131.19) zlabbin = 2;
          else if(-131.19<=zlab && zlab<=-71.19) zlabbin = 3;
          if(-311.19<=zlab_MC && zlab_MC<-251.19) zlabbin_MC = 0;
          else if(-251.19<=zlab_MC && zlab_MC<-191.19) zlabbin_MC = 1;
          else if(-191.19<=zlab_MC && zlab_MC<-131.19) zlabbin_MC = 2;
          else if(-131.19<=zlab_MC && zlab_MC<=-71.19) zlabbin_MC = 3;
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

        // -----------------------------------------------------------------------
        //  Data -----------------------------------------------------------------
        // -----------------------------------------------------------------------

        fAllDISflag = 0;
        int DIS_rec[3][12];

        // Best Primary Vertex
        fBP++;

        // Reconstructed muon
        if((0<isMuPrim->GetLeaf("isMuPrim")->GetValue()))
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
                if(InTarget(x->GetLeaf("x")->GetValue(),y->GetLeaf("y")->GetValue(),z->GetLeaf("z")->GetValue()))
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
                if(InTarget(x->GetLeaf("x")->GetValue(),y->GetLeaf("y")->GetValue(),z->GetLeaf("z")->GetValue())
                    && (-325<z->GetLeaf("z")->GetValue() && z->GetLeaf("z")->GetValue()<-71))
                {
                  fTarg++;

                  if((beam_chi2->GetLeaf("beam_chi2")->GetValue()<10))
                  {
                    fMuchi2++;

                    // Cells crossing
                    if((cellsCrossed->GetLeaf("cellsCrossed")->GetValue()))
                    {
                      fCell++;

                      if((mu_prim_chi2->GetLeaf("mu_prim_chi2")->GetValue()<10))
                      {
                        fMupchi2++;

                        if((MZfirst->GetLeaf("MZfirst")->GetValue()<350))
                        {
                          fMZfirst++;

                          if((trig&2 || trig&4 || trig&8 || trig&512))
                          {
                            fTrig++;

                            // Q2 cut
                            if((Q2>1))
                            {
                              fQ2test++;

                              // y cut
                              if((0.1<yBj && yBj<0.9))
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
                          fMCDIS++;
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
              if(InTarget(x->GetLeaf("x")->GetValue(),y->GetLeaf("y")->GetValue(),z->GetLeaf("z")->GetValue()))
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
              if(InTarget(MC_vx->GetLeaf("MC_vx")->GetValue(),MC_vy->GetLeaf("MC_vy")->GetValue(),MC_vz->GetLeaf("MC_vz")->GetValue())
                  && (-325<MC_vz->GetLeaf("MC_vz")->GetValue() && MC_vz->GetLeaf("MC_vz")->GetValue()<-71))
              {
                // Q2 cut
                if((Q2_MC>1))
                {
                  // y cut
                  if((0.1<yBj_MC && yBj_MC<0.9))
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

        if(kin_flag)
        {
          fTCx.push_back(MC_TCx->GetLeaf("MC_TCx")->GetValue());
          fTCy.push_back(MC_TCy->GetLeaf("MC_TCy")->GetValue());
          if(fAllDISflag)
          {
            fQ2kin.push_back(Q2);
            fXBjkin.push_back(xBj);
            fYBjkin.push_back(yBj);
            fWBjkin.push_back(sqrt(wBj));
            fNukin.push_back(nu);
            fMu.push_back(E_beam->GetLeaf("E_beam")->GetValue());
            fX.push_back(x->GetLeaf("x")->GetValue());
            fY.push_back(y->GetLeaf("y")->GetValue());

            if(int(trig&2) && !int(trig&4) && !int(trig&8) && !int(trig&512))
            {
              fHM04x.push_back(HM04x->GetLeaf("HM04x")->GetValue());
              fHM04y.push_back(HM04y->GetLeaf("HM04y")->GetValue());
              fHM05x.push_back(HM05x->GetLeaf("HM05x")->GetValue());
              fHM05y.push_back(HM05y->GetLeaf("HM05y")->GetValue());
              fHM04MCx.push_back(MC_HM04x->GetLeaf("MC_HM04x")->GetValue());
              fHM04MCy.push_back(MC_HM04y->GetLeaf("MC_HM04y")->GetValue());
              fHM05MCx.push_back(MC_HM05x->GetLeaf("MC_HM05x")->GetValue());
              fHM05MCy.push_back(MC_HM05y->GetLeaf("MC_HM05y")->GetValue());
            }
            if(int(trig&4) && !int(trig&2) && !int(trig&8) && !int(trig&512))
            {
              fHL04x.push_back(HL04x->GetLeaf("HL04x")->GetValue());
              fHL04y.push_back(HL04y->GetLeaf("HL04y")->GetValue());
              fHL05x.push_back(HL05x->GetLeaf("HL05x")->GetValue());
              fHL05y.push_back(HL05y->GetLeaf("HL05y")->GetValue());
              fHL04MCx.push_back(MC_HL04x->GetLeaf("MC_HL04x")->GetValue());
              fHL04MCy.push_back(MC_HL04y->GetLeaf("MC_HL04y")->GetValue());
              fHL05MCx.push_back(MC_HL05x->GetLeaf("MC_HL05x")->GetValue());
              fHL05MCy.push_back(MC_HL05y->GetLeaf("MC_HL05y")->GetValue());
            }
            if(int(trig&8) && !int(trig&4) && !int(trig&2) && !int(trig&512))
            {
              fHO03x.push_back(HO03x->GetLeaf("HO03x")->GetValue());
              fHO03y.push_back(HO03y->GetLeaf("HO03y")->GetValue());
              fHO04x.push_back(HO04x->GetLeaf("HO04x")->GetValue());
              fHO04y.push_back(HO04y->GetLeaf("HO04y")->GetValue());
              fHO03MCx.push_back(MC_HO03x->GetLeaf("MC_HO03x")->GetValue());
              fHO03MCy.push_back(MC_HO03y->GetLeaf("MC_HO03y")->GetValue());
              fHO04MCx.push_back(MC_HO04x->GetLeaf("MC_HO04x")->GetValue());
              fHO04MCy.push_back(MC_HO04y->GetLeaf("MC_HO04y")->GetValue());
            }
            if(int(trig&512) && !int(trig&4) && !int(trig&8) && !int(trig&2))
            {
              fHG01x.push_back(HG01x->GetLeaf("HG01x")->GetValue());
              fHG01y.push_back(HG01y->GetLeaf("HG01y")->GetValue());
              fHG021x.push_back(HG021x->GetLeaf("HG021x")->GetValue());
              fHG021y.push_back(HG021y->GetLeaf("HG021y")->GetValue());
              fHG022x.push_back(HG022x->GetLeaf("HG022x")->GetValue());
              fHG022y.push_back(HG022y->GetLeaf("HG022y")->GetValue());
              fHG01MCx.push_back(MC_HG01x->GetLeaf("MC_HG01x")->GetValue());
              fHG01MCy.push_back(MC_HG01y->GetLeaf("MC_HG01y")->GetValue());
              fHG021MCx.push_back(MC_HG021x->GetLeaf("MC_HG021x")->GetValue());
              fHG021MCy.push_back(MC_HG021y->GetLeaf("MC_HG021y")->GetValue());
              fHG022MCx.push_back(MC_HG022x->GetLeaf("MC_HG022x")->GetValue());
              fHG022MCy.push_back(MC_HG022y->GetLeaf("MC_HG022y")->GetValue());
            }
          }
          if(fAllDISflag_MC)
          {
            fQ2kinMC.push_back(Q2_MC);
            fXBjkinMC.push_back(xBj_MC);
            fYBjkinMC.push_back(yBj_MC);
            fWBjkinMC.push_back(sqrt(wBj_MC));
            fNukinMC.push_back(nu_MC);
            fMuMC.push_back(MCE0);
            fXMC.push_back(MC_vx->GetLeaf("MC_vx")->GetValue());
            fYMC.push_back(MC_vy->GetLeaf("MC_vy")->GetValue());
          }
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
        else if(0.5<=yBj_MC && yBj_MC<0.7) ybin_MC = 4;
        else ybin_MC = 5;

        if(0.1<yBj && yBj<0.15) ybin = 0;
        else if(0.15<=yBj && yBj<0.2) ybin = 1;
        else if(0.2<=yBj && yBj<0.3) ybin = 2;
        else if(0.3<=yBj && yBj<0.5) ybin = 3;
        else if(0.5<=yBj && yBj<0.7) ybin = 4;
        else ybin = 5;

        if(fAllDISflag_MC)
        {
          // z Binnig

          for(int i=0; i<12; i++)
          {
            fNDIS_evt_MC[0][xbin_MC][ybin_MC][i]++;
            fNDIS_evt_MC[1][xbin_MC][ybin_MC][i]++;
            fNDIS_evt_MC[2][xbin_MC][ybin_MC][i]++;
            fNDIS_evt_MC_zvtx[0][xbin_MC][ybin_MC][i][zlabbin_MC]++;
            fNDIS_evt_MC_zvtx[1][xbin_MC][ybin_MC][i][zlabbin_MC]++;
            fNDIS_evt_MC_zvtx[2][xbin_MC][ybin_MC][i][zlabbin_MC]++;

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
              fNDIS_evt_MC[0][xbin_MC][ybin_MC][i]--;
              fNDIS_evt_MC_zvtx[0][xbin_MC][ybin_MC][i][zlabbin_MC]--;
              DIS_MC[0][i] = 0;
            }
            if(fFlag_MC[1][xbin_MC][ybin_MC][i])
            {
              fNDIS_evt_MC[1][xbin_MC][ybin_MC][i]--;
              fNDIS_evt_MC_zvtx[1][xbin_MC][ybin_MC][i][zlabbin_MC]--;
              DIS_MC[1][i] = 0;
            }
            if(fFlag_MC[2][xbin_MC][ybin_MC][i])
            {
              fNDIS_evt_MC[2][xbin_MC][ybin_MC][i]--;
              fNDIS_evt_MC_zvtx[2][xbin_MC][ybin_MC][i][zlabbin_MC]--;
              DIS_MC[2][i] = 0;
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
            fNDIS_evt_zvtx[0][xbin][ybin][i][zlabbin]++;
            fNDIS_evt_zvtx[1][xbin][ybin][i][zlabbin]++;
            fNDIS_evt_zvtx[2][xbin][ybin][i][zlabbin]++;

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
              fNDIS_evt[0][xbin][ybin][i]--;
              fNDIS_evt_zvtx[0][xbin][ybin][i][zlabbin]--;
              DIS_rec[0][i] = 0;
            }
            if(fFlag[1][xbin][ybin][i])
            {
              fNDIS_evt[1][xbin][ybin][i]--;
              fNDIS_evt_zvtx[1][xbin][ybin][i][zlabbin]--;
              DIS_rec[1][i] = 0;
            }
            if(fFlag[2][xbin][ybin][i])
            {
              fNDIS_evt[2][xbin][ybin][i]--;
              fNDIS_evt_zvtx[2][xbin][ybin][i][zlabbin]--;
              DIS_rec[2][i] = 0;
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

            fVertexHadronMC[0]->Fill(MC_vz->GetLeaf("MC_vz")->GetValue());
            if(fId == 8 || fId == 9) fVertexHadronMC[1]->Fill(MC_vz->GetLeaf("MC_vz")->GetValue());
            else fVertexHadronMC[2]->Fill(MC_vz->GetLeaf("MC_vz")->GetValue());

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
                fGnrt_zvtx[xbin_MC][ybin_MC][zbin_u][zlabbin_MC].tab[1][0][3] += 1;
                fMCHplus++;
                idMCrec[1][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zbin_u));
                prevMCrec[1][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),i));
                hidMCrec[1][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),fId));
                pMCrec[1][3].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),MC_p->GetLeaf("MCHadrons.P")->GetValue(i)));
                zMCrec[1][3].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zBj_MC_unid));
                if(kin_flag)
                {
                  fKinematicsMC[3]->Fill(zBj_MC);
                }
              }
              if(fFlag_MC[0][xbin_MC][ybin_MC][zbin]) continue;
              idMCrec[1][0].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zbin));
              prevMCrec[1][0].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),i));
              hidMCrec[1][0].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),fId));
              pMCrec[1][0].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),MC_p->GetLeaf("MCHadrons.P")->GetValue(i)));
              zMCrec[1][0].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zBj_MC));
              fMCPiplus++;
              fGnrt[xbin_MC][ybin_MC][zbin].tab[1][0][0] += 1;
              fGnrt_zvtx[xbin_MC][ybin_MC][zbin][zlabbin_MC].tab[1][0][0] += 1;
            }
            else if(fId==1)
            {
              if(!fFlag_MC[0][xbin_MC][ybin_MC][zbin_u])
              {
                fGnrt[xbin_MC][ybin_MC][zbin_u].tab[0][0][3] += 1;
                fGnrt_zvtx[xbin_MC][ybin_MC][zbin][zlabbin_MC].tab[0][0][3] += 1;
                fMCHminus++;
                idMCrec[0][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zbin_u));
                prevMCrec[0][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),i));
                hidMCrec[0][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),fId));
                pMCrec[0][3].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),MC_p->GetLeaf("MCHadrons.P")->GetValue(i)));
                zMCrec[0][3].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zBj_MC_unid));
                if(kin_flag)
                {
                  fKinematicsMC[3]->Fill(zBj_MC);
                }
              }
              if(fFlag_MC[0][xbin_MC][ybin_MC][zbin]) continue;
              idMCrec[0][0].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zbin));
              prevMCrec[0][0].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),i));
              hidMCrec[0][0].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),fId));
              pMCrec[0][0].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),MC_p->GetLeaf("MCHadrons.P")->GetValue(i)));
              zMCrec[0][0].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zBj_MC));
              fMCPiminus++;
              fGnrt[xbin_MC][ybin_MC][zbin].tab[0][0][0] += 1;
              fGnrt_zvtx[xbin_MC][ybin_MC][zbin][zlabbin_MC].tab[0][0][0] += 1;
            }
            else if(fId==2)
            {
              if(!fFlag_MC[0][xbin_MC][ybin_MC][zbin_u])
              {
                fGnrt[xbin_MC][ybin_MC][zbin_u].tab[1][0][3] += 1;
                fGnrt_zvtx[xbin_MC][ybin_MC][zbin][zlabbin_MC].tab[1][0][3] += 1;
                fMCHplus++;
                idMCrec[1][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zbin_u));
                prevMCrec[1][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),i));
                hidMCrec[1][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),fId));
                pMCrec[1][3].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),MC_p->GetLeaf("MCHadrons.P")->GetValue(i)));
                zMCrec[1][3].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zBj_MC_unid));
                if(kin_flag)
                {
                  fKinematicsMC[3]->Fill(zBj_MC);
                }
              }
              if(fFlag_MC[1][xbin_MC][ybin_MC][zbin]) continue;
              idMCrec[1][1].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zbin));
              prevMCrec[1][1].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),i));
              hidMCrec[1][1].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),fId));
              pMCrec[1][1].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),MC_p->GetLeaf("MCHadrons.P")->GetValue(i)));
              zMCrec[1][1].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zBj_MC));
              fMCKplus++;
              fGnrt[xbin_MC][ybin_MC][zbin].tab[1][0][1] += 1;
              fGnrt_zvtx[xbin_MC][ybin_MC][zbin][zlabbin_MC].tab[1][0][1] += 1;
            }
            else if(fId==3)
            {
              if(!fFlag_MC[0][xbin_MC][ybin_MC][zbin_u])
              {
                fGnrt[xbin_MC][ybin_MC][zbin_u].tab[0][0][3] += 1;
                fGnrt_zvtx[xbin_MC][ybin_MC][zbin][zlabbin_MC].tab[0][0][3] += 1;
                fMCHminus++;
                idMCrec[0][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zbin_u));
                prevMCrec[0][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),i));
                hidMCrec[0][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),fId));
                pMCrec[0][3].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),MC_p->GetLeaf("MCHadrons.P")->GetValue(i)));
                zMCrec[0][3].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zBj_MC_unid));
                if(kin_flag)
                {
                  fKinematicsMC[3]->Fill(zBj_MC);
                }
              }
              if(fFlag_MC[1][xbin_MC][ybin_MC][zbin]) continue;
              idMCrec[0][1].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zbin));
              prevMCrec[0][1].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),i));
              hidMCrec[0][1].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),fId));
              pMCrec[0][1].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),MC_p->GetLeaf("MCHadrons.P")->GetValue(i)));
              zMCrec[0][1].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zBj_MC));
              fMCKminus++;
              fGnrt[xbin_MC][ybin_MC][zbin].tab[0][0][1] += 1;
              fGnrt_zvtx[xbin_MC][ybin_MC][zbin][zlabbin_MC].tab[0][0][1] += 1;
            }
            else if(fId==4)
            {
              if(!fFlag_MC[0][xbin_MC][ybin_MC][zbin_u])
              {
                fGnrt[xbin_MC][ybin_MC][zbin_u].tab[1][0][3] += 1;
                fGnrt_zvtx[xbin_MC][ybin_MC][zbin][zlabbin_MC].tab[1][0][3] += 1;
                fMCHplus++;
                idMCrec[1][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zbin_u));
                prevMCrec[1][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),i));
                hidMCrec[1][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),fId));
                pMCrec[1][3].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),MC_p->GetLeaf("MCHadrons.P")->GetValue(i)));
                zMCrec[1][3].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zBj_MC_unid));
                if(kin_flag)
                {
                  fKinematicsMC[3]->Fill(zBj_MC);
                }
              }
              if(fFlag_MC[2][xbin_MC][ybin_MC][zbin]) continue;
              idMCrec[1][2].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zbin));
              prevMCrec[1][2].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),i));
              hidMCrec[1][2].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),fId));
              pMCrec[1][2].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),MC_p->GetLeaf("MCHadrons.P")->GetValue(i)));
              zMCrec[1][2].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zBj_MC));
              fMCPplus++;
              fGnrt[xbin_MC][ybin_MC][zbin].tab[1][0][2] += 1;
              fGnrt_zvtx[xbin_MC][ybin_MC][zbin][zlabbin_MC].tab[1][0][2] += 1;
            }
            else if(fId==5)
            {
              if(!fFlag_MC[0][xbin_MC][ybin_MC][zbin_u])
              {
                fGnrt[xbin_MC][ybin_MC][zbin_u].tab[0][0][3] += 1;
                fGnrt_zvtx[xbin_MC][ybin_MC][zbin][zlabbin_MC].tab[0][0][3] += 1;
                fMCHminus++;
                idMCrec[0][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zbin_u));
                prevMCrec[0][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),i));
                hidMCrec[0][3].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),fId));
                pMCrec[0][3].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),MC_p->GetLeaf("MCHadrons.P")->GetValue(i)));
                zMCrec[0][3].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zBj_MC_unid));
                if(kin_flag)
                {
                  fKinematicsMC[3]->Fill(zBj_MC);
                }
              }
              if(fFlag_MC[2][xbin_MC][ybin_MC][zbin]) continue;
              idMCrec[0][2].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zbin));
              prevMCrec[0][2].insert(pair<int,int>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),i));
              pMCrec[0][2].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),MC_p->GetLeaf("MCHadrons.P")->GetValue(i)));
              zMCrec[0][2].insert(pair<int,Double_t>(int(MC_recHadIdx->GetLeaf("MCHadrons.recHadIdx")->GetValue(i)),zBj_MC));
              fMCPminus++;
              fGnrt[xbin_MC][ybin_MC][zbin].tab[0][0][2] += 1;
              fGnrt_zvtx[xbin_MC][ybin_MC][zbin][zlabbin_MC].tab[0][0][2] += 1;
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
              fGnrt_zvtx[xbin_MC][ybin_MC][zbin][zlabbin_MC].tab[1][0][3] += 1;
              if(kin_flag)
              {
                fKinematicsMC[3]->Fill(zBj_MC);
              }
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
              fGnrt_zvtx[xbin_MC][ybin_MC][zbin][zlabbin_MC].tab[0][0][3] += 1;
              if(kin_flag)
              {
                fKinematicsMC[3]->Fill(zBj_MC);
              }
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

            fHadrons++;
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

            // Maximum radiation length cumulated
            if(!(hXX0->GetLeaf("Hadrons.XX0")->GetValue(i) < 15)) continue;
            fXX0test++;

            // Chi2/ndf
            if(!(chi2_hadron->GetLeaf("Hadrons.chi2_hadron")->GetValue(i) < 10)) continue;
            fChi2Hadron++;

            // Zfirst
            if(!(HZfirst->GetLeaf("Hadrons.HZfirst")->GetValue(i)<350)) continue;
            fHZfirst++;

            // Zlast
            if(!(350<HZlast->GetLeaf("Hadrons.HZlast")->GetValue(i))) continue;
            fHZlast++;

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
            fZtest++;

            fVertexHadron[0]->Fill(z->GetLeaf("z")->GetValue());
            if(fId == 8 || fId == 9) fVertexHadron[1]->Fill(z->GetLeaf("z")->GetValue());
            else fVertexHadron[2]->Fill(z->GetLeaf("z")->GetValue());

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
                fRcstr_zvtx[xbin][ybin][zbin_u][zlabbin].tab[1][0][3] += 1;
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
              fRcstr_zvtx[xbin][ybin][zbin_u][zlabbin].tab[1][0][0] += 1;
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
                fRcstr_zvtx[xbin][ybin][zbin_u][zlabbin].tab[0][0][3] += 1;
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
              fRcstr_zvtx[xbin][ybin][zbin_u][zlabbin].tab[0][0][0] += 1;
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
                fRcstr_zvtx[xbin][ybin][zbin_u][zlabbin].tab[1][0][3] += 1;
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
              fRcstr_zvtx[xbin][ybin][zbin_u][zlabbin].tab[1][0][1] += 1;
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
                fRcstr_zvtx[xbin][ybin][zbin_u][zlabbin].tab[0][0][3] += 1;
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
              fRcstr_zvtx[xbin][ybin][zbin_u][zlabbin].tab[0][0][1] += 1;
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
                fRcstr_zvtx[xbin][ybin][zbin_u][zlabbin].tab[1][0][3] += 1;
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
              fRcstr_zvtx[xbin][ybin][zbin_u][zlabbin].tab[1][0][2] += 1;
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
                fRcstr_zvtx[xbin][ybin][zbin_u][zlabbin].tab[0][0][3] += 1;
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
              fRcstr_zvtx[xbin][ybin][zbin_u][zlabbin].tab[0][0][2] += 1;
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
              fRcstr_zvtx[xbin][ybin][zbin_u][zlabbin].tab[1][0][3] += 1;
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
              fRcstr_zvtx[xbin][ybin][zbin_u][zlabbin].tab[0][0][3] += 1;
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
              fRcstr_zvtx[xbin][ybin][zbin_u][zlabbin].tab[0][0][0] += 1;
              fRcstr_zvtx[xbin][ybin][zbin_u][zlabbin].tab[0][0][3] += 1;
            }
            else if(fId==9)
            {
              if(fFlag[0][xbin][ybin][zbin]) continue;
              fHplus++; fPiplus++;
              fRcstr[xbin][ybin][zbin].tab[1][0][0] += 1;
              fRcstr[xbin][ybin][zbin].tab[1][0][3] += 1;
              fRcstr_zvtx[xbin][ybin][zbin_u][zlabbin].tab[1][0][0] += 1;
              fRcstr_zvtx[xbin][ybin][zbin_u][zlabbin].tab[1][0][3] += 1;
            }
            else
            {
              continue;
            }

            if(kin_flag)
            {
              fKinematics[3]->Fill(zBj);
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

    ofstream ofs_h(Form("acceptance/%d/hadron/hadron_%s.txt",year,periodName.c_str()), std::ofstream::out | std::ofstream::trunc);
    ofstream ofs_hzvtx(Form("acceptance/%d/hadron/hadron_zvtx_%s.txt",year,periodName.c_str()), std::ofstream::out | std::ofstream::trunc);
    ofstream ofs_d(Form("acceptance/%d/DIS/DIS_%s.txt",year,periodName.c_str()), std::ofstream::out | std::ofstream::trunc);
    ofstream ofs_dzvtx(Form("acceptance/%d/DIS/DIS_zvtx_%s.txt",year,periodName.c_str()), std::ofstream::out | std::ofstream::trunc);

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
              ofs_d << fNDIS_evt[0][i][j][k] << " " << fNDIS_evt_c[0][i][j][k] << " " << fNDIS_evt_MC[0][i][j][k] << " " <<
                       fNDIS_evt[1][i][j][k] << " " << fNDIS_evt_c[1][i][j][k] << " " << fNDIS_evt_MC[1][i][j][k] << " " <<
                       fNDIS_evt[2][i][j][k] << " " << fNDIS_evt_c[2][i][j][k] << " " << fNDIS_evt_MC[2][i][j][k] << endl;

              ofs_dzvtx << fNDIS_evt_zvtx[0][i][j][k][0] << " " << fNDIS_evt_zvtx[0][i][j][k][1] << " " << fNDIS_evt_zvtx[0][i][j][k][2] << " " << fNDIS_evt_zvtx[0][i][j][k][3] << " " <<
                           fNDIS_evt_MC_zvtx[0][i][j][k][0] << " " << fNDIS_evt_MC_zvtx[0][i][j][k][1] << " " << fNDIS_evt_MC_zvtx[0][i][j][k][2] << " " << fNDIS_evt_MC_zvtx[0][i][j][k][3] << " " <<
                           fNDIS_evt_zvtx[1][i][j][k][0] << " " << fNDIS_evt_zvtx[1][i][j][k][1] << " " << fNDIS_evt_zvtx[1][i][j][k][2] << " " << fNDIS_evt_zvtx[1][i][j][k][3] << " " <<
                           fNDIS_evt_MC_zvtx[1][i][j][k][0] << " " << fNDIS_evt_MC_zvtx[1][i][j][k][1] << " " << fNDIS_evt_MC_zvtx[1][i][j][k][2] << " " << fNDIS_evt_MC_zvtx[1][i][j][k][3] << " " <<
                           fNDIS_evt_zvtx[2][i][j][k][0] << " " << fNDIS_evt_zvtx[2][i][j][k][1] << " " << fNDIS_evt_zvtx[2][i][j][k][2] << " " << fNDIS_evt_zvtx[2][i][j][k][3] << " " <<
                           fNDIS_evt_MC_zvtx[2][i][j][k][0] << " " << fNDIS_evt_MC_zvtx[2][i][j][k][1] << " " << fNDIS_evt_MC_zvtx[2][i][j][k][2] << " " << fNDIS_evt_MC_zvtx[2][i][j][k][3] << endl;
            }

            ofs_h << fRcstr[i][j][k].tab[c][0][0] << " " << fRcstr_c[i][j][k].tab[c][0][0] << " " << fGnrt[i][j][k].tab[c][0][0] << " " <<
                     fRcstr[i][j][k].tab[c][0][1] << " " << fRcstr_c[i][j][k].tab[c][0][1] << " " << fGnrt[i][j][k].tab[c][0][1] << " " <<
                     fRcstr[i][j][k].tab[c][0][2] << " " << fRcstr_c[i][j][k].tab[c][0][2] << " " << fGnrt[i][j][k].tab[c][0][2] << " " <<
                     fRcstr[i][j][k].tab[c][0][3] << " " << fRcstr_c[i][j][k].tab[c][0][3] << " " << fGnrt[i][j][k].tab[c][0][3] << " " << endl;

            ofs_hzvtx << fRcstr_zvtx[i][j][k][0].tab[c][0][0] << " " << fRcstr_zvtx[i][j][k][1].tab[c][0][0] << " " << fRcstr_zvtx[i][j][k][2].tab[c][0][0] << " " << fRcstr_zvtx[i][j][k][3].tab[c][0][0] << " " <<
                         fGnrt_zvtx[i][j][k][0].tab[c][0][0] << " " << fGnrt_zvtx[i][j][k][1].tab[c][0][0] << " " << fGnrt_zvtx[i][j][k][2].tab[c][0][0] << " " << fGnrt_zvtx[i][j][k][3].tab[c][0][0] << " " <<
                         fRcstr_zvtx[i][j][k][0].tab[c][0][1] << " " << fRcstr_zvtx[i][j][k][1].tab[c][0][1] << " " << fRcstr_zvtx[i][j][k][2].tab[c][0][1] << " " << fRcstr_zvtx[i][j][k][3].tab[c][0][1] << " " <<
                         fGnrt_zvtx[i][j][k][0].tab[c][0][1] << " " << fGnrt_zvtx[i][j][k][1].tab[c][0][1] << " " << fGnrt_zvtx[i][j][k][2].tab[c][0][1] << " " << fGnrt_zvtx[i][j][k][3].tab[c][0][1] << " " <<
                         fRcstr_zvtx[i][j][k][0].tab[c][0][2] << " " << fRcstr_zvtx[i][j][k][1].tab[c][0][2] << " " << fRcstr_zvtx[i][j][k][2].tab[c][0][2] << " " << fRcstr_zvtx[i][j][k][3].tab[c][0][2] << " " <<
                         fGnrt_zvtx[i][j][k][0].tab[c][0][2] << " " << fGnrt_zvtx[i][j][k][1].tab[c][0][2] << " " << fGnrt_zvtx[i][j][k][2].tab[c][0][2] << " " << fGnrt_zvtx[i][j][k][3].tab[c][0][2] << " " <<
                         fRcstr_zvtx[i][j][k][0].tab[c][0][3] << " " << fRcstr_zvtx[i][j][k][1].tab[c][0][3] << " " << fRcstr_zvtx[i][j][k][2].tab[c][0][3] << " " << fRcstr_zvtx[i][j][k][3].tab[c][0][3] << " " <<
                         fGnrt_zvtx[i][j][k][0].tab[c][0][3] << " " << fGnrt_zvtx[i][j][k][1].tab[c][0][3] << " " << fGnrt_zvtx[i][j][k][2].tab[c][0][3] << " " << fGnrt_zvtx[i][j][k][3].tab[c][0][3] << " " <<  endl;
          }
        }
      }
    }

    ofs_h.close();
    ofs_d.close();
    ofs_hzvtx.close();
    ofs_dzvtx.close();

    resetValues();
  }

  if(kin_flag)
  {
    for(int i=0; i<int(fTCx.size()); i++)
    {
      fTrigCov->Fill(fTCx[i],fTCy[i]);
    }
    for(int i=0; i<int(fQ2kin.size()); i++)
    {
      fKinematics[0]->Fill(fQ2kin[i]);
      fKinematics[1]->Fill(fXBjkin[i]);
      fKinematics[2]->Fill(fYBjkin[i]);
      fKinematics[4]->Fill(fWBjkin[i]);
      fKinematics[5]->Fill(fNukin[i]);
      fKinematics[6]->Fill(fMu[i]);
      fKinematics2D->Fill(fXBjkin[i],fYBjkin[i]);
      fTarget2D->Fill(fX[i],fY[i]);
    }
    for(int i=0; i<int(fHM04x.size()); i++)
    {
      fHM04->Fill(fHM04x[i],fHM04y[i]);
      fHM05->Fill(fHM05x[i],fHM05y[i]);
      fHM04MC->Fill(fHM04MCx[i],fHM04MCy[i]);
      fHM05MC->Fill(fHM05MCx[i],fHM05MCy[i]);
    }
    for(int i=0; i<int(fHL04x.size()); i++)
    {
      fHL04->Fill(fHL04x[i],fHL04y[i]);
      fHL05->Fill(fHL05x[i],fHL05y[i]);
      fHL04x1D->Fill(fHL04x[i]);
      fHL05x1D->Fill(fHL05x[i]);
      fHL04MCx2D->Fill(fHL04x[i]-fHL04MCx[i],fHL04y[i]-fHL04MCy[i]);
      fHL05MCx2D->Fill(fHL05x[i]-fHL05MCx[i],fHL05y[i]-fHL05MCy[i]);
      fHLx2D->Fill(fHL04x[i],fHL05x[i]);
    }
    for(int i=0; i<int(fHL04MCx.size()); i++)
    {
      fHL04MC->Fill(fHL04MCx[i],fHL04MCy[i]);
      fHL05MC->Fill(fHL05MCx[i],fHL05MCy[i]);
      fHL04MCx1D->Fill(fHL04MCx[i]);
      fHL05MCx1D->Fill(fHL05MCx[i]);
    }
    for(int i=0; i<int(fHO03x.size()); i++)
    {
      fHO03->Fill(fHO03x[i],fHO03y[i]);
      fHO04->Fill(fHO04x[i],fHO04y[i]);
      fHO03MC->Fill(fHO03MCx[i],fHO03MCy[i]);
      fHO04MC->Fill(fHO04MCx[i],fHO04MCy[i]);
    }
    for(int i=0; i<int(fHG01x.size()); i++)
    {
      fHG01->Fill(fHG01x[i],fHG01y[i]);
      fHG021->Fill(fHG021x[i],fHG021y[i]);
      fHG022->Fill(fHG022x[i],fHG022y[i]);
      fHG01MC->Fill(fHG01MCx[i],fHG01MCy[i]);
      fHG021MC->Fill(fHG021MCx[i],fHG021MCy[i]);
      fHG022MC->Fill(fHG022MCx[i],fHG022MCy[i]);
    }
    for(int i=0; i<int(fQ2kinMC.size()); i++)
    {
      fKinematicsMC[0]->Fill(fQ2kinMC[i]);
      fKinematicsMC[1]->Fill(fXBjkinMC[i]);
      fKinematicsMC[2]->Fill(fYBjkinMC[i]);
      fKinematicsMC[4]->Fill(fWBjkinMC[i]);
      fKinematicsMC[5]->Fill(fNukinMC[i]);
      fKinematicsMC[6]->Fill(fMuMC[i]);
      fKinematics2DMC->Fill(fXBjkinMC[i],fYBjkinMC[i]);
      fTarget2DMC->Fill(fXMC[i],fYMC[i]);

    }
    save_kin_plots();
  }

  cout << "\n\n";
  cout << "             ********* Cut flow for Reconstructed DIS events after cuts ********* " << endl;
  cout << "             -------------------------------------------------------------------- " << endl;

  cout << '|' << setw(30) << "Cut" << '|' << setw(15) << "Events" << '|' << setw(15) << "Abs." << '|' << setw(15) << "Rel." << '|' << endl;
  cout << '|' << setw(30) << "Best Primary Vertex" << '|' << setw(15) << fBP << '|' << setw(15) << float(fBP)/float(fBP)*100 << '|' << setw(15) << float(fBP)/float(fBP)*100 << '|' << endl;
  cout << '|' << setw(30) << "Mu' found (0,1,1,30,0)" << '|' << setw(15) << fRmu << '|' << setw(15) << float(fRmu)/float(fBP)*100 << '|' << setw(15) << float(fRmu)/float(fBP)*100 << '|' << endl;
  cout << '|' << setw(30) << "BMS" << '|' << setw(15) << fBMS << '|' << setw(15) << float(fBMS)/float(fBP)*100 << '|' << setw(15) << float(fBMS)/float(fRmu)*100 << '|' << endl;
  cout << '|' << setw(30) << "140 < E_mu < 180" << '|' << setw(15) << fBEC << '|' << setw(15) << float(fBEC)/float(fBP)*100 << '|' << setw(15) << float(fBEC)/float(fBMS)*100 << '|' << endl;
  cout << '|' << setw(30) << "Vertex in Target" << '|' << setw(15) << fTarg << '|' << setw(15) << float(fTarg)/float(fBP)*100 << '|' << setw(15) << float(fTarg)/float(fBEC)*100 << '|' << endl;
  cout << '|' << setw(30) << "Mu chi2/ndf < 10" << '|' << setw(15) << fMuchi2 << '|' << setw(15) << float(fMuchi2)/float(fBP)*100 << '|' << setw(15) << float(fMuchi2)/float(fTarg)*100 << '|' << endl;
  cout << '|' << setw(30) << "Beam tarck X Cell" << '|' << setw(15) << fCell << '|' << setw(15) << float(fCell)/float(fBP)*100 << '|' << setw(15) << float(fCell)/float(fMuchi2)*100 << '|' << endl;
  cout << '|' << setw(30) << "Mu' chi2/ndf < 10" << '|' << setw(15) << fMupchi2 << '|' << setw(15) << float(fMupchi2)/float(fBP)*100 << '|' << setw(15) << float(fMupchi2)/float(fCell)*100 << '|' << endl;
  cout << '|' << setw(30) << "Mu' Zfirst < 350" << '|' << setw(15) << fMZfirst << '|' << setw(15) << float(fMZfirst)/float(fBP)*100 << '|' << setw(15) << float(fMZfirst)/float(fMupchi2)*100 << '|' << endl;
  cout << '|' << setw(30) << "Triggers MT/LT/OT/LAST" << '|' << setw(15) << fTrig << '|' << setw(15) << float(fTrig)/float(fBP)*100 << '|' << setw(15) << float(fTrig)/float(fMZfirst)*100 << '|' << endl;
  cout << '|' << setw(30) << "Q2 > 1" << '|' << setw(15) << fQ2test << '|' << setw(15) << float(fQ2test)/float(fBP)*100 << '|' << setw(15) << float(fQ2test)/float(fTrig)*100 << '|' << endl;
  cout << '|' << setw(30) << "0.1 < y < 0.7" << '|' << setw(15) << fYBjtest << '|' << setw(15) << float(fYBjtest)/float(fBP)*100 << '|' << setw(15) << float(fYBjtest)/float(fQ2test)*100 << '|' << endl;
  cout << '|' << setw(30) << "5 < W < 17" << '|' << setw(15) << fWBjtest << '|' << setw(15) << float(fWBjtest)/float(fBP)*100 << '|' << setw(15) << float(fWBjtest)/float(fYBjtest)*100 << '|' << endl;
  cout << '|' << setw(30) << "0.004 < x < 0.4" << '|' << setw(15) << fXBjtest << '|' << setw(15) << float(fXBjtest)/float(fBP)*100 << '|' << setw(15) << float(fXBjtest)/float(fWBjtest)*100 << '|' << endl;

  cout << "\n\n";
  cout << "             ********* Cut flow for Reconstructed hadrons after cuts ********* " << endl;
  cout << "             ----------------------------------------------------------------- " << endl;

  cout << '|' << setw(30) << "Cut" << '|' << setw(15) << "Events" << '|' << setw(15) << "Abs." << '|' << setw(15) << "Rel." << endl;
  cout << '|' << setw(30) << "Hadrons" << '|' << setw(15) << fHadrons << '|' << setw(15) << float(fHadrons)/float(fHadrons)*100 << '|' << setw(15) << float(fHadrons)/float(fHadrons)*100 << endl;
  cout << '|' << setw(30) << "XX0 < 15" << '|' << setw(15) << fXX0test << '|' << setw(15) << float(fXX0test)/float(fHadrons)*100 << '|' << setw(15) << float(fXX0test)/float(fHadrons)*100 << endl;
  cout << '|' << setw(30) << "Chi2/ndf > 10" << '|' << setw(15) << fChi2Hadron << '|' << setw(15) << float(fChi2Hadron)/float(fHadrons)*100 << '|' << setw(15) << float(fChi2Hadron)/float(fXX0test)*100 << endl;
  cout << '|' << setw(30) << "Zfirst < 350 cm" << '|' << setw(15) << fHZfirst << '|' << setw(15) << float(fHZfirst)/float(fHadrons)*100 << '|' << setw(15) << float(fHZfirst)/float(fChi2Hadron)*100 << endl;
  cout << '|' << setw(30) << "Zlast > 350 cm" << '|' << setw(15) << fHZlast << '|' << setw(15) << float(fHZlast)/float(fHadrons)*100 << '|' << setw(15) << float(fHZlast)/float(fHZfirst)*100 << endl;
  cout << '|' << setw(30) << "12 < p_h < 40" << '|' << setw(15) << fMom << '|' << setw(15) << float(fMom)/float(fHadrons)*100 << '|' << setw(15) << float(fMom)/float(fHZlast)*100 << endl;
  cout << '|' << setw(30) << "0.01 < theta_RICH < 0.12" << '|' << setw(15) << fTRICH << '|' << setw(15) << float(fTRICH)/float(fHadrons)*100 << '|' << setw(15) << float(fTRICH)/float(fMom)*100 << endl;
  cout << '|' << setw(30) << "Rich Pipe" << '|' << setw(15) << fPosRICH << '|' << setw(15) << float(fPosRICH)/float(fHadrons)*100 << '|' << setw(15) << float(fPosRICH)/float(fTRICH)*100 << endl;
  cout << '|' << setw(30) << "0.2 < z < 0.85" << '|' << setw(15) << fZtest << '|' << setw(15) << float(fZtest)/float(fHadrons)*100 << '|' << setw(15) << float(fZtest)/float(fPosRICH)*100 << endl;

  cout << "\n\n";
  cout << "             ********* Hadron Content (h,pi,K,p) ********* " << endl;
  cout << "             --------------------------------------------- " << endl;

  cout << '|' << setw(15) << "Hadron"
       << '|' << setw(15) << "h+" << '|' << setw(15) << "h-"
       << '|' << setw(15) << "pi+" << '|' << setw(15) << "pi-"
       << '|' << setw(15) << "K+" << '|' << setw(15) << "K-"
       << '|' << setw(15) << "p+" << '|' << setw(15) << "p-" << endl;
  cout << '|' << setw(15) << "ID"
       << '|' << setw(15) << fHplus << '|' << setw(15) << fHminus
       << '|' << setw(15) << fPiplus << '|' << setw(15) << fPiminus
       << '|' << setw(15) << fKplus << '|' << setw(15) << fKminus
       << '|' << setw(15) << fPplus << '|' << setw(15) << fPminus << endl;
  cout << '|' << setw(15) << "True ID"
       << '|' << setw(15) << fHplus << '|' << setw(15) << fHminus
       << '|' << setw(15) << fPiplus_true << '|' << setw(15) << fPiminus_true
       << '|' << setw(15) << fKplus_true << '|' << setw(15) << fKminus_true
       << '|' << setw(15) << fPplus_true << '|' << setw(15) << fPminus_true << endl;

  trigmaskout.close();

  return 0;
}

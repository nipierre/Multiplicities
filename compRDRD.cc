#include "compRDRD.h"

//Inputs
#define target_file_2012 "data/target-107924-109081.dat"
#define target_file_2016 "data/target-274508-274901.dat"

// Flags
#define Y2006 0
#define Y2012 0
#define Y2016 1

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
    fKinematicsRD[i][0] = new TH1F(Form("Q2 %d",i), "", 50, -1, 2);
    fKinematicsRD[i][1] = new TH1F(Form("xBj %d",i), "", 50, -3, 0);
    fKinematicsRD[i][2] = new TH1F(Form("y %d",i), "", 50, 0, 1);
    fKinematicsRD[i][3] = new TH1F(Form("z %d",i), "", 50, 0, 1);
    fKinematicsRD[i][4] = new TH1F(Form("W %d",i), "", 50, 2, 18);
    fKinematicsRD[i][5] = new TH1F(Form("nu %d",i), "", 50, 0, 160);
    fKinematicsRD[i][6] = new TH1F(Form("Emu %d",i), "", 50, 140, 180);
    fKinematicsRD[i][7] = new TH1F(Form("Emu' %d",i), "", 50, 0, 160);
    fKinematicsRD[i][8] = new TH1F(Form("theta %d",i), "", 50, 0, 0.05);
    fKinematicsRD[i][9] = new TH1F(Form("phi %d",i), "", 50, -1.7, 1.7);
    fKinematicsRD[i][10] = new TH1F(Form("Vertex %d",i), "", 50, -320, -70);
    fKinematicsRD[i][11] = new TH1F(Form("pT %d",i), "", 50, 0, 3);
    fKinematicsRD[i][12] = new TH1F(Form("phadron+e %d",i), "", 50, 0, 40);
    fKinematicsRD2[i][0] = new TH1F(Form("Q2 Ratio %d",i), "", 50, -1, 2);
    fKinematicsRD2[i][1] = new TH1F(Form("xBj Ratio %d",i), "", 50, -3, 0);
    fKinematicsRD2[i][2] = new TH1F(Form("y Ratio %d",i), "", 50, 0, 1);
    fKinematicsRD2[i][3] = new TH1F(Form("z Ratio %d",i), "", 50, 0, 1);
    fKinematicsRD2[i][4] = new TH1F(Form("W Ratio %d",i), "", 50, 2, 18);
    fKinematicsRD2[i][5] = new TH1F(Form("nu Ratio %d",i), "", 50, 0, 160);
    fKinematicsRD2[i][6] = new TH1F(Form("Emu Ratio %d",i), "", 50, 140, 180);
    fKinematicsRD2[i][7] = new TH1F(Form("Emu' Ratio %d",i), "", 50, 0, 160);
    fKinematicsRD2[i][8] = new TH1F(Form("theta Ratio %d",i), "", 50, 0, 0.05);
    fKinematicsRD2[i][9] = new TH1F(Form("phi Ratio %d",i), "", 50, -1.7, 1.7);
    fKinematicsRD2[i][10] = new TH1F(Form("Vertex Ratio %d",i), "", 50, -320, -70);
    fKinematicsRD2[i][11] = new TH1F(Form("pT Ratio %d",i), "", 50, 0, 3);
    fKinematicsRD2[i][12] = new TH1F(Form("phadron+e Ratio %d",i), "", 50, 0, 40);
    BinLogX(fKinematicsRD[i][0]);
    BinLogX(fKinematicsRD2[i][0]);
    BinLogX(fKinematicsRD[i][1]);
    BinLogX(fKinematicsRD2[i][1]);
  }

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
    l1[12][i] = new TLine(0,0.4+i*0.2,40,0.4+i*0.2);
    l1[11][i] = new TLine(0,0.4+i*0.2,3,0.4+i*0.2);
    for(int j=0; j<13; j++)
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
  fKinematicsRD2[i][j]->Sumw2();
  fCountingRD2[i][j] = fKinematicsRD2[i][j]->GetEntries();
  fCountingRD1[i][j] = fKinematicsRD1[i][j]->GetEntries();
  fKinematicsRD2[i][j]->Scale(1/fKinematicsRD2[i][j]->GetEntries());
  fKinematicsRD[i][j]->Scale(1/fKinematicsRD[i][j]->GetEntries());
  fKinematicsRatio[i][j] = (TH1F*)fKinematicsRD[i][j]->Clone();
  fKinematicsRatio[i][j]->SetStats(0);
  fKinematicsRatio[i][j]->Divide(fKinematicsRD2[i][j]);
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
  fKinematicsRD2[i][j]->SetLineColor(kBlue);
  fKinematicsRD2[i][j]->SetFillColor(kBlue);
  fKinematicsRD2[i][j]->SetStats(0);
  fKinematicsRD2[i][j]->SetMinimum(0.);
  fKinematicsRD2[i][j]->Draw();
  fKinematicsRD2[i][j]->GetXaxis()->SetLabelSize(0.08);
  fKinematicsRD2[i][j]->GetYaxis()->SetLabelSize(0.08);
  fKinematicsRD2[i][j]->GetYaxis()->SetNdivisions(304,kTRUE);
  fKinematicsRD[i][j]->SetLineColor(kRed);
  fKinematicsRD[i][j]->Draw("SAME");
}

void save_kin_plots()
{
  TCanvas c1("c1","",3200,1600);
  TCanvas c2("c2","",3200,1600);
  TCanvas c3("c3","",3200,1600);
  TCanvas c4("c4","",3200,1600);
  TCanvas c5("c5","",3200,1600);
  TCanvas c6("c6","",3200,1600);
  TCanvas c8("c8","",1600,1600);
  TCanvas c9("c9","",1600,1600);
  TCanvas c10("c10","",1600,1600);
  TCanvas c11("c11","",1600,1600);
  TCanvas c12("c12","",1600,1600);
  TCanvas c13("c13","",1600,1600);
  TCanvas c14("c14","",3200,1600);
  TCanvas c15("c15","",3200,1600);
  TCanvas c16("c16","",3200,1600);
  TCanvas c17("c17","",3200,1600);
  TCanvas c18("c18","",3200,1600);
  TCanvas c19("c19","",3200,1600);
  TCanvas c20("c20","",1600,1600);
  TCanvas c21("c21","",3200,1600);
  TCanvas c22("c22","",1600,1600);

  c1.Divide(2,4);
  c2.Divide(2,4);
  c3.Divide(2,4);
  c4.Divide(2,4);
  c5.Divide(2,4);
  c6.Divide(2,4);
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
  c20.Divide(1,2);
  c21.Divide(2,4);
  c22.Divide(1,2);

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
    plotting_ratio(i,11);
    c19.Update();
    c19.cd(i+offset+1);
    plotting_device(i,11);
    c19.Update();

    c21.cd(i+offset+1+2);
    plotting_ratio(i,12);
    c21.Update();
    c21.cd(i+offset+1);
    plotting_device(i,12);
    c21.Update();
  }

  c8.cd(2);
  plotting_ratio(4,0);
  gPad->SetLogx();
  c8.Update();
  c8.cd(1);
  plotting_device(4,0);
  gPad->SetLogx();
  c8.Update();

  c9.cd(2);
  plotting_ratio(4,1);
  gPad->SetLogx();
  c9.Update();
  c9.cd(1);
  plotting_device(4,1);
  gPad->SetLogx();
  c9.Update();

  c10.cd(2);
  plotting_ratio(4,2);
  c10.Update();
  c10.cd(1);
  plotting_device(4,2);
  c10.Update();

  c11.cd(2);
  plotting_ratio(4,3);
  c11.Update();
  c11.cd(1);
  plotting_device(4,3);
  c11.Update();

  c12.cd(2);
  plotting_ratio(4,4);
  c12.Update();
  c12.cd(1);
  plotting_device(4,4);
  c12.Update();

  c13.cd(2);
  plotting_ratio(4,5);
  c13.Update();
  c13.cd(1);
  plotting_device(4,5);
  c13.Update();

  c20.cd(2);
  plotting_ratio(4,11);
  c20.Update();
  c20.cd(1);
  plotting_device(4,11);
  c20.Update();

  c22.cd(2);
  plotting_ratio(4,12);
  c22.Update();
  c22.cd(1);
  plotting_device(4,12);
  c22.Update();

  c1.SaveAs("kinRDRD.pdf(","pdf");
  c2.SaveAs("kinRDRD.pdf","pdf");
  c3.SaveAs("kinRDRD.pdf","pdf");
  c4.SaveAs("kinRDRD.pdf","pdf");
  c5.SaveAs("kinRDRD.pdf","pdf");
  c6.SaveAs("kinRDRD.pdf","pdf");
  c8.SaveAs("kinRDRD.pdf","pdf");
  c9.SaveAs("kinRDRD.pdf","pdf");
  c10.SaveAs("kinRDRD.pdf","pdf");
  c11.SaveAs("kinRDRD.pdf","pdf");
  c12.SaveAs("kinRDRD.pdf","pdf");
  c13.SaveAs("kinRDRD.pdf","pdf");
  c14.SaveAs("kinRDRD.pdf","pdf");
  c15.SaveAs("kinRDRD.pdf","pdf");
  c16.SaveAs("kinRDRD.pdf","pdf");
  c17.SaveAs("kinRDRD.pdf","pdf");
  c18.SaveAs("kinRDRD.pdf","pdf");
  c19.SaveAs("kinRDRD.pdf","pdf");
  c20.SaveAs("kinRDRD.pdf","pdf");
  c21.SaveAs("kinRDRD.pdf","pdf");
  c22.SaveAs("kinRDRD.pdf)","pdf");


}

void RD2extraction(string pFilelist)
{
  //Kinematics
  Double_t Q2 = 0;
  Double_t xBj = 0;
  Double_t yBj = 0;
  Double_t zBj = 0;
  Double_t wBj = 0;
  Double_t nu = 0;

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
        if(!InTarget(x->GetLeaf("x")->GetValue(),y->GetLeaf("y")->GetValue(),z->GetLeaf("z")->GetValue())) continue;
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
      fMuRD2[4].push_back(E_beam->GetLeaf("E_beam")->GetValue());

      // y cut
      if(!(fYmin<yBj && yBj<fYmax)) continue;

      // W cut
      if(!(fWmin<sqrt(wBj) && sqrt(wBj)<fWmax)) continue;

      // x cut
      if(!(fXmin<xBj && xBj<fXmax)) continue;

      // MT
      if(int(trig&2) && !int(trig&4) && !int(trig&8) && !int(trig&512))
      {
        fQ2kinRD2[0].push_back(Q2);
        fXBjkinRD2[0].push_back(xBj);
        fYBjkinRD2[0].push_back(yBj);
        fWBjkinRD2[0].push_back(sqrt(wBj));
        fNukinRD2[0].push_back(nu);
        fMuRD2[0].push_back(E_beam->GetLeaf("E_beam")->GetValue());
        fMupRD2[0].push_back(E_mu_prim->GetLeaf("E_mu_prim")->GetValue());
        fThetaRD2[0].push_back(theta_m);
        fPhiRD2[0].push_back(phi_m);
        fVertexRD2[0].push_back(z->GetLeaf("z")->GetValue());
      }
      // LT
      if(int(trig&4) && !int(trig&2) && !int(trig&8)&& !int(trig&512))
      {
        fQ2kinRD2[1].push_back(Q2);
        fXBjkinRD2[1].push_back(xBj);
        fYBjkinRD2[1].push_back(yBj);
        fWBjkinRD2[1].push_back(sqrt(wBj));
        fNukinRD2[1].push_back(nu);
        fMuRD2[1].push_back(E_beam->GetLeaf("E_beam")->GetValue());
        fMupRD2[1].push_back(E_mu_prim->GetLeaf("E_mu_prim")->GetValue());
        fThetaRD2[1].push_back(theta_m);
        fPhiRD2[1].push_back(phi_m);
        fVertexRD2[1].push_back(z->GetLeaf("z")->GetValue());
      }
      // OT
      if(int(trig&8) && !int(trig&2) && !int(trig&4) && !int(trig&512))
      {
        fQ2kinRD2[2].push_back(Q2);
        fXBjkinRD2[2].push_back(xBj);
        fYBjkinRD2[2].push_back(yBj);
        fWBjkinRD2[2].push_back(sqrt(wBj));
        fNukinRD2[2].push_back(nu);
        fMuRD2[2].push_back(E_beam->GetLeaf("E_beam")->GetValue());
        fMupRD2[2].push_back(E_mu_prim->GetLeaf("E_mu_prim")->GetValue());
        fThetaRD2[2].push_back(theta_m);
        fPhiRD2[2].push_back(phi_m);
        fVertexRD2[2].push_back(z->GetLeaf("z")->GetValue());
      }
      // LAST
      if(int(trig&512) && !int(trig&4) && !int(trig&8) && !int(trig&2))
      {
        fQ2kinRD2[3].push_back(Q2);
        fXBjkinRD2[3].push_back(xBj);
        fYBjkinRD2[3].push_back(yBj);
        fWBjkinRD2[3].push_back(sqrt(wBj));
        fNukinRD2[3].push_back(nu);
        fMuRD2[3].push_back(E_beam->GetLeaf("E_beam")->GetValue());
        fMupRD2[3].push_back(E_mu_prim->GetLeaf("E_mu_prim")->GetValue());
        fThetaRD2[3].push_back(theta_m);
        fPhiRD2[3].push_back(phi_m);
        fVertexRD2[3].push_back(z->GetLeaf("z")->GetValue());
      }

      // ALL TRIGGERS
      // if(trig&2 || trig&4 || trig&8)
      if(int(trig&2) || int(trig&4) || int(trig&8) || int(trig&512))
      {
        fQ2kinRD2[4].push_back(Q2);
        fXBjkinRD2[4].push_back(xBj);
        fYBjkinRD2[4].push_back(yBj);
        fWBjkinRD2[4].push_back(sqrt(wBj));
        fNukinRD2[4].push_back(nu);
        // fMu[4].push_back(E_beam->GetLeaf("E_beam")->GetValue());
        fMupRD2[4].push_back(E_mu_prim->GetLeaf("E_mu_prim")->GetValue());
        fThetaRD2[4].push_back(theta_m);
        fPhiRD2[4].push_back(phi_m);
        fVertexRD2[4].push_back(z->GetLeaf("z")->GetValue());
      }
      fXRD2.push_back(x->GetLeaf("x")->GetValue());
      fYRD2.push_back(y->GetLeaf("y")->GetValue());


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
        // if(0.1<zBj && (LH->GetLeaf("Hadrons.LH")->GetValue(3+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(4+6*i)) && (LH->GetLeaf("Hadrons.LH")->GetValue(3+6*i)>fLHsec_tab[3]))
        //   fKinematicsRD2[0][11]->Fill(abs(ph_pl->GetLeaf("Hadrons.ph_pl")->GetValue(i)));

        // Maximum radiation length cumulated
        if(!(hXX0->GetLeaf("Hadrons.XX0")->GetValue(i) < 15)) continue;

        // Chi2/ndf
        if(!(chi2_hadron->GetLeaf("Hadrons.chi2_hadron")->GetValue(i) < 10)) continue;

        // Zfirst
        if(!(HZfirst->GetLeaf("Hadrons.HZfirst")->GetValue(i)<350)) continue;

        // Zlast
        if(!(350<HZlast->GetLeaf("Hadrons.HZlast")->GetValue(i))) continue;

        // Theta cut
        if(!(0.01<thRICH->GetLeaf("Hadrons.thRICH")->GetValue(i) && thRICH->GetLeaf("Hadrons.thRICH")->GetValue(i)<0.12)) continue;

        // RICH position cut
        if(!(pow(RICHx->GetLeaf("Hadrons.RICHx")->GetValue(i),2)+pow(RICHy->GetLeaf("Hadrons.RICHy")->GetValue(i),2)>25)) continue;

        // Momentum cut (12 GeV to 40 GeV, increasing to 3 GeV to 40 GeV)
        if(!(fPmin<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<fPmax)) continue;

        // z cut
        if(!(0.2<zBj && zBj<0.85)) continue;

        // Non null charge
        if(!charge->GetLeaf("Hadrons.charge")->GetValue(i)) continue;

        if(int(trig&2) && !int(trig&4) && !int(trig&8) && !int(trig&512))
        {
          fKinematicsRD2[0][3]->Fill(zBj);
          fKinematicsRD2[0][12]->Fill(p->GetLeaf("Hadrons.P")->GetValue(i));
          fKinematicsRD2[0][11]->Fill(pt->GetLeaf("Hadrons.pt")->GetValue(i));
        }
        if(int(trig&4) && !int(trig&2) && !int(trig&8)&& !int(trig&512))
        {
          fKinematicsRD2[1][3]->Fill(zBj);
          fKinematicsRD2[1][12]->Fill(p->GetLeaf("Hadrons.P")->GetValue(i));
          fKinematicsRD2[1][11]->Fill(pt->GetLeaf("Hadrons.pt")->GetValue(i));
        }
        if(int(trig&8) && !int(trig&2) && !int(trig&4) && !int(trig&512))
        {
          fKinematicsRD2[2][3]->Fill(zBj);
          fKinematicsRD2[2][12]->Fill(p->GetLeaf("Hadrons.P")->GetValue(i));
          fKinematicsRD2[2][11]->Fill(pt->GetLeaf("Hadrons.pt")->GetValue(i));
        }
        if(int(trig&512) && !int(trig&4) && !int(trig&8) && !int(trig&2))
        {
          fKinematicsRD2[3][3]->Fill(zBj);
          fKinematicsRD2[3][12]->Fill(p->GetLeaf("Hadrons.P")->GetValue(i));
          fKinematicsRD2[3][11]->Fill(pt->GetLeaf("Hadrons.pt")->GetValue(i));
        }

        // if(trig&2 || trig&4 || trig&8)
        if(int(trig&2) || int(trig&4) || int(trig&8) || int(trig&512))
        {
          fKinematicsRD2[4][3]->Fill(zBj);
          fKinematicsRD2[4][12]->Fill(p->GetLeaf("Hadrons.P")->GetValue(i));
          fKinematicsRD2[4][11]->Fill(pt->GetLeaf("Hadrons.pt")->GetValue(i));
        }
      }

      //Misc
      fQ2RD2.push_back(Q2);
      fXBjRD2.push_back(xBj);
      fYBjRD2.push_back(yBj);
      fWBjRD2.push_back(wBj);
      fNuRD2.push_back(nu);
    }

    cout << "\n" << endl;
  }

  for(int i=0; i<int(fQ2kinRD2[0].size()); i++)
  {
      fKinematicsRD2[0][0]->Fill(fQ2kinRD2[0][i]);
      fKinematicsRD2[0][1]->Fill(fXBjkinRD2[0][i]);
      fKinematicsRD2[0][2]->Fill(fYBjkinRD2[0][i]);
      fKinematicsRD2[0][4]->Fill(fWBjkinRD2[0][i]);
      fKinematicsRD2[0][5]->Fill(fNukinRD2[0][i]);
      fKinematicsRD2[0][6]->Fill(fMuRD2[0][i]);
      fKinematicsRD2[0][7]->Fill(fMupRD2[0][i]);
      fKinematicsRD2[0][8]->Fill(fThetaRD2[0][i]);
      fKinematicsRD2[0][9]->Fill(fPhiRD2[0][i]);
      fKinematicsRD2[0][10]->Fill(fVertexRD2[0][i]);
  }
  for(int i=0; i<int(fQ2kinRD2[1].size()); i++)
  {
      fKinematicsRD2[1][0]->Fill(fQ2kinRD2[1][i]);
      fKinematicsRD2[1][1]->Fill(fXBjkinRD2[1][i]);
      fKinematicsRD2[1][2]->Fill(fYBjkinRD2[1][i]);
      fKinematicsRD2[1][4]->Fill(fWBjkinRD2[1][i]);
      fKinematicsRD2[1][5]->Fill(fNukinRD2[1][i]);
      fKinematicsRD2[1][6]->Fill(fMuRD2[1][i]);
      fKinematicsRD2[1][7]->Fill(fMupRD2[1][i]);
      fKinematicsRD2[1][8]->Fill(fThetaRD2[1][i]);
      fKinematicsRD2[1][9]->Fill(fPhiRD2[1][i]);
      fKinematicsRD2[1][10]->Fill(fVertexRD2[1][i]);
  }
  for(int i=0; i<int(fQ2kinRD2[2].size()); i++)
  {
      fKinematicsRD2[2][0]->Fill(fQ2kinRD2[2][i]);
      fKinematicsRD2[2][1]->Fill(fXBjkinRD2[2][i]);
      fKinematicsRD2[2][2]->Fill(fYBjkinRD2[2][i]);
      fKinematicsRD2[2][4]->Fill(fWBjkinRD2[2][i]);
      fKinematicsRD2[2][5]->Fill(fNukinRD2[2][i]);
      fKinematicsRD2[2][6]->Fill(fMuRD2[2][i]);
      fKinematicsRD2[2][7]->Fill(fMupRD2[2][i]);
      fKinematicsRD2[2][8]->Fill(fThetaRD2[2][i]);
      fKinematicsRD2[2][9]->Fill(fPhiRD2[2][i]);
      fKinematicsRD2[2][10]->Fill(fVertexRD2[2][i]);
  }
  for(int i=0; i<int(fQ2kinRD2[3].size()); i++)
  {
      fKinematicsRD2[3][0]->Fill(fQ2kinRD2[3][i]);
      fKinematicsRD2[3][1]->Fill(fXBjkinRD2[3][i]);
      fKinematicsRD2[3][2]->Fill(fYBjkinRD2[3][i]);
      fKinematicsRD2[3][4]->Fill(fWBjkinRD2[3][i]);
      fKinematicsRD2[3][5]->Fill(fNukinRD2[3][i]);
      fKinematicsRD2[3][6]->Fill(fMuRD2[3][i]);
      fKinematicsRD2[3][7]->Fill(fMupRD2[3][i]);
      fKinematicsRD2[3][8]->Fill(fThetaRD2[3][i]);
      fKinematicsRD2[3][9]->Fill(fPhiRD2[3][i]);
      fKinematicsRD2[3][10]->Fill(fVertexRD2[3][i]);
  }
  for(int i=0; i<int(fQ2kinRD2[4].size()); i++)
  {
      fKinematicsRD2[4][0]->Fill(fQ2kinRD2[4][i]);
      fKinematicsRD2[4][1]->Fill(fXBjkinRD2[4][i]);
      fKinematicsRD2[4][2]->Fill(fYBjkinRD2[4][i]);
      fKinematicsRD2[4][4]->Fill(fWBjkinRD2[4][i]);
      fKinematicsRD2[4][5]->Fill(fNukinRD2[4][i]);
      // fKinematicsRD2[4][6]->Fill(fMu[4][i]);
      fKinematicsRD2[4][7]->Fill(fMupRD2[4][i]);
      fKinematicsRD2[4][8]->Fill(fThetaRD2[4][i]);
      fKinematicsRD2[4][9]->Fill(fPhiRD2[4][i]);
      fKinematicsRD2[4][10]->Fill(fVertexRD2[4][i]);
  }
  for(int i=0; i<int(fMuRD2[4].size()); i++)
  {
    fKinematicsRD2[4][6]->Fill(fMuRD2[4][i]);
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
        if(!InTarget(x->GetLeaf("x")->GetValue(),y->GetLeaf("y")->GetValue(),z->GetLeaf("z")->GetValue())) continue;
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

      // y cut
      if(!(fYmin<yBj && yBj<fYmax)) continue;

      // W cut
      if(!(fWmin<sqrt(wBj) && sqrt(wBj)<fWmax)) continue;

      // x cut
      if(!(fXmin<xBj && xBj<fXmax)) continue;

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
        // if(0.1<zBj && (LH->GetLeaf("Hadrons.LH")->GetValue(3+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(4+6*i)) && (LH->GetLeaf("Hadrons.LH")->GetValue(3+6*i)>fLHsec_tab[3]))
        //   fKinematicsRD[0][11]->Fill(abs(ph_pl->GetLeaf("Hadrons.ph_pl")->GetValue(i)));

        // Maximum radiation length cumulated
        if(!(hXX0->GetLeaf("Hadrons.XX0")->GetValue(i) < 15)) continue;

        // Chi2/ndf
        if(!(chi2_hadron->GetLeaf("Hadrons.chi2_hadron")->GetValue(i) < 10)) continue;

        // Zfirst
        if(!(HZfirst->GetLeaf("Hadrons.HZfirst")->GetValue(i)<350)) continue;

        // Zlast
        if(!(350<HZlast->GetLeaf("Hadrons.HZlast")->GetValue(i))) continue;

        // Theta cut
        if(!(0.01<thRICH->GetLeaf("Hadrons.thRICH")->GetValue(i) && thRICH->GetLeaf("Hadrons.thRICH")->GetValue(i)<0.12)) continue;

        // RICH position cut
        if(!(pow(RICHx->GetLeaf("Hadrons.RICHx")->GetValue(i),2)+pow(RICHy->GetLeaf("Hadrons.RICHy")->GetValue(i),2)>25)) continue;

        // Momentum cut (12 GeV to 40 GeV, increasing to 3 GeV to 40 GeV)
        if(!(fPmin<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<fPmax)) continue;

        // Non null charge
        if(!charge->GetLeaf("Hadrons.charge")->GetValue(i)) continue;

        // z cut
        if(!(0.2<zBj && zBj<0.85)) continue;

        if(int(trig&2) && !int(trig&4) && !int(trig&8) && !int(trig&512))
        {
          fKinematicsRD[0][3]->Fill(zBj);
          fKinematicsRD[0][12]->Fill(p->GetLeaf("Hadrons.P")->GetValue(i));
          fKinematicsRD[0][11]->Fill(pt->GetLeaf("Hadrons.pt")->GetValue(i));
        }
        if(int(trig&4) && !int(trig&2) && !int(trig&8)&& !int(trig&512))
        {
          fKinematicsRD[1][3]->Fill(zBj);
          fKinematicsRD[1][12]->Fill(p->GetLeaf("Hadrons.P")->GetValue(i));
          fKinematicsRD[1][11]->Fill(pt->GetLeaf("Hadrons.pt")->GetValue(i));
        }
        if(int(trig&8) && !int(trig&2) && !int(trig&4) && !int(trig&512))
        {
          fKinematicsRD[2][3]->Fill(zBj);
          fKinematicsRD[2][12]->Fill(p->GetLeaf("Hadrons.P")->GetValue(i));
          fKinematicsRD[2][11]->Fill(pt->GetLeaf("Hadrons.pt")->GetValue(i));
        }
        if(int(trig&512) && !int(trig&4) && !int(trig&8) && !int(trig&2))
        {
          fKinematicsRD[3][3]->Fill(zBj);
          fKinematicsRD[3][12]->Fill(p->GetLeaf("Hadrons.P")->GetValue(i));
          fKinematicsRD[3][11]->Fill(pt->GetLeaf("Hadrons.pt")->GetValue(i));
        }

        // if(trig&2 || trig&4 || trig&8)
        if(int(trig&2) || int(trig&4) || int(trig&8) || int(trig&512))
        {
          fKinematicsRD[4][3]->Fill(zBj);
          fKinematicsRD[4][12]->Fill(p->GetLeaf("Hadrons.P")->GetValue(i));
          fKinematicsRD[4][11]->Fill(pt->GetLeaf("Hadrons.pt")->GetValue(i));
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

  if(Y2012) InitTargetFile(target_file_2012);
  else if(Y2016) InitTargetFile(target_file_2016);
  create_kin_plots();
  readKinCuts(argv[3]);
  cout << "... Real Data 1 treatment ..." << endl;
  RDextraction(argv[1]);
  cout << "... Real Data 2 treatment ..." << endl;
  RD2extraction(argv[2]);
  cout << "... Saving plots ..." << endl;
  save_kin_plots();

  cout << "\n\n";
  cout << "             ********* Event distribution within MT/LT/OT/LAST in percentage of total ********* " << endl;
  cout << "             ---------------------------------------------------------------------------------- " << endl;

  cout << "\n ==> Real Data 1 <==" << endl;

  cout <<  '|' << setw(15) << "All" << '|' << setw(15) << "MT" << '|' << setw(15) << "LT" << '|' << setw(15) << "OT" << '|' << setw(15) << "LAST" << '|' << endl;
  for(int i=0; i<13; i++)
  {
    cout <<  '|' << setw(15) << fCountingRD1[4][0] << '|' << setw(15) << float(fCountingRD1[0][0])/float(fCountingRD1[4][0])*100
                                                                      << '|' << setw(15) << float(fCountingRD1[1][0])/float(fCountingRD1[4][0])*100
                                                                      << '|' << setw(15) << float(fCountingRD1[2][0])/float(fCountingRD1[4][0])*100
                                                                      << '|' << setw(15) << float(fCountingRD1[3][0])/float(fCountingRD1[4][0])*100 << '|' << endl;
  }

  cout << "\n ==> Real Data 2 <==" << endl;

  cout <<  '|' << setw(15) << "All" << '|' << setw(15) << "MT" << '|' << setw(15) << "LT" << '|' << setw(15) << "OT" << '|' << setw(15) << "LAST" << '|' << endl;
  for(int i=0; i<13; i++)
  {
    cout <<  '|' << setw(15) << fCountingRD2[4][i] << '|' << setw(15) << float(fCountingRD2[0][0])/float(fCountingRD2[4][0])*100
                                                                      << '|' << setw(15) << float(fCountingRD2[1][0])/float(fCountingRD2[4][0])*100
                                                                      << '|' << setw(15) << float(fCountingRD2[2][0])/float(fCountingRD2[4][0])*100
                                                                      << '|' << setw(15) << float(fCountingRD2[3][0])/float(fCountingRD2[4][0])*100 << '|' << endl;
  }


  return 0;
}

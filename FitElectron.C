void FitElectron(double x1, double x2, double y1, double y2)
{
  TFile *_file0 = TFile::Open("RDMC.root");
  TTree *Hadron = (TTree*) _file0->Get("Hadron");

  double par[7] = {0,0,0,0,0,0.95,0.1};

  TH1F *h1 = new TH1F();
  Hadron->Draw("(EECAL1+EECAL2)/phad>>h1",Form("(PID==0 || PID==1) && isinECAL==1 && 0.1<(EECAL1+EECAL2)/phad && (EECAL1+EECAL2)/phad<1.25 && %f<xh && xh<%f && %f<yh && yh<%f",x1,x2,y1,y2));
  h1 = (TH1F*) gDirectory->Get("h1");

  TF1 *f1 = new TF1("f1","pol3+gaus(4)",0.3,1.25);
  f1->SetParameters(par);

  h1->Fit("f1","R");
}

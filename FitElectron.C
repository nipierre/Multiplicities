void FitElectron(double x1, double x2, double y1, double y2)
{
  TFile *_file0 = TFile::Open("RDMC.root");
  TTree *Hadron = (TTree*) _file0->Get("Hadron");

  TH1F *h1 = new TH1F();
  Hadron->Draw("(EECAL1+EECAL2)/phad>>h1",Form("(PID==0 || PID==1) && isinECAL==1 && 0.1<(EECAL1+EECAL2)/phad && (EECAL1+EECAL2)/phad<1.25 && %f<xh && xh<%f && %f<yh && yh<%f",x1,x2,y1,y2));
  h1 = (TH1F*) gDirectory->Get("h1");

  TF1 *f1 = new TF1("f1","pol3+gaus(4)",0.3,1.25);

  double gausCenter = f1->GetParameter(5);
  while(!(0.9<gausCenter && gausCenter<1))
  {
    f1->SetParameter(5,0.95);
    f1->SetParameter(6,0.1);

    h1->Fit("f1","R");
  }
}

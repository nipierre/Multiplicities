double FitElectron(double x1, double x2)
{
  TFile *_file0 = TFile::Open("RDMC.root");
  TTree *Hadron = (TTree*) _file0->Get("Hadron");

  TH1F *h1 = new TH1F();
  Hadron->Draw("(EECAL1+EECAL2)/phad>>h1",Form("(PID==0 || PID==1) && isinECAL==1 && %f<xh && xh<%f",x1,x2));
  h1 = (TH1F*) gDirectory->Get("h1");

  TF1 *f1 = new TF1("f1","pol3+gaus(4)",0.3,1.25);

  h1->Fit("f1","R");

  f1->FixParameter(5,0.95);
  f1->FixParameter(6,0.1);

  h1->Fit("f1","R");

  f1->ReleaseParameter(5);
  f1->ReleaseParameter(6);
  f1->SetParameter(5,0.95);
  f1->SetParameter(6,0.1);

  h1->Fit("f1","R");

  double sigma = f1->GetParameter(6);
  double norm = f1->GetParameter(4);
  double integral = 2*sigma*norm*sqrt(TMath::Pi());
  double econt = integral/h1->GetEntries();

  return econt;
}

double ElectronTable()
{
  double fXrange[10] = {.004,.01,.02,.03,.04,.06,.1,.14,.18,.4};
  double cont[9];

  for(int i=0; i<9; i++)
  {
      cont[i] = FitElectron(fXrange[i], fXrange[i+1]);
  }
  for(int i=0; i<9; i++)
  {
      cout << cont[i] << " ";
  }
}

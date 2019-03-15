double contamination(double x1, double x2)
{

  double RK = 0;
  double picount = 0;
  double piecount = 0;
  double econt = 0;

  TFile *_file0 = TFile::Open("RDMC.root");
  TTree *Hadron = (TTree*) _file0->Get("Hadron");

  TH1F *h1, *h2, *h3, *h4;

  Hadron->Draw("(EECAL1+EECAL2)/phad>>h1",Form("(PID==2 || PID==3) && isinECAL==1 && %f<xh && xh<%f",x1,x2));
  // RK = h1->GetEntries();
  // Hadron->Draw("(EECAL1+EECAL2)/phad>>h2",Form("(PID==2 || PID==3) && isinECAL==1 && 0.8<(EECAL1+EECAL2)/phad && (EECAL1+EECAL2)/phad<1.25 && %f<xh && xh<%f",x1,x2));
  // RK /= h2->GetEntries();
  // Hadron->Draw("(EECAL1+EECAL2)/phad>>h3",Form("(PID==0 || PID==1) && isinECAL==1 && %f<xh && xh<%f",x1,x2));
  // picount = RK*h3->GetEntries();
  // Hadron->Draw("(EECAL1+EECAL2)/phad>>h4",Form("(PID==0 || PID==1) && isinECAL==1 && 0.8<(EECAL1+EECAL2)/phad && (EECAL1+EECAL2)/phad<1.25 && %f<xh && xh<%f",x1,x2));
  // piecount = h4->GetEntries();
  //
  // econt = (piecount - picount)/h3->GetEntries();

  return econt;
}

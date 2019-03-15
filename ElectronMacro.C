double contamination(double x1, double x2, double y1, double y2)
{

  double RK = 0;
  double picount = 0;
  double piecount = 0;
  double econt = 0;

  TFile *_file0 = TFile::Open("RDMC.root");
  TTree *Hadron = (TTree*) _file0->Get("Hadron");

  TH1F *h1 = new TH1F();
  Hadron->Draw("(EECAL1+EECAL2)/phad>>h1",Form("(PID==2 || PID==3) && isinECAL==1 && %f<xh && xh<%f && %f<yh && yh<%f",x1,x2,y1,y2));
  h1 = (TH1F*) gDirectory->Get("h1");

  TH1F *h2 = new TH1F();
  Hadron->Draw("(EECAL1+EECAL2)/phad>>h2",Form("(PID==2 || PID==3) && isinECAL==1 && 0.8<(EECAL1+EECAL2)/phad && (EECAL1+EECAL2)/phad<1.25 && %f<xh && xh<%f && %f<yh && yh<%f",x1,x2,y1,y2));
  h2 = (TH1F*) gDirectory->Get("h2");
  RK = h2->GetEntries()/h1->GetEntries();;

  TH1F *h3 = new TH1F();
  Hadron->Draw("(EECAL1+EECAL2)/phad>>h3",Form("(PID==0 || PID==1) && isinECAL==1 && %f<xh && xh<%f && %f<yh && yh<%f",x1,x2,y1,y2));
  h3 = (TH1F*) gDirectory->Get("h3");
  picount = RK*h3->GetEntries();

  TH1F *h4 = new TH1F();
  Hadron->Draw("(EECAL1+EECAL2)/phad>>h4",Form("(PID==0 || PID==1) && isinECAL==1 && 0.8<(EECAL1+EECAL2)/phad && (EECAL1+EECAL2)/phad<1.25 && %f<xh && xh<%f && %f<yh && yh<%f",x1,x2,y1,y2));
  h4 = (TH1F*) gDirectory->Get("h4");
  piecount = h4->GetEntries();

  econt = (piecount - picount)/h3->GetEntries();

  return econt;
}

void contaminationTable()
{

  double fXrange[10] = {.004,.01,.02,.03,.04,.06,.1,.14,.18,.4};
  double fYrange[6] = {.1,.15,.2,.3,.5,.7};

  for(int i=0; i<9; i++)
  {
    for(int j=0; j<5; j++)
    {
      cout << contamination(fXrange[i], fXrange[i+1], fYrange[j], fYrange[j+1]) << " ";
    }
    cout << endl;
  }
}

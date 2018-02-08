#include "FFPlotter.h"

void whichBins(double px, double py, double pz, int &pxbin, int &pybin, int &pzbin)
{
  pxbin=-1; pybin=-1; pzbin=-1;

  if(0.004<px && px<0.01) pxbin=0;
  else if(0.01<px && px<0.02) pxbin=1;
  else if(0.02<px && px<0.03) pxbin=2;
  else if(0.03<px && px<0.04) pxbin=3;
  else if(0.04<px && px<0.06) pxbin=4;
  else if(0.06<px && px<0.1) pxbin=5;
  else if(0.1<px && px<0.14) pxbin=6;
  else if(0.14<px && px<0.18) pxbin=7;
  else if(0.18<px && px<0.4) pxbin=8;

  if(0.1<py && py<0.15) pybin=0;
  else if(0.15<py && py<0.2) pybin=1;
  else if(0.2<py && py<0.3) pybin=2;
  else if(0.3<py && py<0.5) pybin=3;
  else if(0.5<py && py<0.7) pybin=4;

  if(0.20<pz && pz<0.25) pzbin=0;
  else if(0.25<pz && pz<0.30) pzbin=1;
  else if(0.30<pz && pz<0.35) pzbin=2;
  else if(0.35<pz && pz<0.40) pzbin=3;
  else if(0.40<pz && pz<0.45) pzbin=4;
  else if(0.45<pz && pz<0.50) pzbin=5;
  else if(0.50<pz && pz<0.55) pzbin=6;
  else if(0.55<pz && pz<0.60) pzbin=7;
  else if(0.60<pz && pz<0.65) pzbin=8;
  else if(0.65<pz && pz<0.70) pzbin=9;
  else if(0.70<pz && pz<0.75) pzbin=10;
  else if(0.75<pz && pz<0.85) pzbin=11;
}

void readDataFile(string pF, string option="Pi")
{
  double x, y, z, Q2, Dfav, Dunf;
  int xbin, ybin, zbin;

  cout << "Reading file " << pF << "..." << endl;
  ifstream f(pF);

  if(!f)
  {
    return;
  }

  f >> x;
  do
  {
    f >> y >> Q2 >> z;
// #ifdef DEBUG
    cout << x << " " << y << " " << Q2 << " " << z << endl;
// #endif
    whichBins(x,y,z,xbin,ybin,zbin);
// #ifdef DEBUG
    cout << xbin << " " << ybin << " " << zbin << endl;
// #endif
    fQ2mean[xbin][ybin].push_back(Q2);
    fZ[xbin][ybin].push_back(z);
    if(option=="Pi")
    {
      f >> Dfav >> Dunf;
      cout << Dfav << " " << Dunf << endl;
      fDfav[xbin][ybin].push_back(Dfav*z);
      fDunf[xbin][ybin].push_back(Dunf*z);
    }
  } while(f >> x);
}

double MeanValue(vector<double> pVec)
{
  double mean=0;
  for(int i=0; i<int(pVec.size()); i++)
  {
    mean += pVec[i];
  }
  mean /= pVec.size();

  return mean;
}

void PionPlotter(string pF)
{
  readDataFile(pF,"Pi");
  TCanvas c0("dummy","dummy",3200,1600);
  c0.Print("kinMCRD.pdf(","pdf");
  for(int i=0; i<9; i++)
  {
    for(int j=0; j<5; j++)
    {
      if(int(fZ[i][j].size()))
      {
        double mean;
        mean = MeanValue(fQ2mean[i][j]);
        TCanvas c1("FF","FF",3200,1600);
        c1.Divide(2,1);
        fDfavG[i][j] = new TGraph(int(fZ[i][j].size()),&(fZ[i][j][0]),&(fDfav[i][j][0]));
        fDunfG[i][j] = new TGraph(int(fZ[i][j].size()),&(fZ[i][j][0]),&(fDunf[i][j][0]));
        c1.cd(1);
        fDfavG[i][j]->SetTitle(Form("D^{#Pi}_{fav} @ <Q^{2}>=%f;z;D^{#Pi}_{fav}(z)",mean));
        fDfavG[i][j]->Draw("AC*");
        c1.Update();
        c1.cd(2);
        fDfavG[i][j]->SetTitle(Form("D^{#Pi}_{unf} @ <Q^{2}>=%f;z;D^{#Pi}_{unf}(z)",mean));
        fDunfG[i][j]->Draw("AC*");
        c1.Update();
        c1.Print("kinMCRD.pdf","pdf");
      }
    }
  }
  c0.Print("kinMCRD.pdf)","pdf");
}

int main(int argc, char **argv)
{

  string fileFlag, FileI;

  for (int i = 1; i < argc; i++)
  {
    if (string(argv[i]) == "-h")
      {
        cout << FCYN("HELP : available flags :") << endl;
        cout << FCYN("-pion [Pi file]") << endl;
        cout << FCYN("-pion-next [Pi file]") << endl;
        cout << FCYN("-kaon [K file]") << endl;
        cout << FCYN("-kaon-next [K file]") << endl;
        return 0;
      }
    if(i+1 < argc)
    {
      if (string(argv[i]) == "-pion")
      {
        FileI = argv[i+1];
        fileFlag = "-pion";
        break;
      }
      if (string(argv[i]) == "-pion-next")
      {
        FileI = argv[i+1];
        fileFlag = "-pion-next";
        break;
      }
      if (string(argv[i]) == "-kaon")
      {
        FileI = argv[i+1];
        fileFlag = "-kaon";
        break;
      }
      if (string(argv[i]) == "-kaon-next")
      {
        FileI = argv[i+1];
        fileFlag = "-kaon-next";
        break;
      }
    }
  }

  if(fileFlag!="-pion"
      && fileFlag!="-pion-next"
      && fileFlag!="-kaon"
      && fileFlag!="-kaon-next")
  {
    cout << FRED("ERROR : No valid options specified.") << endl;
    return 1;
  }

  ifstream f1(FileI);

  if(!f1)
  {
      cout << "ERROR : Filename '" << f1 <<
      "' not valid." << endl;
      return 1;
  };

  f1.close();

  if(fileFlag=="-pion")
  {
    PionPlotter(FileI);
  }
  else if(fileFlag=="-pion-next")
  {

  }
  else if(fileFlag=="-kaon")
  {

  }
  else if(fileFlag=="-kaon-next")
  {

  }

  return 0;

}

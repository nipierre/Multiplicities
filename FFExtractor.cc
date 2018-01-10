#include "FFExtractor.h"

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

void readDataFile(string pF, double pMult[9][5][12], int kin_storage=0)
{
  string sdummy;
  double ddummy;
  double x, y, z, Q2;
  int xbin, ybin, zbin;

  cout << "Reading file " << pF << "..." << endl;
  ifstream f(pF);

  if(!f)
  {
    return;
  }

  // for(int i=0; i<15; i++)
  // {
  //   f >> sdummy;
  //   cout << sdummy << " ";
  // }

  f >> x;
  do
  {
    f >> ddummy >> ddummy;
    cout << x << " " << ddummy << " " << ddummy << " ";
    f >> y >> ddummy >> ddummy;
    cout << y << " " << ddummy << " " << ddummy << " ";
    f >> Q2;
    f >> z >> ddummy >> ddummy;
    cout << Q2 << " " << z << " " << ddummy << " " << ddummy << " " << ddummy << " ";
    whichBins(x,y,z,xbin,ybin,zbin);
    if(xbin==-1 || ybin==-1 || zbin==-1) continue;
    // cout << xbin << " " << ybin << " " << zbin << endl;
    f >> pMult[xbin][ybin][zbin] >> ddummy >> ddummy >> ddummy >> ddummy;
    cout << pMult[xbin][ybin][zbin] << " " << ddummy << " " << ddummy << " " << ddummy << " " << ddummy << endl;
    if(kin_storage)
    {
      fQ2[xbin][ybin][zbin] = Q2;
      fX[xbin][ybin][zbin] = x;
      fY[xbin][ybin][zbin] = y;
      fZ[xbin][ybin][zbin] = z;
    }
  } while(f >> x);
}

void PionExtractionDeut(string pf1, string pf2)
{
  const LHAPDF::PDF* basepdf = LHAPDF::mkPDF(fLHGrid);
  const LHAPDF::GridPDF& pdf = * dynamic_cast<const LHAPDF::GridPDF*>(basepdf);
  ofstream ofs_D("Pi_FF_deut.txt", std::ofstream::out | std::ofstream::trunc);

  readDataFile(pf1,fPip_d,1);
  readDataFile(pf2,fPim_d);

  for(int i=0; i<9 ; i++) //x
  {
    for(int j=0; j<5 ; j++) //y
    {
      for(int k=0; k<12 ; k++) //z
      {
        if(!fX[i][j][k]) continue;

        const double u = pdf.xfxQ2(2,fX[i][j][k],fQ2[i][j][k]);
        const double ub = pdf.xfxQ2(-2,fX[i][j][k],fQ2[i][j][k]);
        const double d = pdf.xfxQ2(1,fX[i][j][k],fQ2[i][j][k]);
        const double db = pdf.xfxQ2(-1,fX[i][j][k],fQ2[i][j][k]);
        const double s = pdf.xfxQ2(3,fX[i][j][k],fQ2[i][j][k]);
        const double sb = pdf.xfxQ2(-3,fX[i][j][k],fQ2[i][j][k]);

        const double norm_d = 5*(u+ub+d+db)+2*(s+sb);

        fCoeff2E[0][0] = (4*u*(u+d)+db*(ub+db))/norm_d;
        fCoeff2E[0][1] = (u+d+ub+db+2*(s+sb))/norm_d;

        fCoeff2E[1][0] = (4*ub*(ub+db)+d*(u+d))/norm_d;
        fCoeff2E[1][1] = (u+d+ub+db+2*(s+sb))/norm_d;

        TMatrixD cSol = fCoeff2E.InvertFast();

        fMult2E[0][0] = fPip_d[i][j][k];
        fMult2E[1][0] = fPim_d[i][j][k];

        TMatrixD cRes = cSol*fMult2E;

        fDfav[i][j][k] = cRes[0][0];
        fDunf[i][j][k] = cRes[1][0];

        ofs_D << fX[i][j][k] << " " << fY[i][j][k] << " " << fQ2[i][j][k] << " " << fZ[i][j][k] << " " << fDfav[i][j][k]
        << " " << fDunf[i][j][k] << endl;
      }
    }
  }
  ofs_D.close();
}

void PionExtractionProt(string pf1, string pf2)
{
  const LHAPDF::PDF* basepdf = LHAPDF::mkPDF(fLHGrid);
  const LHAPDF::GridPDF& pdf = * dynamic_cast<const LHAPDF::GridPDF*>(basepdf);
  ofstream ofs_D("Pi_FF_prot.txt", std::ofstream::out | std::ofstream::trunc);

  readDataFile(pf1,fPip_p,1);
  readDataFile(pf2,fPim_p);

  for(int i=0; i<9 ; i++) //x
  {
    for(int j=0; j<5 ; j++) //y
    {
      for(int k=0; k<12 ; k++) //z
      {
        if(!fX[i][j][k]) continue;

        const double u = pdf.xfxQ2(2,fX[i][j][k],fQ2[i][j][k]);
        const double ub = pdf.xfxQ2(-2,fX[i][j][k],fQ2[i][j][k]);
        const double d = pdf.xfxQ2(1,fX[i][j][k],fQ2[i][j][k]);
        const double db = pdf.xfxQ2(-1,fX[i][j][k],fQ2[i][j][k]);
        const double s = pdf.xfxQ2(3,fX[i][j][k],fQ2[i][j][k]);
        const double sb = pdf.xfxQ2(-3,fX[i][j][k],fQ2[i][j][k]);

        const double norm_p = 4*(u+ub)+(d+db)+(s+sb);

        fCoeff2E[0][0] = (4*u+db)/norm_p;
        fCoeff2E[0][1] = (4*ub+d+s+sb)/norm_p;

        fCoeff2E[1][0] = (4*ub+d)/norm_p;
        fCoeff2E[1][1] = (4*u+db+s+sb)/norm_p;

        TMatrixD cSol = fCoeff2E.InvertFast();

        fMult2E[0][0] = fPip_p[i][j][k];
        fMult2E[1][0] = fPim_p[i][j][k];

        TMatrixD cRes = cSol*fMult2E;

        fDfav[i][j][k] = cRes[0][0];
        fDunf[i][j][k] = cRes[1][0];

        ofs_D << fX[i][j][k] << " " << fY[i][j][k] << " " << fQ2[i][j][k] << " " << fZ[i][j][k] << " " << fDfav[i][j][k]
        << " " << fDunf[i][j][k] << endl;
      }
    }
  }
  ofs_D.close();
}

void KaonExtraction3E(string pf1, string pf2, string pf3, string pf4)
{
  const LHAPDF::PDF* basepdf = LHAPDF::mkPDF(fLHGrid);
  const LHAPDF::GridPDF& pdf = * dynamic_cast<const LHAPDF::GridPDF*>(basepdf);
  ofstream ofs_D("K_3FF.txt", std::ofstream::out | std::ofstream::trunc);

  readDataFile(pf1,fKp_p,1);
  readDataFile(pf2,fKm_p);
  readDataFile(pf3,fKp_d);
  readDataFile(pf4,fKm_d);

  for(int i=0; i<9 ; i++) //x
  {
    for(int j=0; j<5 ; j++) //y
    {
      for(int k=0; k<12 ; k++) //z
      {
        if(!fX[i][j][k]) continue;

        const double u = pdf.xfxQ2(2,fX[i][j][k],fQ2[i][j][k]);
        const double ub = pdf.xfxQ2(-2,fX[i][j][k],fQ2[i][j][k]);
        const double d = pdf.xfxQ2(1,fX[i][j][k],fQ2[i][j][k]);
        const double db = pdf.xfxQ2(-1,fX[i][j][k],fQ2[i][j][k]);
        const double s = pdf.xfxQ2(3,fX[i][j][k],fQ2[i][j][k]);
        const double sb = pdf.xfxQ2(-3,fX[i][j][k],fQ2[i][j][k]);

        const double norm_p = 4*(u+ub)+d+db+s+sb;
        const double norm_d = 5*(u+ub+d+db)+2*(s+sb);

        fCoeff[0][0] = 4*u/norm_p;
        fCoeff[0][1] = sb/norm_p;
        fCoeff[0][2] = (d+s+db+4*ub)/norm_p;

        fCoeff[1][0] = 4*ub/norm_p;
        fCoeff[1][1] = s/norm_p;
        fCoeff[1][2] = (d+sb+db+4*u)/norm_p;

        fCoeff[2][0] = 4*(u+ub+d+db)/norm_d;
        fCoeff[2][1] = 2*(s+sb)/norm_d;
        fCoeff[2][2] = (6*(u+ub+d+db)+2*(s+sb))/norm_d;

        TMatrixD cSol = fCoeff.InvertFast();

        fKpm_d[i][j][k] = fKp_d[i][j][k]+fKm_d[i][j][k];

        fMult[0][0] = fKp_p[i][j][k];
        fMult[1][0] = fKm_p[i][j][k];
        fMult[2][0] = fKpm_d[i][j][k];

        TMatrixD cRes = cSol*fMult;

        fDfav[i][j][k] = cRes[0][0];
        fDstr[i][j][k] = cRes[1][0];
        fDunf[i][j][k] = cRes[2][0];

        ofs_D << fX[i][j][k] << " " << fY[i][j][k] << " " << fQ2[i][j][k] << " " << fZ[i][j][k] << " " << fDfav[i][j][k]
        << " " << fDstr[i][j][k] << " " << fDunf[i][j][k] << endl;
      }
    }
  }
  ofs_D.close();
}

void KaonExtraction4E(string pf1, string pf2, string pf3, string pf4)
{
  const LHAPDF::PDF* basepdf = LHAPDF::mkPDF(fLHGrid);
  const LHAPDF::GridPDF& pdf = * dynamic_cast<const LHAPDF::GridPDF*>(basepdf);
  ofstream ofs_D("K_4FF.txt", std::ofstream::out | std::ofstream::trunc);

  readDataFile(pf1,fKp_p,1);
  readDataFile(pf2,fKm_p);
  readDataFile(pf3,fKp_d);
  readDataFile(pf4,fKm_d);

  for(int i=0; i<9 ; i++) //x
  {
    for(int j=0; j<5 ; j++) //y
    {
      for(int k=0; k<12 ; k++) //z
      {
        if(!fX[i][j][k]) continue;

        const double u = pdf.xfxQ2(2,fX[i][j][k],fQ2[i][j][k]);
        const double ub = pdf.xfxQ2(-2,fX[i][j][k],fQ2[i][j][k]);
        const double d = pdf.xfxQ2(1,fX[i][j][k],fQ2[i][j][k]);
        const double db = pdf.xfxQ2(-1,fX[i][j][k],fQ2[i][j][k]);
        const double s = pdf.xfxQ2(3,fX[i][j][k],fQ2[i][j][k]);
        const double sb = pdf.xfxQ2(-3,fX[i][j][k],fQ2[i][j][k]);

        const double norm_p = 4*(u+ub)+d+db+s+sb;
        const double norm_d = 5*(u+ub+d+db)+2*(s+sb);

        fCoeff4E[0][0] = 4*u/norm_p;
        fCoeff4E[0][1] = sb/norm_p;
        fCoeff4E[0][2] = (d+s)/norm_p;
        fCoeff4E[0][3] = (db+4*ub)/norm_p;

        fCoeff4E[1][0] = 4*ub/norm_p;
        fCoeff4E[1][1] = s/norm_p;
        fCoeff4E[1][2] = (db+sb)/norm_p;
        fCoeff4E[1][3] = (d+4*u)/norm_p;

        fCoeff4E[2][0] = 4*(u+d)/norm_d; //Du
        fCoeff4E[2][1] = 2*sb/norm_d; //Dsb
        fCoeff4E[2][2] = (u+d+2*s)/norm_d; //Ds+Dd
        fCoeff4E[2][3] = (5*(ub+db))/norm_d; //Dub+Ddb

        fCoeff4E[3][0] = 4*(ub+db)/norm_d; //Dub
        fCoeff4E[3][1] = 2*s/norm_d; //Ds
        fCoeff4E[3][2] = (ub+db+2*sb)/norm_d; //Dsb+Ddb
        fCoeff4E[3][3] = (5*(u+d))/norm_d; //Du+Dd

        TMatrixD cSol = fCoeff4E.InvertFast();

        fMult4E[0][0] = fKp_p[i][j][k];
        fMult4E[1][0] = fKm_p[i][j][k];
        fMult4E[2][0] = fKp_d[i][j][k];
        fMult4E[3][0] = fKm_d[i][j][k];

        TMatrixD cRes = cSol*fMult4E;

        fDfav[i][j][k] = cRes[0][0];
        fDstr[i][j][k] = cRes[1][0];
        fDunf1[i][j][k] = cRes[2][0];
        fDunf2[i][j][k] = cRes[3][0];

        ofs_D << fX[i][j][k] << " " << fY[i][j][k] << " " << fQ2[i][j][k] << " " << fZ[i][j][k] << " " << fDfav[i][j][k]
        << " " << fDstr[i][j][k] << " " << fDunf1[i][j][k] << " " << fDunf2[i][j][k] << endl;
      }
    }
  }
  ofs_D.close();
}

void createDummyData(string pf1, string pf2)
{
  readDataFile(pf1,fKp_d,1);
  readDataFile(pf2,fKm_d);

  ofstream ofs_p("K+_prot.txt", std::ofstream::out | std::ofstream::trunc);
  ofstream ofs_m("K-_prot.txt", std::ofstream::out | std::ofstream::trunc);

  for(int i=0; i<9 ; i++) //x
  {
    for(int j=0; j<5 ; j++) //y
    {
      for(int k=0; k<12 ; k++) //z
      {
        if(fKp_d[i][j][k])
        {
          fKp_p[i][j][k] = fKp_d[i][j][k]*(1+((-1)^(rand()%2))*(rand()%10)/100);
          ofs_p << fX[i][j][k] << " 0 0 " << fY[i][j][k] << " 0 0 " << fQ2[i][j][k] << " " << fZ[i][j][k] << " 0 0 " << fKp_p[i][j][k] << " 0 0 0 0 " << endl;
        }
        if(fKm_d[i][j][k])
        {
          fKm_p[i][j][k] = fKm_d[i][j][k]*(1+((-1)^(rand()%2))*(rand()%10)/100);
          ofs_m << fX[i][j][k] << " 0 0 " << fY[i][j][k] << " 0 0 " << fQ2[i][j][k] << " " << fZ[i][j][k] << " 0 0 " << fKm_p[i][j][k] << " 0 0 0 0 " << endl;
        }
      }
    }
  }

  ofs_p.close();
  ofs_m.close();
}

int main(int argc, char **argv)
{

  for (int i = 1; i < argc; i++)
  {
    if (string(argv[i]) == "-h")
      {
        cout << FCYN("HELP : available flags :") << endl;
        cout << FCYN("-pion-deut [Pi+ file] [Pi- file]") << endl;
        cout << FCYN("-pion-prot [Pi+ file] [Pi- file]") << endl;
        cout << FCYN("-kaon-3 [K+ prot] [K- prot] [K+ deut] [K- deut]") << endl;
        cout << FCYN("-kaon-4 [K+ prot] [K- prot] [K+ deut] [K- deut]") << endl;
        cout << FCYN("-dummy-data [Mult base file]") << endl;
        return 0;
      }
    if(i+2 < argc)
    {
      if (string(argv[i]) == "-dummy-data")
      {
        MultFileP = argv[i+1];
        MultFileM = argv[i+1];
        fileFlag = "-dummy-data";
        break;
      }
      if (string(argv[i]) == "-pion-deut")
      {
        FilePiP = argv[i+1];
        FilePiM = argv[i+2];
        fileFlag = "-pion-deut";
        break;
      }
      if (string(argv[i]) == "-pion-prot")
      {
        FilePiP = argv[i+1];
        FilePiM = argv[i+2];
        fileFlag = "-pion-prot";
        break;
      }
    }
    if(i+4 < argc)
    {
      if (string(argv[i]) == "-kaon-3")
      {
        FileKPP = argv[i+1];
        FileKMP = argv[i+2];
        FileKPD = argv[i+1];
        FileKMD = argv[i+2];
        fileFlag = "-kaon-3";
        break;
      }
      if (string(argv[i]) == "-kaon-4")
      {
        FileKPP = argv[i+1];
        FileKMP = argv[i+2];
        FileKPD = argv[i+1];
        FileKMD = argv[i+2];
        fileFlag = "-kaon-4";
        break;
      }
    }
  }

  if(fileFlag!="-dummy-data"
      && fileFlag!="-pion-prot"
      && fileFlag!="-pion-deut"
      && fileFlag!="-kaon-3"
      && fileFlag!="-kaon-4")
  {
    cout << "ERROR : No valid options specified." << endl;
    return 1;
  }

  fLHGrid = "MSTW2008lo68cl";

  if(fileFlag=="-dummy-data")
  {
    ifstream f1(MultFileP);
    ifstream f2(MultFileM);

    if(!f1 || !f2)
    {
        cout << "ERROR : Filename '" << f1 <<
        "' or '" << f2 <<
        "' not valid." << endl;
        return 1;
    };

    f1.close(); f2.close();
  }
  else if(fileFlag=="-pion-prot" || fileFlag=="-pion-deut")
  {
    ifstream f1(FilePiP);
    ifstream f2(FilePiM);

    if(!f1 || !f2)
    {
        cout << "ERROR : Filename '" << f1 <<
        "' or '" << f2 <<
        "' not valid." << endl;
        return 1;
    };

    f1.close(); f2.close();
  }
  else if(fileFlag=="-kaon-3" || fileFlag=="-kaon-4")
  {
    ifstream f1(FileKPP);
    ifstream f2(FileKMP);
    ifstream f3(FileKPD);
    ifstream f4(FileKMD);

    if(!f1 || !f2 || !f3 || !f4)
    {
        cout << "ERROR : Filename '" << f1 <<
        "' or '" << f2 <<
        "' or '" << f3 <<
        "' or '" << f4 <<
        "' not valid." << endl;
        return 1;
    };

    f1.close(); f2.close(); f3.close(); f4.close();
  }

  if(fileFlag=="-pion-deut")
  {
    PionExtractionDeut(FilePiP,FilePiM);
  }
  else if(fileFlag=="-kaon-3")
  {
    KaonExtraction3E(FileKPP,FileKMP,FileKPD,FileKMD);
  }
  else if(fileFlag=="-kaon-4")
  {
    KaonExtraction4E(FileKPP,FileKMP,FileKPD,FileKMD);
  }
  else if(fileFlag=="-dummy-data")
  {
    createDummyData(argv[1],argv[2]);
  }

  return 0;
}

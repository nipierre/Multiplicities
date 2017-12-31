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
}

void readDataFile(string pF, double pMult[9][5][10], int kin_storage=0)
{
  string sdummy;
  double ddummy;
  double x, y, z, Q2;
  int xbin, ybin, zbin;

  cout << "Reading file " << pF << "..." << endl;
  ifstream f(pF);

  for(int i=0; i<15; i++) f >> sdummy;

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
    f >> pMult[xbin][ybin][zbin] >> ddummy >> ddummy >> ddummy >> ddummy;
    cout << pMult[xbin][ybin][zbin] << " " << ddummy << " " << ddummy << " " << ddummy << " " << ddummy << endl;
    if(kin_storage)
    {
      fQ2[xbin][ybin][zbin] = Q2;
      fX[xbin][ybin][zbin] = x;
      fY[xbin][ybin][zbin] = y;
      fZ[xbin][ybin][zbin] = z;
    }
  } while(!(f >> x));
}

void PionExtraction(string pf1, string pf2)
{

}

void KaonExtraction3E(string pf1, string pf2, string pf3, string pf4)
{
  const LHAPDF::PDF* basepdf = LHAPDF::mkPDF(fLHGrid);
  const LHAPDF::GridPDF& pdf = * dynamic_cast<const LHAPDF::GridPDF*>(basepdf);

  readDataFile(pf1,fKp_p,1);
  readDataFile(pf2,fKm_p);
  readDataFile(pf3,fKp_d);
  readDataFile(pf4,fKm_d);

  for(int i=0; i<9 ; i++) //x
  {
    for(int j=0; j<5 ; j++) //y
    {
      for(int k=0; k<10 ; k++) //z
      {
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
        fDunf[i][j][k] = cRes[3][0];
      }
    }
  }
}

void KaonExtraction4E(string pf1, string pf2, string pf3, string pf4)
{
  const LHAPDF::PDF* basepdf = LHAPDF::mkPDF(fLHGrid);
  const LHAPDF::GridPDF& pdf = * dynamic_cast<const LHAPDF::GridPDF*>(basepdf);

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
      }
    }
  }
}

void createDummyData(string pf1, string pf2)
{
  readDataFile(pf1,fKp_d,1);
  readDataFile(pf2,fKm_d);

  ofstream ofs_p("K+_prot.txt", std::ofstream::out | std::ofstream::trunc);
  ofstream ofs_m("K-_prot.txt", std::ofstream::out | std::ofstream::trunc);

  ofs_p << "'$x$' '$x$ LOW' '$x$ HIGH' '$y$' '$y$ LOW' '$y$ HIGH' '$<Q^2> (GeV/c)^2$' '$z^{K^{+}}$' '$z^{K^{+}}$ LOW' '$z^{K^{+}}$ HIGH' '$M^{K^{+}}$' 'stat +' 'stat -' 'sys +' 'sys -'" << endl;
  ofs_m << "'$x$' '$x$ LOW' '$x$ HIGH' '$y$' '$y$ LOW' '$y$ HIGH' '$<Q^2> (GeV/c)^2$' '$z^{K^{-}}$' '$z^{K^{-}}$ LOW' '$z^{K^{-}}$ HIGH' '$M^{K^{-}}$' 'stat +' 'stat -' 'sys +' 'sys -'" << endl;

  for(int i=0; i<9 ; i++) //x
  {
    for(int j=0; j<5 ; j++) //y
    {
      for(int k=0; k<12 ; k++) //z
      {
        fKp_p[i][j][k] = fKp_d[i][j][k]*(1+(-1)^(rand()%2)*(rand()%10)/100);
        fKm_p[i][j][k] = fKm_d[i][j][k]*(1+(-1)^(rand()%2)*(rand()%10)/100);
        ofs_p << fX[i][j][k] << " 0 0 " << fY[i][j][k] << " 0 0 " << fQ2[i][j][k] << fZ[i][j][k] << " 0 0 " << fKp_p[i][j][k] << " 0 0 0 0 " << endl;
        ofs_m << fX[i][j][k] << " 0 0 " << fY[i][j][k] << " 0 0 " << fQ2[i][j][k] << fZ[i][j][k] << " 0 0 " << fKm_p[i][j][k] << " 0 0 0 0 " << endl;
      }
    }
  }

  ofs_p.close();
  ofs_m.close();
}

int main(int argc, char **argv)
{

  if(argc!=4 && argc!=6)
  {
    cout << "ERROR : Number of arguments not valid." << endl;
    return 1;
  }

  fLHGrid = "MSTW08";

  if(argc==4)
  {
    ifstream f1(argv[1]);
    ifstream f2(argv[2]);

    if(!f1 || !f2)
    {
        cout << "ERROR : Filename '" << f1 <<
        "' or '" << f2 <<
        "' not valid." << endl;
        return 1;
    };

    f1.close(); f2.close();
  }
  else if(argc==6)
  {
    ifstream f1(argv[1]);
    ifstream f2(argv[2]);
    ifstream f3(argv[3]);
    ifstream f4(argv[4]);

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

  if(atoi(argv[3])==1)
  {
    PionExtraction(argv[1],argv[2]);
  }
  else if(atoi(argv[5])==2)
  {
    KaonExtraction3E(argv[1],argv[2],argv[3],argv[4]);
  }
  else if(atoi(argv[5])==3)
  {
    KaonExtraction4E(argv[1],argv[2],argv[3],argv[4]);
  }
  else if(atoi(argv[3])==4)
  {
    createDummyData(argv[1],argv[2]);
  }
  else
  {
    cout << "Error : No valid option." << endl;
    return 1;
  }

  return 0;
}

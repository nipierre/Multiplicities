#include "analySIDIS_split.h"

using namespace std;

//******************************************************************************
//Inputs
//******************************************************************************

#define mat_RICH_2006_name "data/rich_mat_2006.txt"
#define mat_RICH_2016_name "data/rich_mat_2016.txt"
#define err_RICH_name "data/rich_mat_error.txt"
#define target_file_2012 "data/target-107924-109081.dat"
#define target_file_2016 "data/target-274508-274901.dat"
#define proton_sirc "data/proton_semi_inclusive_RC.txt"
#define proton_irc "data/hh160_r1998_f2tulay_compass_grv.asy_hcorr.txt"
#define ElectronPi "data/electron_pion_contamination.txt"
#define ElectronPiVtx "data/electron_pion_contamination_vtx.txt"
#define ElectronPiTheta "data/electron_pion_contamination_Theta.txt"
#define ElectronPipT "data/electron_pion_contamination_pT.txt"

//******************************************************************************
//User dependant input
//******************************************************************************

#define data_path "/sps/compass/npierre"

//******************************************************************************
// Flags
//******************************************************************************

#define Y2006 0
#define Y2012 0
#define Y2016 1
#define MOMENTUM_DOWN 12
#define MOMENTUM_UP 40
#define XMIN 0.004
#define XMAX 0.4
#define YMIN 0.1
#define YMAX 0.7
#define HXX0LIMIT 15
#define MUCHARGE_SEPARATION 0
#define MUCHARGE 1

#define IRC 0
#define SIRC 0
#define RICH 1

//******************************************************************************
// Progress bar
//******************************************************************************

# define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
# define PBWIDTH 60

void printProgress(int event, int total)
{
    string points[6] = {"   ",".  ",".. ","..."," ..","  ."};
    float percentage = float(event)/float(total);
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf ("\r Progress%s %3d%% [%.*s%*s] (%d/%d)", points[int(event/24)%6].c_str(), val, lpad, PBSTR, rpad, "", event, total);
    fflush (stdout);
}

//******************************************************************************
// Fusion sort
//******************************************************************************

void fusion(Float_t* tab, Int_t beg1 , Int_t end1, Int_t end2)
{
   Float_t* tab2 = new Float_t[end1-beg1+1];
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

void fusionSort2(Float_t* tab, Int_t begin, Int_t end)
{
   if(begin!=end)
   {
      Int_t mid = (begin+end)/2;
      fusionSort2(tab, begin, mid);
      fusionSort2(tab, mid+1, end);
      fusion(tab, begin, mid, end);
   }
}

void fusionSort(Float_t* tab, Int_t len)
{
   if(len>0)
      fusionSort2(tab, 0, len-1);
}

//******************************************************************************
// Target Management
//******************************************************************************

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

void CellCenter(Float_t z, Float_t& xc, Float_t& yc, Float_t& R)
{
  xc = 1000000;
  yc = 1000000;

  for(Int_t i = 0; i < int(fZv.size()-1); i++)
  {
    Float_t z1 = fZv[i];
    Float_t z2 = fZv[i+1];

    if( z2 < z ) continue;
    if( z1 > z ) continue;

    Float_t xc1 = fXv[i];
    Float_t xc2 = fXv[i+1];

    Float_t yc1 = fYv[i];
    Float_t yc2 = fYv[i+1];

    Float_t rc1 = fRv[i];
    Float_t rc2 = fRv[i+1];

    Float_t dxcdz = (xc2-xc1)/(z2-z1);
    Float_t dycdz = (yc2-yc1)/(z2-z1);
    Float_t drcdz = (rc2-rc1)/(z2-z1);

    Float_t dz = z-z1;
    xc = xc1 + dxcdz*dz;
    yc = yc1 + dycdz*dz;
    R = rc1 + drcdz*dz;

    break;
  }
}

bool InTarget(Float_t xvtx, Float_t yvtx, Float_t zvtx)
{
  Float_t xc, yc, R;
  CellCenter(zvtx, xc, yc, R);
  Float_t dx = xvtx-xc;
  Float_t dy = yvtx-yc;
  Float_t r = sqrt(dx*dx + dy*dy);

  return( r < 1.9 && yvtx < 1.2 );
}

void LoadInclusiveRadiativeCorrection()
{
  string sdum;

  ifstream proton(proton_irc);
  for(int i=0; i<19; i++)
  {
    for(int j=0; j<5; j++)
    {
      proton >> sdum;
#ifdef DEBUG
      cout << sdum << "\t";
#endif

      for(int k=0; k<6; k++)
      {
        proton >> fInclusiveRCproton[k+j*6][i] >> sdum;
#ifdef DEBUG
        cout << " " << fInclusiveRCproton[k+j*6][i] << sdum;
#endif
      }

#ifdef DEBUG
      cout << endl;
#endif
    }
  }
  proton.close();
}

//******************************************************************************
// RADIATIVE CORRECTION
//******************************************************************************

void LoadSemiInclusiveRadiativeCorrection()
{
  string sdum;

  ifstream proton(proton_sirc);

  for(int i=0; i<9; i++)
  {
    for(int j=0; j<6; j++)
    {
      for(int k=0; k<14; k++)
      {
        for(int l=0; l<7; l++)
        {
          proton >> sdum;
#ifdef DEBUG
          cout << sdum << "\t";
#endif
        }
        proton >> fSemiInclusiveRCproton[i][j][k];
#ifdef DEBUG
        cout << fSemiInclusiveRCproton[i][j][k] << "\t";
#endif
        proton >> sdum;
#ifdef DEBUG
        cout << sdum << endl;
#endif
      }
    }
  }
  proton.close();
}

Float_t GetInclusiveRadiativeCorrection(Float_t x, Float_t y)
{
  int xb=-1;
  int yb=-1;

  if(0.0<x && x<0.00006) xb = 0;
  else if(0.00006<=x && x<0.00085) xb = 1;
  else if(0.00085<=x && x<0.00015) xb = 2;
  else if(0.00015<=x && x<0.00025) xb = 3;
  else if(0.00025<=x && x<0.0004) xb = 4;
  else if(0.0004<=x && x<0.0006) xb = 5;
  else if(0.0006<=x && x<0.0085) xb = 6;
  else if(0.0085<=x && x<0.0015) xb = 7;
  else if(0.0015<=x && x<0.003) xb = 8;
  else if(0.003<=x && x<0.005) xb = 9;
  else if(0.005<=x && x<0.007) xb = 10;
  else if(0.007<=x && x<0.009) xb = 11;
  else if(0.009<=x && x<0.0115) xb = 12;
  else if(0.0115<=x && x<0.0145) xb = 13;
  else if(0.0145<=x && x<0.018) xb = 14;
  else if(0.018<=x && x<0.025) xb = 15;
  else if(0.025<=x && x<0.035) xb = 16;
  else if(0.035<=x && x<0.05) xb = 17;
  else if(0.05<=x && x<0.07) xb = 18;
  else if(0.07<=x && x<0.09) xb = 19;
  else if(0.09<=x && x<0.125) xb = 20;
  else if(0.125<=x && x<0.175) xb = 21;
  else if(0.175<=x && x<0.25) xb = 22;
  else if(0.25<=x && x<0.35) xb = 23;
  else if(0.35<=x && x<0.45) xb = 24;
  else if(0.45<=x && x<0.55) xb = 25;
  else if(0.55<=x && x<0.65) xb = 26;
  else if(0.65<=x && x<0.75) xb = 27;
  else if(0.75<=x && x<0.85) xb = 28;
  else xb = 29;

  if(0.0<y && y<0.7) yb = 0;
  else if(0.7<=y && y<0.125) yb = 1;
  else if(0.125<=y && y<0.175) yb = 2;
  else if(0.175<=y && y<0.225) yb = 3;
  else if(0.225<=y && y<0.275) yb = 4;
  else if(0.275<=y && y<0.325) yb = 5;
  else if(0.325<=y && y<0.375) yb = 6;
  else if(0.375<=y && y<0.425) yb = 7;
  else if(0.425<=y && y<0.475) yb = 8;
  else if(0.475<=y && y<0.525) yb = 9;
  else if(0.525<=y && y<0.575) yb = 10;
  else if(0.575<=y && y<0.625) yb = 11;
  else if(0.625<=y && y<0.675) yb = 12;
  else if(0.675<=y && y<0.725) yb = 13;
  else if(0.725<=y && y<0.775) yb = 14;
  else if(0.775<=y && y<0.825) yb = 15;
  else if(0.825<=y && y<0.875) yb = 16;
  else if(0.875<=y && y<0.925) yb = 17;
  else yb = 18;

  if(Y2006 || !IRC)
  {
    return 1;
  }
  else if(Y2012 || Y2016)
  {
    return 1;
  }
  else
  {
    cout << "ERROR in GetInclusiveRadiativeCorrection : Year not recognized. No correction applied." << endl;
    return 1;
  }
}

Float_t GetSemiInclusiveRadiativeCorrection(Float_t x, Float_t y, Float_t z)
{
  int xb=-1;
  int yb=-1;
  int zb=-1;

  if(0.004<x && x<0.01) xb = 0;
  else if(0.01<=x && x<0.02) xb = 1;
  else if(0.02<=x && x<0.03) xb = 2;
  else if(0.03<=x && x<0.04) xb = 3;
  else if(0.04<=x && x<0.06) xb = 4;
  else if(0.06<=x && x<0.1) xb = 5;
  else if(0.1<=x && x<0.14) xb = 6;
  else if(0.14<=x && x<0.18) xb = 7;
  else xb = 8;

  if(0.1<y && y<0.15) yb = 0;
  else if(0.15<=y && y<0.2) yb = 1;
  else if(0.2<=y && y<0.3) yb = 2;
  else if(0.3<=y && y<0.5) yb = 3;
  else if(0.5<=y && y<0.7) yb = 4;
  else yb = 5;

  if(0<z && z<0.2) zb = 0;
  else if(0.2<z && z<0.25) zb = 1;
  else if(0.25<=z && z<0.30) zb = 2;
  else if(0.30<=z && z<0.35) zb = 3;
  else if(0.35<=z && z<0.40) zb = 4;
  else if(0.40<=z && z<0.45) zb = 5;
  else if(0.45<=z && z<0.50) zb = 6;
  else if(0.50<=z && z<0.55) zb = 7;
  else if(0.55<=z && z<0.60) zb = 8;
  else if(0.60<=z && z<0.65) zb = 9;
  else if(0.65<=z && z<0.70) zb = 10;
  else if(0.70<=z && z<0.75) zb = 11;
  else if(0.75<=z && z<0.85) zb = 12;
  else zb = 13;

  if(Y2006 || !SIRC)
  {
    return 1;
  }
  else if(Y2012 || Y2016)
  {
    return 1;
  }
  else
  {
    cout << "ERROR in GetSemiInclusiveRadiativeCorrection : Year not recognized. No correction applied." << endl;
    return 1;
  }
}

//******************************************************************************
// ELECTRON CORRECTION
//******************************************************************************

void LoadElectronCorrection()
{
  ifstream epi(ElectronPi);

  for(int c=0; c<2; c++)
  {
    for(int i=0; i<9; i++)
    {
      for(int j=0; j<6; j++)
      {
        for(int k=0; k<12; k++)
        {
          epi >> fCepi[c][1][i][j][k] >> fCepi[c][0][i][j][k];
          // cout << fCepi[c][1][i][j][k] << " " << fCepi[c][0][i][j][k] << endl;
        }
      }
    }
  }

  epi.close();

  ifstream epiVtx(ElectronPiVtx);

  for(int c=0; c<2; c++)
  {
    for(int i=0; i<9; i++)
    {
      for(int j=0; j<6; j++)
      {
        for(int k=0; k<12; k++)
        {
          for(int zv=0; zv<4; zv++)
          {
            epiVtx >> fCepiVtx[c][1][i][j][k][zv] >> fCepiVtx[c][0][i][j][k][zv];
            // cout << fCepiVtx[c][1][i][j][k][zv] << " " << fCepiVtx[c][0][i][j][k][zv] << endl;
          }
        }
      }
    }
  }

  epiVtx.close();

  ifstream epiTheta(ElectronPiTheta);

  for(int c=0; c<2; c++)
  {
    for(int i=0; i<9; i++)
    {
      for(int j=0; j<6; j++)
      {
        for(int k=0; k<12; k++)
        {
          for(int th=0; th<8; th++)
          {
            epiTheta >> fCepiTh[c][1][i][j][k][th] >> fCepiTh[c][0][i][j][k][th];
          }
        }
      }
    }
  }

  epiTheta.close();

  ifstream epipT(ElectronPipT);

  for(int c=0; c<2; c++)
  {
    for(int i=0; i<9; i++)
    {
      for(int j=0; j<6; j++)
      {
        for(int k=0; k<12; k++)
        {
          for(int pt=0; pt<10; pt++)
          {
            epipT >> fCepipT[c][1][i][j][k][pt] >> fCepipT[c][0][i][j][k][pt];
          }
        }
      }
    }
  }

  epipT.close();

}

//******************************************************************************
// RICH MATRICES
//******************************************************************************

void load_rich_mat_2006(string prich, string prich_err)
{

  pi_sigma_uni[0][0] = 1;
  k_sigma_uni[1][1] = 1;
  p_sigma_uni[2][2] = 1;
  pi_vect[0][0] = 1;
  k_vect[1][0] = 1;
  p_vect[2][0] = 1;

  for(int i=0; i<2; i++)
  {
    for(int j=0; j<10; j++)
    {
      rich_mat_p[i][j].ResizeTo(3,3);
      rich_mat_m[i][j].ResizeTo(3,3);
      inv_rich_p[i][j].ResizeTo(3,3);
      inv_rich_m[i][j].ResizeTo(3,3);
    }
  }

  for(int i=0; i<3; i++)
  {
    for(int j=0; j<3; j++)
    {
      err_rich_p[i][j].ResizeTo(3,3);
      err_rich_m[i][j].ResizeTo(3,3);
    }
  }

  ifstream matRICH(prich);

  for(int i=0; i<24; i++)
  {
    if(i < 4)
    {
      for(int j=0; j<20; j++)
        matRICH >> dummy;
    }
    else
    {
      matRICH >> mat_bin[0][i] >> mat_bin[1][i];
      for(int j=0; j<9; j++)
      {
        matRICH >> rich_mat_p[i%2][(i-4)/2][(int)j/3][j%3];
      }
      for(int j=0; j<9; j++)
      {
        matRICH >> rich_mat_m[i%2][(i-4)/2][(int)j/3][j%3];
      }
    }
  }

  matRICH.close();

  for(int i=0; i<10; i++)
  {
    cout << rich_mat_p[0][i][0][0] << " " << rich_mat_p[0][i][0][1] << " " << rich_mat_p[0][i][0][2] << endl;
    cout << rich_mat_p[0][i][1][0] << " " << rich_mat_p[0][i][1][1] << " " << rich_mat_p[0][i][1][2] << endl;
    cout << rich_mat_p[0][i][2][0] << " " << rich_mat_p[0][i][2][1] << " " << rich_mat_p[0][i][2][2] << endl;

    cout << endl;

    cout << rich_mat_p[1][i][0][0] << " " << rich_mat_p[1][i][0][1] << " " << rich_mat_p[1][i][0][2] << endl;
    cout << rich_mat_p[1][i][1][0] << " " << rich_mat_p[1][i][1][1] << " " << rich_mat_p[1][i][1][2] << endl;
    cout << rich_mat_p[1][i][2][0] << " " << rich_mat_p[1][i][2][1] << " " << rich_mat_p[1][i][2][2] << endl;

    cout << endl;

    cout << rich_mat_m[0][i][0][0] << " " << rich_mat_m[0][i][0][1] << " " << rich_mat_m[0][i][0][2] << endl;
    cout << rich_mat_m[0][i][1][0] << " " << rich_mat_m[0][i][1][1] << " " << rich_mat_m[0][i][1][2] << endl;
    cout << rich_mat_m[0][i][2][0] << " " << rich_mat_m[0][i][2][1] << " " << rich_mat_m[0][i][2][2] << endl;

    cout << endl;

    cout << inv_rich_m[1][i][0][0] << " " << inv_rich_m[1][i][0][1] << " " << inv_rich_m[1][i][0][2] << endl;
    cout << inv_rich_m[1][i][1][0] << " " << inv_rich_m[1][i][1][1] << " " << inv_rich_m[1][i][1][2] << endl;
    cout << inv_rich_m[1][i][2][0] << " " << inv_rich_m[1][i][2][1] << " " << inv_rich_m[1][i][2][2] << endl;

    cout << endl;

    inv_rich_p[0][i] = rich_mat_p[0][i].InvertFast();
    inv_rich_p[1][i] = rich_mat_p[1][i].InvertFast();
    inv_rich_m[0][i] = rich_mat_m[0][i].InvertFast();
    inv_rich_m[1][i] = rich_mat_m[1][i].InvertFast();

    cout << inv_rich_p[0][i][0][0] << " " << inv_rich_p[0][i][0][1] << " " << inv_rich_p[0][i][0][2] << endl;
    cout << inv_rich_p[0][i][1][0] << " " << inv_rich_p[0][i][1][1] << " " << inv_rich_p[0][i][1][2] << endl;
    cout << inv_rich_p[0][i][2][0] << " " << inv_rich_p[0][i][2][1] << " " << inv_rich_p[0][i][2][2] << endl;

    cout << endl;


    cout << inv_rich_p[1][i][0][0] << " " << inv_rich_p[1][i][0][1] << " " << inv_rich_p[1][i][0][2] << endl;
    cout << inv_rich_p[1][i][1][0] << " " << inv_rich_p[1][i][1][1] << " " << inv_rich_p[1][i][1][2] << endl;
    cout << inv_rich_p[1][i][2][0] << " " << inv_rich_p[1][i][2][1] << " " << inv_rich_p[1][i][2][2] << endl;

    cout << endl;


    cout << inv_rich_m[0][i][0][0] << " " << inv_rich_m[0][i][0][1] << " " << inv_rich_m[0][i][0][2] << endl;
    cout << inv_rich_m[0][i][1][0] << " " << inv_rich_m[0][i][1][1] << " " << inv_rich_m[0][i][1][2] << endl;
    cout << inv_rich_m[0][i][2][0] << " " << inv_rich_m[0][i][2][1] << " " << inv_rich_m[0][i][2][2] << endl;

    cout << endl;

    cout << rich_mat_m[1][i][0][0] << " " << rich_mat_m[1][i][0][1] << " " << rich_mat_m[1][i][0][2] << endl;
    cout << rich_mat_m[1][i][1][0] << " " << rich_mat_m[1][i][1][1] << " " << rich_mat_m[1][i][1][2] << endl;
    cout << rich_mat_m[1][i][2][0] << " " << rich_mat_m[1][i][2][1] << " " << rich_mat_m[1][i][2][2] << endl;

    cout << endl;

  }

  // Errors YODO

  ifstream errRICH(prich_err);


  for(int loop=0; loop<24; loop++)
  {
    if(loop < 4)
    {
      for(int j=0; j<38; j++)
        errRICH >> dummy;
    }
    else
    {
      errRICH >> err_bin[0][loop-4] >> err_bin[1][loop-4];

      errRICH >> err_rich_p[0][0][0][0];
      errRICH >> err_rich_p[0][1][0][1];
      errRICH >> err_rich_p[0][2][0][2];//3
      errRICH >> err_rich_p[0][0][1][0];
      errRICH >> err_rich_p[0][0][2][0];
      errRICH >> err_rich_p[1][0][2][0];//6
      errRICH >> err_rich_p[1][0][1][0];
      errRICH >> err_rich_p[1][1][1][1];
      errRICH >> err_rich_p[1][2][1][2];//9
      errRICH >> err_rich_p[0][1][1][1];
      errRICH >> err_rich_p[0][1][2][1];
      errRICH >> err_rich_p[1][1][2][1];//12
      errRICH >> err_rich_p[2][0][2][0];
      errRICH >> err_rich_p[2][1][2][1];
      errRICH >> err_rich_p[2][2][2][2];//15
      errRICH >> err_rich_p[0][2][1][2];
      errRICH >> err_rich_p[0][2][2][2];
      errRICH >> err_rich_p[1][2][2][2];//18

      err_rich_p[0][0][0][0] = pow(err_rich_p[0][0][0][0],2);
      err_rich_p[0][1][0][1] = pow(err_rich_p[0][1][0][1],2);
      err_rich_p[0][2][0][2] = pow(err_rich_p[0][2][0][2],2);
      err_rich_p[1][0][1][0] = pow(err_rich_p[1][0][1][0],2);
      err_rich_p[1][1][1][1] = pow(err_rich_p[1][1][1][1],2);
      err_rich_p[1][2][1][2] = pow(err_rich_p[1][2][1][2],2);
      err_rich_p[2][0][2][0] = pow(err_rich_p[2][0][2][0],2);
      err_rich_p[2][1][2][1] = pow(err_rich_p[2][1][2][1],2);
      err_rich_p[2][2][2][2] = pow(err_rich_p[2][2][2][2],2);

      err_rich_p[0][0][1][0] = err_rich_p[1][0][0][0];
      err_rich_p[0][0][2][0] = err_rich_p[2][0][0][0];
      err_rich_p[1][0][2][0] = err_rich_p[2][0][1][0];

      err_rich_p[0][1][1][1] = err_rich_p[1][1][0][1];
      err_rich_p[0][1][2][1] = err_rich_p[2][1][0][1];
      err_rich_p[1][1][2][1] = err_rich_p[2][1][1][1];

      err_rich_p[0][2][1][2] = err_rich_p[1][2][0][2];
      err_rich_p[0][2][2][2] = err_rich_p[2][2][0][2];
      err_rich_p[1][2][2][2] = err_rich_p[2][2][1][2];

      errRICH >> err_rich_m[0][0][0][0];
      errRICH >> err_rich_m[0][1][0][1];
      errRICH >> err_rich_m[0][2][0][2];//3
      errRICH >> err_rich_m[0][0][1][0];
      errRICH >> err_rich_m[0][0][2][0];
      errRICH >> err_rich_m[1][0][2][0];//6
      errRICH >> err_rich_m[1][0][1][0];
      errRICH >> err_rich_m[1][1][1][1];
      errRICH >> err_rich_m[1][2][1][2];//9
      errRICH >> err_rich_m[0][1][1][1];
      errRICH >> err_rich_m[0][1][2][1];
      errRICH >> err_rich_m[1][1][2][1];//12
      errRICH >> err_rich_m[2][0][2][0];
      errRICH >> err_rich_m[2][1][2][1];
      errRICH >> err_rich_m[2][2][2][2];//15
      errRICH >> err_rich_m[0][2][1][2];
      errRICH >> err_rich_m[0][2][2][2];
      errRICH >> err_rich_m[1][2][2][2];//18

      err_rich_m[0][0][0][0] = pow(err_rich_m[0][0][0][0],2);
      err_rich_m[0][1][0][1] = pow(err_rich_m[0][1][0][1],2);
      err_rich_m[0][2][0][2] = pow(err_rich_m[0][2][0][2],2);
      err_rich_m[1][0][1][0] = pow(err_rich_m[1][0][1][0],2);
      err_rich_m[1][1][1][1] = pow(err_rich_m[1][1][1][1],2);
      err_rich_m[1][2][1][2] = pow(err_rich_m[1][2][1][2],2);
      err_rich_m[2][0][2][0] = pow(err_rich_m[2][0][2][0],2);
      err_rich_m[2][1][2][1] = pow(err_rich_m[2][1][2][1],2);
      err_rich_m[2][2][2][2] = pow(err_rich_m[2][2][2][2],2);

      err_rich_m[0][0][1][0] = err_rich_m[1][0][0][0];
      err_rich_m[0][0][2][0] = err_rich_m[2][0][0][0];
      err_rich_m[1][0][2][0] = err_rich_m[2][0][1][0];

      err_rich_m[0][1][1][1] = err_rich_m[1][1][0][1];
      err_rich_m[0][1][2][1] = err_rich_m[2][1][0][1];
      err_rich_m[1][1][2][1] = err_rich_m[2][1][1][1];

      err_rich_m[0][2][1][2] = err_rich_m[1][2][0][2];
      err_rich_m[0][2][2][2] = err_rich_m[2][2][0][2];
      err_rich_m[1][2][2][2] = err_rich_m[2][2][1][2];

#ifdef DEBUG
      for(int i=0; i<9; i++)
      {
        for(int j=0; j<9; j++)
        {
          cout << err_rich_p[i][j] << " ";
        }
        cout << endl;
      }
#endif

      cout << "\n" << endl;

      for(int i=0; i<3; i++)
      {
        cov1_pi[0][i] = 0;
        cov1_pi[1][i] = 0;
        cov1_k[0][i] = 0;
        cov1_k[1][i] = 0;
        cov1_p[0][i] = 0;
        cov1_p[1][i] = 0;
        cov2[0][i] = 0;
        cov2[1][i] = 0;
      }

      for(int i=0; i<3; i++)
      {
        for(int j=0; j<3; j++)
        {
          cov1_pi[0][i] += pow(inv_rich_p[loop%2][(loop-4)/2][i][j]*pi_sigma_uni[j][j],2);
          cov1_pi[1][i] += pow(inv_rich_m[loop%2][(loop-4)/2][i][j]*pi_sigma_uni[j][j],2);
          cov1_k[0][i] += pow(inv_rich_p[loop%2][(loop-4)/2][i][j]*k_sigma_uni[j][j],2);
          cov1_k[1][i] += pow(inv_rich_m[loop%2][(loop-4)/2][i][j]*k_sigma_uni[j][j],2);
          cov1_p[0][i] += pow(inv_rich_p[loop%2][(loop-4)/2][i][j]*p_sigma_uni[j][j],2);
          cov1_p[1][i] += pow(inv_rich_m[loop%2][(loop-4)/2][i][j]*p_sigma_uni[j][j],2);

          for(int k=0; k<3; k++)
          {
            for(int l=0; l<3; l++)
            {
              for(int m=0; m<3; m++)
              {
                    cov2[0][i] += inv_rich_p[loop%2][(loop-4)/2][i][j]
                                   *inv_rich_p[loop%2][(loop-4)/2][i][l]
                                   *inv_rich_p[loop%2][(loop-4)/2][k][i]
                                   *inv_rich_p[loop%2][(loop-4)/2][m][i]
                                   *err_rich_p[j][k][l][m];
#ifdef DEBUG
                                  cout << err_rich_p[j*3+k][l*3+m] << endl;
                                  cout << j*3+k << " " << l*3+m << endl;
#endif
                    cov2[1][i] += inv_rich_m[loop%2][(loop-4)/2][i][j]
                                   *inv_rich_m[loop%2][(loop-4)/2][i][l]
                                   *inv_rich_m[loop%2][(loop-4)/2][k][i]
                                   *inv_rich_m[loop%2][(loop-4)/2][m][i]
                                   *err_rich_m[j][k][l][m];
              }
            }
          }
        }
        pi_unfolding_err_p[loop%2][(loop-4)/2][i] = cov1_pi[0][i] + cov2[0][i];
        pi_unfolding_err_m[loop%2][(loop-4)/2][i] = cov1_pi[1][i] + cov2[1][i];
        k_unfolding_err_p[loop%2][(loop-4)/2][i] = cov1_k[0][i] + cov2[0][i];
        k_unfolding_err_m[loop%2][(loop-4)/2][i] = cov1_k[1][i] + cov2[1][i];
        p_unfolding_err_p[loop%2][(loop-4)/2][i] = cov1_p[0][i] + cov2[0][i];
        p_unfolding_err_m[loop%2][(loop-4)/2][i] = cov1_p[1][i] + cov2[1][i];
      }
    }
  }

  errRICH.close();
}

void load_rich_mat_2016(string prich, string prich_err)
{

  pi_sigma_uni[0][0] = 1;
  k_sigma_uni[1][1] = 1;
  p_sigma_uni[2][2] = 1;
  pi_vect[0][0] = 1;
  k_vect[1][0] = 1;
  p_vect[2][0] = 1;

  for(int i=0; i<2; i++)
  {
    for(int j=0; j<10; j++)
    {
      rich_mat_p[i][j].ResizeTo(3,3);
      rich_mat_m[i][j].ResizeTo(3,3);
      inv_rich_p[i][j].ResizeTo(3,3);
      inv_rich_m[i][j].ResizeTo(3,3);
    }
  }

  for(int i=0; i<3; i++)
  {
    for(int j=0; j<3; j++)
    {
      err_rich_p[i][j].ResizeTo(3,3);
      err_rich_m[i][j].ResizeTo(3,3);
    }
  }

  ifstream matRICH(prich);

  for(int i=0; i<28; i++)
  {
    if(i < 8)
    {
      for(int j=0; j<20; j++)
        matRICH >> dummy;
    }
    else
    {
      matRICH >> mat_bin[0][i] >> mat_bin[1][i];
      for(int k=0; k<3; k++)
      {
        for(int j=0; j<3; j++)
        {
          matRICH >> rich_mat_m[i%2][(i-8)/2][j][k];
        }
        for(int j=0; j<3; j++)
        {
          matRICH >> rich_mat_p[i%2][(i-8)/2][j][k];
        }
      }
    }
  }

  matRICH.close();

  for(int i=0; i<10; i++)
  {

    cout << rich_mat_m[0][i][0][0] << " " << rich_mat_m[0][i][0][1] << " " << rich_mat_m[0][i][0][2] << endl;
    cout << rich_mat_m[0][i][1][0] << " " << rich_mat_m[0][i][1][1] << " " << rich_mat_m[0][i][1][2] << endl;
    cout << rich_mat_m[0][i][2][0] << " " << rich_mat_m[0][i][2][1] << " " << rich_mat_m[0][i][2][2] << endl;

    cout << endl;

    cout << rich_mat_p[0][i][0][0] << " " << rich_mat_p[0][i][0][1] << " " << rich_mat_p[0][i][0][2] << endl;
    cout << rich_mat_p[0][i][1][0] << " " << rich_mat_p[0][i][1][1] << " " << rich_mat_p[0][i][1][2] << endl;
    cout << rich_mat_p[0][i][2][0] << " " << rich_mat_p[0][i][2][1] << " " << rich_mat_p[0][i][2][2] << endl;

    cout << endl;

    cout << rich_mat_m[1][i][0][0] << " " << rich_mat_m[1][i][0][1] << " " << rich_mat_m[1][i][0][2] << endl;
    cout << rich_mat_m[1][i][1][0] << " " << rich_mat_m[1][i][1][1] << " " << rich_mat_m[1][i][1][2] << endl;
    cout << rich_mat_m[1][i][2][0] << " " << rich_mat_m[1][i][2][1] << " " << rich_mat_m[1][i][2][2] << endl;

    cout << endl;

    cout << rich_mat_p[1][i][0][0] << " " << rich_mat_p[1][i][0][1] << " " << rich_mat_p[1][i][0][2] << endl;
    cout << rich_mat_p[1][i][1][0] << " " << rich_mat_p[1][i][1][1] << " " << rich_mat_p[1][i][1][2] << endl;
    cout << rich_mat_p[1][i][2][0] << " " << rich_mat_p[1][i][2][1] << " " << rich_mat_p[1][i][2][2] << endl;

    cout << endl;

    inv_rich_p[0][i] = rich_mat_p[0][i].InvertFast();
    inv_rich_p[1][i] = rich_mat_p[1][i].InvertFast();
    inv_rich_m[0][i] = rich_mat_m[0][i].InvertFast();
    inv_rich_m[1][i] = rich_mat_m[1][i].InvertFast();

    // cout << inv_rich_p[0][i][0][0] << " " << inv_rich_p[0][i][0][1] << " " << inv_rich_p[0][i][0][2] << endl;
    // cout << inv_rich_p[0][i][1][0] << " " << inv_rich_p[0][i][1][1] << " " << inv_rich_p[0][i][1][2] << endl;
    // cout << inv_rich_p[0][i][2][0] << " " << inv_rich_p[0][i][2][1] << " " << inv_rich_p[0][i][2][2] << endl;
    //
    // cout << endl;
    //
    //
    // cout << inv_rich_p[1][i][0][0] << " " << inv_rich_p[1][i][0][1] << " " << inv_rich_p[1][i][0][2] << endl;
    // cout << inv_rich_p[1][i][1][0] << " " << inv_rich_p[1][i][1][1] << " " << inv_rich_p[1][i][1][2] << endl;
    // cout << inv_rich_p[1][i][2][0] << " " << inv_rich_p[1][i][2][1] << " " << inv_rich_p[1][i][2][2] << endl;
    //
    // cout << endl;
    //
    //
    // cout << inv_rich_m[0][i][0][0] << " " << inv_rich_m[0][i][0][1] << " " << inv_rich_m[0][i][0][2] << endl;
    // cout << inv_rich_m[0][i][1][0] << " " << inv_rich_m[0][i][1][1] << " " << inv_rich_m[0][i][1][2] << endl;
    // cout << inv_rich_m[0][i][2][0] << " " << inv_rich_m[0][i][2][1] << " " << inv_rich_m[0][i][2][2] << endl;
    //
    // cout << endl;
    //
    // cout << inv_rich_m[1][i][0][0] << " " << inv_rich_m[1][i][0][1] << " " << inv_rich_m[1][i][0][2] << endl;
    // cout << inv_rich_m[1][i][1][0] << " " << inv_rich_m[1][i][1][1] << " " << inv_rich_m[1][i][1][2] << endl;
    // cout << inv_rich_m[1][i][2][0] << " " << inv_rich_m[1][i][2][1] << " " << inv_rich_m[1][i][2][2] << endl;
    //
    // cout << endl;
  }

  // Errors YODO

  ifstream errRICH(prich_err);


  for(int loop=0; loop<24; loop++)
  {
    if(loop < 4)
    {
      for(int j=0; j<38; j++)
        errRICH >> dummy;
    }
    else
    {
      errRICH >> err_bin[0][loop-4] >> err_bin[1][loop-4];

      errRICH >> err_rich_p[0][0][0][0];
      errRICH >> err_rich_p[0][1][0][1];
      errRICH >> err_rich_p[0][2][0][2];//3
      errRICH >> err_rich_p[0][0][1][0];
      errRICH >> err_rich_p[0][0][2][0];
      errRICH >> err_rich_p[1][0][2][0];//6
      errRICH >> err_rich_p[1][0][1][0];
      errRICH >> err_rich_p[1][1][1][1];
      errRICH >> err_rich_p[1][2][1][2];//9
      errRICH >> err_rich_p[0][1][1][1];
      errRICH >> err_rich_p[0][1][2][1];
      errRICH >> err_rich_p[1][1][2][1];//12
      errRICH >> err_rich_p[2][0][2][0];
      errRICH >> err_rich_p[2][1][2][1];
      errRICH >> err_rich_p[2][2][2][2];//15
      errRICH >> err_rich_p[0][2][1][2];
      errRICH >> err_rich_p[0][2][2][2];
      errRICH >> err_rich_p[1][2][2][2];//18

      err_rich_p[0][0][0][0] = pow(err_rich_p[0][0][0][0],2);
      err_rich_p[0][1][0][1] = pow(err_rich_p[0][1][0][1],2);
      err_rich_p[0][2][0][2] = pow(err_rich_p[0][2][0][2],2);
      err_rich_p[1][0][1][0] = pow(err_rich_p[1][0][1][0],2);
      err_rich_p[1][1][1][1] = pow(err_rich_p[1][1][1][1],2);
      err_rich_p[1][2][1][2] = pow(err_rich_p[1][2][1][2],2);
      err_rich_p[2][0][2][0] = pow(err_rich_p[2][0][2][0],2);
      err_rich_p[2][1][2][1] = pow(err_rich_p[2][1][2][1],2);
      err_rich_p[2][2][2][2] = pow(err_rich_p[2][2][2][2],2);

      err_rich_p[0][0][1][0] = err_rich_p[1][0][0][0];
      err_rich_p[0][0][2][0] = err_rich_p[2][0][0][0];
      err_rich_p[1][0][2][0] = err_rich_p[2][0][1][0];

      err_rich_p[0][1][1][1] = err_rich_p[1][1][0][1];
      err_rich_p[0][1][2][1] = err_rich_p[2][1][0][1];
      err_rich_p[1][1][2][1] = err_rich_p[2][1][1][1];

      err_rich_p[0][2][1][2] = err_rich_p[1][2][0][2];
      err_rich_p[0][2][2][2] = err_rich_p[2][2][0][2];
      err_rich_p[1][2][2][2] = err_rich_p[2][2][1][2];

      errRICH >> err_rich_m[0][0][0][0];
      errRICH >> err_rich_m[0][1][0][1];
      errRICH >> err_rich_m[0][2][0][2];//3
      errRICH >> err_rich_m[0][0][1][0];
      errRICH >> err_rich_m[0][0][2][0];
      errRICH >> err_rich_m[1][0][2][0];//6
      errRICH >> err_rich_m[1][0][1][0];
      errRICH >> err_rich_m[1][1][1][1];
      errRICH >> err_rich_m[1][2][1][2];//9
      errRICH >> err_rich_m[0][1][1][1];
      errRICH >> err_rich_m[0][1][2][1];
      errRICH >> err_rich_m[1][1][2][1];//12
      errRICH >> err_rich_m[2][0][2][0];
      errRICH >> err_rich_m[2][1][2][1];
      errRICH >> err_rich_m[2][2][2][2];//15
      errRICH >> err_rich_m[0][2][1][2];
      errRICH >> err_rich_m[0][2][2][2];
      errRICH >> err_rich_m[1][2][2][2];//18

      err_rich_m[0][0][0][0] = pow(err_rich_m[0][0][0][0],2);
      err_rich_m[0][1][0][1] = pow(err_rich_m[0][1][0][1],2);
      err_rich_m[0][2][0][2] = pow(err_rich_m[0][2][0][2],2);
      err_rich_m[1][0][1][0] = pow(err_rich_m[1][0][1][0],2);
      err_rich_m[1][1][1][1] = pow(err_rich_m[1][1][1][1],2);
      err_rich_m[1][2][1][2] = pow(err_rich_m[1][2][1][2],2);
      err_rich_m[2][0][2][0] = pow(err_rich_m[2][0][2][0],2);
      err_rich_m[2][1][2][1] = pow(err_rich_m[2][1][2][1],2);
      err_rich_m[2][2][2][2] = pow(err_rich_m[2][2][2][2],2);

      err_rich_m[0][0][1][0] = err_rich_m[1][0][0][0];
      err_rich_m[0][0][2][0] = err_rich_m[2][0][0][0];
      err_rich_m[1][0][2][0] = err_rich_m[2][0][1][0];

      err_rich_m[0][1][1][1] = err_rich_m[1][1][0][1];
      err_rich_m[0][1][2][1] = err_rich_m[2][1][0][1];
      err_rich_m[1][1][2][1] = err_rich_m[2][1][1][1];

      err_rich_m[0][2][1][2] = err_rich_m[1][2][0][2];
      err_rich_m[0][2][2][2] = err_rich_m[2][2][0][2];
      err_rich_m[1][2][2][2] = err_rich_m[2][2][1][2];

#ifdef DEBUG
      for(int i=0; i<9; i++)
      {
        for(int j=0; j<9; j++)
        {
          cout << err_rich_p[i][j] << " ";
        }
        cout << endl;
      }
#endif

      cout << "\n" << endl;

      for(int i=0; i<3; i++)
      {
        cov1_pi[0][i] = 0;
        cov1_pi[1][i] = 0;
        cov1_k[0][i] = 0;
        cov1_k[1][i] = 0;
        cov1_p[0][i] = 0;
        cov1_p[1][i] = 0;
        cov2[0][i] = 0;
        cov2[1][i] = 0;
      }

      for(int i=0; i<3; i++)
      {
        for(int j=0; j<3; j++)
        {
          cov1_pi[0][i] += pow(inv_rich_p[loop%2][(loop-4)/2][i][j]*pi_sigma_uni[j][j],2);
          cov1_pi[1][i] += pow(inv_rich_m[loop%2][(loop-4)/2][i][j]*pi_sigma_uni[j][j],2);
          cov1_k[0][i] += pow(inv_rich_p[loop%2][(loop-4)/2][i][j]*k_sigma_uni[j][j],2);
          cov1_k[1][i] += pow(inv_rich_m[loop%2][(loop-4)/2][i][j]*k_sigma_uni[j][j],2);
          cov1_p[0][i] += pow(inv_rich_p[loop%2][(loop-4)/2][i][j]*p_sigma_uni[j][j],2);
          cov1_p[1][i] += pow(inv_rich_m[loop%2][(loop-4)/2][i][j]*p_sigma_uni[j][j],2);

          for(int k=0; k<3; k++)
          {
            for(int l=0; l<3; l++)
            {
              for(int m=0; m<3; m++)
              {
                    cov2[0][i] += inv_rich_p[loop%2][(loop-4)/2][i][j]
                                   *inv_rich_p[loop%2][(loop-4)/2][i][l]
                                   *inv_rich_p[loop%2][(loop-4)/2][k][i]
                                   *inv_rich_p[loop%2][(loop-4)/2][m][i]
                                   *err_rich_p[j][k][l][m];
#ifdef DEBUG
                                  cout << err_rich_p[j*3+k][l*3+m] << endl;
                                  cout << j*3+k << " " << l*3+m << endl;
#endif
                    cov2[1][i] += inv_rich_m[loop%2][(loop-4)/2][i][j]
                                   *inv_rich_m[loop%2][(loop-4)/2][i][l]
                                   *inv_rich_m[loop%2][(loop-4)/2][k][i]
                                   *inv_rich_m[loop%2][(loop-4)/2][m][i]
                                   *err_rich_m[j][k][l][m];
              }
            }
          }
        }
        pi_unfolding_err_p[loop%2][(loop-4)/2][i] = cov1_pi[0][i] + cov2[0][i];
        pi_unfolding_err_m[loop%2][(loop-4)/2][i] = cov1_pi[1][i] + cov2[1][i];
        k_unfolding_err_p[loop%2][(loop-4)/2][i] = cov1_k[0][i] + cov2[0][i];
        k_unfolding_err_m[loop%2][(loop-4)/2][i] = cov1_k[1][i] + cov2[1][i];
        p_unfolding_err_p[loop%2][(loop-4)/2][i] = cov1_p[0][i] + cov2[0][i];
        p_unfolding_err_m[loop%2][(loop-4)/2][i] = cov1_p[1][i] + cov2[1][i];
      }
    }
  }

  errRICH.close();
}

void load_rich_mat_dummy(string prich, string prich_err)
{

  pi_sigma_uni[0][0] = 1;
  k_sigma_uni[1][1] = 1;
  p_sigma_uni[2][2] = 1;
  pi_vect[0][0] = 1;
  k_vect[1][0] = 1;
  p_vect[2][0] = 1;

  for(int i=0; i<2; i++)
  {
    for(int j=0; j<10; j++)
    {
      rich_mat_p[i][j].ResizeTo(3,3);
      rich_mat_m[i][j].ResizeTo(3,3);
      inv_rich_p[i][j].ResizeTo(3,3);
      inv_rich_m[i][j].ResizeTo(3,3);
    }
  }

  for(int i=0; i<3; i++)
  {
    for(int j=0; j<3; j++)
    {
      err_rich_p[i][j].ResizeTo(3,3);
      err_rich_m[i][j].ResizeTo(3,3);
    }
  }

  for(int i=0; i<20; i++)
  {
        rich_mat_m[i%2][(i)/2][0][0]=1;
        rich_mat_m[i%2][(i)/2][1][1]=1;
        rich_mat_m[i%2][(i)/2][2][2]=1;
        rich_mat_p[i%2][(i)/2][0][0]=1;
        rich_mat_p[i%2][(i)/2][1][1]=1;
        rich_mat_p[i%2][(i)/2][2][2]=1;
  }

  for(int i=0; i<10; i++)
  {
    cout << rich_mat_p[0][i][0][0] << " " << rich_mat_p[0][i][0][1] << " " << rich_mat_p[0][i][0][2] << endl;
    cout << rich_mat_p[0][i][1][0] << " " << rich_mat_p[0][i][1][1] << " " << rich_mat_p[0][i][1][2] << endl;
    cout << rich_mat_p[0][i][2][0] << " " << rich_mat_p[0][i][2][1] << " " << rich_mat_p[0][i][2][2] << endl;

    cout << endl;

    cout << rich_mat_p[1][i][0][0] << " " << rich_mat_p[1][i][0][1] << " " << rich_mat_p[1][i][0][2] << endl;
    cout << rich_mat_p[1][i][1][0] << " " << rich_mat_p[1][i][1][1] << " " << rich_mat_p[1][i][1][2] << endl;
    cout << rich_mat_p[1][i][2][0] << " " << rich_mat_p[1][i][2][1] << " " << rich_mat_p[1][i][2][2] << endl;

    cout << endl;

    cout << rich_mat_m[0][i][0][0] << " " << rich_mat_m[0][i][0][1] << " " << rich_mat_m[0][i][0][2] << endl;
    cout << rich_mat_m[0][i][1][0] << " " << rich_mat_m[0][i][1][1] << " " << rich_mat_m[0][i][1][2] << endl;
    cout << rich_mat_m[0][i][2][0] << " " << rich_mat_m[0][i][2][1] << " " << rich_mat_m[0][i][2][2] << endl;

    cout << endl;

    cout << rich_mat_m[1][i][0][0] << " " << rich_mat_m[1][i][0][1] << " " << rich_mat_m[1][i][0][2] << endl;
    cout << rich_mat_m[1][i][1][0] << " " << rich_mat_m[1][i][1][1] << " " << rich_mat_m[1][i][1][2] << endl;
    cout << rich_mat_m[1][i][2][0] << " " << rich_mat_m[1][i][2][1] << " " << rich_mat_m[1][i][2][2] << endl;

    cout << endl;

    inv_rich_p[0][i] = rich_mat_p[0][i].InvertFast();
    inv_rich_p[1][i] = rich_mat_p[1][i].InvertFast();
    inv_rich_m[0][i] = rich_mat_m[0][i].InvertFast();
    inv_rich_m[1][i] = rich_mat_m[1][i].InvertFast();

    cout << inv_rich_p[0][i][0][0] << " " << inv_rich_p[0][i][0][1] << " " << inv_rich_p[0][i][0][2] << endl;
    cout << inv_rich_p[0][i][1][0] << " " << inv_rich_p[0][i][1][1] << " " << inv_rich_p[0][i][1][2] << endl;
    cout << inv_rich_p[0][i][2][0] << " " << inv_rich_p[0][i][2][1] << " " << inv_rich_p[0][i][2][2] << endl;

    cout << endl;


    cout << inv_rich_p[1][i][0][0] << " " << inv_rich_p[1][i][0][1] << " " << inv_rich_p[1][i][0][2] << endl;
    cout << inv_rich_p[1][i][1][0] << " " << inv_rich_p[1][i][1][1] << " " << inv_rich_p[1][i][1][2] << endl;
    cout << inv_rich_p[1][i][2][0] << " " << inv_rich_p[1][i][2][1] << " " << inv_rich_p[1][i][2][2] << endl;

    cout << endl;


    cout << inv_rich_m[0][i][0][0] << " " << inv_rich_m[0][i][0][1] << " " << inv_rich_m[0][i][0][2] << endl;
    cout << inv_rich_m[0][i][1][0] << " " << inv_rich_m[0][i][1][1] << " " << inv_rich_m[0][i][1][2] << endl;
    cout << inv_rich_m[0][i][2][0] << " " << inv_rich_m[0][i][2][1] << " " << inv_rich_m[0][i][2][2] << endl;

    cout << endl;

    cout << inv_rich_m[1][i][0][0] << " " << inv_rich_m[1][i][0][1] << " " << inv_rich_m[1][i][0][2] << endl;
    cout << inv_rich_m[1][i][1][0] << " " << inv_rich_m[1][i][1][1] << " " << inv_rich_m[1][i][1][2] << endl;
    cout << inv_rich_m[1][i][2][0] << " " << inv_rich_m[1][i][2][1] << " " << inv_rich_m[1][i][2][2] << endl;

    cout << endl;
  }

  // Errors YODO

  ifstream errRICH(prich_err);


  for(int loop=0; loop<24; loop++)
  {
    if(loop < 4)
    {
      for(int j=0; j<38; j++)
        errRICH >> dummy;
    }
    else
    {
      errRICH >> err_bin[0][loop-4] >> err_bin[1][loop-4];

      errRICH >> err_rich_p[0][0][0][0];
      errRICH >> err_rich_p[0][1][0][1];
      errRICH >> err_rich_p[0][2][0][2];//3
      errRICH >> err_rich_p[0][0][1][0];
      errRICH >> err_rich_p[0][0][2][0];
      errRICH >> err_rich_p[1][0][2][0];//6
      errRICH >> err_rich_p[1][0][1][0];
      errRICH >> err_rich_p[1][1][1][1];
      errRICH >> err_rich_p[1][2][1][2];//9
      errRICH >> err_rich_p[0][1][1][1];
      errRICH >> err_rich_p[0][1][2][1];
      errRICH >> err_rich_p[1][1][2][1];//12
      errRICH >> err_rich_p[2][0][2][0];
      errRICH >> err_rich_p[2][1][2][1];
      errRICH >> err_rich_p[2][2][2][2];//15
      errRICH >> err_rich_p[0][2][1][2];
      errRICH >> err_rich_p[0][2][2][2];
      errRICH >> err_rich_p[1][2][2][2];//18

      err_rich_p[0][0][0][0] = pow(err_rich_p[0][0][0][0],2);
      err_rich_p[0][1][0][1] = pow(err_rich_p[0][1][0][1],2);
      err_rich_p[0][2][0][2] = pow(err_rich_p[0][2][0][2],2);
      err_rich_p[1][0][1][0] = pow(err_rich_p[1][0][1][0],2);
      err_rich_p[1][1][1][1] = pow(err_rich_p[1][1][1][1],2);
      err_rich_p[1][2][1][2] = pow(err_rich_p[1][2][1][2],2);
      err_rich_p[2][0][2][0] = pow(err_rich_p[2][0][2][0],2);
      err_rich_p[2][1][2][1] = pow(err_rich_p[2][1][2][1],2);
      err_rich_p[2][2][2][2] = pow(err_rich_p[2][2][2][2],2);

      err_rich_p[0][0][1][0] = err_rich_p[1][0][0][0];
      err_rich_p[0][0][2][0] = err_rich_p[2][0][0][0];
      err_rich_p[1][0][2][0] = err_rich_p[2][0][1][0];

      err_rich_p[0][1][1][1] = err_rich_p[1][1][0][1];
      err_rich_p[0][1][2][1] = err_rich_p[2][1][0][1];
      err_rich_p[1][1][2][1] = err_rich_p[2][1][1][1];

      err_rich_p[0][2][1][2] = err_rich_p[1][2][0][2];
      err_rich_p[0][2][2][2] = err_rich_p[2][2][0][2];
      err_rich_p[1][2][2][2] = err_rich_p[2][2][1][2];

      errRICH >> err_rich_m[0][0][0][0];
      errRICH >> err_rich_m[0][1][0][1];
      errRICH >> err_rich_m[0][2][0][2];//3
      errRICH >> err_rich_m[0][0][1][0];
      errRICH >> err_rich_m[0][0][2][0];
      errRICH >> err_rich_m[1][0][2][0];//6
      errRICH >> err_rich_m[1][0][1][0];
      errRICH >> err_rich_m[1][1][1][1];
      errRICH >> err_rich_m[1][2][1][2];//9
      errRICH >> err_rich_m[0][1][1][1];
      errRICH >> err_rich_m[0][1][2][1];
      errRICH >> err_rich_m[1][1][2][1];//12
      errRICH >> err_rich_m[2][0][2][0];
      errRICH >> err_rich_m[2][1][2][1];
      errRICH >> err_rich_m[2][2][2][2];//15
      errRICH >> err_rich_m[0][2][1][2];
      errRICH >> err_rich_m[0][2][2][2];
      errRICH >> err_rich_m[1][2][2][2];//18

      err_rich_m[0][0][0][0] = pow(err_rich_m[0][0][0][0],2);
      err_rich_m[0][1][0][1] = pow(err_rich_m[0][1][0][1],2);
      err_rich_m[0][2][0][2] = pow(err_rich_m[0][2][0][2],2);
      err_rich_m[1][0][1][0] = pow(err_rich_m[1][0][1][0],2);
      err_rich_m[1][1][1][1] = pow(err_rich_m[1][1][1][1],2);
      err_rich_m[1][2][1][2] = pow(err_rich_m[1][2][1][2],2);
      err_rich_m[2][0][2][0] = pow(err_rich_m[2][0][2][0],2);
      err_rich_m[2][1][2][1] = pow(err_rich_m[2][1][2][1],2);
      err_rich_m[2][2][2][2] = pow(err_rich_m[2][2][2][2],2);

      err_rich_m[0][0][1][0] = err_rich_m[1][0][0][0];
      err_rich_m[0][0][2][0] = err_rich_m[2][0][0][0];
      err_rich_m[1][0][2][0] = err_rich_m[2][0][1][0];

      err_rich_m[0][1][1][1] = err_rich_m[1][1][0][1];
      err_rich_m[0][1][2][1] = err_rich_m[2][1][0][1];
      err_rich_m[1][1][2][1] = err_rich_m[2][1][1][1];

      err_rich_m[0][2][1][2] = err_rich_m[1][2][0][2];
      err_rich_m[0][2][2][2] = err_rich_m[2][2][0][2];
      err_rich_m[1][2][2][2] = err_rich_m[2][2][1][2];

#ifdef DEBUG
      for(int i=0; i<9; i++)
      {
        for(int j=0; j<9; j++)
        {
          cout << err_rich_p[i][j] << " ";
        }
        cout << endl;
      }
#endif

      cout << "\n" << endl;

      for(int i=0; i<3; i++)
      {
        cov1_pi[0][i] = 0;
        cov1_pi[1][i] = 0;
        cov1_k[0][i] = 0;
        cov1_k[1][i] = 0;
        cov1_p[0][i] = 0;
        cov1_p[1][i] = 0;
        cov2[0][i] = 0;
        cov2[1][i] = 0;
      }

      for(int i=0; i<3; i++)
      {
        for(int j=0; j<3; j++)
        {
          cov1_pi[0][i] += pow(inv_rich_p[loop%2][(loop-4)/2][i][j]*pi_sigma_uni[j][j],2);
          cov1_pi[1][i] += pow(inv_rich_m[loop%2][(loop-4)/2][i][j]*pi_sigma_uni[j][j],2);
          cov1_k[0][i] += pow(inv_rich_p[loop%2][(loop-4)/2][i][j]*k_sigma_uni[j][j],2);
          cov1_k[1][i] += pow(inv_rich_m[loop%2][(loop-4)/2][i][j]*k_sigma_uni[j][j],2);
          cov1_p[0][i] += pow(inv_rich_p[loop%2][(loop-4)/2][i][j]*p_sigma_uni[j][j],2);
          cov1_p[1][i] += pow(inv_rich_m[loop%2][(loop-4)/2][i][j]*p_sigma_uni[j][j],2);

          for(int k=0; k<3; k++)
          {
            for(int l=0; l<3; l++)
            {
              for(int m=0; m<3; m++)
              {
                    cov2[0][i] += inv_rich_p[loop%2][(loop-4)/2][i][j]
                                   *inv_rich_p[loop%2][(loop-4)/2][i][l]
                                   *inv_rich_p[loop%2][(loop-4)/2][k][i]
                                   *inv_rich_p[loop%2][(loop-4)/2][m][i]
                                   *err_rich_p[j][k][l][m];
#ifdef DEBUG
                                  cout << err_rich_p[j*3+k][l*3+m] << endl;
                                  cout << j*3+k << " " << l*3+m << endl;
#endif
                    cov2[1][i] += inv_rich_m[loop%2][(loop-4)/2][i][j]
                                   *inv_rich_m[loop%2][(loop-4)/2][i][l]
                                   *inv_rich_m[loop%2][(loop-4)/2][k][i]
                                   *inv_rich_m[loop%2][(loop-4)/2][m][i]
                                   *err_rich_m[j][k][l][m];
              }
            }
          }
        }
        pi_unfolding_err_p[loop%2][(loop-4)/2][i] = cov1_pi[0][i] + cov2[0][i];
        pi_unfolding_err_m[loop%2][(loop-4)/2][i] = cov1_pi[1][i] + cov2[1][i];
        k_unfolding_err_p[loop%2][(loop-4)/2][i] = cov1_k[0][i] + cov2[0][i];
        k_unfolding_err_m[loop%2][(loop-4)/2][i] = cov1_k[1][i] + cov2[1][i];
        p_unfolding_err_p[loop%2][(loop-4)/2][i] = cov1_p[0][i] + cov2[0][i];
        p_unfolding_err_m[loop%2][(loop-4)/2][i] = cov1_p[1][i] + cov2[1][i];
      }
    }
  }

  errRICH.close();
}

//******************************************************************************
// LOG BINNING
//******************************************************************************

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

void BinLogY(TH1*h)
{
   TAxis *axis = h->GetYaxis();
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

//******************************************************************************
// PLOT HANDLING
//******************************************************************************

void create_kin_plots()
{
  fKinematics[0] = new TH1F("Q^{2}", "Q^{2}", 100, 0, 2);
  fKinematics[1] = new TH1F("x_{Bj}", "x_{Bj}", 100, -3, 0);
  fKinematics[2] = new TH1F("y", "y", 100, 0, 1);
  fKinematics[3] = new TH1F("z", "z", 100, 0, 1);
  fKinematics[4] = new TH1F("W", "W", 100, 2, 18);
  fKinematics[5] = new TH1F("#nu", "#nu", 100, 0, 160);
  fKinematics2D[0] = new TH2F("DIS kin space (x,y)", "DIS kin space (x,y)", 300, -3, 0, 300, 0, 0.9);
  fKinematics2D[1] = new TH2F("DIS kin space (x,Q2)", "DIS kin space (x,Q2)", 300, -3, 0, 300, -0.097, 2);
  fKinematicsRICH = new TH2F("RICH spectrum", "RICH spectrum", 500, 0, 60, 500, 20, 60);
  fTarget2D = new TH2F("Target xy", "Target xy", 100, -3, 3, 100, -3, 3);
  fHO03 = new TH2F("HO03", "HO03", 100, -120, 120, 100, -60, 60);
  fHO04 = new TH2F("HO04", "HO04", 100, -250, 250, 100, -100, 100);
  fRICHLH = new TH2F("RICH LH", "RICH LH", 100, -2, 1.5, 100, -2, 1.5);
  fInTarget[0] = new TH2F("Target XZ", "Target XZ", 100, -350, -50, 100, -3, 3);
  fInTarget[1] = new TH2F("Target YZ", "Target YZ", 100, -350, -50, 100, -3, 3);
  fAllTarget[0] = new TH2F("Target XZ2", "Target XZ2", 100, -350, -50, 100, -3, 3);
  fAllTarget[1] = new TH2F("Target YZ2", "Target YZ2", 100, -350, -50, 100, -3, 3);
  fxMT = new TH1F("xMT","XMT", 100, -3, 0);
  fxLT = new TH1F("xLT","xLT", 100, -3, 0);
  fxOT = new TH1F("xOT","xOT", 100, -3, 0);
  fxLAST = new TH1F("xLAST","xLAST", 100, -3, 0);
  fQ2MT = new TH1F("Q2MT","Q2MT", 100, -0.09, 2);
  fQ2LT = new TH1F("Q2LT","Q2LT", 100, -0.09, 2);
  fQ2OT = new TH1F("Q2OT","Q2OT", 100, -0.09, 2);
  fQ2LAST = new TH1F("Q2LAST","Q2LAST", 100, -0.09, 2);
  for(int i=0; i<10; i++)
  {
    fInTarget[i+2] = new TH2F(Form("Target XY z%d",i), Form("Target YZ z%d",i), 500, -3, 3, 500, -3, 3);
    fAllTarget[i+2] = new TH2F(Form("Target XY2 z%d",i), Form("Target XZ2 z%d",i), 200, -3, 3, 200, -3, 3);
  }
  fZvtx[0] = new TH1F("Zvtx1","Zvtx1", 200, -500, 500);
  fZvtx[1] = new TH1F("Zvtx2","Zvtx2", 200, -500, 500);
  fQ2k[0] = new TH1F("Q^{2}1","Q^{2}1", 200, -1*log10(8), 2);
  fQ2k[1] = new TH1F("Q^{2}2","Q^{2}2", 200, -1*log10(8), 2);
  fYk[0] = new TH1F("y1","y1", 200, 0, 1);
  fYk[1] = new TH1F("y2","y2", 200, 0, 1);
  fThRich[0] = new TH1F("#theta_{RICH}1","#theta_{RICH}1", 200, 0, 0.8);
  fThRich[1] = new TH1F("#theta_{RICH}2","#theta_{RICH}2", 200, 0, 0.8);
  fZk[0] = new TH1F("z1","z1", 200, 0, 1);
  fZk[1] = new TH1F("z2","z2", 200, 0, 1);
  fPk[0] = new TH1F("p_{h}1","p_{h}1", 200, 0, 80);
  fPk[1] = new TH1F("p_{h}2","p_{h}2", 200, 0, 80);
  BinLogX(fKinematics[0]);
  BinLogX(fKinematics[1]);
  BinLogX(fKinematics2D[0]);
  fKinematics2D[0]->SetStats(0);
  BinLogX(fKinematics2D[1]);
  BinLogY(fKinematics2D[1]);
  fKinematics2D[1]->SetStats(0);
  BinLogX(fQ2k[0]);
  BinLogX(fQ2k[1]);
  BinLogX(fxMT);
  BinLogX(fxLT);
  BinLogX(fxOT);
  BinLogX(fxLAST);
  BinLogX(fQ2MT);
  BinLogX(fQ2LT);
  BinLogX(fQ2OT);
  BinLogX(fQ2LAST);
}

void save_kin_plots()
{
  // gStyle->SetPalette(kRainBow);
  gStyle->SetPalette(kRainBow);
  c1.Divide(1,1);
  c2.Divide(1,1);
  c3.Divide(1,1);
  c4.Divide(1,1);
  c5.Divide(1,1);
  c6.Divide(1,1);
  c7.Divide(1,1);
  c71.Divide(1,1);
  c8.Divide(1,1);
  c9.Divide(1,1);
  c10.Divide(1,1);
  c11.Divide(1,1);
  c12.Divide(1,1);
  c13.Divide(1,2);
  c14.Divide(5,2);
  c15.Divide(3,2);
  c16.Divide(1,1);
  c17.Divide(1,1);
  c1.cd(1);
  fKinematics[0]->SetStats(0);
  fKinematics[0]->Draw();
  gPad->SetLogx();
  c1.Update();
  c2.cd(1);
  fKinematics[1]->SetStats(0);
  fKinematics[1]->Draw();
  gPad->SetLogx();
  c2.Update();
  c3.cd(1);
  fKinematics[2]->SetStats(0);
  fKinematics[2]->Draw();
  c3.Update();
  c4.cd(1);
  fKinematics[3]->SetStats(0);
  fKinematics[3]->Draw();
  c4.Update();
  c5.cd(1);
  fKinematics[4]->SetStats(0);
  fKinematics[4]->Draw();
  c5.Update();
  c6.cd(1);
  fKinematics[5]->SetStats(0);
  fKinematics[5]->Draw();
  c6.Update();
  c7.cd(1);
  fKinematics2D[0]->SetStats(0);
  fKinematics2D[0]->Draw("COLZ");
  gPad->SetLogx();
  c7.Update();
  c71.cd(1);
  fKinematics2D[1]->SetStats(0);
  fKinematics2D[1]->Draw("COLZ");
  gPad->SetLogx();
  gPad->SetLogy();
  gPad->SetLogz();
  c71.Update();
  c8.cd(1);
  fTarget2D->SetStats(0);
  fTarget2D->Draw("COLZ");
  c8.Update();
  c9.cd(1);
  fRICHLH->SetStats(0);
  fRICHLH->Draw("COLZ");
  c9.Update();
  c10.cd(1);
  fHO03->SetStats(0);
  fHO03->Draw("COLZ");
  c10.Update();
  c11.cd(1);
  fHO04->SetStats(0);
  fHO04->Draw("COLZ");
  c11.Update();
  gStyle->SetPalette(kColorPrintableOnGrey);
  c12.cd(1);
  fKinematicsRICH->SetStats(0);
  fKinematicsRICH->Draw("COLZ");
  gPad->SetLogz();
  TF1 *fe = new TF1("fe","acos((1/1.00142)*sqrt((pow(0.000511,2)/pow(x,2))+1))*1000",0,60);
  fe->SetLineColor(10);
  fe->SetLineStyle(9);
  fe->Draw("SAME");
  TF1 *fmu = new TF1("fmu","acos((1/1.00142)*sqrt((pow(0.10566,2)/pow(x,2))+1))*1000",0,60);
  fmu->SetLineColor(10);
  fmu->SetLineStyle(9);
  fmu->Draw("SAME");
  TF1 *fpi = new TF1("fpi","acos((1/1.00142)*sqrt((pow(0.13957,2)/pow(x,2))+1))*1000",0,60);
  fpi->SetLineColor(10);
  fpi->SetLineStyle(9);
  fpi->Draw("SAME");
  TF1 *fk = new TF1("fk","acos((1/1.00142)*sqrt((pow(0.4937,2)/pow(x,2))+1))*1000",0,60);
  fk->SetLineColor(10);
  fk->SetLineStyle(9);
  fk->Draw("SAME");
  TF1 *fp = new TF1("fp","acos((1/1.00142)*sqrt((pow(0.93827,2)/pow(x,2))+1))*1000",0,60);
  fp->SetLineColor(10);
  fp->SetLineStyle(9);
  fp->Draw("SAME");
  c12.Update();
  // c13.cd(1);
  // fAllTarget[0]->SetStats(0);
  // fInTarget[0]->SetStats(0);
  // fInTarget[0]->SetMarkerColor(kRed);
  // fAllTarget[0]->Draw("");
  // fInTarget[0]->Draw("SAME");
  // c13.Update();
  // c13.cd(2);
  // fAllTarget[1]->SetStats(0);
  // fInTarget[1]->SetStats(0);
  // fInTarget[1]->SetMarkerColor(kRed);
  // fAllTarget[1]->Draw("");
  // fInTarget[1]->Draw("SAME");
  // c13.Update();
  // for(int i=0; i<10; i++)
  // {
  //   c14.cd(i+1);
  //   fAllTarget[i+2]->SetStats(0);
  //   fInTarget[i+2]->SetStats(0);
  //   fInTarget[i+2]->SetMarkerColor(kRed);
  //   fAllTarget[i+2]->Draw("");
  //   fInTarget[i+2]->Draw("SAME");
  //   c14.Update();
  // }
  gStyle->SetPalette(kRainBow);
  c15.cd(1);
  fZvtx[0]->SetStats(0);
  fZvtx[0]->SetTitle("z_{vtx}");
  fZvtx[0]->GetXaxis()->SetTitle("z_{vtx}");
  fZvtx[0]->GetYaxis()->SetTitle("Entries");
  fZvtx[1]->SetStats(0);
  fZvtx[1]->SetFillColor(kYellow);
  fZvtx[0]->Draw();
  fZvtx[1]->Draw("SAME");
  c15.Update();
  c15.cd(2);
  fQ2k[0]->SetStats(0);
  fQ2k[0]->SetTitle("Q^{2}");
  fQ2k[0]->GetXaxis()->SetTitle("Q^{2}");
  fQ2k[0]->GetYaxis()->SetTitle("Entries");
  fQ2k[1]->SetStats(0);
  fQ2k[1]->SetFillColor(kYellow);
  fQ2k[0]->Draw();
  fQ2k[1]->Draw("SAME");
  gPad->SetLogx();
  c15.Update();
  c15.cd(3);
  fYk[0]->SetStats(0);
  fYk[0]->SetTitle("y");
  fYk[0]->GetXaxis()->SetTitle("y");
  fYk[0]->GetYaxis()->SetTitle("Entries");
  fYk[1]->SetStats(0);
  fYk[1]->SetFillColor(kYellow);
  fYk[0]->Draw();
  fYk[1]->Draw("SAME");
  c15.Update();
  c15.cd(4);
  fThRich[0]->SetStats(0);
  fThRich[0]->SetTitle("#theta_{RICH}");
  fThRich[0]->GetXaxis()->SetTitle("#theta_{RICH}");
  fThRich[0]->GetYaxis()->SetTitle("Entries");
  fThRich[1]->SetStats(0);
  fThRich[1]->SetFillColor(kYellow);
  fThRich[0]->Draw();
  fThRich[1]->Draw("SAME");
  c15.Update();
  c15.cd(5);
  fPk[0]->SetStats(0);
  fPk[0]->SetTitle("p_{h}");
  fPk[0]->GetXaxis()->SetTitle("p_{h}");
  fPk[0]->GetYaxis()->SetTitle("Entries");
  fPk[1]->SetStats(0);
  fPk[1]->SetFillColor(kYellow);
  fPk[0]->Draw();
  fPk[1]->Draw("SAME");
  c15.Update();
  c15.cd(6);
  fZk[0]->SetStats(0);
  fZk[0]->SetTitle("z");
  fZk[0]->GetXaxis()->SetTitle("z");
  fZk[0]->GetYaxis()->SetTitle("Entries");
  fZk[1]->SetStats(0);
  fZk[1]->SetFillColor(kYellow);
  fZk[0]->Draw();
  fZk[1]->Draw("SAME");
  c15.Update();
  c16.cd(1);
  fQ2MT->SetStats(0);
  fQ2MT->GetXaxis()->SetTitle("Q^{2} (GeV/c)^{2}");
  fQ2MT->GetYaxis()->SetTitle("Entries");
  fQ2MT->SetLineColor(kRed);
  fQ2LT->SetStats(0);
  fQ2LT->SetLineColor(kMagenta);
  fQ2OT->SetStats(0);
  fQ2OT->SetLineColor(kOrange);
  fQ2LAST->SetStats(0);
  fQ2LAST->SetLineColor(kBlue);
  gPad->SetLogx();
  fQ2MT->Draw();
  fQ2LT->Draw("SAME");
  fQ2OT->Draw("SAME");
  fQ2LAST->Draw("SAME");
  c16.Update();
  c17.cd(1);
  fxMT->SetStats(0);
  fxMT->GetXaxis()->SetTitle("x");
  fxMT->GetYaxis()->SetTitle("Entries");
  fxMT->SetLineColor(kRed);
  fxLT->SetStats(0);
  fxLT->SetLineColor(kMagenta);
  fxOT->SetStats(0);
  fxOT->SetLineColor(kOrange);
  fxLAST->SetStats(0);
  fxLAST->SetLineColor(kBlue);
  gPad->SetLogx();
  fxMT->Draw();
  fxLT->Draw("SAME");
  fxOT->Draw("SAME");
  fxLAST->Draw("SAME");
  c17.Update();

  c1.Print("kinSIDIS.pdf(","pdf");
  c2.Print("kinSIDIS.pdf","pdf");
  c3.Print("kinSIDIS.pdf","pdf");
  c4.Print("kinSIDIS.pdf","pdf");
  c5.Print("kinSIDIS.pdf","pdf");
  c6.Print("kinSIDIS.pdf","pdf");
  c7.Print("kinSIDIS.pdf","pdf");
  c71.Print("kinSIDIS.pdf","pdf");
  c8.Print("kinSIDIS.pdf","pdf");
  c9.Print("kinSIDIS.pdf","pdf");
  c10.Print("kinSIDIS.pdf","pdf");
  c11.Print("kinSIDIS.pdf","pdf");
  c12.Print("kinSIDIS.pdf","pdf");
  c13.Print("kinSIDIS.pdf","pdf");
  c14.Print("kinSIDIS.pdf","pdf");
  c15.Print("kinSIDIS.pdf","pdf");
  c16.Print("kinSIDIS.pdf","pdf");
  c17.Print("kinSIDIS.pdf)","pdf");
}

//******************************************************************************
// RESETTING CONTAINERS
//******************************************************************************

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
          for(int m=0; m<2; m++)
          {
            fNDIS_evt[0][m][i][j][k]=0; fNDIS_evt_err[0][m][i][j][k]=0;
            fNDIS_evt[1][m][i][j][k]=0; fNDIS_evt_err[1][m][i][j][k]=0;
            fNDIS_evt[2][m][i][j][k]=0; fNDIS_evt_err[2][m][i][j][k]=0;
            fBinning[i][j][k].tab[c][m][0][0] = 0; fBinning[i][j][k].tab[c][m][1][0] = 0;
            fBinning[i][j][k].tab[c][m][0][1] = 0; fBinning[i][j][k].tab[c][m][1][1] = 0;
            fBinning[i][j][k].tab[c][m][0][2] = 0; fBinning[i][j][k].tab[c][m][1][2] = 0;
            fBinning[i][j][k].tab[c][m][0][3] = 0; fBinning[i][j][k].tab[c][m][1][3] = 0;
            fBinning_loose[i][j][k].tab[c][m][0][0] = 0; fBinning_severe[i][j][k].tab[c][m][0][0] = 0;
            fBinning_loose[i][j][k].tab[c][m][0][1] = 0; fBinning_severe[i][j][k].tab[c][m][0][1] = 0;
            fBinning_loose[i][j][k].tab[c][m][0][2] = 0; fBinning_severe[i][j][k].tab[c][m][0][2] = 0;
            fBinning_loose[i][j][k].tab[c][m][0][3] = 0; fBinning_severe[i][j][k].tab[c][m][0][3] = 0;
          }
          for(int ll=0; ll<4; ll++)
          {
            fMeanvalues_data[i][j][k].tab[c][ll][0]=0;
            fMeanvalues_data[i][j][k].tab[c][ll][1]=0;
            fMeanvalues_data[i][j][k].tab[c][ll][2]=0;
            fMeanvalues_data[i][j][k].tab[c][ll][3]=0;
            fMeanvalues_size[i][j][k].tab[c][ll][0]=0;
            fMeanvalues_size[i][j][k].tab[c][ll][1]=0;
            fMeanvalues_size[i][j][k].tab[c][ll][2]=0;
            fMeanvalues_size[i][j][k].tab[c][ll][3]=0;
            fMeanvalues[i][j][k].vec[c][ll][0].clear();
            fMeanvalues[i][j][k].vec[c][ll][1].clear();
            fMeanvalues[i][j][k].vec[c][ll][2].clear();
            fMeanvalues[i][j][k].vec[c][ll][3].clear();
          }

        }
      }
    }
  }
}

//******************************************************************************
// MAIN FUNCTION
//******************************************************************************

int main(int argc, char **argv)
{

  if(argc < 2)
  {
    cout << "ERROR : Not enough arguments." << endl;
    cout << "Asked : at least 1 *** Received : " << argc-1 << endl;
    cout << "./analySIDIS_split periodFile" << endl;

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
  Float_t Q2 = 0;
  Float_t xBj = 0;
  Float_t yBj = 0;
  Float_t zBj = 0;
  Float_t wBj = 0;
  Float_t nu = 0;

  if(kin_flag) create_kin_plots();

  if(RICH)
  {
    if(Y2006 || Y2012) load_rich_mat_2006(mat_RICH_2006_name, err_RICH_name);
    if(Y2016) load_rich_mat_2016(mat_RICH_2016_name, err_RICH_name);
  }
  else
  {
    load_rich_mat_dummy(mat_RICH_2016_name, err_RICH_name);
  }

  //cout << pi_sigma_uni[0][0] << " " << pi_sigma_uni[1][1] << " " << pi_sigma_uni[2][2] << endl;

  // Target cells
  if(Y2012) InitTargetFile(target_file_2012);
  else if(Y2016) InitTargetFile(target_file_2016);

  LoadInclusiveRadiativeCorrection();
  LoadSemiInclusiveRadiativeCorrection();
  LoadElectronCorrection();

  // ofstream test("test.txt", std::ofstream::out | std::ofstream::trunc);
  ofstream count("count.txt", std::ofstream::out | std::ofstream::trunc);

  //----------------------------------------------------------------------------
  //--------- nu cut prep ------------------------------------------------------
  //----------------------------------------------------------------------------

  for(int i=0; i<12; i++)
  {
    fNu_max[1][i] = sqrt(pow(MOMENTUM_UP,2)+pow(fM_K,2))/fZrange[i+1];
    fNu_min[1][i] = sqrt(pow(MOMENTUM_DOWN,2)+pow(fM_K,2))/fZrange[i];

    fNu_max[2][i] = sqrt(pow(MOMENTUM_UP,2)+pow(fM_p,2))/fZrange[i+1];
    fNu_min[2][i] = sqrt(pow(MOMENTUM_DOWN,2)+pow(fM_p,2))/fZrange[i];

    fNu_max[0][i] = sqrt(pow(MOMENTUM_UP,2)+pow(fM_pi,2))/fZrange[i+1];
    fNu_min[0][i] = sqrt(pow(MOMENTUM_DOWN,2)+pow(fM_pi,2))/fZrange[i];
  }

  // List of files
  ifstream periods(argv[1]);
  string filelist, periodName;
  int periodBit;
  while(periods >> periodName)
  {
    periods >> periodBit;
    fPeriodName.push_back(periodName);
    fPeriodBit.push_back(periodBit);
    if(!periodBit) continue;
    filelist = Form("%s/%s/filelist.txt",data_path,periodName.c_str());
    ifstream list(filelist);
    string filename;

    while(list >> filename)
    {
      TFile *f;

      cout << ".. Processing file " << filename << " .." << endl;
      f = TFile::Open(filename.c_str());

      if(!f) continue;

      fFilesNumber++;

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
      TBranch *inTarget = (TBranch*) tree->FindBranch("inTarget");
      TBranch *BMS = (TBranch*) tree->FindBranch("BMS");

      //Hadrons
      TBranch *p = (TBranch*) tree->FindBranch("Hadrons.P");
      TBranch *th = (TBranch*) tree->FindBranch("Hadrons.th");
      TBranch *pt = (TBranch*) tree->FindBranch("Hadrons.pt");
      TBranch *ph = (TBranch*) tree->FindBranch("Hadrons.ph");
      TBranch *hXX0 = (TBranch*) tree->FindBranch("Hadrons.XX0");
      TBranch *inHCALacc = (TBranch*) tree->FindBranch("Hadrons.inHCALacc");
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

      vector<Pvsz> Pvszlocal;
      vector<Pvsz> Pvszloose;
      vector<Pvsz> Pvszsevere;
      vector<Pvsz> Pvsz_errlocal;
      vector<Float_t> XBjlocal;
      vector<Float_t> YBjlocal;
      vector<Float_t> Q2local;
      vector<Float_t> Zvtxlocal;
      vector<Float_t> XBjloose;
      vector<Float_t> YBjloose;
      vector<Float_t> Q2loose;
      vector<Float_t> XBjsevere;
      vector<Float_t> YBjsevere;
      vector<Float_t> Q2severe;

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
        inTarget->GetEntry(ip);
        BMS->GetEntry(ip);

        //Hadrons
        p->GetEntry(ip);
        th->GetEntry(ip);
        pt->GetEntry(ip);
        ph->GetEntry(ip);
        hXX0->GetEntry(ip);
        inHCALacc->GetEntry(ip);
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
        //--------- Vertex Study ---------------------------------------------------
        //--------------------------------------------------------------------------

        Float_t zlab = z->GetLeaf("z")->GetValue();
        zlabbin=-1;

        if(Y2016)
        {
          if(-325<=zlab && zlab<-261.5) zlabbin = 0;
          else if(-261.5<=zlab && zlab<-198) zlabbin = 1;
          else if(-198<=zlab && zlab<-134.5) zlabbin = 2;
          else if(-134.5<=zlab && zlab<=-71) zlabbin = 3;
        }

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
        static const Float_t dz = 2;

        static const Float_t mcxU = -0.085;
        static const Float_t mcyU = 0.33;
        static const Float_t mczU_1 = -65+dz+4;
        //static const float mczU_2 = -35+dz;
        //static const float mczC_1 = -30+dz+8;
        //static const float mczC_2 = 30+dz;
        static const Float_t mcxD = -0.085;
        static const Float_t mcyD = 0.33;
        //static const float mczD_1 = 35+dz+2;
        static const Float_t mczD_2 = 65+dz;

        float mcR    = 1.4;
        //float mcyCUT = 1.4;

        //target position data 2006
        static const Float_t xU = -0.1;
        static const Float_t yU = 0.33;
        static const Float_t zU_1 = -65+dz+4;
        //static const float zU_2 = -35+dz;
        //static const float zC_1 = -30+dz+8;
        //static const float zC_2 =  30+dz;
        static const Float_t xD = -0.07;
        static const Float_t yD = 0.33;
        //static const float zD_1 =  35+dz+2;
        static const Float_t zD_2 =  65+dz;

        Float_t R    = 1.4;//1.4;
        Float_t yCUT = 1.4;

        Float_t mcxC = (mcxD-mcxU) * (mczU_1-z->GetLeaf("z")->GetValue()) / (mczU_1-mczD_2) + mcxU;
        Float_t mcyC = (mcyD-mcyU) * (mczU_1-z->GetLeaf("z")->GetValue()) / (mczU_1-mczD_2) + mcyU;
        Float_t mcr = sqrt( (x->GetLeaf("x")->GetValue()-mcxC)*(x->GetLeaf("x")->GetValue()-mcxC)
                      + (y->GetLeaf("y")->GetValue()-mcyC)*(y->GetLeaf("y")->GetValue()-mcyC) );
        Float_t xC = (xD-xU) * (zU_1-z->GetLeaf("z")->GetValue()) / (zU_1-zD_2) + xU;
        Float_t yC = (yD-yU) * (zU_1-z->GetLeaf("z")->GetValue()) / (zU_1-zD_2) + yU;
        Float_t r = sqrt( (x->GetLeaf("x")->GetValue()-xC)*(x->GetLeaf("x")->GetValue()-xC)
                      + (y->GetLeaf("y")->GetValue()-yC)*(y->GetLeaf("y")->GetValue()-yC) );

        // -------------------------------------------------------------------------
        // --------- DIS Selection -------------------------------------------------
        // -------------------------------------------------------------------------

        // Mu Charge
        if(MUCHARGE_SEPARATION)
        {
          if(!(Charge->GetLeaf("Charge")->GetValue()==MUCHARGE)) continue;
        }
        fMuCharge = Charge->GetLeaf("Charge")->GetValue()==1 ? 1 : 0;

        // Best Primary Vertex
        fBP++;

        // IsMuPrim
        if(!(0<isMuPrim->GetLeaf("isMuPrim")->GetValue())) continue;
        fRmu++;

        if(kin_flag)
        {
          int bin=-1;
          if(-325<z->GetLeaf("z")->GetValue() && z->GetLeaf("z")->GetValue()<-299.6) bin=0;
          else if(-299.6<z->GetLeaf("z")->GetValue() && z->GetLeaf("z")->GetValue()<-274.2) bin=1;
          else if(-274.2<z->GetLeaf("z")->GetValue() && z->GetLeaf("z")->GetValue()<-248.8) bin=2;
          else if(-248.8<z->GetLeaf("z")->GetValue() && z->GetLeaf("z")->GetValue()<-223.4) bin=3;
          else if(-223.4<z->GetLeaf("z")->GetValue() && z->GetLeaf("z")->GetValue()<-198) bin=4;
          else if(-198<z->GetLeaf("z")->GetValue() && z->GetLeaf("z")->GetValue()<-172.6) bin=5;
          else if(-172.6<z->GetLeaf("z")->GetValue() && z->GetLeaf("z")->GetValue()<-147.2) bin=6;
          else if(-147.2<z->GetLeaf("z")->GetValue() && z->GetLeaf("z")->GetValue()<-121.8) bin=7;
          else if(-121.8<z->GetLeaf("z")->GetValue() && z->GetLeaf("z")->GetValue()<-96.4) bin=8;
          else if(-96.4<z->GetLeaf("z")->GetValue() && z->GetLeaf("z")->GetValue()<-71) bin=9;
          fAllTarget[0]->Fill(z->GetLeaf("z")->GetValue(),x->GetLeaf("x")->GetValue());
          fAllTarget[1]->Fill(z->GetLeaf("z")->GetValue(),y->GetLeaf("y")->GetValue());
          if(bin!=-1) fAllTarget[bin+2]->Fill(x->GetLeaf("x")->GetValue(),y->GetLeaf("y")->GetValue());
          fZvtx[0]->Fill(z->GetLeaf("z")->GetValue());
          if(-325<z->GetLeaf("z")->GetValue() && z->GetLeaf("z")->GetValue()<-71) fZvtx[1]->Fill(z->GetLeaf("z")->GetValue());
        }

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
          if(!inTarget->GetLeaf("inTarget")->GetValue()) continue;
          if(!(-325<z->GetLeaf("z")->GetValue() && z->GetLeaf("z")->GetValue()<-71)) continue;
        }
        //2016 ---
        fTarg++;

        if(kin_flag)
        {
          int bin=-1;
          if(-325<z->GetLeaf("z")->GetValue() && z->GetLeaf("z")->GetValue()<-299.6) bin=0;
          else if(-299.6<z->GetLeaf("z")->GetValue() && z->GetLeaf("z")->GetValue()<-274.2) bin=1;
          else if(-274.2<z->GetLeaf("z")->GetValue() && z->GetLeaf("z")->GetValue()<-248.8) bin=2;
          else if(-248.8<z->GetLeaf("z")->GetValue() && z->GetLeaf("z")->GetValue()<-223.4) bin=3;
          else if(-223.4<z->GetLeaf("z")->GetValue() && z->GetLeaf("z")->GetValue()<-198) bin=4;
          else if(-198<z->GetLeaf("z")->GetValue() && z->GetLeaf("z")->GetValue()<-172.6) bin=5;
          else if(-172.6<z->GetLeaf("z")->GetValue() && z->GetLeaf("z")->GetValue()<-147.2) bin=6;
          else if(-147.2<z->GetLeaf("z")->GetValue() && z->GetLeaf("z")->GetValue()<-121.8) bin=7;
          else if(-121.8<z->GetLeaf("z")->GetValue() && z->GetLeaf("z")->GetValue()<-96.4) bin=8;
          else if(-96.4<z->GetLeaf("z")->GetValue() && z->GetLeaf("z")->GetValue()<-71) bin=9;
          fInTarget[0]->Fill(z->GetLeaf("z")->GetValue(),x->GetLeaf("x")->GetValue());
          fInTarget[1]->Fill(z->GetLeaf("z")->GetValue(),y->GetLeaf("y")->GetValue());
          fInTarget[bin+2]->Fill(x->GetLeaf("x")->GetValue(),y->GetLeaf("y")->GetValue());
          fZvtx[1]->Fill(z->GetLeaf("z")->GetValue());
        }

        // Energy of the muon beam
        if(!(140<E_beam->GetLeaf("E_beam")->GetValue() && E_beam->GetLeaf("E_beam")->GetValue()<180)) continue;
        fBEC++;

        // BMS
        if(!(BMS->GetLeaf("BMS")->GetValue()>3)) continue;
        fBMS++;

        // Chi2 beam
        if(!(beam_chi2->GetLeaf("beam_chi2")->GetValue()<10)) continue;
        fMuchi2++;

        // Cells crossing
        if(!(cellsCrossed->GetLeaf("cellsCrossed")->GetValue())) continue;
        fCell++;

        if(!(mu_prim_chi2->GetLeaf("mu_prim_chi2")->GetValue()<10)) continue;
        fMupchi2++;

        if(!(MZfirst->GetLeaf("MZfirst")->GetValue()<350)) continue;
        fMZfirst++;

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
          if(!(trig&2 || trig&4 || trig&8 || trig&512)) continue;
        }
        //2016 ---
        fTrig++;

        if(kin_flag) fQ2k[0]->Fill(Q2);

        // Q2 cut
        if(!(Q2>1)) continue;
        // if(!(Q2>0.85)) continue;
        fQ2test++;

        if(kin_flag)
        {
          fQ2k[1]->Fill(Q2);
          fYk[0]->Fill(yBj);
        }

        // y cut
        if(!(YMIN<yBj && yBj<YMAX)) continue;
        fYBjtest++;

        if(kin_flag) fYk[1]->Fill(yBj);

        // W cut
        if(!(5<sqrt(wBj) && sqrt(wBj)<17)) continue;
        fWBjtest++;

        // x cut
        if(!(XMIN<xBj && xBj<XMAX)) continue;
        fXBjtest++;

        if(kin_flag)
        {
          fQ2kin.push_back(Q2);
          fXBjkin.push_back(xBj);
          fYBjkin.push_back(yBj);
          fWBjkin.push_back(sqrt(wBj));
          fNukin.push_back(nu);
          fX.push_back(x->GetLeaf("x")->GetValue());
          fY.push_back(y->GetLeaf("y")->GetValue());
          fHO03x.push_back(HO03x->GetLeaf("HO03x")->GetValue());
          fHO03y.push_back(HO03y->GetLeaf("HO03y")->GetValue());
          fHO04x.push_back(HO04x->GetLeaf("HO04x")->GetValue());
          fHO04y.push_back(HO04y->GetLeaf("HO04y")->GetValue());
          if(trig&2)
          {
            fQ2MT->Fill(Q2);
            fxMT->Fill(xBj);
          }
          if(trig&4)
          {
            fQ2LT->Fill(Q2);
            fxLT->Fill(xBj);
          }
          if(trig&8)
          {
            fQ2OT->Fill(Q2);
            fxOT->Fill(xBj);
          }
          if(trig&512)
          {
            fQ2LAST->Fill(Q2);
            fxLAST->Fill(xBj);
          }
        }

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
        else if(0.5<yBj && yBj<0.7) ybin = 4;
        else ybin = 5;


        // z binning

        for(int i=0; i<12; i++)
        {
          fNDIS_evt[0][fMuCharge][xbin][ybin][i] += 1*GetInclusiveRadiativeCorrection(xBj,yBj);
          fNDIS_evt[1][fMuCharge][xbin][ybin][i] += 1*GetInclusiveRadiativeCorrection(xBj,yBj);
          fNDIS_evt[2][fMuCharge][xbin][ybin][i] += 1*GetInclusiveRadiativeCorrection(xBj,yBj);
          fNDIS_evt_zvtx[0][fMuCharge][xbin][ybin][i][zlabbin] += 1*GetInclusiveRadiativeCorrection(xBj,yBj);
          fNDIS_evt_zvtx[1][fMuCharge][xbin][ybin][i][zlabbin] += 1*GetInclusiveRadiativeCorrection(xBj,yBj);
          fNDIS_evt_zvtx[2][fMuCharge][xbin][ybin][i][zlabbin] += 1*GetInclusiveRadiativeCorrection(xBj,yBj);

          fNDIS_evt_err[0][fMuCharge][xbin][ybin][i] += pow(GetInclusiveRadiativeCorrection(xBj,yBj),2);
          fNDIS_evt_err[1][fMuCharge][xbin][ybin][i] += pow(GetInclusiveRadiativeCorrection(xBj,yBj),2);
          fNDIS_evt_err[2][fMuCharge][xbin][ybin][i] += pow(GetInclusiveRadiativeCorrection(xBj,yBj),2);
          fNDIS_evt_err_zvtx[0][fMuCharge][xbin][ybin][i][zlabbin] += pow(GetInclusiveRadiativeCorrection(xBj,yBj),2);
          fNDIS_evt_err_zvtx[1][fMuCharge][xbin][ybin][i][zlabbin] += pow(GetInclusiveRadiativeCorrection(xBj,yBj),2);
          fNDIS_evt_err_zvtx[2][fMuCharge][xbin][ybin][i][zlabbin] += pow(GetInclusiveRadiativeCorrection(xBj,yBj),2);

          fFlag[0][xbin][ybin][i]=0;
          fFlag[1][xbin][ybin][i]=0;
          fFlag[2][xbin][ybin][i]=0;

          // nu cut
          if(!(fNu_min[0][i]<nu && nu<fNu_max[0][i]))
          {
            fFlag[0][xbin][ybin][i]=1;
            fNDIS_evt[0][fMuCharge][xbin][ybin][i] -= GetInclusiveRadiativeCorrection(xBj,yBj);
            fNDIS_evt_zvtx[0][fMuCharge][xbin][ybin][i][zlabbin] -= GetInclusiveRadiativeCorrection(xBj,yBj);
            fNDIS_evt_err[0][fMuCharge][xbin][ybin][i] -= pow(GetInclusiveRadiativeCorrection(xBj,yBj),2);
            fNDIS_evt_err_zvtx[0][fMuCharge][xbin][ybin][i][zlabbin] -= pow(GetInclusiveRadiativeCorrection(xBj,yBj),2);
          }
          if(!(fNu_min[1][i]<nu && nu<fNu_max[1][i]))
          {
            fFlag[1][xbin][ybin][i]=1;
            fNDIS_evt[1][fMuCharge][xbin][ybin][i] -= GetInclusiveRadiativeCorrection(xBj,yBj);
            fNDIS_evt_zvtx[1][fMuCharge][xbin][ybin][i][zlabbin] -= GetInclusiveRadiativeCorrection(xBj,yBj);
            fNDIS_evt_err[1][fMuCharge][xbin][ybin][i] -= pow(GetInclusiveRadiativeCorrection(xBj,yBj),2);
            fNDIS_evt_err_zvtx[1][fMuCharge][xbin][ybin][i][zlabbin] -= pow(GetInclusiveRadiativeCorrection(xBj,yBj),2);
          }
          if(!(fNu_min[2][i]<nu && nu<fNu_max[2][i]))
          {
            fFlag[2][xbin][ybin][i]=1;
            fNDIS_evt[2][fMuCharge][xbin][ybin][i] -= GetInclusiveRadiativeCorrection(xBj,yBj);
            fNDIS_evt_zvtx[2][fMuCharge][xbin][ybin][i][zlabbin] -= GetInclusiveRadiativeCorrection(xBj,yBj);
            fNDIS_evt_err[2][fMuCharge][xbin][ybin][i] -= pow(GetInclusiveRadiativeCorrection(xBj,yBj),2);
            fNDIS_evt_err_zvtx[2][fMuCharge][xbin][ybin][i][zlabbin] -= pow(GetInclusiveRadiativeCorrection(xBj,yBj),2);
          }
        }

        // -------------------------------------------------------------------------
        // --------- Hadrons Selection ---------------------------------------------
        // -------------------------------------------------------------------------

        Pvsz pzcontainer;
        Pvsz pzcontainer_loose;
        Pvsz pzcontainer_severe;
        Pvsz pzcontainer_err;
        vector<Float_t> thlocal;
        vector<Float_t> ptlocal;
        hadiden hadcontainer;

        for(int i=0; i<p->GetLeaf("Hadrons.P")->GetLen(); i++)
        {

          fHadrons++;

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

            fLHsec_tab = new Float_t[5];
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

            fLHsec_tab = new Float_t[4];
            fLHsec_tab[0] = LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i);
            fLHsec_tab[1] = LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i);
            fLHsec_tab[2] = LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i);
            fLHsec_tab[3] = LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i);
            fusionSort(fLHsec_tab,4);
            fLHsec = fLHsec_tab[2];
          }

          set<Float_t>::iterator it = fLHsec_set.begin();
          advance(it, fLHsec_set.size()-2);
  #ifdef DEBUG
          if(*it != fLHsec) {cout << i << " : "
          <<  LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)  << " "
          << LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)  << " "
          << LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)  << " "
          << LH->GetLeaf("Hadrons.LH")->GetValue(3+6*i)  << " "
          << LH->GetLeaf("Hadrons.LH")->GetValue(4+6*i)  << " "
          << LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)  << " "
          << endl; cout << *it << " " << fLHsec << endl;}
  #endif

          //**********************************************************************

          // Hadron identification cuts ------------------------------------------

          // Charge +

          if(charge->GetLeaf("Hadrons.charge")->GetValue(i) == 1)
          {
            // Normal cuts ---

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

            // Normal cuts ---


            //**********************************************************************

            // Error associated to RICH unfolding ----------------------------------

            // Loose cuts ---

            if((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)>0)
               && (LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i))
               && (LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i))
               && (LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)/fLHsec>1.00)
               && (LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)>2.00)) fId_loose = 0;

            else if((LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)>0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i))
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i))
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)/fLHsec>1.06)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)>2.00)) fId_loose = 2;

            else if((8.9<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<=17.95-5)
                    && (((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)<2.3)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)<3.0))
                      || ((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i) == 0)
                      && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i) == 0)
                      && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i) == 0)
                      && (LH->GetLeaf("Hadrons.LH")->GetValue(3+6*i) == 0)
                      && (LH->GetLeaf("Hadrons.LH")->GetValue(4+6*i) == 0)
                      && (LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i) == 0)))) fId_loose = 4;

            else if((p->GetLeaf("Hadrons.P")->GetValue(i)>(17.95+5))
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i))
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i))
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i))
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)/fLHsec>0.98)) fId_loose = 4;

            else if(((17.95-5)<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<(17.95+5))
                    && (((LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i))
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i))
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i))
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)/fLHsec>0.98))
                    || (((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)<2.3)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)<3.0))
                    || ((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i) == 0)
                      && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i) == 0)
                      && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i) == 0)
                      && (LH->GetLeaf("Hadrons.LH")->GetValue(3+6*i) == 0)
                      && (LH->GetLeaf("Hadrons.LH")->GetValue(4+6*i) == 0)
                      && (LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i) == 0))))) fId_loose = 4;
            else fId_loose = 6;

            // Loose cuts ---

            // Severe cuts ---

            if((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)>0)
               && (LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i))
               && (LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i))
               && (LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)/fLHsec>1.06)
               && (LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)>2.04)) fId_severe = 0;

            else if((LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)>0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i))
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i))
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)/fLHsec>1.10)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)>2.16)) fId_severe = 2;

            else if((8.9<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<=17.95-5)
                    && (((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)<2.0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)<2.7))
                      || ((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i) == 0)
                      && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i) == 0)
                      && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i) == 0)
                      && (LH->GetLeaf("Hadrons.LH")->GetValue(3+6*i) == 0)
                      && (LH->GetLeaf("Hadrons.LH")->GetValue(4+6*i) == 0)
                      && (LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i) == 0)))) fId_severe = 4;

            else if((p->GetLeaf("Hadrons.P")->GetValue(i)>(17.95+5))
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i))
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i))
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i))
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)/fLHsec>1.06)) fId_severe = 4;

            else if(((17.95-5)<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<(17.95+5))
                    && (((LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i))
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i))
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i))
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)/fLHsec>1.06))
                    || (((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)<2.0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)<2.7))
                    || ((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i) == 0)
                      && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i) == 0)
                      && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i) == 0)
                      && (LH->GetLeaf("Hadrons.LH")->GetValue(3+6*i) == 0)
                      && (LH->GetLeaf("Hadrons.LH")->GetValue(4+6*i) == 0)
                      && (LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i) == 0))))) fId_severe = 4;
            else fId_severe = 6;

            // Severe cuts ---
          }

          // Charge -

          else if(charge->GetLeaf("Hadrons.charge")->GetValue(i) == -1)
          {
            // Normal cuts ---

            if((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)>0)
                && (LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i))
                && (LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i))
                && (LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)/fLHsec>1.02)
                && (LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)>2.02)) fId = 1;

            else if((LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)>0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i))
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i))
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

            // Normal cuts ---


            //**********************************************************************

            // Error associated to RICH unfolding ----------------------------------

            // Loose cuts ---

            if((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)>0)
               && (LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i))
               && (LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i))
               && (LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)/fLHsec>1.00)
               && (LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)>2.00)) fId_loose = 1;

            else if((LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)>0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i))
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i))
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)/fLHsec>1.06)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)>2.00)) fId_loose = 3;

            else if((8.9<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<=17.95-5)
                    && (((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)<2.3)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)<3.0))
                      || ((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i) == 0)
                      && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i) == 0)
                      && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i) == 0)
                      && (LH->GetLeaf("Hadrons.LH")->GetValue(3+6*i) == 0)
                      && (LH->GetLeaf("Hadrons.LH")->GetValue(4+6*i) == 0)
                      && (LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i) == 0)))) fId_loose = 5;

            else if((p->GetLeaf("Hadrons.P")->GetValue(i)>(17.95+5))
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i))
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i))
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i))
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)/fLHsec>0.98)) fId_loose = 5;

            else if(((17.95-5)<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<(17.95+5))
                    && (((LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i))
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i))
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i))
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)/fLHsec>0.98))
                    || (((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)<2.3)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)<3.0))
                    || ((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i) == 0)
                      && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i) == 0)
                      && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i) == 0)
                      && (LH->GetLeaf("Hadrons.LH")->GetValue(3+6*i) == 0)
                      && (LH->GetLeaf("Hadrons.LH")->GetValue(4+6*i) == 0)
                      && (LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i) == 0))))) fId_loose = 5;
            else fId_loose = 7;

            // Loose cuts ---

            // Severe cuts ---

            if((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)>0)
               && (LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i))
               && (LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i))
               && (LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)/fLHsec>1.06)
               && (LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)>2.04)) fId_severe = 1;

            else if((LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)>0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i))
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i))
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)/fLHsec>1.10)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)>2.16)) fId_severe = 3;

            else if((8.9<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<=17.95-5)
                    && (((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)<2.0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)<2.7))
                      || ((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i) == 0)
                      && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i) == 0)
                      && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i) == 0)
                      && (LH->GetLeaf("Hadrons.LH")->GetValue(3+6*i) == 0)
                      && (LH->GetLeaf("Hadrons.LH")->GetValue(4+6*i) == 0)
                      && (LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i) == 0)))) fId_severe = 5;

            else if((p->GetLeaf("Hadrons.P")->GetValue(i)>(17.95+5))
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i))
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i))
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i))
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)/fLHsec>1.06)) fId_severe = 5;

            else if(((17.95-5)<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<(17.95+5))
                    && (((LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i))
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i))
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)>LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i))
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i)/fLHsec>1.06))
                    || (((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)<2.0)
                    && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)/LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i)<2.7))
                    || ((LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i) == 0)
                      && (LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i) == 0)
                      && (LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i) == 0)
                      && (LH->GetLeaf("Hadrons.LH")->GetValue(3+6*i) == 0)
                      && (LH->GetLeaf("Hadrons.LH")->GetValue(4+6*i) == 0)
                      && (LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i) == 0))))) fId_severe = 5;
            else fId_severe = 7;

            // Severe cuts ---
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

          // Maximum radiation length cumulated
          if(!(hXX0->GetLeaf("Hadrons.XX0")->GetValue(i) < HXX0LIMIT)) continue;
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

          if(kin_flag) fThRich[0]->Fill(thRICH->GetLeaf("Hadrons.thRICH")->GetValue(i));

          // Theta cut
          if(!(0.01<thRICH->GetLeaf("Hadrons.thRICH")->GetValue(i) && thRICH->GetLeaf("Hadrons.thRICH")->GetValue(i)<0.12)) continue;
          fTRICH++;

          if(kin_flag) fThRich[1]->Fill(thRICH->GetLeaf("Hadrons.thRICH")->GetValue(i));

          // RICH position cut
          if(!(pow(RICHx->GetLeaf("Hadrons.RICHx")->GetValue(i),2)+pow(RICHy->GetLeaf("Hadrons.RICHy")->GetValue(i),2)>25)) continue;
          fPosRICH++;

          if(kin_flag)
          {
            fKinematicsRICH->Fill(p->GetLeaf("Hadrons.P")->GetValue(i),thC->GetLeaf("Hadrons.thC")->GetValue(i)*1000);
            fPk[0]->Fill(p->GetLeaf("Hadrons.P")->GetValue(i));
          }

          // Momentum cut (12 GeV to 40 GeV, increasing to 3 GeV to 40 GeV)
          if(!(MOMENTUM_DOWN<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<MOMENTUM_UP)) continue;
          fMom++;

          if(kin_flag) fPk[1]->Fill(p->GetLeaf("Hadrons.P")->GetValue(i));

          // Non null charge
          if(!charge->GetLeaf("Hadrons.charge")->GetValue(i)) continue;

          Int_t theta_bin, mom_bin;
          TMatrixD res_vect(3,1);
          Float_t res_vect_err[3];
          Float_t hadron_nb;

          // Theta and momentum binning

          if(0.01<thRICH->GetLeaf("Hadrons.thRICH")->GetValue(i) && thRICH->GetLeaf("Hadrons.thRICH")->GetValue(i)<0.04)
          {
            theta_bin = 0;
            if(MOMENTUM_DOWN<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<13) mom_bin = 0; // Here from 3 to 13 GeV
            else if(13<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<15) mom_bin = 1;
            else if(15<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<17) mom_bin = 2;
            else if(17<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<19) mom_bin = 3;
            else if(19<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<22) mom_bin = 4;
            else if(22<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<25) mom_bin = 5;
            else if(25<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<27) mom_bin = 6;
            else if(27<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<30) mom_bin = 7;
            else if(30<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<35) mom_bin = 8;
            else mom_bin = 9;
          }
          else
          {
            theta_bin = 1;
            if(MOMENTUM_DOWN<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<13) mom_bin = 0; // Here from 3 to 13 GeV
            else if(13<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<15) mom_bin = 1;
            else if(15<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<17) mom_bin = 2;
            else if(17<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<19) mom_bin = 3;
            else if(19<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<22) mom_bin = 4;
            else if(22<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<25) mom_bin = 5;
            else if(25<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<27) mom_bin = 6;
            else if(27<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<30) mom_bin = 7;
            else if(30<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<35) mom_bin = 8;
            else mom_bin = 9;
          }

          if(kin_flag) fZk[0]->Fill(zBj);

          // z cut
          if(!(0.2<zBj && zBj<0.85)) continue;
          fZtest++;

          int dz = abs(z->GetLeaf("z")->GetValue()-70);
          int ydy = y->GetLeaf("y")->GetValue()+dz*tan(th->GetLeaf("Hadrons.th")->GetValue(i))*sin(ph->GetLeaf("Hadrons.ph")->GetValue(i));
          int xdx = x->GetLeaf("x")->GetValue()+dz*tan(th->GetLeaf("Hadrons.th")->GetValue(i))*cos(ph->GetLeaf("Hadrons.ph")->GetValue(i));
          // if(!( ( -35 < xdx && xdx < 35 ) && ( -25 < ydy && ydy < 25 ) )) continue;

          // test << xBj << " " << Q2 << " " << yBj << " " << p->GetLeaf("Hadrons.P")->GetValue(i) << " "
          //      << thRICH->GetLeaf("Hadrons.thRICH")->GetValue(i) << " " << zBj << " " << fId << " "
          //      << LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i) << " "
          //      << LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i) << " "
          //      << LH->GetLeaf("Hadrons.LH")->GetValue(2+6*i) << " "
          //      << LH->GetLeaf("Hadrons.LH")->GetValue(3+6*i) << " "
          //      << LH->GetLeaf("Hadrons.LH")->GetValue(4+6*i) << " "
          //      << LH->GetLeaf("Hadrons.LH")->GetValue(5+6*i) << " " << endl;


          if(kin_flag)
          {
            fKinematics[3]->Fill(zBj);
            fZk[1]->Fill(zBj);
          }

          if(0.2<zBj && zBj<0.25) zbin = 0;
          else if(0.25<zBj && zBj<0.30) zbin = 1;
          else if(0.30<zBj && zBj<0.35) zbin = 2;
          else if(0.35<zBj && zBj<0.40) zbin = 3;
          else if(0.40<zBj && zBj<0.45) zbin = 4;
          else if(0.45<zBj && zBj<0.50) zbin = 5;
          else if(0.50<zBj && zBj<0.55) zbin = 6;
          else if(0.55<zBj && zBj<0.60) zbin = 7;
          else if(0.60<zBj && zBj<0.65) zbin = 8;
          else if(0.65<zBj && zBj<0.70) zbin = 9;
          else if(0.70<zBj && zBj<0.75) zbin = 10;
          else zbin = 11;


          if(kin_flag && zBj > 0.70 && 35<p->GetLeaf("Hadrons.P")->GetValue(i))
          {
            fLHpi.push_back(log10(LH->GetLeaf("Hadrons.LH")->GetValue(0+6*i)));
            fLHK.push_back(log10(LH->GetLeaf("Hadrons.LH")->GetValue(1+6*i)));
          }

          //**********************************************************************

          // Save of hadrons

          int hadron_flag = 0;

          if(fId==0)
          {
            if(!fFlag[0][xbin][ybin][zbin])
            {
              fHplus++; fPiplus++;
              res_vect = inv_rich_p[theta_bin][mom_bin]*pi_vect;
              for(int rce=0; rce<3; rce++) res_vect_err[rce] = pi_unfolding_err_p[theta_bin][mom_bin][rce];
              hadron_nb = 1;
              fPiplus_true += res_vect[0][0]; fKplus_true += res_vect[1][0]; fPplus_true += res_vect[2][0];
              fPiplus_err += pow(res_vect_err[0],2);
              fKplus_err += pow(res_vect_err[1],2);
              fPplus_err += pow(res_vect_err[2],2);
              res_vect[0][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect[1][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect[2][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect_err[0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect_err[1] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect_err[2] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              hadron_nb *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);

              pzcontainer.vec[1][0].push_back(zBj);
              pzcontainer.vec[1][1].push_back(res_vect[0][0]);
              pzcontainer.vec[1][2].push_back(res_vect[1][0]);
              pzcontainer.vec[1][3].push_back(res_vect[2][0]);
              thlocal.push_back(th->GetLeaf("Hadrons.th")->GetValue(i));
              ptlocal.push_back(pow(pt->GetLeaf("Hadrons.pt")->GetValue(i),2));

              pzcontainer_err.vec[1][0].push_back(zBj);
              pzcontainer_err.vec[1][1].push_back(pow(res_vect_err[0],2));
              pzcontainer_err.vec[1][2].push_back(pow(res_vect_err[1],2));
              pzcontainer_err.vec[1][3].push_back(pow(res_vect_err[2],2));

              pzcontainer.vec[1][4].push_back(hadron_nb);
              pzcontainer_err.vec[1][4].push_back(pow(hadron_nb,2));

              hadcontainer.vec.push_back(0);
  #ifdef DEBUG
              cout << res_vect[0][0] << " " << res_vect[1][0] << " " << res_vect[2][0] << endl;
  #endif
            }
          }
          else if(fId==1)
          {
            if(!fFlag[0][xbin][ybin][zbin])
            {
              fHminus++; fPiminus++;
              res_vect = inv_rich_m[theta_bin][mom_bin]*pi_vect;
              for(int rce=0; rce<3; rce++) res_vect_err[rce] = pi_unfolding_err_m[theta_bin][mom_bin][rce];
              hadron_nb = 1;
              fPiminus_true += res_vect[0][0]; fKminus_true += res_vect[1][0]; fPminus_true += res_vect[2][0];
              fPiminus_err += pow(res_vect_err[0],2);
              fKminus_err += pow(res_vect_err[1],2);
              fPminus_err += pow(res_vect_err[2],2);
              res_vect[0][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect[1][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect[2][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect_err[0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect_err[1] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect_err[2] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              hadron_nb *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);

              pzcontainer.vec[0][0].push_back(zBj);
              pzcontainer.vec[0][1].push_back(res_vect[0][0]);
              pzcontainer.vec[0][2].push_back(res_vect[1][0]);
              pzcontainer.vec[0][3].push_back(res_vect[2][0]);
              thlocal.push_back(th->GetLeaf("Hadrons.th")->GetValue(i));
              ptlocal.push_back(pow(pt->GetLeaf("Hadrons.pt")->GetValue(i),2));

              pzcontainer_err.vec[0][0].push_back(zBj);
              pzcontainer_err.vec[0][1].push_back(pow(res_vect_err[0],2));
              pzcontainer_err.vec[0][2].push_back(pow(res_vect_err[1],2));
              pzcontainer_err.vec[0][3].push_back(pow(res_vect_err[2],2));

              pzcontainer.vec[0][4].push_back(hadron_nb);
              pzcontainer_err.vec[0][4].push_back(pow(hadron_nb,2));

              hadcontainer.vec.push_back(1);
  #ifdef DEBUG
              cout << res_vect[0][0] << " " << res_vect[1][0] << " " << res_vect[2][0] << endl;
  #endif
            }
          }
          else if(fId==2)
          {
            if(!fFlag[0][xbin][ybin][zbin])
            {
              fHplus++;
              pzcontainer.vec[1][4].push_back(1*GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj));
              pzcontainer_err.vec[1][4].push_back(pow(GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj),2));
              thlocal.push_back(th->GetLeaf("Hadrons.th")->GetValue(i));
              ptlocal.push_back(pow(pt->GetLeaf("Hadrons.pt")->GetValue(i),2));
              hadron_flag = 1;
            }
            if(!fFlag[1][xbin][ybin][zbin])
            {
              fKplus++;
              res_vect = inv_rich_p[theta_bin][mom_bin]*k_vect;
              for(int rce=0; rce<3; rce++) res_vect_err[rce] = k_unfolding_err_p[theta_bin][mom_bin][rce];
              hadron_nb = 1;
              fPiplus_true += res_vect[0][0]; fKplus_true += res_vect[1][0]; fPplus_true += res_vect[2][0];
              fPiplus_err += pow(res_vect_err[0],2);
              fKplus_err += pow(res_vect_err[1],2);
              fPplus_err += pow(res_vect_err[2],2);
              res_vect[0][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect[1][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect[2][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect_err[0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect_err[1] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect_err[2] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              hadron_nb *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);

              pzcontainer.vec[1][0].push_back(zBj);
              pzcontainer.vec[1][1].push_back(res_vect[0][0]);
              pzcontainer.vec[1][2].push_back(res_vect[1][0]);
              pzcontainer.vec[1][3].push_back(res_vect[2][0]);

              pzcontainer_err.vec[1][0].push_back(zBj);
              pzcontainer_err.vec[1][1].push_back(pow(res_vect_err[0],2));
              pzcontainer_err.vec[1][2].push_back(pow(res_vect_err[1],2));
              pzcontainer_err.vec[1][3].push_back(pow(res_vect_err[2],2));

              if(!hadron_flag)
              {
                pzcontainer.vec[1][4].push_back(0);
                pzcontainer_err.vec[1][4].push_back(0);
                thlocal.push_back(th->GetLeaf("Hadrons.th")->GetValue(i));
                ptlocal.push_back(pow(pt->GetLeaf("Hadrons.pt")->GetValue(i),2));
              }

              hadcontainer.vec.push_back(2);
  #ifdef DEBUG
              cout << res_vect[0][0] << " " << res_vect[1][0] << " " << res_vect[2][0] << endl;
  #endif
            }
          }
          else if(fId==3)
          {
            if(!fFlag[0][xbin][ybin][zbin])
            {
              fHminus++;
              pzcontainer.vec[0][4].push_back(1*GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj));
              pzcontainer_err.vec[0][4].push_back(pow(GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj),2));
              thlocal.push_back(th->GetLeaf("Hadrons.th")->GetValue(i));
              ptlocal.push_back(pow(pt->GetLeaf("Hadrons.pt")->GetValue(i),2));
              hadron_flag = 1;
            }
            if(!fFlag[1][xbin][ybin][zbin])
            {
              fKminus++;
              res_vect = inv_rich_m[theta_bin][mom_bin]*k_vect;
              for(int rce=0; rce<3; rce++) res_vect_err[rce] = k_unfolding_err_m[theta_bin][mom_bin][rce];
              hadron_nb = 1;
              fPiminus_true += res_vect[0][0]; fKminus_true += res_vect[1][0]; fPminus_true += res_vect[2][0];
              fPiminus_err += pow(res_vect_err[0],2);
              fKminus_err += pow(res_vect_err[1],2);
              fPminus_err += pow(res_vect_err[2],2);
              res_vect[0][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect[1][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect[2][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect_err[0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect_err[1] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect_err[2] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              hadron_nb *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);

              pzcontainer.vec[0][0].push_back(zBj);
              pzcontainer.vec[0][1].push_back(res_vect[0][0]);
              pzcontainer.vec[0][2].push_back(res_vect[1][0]);
              pzcontainer.vec[0][3].push_back(res_vect[2][0]);

              pzcontainer_err.vec[0][0].push_back(zBj);
              pzcontainer_err.vec[0][1].push_back(pow(res_vect_err[0],2));
              pzcontainer_err.vec[0][2].push_back(pow(res_vect_err[1],2));
              pzcontainer_err.vec[0][3].push_back(pow(res_vect_err[2],2));

              if(!hadron_flag)
              {
                pzcontainer.vec[0][4].push_back(0);
                pzcontainer_err.vec[0][4].push_back(0);
                thlocal.push_back(th->GetLeaf("Hadrons.th")->GetValue(i));
                ptlocal.push_back(pow(pt->GetLeaf("Hadrons.pt")->GetValue(i),2));
              }

              hadcontainer.vec.push_back(3);
  #ifdef DEBUG
              cout << res_vect[0][0] << " " << res_vect[1][0] << " " << res_vect[2][0] << endl;
  #endif
            }
          }
          else if(fId==4)
          {
            if(!fFlag[0][xbin][ybin][zbin])
            {
              fHplus++;
              pzcontainer.vec[1][4].push_back(1*GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj));
              pzcontainer_err.vec[1][4].push_back(pow(GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj),2));
              thlocal.push_back(th->GetLeaf("Hadrons.th")->GetValue(i));
              ptlocal.push_back(pow(pt->GetLeaf("Hadrons.pt")->GetValue(i),2));
              hadron_flag = 1;
            }
            if(!fFlag[2][xbin][ybin][zbin])
            {
              fPplus++;
              res_vect = inv_rich_p[theta_bin][mom_bin]*p_vect;
              for(int rce=0; rce<3; rce++) res_vect_err[rce] = p_unfolding_err_p[theta_bin][mom_bin][rce];
              hadron_nb = 1;
              fPiplus_true += res_vect[0][0]; fKplus_true += res_vect[1][0]; fPplus_true += res_vect[2][0];
              fPiplus_err += pow(res_vect_err[0],2);
              fKplus_err += pow(res_vect_err[1],2);
              fPplus_err += pow(res_vect_err[2],2);
              res_vect[0][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect[1][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect[2][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect_err[0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect_err[1] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect_err[2] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              hadron_nb *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);

              pzcontainer.vec[1][0].push_back(zBj);
              pzcontainer.vec[1][1].push_back(res_vect[0][0]);
              pzcontainer.vec[1][2].push_back(res_vect[1][0]);
              pzcontainer.vec[1][3].push_back(res_vect[2][0]);

              pzcontainer_err.vec[1][0].push_back(zBj);
              pzcontainer_err.vec[1][1].push_back(pow(res_vect_err[0],2));
              pzcontainer_err.vec[1][2].push_back(pow(res_vect_err[1],2));
              pzcontainer_err.vec[1][3].push_back(pow(res_vect_err[2],2));

              if(!hadron_flag)
              {
                pzcontainer.vec[1][4].push_back(0);
                pzcontainer_err.vec[1][4].push_back(0);
                thlocal.push_back(th->GetLeaf("Hadrons.th")->GetValue(i));
                ptlocal.push_back(pow(pt->GetLeaf("Hadrons.pt")->GetValue(i),2));
              }

              hadcontainer.vec.push_back(0);
  #ifdef DEBUG
              cout << res_vect[0][0] << " " << res_vect[1][0] << " " << res_vect[2][0] << endl;
  #endif
            }
          }
          else if(fId==5)
          {
            if(!fFlag[0][xbin][ybin][zbin])
            {
              fHminus++;
              pzcontainer.vec[0][4].push_back(1*GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj));
              pzcontainer_err.vec[0][4].push_back(pow(GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj),2));
              thlocal.push_back(th->GetLeaf("Hadrons.th")->GetValue(i));
              ptlocal.push_back(pow(pt->GetLeaf("Hadrons.pt")->GetValue(i),2));
              hadron_flag = 1;
            }
            if(!fFlag[2][xbin][ybin][zbin])
            {
              fPminus++;
              res_vect = inv_rich_m[theta_bin][mom_bin]*p_vect;
              for(int rce=0; rce<3; rce++) res_vect_err[rce] = p_unfolding_err_m[theta_bin][mom_bin][rce];
              hadron_nb = 1;
              fPiminus_true += res_vect[0][0]; fKminus_true += res_vect[1][0]; fPminus_true += res_vect[2][0];
              fPiminus_err += pow(res_vect_err[0],2);
              fKminus_err += pow(res_vect_err[1],2);
              fPminus_err += pow(res_vect_err[2],2);
              res_vect[0][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect[1][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect[2][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect_err[0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect_err[1] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect_err[2] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              hadron_nb *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);

              pzcontainer.vec[0][0].push_back(zBj);
              pzcontainer.vec[0][1].push_back(res_vect[0][0]);
              pzcontainer.vec[0][2].push_back(res_vect[1][0]);
              pzcontainer.vec[0][3].push_back(res_vect[2][0]);

              pzcontainer_err.vec[0][0].push_back(zBj);
              pzcontainer_err.vec[0][1].push_back(pow(res_vect_err[0],2));
              pzcontainer_err.vec[0][2].push_back(pow(res_vect_err[1],2));
              pzcontainer_err.vec[0][3].push_back(pow(res_vect_err[2],2));

              if(!hadron_flag)
              {
                pzcontainer.vec[0][4].push_back(0);
                pzcontainer_err.vec[0][4].push_back(0);
                thlocal.push_back(th->GetLeaf("Hadrons.th")->GetValue(i));
                ptlocal.push_back(pow(pt->GetLeaf("Hadrons.pt")->GetValue(i),2));
              }

              hadcontainer.vec.push_back(5);
  #ifdef DEBUG
              cout << res_vect[0][0] << " " << res_vect[1][0] << " " << res_vect[2][0] << endl;
  #endif
            }
          }
          else if(fId==6)
          {
            if(!fFlag[0][xbin][ybin][zbin])
            {
              fHplus++;
              pzcontainer.vec[1][0].push_back(zBj); pzcontainer_err.vec[1][0].push_back(zBj);
              pzcontainer.vec[1][1].push_back(0); pzcontainer_err.vec[1][1].push_back(0);
              pzcontainer.vec[1][2].push_back(0); pzcontainer_err.vec[1][2].push_back(0);
              pzcontainer.vec[1][3].push_back(0); pzcontainer_err.vec[1][3].push_back(0);
              pzcontainer.vec[1][4].push_back(1*GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj));
              pzcontainer_err.vec[1][4].push_back(pow(GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj),2));
              thlocal.push_back(th->GetLeaf("Hadrons.th")->GetValue(i));
              ptlocal.push_back(pow(pt->GetLeaf("Hadrons.pt")->GetValue(i),2));
              hadcontainer.vec.push_back(6);
            }
          }
          else if(fId==7)
          {
            if(!fFlag[0][xbin][ybin][zbin])
            {
              fHminus++;
              pzcontainer.vec[0][0].push_back(zBj); pzcontainer_err.vec[0][0].push_back(zBj);
              pzcontainer.vec[0][1].push_back(0); pzcontainer_err.vec[0][1].push_back(0);
              pzcontainer.vec[0][2].push_back(0); pzcontainer_err.vec[0][2].push_back(0);
              pzcontainer.vec[0][3].push_back(0); pzcontainer_err.vec[0][3].push_back(0);
              pzcontainer.vec[0][4].push_back(1*GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj));
              pzcontainer_err.vec[0][4].push_back(pow(GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj),2));
              thlocal.push_back(th->GetLeaf("Hadrons.th")->GetValue(i));
              ptlocal.push_back(pow(pt->GetLeaf("Hadrons.pt")->GetValue(i),2));
              hadcontainer.vec.push_back(7);
            }
          }
          else
          {

          }


          //**********************************************************************

          // Loose cut

          hadron_flag = 0;

          if(fId_loose==0)
          {
            if(!fFlag[0][xbin][ybin][zbin])
            {
              res_vect = inv_rich_p[theta_bin][mom_bin]*pi_vect;
              hadron_nb = 1;
              res_vect[0][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect[1][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect[2][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              hadron_nb *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              pzcontainer_loose.vec[1][0].push_back(zBj);
              pzcontainer_loose.vec[1][1].push_back(res_vect[0][0]);
              pzcontainer_loose.vec[1][2].push_back(res_vect[1][0]);
              pzcontainer_loose.vec[1][3].push_back(res_vect[2][0]);
              pzcontainer_loose.vec[1][4].push_back(hadron_nb);
            }
          }
          else if(fId_loose==1)
          {
            if(!fFlag[0][xbin][ybin][zbin])
            {
              res_vect = inv_rich_m[theta_bin][mom_bin]*pi_vect;
              hadron_nb = 1;
              res_vect[0][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect[1][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect[2][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              hadron_nb *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              pzcontainer_loose.vec[0][0].push_back(zBj);
              pzcontainer_loose.vec[0][1].push_back(res_vect[0][0]);
              pzcontainer_loose.vec[0][2].push_back(res_vect[1][0]);
              pzcontainer_loose.vec[0][3].push_back(res_vect[2][0]);
              pzcontainer_loose.vec[0][4].push_back(hadron_nb);
            }
          }
          else if(fId_loose==2)
          {
            if(!fFlag[0][xbin][ybin][zbin])
            {
              pzcontainer_loose.vec[1][4].push_back(1*GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj));
              hadron_flag = 1;
            }
            if(!fFlag[1][xbin][ybin][zbin])
            {
              res_vect = inv_rich_p[theta_bin][mom_bin]*k_vect;
              res_vect[0][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect[1][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect[2][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              pzcontainer_loose.vec[1][0].push_back(zBj);
              pzcontainer_loose.vec[1][1].push_back(res_vect[0][0]);
              pzcontainer_loose.vec[1][2].push_back(res_vect[1][0]);
              pzcontainer_loose.vec[1][3].push_back(res_vect[2][0]);
              if(!hadron_flag) pzcontainer_loose.vec[1][4].push_back(0);
            }
          }
          else if(fId_loose==3)
          {
            if(!fFlag[0][xbin][ybin][zbin])
            {
              pzcontainer_loose.vec[0][4].push_back(1*GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj));
              hadron_flag = 1;
            }
            if(!fFlag[1][xbin][ybin][zbin])
            {
              res_vect = inv_rich_m[theta_bin][mom_bin]*k_vect;
              res_vect[0][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect[1][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect[2][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              pzcontainer_loose.vec[0][0].push_back(zBj);
              pzcontainer_loose.vec[0][1].push_back(res_vect[0][0]);
              pzcontainer_loose.vec[0][2].push_back(res_vect[1][0]);
              pzcontainer_loose.vec[0][3].push_back(res_vect[2][0]);
              if(!hadron_flag) pzcontainer_loose.vec[0][4].push_back(0);
            }
          }
          else if(fId_loose==4)
          {
            if(!fFlag[0][xbin][ybin][zbin])
            {
              pzcontainer_loose.vec[1][4].push_back(1*GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj));
              hadron_flag = 1;
            }
            if(!fFlag[2][xbin][ybin][zbin])
            {
              res_vect = inv_rich_p[theta_bin][mom_bin]*p_vect;
              res_vect[0][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect[1][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect[2][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              pzcontainer_loose.vec[1][0].push_back(zBj);
              pzcontainer_loose.vec[1][1].push_back(res_vect[0][0]);
              pzcontainer_loose.vec[1][2].push_back(res_vect[1][0]);
              pzcontainer_loose.vec[1][3].push_back(res_vect[2][0]);
              if(!hadron_flag) pzcontainer_loose.vec[1][4].push_back(0);
            }
          }
          else if(fId_loose==5)
          {
            if(!fFlag[0][xbin][ybin][zbin])
            {
              pzcontainer_loose.vec[0][4].push_back(1*GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj));
              hadron_flag = 1;
            }
            if(!fFlag[2][xbin][ybin][zbin])
            {
              res_vect = inv_rich_m[theta_bin][mom_bin]*p_vect;
              res_vect[0][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect[1][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect[2][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              pzcontainer_loose.vec[0][0].push_back(zBj);
              pzcontainer_loose.vec[0][1].push_back(res_vect[0][0]);
              pzcontainer_loose.vec[0][2].push_back(res_vect[1][0]);
              pzcontainer_loose.vec[0][3].push_back(res_vect[2][0]);
              if(!hadron_flag) pzcontainer_loose.vec[0][4].push_back(0);
            }
          }
          else if(fId_loose==6)
          {
            if(!fFlag[0][xbin][ybin][zbin])
            {
              pzcontainer_loose.vec[1][0].push_back(zBj);
              pzcontainer_loose.vec[1][1].push_back(0);
              pzcontainer_loose.vec[1][2].push_back(0);
              pzcontainer_loose.vec[1][3].push_back(0);
              hadron_nb = 1;
              hadron_nb *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              pzcontainer_loose.vec[1][4].push_back(hadron_nb);
            }
          }
          else if(fId_loose==7)
          {
            if(!fFlag[0][xbin][ybin][zbin])
            {
              pzcontainer_loose.vec[0][0].push_back(zBj);
              pzcontainer_loose.vec[0][1].push_back(0);
              pzcontainer_loose.vec[0][2].push_back(0);
              pzcontainer_loose.vec[0][3].push_back(0);
              hadron_nb = 1;
              hadron_nb *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              pzcontainer_loose.vec[0][4].push_back(hadron_nb);
            }
          }
          else
          {

          }


          // Severe cut

          hadron_flag = 0;

          if(fId_severe==0)
          {
            if(!fFlag[0][xbin][ybin][zbin])
            {
              res_vect = inv_rich_p[theta_bin][mom_bin]*pi_vect;
              hadron_nb = 1;
              res_vect[0][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect[1][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect[2][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              hadron_nb *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              pzcontainer_severe.vec[1][0].push_back(zBj);
              pzcontainer_severe.vec[1][1].push_back(res_vect[0][0]);
              pzcontainer_severe.vec[1][2].push_back(res_vect[1][0]);
              pzcontainer_severe.vec[1][3].push_back(res_vect[2][0]);
              pzcontainer_severe.vec[1][4].push_back(hadron_nb);
            }
          }
          else if(fId_severe==1)
          {
            if(!fFlag[0][xbin][ybin][zbin])
            {
              res_vect = inv_rich_m[theta_bin][mom_bin]*pi_vect;
              hadron_nb = 1;
              res_vect[0][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect[1][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect[2][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              hadron_nb *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              pzcontainer_severe.vec[0][0].push_back(zBj);
              pzcontainer_severe.vec[0][1].push_back(res_vect[0][0]);
              pzcontainer_severe.vec[0][2].push_back(res_vect[1][0]);
              pzcontainer_severe.vec[0][3].push_back(res_vect[2][0]);
              pzcontainer_severe.vec[0][4].push_back(hadron_nb);
            }
          }
          else if(fId_severe==2)
          {
            if(!fFlag[0][xbin][ybin][zbin])
            {
              pzcontainer_severe.vec[1][4].push_back(1*GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj));
              hadron_flag = 1;
            }
            if(!fFlag[1][xbin][ybin][zbin])
            {
              res_vect = inv_rich_p[theta_bin][mom_bin]*k_vect;
              res_vect[0][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect[1][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect[2][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              pzcontainer_severe.vec[1][0].push_back(zBj);
              pzcontainer_severe.vec[1][1].push_back(res_vect[0][0]);
              pzcontainer_severe.vec[1][2].push_back(res_vect[1][0]);
              pzcontainer_severe.vec[1][3].push_back(res_vect[2][0]);
              if(!hadron_flag) pzcontainer_severe.vec[1][4].push_back(0);
            }
          }
          else if(fId_severe==3)
          {
            if(!fFlag[0][xbin][ybin][zbin])
            {
              pzcontainer_severe.vec[0][4].push_back(1*GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj));
              hadron_flag = 1;
            }
            if(!fFlag[1][xbin][ybin][zbin])
            {
              res_vect = inv_rich_m[theta_bin][mom_bin]*k_vect;
              res_vect[0][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect[1][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect[2][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              pzcontainer_severe.vec[0][0].push_back(zBj);
              pzcontainer_severe.vec[0][1].push_back(res_vect[0][0]);
              pzcontainer_severe.vec[0][2].push_back(res_vect[1][0]);
              pzcontainer_severe.vec[0][3].push_back(res_vect[2][0]);
              if(!hadron_flag) pzcontainer_severe.vec[0][4].push_back(0);
            }
          }
          else if(fId_severe==4)
          {
            if(!fFlag[0][xbin][ybin][zbin])
            {
              pzcontainer_severe.vec[1][4].push_back(1*GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj));
              hadron_flag = 1;
            }
            if(!fFlag[2][xbin][ybin][zbin])
            {
              res_vect = inv_rich_p[theta_bin][mom_bin]*p_vect;
              res_vect[0][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect[1][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect[2][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              pzcontainer_severe.vec[1][0].push_back(zBj);
              pzcontainer_severe.vec[1][1].push_back(res_vect[0][0]);
              pzcontainer_severe.vec[1][2].push_back(res_vect[1][0]);
              pzcontainer_severe.vec[1][3].push_back(res_vect[2][0]);
              if(!hadron_flag) pzcontainer_severe.vec[1][4].push_back(0);
            }
          }
          else if(fId_severe==5)
          {
            if(!fFlag[0][xbin][ybin][zbin])
            {
              pzcontainer_severe.vec[0][4].push_back(1*GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj));
              hadron_flag = 1;
            }
            if(!fFlag[2][xbin][ybin][zbin])
            {
              res_vect = inv_rich_m[theta_bin][mom_bin]*p_vect;
              res_vect[0][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect[1][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              res_vect[2][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              pzcontainer_severe.vec[0][0].push_back(zBj);
              pzcontainer_severe.vec[0][1].push_back(res_vect[0][0]);
              pzcontainer_severe.vec[0][2].push_back(res_vect[1][0]);
              pzcontainer_severe.vec[0][3].push_back(res_vect[2][0]);
              if(!hadron_flag) pzcontainer_severe.vec[0][4].push_back(0);
            }
          }
          else if(fId_severe==6)
          {
            if(!fFlag[0][xbin][ybin][zbin])
            {
              pzcontainer_severe.vec[1][0].push_back(zBj);
              pzcontainer_severe.vec[1][1].push_back(0);
              pzcontainer_severe.vec[1][2].push_back(0);
              pzcontainer_severe.vec[1][3].push_back(0);
              hadron_nb = 1;
              hadron_nb *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              pzcontainer_severe.vec[1][4].push_back(hadron_nb);
            }
          }
          else if(fId_severe==7)
          {
            if(!fFlag[0][xbin][ybin][zbin])
            {
              pzcontainer_severe.vec[0][0].push_back(zBj);
              pzcontainer_severe.vec[0][1].push_back(0);
              pzcontainer_severe.vec[0][2].push_back(0);
              pzcontainer_severe.vec[0][3].push_back(0);
              hadron_nb = 1;
              hadron_nb *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
              pzcontainer_severe.vec[0][4].push_back(hadron_nb);
            }
          }
          else
          {

          }
        }

        //Misc
        fQ2.push_back(Q2);
        fXBj.push_back(xBj);
        fYBj.push_back(yBj);
        fWBj.push_back(wBj);
        fNu.push_back(nu);

        fPvsz.push_back(pzcontainer);
        fPvsz_err.push_back(pzcontainer_err);
        fHadiden.push_back(hadcontainer);
        fTheta.push_back(thlocal);
        fpT.push_back(ptlocal);

        Q2local.push_back(Q2);
        Pvszlocal.push_back(pzcontainer);
        Pvsz_errlocal.push_back(pzcontainer_err);
        XBjlocal.push_back(xBj);
        YBjlocal.push_back(yBj);
        Zvtxlocal.push_back(zlab);

        Q2loose.push_back(Q2);
        Pvszloose.push_back(pzcontainer_loose);
        XBjloose.push_back(xBj);
        YBjloose.push_back(yBj);

        Q2severe.push_back(Q2);
        Pvszsevere.push_back(pzcontainer_severe);
        XBjsevere.push_back(xBj);
        YBjsevere.push_back(yBj);

      }
      cout << "\n" << endl;

      // Loose cut

      for(int i=0; i<int(Q2loose.size()); i++)
      {
        if(0.004<=XBjloose[i] && XBjloose[i]<0.01) xbin = 0;
        else if(0.01<=XBjloose[i] && XBjloose[i]<0.02) xbin = 1;
        else if(0.02<=XBjloose[i] && XBjloose[i]<0.03) xbin = 2;
        else if(0.03<=XBjloose[i] && XBjloose[i]<0.04) xbin = 3;
        else if(0.04<=XBjloose[i] && XBjloose[i]<0.06) xbin = 4;
        else if(0.06<=XBjloose[i] && XBjloose[i]<0.1) xbin = 5;
        else if(0.1<=XBjloose[i] && XBjloose[i]<0.14) xbin = 6;
        else if(0.1<=XBjloose[i] && XBjloose[i]<0.18) xbin = 7;
        else xbin = 8;

        if(0.1<YBjloose[i] && YBjloose[i]<0.15) ybin = 0;
        else if(0.15<YBjloose[i] && YBjloose[i]<0.2) ybin = 1;
        else if(0.2<YBjloose[i] && YBjloose[i]<0.3) ybin = 2;
        else if(0.3<YBjloose[i] && YBjloose[i]<0.5) ybin = 3;
        else if(0.5<YBjloose[i] && YBjloose[i]<0.7) ybin = 4;
        else ybin = 5;

        for(int j=0; j<2; j++)
        {
          for(int l=0; l<int(Pvszloose[i].vec[j][0].size()); l++)
          {
            if(0.2<Pvszloose[i].vec[j][0][l] && Pvszloose[i].vec[j][0][l]<0.25) zbin = 0;
            else if(0.25<Pvszloose[i].vec[j][0][l] && Pvszloose[i].vec[j][0][l]<0.30) zbin = 1;
            else if(0.30<Pvszloose[i].vec[j][0][l] && Pvszloose[i].vec[j][0][l]<0.35) zbin = 2;
            else if(0.35<Pvszloose[i].vec[j][0][l] && Pvszloose[i].vec[j][0][l]<0.40) zbin = 3;
            else if(0.40<Pvszloose[i].vec[j][0][l] && Pvszloose[i].vec[j][0][l]<0.45) zbin = 4;
            else if(0.45<Pvszloose[i].vec[j][0][l] && Pvszloose[i].vec[j][0][l]<0.50) zbin = 5;
            else if(0.50<Pvszloose[i].vec[j][0][l] && Pvszloose[i].vec[j][0][l]<0.55) zbin = 6;
            else if(0.55<Pvszloose[i].vec[j][0][l] && Pvszloose[i].vec[j][0][l]<0.60) zbin = 7;
            else if(0.60<Pvszloose[i].vec[j][0][l] && Pvszloose[i].vec[j][0][l]<0.65) zbin = 8;
            else if(0.65<Pvszloose[i].vec[j][0][l] && Pvszloose[i].vec[j][0][l]<0.70) zbin = 9;
            else if(0.70<Pvszloose[i].vec[j][0][l] && Pvszloose[i].vec[j][0][l]<0.75) zbin = 10;
            else zbin = 11;

            for(int ll=0; ll<4; ll++)
            {
              fBinning_loose[xbin][ybin][zbin].tab[j][fMuCharge][0][ll] += Pvszloose[i].vec[j][ll+1][l];
            }
          }
        }

      }

      // Severe cut

      for(int i=0; i<int(Q2severe.size()); i++)
      {
        if(0.004<=XBjsevere[i] && XBjsevere[i]<0.01) xbin = 0;
        else if(0.01<=XBjsevere[i] && XBjsevere[i]<0.02) xbin = 1;
        else if(0.02<=XBjsevere[i] && XBjsevere[i]<0.03) xbin = 2;
        else if(0.03<=XBjsevere[i] && XBjsevere[i]<0.04) xbin = 3;
        else if(0.04<=XBjsevere[i] && XBjsevere[i]<0.06) xbin = 4;
        else if(0.06<=XBjsevere[i] && XBjsevere[i]<0.1) xbin = 5;
        else if(0.1<=XBjsevere[i] && XBjsevere[i]<0.14) xbin = 6;
        else if(0.1<=XBjsevere[i] && XBjsevere[i]<0.18) xbin = 7;
        else xbin = 8;

        if(0.1<YBjsevere[i] && YBjsevere[i]<0.15) ybin = 0;
        else if(0.15<YBjsevere[i] && YBjsevere[i]<0.2) ybin = 1;
        else if(0.2<YBjsevere[i] && YBjsevere[i]<0.3) ybin = 2;
        else if(0.3<YBjsevere[i] && YBjsevere[i]<0.5) ybin = 3;
        else if(0.5<YBjsevere[i] && YBjsevere[i]<0.7) ybin = 4;
        else ybin = 5;

        for(int j=0; j<2; j++)
        {
          for(int l=0; l<int(Pvszsevere[i].vec[j][0].size()); l++)
          {
            if(0.2<Pvszsevere[i].vec[j][0][l] && Pvszsevere[i].vec[j][0][l]<0.25) zbin = 0;
            else if(0.25<Pvszsevere[i].vec[j][0][l] && Pvszsevere[i].vec[j][0][l]<0.30) zbin = 1;
            else if(0.30<Pvszsevere[i].vec[j][0][l] && Pvszsevere[i].vec[j][0][l]<0.35) zbin = 2;
            else if(0.35<Pvszsevere[i].vec[j][0][l] && Pvszsevere[i].vec[j][0][l]<0.40) zbin = 3;
            else if(0.40<Pvszsevere[i].vec[j][0][l] && Pvszsevere[i].vec[j][0][l]<0.45) zbin = 4;
            else if(0.45<Pvszsevere[i].vec[j][0][l] && Pvszsevere[i].vec[j][0][l]<0.50) zbin = 5;
            else if(0.50<Pvszsevere[i].vec[j][0][l] && Pvszsevere[i].vec[j][0][l]<0.55) zbin = 6;
            else if(0.55<Pvszsevere[i].vec[j][0][l] && Pvszsevere[i].vec[j][0][l]<0.60) zbin = 7;
            else if(0.60<Pvszsevere[i].vec[j][0][l] && Pvszsevere[i].vec[j][0][l]<0.65) zbin = 8;
            else if(0.65<Pvszsevere[i].vec[j][0][l] && Pvszsevere[i].vec[j][0][l]<0.70) zbin = 9;
            else if(0.70<Pvszsevere[i].vec[j][0][l] && Pvszsevere[i].vec[j][0][l]<0.75) zbin = 10;
            else zbin = 11;

            for(int ll=0; ll<4; ll++)
            {
              fBinning_severe[xbin][ybin][zbin].tab[j][fMuCharge][0][ll] += Pvszsevere[i].vec[j][ll+1][l];
            }
          }
        }

      }

      for(int i=0; i<int(Q2local.size()); i++)
      {
        if(0.004<=XBjlocal[i] && XBjlocal[i]<0.01) xbin = 0;
        else if(0.01<=XBjlocal[i] && XBjlocal[i]<0.02) xbin = 1;
        else if(0.02<=XBjlocal[i] && XBjlocal[i]<0.03) xbin = 2;
        else if(0.03<=XBjlocal[i] && XBjlocal[i]<0.04) xbin = 3;
        else if(0.04<=XBjlocal[i] && XBjlocal[i]<0.06) xbin = 4;
        else if(0.06<=XBjlocal[i] && XBjlocal[i]<0.1) xbin = 5;
        else if(0.1<=XBjlocal[i] && XBjlocal[i]<0.14) xbin = 6;
        else if(0.1<=XBjlocal[i] && XBjlocal[i]<0.18) xbin = 7;
        else xbin = 8;

        if(0.1<YBjlocal[i] && YBjlocal[i]<0.15) ybin = 0;
        else if(0.15<YBjlocal[i] && YBjlocal[i]<0.2) ybin = 1;
        else if(0.2<YBjlocal[i] && YBjlocal[i]<0.3) ybin = 2;
        else if(0.3<YBjlocal[i] && YBjlocal[i]<0.5) ybin = 3;
        else if(0.5<YBjlocal[i] && YBjlocal[i]<0.7) ybin = 4;
        else ybin = 5;

        if(-325<=Zvtxlocal[i] && Zvtxlocal[i]<-261.5) zlabbin = 0;
        else if(-261.5<=Zvtxlocal[i] && Zvtxlocal[i]<-198) zlabbin = 1;
        else if(-198<=Zvtxlocal[i] && Zvtxlocal[i]<-134.5) zlabbin = 2;
        else if(-134.5<=Zvtxlocal[i] && Zvtxlocal[i]<=-71) zlabbin = 3;

        for(int j=0; j<2; j++)
        {
          for(int l=0; l<int(Pvszlocal[i].vec[j][0].size()); l++)
          {
            if(0.2<Pvszlocal[i].vec[j][0][l] && Pvszlocal[i].vec[j][0][l]<0.25) zbin = 0;
            else if(0.25<Pvszlocal[i].vec[j][0][l] && Pvszlocal[i].vec[j][0][l]<0.30) zbin = 1;
            else if(0.30<Pvszlocal[i].vec[j][0][l] && Pvszlocal[i].vec[j][0][l]<0.35) zbin = 2;
            else if(0.35<Pvszlocal[i].vec[j][0][l] && Pvszlocal[i].vec[j][0][l]<0.40) zbin = 3;
            else if(0.40<Pvszlocal[i].vec[j][0][l] && Pvszlocal[i].vec[j][0][l]<0.45) zbin = 4;
            else if(0.45<Pvszlocal[i].vec[j][0][l] && Pvszlocal[i].vec[j][0][l]<0.50) zbin = 5;
            else if(0.50<Pvszlocal[i].vec[j][0][l] && Pvszlocal[i].vec[j][0][l]<0.55) zbin = 6;
            else if(0.55<Pvszlocal[i].vec[j][0][l] && Pvszlocal[i].vec[j][0][l]<0.60) zbin = 7;
            else if(0.60<Pvszlocal[i].vec[j][0][l] && Pvszlocal[i].vec[j][0][l]<0.65) zbin = 8;
            else if(0.65<Pvszlocal[i].vec[j][0][l] && Pvszlocal[i].vec[j][0][l]<0.70) zbin = 9;
            else if(0.70<Pvszlocal[i].vec[j][0][l] && Pvszlocal[i].vec[j][0][l]<0.75) zbin = 10;
            else zbin = 11;

            if(0<=fTheta[i][l] && fTheta[i][l]<0.015) thbin = 0;
            else if(0.015<=fTheta[i][l] && fTheta[i][l]<0.025) thbin = 1;
            else if(0.025<=fTheta[i][l] && fTheta[i][l]<0.035) thbin = 2;
            else if(0.035<=fTheta[i][l] && fTheta[i][l]<0.045) thbin = 3;
            else if(0.045<=fTheta[i][l] && fTheta[i][l]<0.058) thbin = 4;
            else if(0.058<=fTheta[i][l] && fTheta[i][l]<0.072) thbin = 5;
            else if(0.072<=fTheta[i][l] && fTheta[i][l]<0.088) thbin = 6;
            else if(0.088<=fTheta[i][l] && fTheta[i][l]<0.2) thbin = 7;
            else thbin = -1;

            if(0.02<=fpT[i][l] && fpT[i][l]<0.08) ptbin = 0;
            else if(0.08<=fpT[i][l] && fpT[i][l]<0.14) ptbin = 1;
            else if(0.14<=fpT[i][l] && fpT[i][l]<0.23) ptbin = 2;
            else if(0.23<=fpT[i][l] && fpT[i][l]<0.35) ptbin = 3;
            else if(0.35<=fpT[i][l] && fpT[i][l]<0.52) ptbin = 4;
            else if(0.52<=fpT[i][l] && fpT[i][l]<0.76) ptbin = 5;
            else if(0.76<=fpT[i][l] && fpT[i][l]<1.12) ptbin = 6;
            else if(1.12<=fpT[i][l] && fpT[i][l]<1.52) ptbin = 7;
            else if(1.52<=fpT[i][l] && fpT[i][l]<2.05) ptbin = 8;
            else if(2.05<=fpT[i][l] && fpT[i][l]<3.0) ptbin = 9;
            else ptbin = -1;

            for(int ll=0; ll<4; ll++)
            {
              fBinning[xbin][ybin][zbin].tab[j][fMuCharge][0][ll] += Pvszlocal[i].vec[j][ll+1][l];
              fBinning_zvtx[xbin][ybin][zbin][zlabbin].tab[j][fMuCharge][0][ll] += Pvszlocal[i].vec[j][ll+1][l];
              if(thbin!=-1) fBinning_theta[xbin][ybin][zbin][thbin].tab[j][fMuCharge][0][ll] += Pvszlocal[i].vec[j][ll+1][l];
              if(ptbin!=-1) fBinning_pt[xbin][ybin][zbin][ptbin].tab[j][fMuCharge][0][ll] += Pvszlocal[i].vec[j][ll+1][l];
              fBinning[xbin][ybin][zbin].tab[j][fMuCharge][1][ll] += Pvsz_errlocal[i].vec[j][ll+1][l];
              fBinning_zvtx[xbin][ybin][zbin][zlabbin].tab[j][fMuCharge][1][ll] += Pvsz_errlocal[i].vec[j][ll+1][l];
              if(thbin!=-1) fBinning_theta[xbin][ybin][zbin][thbin].tab[j][fMuCharge][1][ll] += Pvsz_errlocal[i].vec[j][ll+1][l];
              if(ptbin!=-1) fBinning_pt[xbin][ybin][zbin][ptbin].tab[j][fMuCharge][1][ll] += Pvsz_errlocal[i].vec[j][ll+1][l];
              fMeanvalues[xbin][ybin][zbin].vec[j][ll][2].push_back(Q2local[i]);
              fMeanvalues[xbin][ybin][zbin].vec[j][ll][0].push_back(XBjlocal[i]);
              fMeanvalues[xbin][ybin][zbin].vec[j][ll][1].push_back(YBjlocal[i]);
              fMeanvalues[xbin][ybin][zbin].vec[j][ll][3].push_back(Pvszlocal[i].vec[j][0][l]);
            }
          }
        }
      }

      Pvszlocal.clear();
      Pvsz_errlocal.clear();
      XBjlocal.clear();
      YBjlocal.clear();
      Q2local.clear();
      Zvtxlocal.clear();
      Pvszloose.clear();
      XBjloose.clear();
      YBjloose.clear();
      Q2loose.clear();
      Pvszsevere.clear();
      XBjsevere.clear();
      YBjsevere.clear();
      Q2severe.clear();
      fTheta.clear();
      fpT.clear();
      f->Close();
    }

    for(int i=0; i<9; i++)
    {
      for(int j=0; j<6; j++)
      {
        for(int k=0; k<12; k++)
        {
          for(int c=0; c<2; c++)
          {
            for(int ll=0; ll<4; ll++)
            {
              for(int sv=0; sv<int(fMeanvalues[i][j][k].vec[c][ll][0].size()); sv++)
              {
                fMeanvalues_data[i][j][k].tab[c][ll][0] += fMeanvalues[i][j][k].vec[c][ll][0][sv];
                fMeanvalues_data[i][j][k].tab[c][ll][1] += fMeanvalues[i][j][k].vec[c][ll][1][sv];
                fMeanvalues_data[i][j][k].tab[c][ll][2] += fMeanvalues[i][j][k].vec[c][ll][2][sv];
                fMeanvalues_data[i][j][k].tab[c][ll][3] += fMeanvalues[i][j][k].vec[c][ll][3][sv];
              }

              if(int(fMeanvalues[i][j][k].vec[c][ll][0].size()))
              {
                fMeanvalues_size[i][j][k].tab[c][ll][0] = int(fMeanvalues[i][j][k].vec[c][ll][0].size());
                fMeanvalues_size[i][j][k].tab[c][ll][1] = int(fMeanvalues[i][j][k].vec[c][ll][1].size());
                fMeanvalues_size[i][j][k].tab[c][ll][2] = int(fMeanvalues[i][j][k].vec[c][ll][2].size());
                fMeanvalues_size[i][j][k].tab[c][ll][3] = int(fMeanvalues[i][j][k].vec[c][ll][3].size());
              }
            }
          }
        }
      }
    }

    ofstream ofs_h(Form("rawmult/%d/hadron_%s.txt",year,periodName.c_str()), std::ofstream::out | std::ofstream::trunc);
    ofstream ofs_hzvtx(Form("rawmult/%d/hadron_zvtx_%s.txt",year,periodName.c_str()), std::ofstream::out | std::ofstream::trunc);
    ofstream ofs_hth(Form("rawmult/%d/hadron_theta_%s.txt",year,periodName.c_str()), std::ofstream::out | std::ofstream::trunc);
    ofstream ofs_hpt(Form("rawmult/%d/hadron_pT_%s.txt",year,periodName.c_str()), std::ofstream::out | std::ofstream::trunc);
    ofstream ofs_d(Form("rawmult/%d/DIS_%s.txt",year,periodName.c_str()), std::ofstream::out | std::ofstream::trunc);
    ofstream ofs_dzvtx(Form("rawmult/%d/DIS_zvtx_%s.txt",year,periodName.c_str()), std::ofstream::out | std::ofstream::trunc);
    ofstream xc(Form("rawmult/%d/xcheck_%s.txt",year,periodName.c_str()), std::ofstream::out | std::ofstream::trunc);

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
              ofs_d << fNDIS_evt[0][1][i][j][k] << " " << fNDIS_evt_err[0][1][i][j][k] << " " << fNDIS_evt[1][1][i][j][k] << " " << fNDIS_evt_err[1][1][i][j][k] << " " << fNDIS_evt[2][1][i][j][k] << " " << fNDIS_evt_err[2][1][i][j][k] << " " <<
                       fNDIS_evt[0][0][i][j][k] << " " << fNDIS_evt_err[0][0][i][j][k] << " " << fNDIS_evt[1][0][i][j][k] << " " << fNDIS_evt_err[1][0][i][j][k] << " " << fNDIS_evt[2][0][i][j][k] << " " << fNDIS_evt_err[2][0][i][j][k];
              for(int zv=0; zv<4; zv++)
                ofs_dzvtx << fNDIS_evt_zvtx[0][1][i][j][k][zv] << " " << fNDIS_evt_err_zvtx[0][1][i][j][k][zv] << " " << fNDIS_evt_zvtx[1][1][i][j][k][zv] << " " << fNDIS_evt_err_zvtx[1][1][i][j][k][zv] << " " << fNDIS_evt_zvtx[2][1][i][j][k][zv] << " " << fNDIS_evt_err_zvtx[2][1][i][j][k][zv] << " " <<
                             fNDIS_evt_zvtx[0][0][i][j][k][zv] << " " << fNDIS_evt_err_zvtx[0][0][i][j][k][zv] << " " << fNDIS_evt_zvtx[1][0][i][j][k][zv] << " " << fNDIS_evt_err_zvtx[1][0][i][j][k][zv] << " " << fNDIS_evt_zvtx[2][0][i][j][k][zv] << " " << fNDIS_evt_err_zvtx[2][0][i][j][k][zv] << " ";
            }

            for(int ll=0; ll<4; ll++)
            {
              ofs_d << " " << fMeanvalues_data[i][j][k].tab[c][ll][0] << " " <<
                              fMeanvalues_data[i][j][k].tab[c][ll][1] << " " <<
                              fMeanvalues_data[i][j][k].tab[c][ll][2] << " " <<
                              fMeanvalues_data[i][j][k].tab[c][ll][3] << " " << fMeanvalues_size[i][j][k].tab[c][ll][0];
            }
            ofs_d << endl;
            ofs_dzvtx << endl;

            ofs_h << fBinning[i][j][k].tab[c][1][0][0]/**fCepi[c][1][i][j][k]*/ << " " << fBinning[i][j][k].tab[c][1][1][0]/**fCepi[c][1][i][j][k]*/ << " " << fBinning_loose[i][j][k].tab[c][1][0][0]/**fCepi[c][1][i][j][k]*/ << " " << fBinning_severe[i][j][k].tab[c][1][0][0]/**fCepi[c][1][i][j][k]*/ << " " <<
                     fBinning[i][j][k].tab[c][1][0][1] << " " << fBinning[i][j][k].tab[c][1][1][1] << " " << fBinning_loose[i][j][k].tab[c][1][0][1] << " " << fBinning_severe[i][j][k].tab[c][1][0][1] << " " <<
                     fBinning[i][j][k].tab[c][1][0][2] << " " << fBinning[i][j][k].tab[c][1][1][2] << " " << fBinning_loose[i][j][k].tab[c][1][0][2] << " " << fBinning_severe[i][j][k].tab[c][1][0][2] << " " <<
                     fBinning[i][j][k].tab[c][1][0][3]/*-fBinning[i][j][k].tab[c][1][0][0]*(1-fCepi[c][1][i][j][k])*/ << " " << fBinning[i][j][k].tab[c][1][1][3]/*-fBinning[i][j][k].tab[c][1][1][0]*(1-fCepi[c][1][i][j][k])*/ << " " << fBinning_loose[i][j][k].tab[c][1][0][3]/*-fBinning_loose[i][j][k].tab[c][1][0][0]*(1-fCepi[c][1][i][j][k])*/ << " " << fBinning_severe[i][j][k].tab[c][1][0][3]/*-fBinning_severe[i][j][k].tab[c][1][0][0]*(1-fCepi[c][1][i][j][k])*/ << " " <<
                     fBinning[i][j][k].tab[c][0][0][0]/**fCepi[c][0][i][j][k]*/ << " " << fBinning[i][j][k].tab[c][0][1][0]/**fCepi[c][0][i][j][k]*/ << " " << fBinning_loose[i][j][k].tab[c][0][0][0]/**fCepi[c][0][i][j][k]*/ << " " << fBinning_severe[i][j][k].tab[c][0][0][0]/**fCepi[c][0][i][j][k]*/ << " " <<
                     fBinning[i][j][k].tab[c][0][0][1] << " " << fBinning[i][j][k].tab[c][0][1][1] << " " << fBinning_loose[i][j][k].tab[c][0][0][1] << " " << fBinning_severe[i][j][k].tab[c][0][0][1] << " " <<
                     fBinning[i][j][k].tab[c][0][0][2] << " " << fBinning[i][j][k].tab[c][0][1][2] << " " << fBinning_loose[i][j][k].tab[c][0][0][2] << " " << fBinning_severe[i][j][k].tab[c][0][0][2] << " " <<
                     fBinning[i][j][k].tab[c][0][0][3]/*-fBinning[i][j][k].tab[c][0][0][0]*(1-fCepi[c][1][i][j][k])*/ << " " << fBinning[i][j][k].tab[c][0][1][3]/*-fBinning[i][j][k].tab[c][0][1][0]*(1-fCepi[c][1][i][j][k])*/ << " " << fBinning_loose[i][j][k].tab[c][0][0][3]/*-fBinning_loose[i][j][k].tab[c][0][0][0]*(1-fCepi[c][1][i][j][k])*/ << " " << fBinning_severe[i][j][k].tab[c][0][0][3]/*-fBinning_severe[i][j][k].tab[c][0][0][0]*(1-fCepi[c][1][i][j][k])*/ << " " << endl;

            for(int zv=0; zv<4; zv++)
            {
              ofs_hzvtx << fBinning_zvtx[i][j][k][zv].tab[c][1][0][0]/**fCepiVtx[c][1][i][j][k][zv]*/ << " " << fBinning_zvtx[i][j][k][zv].tab[c][1][1][0]/**fCepiVtx[c][1][i][j][k][zv]*/ << " " <<
                           fBinning_zvtx[i][j][k][zv].tab[c][1][0][1] << " " << fBinning_zvtx[i][j][k][zv].tab[c][1][1][1] << " " <<
                           fBinning_zvtx[i][j][k][zv].tab[c][1][0][2] << " " << fBinning_zvtx[i][j][k][zv].tab[c][1][1][2] << " " <<
                           fBinning_zvtx[i][j][k][zv].tab[c][1][0][3]/*-fBinning_zvtx[i][j][k][zv].tab[c][1][0][0]*(1-fCepiVtx[c][1][i][j][k][zv])*/ << " " << fBinning_zvtx[i][j][k][zv].tab[c][1][1][3]/*-fBinning_zvtx[i][j][k][zv].tab[c][1][0][0]*(1-fCepiVtx[c][1][i][j][k][zv])*/ << " " <<
                           fBinning_zvtx[i][j][k][zv].tab[c][0][0][0]/**fCepiVtx[c][0][i][j][k][zv]*/ << " " << fBinning_zvtx[i][j][k][zv].tab[c][0][1][0]/**fCepiVtx[c][0][i][j][k][zv]*/ << " " <<
                           fBinning_zvtx[i][j][k][zv].tab[c][0][0][1] << " " << fBinning_zvtx[i][j][k][zv].tab[c][0][1][1] << " " <<
                           fBinning_zvtx[i][j][k][zv].tab[c][0][0][2] << " " << fBinning_zvtx[i][j][k][zv].tab[c][0][1][2] << " " <<
                           fBinning_zvtx[i][j][k][zv].tab[c][0][0][3]/*-fBinning_zvtx[i][j][k][zv].tab[c][0][0][0]*(1-fCepiVtx[c][0][i][j][k][zv])*/ << " " << fBinning_zvtx[i][j][k][zv].tab[c][0][1][3]/*-fBinning_zvtx[i][j][k][zv].tab[c][0][0][0]*(1-fCepiVtx[c][0][i][j][k][zv])*/ << " " << endl;
            }
            for(int th=0; th<8; th++)
            {
              ofs_hth << fBinning_theta[i][j][k][th].tab[c][1][0][0]/**fCepiTh[c][1][i][j][k][th]*/ << " " << fBinning_theta[i][j][k][th].tab[c][1][1][0]/**fCepiTh[c][1][i][j][k][th]*/ << " " <<
                           fBinning_theta[i][j][k][th].tab[c][1][0][1] << " " << fBinning_theta[i][j][k][th].tab[c][1][1][1] << " " <<
                           fBinning_theta[i][j][k][th].tab[c][1][0][2] << " " << fBinning_theta[i][j][k][th].tab[c][1][1][2] << " " <<
                           fBinning_theta[i][j][k][th].tab[c][1][0][3]/*-fBinning_theta[i][j][k][th].tab[c][1][0][0]*(1-fCepiTh[c][1][i][j][k][th])*/ << " " << fBinning_theta[i][j][k][th].tab[c][1][1][3]/*-fBinning_theta[i][j][k][th].tab[c][1][0][0]*(1-fCepiTh[c][1][i][j][k][th])*/ << " " <<
                           fBinning_theta[i][j][k][th].tab[c][0][0][0]/**fCepiTh[c][0][i][j][k][th]*/ << " " << fBinning_theta[i][j][k][th].tab[c][0][1][0]/**fCepiTh[c][0][i][j][k][th]*/ << " " <<
                           fBinning_theta[i][j][k][th].tab[c][0][0][1] << " " << fBinning_theta[i][j][k][th].tab[c][0][1][1] << " " <<
                           fBinning_theta[i][j][k][th].tab[c][0][0][2] << " " << fBinning_theta[i][j][k][th].tab[c][0][1][2] << " " <<
                           fBinning_theta[i][j][k][th].tab[c][0][0][3]/*-fBinning_theta[i][j][k][th].tab[c][0][0][0]*(1-fCepiTh[c][0][i][j][k][th])*/ << " " << fBinning_theta[i][j][k][th].tab[c][0][1][3]/*-fBinning_theta[i][j][k][th].tab[c][0][0][0]*(1-fCepiTh[c][0][i][j][k][th])*/ << " " << endl;
            }
            for(int pt=0; pt<10; pt++)
            {
              ofs_hpt << fBinning_pt[i][j][k][pt].tab[c][1][0][0]/**fCepipT[c][1][i][j][k][pt]*/ << " " << fBinning_pt[i][j][k][pt].tab[c][1][1][0]/**fCepipT[c][1][i][j][k][pt]*/ << " " <<
                           fBinning_pt[i][j][k][pt].tab[c][1][0][1] << " " << fBinning_pt[i][j][k][pt].tab[c][1][1][1] << " " <<
                           fBinning_pt[i][j][k][pt].tab[c][1][0][2] << " " << fBinning_pt[i][j][k][pt].tab[c][1][1][2] << " " <<
                           fBinning_pt[i][j][k][pt].tab[c][1][0][3]/*-fBinning_pt[i][j][k][pt].tab[c][1][0][0]*(1-fCepipT[c][1][i][j][k][pt])*/ << " " << fBinning_pt[i][j][k][pt].tab[c][1][1][3]/*-fBinning_pt[i][j][k][pt].tab[c][1][0][0]*(1-fCepipT[c][1][i][j][k][pt])*/ << " " <<
                           fBinning_pt[i][j][k][pt].tab[c][0][0][0]/**fCepipT[c][0][i][j][k][pt]*/ << " " << fBinning_pt[i][j][k][pt].tab[c][0][1][0]/**fCepipT[c][0][i][j][k][pt]*/ << " " <<
                           fBinning_pt[i][j][k][pt].tab[c][0][0][1] << " " << fBinning_pt[i][j][k][pt].tab[c][0][1][1] << " " <<
                           fBinning_pt[i][j][k][pt].tab[c][0][0][2] << " " << fBinning_pt[i][j][k][pt].tab[c][0][1][2] << " " <<
                           fBinning_pt[i][j][k][pt].tab[c][0][0][3]/*-fBinning_pt[i][j][k][pt].tab[c][0][0][0]*(1-fCepipT[c][0][i][j][k][pt])*/ << " " << fBinning_pt[i][j][k][pt].tab[c][0][1][3]/*-fBinning_pt[i][j][k][pt].tab[c][0][0][0]*(1-fCepipT[c][0][i][j][k][pt])*/ << " " << endl;
            }

        	  xc << c << " " << fXrange[i] << " " << fYrange[j] << " " << fZrange[k] << " " <<
        		fNDIS_evt[0][0][i][j][k] << " " << fBinning[i][j][k].tab[c][0][0][0] << " " << fBinning[i][j][k].tab[c][0][0][1] << " " <<
        		fBinning[i][j][k].tab[c][0][0][2] << " " << fBinning[i][j][k].tab[c][0][0][3] << endl;

  	      }
        }
      }
    }
    for(int i=0; i<9; i++)
    {
      for(int j=0; j<5; j++)
      {
        for(int k=0; k<12; k++)
        {
          count << fXrange[i] << " " << fYrange[j] << " " << fZrange[k] << " "
               << fBinning[i][j][k].tab[1][1][0][0]+fBinning[i][j][k].tab[1][0][0][0] << " "
               << fBinning[i][j][k].tab[0][1][0][0]+fBinning[i][j][k].tab[0][0][0][0] << " "
               << fNDIS_evt[0][0][i][j][k]+fNDIS_evt[0][1][i][j][k] << " "
               << fBinning[i][j][k].tab[1][1][0][1]+fBinning[i][j][k].tab[1][0][0][1] << " "
               << fBinning[i][j][k].tab[0][1][0][1]+fBinning[i][j][k].tab[0][0][0][1] << " "
               << fNDIS_evt[1][0][i][j][k]+fNDIS_evt[1][1][i][j][k] << " "
               << fBinning[i][j][k].tab[1][1][0][2]+fBinning[i][j][k].tab[1][0][0][2] << " "
               << fBinning[i][j][k].tab[0][1][0][2]+fBinning[i][j][k].tab[0][0][0][2] << " "
               << fNDIS_evt[2][0][i][j][k]+fNDIS_evt[2][1][i][j][k] << " "
               << fBinning[i][j][k].tab[1][1][0][3]+fBinning[i][j][k].tab[1][0][0][3] << " "
               << fBinning[i][j][k].tab[0][1][0][3]+fBinning[i][j][k].tab[0][0][0][3] << " "
               << fNDIS_evt[0][0][i][j][k]+fNDIS_evt[0][1][i][j][k] << endl;
        }
      }
    }

    ofs_h.close();
    ofs_d.close();
    ofs_hzvtx.close();
    ofs_dzvtx.close();
    xc.close();

    resetValues();
  }

  if(kin_flag)
  {
    for(int i=0; i<int(fQ2.size()); i++)
    {
      fKinematics[0]->Fill(fQ2kin[i]);
      fKinematics[1]->Fill(fXBjkin[i]);
      fKinematics[2]->Fill(fYBjkin[i]);
      fKinematics[4]->Fill(fWBjkin[i]);
      fKinematics[5]->Fill(fNukin[i]);
      fKinematics2D[0]->Fill(fXBjkin[i],fYBjkin[i]);
      fKinematics2D[1]->Fill(fXBjkin[i],fQ2kin[i]);
      fTarget2D->Fill(fX[i],fY[i]);
      fHO03->Fill(fHO03x[i],fHO03y[i]);
      fHO04->Fill(fHO04x[i],fHO04y[i]);
    }
    for(int i=0; i<int(fLHpi.size()); i++)
    {
      fRICHLH->Fill(fLHK[i],fLHpi[i]);
    }
    save_kin_plots();
  }

  cout << "\n\n";
  cout << "             ********* Cut flow for Reconstructed DIS events after cuts ********* " << endl;
  cout << "             -------------------------------------------------------------------- " << endl;
  cout << "             from " << fFilesNumber << " file(s)" << endl;

  for(int j=0; j<2; j++)
  {
    for(int i=0; i<int(fPeriodBit.size()); i++)
    {
      if(!j)
        cout << setw(6) << fPeriodName[i];
      else
        cout << setw(6) << fPeriodBit[i];
    }
    cout << endl;
  }

  cout.precision(16);

  cout << '|' << setw(30) << "Cut" << '|' << setw(15) << "Events" << '|' << setw(15) << "Abs." << '|' << setw(15) << "Rel." << '|' << endl;
  cout << '|' << setw(30) << "Best Primary Vertex" << '|' << setw(15) << fBP << '|' << setw(15) << float(fBP)/float(fBP)*100 << '|' << setw(15) << float(fBP)/float(fBP)*100 << '|' << endl;
  cout << '|' << setw(30) << "Mu' found (0,1,1,1,30)" << '|' << setw(15) << fRmu << '|' << setw(15) << float(fRmu)/float(fBP)*100 << '|' << setw(15) << float(fRmu)/float(fBP)*100 << '|' << endl;
  cout << '|' << setw(30) << "Vertex in Target" << '|' << setw(15) << fTarg << '|' << setw(15) << float(fTarg)/float(fBP)*100 << '|' << setw(15) << float(fTarg)/float(fRmu)*100 << '|' << endl;
  cout << '|' << setw(30) << "140 < E_mu < 180" << '|' << setw(15) << fBEC << '|' << setw(15) << float(fBEC)/float(fBP)*100 << '|' << setw(15) << float(fBEC)/float(fTarg)*100 << '|' << endl;
  cout << '|' << setw(30) << "BMS" << '|' << setw(15) << fBMS << '|' << setw(15) << float(fBMS)/float(fBP)*100 << '|' << setw(15) << float(fBMS)/float(fBEC)*100 << '|' << endl;
  cout << '|' << setw(30) << "Mu chi2/ndf < 10" << '|' << setw(15) << fMuchi2 << '|' << setw(15) << float(fMuchi2)/float(fBP)*100 << '|' << setw(15) << float(fMuchi2)/float(fBMS)*100 << '|' << endl;
  cout << '|' << setw(30) << "Beam track X Cell" << '|' << setw(15) << fCell << '|' << setw(15) << float(fCell)/float(fBP)*100 << '|' << setw(15) << float(fCell)/float(fMuchi2)*100 << '|' << endl;
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
  cout << "             from " << fFilesNumber << " file(s)" << endl;

  for(int j=0; j<2; j++)
  {
    for(int i=0; i<int(fPeriodBit.size()); i++)
    {
      if(!j)
        cout << setw(6) << fPeriodName[i];
      else
        cout << setw(6) << fPeriodBit[i];
    }
    cout << endl;
  }

  cout << '|' << setw(30) << "Cut" << '|' << setw(15) << "Events" << '|' << setw(15) << "Abs." << '|' << setw(15) << "Rel." << endl;
  cout << '|' << setw(30) << "Hadrons" << '|' << setw(15) << fHadrons << '|' << setw(15) << float(fHadrons)/float(fHadrons)*100 << '|' << setw(15) << float(fHadrons)/float(fHadrons)*100 << endl;
  cout << '|' << setw(30) << "XX0 < 15" << '|' << setw(15) << fXX0test << '|' << setw(15) << float(fXX0test)/float(fHadrons)*100 << '|' << setw(15) << float(fXX0test)/float(fHadrons)*100 << endl;
  cout << '|' << setw(30) << "Chi2/ndf < 10" << '|' << setw(15) << fChi2Hadron << '|' << setw(15) << float(fChi2Hadron)/float(fHadrons)*100 << '|' << setw(15) << float(fChi2Hadron)/float(fXX0test)*100 << endl;
  cout << '|' << setw(30) << "Zfirst < 350 cm" << '|' << setw(15) << fHZfirst << '|' << setw(15) << float(fHZfirst)/float(fHadrons)*100 << '|' << setw(15) << float(fHZfirst)/float(fChi2Hadron)*100 << endl;
  cout << '|' << setw(30) << "Zlast > 350 cm" << '|' << setw(15) << fHZlast << '|' << setw(15) << float(fHZlast)/float(fHadrons)*100 << '|' << setw(15) << float(fHZlast)/float(fHZfirst)*100 << endl;
  cout << '|' << setw(30) << "0.01 < theta_RICH < 0.12" << '|' << setw(15) << fTRICH << '|' << setw(15) << float(fTRICH)/float(fHadrons)*100 << '|' << setw(15) << float(fTRICH)/float(fHZlast)*100 << endl;
  cout << '|' << setw(30) << "Rich Pipe" << '|' << setw(15) << fPosRICH << '|' << setw(15) << float(fPosRICH)/float(fHadrons)*100 << '|' << setw(15) << float(fPosRICH)/float(fTRICH)*100 << endl;
  cout << '|' << setw(30) << "12 < p_h < 40" << '|' << setw(15) << fMom << '|' << setw(15) << float(fMom)/float(fHadrons)*100 << '|' << setw(15) << float(fMom)/float(fPosRICH)*100 << endl;
  cout << '|' << setw(30) << "0.2 < z < 0.85" << '|' << setw(15) << fZtest << '|' << setw(15) << float(fZtest)/float(fHadrons)*100 << '|' << setw(15) << float(fZtest)/float(fMom)*100 << endl;

  cout << "\n\n";
  cout << "             ********* Hadron Content (h,pi,K,p) ********* " << endl;
  cout << "             --------------------------------------------- " << endl;

  for(int j=0; j<2; j++)
  {
    for(int i=0; i<int(fPeriodBit.size()); i++)
    {
      if(!j)
        cout << setw(6) << fPeriodName[i];
      else
        cout << setw(6) << fPeriodBit[i];
    }
    cout << endl;
  }

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


  // test.close();
  count.close();

  return 0;
}

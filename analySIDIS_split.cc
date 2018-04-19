#include "analySIDIS_split.h"

using namespace std;

//Inputs
#define mat_RICH_name "data/rich_mat.txt"
#define err_RICH_name "data/rich_mat_error.txt"
#define target_file_2012 "data/target-107924-109081.dat"
#define target_file_2016 "data/target-274508-274901.dat"
#define proton_sirc "data/proton_semi_inclusive_RC.txt"
#define proton_irc "data/hh160_r1998_f2tulay_compass_grv.asy_hcorr.txt"

// Flags
#define Y2006 0
#define Y2012 0
#define Y2016 1
#define RCUTSTUDY_ON 0
#define MOMENTUM 3

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
    float z, x, y, dummy;
    fin >> z >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> x >> y;
    fZv.push_back(z);
    fXv.push_back(x);
    fYv.push_back(y);
  }
  cout<<"INFO : Target cell description loaded"<<endl;
}

void CellCenter(Double_t z, Double_t& xc, Double_t& yc)
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

    Double_t dxcdz = (xc2-xc1)/(z2-z1);
    Double_t dycdz = (yc2-yc1)/(z2-z1);

    Double_t dz = z-z1;
    xc = xc1 + dxcdz*dz;
    yc = yc1 + dycdz*dz;

    break;
  }
}

bool InTarget(Double_t xvtx, Double_t yvtx, Double_t zvtx, Double_t R)
{
  Double_t xc, yc;
  CellCenter(zvtx, xc, yc);
  Double_t dx = xvtx-xc;
  Double_t dy = yvtx-yc;
  Double_t r = sqrt(dx*dx + dy*dy);

  return( r <= R );
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
        proton >> fInclusiveRCproton[i][k+j*6] >> sdum;
#ifdef DEBUG
        cout << " " << fInclusiveRCproton[i][k+j*6] << sdum;
#endif
      }

#ifdef DEBUG
      cout << endl;
#endif
    }
  }
  proton.close();
}

void LoadSemiInclusiveRadiativeCorrection()
{
  string sdum;

  ifstream proton(proton_sirc);
  for(int i=0; i<9; i++)
  {
    proton >> sdum;
#ifdef DEBUG
    cout << " " << sdum;
#endif
  }
  for(int i=0; i<9; i++)
  {
    for(int j=0; j<6; j++)
    {
      proton >> sdum;
#ifdef DEBUG
      cout << sdum << "\t";
#endif
      proton >> sdum;
#ifdef DEBUG
      cout << sdum << "\t";
#endif

      for(int k=0; k<14; k++)
      {
        proton >> fSemiInclusiveRCproton[i][j][k];
#ifdef DEBUG
        cout << " " << fSemiInclusiveRCproton[i][j][k];
#endif
      }

#ifdef DEBUG
      cout << endl;
#endif
    }
  }
  proton.close();
}

Double_t GetInclusiveRadiativeCorrection(Double_t x, Double_t y)
{
  int xb, yb;

  if(0.00005<x && x<0.00007) xb = 0;
  else if(0.00007<=x && x<0.0001) xb = 1;
  else if(0.0001<=x && x<0.0002) xb = 2;
  else if(0.0002<=x && x<0.0003) xb = 3;
  else if(0.0003<=x && x<0.0005) xb = 4;
  else if(0.0005<=x && x<0.0007) xb = 5;
  else if(0.0007<=x && x<0.001) xb = 6;
  else if(0.001<=x && x<0.002) xb = 7;
  else if(0.002<=x && x<0.004) xb = 8;
  else if(0.004<=x && x<0.006) xb = 9;
  else if(0.006<=x && x<0.008) xb = 10;
  else if(0.008<=x && x<0.01) xb = 11;
  else if(0.01<=x && x<0.013) xb = 12;
  else if(0.013<=x && x<0.016) xb = 13;
  else if(0.016<=x && x<0.02) xb = 14;
  else if(0.02<=x && x<0.03) xb = 15;
  else if(0.03<=x && x<0.04) xb = 16;
  else if(0.04<=x && x<0.06) xb = 17;
  else if(0.06<=x && x<0.08) xb = 18;
  else if(0.08<=x && x<0.1) xb = 19;
  else if(0.1<=x && x<0.15) xb = 20;
  else if(0.15<=x && x<0.2) xb = 21;
  else if(0.2<=x && x<0.3) xb = 22;
  else if(0.3<=x && x<0.4) xb = 23;
  else if(0.4<=x && x<0.5) xb = 24;
  else if(0.5<=x && x<0.6) xb = 25;
  else if(0.6<=x && x<0.7) xb = 26;
  else if(0.7<=x && x<0.8) xb = 27;
  else xb = 28;

  if(0.05<y && y<0.1) yb = 0;
  else if(0.1<=y && y<0.15) yb = 1;
  else if(0.15<=y && y<0.2) yb = 2;
  else if(0.2<=y && y<0.25) yb = 3;
  else if(0.25<=y && y<0.3) yb = 4;
  else if(0.3<=y && y<0.35) yb = 5;
  else if(0.3<=y && y<0.35) yb = 6;
  else if(0.35<=y && y<0.4) yb = 7;
  else if(0.4<=y && y<0.45) yb = 8;
  else if(0.45<=y && y<0.5) yb = 9;
  else if(0.5<=y && y<0.55) yb = 10;
  else if(0.55<=y && y<0.6) yb = 11;
  else if(0.6<=y && y<0.65) yb = 12;
  else if(0.65<=y && y<0.7) yb = 13;
  else if(0.7<=y && y<0.75) yb = 14;
  else if(0.75<=y && y<0.8) yb = 15;
  else if(0.8<=y && y<0.85) yb = 16;
  else yb = 17;

  if(2006)
  {
    return 1;
  }
  else if(2012 || 2016)
  {
    return fInclusiveRCproton[xb][yb];
  }
  else
  {
    cout << "ERROR in GetInclusiveRadiativeCorrection : Year not recognized. No correction applied." << endl;
    return 1;
  }
}

Double_t GetSemiInclusiveRadiativeCorrection(Double_t x, Double_t y, Double_t z)
{
  int xb, yb, zb;

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

  if(2006)
  {
    return 1;
  }
  else if(2012 || 2016)
  {
    return fSemiInclusiveRCproton[xb][yb][zb];
  }
  else
  {
    cout << "ERROR in GetSemiInclusiveRadiativeCorrection : Year not recognized. No correction applied." << endl;
    return 1;
  }
}

void load_rich_mat(string prich, string prich_err)
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
    inv_rich_p[0][i] = rich_mat_p[0][i].InvertFast();
    inv_rich_p[1][i] = rich_mat_p[1][i].InvertFast();
    inv_rich_m[0][i] = rich_mat_m[0][i].InvertFast();
    inv_rich_m[1][i] = rich_mat_m[1][i].InvertFast();
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

void create_kin_plots()
{
  fKinematics[0] = new TH1F("Q^{2}", "Q^{2}", 100, 0, 2);
  fKinematics[1] = new TH1F("x_{Bj}", "x_{Bj}", 100, -3, 0);
  fKinematics[2] = new TH1F("y", "y", 100, 0, 1);
  fKinematics[3] = new TH1F("z", "z", 100, 0, 1);
  fKinematics[4] = new TH1F("W", "W", 100, 2, 18);
  fKinematics[5] = new TH1F("#nu", "#nu", 100, 0, 160);
  fKinematics2D = new TH2F("DIS kin space", "DIS kin space", 100, -3, 0, 100, 0.1, 0.7);
  fTarget2D = new TH2F("Target xy", "Target xy", 100, -3, 3, 100, -3, 3);
  fHO03 = new TH2F("HO03", "HO03", 100, -120, 120, 100, -60, 60);
  fHO04 = new TH2F("HO04", "HO04", 100, -250, 250, 100, -100, 100);
  fRICHLH = new TH2F("RICH LH", "RICH LH", 100, -2, 2, 100, -2, 2);
  BinLogX(fKinematics[0]);
  BinLogX(fKinematics[1]);
  BinLogX(fKinematics2D);
}

void save_kin_plots()
{
  c1.Divide(1,1);
  c2.Divide(1,1);
  c3.Divide(1,1);
  c4.Divide(1,1);
  c5.Divide(1,1);
  c6.Divide(1,1);
  c7.Divide(1,1);
  c8.Divide(1,1);
  c9.Divide(1,1);
  c10.Divide(1,1);
  c11.Divide(1,1);
  c1.cd(1);
  fKinematics[0]->Draw();
  gPad->SetLogx();
  c1.Update();
  c2.cd(1);
  fKinematics[1]->Draw();
  gPad->SetLogx();
  c2.Update();
  c3.cd(1);
  fKinematics[2]->Draw();
  c3.Update();
  c4.cd(1);
  fKinematics[3]->Draw();
  c4.Update();
  c5.cd(1);
  fKinematics[4]->Draw();
  c5.Update();
  c6.cd(1);
  fKinematics[5]->Draw();
  c6.Update();
  c7.cd(1);
  fKinematics2D->Draw("COLZ");
  gPad->SetLogx();
  c7.Update();
  c8.cd(1);
  fTarget2D->Draw("COLZ");
  c8.Update();
  c9.cd(1);
  fRICHLH->Draw();
  c9.Update();
  c10.cd(1);
  fHO03->Draw("COLZ");
  c10.Update();
  c11.cd(1);
  fHO04->Draw("COLZ");
  c11.Update();

  c1.Print("kinSIDIS.pdf(","pdf");
  c2.Print("kinSIDIS.pdf","pdf");
  c3.Print("kinSIDIS.pdf","pdf");
  c4.Print("kinSIDIS.pdf","pdf");
  c5.Print("kinSIDIS.pdf","pdf");
  c6.Print("kinSIDIS.pdf","pdf");
  c7.Print("kinSIDIS.pdf","pdf");
  c8.Print("kinSIDIS.pdf","pdf");
  c9.Print("kinSIDIS.pdf","pdf");
  c10.Print("kinSIDIS.pdf","pdf");
  c11.Print("kinSIDIS.pdf)","pdf");
}

int main(int argc, char **argv)
{

  if(argc < 2)
  {
    cout << "ERROR : Not enough arguments." << endl;
    cout << "Asked : at least 1 *** Received : " << argc-1 << endl;
    cout << "./analySIDIS_split filelist" << endl;

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

  //Kinematics
  Double_t Q2 = 0;
  Double_t xBj = 0;
  Double_t yBj = 0;
  Double_t zBj = 0;
  Double_t wBj = 0;
  Double_t nu = 0;

  if(kin_flag) create_kin_plots();

  load_rich_mat(mat_RICH_name, err_RICH_name);

  //cout << pi_sigma_uni[0][0] << " " << pi_sigma_uni[1][1] << " " << pi_sigma_uni[2][2] << endl;

  // Target cells
  if(Y2012) InitTargetFile(target_file_2012);
  else if(Y2016) InitTargetFile(target_file_2016);

  //----------------------------------------------------------------------------
  //--------- nu cut prep ------------------------------------------------------
  //----------------------------------------------------------------------------

  for(int i=0; i<12; i++)
  {
    fNu_max[1][i] = sqrt(pow(40,2)+pow(fM_K,2))/fZrange[i+1];
    fNu_min[1][i] = sqrt(pow(MOMENTUM,2)+pow(fM_K,2))/fZrange[i];

    fNu_max[2][i] = sqrt(pow(40,2)+pow(fM_p,2))/fZrange[i+1];
    fNu_min[2][i] = sqrt(pow(MOMENTUM,2)+pow(fM_p,2))/fZrange[i];

    fNu_max[0][i] = sqrt(pow(40,2)+pow(fM_pi,2))/fZrange[i+1];
    fNu_min[0][i] = sqrt(pow(MOMENTUM,2)+pow(fM_pi,2))/fZrange[i];
  }

  // List of files

  ifstream list(argv[1]);
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
    TBranch *cellsCrossed = (TBranch*) tree->FindBranch("cellsCrossed");
    TBranch *backPropFlag = (TBranch*) tree->FindBranch("backPropFlag");

    //Hadrons
    TBranch *p = (TBranch*) tree->FindBranch("Hadrons.P");
    TBranch *th = (TBranch*) tree->FindBranch("Hadrons.th");
    TBranch *ph = (TBranch*) tree->FindBranch("Hadrons.ph");
    TBranch *hXX0 = (TBranch*) tree->FindBranch("Hadrons.XX0");
    TBranch *inHCALacc = (TBranch*) tree->FindBranch("Hadrons.inHCALacc");
    TBranch *HCAL = (TBranch*) tree->FindBranch("Hadrons.HCAL");
    TBranch *charge = (TBranch*) tree->FindBranch("Hadrons.charge");
    TBranch *thRICH = (TBranch*) tree->FindBranch("Hadrons.thRICH");
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

    // Loopy loop over the events
    Int_t N = (Int_t) tree->GetEntries();

    vector<Pvsz> Pvszlocal;
    vector<Pvsz> Pvszloose;
    vector<Pvsz> Pvszsevere;
    vector<Pvsz> Pvsz_errlocal;
    vector<Double_t> XBjlocal;
    vector<Double_t> YBjlocal;
    vector<Double_t> Q2local;
    vector<Double_t> XBjloose;
    vector<Double_t> YBjloose;
    vector<Double_t> Q2loose;
    vector<Double_t> XBjsevere;
    vector<Double_t> YBjsevere;
    vector<Double_t> Q2severe;

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
      cellsCrossed->GetEntry(ip);
      backPropFlag->GetEntry(ip);

      //Hadrons
      p->GetEntry(ip);
      th->GetEntry(ip);
      ph->GetEntry(ip);
      hXX0->GetEntry(ip);
      inHCALacc->GetEntry(ip);
      HCAL->GetEntry(ip);
      charge->GetEntry(ip);
      thRICH->GetEntry(ip);
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

      // -------------------------------------------------------------------------
      // --------- Calculation ---------------------------------------------------
      // -------------------------------------------------------------------------

      Double_t zlab = z->GetLeaf("z")->GetValue();
      Int_t zlabbin;

      if(Y2012 || RCUTSTUDY_ON)
      {
        if(!(-311.2<zlab && zlab<-71.2)) continue;

        if(-311.2<=zlab && zlab<-301.2) zlabbin = 0;
        else if(-301.2<=zlab && zlab<-291.2) zlabbin = 1;
        else if(-291.2<=zlab && zlab<-281.2) zlabbin = 2;
        else if(-281.2<=zlab && zlab<-271.2) zlabbin = 3;
        else if(-271.2<=zlab && zlab<-261.2) zlabbin = 4;
        else if(-261.2<=zlab && zlab<-251.2) zlabbin = 5;
        else if(-251.2<=zlab && zlab<-241.2) zlabbin = 6;
        else if(-241.2<=zlab && zlab<-231.2) zlabbin = 7;
        else if(-231.2<=zlab && zlab<-221.2) zlabbin = 8;
        else if(-221.2<=zlab && zlab<-211.2) zlabbin = 9;
        else if(-211.2<=zlab && zlab<-201.2) zlabbin = 10;
        else if(-201.2<=zlab && zlab<-191.2) zlabbin = 11;
        else if(-191.2<=zlab && zlab<-181.2) zlabbin = 12;
        else if(-181.2<=zlab && zlab<-171.2) zlabbin = 13;
        else if(-171.2<=zlab && zlab<-161.2) zlabbin = 14;
        else if(-161.2<=zlab && zlab<-151.2) zlabbin = 15;
        else if(-151.2<=zlab && zlab<-141.2) zlabbin = 16;
        else if(-141.2<=zlab && zlab<-131.2) zlabbin = 17;
        else if(-131.2<=zlab && zlab<-121.2) zlabbin = 18;
        else if(-121.2<=zlab && zlab<-111.2) zlabbin = 19;
        else if(-111.2<=zlab && zlab<-101.2) zlabbin = 20;
        else if(-101.2<=zlab && zlab<-91.2) zlabbin = 21;
        else if(-91.2<=zlab && zlab<-81.2) zlabbin = 22;
        else zlabbin = 23;
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

      // Best Primary Vertex
      fBP++;

      // Reconstructed muon
      if(!(0<E_beam->GetLeaf("E_beam")->GetValue())) continue;
      fRmu++;

      Double_t mxc, myc;

      if(Y2012 || RCUTSTUDY_ON)
      {
        CellCenter(z->GetLeaf("z")->GetValue(), mxc, myc);
      }

      //Rcut study ---
      if(RCUTSTUDY_ON)
      {
        fRstudy[zlabbin].vec.push_back(sqrt((pow(x->GetLeaf("x")->GetValue()-mxc,2)+pow(y->GetLeaf("y")->GetValue()-myc,2))));
        fRstudy_xy[zlabbin].vec[0].push_back(x->GetLeaf("x")->GetValue());
        fRstudy_xy[zlabbin].vec[1].push_back(y->GetLeaf("y")->GetValue());

        if(!InTarget(x->GetLeaf("x")->GetValue(),y->GetLeaf("y")->GetValue(),z->GetLeaf("z")->GetValue(),fRcutval[zlabbin])) continue;
        fTarg++;

        fR_xy[zlabbin].vec[0].push_back(x->GetLeaf("x")->GetValue());
        fR_xy[zlabbin].vec[1].push_back(y->GetLeaf("y")->GetValue());
      }
      //Rcut study ---


      //BMS (reconstructed beam track)
      if((backPropFlag->GetLeaf("backPropFlag")->GetValue())) continue;
      fBMS++;

      // Energy of the muon beam
      if(!(140<E_beam->GetLeaf("E_beam")->GetValue() && E_beam->GetLeaf("E_beam")->GetValue()<180)) continue;
      fBEC++;

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
        if(!InTarget(x->GetLeaf("x")->GetValue(),y->GetLeaf("y")->GetValue(),z->GetLeaf("z")->GetValue(),fRcutval[zlabbin])) continue;
      }
      //2012 ---
      //2016 ---
      else if(Y2016)
      {
        if(!InTarget(x->GetLeaf("x")->GetValue(),y->GetLeaf("y")->GetValue(),z->GetLeaf("z")->GetValue(),1.5)) continue;
      }
      //2016 ---
      fTarg++;

      // Cells crossing
      //if(!(cellsCrossed->GetLeaf("cellsCrossed")->GetValue())) continue;
      //fCell++;

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

      // Q2 cut
      if(!(Q2>1)) continue;
      fQ2test++;

      // y cut
      if(!(0.1<yBj && yBj<0.9)) continue;
      fYBjtest++;

      // W cut
      if(!(5<sqrt(wBj) && sqrt(wBj)<17)) continue;
      fWBjtest++;

      // x cut
      if(!(0.004<xBj && xBj<0.4)) continue;
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
      else ybin = 4;


      // z binning

      for(int i=0; i<12; i++)
      {
        fNDIS_evt[0][xbin][ybin][i] += 1*GetInclusiveRadiativeCorrection(xBj,yBj);
        fNDIS_evt[1][xbin][ybin][i] += 1*GetInclusiveRadiativeCorrection(xBj,yBj);
        fNDIS_evt[2][xbin][ybin][i] += 1*GetInclusiveRadiativeCorrection(xBj,yBj);

        fNDIS_evt_err[0][xbin][ybin][i] += pow(GetInclusiveRadiativeCorrection(xBj,yBj),2);
        fNDIS_evt_err[1][xbin][ybin][i] += pow(GetInclusiveRadiativeCorrection(xBj,yBj),2);
        fNDIS_evt_err[2][xbin][ybin][i] += pow(GetInclusiveRadiativeCorrection(xBj,yBj),2);

        fFlag[0][xbin][ybin][i]=0;
        fFlag[1][xbin][ybin][i]=0;
        fFlag[2][xbin][ybin][i]=0;

        // nu cut
        if(!(fNu_min[0][i]<nu && nu<fNu_max[0][i]))
        {
          fFlag[0][xbin][ybin][i]=1;
        }
        if(!(fNu_min[1][i]<nu && nu<fNu_max[1][i]))
        {
          fFlag[1][xbin][ybin][i]=1;
        }
        if(!(fNu_min[2][i]<nu && nu<fNu_max[2][i]))
        {
          fFlag[2][xbin][ybin][i]=1;
        }
        if(fFlag[0][xbin][ybin][i] /*|| fFlag[1][xbin][ybin][i] || fFlag[2][xbin][ybin][i]*/)
        {
          fNDIS_evt[0][xbin][ybin][i] -= GetInclusiveRadiativeCorrection(xBj,yBj);
         // fNDIS_evt[1][xbin][ybin][i] -= GetInclusiveRadiativeCorrection(xBj,yBj);
         // fNDIS_evt[2][xbin][ybin][i] -= GetInclusiveRadiativeCorrection(xBj,yBj);
          fNDIS_evt_err[0][xbin][ybin][i] -= pow(GetInclusiveRadiativeCorrection(xBj,yBj),2);
         // fNDIS_evt_err[1][xbin][ybin][i] -= GetInclusiveRadiativeCorrection(xBj,yBj),2);
         // fNDIS_evt_err[2][xbin][ybin][i] -= GetInclusiveRadiativeCorrection(xBj,yBj),2);
        }
      }

      // -------------------------------------------------------------------------
      // --------- Hadrons Selection ---------------------------------------------
      // -------------------------------------------------------------------------

      Pvsz pzcontainer;
      Pvsz pzcontainer_loose;
      Pvsz pzcontainer_severe;
      Pvsz pzcontainer_err;
      hadiden hadcontainer;

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
          else fId = 6;

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
          else fId = 7;

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
        if(!(hXX0->GetLeaf("Hadrons.XX0")->GetValue(i) < 15)) continue;
        fXX0test++;

        // Momentum cut (12 GeV to 40 GeV, increasing to 3 GeV to 40 GeV)
        if(!(MOMENTUM<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<40)) continue;
        fMom++;

        // Theta cut
        if(!(0.01<thRICH->GetLeaf("Hadrons.thRICH")->GetValue(i) && thRICH->GetLeaf("Hadrons.thRICH")->GetValue(i)<0.12)) continue;
        fTRICH++;

        // RICH position cut
        if(!(pow(RICHx->GetLeaf("Hadrons.RICHx")->GetValue(i),2)+pow(RICHy->GetLeaf("Hadrons.RICHy")->GetValue(i),2)>25)) continue;
        fPosRICH++;

        // Non null charge
        if(!charge->GetLeaf("Hadrons.charge")->GetValue(i)) continue;

        Int_t theta_bin, mom_bin;
        TMatrixD res_vect(3,1);
        Double_t res_vect_err[3];
        Double_t hadron_nb;

        // Theta and momentum binning

        if(0.01<thRICH->GetLeaf("Hadrons.thRICH")->GetValue(i) && thRICH->GetLeaf("Hadrons.thRICH")->GetValue(i)<0.04)
        {
          theta_bin = 0;
          if(MOMENTUM<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<13) mom_bin = 0; // Here from 3 to 13 GeV
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
          if(MOMENTUM<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<13) mom_bin = 0; // Here from 3 to 13 GeV
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

        // z cut
        if(!(0.2<zBj && zBj<0.85)) continue;

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

            res_vect[0][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
            res_vect[1][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
            res_vect[2][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
            res_vect_err[0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
            res_vect_err[1] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
            res_vect_err[2] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
            hadron_nb *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
            fPiplus_true += res_vect[0][0]; fKplus_true += res_vect[1][0]; fPplus_true += res_vect[2][0];
            fPiplus_err += pow(res_vect_err[0],2);
            fKplus_err += pow(res_vect_err[1],2);
            fPplus_err += pow(res_vect_err[2],2);

            pzcontainer.vec[1][0].push_back(zBj);
            pzcontainer.vec[1][1].push_back(res_vect[0][0]);
            pzcontainer.vec[1][2].push_back(res_vect[1][0]);
            pzcontainer.vec[1][3].push_back(res_vect[2][0]);

            pzcontainer_err.vec[1][0].push_back(zBj);
            pzcontainer_err.vec[1][1].push_back(pow(res_vect_err[0],2));
            pzcontainer_err.vec[1][2].push_back(pow(res_vect_err[1],2));
            pzcontainer_err.vec[1][3].push_back(pow(res_vect_err[2],2));

            pzcontainer.vec[1][4].push_back(hadron_nb);
            pzcontainer_err.vec[1][4].push_back(hadron_nb);

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
            res_vect[0][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
            res_vect[1][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
            res_vect[2][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
            res_vect_err[0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
            res_vect_err[1] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
            res_vect_err[2] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
            hadron_nb *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
            fPiminus_true += res_vect[0][0]; fKminus_true += res_vect[1][0]; fPminus_true += res_vect[2][0];
            fPiminus_err += pow(res_vect_err[0],2);
            fKminus_err += pow(res_vect_err[1],2);
            fPminus_err += pow(res_vect_err[2],2);

            pzcontainer.vec[0][0].push_back(zBj);
            pzcontainer.vec[0][1].push_back(res_vect[0][0]);
            pzcontainer.vec[0][2].push_back(res_vect[1][0]);
            pzcontainer.vec[0][3].push_back(res_vect[2][0]);

            pzcontainer_err.vec[0][0].push_back(zBj);
            pzcontainer_err.vec[0][1].push_back(pow(res_vect_err[0],2));
            pzcontainer_err.vec[0][2].push_back(pow(res_vect_err[1],2));
            pzcontainer_err.vec[0][3].push_back(pow(res_vect_err[2],2));

            pzcontainer.vec[0][4].push_back(hadron_nb);
            pzcontainer_err.vec[0][4].push_back(hadron_nb);

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
            hadron_flag = 1;
          }
          if(!fFlag[1][xbin][ybin][zbin])
          {
            fKplus++;
            res_vect = inv_rich_p[theta_bin][mom_bin]*k_vect;
            for(int rce=0; rce<3; rce++) res_vect_err[rce] = k_unfolding_err_p[theta_bin][mom_bin][rce];
            hadron_nb = 1;
            res_vect[0][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
            res_vect[1][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
            res_vect[2][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
            res_vect_err[0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
            res_vect_err[1] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
            res_vect_err[2] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
            hadron_nb *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
            fPiplus_true += res_vect[0][0]; fKplus_true += res_vect[1][0]; fPplus_true += res_vect[2][0];
            fPiplus_err += pow(res_vect_err[0],2);
            fKplus_err += pow(res_vect_err[1],2);
            fPplus_err += pow(res_vect_err[2],2);

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
            hadron_flag = 1;
          }
          if(!fFlag[1][xbin][ybin][zbin])
          {
            fKminus++;
            res_vect = inv_rich_m[theta_bin][mom_bin]*k_vect;
            for(int rce=0; rce<3; rce++) res_vect_err[rce] = k_unfolding_err_m[theta_bin][mom_bin][rce];
            hadron_nb = 1;
            res_vect[0][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
            res_vect[1][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
            res_vect[2][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
            res_vect_err[0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
            res_vect_err[1] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
            res_vect_err[2] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
            hadron_nb *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
            fPiminus_true += res_vect[0][0]; fKminus_true += res_vect[1][0]; fPminus_true += res_vect[2][0];
            fPiminus_err += pow(res_vect_err[0],2);
            fKminus_err += pow(res_vect_err[1],2);
            fPminus_err += pow(res_vect_err[2],2);

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
            hadron_flag = 1;
          }
          if(!fFlag[2][xbin][ybin][zbin])
          {
            fPplus++;
            res_vect = inv_rich_p[theta_bin][mom_bin]*p_vect;
            for(int rce=0; rce<3; rce++) res_vect_err[rce] = p_unfolding_err_p[theta_bin][mom_bin][rce];
            hadron_nb = 1;
            res_vect[0][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
            res_vect[1][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
            res_vect[2][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
            res_vect_err[0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
            res_vect_err[1] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
            res_vect_err[2] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
            hadron_nb *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
            fPiplus_true += res_vect[0][0]; fKplus_true += res_vect[1][0]; fPplus_true += res_vect[2][0];
            fPiplus_err += pow(res_vect_err[0],2);
            fKplus_err += pow(res_vect_err[1],2);
            fPplus_err += pow(res_vect_err[2],2);

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
            hadron_flag = 1;
          }
          if(!fFlag[2][xbin][ybin][zbin])
          {
            fPminus++;
            res_vect = inv_rich_m[theta_bin][mom_bin]*p_vect;
            for(int rce=0; rce<3; rce++) res_vect_err[rce] = p_unfolding_err_m[theta_bin][mom_bin][rce];
            hadron_nb = 1;
            res_vect[0][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
            res_vect[1][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
            res_vect[2][0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
            res_vect_err[0] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
            res_vect_err[1] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
            res_vect_err[2] *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
            hadron_nb *= GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj);
            fPiminus_true += res_vect[0][0]; fKminus_true += res_vect[1][0]; fPminus_true += res_vect[2][0];
            fPiminus_err += pow(res_vect_err[0],2);
            fKminus_err += pow(res_vect_err[1],2);
            fPminus_err += pow(res_vect_err[2],2);

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
            pzcontainer_err.vec[0][4].push_back(pow(GetSemiInclusiveRadiativeCorrection(xBj,yBj,zBj),2));
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

      Q2local.push_back(Q2);
      Pvszlocal.push_back(pzcontainer);
      Pvsz_errlocal.push_back(pzcontainer_err);
      XBjlocal.push_back(xBj);
      YBjlocal.push_back(yBj);

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
      else ybin = 4;

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
            fBinning_loose[xbin][ybin][zbin].tab[j][0][ll] += Pvszloose[i].vec[j][ll+1][l];
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
      else ybin = 4;

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
            fBinning_severe[xbin][ybin][zbin].tab[j][0][ll] += Pvszsevere[i].vec[j][ll+1][l];
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
      else ybin = 4;

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

          if(kin_flag) fKinematics[3]->Fill(Pvszlocal[i].vec[j][0][l]);

          for(int ll=0; ll<4; ll++)
          {
            fBinning[xbin][ybin][zbin].tab[j][0][ll] += Pvszlocal[i].vec[j][ll+1][l];
            //fBinning[xbin][ybin][zbin].tab[j][1][ll] += Pvsz_errlocal[i].vec[j][ll+1][l];
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
    Pvszloose.clear();
    XBjloose.clear();
    YBjloose.clear();
    Q2loose.clear();
    Pvszsevere.clear();
    XBjsevere.clear();
    YBjsevere.clear();
    Q2severe.clear();
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
      fKinematics2D->Fill(fXBjkin[i],fYBjkin[i]);
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

  for(int i=0; i<9; i++)
  {
    for(int j=0; j<5; j++)
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
              fMeanvalues_size[i][j][k].tab[c][ll][0] = int(fMeanvalues[i][j][k].vec[c][ll][1].size());
              fMeanvalues_size[i][j][k].tab[c][ll][0] = int(fMeanvalues[i][j][k].vec[c][ll][2].size());
              fMeanvalues_size[i][j][k].tab[c][ll][0] = int(fMeanvalues[i][j][k].vec[c][ll][3].size());
            }
          }
        }
      }
    }
  }


  cout <<
  fBP << " Best Primary (entries in disevent.root)\n\n" <<
  fRmu << " Reconstr. Mu (E_Beam>0)\n\n" <<
  fBMS << " BMS\n\n" <<
  fBEC << " Beam Energy Cuts\n\n" <<
  fTarg << " Event in Data Target\n\n" <<
  //fCell << " X Cells\n\n" <<
  fTrig << " O&IM Triggers\n\n" <<
  fQ2test << " Q>1\n\n" <<
  fYBjtest << " 0.1<y<0.7\n\n" <<
  fWBjtest << " 5<W<17\n\n" <<
  fXBjtest << " 0.004<Xbj<0.4\n\n" <<
  fXX0test << " XX0\n\n" <<
  fMom << " Momentum\n\n" <<
  fTRICH << " Theta RICH\n\n" <<
  fPosRICH << " Position RICH\n\n" <<
  fHplus << " h+\n\n" <<
  fHminus << " h-\n\n" <<
  fPiplus << " pi+\n\n" <<
  fPiminus << " pi-\n\n" <<
  fKplus << " K+\n\n" <<
  fKminus << " K-\n\n" <<
  fPplus << " p+\n\n" <<
  fPminus << " p-\n\n" <<
  "true pions : + " << fPiplus_true << " - " << fPiminus_true << "\n\n" <<
  "true kaons : + " << fKplus_true << " - " << fKminus_true << "\n\n" <<
  "true protons : + " << fPplus_true << " - " << fPminus_true << "\n\n";

  ofstream shout(Form("rawmult/%d/shout.txt",year), std::ofstream::out | std::ofstream::trunc);

  shout <<
  fBP << " Best Primary (entries in disevent.root)\n\n" <<
  fRmu << " Reconstr. Mu (E_Beam>0)\n\n" <<
  fBMS << " BMS\n\n" <<
  fBEC << " Beam Energy Cuts\n\n" <<
  fTarg << " Event in Data Target\n\n" <<
  //fCell << " X Cells\n\n" <<
  fTrig << " O&IM Triggers\n\n" <<
  fQ2test << " Q>1\n\n" <<
  fYBjtest << " 0.1<y<0.7\n\n" <<
  fWBjtest << " 5<W<17\n\n" <<
  fXX0test << " XX0\n\n" <<
  fMom << " Momentum\n\n" <<
  fTRICH << " Theta RICH\n\n" <<
  fPosRICH << " Position RICH\n\n" <<
  fHplus << " h+\n\n" <<
  fHminus << " h-\n\n" <<
  fPiplus << " pi+\n\n" <<
  fPiminus << " pi-\n\n" <<
  fKplus << " K+\n\n" <<
  fKminus << " K-\n\n" <<
  fPplus << " p+\n\n" <<
  fPminus << " p-\n\n" <<
  "true pions : + " << fPiplus_true << " - " << fPiminus_true << "\n\n" <<
  "true kaons : + " << fKplus_true << " - " << fKminus_true << "\n\n" <<
  "true protons : + " << fPplus_true << " - " << fPminus_true << "\n\n";

  shout.close();

  ofstream ofs_h(Form("rawmult/%d/hadron_0.txt",year), std::ofstream::out | std::ofstream::trunc);
  ofstream ofs_d("rawmult/%d/DIS_0.txt",year, std::ofstream::out | std::ofstream::trunc);
  ofstream xc("rawmult/%d/xcheck.txt",year, std::ofstream::out | std::ofstream::trunc);

  for(int c=0; c<2; c++)
  {
    for(int i=0; i<9; i++)
    {
      for(int j=0; j<5; j++)
      {
        for(int k=0; k<12; k++)
        {
          if(!c)
          {
            ofs_d << fNDIS_evt[0][i][j][k] << " " << fNDIS_evt_err[0][i][j][k];
          }

          for(int ll=0; ll<4; ll++)
          {
            ofs_d << " " << fMeanvalues_data[i][j][k].tab[c][ll][0] << " " <<
                            fMeanvalues_data[i][j][k].tab[c][ll][1] << " " <<
                            fMeanvalues_data[i][j][k].tab[c][ll][2] << " " <<
                            fMeanvalues_data[i][j][k].tab[c][ll][3] << " " << fMeanvalues_size[i][j][k].tab[c][ll][0];
          }
          ofs_d << endl;

          ofs_h << fBinning[i][j][k].tab[c][0][0] << " " << fBinning[i][j][k].tab[c][1][0] << " " << fBinning_loose[i][j][k].tab[c][0][0] << " " << fBinning_severe[i][j][k].tab[c][0][0] << " " <<
                   fBinning[i][j][k].tab[c][0][1] << " " << fBinning[i][j][k].tab[c][1][1] << " " << fBinning_loose[i][j][k].tab[c][0][1] << " " << fBinning_severe[i][j][k].tab[c][0][1] << " " <<
                   fBinning[i][j][k].tab[c][0][2] << " " << fBinning[i][j][k].tab[c][1][2] << " " << fBinning_loose[i][j][k].tab[c][0][2] << " " << fBinning_severe[i][j][k].tab[c][0][2] << " " <<
                   fBinning[i][j][k].tab[c][0][3] << " " << fBinning[i][j][k].tab[c][1][3] << " " << fBinning_loose[i][j][k].tab[c][0][3] << " " << fBinning_severe[i][j][k].tab[c][0][3] << " " << endl;

	  xc << c << " " << fXrange[i] << " " << fYrange[j] << " " << fZrange[k] << " " <<
		fNDIS_evt[0][i][j][k] << " " << fBinning[i][j][k].tab[c][0][0] << " " << fBinning[i][j][k].tab[c][0][1] << " " <<
		fBinning[i][j][k].tab[c][0][2] << " " << fBinning[i][j][k].tab[c][0][3] << endl;

	}
      }
    }
  }

  ofs_h.close();
  ofs_d.close();
  xc.close();

  return 0;
}

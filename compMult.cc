#include "compMult.h"

using namespace std;

void LoadMultiplicityFiles(string pfile1, strin pfile2)
{
  string sdum;
  Double_t x,y,z;

  ifstream mult1(pfile1);

  for(int i=0; i<9; i++)
  {
    for(int j=0; j<6; j++)
    {
      for(int k=0; k<12; k++)
      {
          mult1 >> x >> y >> z;
          cout << x << "\t" << y << "\t" << z << "\t" << endl;
          for(int l=0; l<4; l++) mult1 >> sdum;
          mult1 >> fMultiplicities1[i][j][k].tab[1][0][0] >> fMultiplicities1[i][j][k].tab[1][1][0] >> fMultiplicities1[i][j][k].tab[1][2][0];
          for(int l=0; l<5; l++) mult1 >> sdum;
          mult1 >> fMultiplicities1[i][j][k].tab[0][0][0] >> fMultiplicities1[i][j][k].tab[0][1][0] >> fMultiplicities1[i][j][k].tab[0][2][0];
          mult1 >> sdum;
      }
    }
  }
  mult1.close();

  ifstream mult2(pfile1);

  for(int i=0; i<9; i++)
  {
    for(int j=0; j<6; j++)
    {
      for(int k=0; k<12; k++)
      {
          mult >> x >> y >> z;
          cout << x << "\t" << y << "\t" << z << "\t" << endl;
          for(int l=0; l<4; l++) mult >> sdum;
          mult >> fMultiplicities2[i][j][k].tab[1][0][0] >> fMultiplicities2[i][j][k].tab[1][1][0] >> fMultiplicities2[i][j][k].tab[1][2][0];
          for(int l=0; l<5; l++) mult >> sdum;
          mult >> fMultiplicities2[i][j][k].tab[0][0][0] >> fMultiplicities2[i][j][k].tab[0][1][0] >> fMultiplicities2[i][j][k].tab[0][2][0];
          mult >> sdum;
      }
    }
  }
  mult2.close();
}

void yweightedavg()
{
  for(int c=0; c<2; c++)
  {
    for(int x=0; x<9; x++)
    {
      for(int z=0; z<12; z++)
      {
        for(int i=0; i<6; i++)
        {
          if(fMultiplicities1[x][i][z].tab[c][0][0])
          {
            fMultiplicities1_yavg[x][z].tab[c][0][0]+=fMultiplicities1[x][i][z].tab[c][0][0]/fMultiplicities1[x][i][z].tab[c][1][0];
            fMultiplicities1_yavg[x][z].tab[c][1][0]+=1/fMultiplicities1[x][i][z].tab[c][1][0];
            fMultiplicities1_yavg[x][z].tab[c][2][0]+=1/fMultiplicities1[x][i][z].tab[c][2][0];
          }
          if(fMultiplicities2[x][i][z].tab[c][0][0])
          {
            fMultiplicities2_yavg[x][z].tab[c][0][0]+=fMultiplicities2[x][i][z].tab[c][0][0]/fMultiplicities2[x][i][z].tab[c][1][0];
            fMultiplicities2_yavg[x][z].tab[c][1][0]+=1/fMultiplicities2[x][i][z].tab[c][1][0];
            fMultiplicities2_yavg[x][z].tab[c][2][0]+=1/fMultiplicities2[x][i][z].tab[c][2][0];
          }
        }
        if(fMultiplicities1_yavg[x][z].tab[c][0][l])
        {
          fMultiplicities1_yavg[x][z].tab[c][1][l]=1/fMultiplicities1_yavg[x][z].tab[c][1][l];
          fMultiplicities1_yavg[x][z].tab[c][2][l]=1/fMultiplicities1_yavg[x][z].tab[c][2][l];
          fMultiplicities1_yavg[x][z].tab[c][0][l]*=fMultiplicities1_yavg[x][z].tab[c][1][l];
        }
        if(fMultiplicities2_yavg[x][z].tab[c][0][l])
        {
          fMultiplicities2_yavg[x][z].tab[c][1][l]=1/fMultiplicities2_yavg[x][z].tab[c][1][l];
          fMultiplicities2_yavg[x][z].tab[c][2][l]=1/fMultiplicities2_yavg[x][z].tab[c][2][l];
          fMultiplicities2_yavg[x][z].tab[c][0][l]*=fMultiplicities2_yavg[x][z].tab[c][1][l];
        }
      }
    }
  }
}

int main(int argc, char **argv)
{

  LoadMultiplicityFiles(string argv[1], strin argv[2])

  r.push_back(fMultiplicities2[i][j][k].tab[0][0][0] ? fMultiplicities1[i][j][k].tab[1][0][0]/fMultiplicities2[i][j][k].tab[0][0][0] : 0);
  r_err.push_back(sqrt((fMultiplicities1[i][j][k].tab[1][1][0]+pow(fMultiplicities2[i][j][k].tab[0][1][0],2)*fMultiplicities1[i][j][k].tab[1][0][0]
                          /pow(fMultiplicities2[i][j][k].tab[0][0][0],2))/pow(fMultiplicities2[i][j][k].tab[0][0][0],2)));

  yweightedavg();

  r_y.push_back(fMultiplicities2_yavg[i][k].tab[0][0][0] ? fMultiplicities1_yavg[i][k].tab[1][0][0]/fMultiplicities2_yavg[i][k].tab[0][0][0] : 0);
  r_y_err.push_back(sqrt((fMultiplicities1_yavg[i][k].tab[1][1][0]+pow(fMultiplicities2_yavg[i][k].tab[0][1][0],2)*fMultiplicities1_yavg[i][k].tab[1][0][0]
                          /pow(fMultiplicities2_yavg[i][k].tab[0][0][0],2))/pow(fMultiplicities2_yavg[i][k].tab[0][0][0],2)));

  return 0;
}

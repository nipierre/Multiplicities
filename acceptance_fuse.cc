#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TMatrixD.h>
#include <TMatrixTUtils.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TGraphErrors.h>
#include <TArrow.h>

#include "acceptance_fuse.h"

using namespace std;

int main(int argc, char **argv)
{

  if(argc < 3)
  {
    cout << "ERROR : Not enough arguments." << endl;
    cout << "Asked : at least 2 *** Received : " << argc-1 << endl;
    cout << "./accfuse period filelist" << endl;

    return 1;
  }

  string cFilelist = argv[2];
  string periodName = argv[1];
  int year = 2016;
  int dummyd;

  ifstream list(cFilelist.c_str());
  string filename;

  while(list >> filename)
  {
    ifstream DIS_file(Form("%s/DIS_%s.txt",filename.c_str(),periodName.c_str()));
    ifstream DIS_zvtx_file(Form("%s/DIS_zvtx_%s.txt",filename.c_str(),periodName.c_str()));
    ifstream had_file(Form("%s/hadron_%s.txt",filename.c_str(),periodName.c_str()));
    ifstream had_zvtx_file(Form("%s/hadron_zvtx_%s.txt",filename.c_str(),periodName.c_str()));

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
              DIS_file >> dummyd;
              fNDIS_evt[0][i][j][k] += dummyd;
              DIS_file >> dummyd;
              fNDIS_evt_c[0][i][j][k] += dummyd;
              DIS_file >> dummyd;
              fNDIS_evt_MC[0][i][j][k] += dummyd;
              DIS_file >> dummyd;
              fNDIS_evt[1][i][j][k] += dummyd;
              DIS_file >> dummyd;
              fNDIS_evt_c[1][i][j][k] += dummyd;
              DIS_file >> dummyd;
              fNDIS_evt_MC[1][i][j][k] += dummyd;
              DIS_file >> dummyd;
              fNDIS_evt[2][i][j][k] += dummyd;
              DIS_file >> dummyd;
              fNDIS_evt_c[2][i][j][k] += dummyd;
              DIS_file >> dummyd;
              fNDIS_evt_MC[2][i][j][k] += dummyd;
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_zvtx[0][i][j][k][0] += dummyd;
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_zvtx[0][i][j][k][1] += dummyd;
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_zvtx[0][i][j][k][2] += dummyd;
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_zvtx[0][i][j][k][3] += dummyd;
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_MC_zvtx[0][i][j][k][0] += dummyd;
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_MC_zvtx[0][i][j][k][1] += dummyd;
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_MC_zvtx[0][i][j][k][2] += dummyd;
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_MC_zvtx[0][i][j][k][3] += dummyd;
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_zvtx[1][i][j][k][0] += dummyd;
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_zvtx[1][i][j][k][1] += dummyd;
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_zvtx[1][i][j][k][2] += dummyd;
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_zvtx[1][i][j][k][3] += dummyd;
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_MC_zvtx[1][i][j][k][0] += dummyd;
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_MC_zvtx[1][i][j][k][1] += dummyd;
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_MC_zvtx[1][i][j][k][2] += dummyd;
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_MC_zvtx[1][i][j][k][3] += dummyd;
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_zvtx[2][i][j][k][0] += dummyd;
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_zvtx[2][i][j][k][1] += dummyd;
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_zvtx[2][i][j][k][2] += dummyd;
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_zvtx[2][i][j][k][3] += dummyd;
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_MC_zvtx[2][i][j][k][0] += dummyd;
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_MC_zvtx[2][i][j][k][1] += dummyd;
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_MC_zvtx[2][i][j][k][2] += dummyd;
              DIS_zvtx_file >> dummyd;
              fNDIS_evt_MC_zvtx[2][i][j][k][3] += dummyd;

            }

            had_file >> dummyd;
            fRcstr[i][j][k].tab[c][0][0] += dummyd;
            had_file >> dummyd;
            fRcstr_c[i][j][k].tab[c][0][0] += dummyd;
            had_file >> dummyd;
            fGnrt[i][j][k].tab[c][0][0] += dummyd;

            had_file >> dummyd;
            fRcstr[i][j][k].tab[c][0][1] += dummyd;
            had_file >> dummyd;
            fRcstr_c[i][j][k].tab[c][0][1] += dummyd;
            had_file >> dummyd;
            fGnrt[i][j][k].tab[c][0][1] += dummyd;

            had_file >> dummyd;
            fRcstr[i][j][k].tab[c][0][2] += dummyd;
            had_file >> dummyd;
            fRcstr_c[i][j][k].tab[c][0][2] += dummyd;
            had_file >> dummyd;
            fGnrt[i][j][k].tab[c][0][2] += dummyd;

            had_file >> dummyd;
            fRcstr[i][j][k].tab[c][0][3] += dummyd;
            had_file >> dummyd;
            fRcstr_c[i][j][k].tab[c][0][3] += dummyd;
            had_file >> dummyd;
            fGnrt[i][j][k].tab[c][0][3] += dummyd;

            had_zvtx_file >> dummyd;
            fRcstr_zvtx[i][j][k][0].tab[c][0][0] += dummyd;
            had_zvtx_file >> dummyd;
            fRcstr_zvtx[i][j][k][1].tab[c][0][0] += dummyd;
            had_zvtx_file >> dummyd;
            fRcstr_zvtx[i][j][k][2].tab[c][0][0] += dummyd;
            had_zvtx_file >> dummyd;
            fRcstr_zvtx[i][j][k][3].tab[c][0][0] += dummyd;
            had_zvtx_file >> dummyd;
            fGnrt_zvtx[i][j][k][0].tab[c][0][0] += dummyd;
            had_zvtx_file >> dummyd;
            fGnrt_zvtx[i][j][k][1].tab[c][0][0] += dummyd;
            had_zvtx_file >> dummyd;
            fGnrt_zvtx[i][j][k][2].tab[c][0][0] += dummyd;
            had_zvtx_file >> dummyd;
            fGnrt_zvtx[i][j][k][3].tab[c][0][0] += dummyd;

            had_zvtx_file >> dummyd;
            fRcstr_zvtx[i][j][k][0].tab[c][0][1] += dummyd;
            had_zvtx_file >> dummyd;
            fRcstr_zvtx[i][j][k][1].tab[c][0][1] += dummyd;
            had_zvtx_file >> dummyd;
            fRcstr_zvtx[i][j][k][2].tab[c][0][1] += dummyd;
            had_zvtx_file >> dummyd;
            fRcstr_zvtx[i][j][k][3].tab[c][0][1] += dummyd;
            had_zvtx_file >> dummyd;
            fGnrt_zvtx[i][j][k][0].tab[c][0][1] += dummyd;
            had_zvtx_file >> dummyd;
            fGnrt_zvtx[i][j][k][1].tab[c][0][1] += dummyd;
            had_zvtx_file >> dummyd;
            fGnrt_zvtx[i][j][k][2].tab[c][0][1] += dummyd;
            had_zvtx_file >> dummyd;
            fGnrt_zvtx[i][j][k][3].tab[c][0][1] += dummyd;

            had_zvtx_file >> dummyd;
            fRcstr_zvtx[i][j][k][0].tab[c][0][2] += dummyd;
            had_zvtx_file >> dummyd;
            fRcstr_zvtx[i][j][k][1].tab[c][0][2] += dummyd;
            had_zvtx_file >> dummyd;
            fRcstr_zvtx[i][j][k][2].tab[c][0][2] += dummyd;
            had_zvtx_file >> dummyd;
            fRcstr_zvtx[i][j][k][3].tab[c][0][2] += dummyd;
            had_zvtx_file >> dummyd;
            fGnrt_zvtx[i][j][k][0].tab[c][0][2] += dummyd;
            had_zvtx_file >> dummyd;
            fGnrt_zvtx[i][j][k][1].tab[c][0][2] += dummyd;
            had_zvtx_file >> dummyd;
            fGnrt_zvtx[i][j][k][2].tab[c][0][2] += dummyd;
            had_zvtx_file >> dummyd;
            fGnrt_zvtx[i][j][k][3].tab[c][0][2] += dummyd;

            had_zvtx_file >> dummyd;
            fRcstr_zvtx[i][j][k][0].tab[c][0][3] += dummyd;
            had_zvtx_file >> dummyd;
            fRcstr_zvtx[i][j][k][1].tab[c][0][3] += dummyd;
            had_zvtx_file >> dummyd;
            fRcstr_zvtx[i][j][k][2].tab[c][0][3] += dummyd;
            had_zvtx_file >> dummyd;
            fRcstr_zvtx[i][j][k][3].tab[c][0][3] += dummyd;
            had_zvtx_file >> dummyd;
            fGnrt_zvtx[i][j][k][0].tab[c][0][3] += dummyd;
            had_zvtx_file >> dummyd;
            fGnrt_zvtx[i][j][k][1].tab[c][0][3] += dummyd;
            had_zvtx_file >> dummyd;
            fGnrt_zvtx[i][j][k][2].tab[c][0][3] += dummyd;
            had_zvtx_file >> dummyd;
            fGnrt_zvtx[i][j][k][3].tab[c][0][3] += dummyd;
          }
        }
      }
    }

    DIS_file.close();
    had_file.close();
    DIS_zvtx_file.close();
    had_zvtx_file.close();
  }

  ofstream ofs_h(Form("acceptance/%d/hadron/hadron_%s.txt",year,periodName.c_str()), std::ofstream::out | std::ofstream::trunc);
  ofstream ofs_hmult(Form("acceptance/%d/hadron/hadron_mult_%s.txt",year,periodName.c_str()), std::ofstream::out | std::ofstream::trunc);
  ofstream ofs_hzvtx(Form("acceptance/%d/hadron/hadron_zvtx_%s.txt",year,periodName.c_str()), std::ofstream::out | std::ofstream::trunc);
  ofstream ofs_hmult_zvtx(Form("acceptance/%d/hadron/hadron_mult_zvtx_%s.txt",year,periodName.c_str()), std::ofstream::out | std::ofstream::trunc);
  ofstream ofs_d(Form("acceptance/%d/DIS/DIS_%s.txt",year,periodName.c_str()), std::ofstream::out | std::ofstream::trunc);
  ofstream ofs_dmult(Form("acceptance/%d/DIS/DIS_mult_%s.txt",year,periodName.c_str()), std::ofstream::out | std::ofstream::trunc);
  ofstream ofs_dzvtx(Form("acceptance/%d/DIS/DIS_zvtx_%s.txt",year,periodName.c_str()), std::ofstream::out | std::ofstream::trunc);
  ofstream ofs_dmult_zvtx(Form("acceptance/%d/DIS/DIS_mult_zvtx_%s.txt",year,periodName.c_str()), std::ofstream::out | std::ofstream::trunc);

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
            ofs_d << fNDIS_evt[0][i][j][k] << " " << fNDIS_evt_c[0][i][j][k] << " " << fNDIS_evt_MC[0][i][j][k] << " " <<
                     fNDIS_evt[1][i][j][k] << " " << fNDIS_evt_c[1][i][j][k] << " " << fNDIS_evt_MC[1][i][j][k] << " " <<
                     fNDIS_evt[2][i][j][k] << " " << fNDIS_evt_c[2][i][j][k] << " " << fNDIS_evt_MC[2][i][j][k] << endl;

            ofs_dmult << fNDIS_evt[0][i][j][k] << " " << fNDIS_evt[0][i][j][k];

            for(int ll=0; ll<4; ll++)
            {
                ofs_dmult << " " << 0 << " " <<
                                    0 << " " <<
                                    0 << " " <<
                                    0 << " " << 0;
            }
            ofs_dmult << endl;

            ofs_dzvtx << fNDIS_evt_zvtx[0][i][j][k][0] << " " << fNDIS_evt_zvtx[0][i][j][k][1] << " " << fNDIS_evt_zvtx[0][i][j][k][2] << " " << fNDIS_evt_zvtx[0][i][j][k][3] << " " <<
                         fNDIS_evt_MC_zvtx[0][i][j][k][0] << " " << fNDIS_evt_MC_zvtx[0][i][j][k][1] << " " << fNDIS_evt_MC_zvtx[0][i][j][k][2] << " " << fNDIS_evt_MC_zvtx[0][i][j][k][3] << " " <<
                         fNDIS_evt_zvtx[1][i][j][k][0] << " " << fNDIS_evt_zvtx[1][i][j][k][1] << " " << fNDIS_evt_zvtx[1][i][j][k][2] << " " << fNDIS_evt_zvtx[1][i][j][k][3] << " " <<
                         fNDIS_evt_MC_zvtx[1][i][j][k][0] << " " << fNDIS_evt_MC_zvtx[1][i][j][k][1] << " " << fNDIS_evt_MC_zvtx[1][i][j][k][2] << " " << fNDIS_evt_MC_zvtx[1][i][j][k][3] << " " <<
                         fNDIS_evt_zvtx[2][i][j][k][0] << " " << fNDIS_evt_zvtx[2][i][j][k][1] << " " << fNDIS_evt_zvtx[2][i][j][k][2] << " " << fNDIS_evt_zvtx[2][i][j][k][3] << " " <<
                         fNDIS_evt_MC_zvtx[2][i][j][k][0] << " " << fNDIS_evt_MC_zvtx[2][i][j][k][1] << " " << fNDIS_evt_MC_zvtx[2][i][j][k][2] << " " << fNDIS_evt_MC_zvtx[2][i][j][k][3] << endl;

            for(int zv=0; zv<4; zv++)
                ofs_dmult_zvtx << fNDIS_evt_zvtx[0][i][j][k][zv] << " " << fNDIS_evt_zvtx[0][i][j][k][zv] << " ";

            ofs_dmult_zvtx << endl;
          }

          ofs_h << fRcstr[i][j][k].tab[c][0][0] << " " << fRcstr_c[i][j][k].tab[c][0][0] << " " << fGnrt[i][j][k].tab[c][0][0] << " " <<
                   fRcstr[i][j][k].tab[c][0][1] << " " << fRcstr_c[i][j][k].tab[c][0][1] << " " << fGnrt[i][j][k].tab[c][0][1] << " " <<
                   fRcstr[i][j][k].tab[c][0][2] << " " << fRcstr_c[i][j][k].tab[c][0][2] << " " << fGnrt[i][j][k].tab[c][0][2] << " " <<
                   fRcstr[i][j][k].tab[c][0][3] << " " << fRcstr_c[i][j][k].tab[c][0][3] << " " << fGnrt[i][j][k].tab[c][0][3] << " " << endl;

          ofs_hmult << fRcstr[i][j][k].tab[c][0][0] << " " << fRcstr[i][j][k].tab[c][0][0] << " " << 0 << " " << 0 << " " <<
                    fRcstr[i][j][k].tab[c][0][1] << " " << fRcstr[i][j][k].tab[c][0][1] << " " << 0 << " " << 0 << " " <<
                    fRcstr[i][j][k].tab[c][0][2] << " " << fRcstr[i][j][k].tab[c][0][2] << " " << 0 << " " << 0 << " " <<
                    fRcstr[i][j][k].tab[c][0][3] << " " << fRcstr[i][j][k].tab[c][0][3] << " " << 0 << " " << 0 << " " <<
          endl;

          ofs_hzvtx << fRcstr_zvtx[i][j][k][0].tab[c][0][0] << " " << fRcstr_zvtx[i][j][k][1].tab[c][0][0] << " " << fRcstr_zvtx[i][j][k][2].tab[c][0][0] << " " << fRcstr_zvtx[i][j][k][3].tab[c][0][0] << " " <<
                       fGnrt_zvtx[i][j][k][0].tab[c][0][0] << " " << fGnrt_zvtx[i][j][k][1].tab[c][0][0] << " " << fGnrt_zvtx[i][j][k][2].tab[c][0][0] << " " << fGnrt_zvtx[i][j][k][3].tab[c][0][0] << " " <<
                       fRcstr_zvtx[i][j][k][0].tab[c][0][1] << " " << fRcstr_zvtx[i][j][k][1].tab[c][0][1] << " " << fRcstr_zvtx[i][j][k][2].tab[c][0][1] << " " << fRcstr_zvtx[i][j][k][3].tab[c][0][1] << " " <<
                       fGnrt_zvtx[i][j][k][0].tab[c][0][1] << " " << fGnrt_zvtx[i][j][k][1].tab[c][0][1] << " " << fGnrt_zvtx[i][j][k][2].tab[c][0][1] << " " << fGnrt_zvtx[i][j][k][3].tab[c][0][1] << " " <<
                       fRcstr_zvtx[i][j][k][0].tab[c][0][2] << " " << fRcstr_zvtx[i][j][k][1].tab[c][0][2] << " " << fRcstr_zvtx[i][j][k][2].tab[c][0][2] << " " << fRcstr_zvtx[i][j][k][3].tab[c][0][2] << " " <<
                       fGnrt_zvtx[i][j][k][0].tab[c][0][2] << " " << fGnrt_zvtx[i][j][k][1].tab[c][0][2] << " " << fGnrt_zvtx[i][j][k][2].tab[c][0][2] << " " << fGnrt_zvtx[i][j][k][3].tab[c][0][2] << " " <<
                       fRcstr_zvtx[i][j][k][0].tab[c][0][3] << " " << fRcstr_zvtx[i][j][k][1].tab[c][0][3] << " " << fRcstr_zvtx[i][j][k][2].tab[c][0][3] << " " << fRcstr_zvtx[i][j][k][3].tab[c][0][3] << " " <<
                       fGnrt_zvtx[i][j][k][0].tab[c][0][3] << " " << fGnrt_zvtx[i][j][k][1].tab[c][0][3] << " " << fGnrt_zvtx[i][j][k][2].tab[c][0][3] << " " << fGnrt_zvtx[i][j][k][3].tab[c][0][3] << " " <<  endl;

         for(int zv=0; zv<4; zv++)
         {
           ofs_hmult_zvtx << fRcstr_zvtx[i][j][k][zv].tab[c][0][0] << " " << fRcstr_zvtx[i][j][k][zv].tab[c][0][0] << " " <<
                    fRcstr_zvtx[i][j][k][zv].tab[c][0][1] << " " << fRcstr_zvtx[i][j][k][zv].tab[c][0][1] << " " <<
                    fRcstr_zvtx[i][j][k][zv].tab[c][0][2] << " " << fRcstr_zvtx[i][j][k][zv].tab[c][0][2] << " " <<
                    fRcstr_zvtx[i][j][k][zv].tab[c][0][3] << " " << fRcstr_zvtx[i][j][k][zv].tab[c][0][3] << " " << endl;
         }
        }
      }
    }
  }

  ofs_h.close();
  ofs_hmult.close();
  ofs_d.close();
  ofs_dmult.close();
  ofs_hzvtx.close();
  ofs_hmult_zvtx.close();
  ofs_dzvtx.close();
  ofs_dmult_zvtx.close();

  return 0;
}

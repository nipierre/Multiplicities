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
    cout << "./accfuse period filelist_mu+ filelist_mu-" << endl;

    return 1;
  }

  string cFilelist1 = argv[2];
  string cFilelist2 = argv[3];
  string periodName = argv[1];
  int year = 2016;
  int dummyd;

  ifstream list1(cFilelist1.c_str());
  string filename;

  while(list1 >> filename)
  {
    ifstream DIS_file1(Form("%s/DIS_%s.txt",filename.c_str(),periodName.c_str()));
    ifstream DIS_zvtx_file1(Form("%s/DIS_zvtx_%s.txt",filename.c_str(),periodName.c_str()));
    ifstream had_file1(Form("%s/hadron_%s.txt",filename.c_str(),periodName.c_str()));
    ifstream had_zvtx_file1(Form("%s/hadron_zvtx_%s.txt",filename.c_str(),periodName.c_str()));
    ifstream electron_file1(Form("%s/electron_%s.txt",filename.c_str(),periodName.c_str()));
    ifstream electron_zvtx_file1(Form("%s/electron_zvtx_%s.txt",filename.c_str(),periodName.c_str()));

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
              DIS_file1 >> dummyd;
              fNDIS_evt[0][1][i][j][k] += dummyd;
              DIS_file1 >> dummyd;
              fNDIS_evt_c[0][1][i][j][k] += dummyd;
              DIS_file1 >> dummyd;
              fNDIS_evt_MC[0][1][i][j][k] += dummyd;
              DIS_file1 >> dummyd;
              fNDIS_evt[1][1][i][j][k] += dummyd;
              DIS_file1 >> dummyd;
              fNDIS_evt_c[1][1][i][j][k] += dummyd;
              DIS_file1 >> dummyd;
              fNDIS_evt_MC[1][1][i][j][k] += dummyd;
              DIS_file1 >> dummyd;
              fNDIS_evt[2][1][i][j][k] += dummyd;
              DIS_file1 >> dummyd;
              fNDIS_evt_c[2][1][i][j][k] += dummyd;
              DIS_file1 >> dummyd;
              fNDIS_evt_MC[2][1][i][j][k] += dummyd;
              DIS_zvtx_file1 >> dummyd;
              fNDIS_evt_zvtx[0][1][i][j][k][0] += dummyd;
              DIS_zvtx_file1 >> dummyd;
              fNDIS_evt_zvtx[0][1][i][j][k][1] += dummyd;
              DIS_zvtx_file1 >> dummyd;
              fNDIS_evt_zvtx[0][1][i][j][k][2] += dummyd;
              DIS_zvtx_file1 >> dummyd;
              fNDIS_evt_zvtx[0][1][i][j][k][3] += dummyd;
              DIS_zvtx_file1 >> dummyd;
              fNDIS_evt_MC_zvtx[0][1][i][j][k][0] += dummyd;
              DIS_zvtx_file1 >> dummyd;
              fNDIS_evt_MC_zvtx[0][1][i][j][k][1] += dummyd;
              DIS_zvtx_file1 >> dummyd;
              fNDIS_evt_MC_zvtx[0][1][i][j][k][2] += dummyd;
              DIS_zvtx_file1 >> dummyd;
              fNDIS_evt_MC_zvtx[0][1][i][j][k][3] += dummyd;
              DIS_zvtx_file1 >> dummyd;
              fNDIS_evt_zvtx[1][1][i][j][k][0] += dummyd;
              DIS_zvtx_file1 >> dummyd;
              fNDIS_evt_zvtx[1][1][i][j][k][1] += dummyd;
              DIS_zvtx_file1 >> dummyd;
              fNDIS_evt_zvtx[1][1][i][j][k][2] += dummyd;
              DIS_zvtx_file1 >> dummyd;
              fNDIS_evt_zvtx[1][1][i][j][k][3] += dummyd;
              DIS_zvtx_file1 >> dummyd;
              fNDIS_evt_MC_zvtx[1][1][i][j][k][0] += dummyd;
              DIS_zvtx_file1 >> dummyd;
              fNDIS_evt_MC_zvtx[1][1][i][j][k][1] += dummyd;
              DIS_zvtx_file1 >> dummyd;
              fNDIS_evt_MC_zvtx[1][1][i][j][k][2] += dummyd;
              DIS_zvtx_file1 >> dummyd;
              fNDIS_evt_MC_zvtx[1][1][i][j][k][3] += dummyd;
              DIS_zvtx_file1 >> dummyd;
              fNDIS_evt_zvtx[2][1][i][j][k][0] += dummyd;
              DIS_zvtx_file1 >> dummyd;
              fNDIS_evt_zvtx[2][1][i][j][k][1] += dummyd;
              DIS_zvtx_file1 >> dummyd;
              fNDIS_evt_zvtx[2][1][i][j][k][2] += dummyd;
              DIS_zvtx_file1 >> dummyd;
              fNDIS_evt_zvtx[2][1][i][j][k][3] += dummyd;
              DIS_zvtx_file1 >> dummyd;
              fNDIS_evt_MC_zvtx[2][1][i][j][k][0] += dummyd;
              DIS_zvtx_file1 >> dummyd;
              fNDIS_evt_MC_zvtx[2][1][i][j][k][1] += dummyd;
              DIS_zvtx_file1 >> dummyd;
              fNDIS_evt_MC_zvtx[2][1][i][j][k][2] += dummyd;
              DIS_zvtx_file1 >> dummyd;
              fNDIS_evt_MC_zvtx[2][1][i][j][k][3] += dummyd;

            }

            had_file1 >> dummyd;
            fRcstr[i][j][k].tab[c][1][0][0] += dummyd;
            had_file1 >> dummyd;
            fRcstr_c[i][j][k].tab[c][1][0][0] += dummyd;
            had_file1 >> dummyd;
            fGnrt[i][j][k].tab[c][1][0][0] += dummyd;

            had_file1 >> dummyd;
            fRcstr[i][j][k].tab[c][1][0][1] += dummyd;
            had_file1 >> dummyd;
            fRcstr_c[i][j][k].tab[c][1][0][1] += dummyd;
            had_file1 >> dummyd;
            fGnrt[i][j][k].tab[c][1][0][1] += dummyd;

            had_file1 >> dummyd;
            fRcstr[i][j][k].tab[c][1][0][2] += dummyd;
            had_file1 >> dummyd;
            fRcstr_c[i][j][k].tab[c][1][0][2] += dummyd;
            had_file1 >> dummyd;
            fGnrt[i][j][k].tab[c][1][0][2] += dummyd;

            had_file1 >> dummyd;
            fRcstr[i][j][k].tab[c][1][0][3] += dummyd;
            had_file1 >> dummyd;
            fRcstr_c[i][j][k].tab[c][1][0][3] += dummyd;
            had_file1 >> dummyd;
            fGnrt[i][j][k].tab[c][1][0][3] += dummyd;

            had_zvtx_file1 >> dummyd;
            fRcstr_zvtx[i][j][k][0].tab[c][1][0][0] += dummyd;
            had_zvtx_file1 >> dummyd;
            fRcstr_zvtx[i][j][k][1].tab[c][1][0][0] += dummyd;
            had_zvtx_file1 >> dummyd;
            fRcstr_zvtx[i][j][k][2].tab[c][1][0][0] += dummyd;
            had_zvtx_file1 >> dummyd;
            fRcstr_zvtx[i][j][k][3].tab[c][1][0][0] += dummyd;
            had_zvtx_file1 >> dummyd;
            fGnrt_zvtx[i][j][k][0].tab[c][1][0][0] += dummyd;
            had_zvtx_file1 >> dummyd;
            fGnrt_zvtx[i][j][k][1].tab[c][1][0][0] += dummyd;
            had_zvtx_file1 >> dummyd;
            fGnrt_zvtx[i][j][k][2].tab[c][1][0][0] += dummyd;
            had_zvtx_file1 >> dummyd;
            fGnrt_zvtx[i][j][k][3].tab[c][1][0][0] += dummyd;

            had_zvtx_file1 >> dummyd;
            fRcstr_zvtx[i][j][k][0].tab[c][1][0][1] += dummyd;
            had_zvtx_file1 >> dummyd;
            fRcstr_zvtx[i][j][k][1].tab[c][1][0][1] += dummyd;
            had_zvtx_file1 >> dummyd;
            fRcstr_zvtx[i][j][k][2].tab[c][1][0][1] += dummyd;
            had_zvtx_file1 >> dummyd;
            fRcstr_zvtx[i][j][k][3].tab[c][1][0][1] += dummyd;
            had_zvtx_file1 >> dummyd;
            fGnrt_zvtx[i][j][k][0].tab[c][1][0][1] += dummyd;
            had_zvtx_file1 >> dummyd;
            fGnrt_zvtx[i][j][k][1].tab[c][1][0][1] += dummyd;
            had_zvtx_file1 >> dummyd;
            fGnrt_zvtx[i][j][k][2].tab[c][1][0][1] += dummyd;
            had_zvtx_file1 >> dummyd;
            fGnrt_zvtx[i][j][k][3].tab[c][1][0][1] += dummyd;

            had_zvtx_file1 >> dummyd;
            fRcstr_zvtx[i][j][k][0].tab[c][1][0][2] += dummyd;
            had_zvtx_file1 >> dummyd;
            fRcstr_zvtx[i][j][k][1].tab[c][1][0][2] += dummyd;
            had_zvtx_file1 >> dummyd;
            fRcstr_zvtx[i][j][k][2].tab[c][1][0][2] += dummyd;
            had_zvtx_file1 >> dummyd;
            fRcstr_zvtx[i][j][k][3].tab[c][1][0][2] += dummyd;
            had_zvtx_file1 >> dummyd;
            fGnrt_zvtx[i][j][k][0].tab[c][1][0][2] += dummyd;
            had_zvtx_file1 >> dummyd;
            fGnrt_zvtx[i][j][k][1].tab[c][1][0][2] += dummyd;
            had_zvtx_file1 >> dummyd;
            fGnrt_zvtx[i][j][k][2].tab[c][1][0][2] += dummyd;
            had_zvtx_file1 >> dummyd;
            fGnrt_zvtx[i][j][k][3].tab[c][1][0][2] += dummyd;

            had_zvtx_file1 >> dummyd;
            fRcstr_zvtx[i][j][k][0].tab[c][1][0][3] += dummyd;
            had_zvtx_file1 >> dummyd;
            fRcstr_zvtx[i][j][k][1].tab[c][1][0][3] += dummyd;
            had_zvtx_file1 >> dummyd;
            fRcstr_zvtx[i][j][k][2].tab[c][1][0][3] += dummyd;
            had_zvtx_file1 >> dummyd;
            fRcstr_zvtx[i][j][k][3].tab[c][1][0][3] += dummyd;
            had_zvtx_file1 >> dummyd;
            fGnrt_zvtx[i][j][k][0].tab[c][1][0][3] += dummyd;
            had_zvtx_file1 >> dummyd;
            fGnrt_zvtx[i][j][k][1].tab[c][1][0][3] += dummyd;
            had_zvtx_file1 >> dummyd;
            fGnrt_zvtx[i][j][k][2].tab[c][1][0][3] += dummyd;
            had_zvtx_file1 >> dummyd;
            fGnrt_zvtx[i][j][k][3].tab[c][1][0][3] += dummyd;

            electron_file1 >> dummyd;
            fRcstr[i][j][k].tab[c][1][0][4] += dummyd;

            electron_zvtx_file1 >> dummyd;
            fRcstr_zvtx[i][j][k][0].tab[c][1][0][4] += dummyd;
            electron_zvtx_file1 >> dummyd;
            fRcstr_zvtx[i][j][k][1].tab[c][1][0][4] += dummyd;
            electron_zvtx_file1 >> dummyd;
            fRcstr_zvtx[i][j][k][2].tab[c][1][0][4] += dummyd;
            electron_zvtx_file1 >> dummyd;
            fRcstr_zvtx[i][j][k][3].tab[c][1][0][4] += dummyd;
          }
        }
      }
    }

    DIS_file1.close();
    had_file1.close();
    electron_file1.close();
    DIS_zvtx_file1.close();
    had_zvtx_file1.close();
    electron_zvtx_file1.close();
  }

  ifstream list2(cFilelist2.c_str());

  while(list2 >> filename)
  {
    ifstream DIS_file2(Form("%s/DIS_%s.txt",filename.c_str(),periodName.c_str()));
    ifstream DIS_zvtx_file2(Form("%s/DIS_zvtx_%s.txt",filename.c_str(),periodName.c_str()));
    ifstream had_file2(Form("%s/hadron_%s.txt",filename.c_str(),periodName.c_str()));
    ifstream had_zvtx_file2(Form("%s/hadron_zvtx_%s.txt",filename.c_str(),periodName.c_str()));
    ifstream electron_file2(Form("%s/electron_%s.txt",filename.c_str(),periodName.c_str()));
    ifstream electron_zvtx_file2(Form("%s/electron_zvtx_%s.txt",filename.c_str(),periodName.c_str()));

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
              DIS_file2 >> dummyd;
              fNDIS_evt[0][0][i][j][k] += dummyd;
              DIS_file2 >> dummyd;
              fNDIS_evt_c[0][0][i][j][k] += dummyd;
              DIS_file2 >> dummyd;
              fNDIS_evt_MC[0][0][i][j][k] += dummyd;
              DIS_file2 >> dummyd;
              fNDIS_evt[1][0][i][j][k] += dummyd;
              DIS_file2 >> dummyd;
              fNDIS_evt_c[1][0][i][j][k] += dummyd;
              DIS_file2 >> dummyd;
              fNDIS_evt_MC[1][0][i][j][k] += dummyd;
              DIS_file2 >> dummyd;
              fNDIS_evt[2][0][i][j][k] += dummyd;
              DIS_file2 >> dummyd;
              fNDIS_evt_c[2][0][i][j][k] += dummyd;
              DIS_file2 >> dummyd;
              fNDIS_evt_MC[2][0][i][j][k] += dummyd;
              DIS_zvtx_file2 >> dummyd;
              fNDIS_evt_zvtx[0][0][i][j][k][0] += dummyd;
              DIS_zvtx_file2 >> dummyd;
              fNDIS_evt_zvtx[0][0][i][j][k][1] += dummyd;
              DIS_zvtx_file2 >> dummyd;
              fNDIS_evt_zvtx[0][0][i][j][k][2] += dummyd;
              DIS_zvtx_file2 >> dummyd;
              fNDIS_evt_zvtx[0][0][i][j][k][3] += dummyd;
              DIS_zvtx_file2 >> dummyd;
              fNDIS_evt_MC_zvtx[0][0][i][j][k][0] += dummyd;
              DIS_zvtx_file2 >> dummyd;
              fNDIS_evt_MC_zvtx[0][0][i][j][k][1] += dummyd;
              DIS_zvtx_file2 >> dummyd;
              fNDIS_evt_MC_zvtx[0][0][i][j][k][2] += dummyd;
              DIS_zvtx_file2 >> dummyd;
              fNDIS_evt_MC_zvtx[0][0][i][j][k][3] += dummyd;
              DIS_zvtx_file2 >> dummyd;
              fNDIS_evt_zvtx[1][0][i][j][k][0] += dummyd;
              DIS_zvtx_file2 >> dummyd;
              fNDIS_evt_zvtx[1][0][i][j][k][1] += dummyd;
              DIS_zvtx_file2 >> dummyd;
              fNDIS_evt_zvtx[1][0][i][j][k][2] += dummyd;
              DIS_zvtx_file2 >> dummyd;
              fNDIS_evt_zvtx[1][0][i][j][k][3] += dummyd;
              DIS_zvtx_file2 >> dummyd;
              fNDIS_evt_MC_zvtx[1][0][i][j][k][0] += dummyd;
              DIS_zvtx_file2 >> dummyd;
              fNDIS_evt_MC_zvtx[1][0][i][j][k][1] += dummyd;
              DIS_zvtx_file2 >> dummyd;
              fNDIS_evt_MC_zvtx[1][0][i][j][k][2] += dummyd;
              DIS_zvtx_file2 >> dummyd;
              fNDIS_evt_MC_zvtx[1][0][i][j][k][3] += dummyd;
              DIS_zvtx_file2 >> dummyd;
              fNDIS_evt_zvtx[2][0][i][j][k][0] += dummyd;
              DIS_zvtx_file2 >> dummyd;
              fNDIS_evt_zvtx[2][0][i][j][k][1] += dummyd;
              DIS_zvtx_file2 >> dummyd;
              fNDIS_evt_zvtx[2][0][i][j][k][2] += dummyd;
              DIS_zvtx_file2 >> dummyd;
              fNDIS_evt_zvtx[2][0][i][j][k][3] += dummyd;
              DIS_zvtx_file2 >> dummyd;
              fNDIS_evt_MC_zvtx[2][0][i][j][k][0] += dummyd;
              DIS_zvtx_file2 >> dummyd;
              fNDIS_evt_MC_zvtx[2][0][i][j][k][1] += dummyd;
              DIS_zvtx_file2 >> dummyd;
              fNDIS_evt_MC_zvtx[2][0][i][j][k][2] += dummyd;
              DIS_zvtx_file2 >> dummyd;
              fNDIS_evt_MC_zvtx[2][0][i][j][k][3] += dummyd;

            }

            had_file2 >> dummyd;
            fRcstr[i][j][k].tab[c][0][0][0] += dummyd;
            had_file2 >> dummyd;
            fRcstr_c[i][j][k].tab[c][0][0][0] += dummyd;
            had_file2 >> dummyd;
            fGnrt[i][j][k].tab[c][0][0][0] += dummyd;

            had_file2 >> dummyd;
            fRcstr[i][j][k].tab[c][0][0][1] += dummyd;
            had_file2 >> dummyd;
            fRcstr_c[i][j][k].tab[c][0][0][1] += dummyd;
            had_file2 >> dummyd;
            fGnrt[i][j][k].tab[c][0][0][1] += dummyd;

            had_file2 >> dummyd;
            fRcstr[i][j][k].tab[c][0][0][2] += dummyd;
            had_file2 >> dummyd;
            fRcstr_c[i][j][k].tab[c][0][0][2] += dummyd;
            had_file2 >> dummyd;
            fGnrt[i][j][k].tab[c][0][0][2] += dummyd;

            had_file2 >> dummyd;
            fRcstr[i][j][k].tab[c][0][0][3] += dummyd;
            had_file2 >> dummyd;
            fRcstr_c[i][j][k].tab[c][0][0][3] += dummyd;
            had_file2 >> dummyd;
            fGnrt[i][j][k].tab[c][0][0][3] += dummyd;

            had_zvtx_file2 >> dummyd;
            fRcstr_zvtx[i][j][k][0].tab[c][0][0][0] += dummyd;
            had_zvtx_file2 >> dummyd;
            fRcstr_zvtx[i][j][k][1].tab[c][0][0][0] += dummyd;
            had_zvtx_file2 >> dummyd;
            fRcstr_zvtx[i][j][k][2].tab[c][0][0][0] += dummyd;
            had_zvtx_file2 >> dummyd;
            fRcstr_zvtx[i][j][k][3].tab[c][0][0][0] += dummyd;
            had_zvtx_file2 >> dummyd;
            fGnrt_zvtx[i][j][k][0].tab[c][0][0][0] += dummyd;
            had_zvtx_file2 >> dummyd;
            fGnrt_zvtx[i][j][k][1].tab[c][0][0][0] += dummyd;
            had_zvtx_file2 >> dummyd;
            fGnrt_zvtx[i][j][k][2].tab[c][0][0][0] += dummyd;
            had_zvtx_file2 >> dummyd;
            fGnrt_zvtx[i][j][k][3].tab[c][0][0][0] += dummyd;

            had_zvtx_file2 >> dummyd;
            fRcstr_zvtx[i][j][k][0].tab[c][0][0][1] += dummyd;
            had_zvtx_file2 >> dummyd;
            fRcstr_zvtx[i][j][k][1].tab[c][0][0][1] += dummyd;
            had_zvtx_file2 >> dummyd;
            fRcstr_zvtx[i][j][k][2].tab[c][0][0][1] += dummyd;
            had_zvtx_file2 >> dummyd;
            fRcstr_zvtx[i][j][k][3].tab[c][0][0][1] += dummyd;
            had_zvtx_file2 >> dummyd;
            fGnrt_zvtx[i][j][k][0].tab[c][0][0][1] += dummyd;
            had_zvtx_file2 >> dummyd;
            fGnrt_zvtx[i][j][k][1].tab[c][0][0][1] += dummyd;
            had_zvtx_file2 >> dummyd;
            fGnrt_zvtx[i][j][k][2].tab[c][0][0][1] += dummyd;
            had_zvtx_file2 >> dummyd;
            fGnrt_zvtx[i][j][k][3].tab[c][0][0][1] += dummyd;

            had_zvtx_file2 >> dummyd;
            fRcstr_zvtx[i][j][k][0].tab[c][0][0][2] += dummyd;
            had_zvtx_file2 >> dummyd;
            fRcstr_zvtx[i][j][k][1].tab[c][0][0][2] += dummyd;
            had_zvtx_file2 >> dummyd;
            fRcstr_zvtx[i][j][k][2].tab[c][0][0][2] += dummyd;
            had_zvtx_file2 >> dummyd;
            fRcstr_zvtx[i][j][k][3].tab[c][0][0][2] += dummyd;
            had_zvtx_file2 >> dummyd;
            fGnrt_zvtx[i][j][k][0].tab[c][0][0][2] += dummyd;
            had_zvtx_file2 >> dummyd;
            fGnrt_zvtx[i][j][k][1].tab[c][0][0][2] += dummyd;
            had_zvtx_file2 >> dummyd;
            fGnrt_zvtx[i][j][k][2].tab[c][0][0][2] += dummyd;
            had_zvtx_file2 >> dummyd;
            fGnrt_zvtx[i][j][k][3].tab[c][0][0][2] += dummyd;

            had_zvtx_file2 >> dummyd;
            fRcstr_zvtx[i][j][k][0].tab[c][0][0][3] += dummyd;
            had_zvtx_file2 >> dummyd;
            fRcstr_zvtx[i][j][k][1].tab[c][0][0][3] += dummyd;
            had_zvtx_file2 >> dummyd;
            fRcstr_zvtx[i][j][k][2].tab[c][0][0][3] += dummyd;
            had_zvtx_file2 >> dummyd;
            fRcstr_zvtx[i][j][k][3].tab[c][0][0][3] += dummyd;
            had_zvtx_file2 >> dummyd;
            fGnrt_zvtx[i][j][k][0].tab[c][0][0][3] += dummyd;
            had_zvtx_file2 >> dummyd;
            fGnrt_zvtx[i][j][k][1].tab[c][0][0][3] += dummyd;
            had_zvtx_file2 >> dummyd;
            fGnrt_zvtx[i][j][k][2].tab[c][0][0][3] += dummyd;
            had_zvtx_file2 >> dummyd;
            fGnrt_zvtx[i][j][k][3].tab[c][0][0][3] += dummyd;

            electron_file2 >> dummyd;
            fRcstr[i][j][k].tab[c][0][0][4] += dummyd;

            electron_zvtx_file2 >> dummyd;
            fRcstr_zvtx[i][j][k][0].tab[c][0][0][4] += dummyd;
            electron_zvtx_file2 >> dummyd;
            fRcstr_zvtx[i][j][k][1].tab[c][0][0][4] += dummyd;
            electron_zvtx_file2 >> dummyd;
            fRcstr_zvtx[i][j][k][2].tab[c][0][0][4] += dummyd;
            electron_zvtx_file2 >> dummyd;
            fRcstr_zvtx[i][j][k][3].tab[c][0][0][4] += dummyd;
          }
        }
      }
    }

    DIS_file2.close();
    had_file2.close();
    electron_file2.close();
    DIS_zvtx_file2.close();
    had_zvtx_file2.close();
    electron_zvtx_file2.close();
  }

  ofstream ofs_h(Form("acceptance/%d/hadron/hadron_%s.txt",year,periodName.c_str()), std::ofstream::out | std::ofstream::trunc);
  ofstream ofs_hmult(Form("acceptance/%d/hadron/hadron_mult_%s.txt",year,periodName.c_str()), std::ofstream::out | std::ofstream::trunc);
  ofstream ofs_hzvtx(Form("acceptance/%d/hadron/hadron_zvtx_%s.txt",year,periodName.c_str()), std::ofstream::out | std::ofstream::trunc);
  ofstream ofs_hmult_zvtx(Form("acceptance/%d/hadron/hadron_mult_zvtx_%s.txt",year,periodName.c_str()), std::ofstream::out | std::ofstream::trunc);
  ofstream ofs_d(Form("acceptance/%d/DIS/DIS_%s.txt",year,periodName.c_str()), std::ofstream::out | std::ofstream::trunc);
  ofstream ofs_dmult(Form("acceptance/%d/DIS/DIS_mult_%s.txt",year,periodName.c_str()), std::ofstream::out | std::ofstream::trunc);
  ofstream ofs_dzvtx(Form("acceptance/%d/DIS/DIS_zvtx_%s.txt",year,periodName.c_str()), std::ofstream::out | std::ofstream::trunc);
  ofstream ofs_dmult_zvtx(Form("acceptance/%d/DIS/DIS_mult_zvtx_%s.txt",year,periodName.c_str()), std::ofstream::out | std::ofstream::trunc);
  ofstream ofs_e(Form("acceptance/%d/electro/electron_%s.txt",year,periodName.c_str()), std::ofstream::out | std::ofstream::trunc);
  ofstream ofs_ezvtx(Form("acceptance/%d/electro/electron_zvtx_%s.txt",year,periodName.c_str()), std::ofstream::out | std::ofstream::trunc);

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
            ofs_d << fNDIS_evt[0][1][i][j][k] << " " << fNDIS_evt_c[0][1][i][j][k] << " " << fNDIS_evt_MC[0][1][i][j][k] << " " <<
                     fNDIS_evt[1][1][i][j][k] << " " << fNDIS_evt_c[1][1][i][j][k] << " " << fNDIS_evt_MC[1][1][i][j][k] << " " <<
                     fNDIS_evt[2][1][i][j][k] << " " << fNDIS_evt_c[2][1][i][j][k] << " " << fNDIS_evt_MC[2][1][i][j][k] << " " <<
                     fNDIS_evt[0][0][i][j][k] << " " << fNDIS_evt_c[0][0][i][j][k] << " " << fNDIS_evt_MC[0][0][i][j][k] << " " <<
                     fNDIS_evt[1][0][i][j][k] << " " << fNDIS_evt_c[1][0][i][j][k] << " " << fNDIS_evt_MC[1][0][i][j][k] << " " <<
                     fNDIS_evt[2][0][i][j][k] << " " << fNDIS_evt_c[2][0][i][j][k] << " " << fNDIS_evt_MC[2][0][i][j][k] << endl;

            ofs_dmult << fNDIS_evt[0][1][i][j][k] << " " << fNDIS_evt[0][1][i][j][k] << " " << fNDIS_evt[0][0][i][j][k] << " " << fNDIS_evt[0][0][i][j][k];

            for(int ll=0; ll<4; ll++)
            {
                ofs_dmult << " " << 0 << " " <<
                                    0 << " " <<
                                    0 << " " <<
                                    0 << " " << 0;
            }
            ofs_dmult << endl;

            ofs_dzvtx << fNDIS_evt_zvtx[0][1][i][j][k][0] << " " << fNDIS_evt_zvtx[0][1][i][j][k][1] << " " << fNDIS_evt_zvtx[0][1][i][j][k][2] << " " << fNDIS_evt_zvtx[0][1][i][j][k][3] << " " <<
                         fNDIS_evt_MC_zvtx[0][1][i][j][k][0] << " " << fNDIS_evt_MC_zvtx[0][1][i][j][k][1] << " " << fNDIS_evt_MC_zvtx[0][1][i][j][k][2] << " " << fNDIS_evt_MC_zvtx[0][1][i][j][k][3] << " " <<
                         fNDIS_evt_zvtx[1][1][i][j][k][0] << " " << fNDIS_evt_zvtx[1][1][i][j][k][1] << " " << fNDIS_evt_zvtx[1][1][i][j][k][2] << " " << fNDIS_evt_zvtx[1][1][i][j][k][3] << " " <<
                         fNDIS_evt_MC_zvtx[1][1][i][j][k][0] << " " << fNDIS_evt_MC_zvtx[1][1][i][j][k][1] << " " << fNDIS_evt_MC_zvtx[1][1][i][j][k][2] << " " << fNDIS_evt_MC_zvtx[1][1][i][j][k][3] << " " <<
                         fNDIS_evt_zvtx[2][1][i][j][k][0] << " " << fNDIS_evt_zvtx[2][1][i][j][k][1] << " " << fNDIS_evt_zvtx[2][1][i][j][k][2] << " " << fNDIS_evt_zvtx[2][1][i][j][k][3] << " " <<
                         fNDIS_evt_MC_zvtx[2][1][i][j][k][0] << " " << fNDIS_evt_MC_zvtx[2][1][i][j][k][1] << " " << fNDIS_evt_MC_zvtx[2][1][i][j][k][2] << " " << fNDIS_evt_MC_zvtx[2][1][i][j][k][3] << " " <<
                         fNDIS_evt_zvtx[0][0][i][j][k][0] << " " << fNDIS_evt_zvtx[0][0][i][j][k][1] << " " << fNDIS_evt_zvtx[0][0][i][j][k][2] << " " << fNDIS_evt_zvtx[0][0][i][j][k][3] << " " <<
                         fNDIS_evt_MC_zvtx[0][0][i][j][k][0] << " " << fNDIS_evt_MC_zvtx[0][0][i][j][k][1] << " " << fNDIS_evt_MC_zvtx[0][0][i][j][k][2] << " " << fNDIS_evt_MC_zvtx[0][0][i][j][k][3] << " " <<
                         fNDIS_evt_zvtx[1][0][i][j][k][0] << " " << fNDIS_evt_zvtx[1][0][i][j][k][1] << " " << fNDIS_evt_zvtx[1][0][i][j][k][2] << " " << fNDIS_evt_zvtx[1][0][i][j][k][3] << " " <<
                         fNDIS_evt_MC_zvtx[1][0][i][j][k][0] << " " << fNDIS_evt_MC_zvtx[1][0][i][j][k][1] << " " << fNDIS_evt_MC_zvtx[1][0][i][j][k][2] << " " << fNDIS_evt_MC_zvtx[1][0][i][j][k][3] << " " <<
                         fNDIS_evt_zvtx[2][0][i][j][k][0] << " " << fNDIS_evt_zvtx[2][0][i][j][k][1] << " " << fNDIS_evt_zvtx[2][0][i][j][k][2] << " " << fNDIS_evt_zvtx[2][0][i][j][k][3] << " " <<
                         fNDIS_evt_MC_zvtx[2][0][i][j][k][0] << " " << fNDIS_evt_MC_zvtx[2][0][i][j][k][1] << " " << fNDIS_evt_MC_zvtx[2][0][i][j][k][2] << " " << fNDIS_evt_MC_zvtx[2][0][i][j][k][3] << endl;

            for(int zv=0; zv<4; zv++)
                ofs_dmult_zvtx << fNDIS_evt_zvtx[0][1][i][j][k][zv] << " " << fNDIS_evt_zvtx[0][1][i][j][k][zv] << " " << fNDIS_evt_zvtx[0][0][i][j][k][zv] << " " << fNDIS_evt_zvtx[0][0][i][j][k][zv] << " ";

            ofs_dmult_zvtx << endl;
          }

          ofs_h << fRcstr[i][j][k].tab[c][1][0][0] << " " << fRcstr_c[i][j][k].tab[c][1][0][0] << " " << fGnrt[i][j][k].tab[c][1][0][0] << " " <<
                   fRcstr[i][j][k].tab[c][1][0][1] << " " << fRcstr_c[i][j][k].tab[c][1][0][1] << " " << fGnrt[i][j][k].tab[c][1][0][1] << " " <<
                   fRcstr[i][j][k].tab[c][1][0][2] << " " << fRcstr_c[i][j][k].tab[c][1][0][2] << " " << fGnrt[i][j][k].tab[c][1][0][2] << " " <<
                   fRcstr[i][j][k].tab[c][1][0][3] << " " << fRcstr_c[i][j][k].tab[c][1][0][3] << " " << fGnrt[i][j][k].tab[c][1][0][3] << " " <<
                   fRcstr[i][j][k].tab[c][0][0][0] << " " << fRcstr_c[i][j][k].tab[c][0][0][0] << " " << fGnrt[i][j][k].tab[c][0][0][0] << " " <<
                   fRcstr[i][j][k].tab[c][0][0][1] << " " << fRcstr_c[i][j][k].tab[c][0][0][1] << " " << fGnrt[i][j][k].tab[c][0][0][1] << " " <<
                   fRcstr[i][j][k].tab[c][0][0][2] << " " << fRcstr_c[i][j][k].tab[c][0][0][2] << " " << fGnrt[i][j][k].tab[c][0][0][2] << " " <<
                   fRcstr[i][j][k].tab[c][0][0][3] << " " << fRcstr_c[i][j][k].tab[c][0][0][3] << " " << fGnrt[i][j][k].tab[c][0][0][3] << " " << endl;

          ofs_hmult << fRcstr[i][j][k].tab[c][1][0][0] << " " << fRcstr[i][j][k].tab[c][1][0][0] << " " << 0 << " " << 0 << " " <<
                       fRcstr[i][j][k].tab[c][1][0][1] << " " << fRcstr[i][j][k].tab[c][1][0][1] << " " << 0 << " " << 0 << " " <<
                       fRcstr[i][j][k].tab[c][1][0][2] << " " << fRcstr[i][j][k].tab[c][1][0][2] << " " << 0 << " " << 0 << " " <<
                       fRcstr[i][j][k].tab[c][1][0][3] << " " << fRcstr[i][j][k].tab[c][1][0][3] << " " << 0 << " " << 0 << " " <<
                       fRcstr[i][j][k].tab[c][0][0][0] << " " << fRcstr[i][j][k].tab[c][0][0][0] << " " << 0 << " " << 0 << " " <<
                       fRcstr[i][j][k].tab[c][0][0][1] << " " << fRcstr[i][j][k].tab[c][0][0][1] << " " << 0 << " " << 0 << " " <<
                       fRcstr[i][j][k].tab[c][0][0][2] << " " << fRcstr[i][j][k].tab[c][0][0][2] << " " << 0 << " " << 0 << " " <<
                       fRcstr[i][j][k].tab[c][0][0][3] << " " << fRcstr[i][j][k].tab[c][0][0][3] << " " << 0 << " " << 0 << " " << endl;

          ofs_e << fRcstr[i][j][k].tab[c][1][0][0]/(fRcstr[i][j][k].tab[c][1][0][0]+fRcstr[i][j][k].tab[c][1][0][4]) << " "
                << fRcstr[i][j][k].tab[c][0][0][0]/(fRcstr[i][j][k].tab[c][0][0][0]+fRcstr[i][j][k].tab[c][0][0][4]) << endl;

          ofs_hzvtx << fRcstr_zvtx[i][j][k][0].tab[c][1][0][0] << " " << fRcstr_zvtx[i][j][k][1].tab[c][1][0][0] << " " << fRcstr_zvtx[i][j][k][2].tab[c][1][0][0] << " " << fRcstr_zvtx[i][j][k][3].tab[c][1][0][0] << " " <<
                       fGnrt_zvtx[i][j][k][0].tab[c][1][0][0] << " " << fGnrt_zvtx[i][j][k][1].tab[c][1][0][0] << " " << fGnrt_zvtx[i][j][k][2].tab[c][1][0][0] << " " << fGnrt_zvtx[i][j][k][3].tab[c][1][0][0] << " " <<
                       fRcstr_zvtx[i][j][k][0].tab[c][1][0][1] << " " << fRcstr_zvtx[i][j][k][1].tab[c][1][0][1] << " " << fRcstr_zvtx[i][j][k][2].tab[c][1][0][1] << " " << fRcstr_zvtx[i][j][k][3].tab[c][1][0][1] << " " <<
                       fGnrt_zvtx[i][j][k][0].tab[c][1][0][1] << " " << fGnrt_zvtx[i][j][k][1].tab[c][1][0][1] << " " << fGnrt_zvtx[i][j][k][2].tab[c][1][0][1] << " " << fGnrt_zvtx[i][j][k][3].tab[c][1][0][1] << " " <<
                       fRcstr_zvtx[i][j][k][0].tab[c][1][0][2] << " " << fRcstr_zvtx[i][j][k][1].tab[c][1][0][2] << " " << fRcstr_zvtx[i][j][k][2].tab[c][1][0][2] << " " << fRcstr_zvtx[i][j][k][3].tab[c][1][0][2] << " " <<
                       fGnrt_zvtx[i][j][k][0].tab[c][1][0][2] << " " << fGnrt_zvtx[i][j][k][1].tab[c][1][0][2] << " " << fGnrt_zvtx[i][j][k][2].tab[c][1][0][2] << " " << fGnrt_zvtx[i][j][k][3].tab[c][1][0][2] << " " <<
                       fRcstr_zvtx[i][j][k][0].tab[c][1][0][3] << " " << fRcstr_zvtx[i][j][k][1].tab[c][1][0][3] << " " << fRcstr_zvtx[i][j][k][2].tab[c][1][0][3] << " " << fRcstr_zvtx[i][j][k][3].tab[c][1][0][3] << " " <<
                       fGnrt_zvtx[i][j][k][0].tab[c][1][0][3] << " " << fGnrt_zvtx[i][j][k][1].tab[c][1][0][3] << " " << fGnrt_zvtx[i][j][k][2].tab[c][1][0][3] << " " << fGnrt_zvtx[i][j][k][3].tab[c][1][0][3] << " " <<
                       fRcstr_zvtx[i][j][k][0].tab[c][0][0][0] << " " << fRcstr_zvtx[i][j][k][1].tab[c][0][0][0] << " " << fRcstr_zvtx[i][j][k][2].tab[c][0][0][0] << " " << fRcstr_zvtx[i][j][k][3].tab[c][0][0][0] << " " <<
                       fGnrt_zvtx[i][j][k][0].tab[c][0][0][0] << " " << fGnrt_zvtx[i][j][k][1].tab[c][0][0][0] << " " << fGnrt_zvtx[i][j][k][2].tab[c][0][0][0] << " " << fGnrt_zvtx[i][j][k][3].tab[c][0][0][0] << " " <<
                       fRcstr_zvtx[i][j][k][0].tab[c][0][0][1] << " " << fRcstr_zvtx[i][j][k][1].tab[c][0][0][1] << " " << fRcstr_zvtx[i][j][k][2].tab[c][0][0][1] << " " << fRcstr_zvtx[i][j][k][3].tab[c][0][0][1] << " " <<
                       fGnrt_zvtx[i][j][k][0].tab[c][0][0][1] << " " << fGnrt_zvtx[i][j][k][1].tab[c][0][0][1] << " " << fGnrt_zvtx[i][j][k][2].tab[c][0][0][1] << " " << fGnrt_zvtx[i][j][k][3].tab[c][0][0][1] << " " <<
                       fRcstr_zvtx[i][j][k][0].tab[c][0][0][2] << " " << fRcstr_zvtx[i][j][k][1].tab[c][0][0][2] << " " << fRcstr_zvtx[i][j][k][2].tab[c][0][0][2] << " " << fRcstr_zvtx[i][j][k][3].tab[c][0][0][2] << " " <<
                       fGnrt_zvtx[i][j][k][0].tab[c][0][0][2] << " " << fGnrt_zvtx[i][j][k][1].tab[c][0][0][2] << " " << fGnrt_zvtx[i][j][k][2].tab[c][0][0][2] << " " << fGnrt_zvtx[i][j][k][3].tab[c][0][0][2] << " " <<
                       fRcstr_zvtx[i][j][k][0].tab[c][0][0][3] << " " << fRcstr_zvtx[i][j][k][1].tab[c][0][0][3] << " " << fRcstr_zvtx[i][j][k][2].tab[c][0][0][3] << " " << fRcstr_zvtx[i][j][k][3].tab[c][0][0][3] << " " <<
                       fGnrt_zvtx[i][j][k][0].tab[c][0][0][3] << " " << fGnrt_zvtx[i][j][k][1].tab[c][0][0][3] << " " << fGnrt_zvtx[i][j][k][2].tab[c][0][0][3] << " " << fGnrt_zvtx[i][j][k][3].tab[c][0][0][3] << " " <<  endl;

         for(int zv=0; zv<4; zv++)
         {
           ofs_hmult_zvtx << fRcstr_zvtx[i][j][k][zv].tab[c][1][0][0] << " " << fRcstr_zvtx[i][j][k][zv].tab[c][1][0][0] << " " <<
                             fRcstr_zvtx[i][j][k][zv].tab[c][1][0][1] << " " << fRcstr_zvtx[i][j][k][zv].tab[c][1][0][1] << " " <<
                             fRcstr_zvtx[i][j][k][zv].tab[c][1][0][2] << " " << fRcstr_zvtx[i][j][k][zv].tab[c][1][0][2] << " " <<
                             fRcstr_zvtx[i][j][k][zv].tab[c][1][0][3] << " " << fRcstr_zvtx[i][j][k][zv].tab[c][1][0][3] << " " <<
                             fRcstr_zvtx[i][j][k][zv].tab[c][0][0][0] << " " << fRcstr_zvtx[i][j][k][zv].tab[c][0][0][0] << " " <<
                             fRcstr_zvtx[i][j][k][zv].tab[c][0][0][1] << " " << fRcstr_zvtx[i][j][k][zv].tab[c][0][0][1] << " " <<
                             fRcstr_zvtx[i][j][k][zv].tab[c][0][0][2] << " " << fRcstr_zvtx[i][j][k][zv].tab[c][0][0][2] << " " <<
                             fRcstr_zvtx[i][j][k][zv].tab[c][0][0][3] << " " << fRcstr_zvtx[i][j][k][zv].tab[c][0][0][3] << " " << endl;

           ofs_ezvtx << fRcstr_zvtx[i][j][k][zv].tab[c][1][0][0]/(fRcstr_zvtx[i][j][k][zv].tab[c][1][0][0]+fRcstr_zvtx[i][j][k][zv].tab[c][1][0][4]) << " "
                     << fRcstr_zvtx[i][j][k][zv].tab[c][0][0][0]/(fRcstr_zvtx[i][j][k][zv].tab[c][0][0][0]+fRcstr_zvtx[i][j][k][zv].tab[c][0][0][4]) << endl;
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

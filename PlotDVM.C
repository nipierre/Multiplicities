void PlotDVM()
{
  int x,y,z;
  float dis,had;
  double fDiffVectorMeson[2][9][6][12][4];
  double fDiffVectorMesonErr[2][9][6][12][4];

  Int_t fMarkerColor[2] = {4,2};
  Int_t fMarkerStyle[2] = {24,20};
  Int_t fMarkerColorAlt[5] = {2,95,209,226,221};
  Int_t fMarkerStyleAlt[5][2] = {{24,20},{26,22},{25,21},{27,33},{28,34}};

  TCanvas c1("Pion_DVM","Pion_DVM",3200,1600);
  TCanvas c2("Kaon_DVM","Kaon_DVM",3200,1600);

  double z_range[12] = {.225,.275,.325,.375,.425,.475,.525,.575,.625,.675,.725,.8};

  c1.SetFillColor(0);
  c2.SetFillColor(0);

  c1.Divide(9,5,0,0);
  c2.Divide(9,5,0,0);

  TGraphErrors* P_d[2][9][6];
  TGraphErrors* K_d[2][9][6];

  std::vector<double> p_d[2][9][6];
  std::vector<double> k_d[2][9][6];
  std::vector<double> p_err[2][9][6];
  std::vector<double> k_err[2][9][6];
  std::vector<double> z_range_p[2][9][6];
  std::vector<double> z_range_k[2][9][6];

  ifstream DVM("data/DVM_2016.dat");

  while(DVM >> x)
  {
    DVM >> y >> z;
    DVM >> had >> dis;
    fDiffVectorMeson[1][x-1][y-1][z-1][0] = fDiffVectorMeson[1][x-1][y-1][z-1][3] = had/dis;
    DVM >> had >> dis;
    fDiffVectorMeson[0][x-1][y-1][z-1][0] = fDiffVectorMeson[0][x-1][y-1][z-1][3] = had/dis;
    DVM >> had >> dis;
    fDiffVectorMeson[1][x-1][y-1][z-1][1] = had/dis;
    DVM >> had >> dis;
    fDiffVectorMeson[0][x-1][y-1][z-1][1] = had/dis;
  }
  DVM.close();

  ifstream DVMerr("data/DVMErr.dat");

  while(DVMerr >> x)
  {
    DVMerr >> y >> z;
    DVMerr >> had >> dis;
    DVMerr >> fDiffVectorMesonErr[1][x-1][y-1][z-1][0];
    // cout << fDiffVectorMesonErr[1][x-1][y-1][z-1][0] << " ";
    DVMerr >> had >> dis;
    DVMerr >> fDiffVectorMesonErr[0][x-1][y-1][z-1][0];
    // cout << fDiffVectorMesonErr[0][x-1][y-1][z-1][0] << " ";
    DVMerr >> had >> dis;
    DVMerr >> fDiffVectorMesonErr[1][x-1][y-1][z-1][1];
    // cout << fDiffVectorMesonErr[1][x-1][y-1][z-1][1] << " ";
    DVMerr >> had >> dis;
    DVMerr >> fDiffVectorMesonErr[1][x-1][y-1][z-1][1];
    // cout << fDiffVectorMesonErr[1][x-1][y-1][z-1][1] << " ";
  }
  DVMerr.close();

  for(int i=0; i<9; i++)
  {
    for(int j=0; j<5; j++)
    {
      for(int c=0; c<2; c++)
      {
        for(int k=0; k<12; k++)
        {
          if((j==4 && k>4)
          || (i==0 && j==5 && k==3)
          || (i>4 && j==4 && k==5)
          || (i==0 && j==3 && k==9)
          || (j==3 && k==11)
          || (i==8 && j==4 && k==4)
          || (i==6 && j==3 && k==10)
          || (i==7 && j==3 && k>7)
          || (i==8 && j==3 && k>6)
          || (i>5 && j==5 && k==2)
          || (j==1 && k==3)
          || (i<7 && j==4 && k>5)
          || (i>6 && j==4)
          || (i>5 && j==3 && k>10)
          || (i>7 && j==3 && k>8)
          || (j==2 && k<2)
          || (j==1 && k<4)
          || (j==0 && k<7)
          || (i==8 && j==0)
          || (i==0 && j<3)
          || (i==1 && j==0))
          {
            fDiffVectorMeson[c][i][j][k][0]=0;
            fDiffVectorMeson[c][i][j][k][1]=0;
          }
          p_d[c][i][j].push_back(fDiffVectorMeson[c][i][j][k][0]);
          k_d[c][i][j].push_back(fDiffVectorMeson[c][i][j][k][1]);
          p_err[c][i][j].push_back(sqrt(fDiffVectorMesonErr[c][i][j][k][0]));
          k_err[c][i][j].push_back(sqrt(fDiffVectorMesonErr[c][i][j][k][1]));

          z_range_p[c][i][j].push_back(z_range[k]);
          z_range_k[c][i][j].push_back(z_range[k]);
        }
        for(int k=12; k>0; k--)
        {
          if(p_d[c][i][j][k-1]==0) {p_d[c][i][j].erase(p_d[c][i][j].begin()+k-1); p_err[c][i][j].erase(p_err[c][i][j].begin()+k-1); /*p_sys[c][i][j].erase(p_sys[c][i][j].begin()+k-1);*/ z_range_p[c][i][j].erase(z_range_p[c][i][j].begin()+k-1);}
          if(k_d[c][i][j][k-1]==0) {k_d[c][i][j].erase(k_d[c][i][j].begin()+k-1); k_err[c][i][j].erase(k_err[c][i][j].begin()+k-1); /*k_sys[c][i][j].erase(k_sys[c][i][j].begin()+k-1);*/ z_range_k[c][i][j].erase(z_range_k[c][i][j].begin()+k-1);}
        }

        bool p_empty = 0;
        bool k_empty = 0;

        if(!(Int_t(p_d[c][i][j].size()))) p_empty = 1;
        if(!(Int_t(k_d[c][i][j].size()))) k_empty = 1;

        P_d[c][i][j] = new TGraphErrors(Int_t(p_d[c][i][j].size()),&(z_range_p[c][i][j][0]),&(p_d[c][i][j][0]),0,&(p_err[c][i][j][0]));
        K_d[c][i][j] = new TGraphErrors(Int_t(k_d[c][i][j].size()),&(z_range_k[c][i][j][0]),&(k_d[c][i][j][0]),0,&(k_err[c][i][j][0]));
        // P_sys[c][i][j] = new TGraphAsymmErrors(Int_t(p_m[c][i][j].size()),&(z_range_p[c][i][j][0]), &p_yoffset[0], &errorx[0], &errorx[0], 0, &(p_sys[c][i][j][0]));
        // K_sys[c][i][j] = new TGraphAsymmErrors(Int_t(k_m[c][i][j].size()),&(z_range_k[c][i][j][0]), &k_yoffset[0], &errorx[0], &errorx[0], 0, &(k_sys[c][i][j][0]));

        P_d[c][i][j]->SetMarkerColor(fMarkerColor[c]);
        K_d[c][i][j]->SetMarkerColor(fMarkerColor[c]);

        P_d[c][i][j]->SetMarkerSize(1);
        K_d[c][i][j]->SetMarkerSize(1);

        P_d[c][i][j]->SetMarkerStyle(fMarkerStyle[c]);
        K_d[c][i][j]->SetMarkerStyle(fMarkerStyle[c]);

        P_d[c][i][j]->SetTitle("");
        K_d[c][i][j]->SetTitle("");

        P_d[c][i][j]->GetYaxis()->SetTitle("");
        K_d[c][i][j]->GetYaxis()->SetTitle("");

        // P_sys[c][i][j]->SetFillColor(fMarkerColor[j]);
        // K_sys[c][i][j]->SetFillColor(fMarkerColor[j]);

        if(!p_empty)
        {
          c1.cd(i+1+9*j);
          if(P_d[c][i][j])
          {
            if(!c)
            {
              P_d[c][i][j]->Draw("SAMEPA");
              P_d[c][i][j]->GetXaxis()->SetLimits(0.1,0.9);
              P_d[c][i][j]->SetMinimum(0.6);
              P_d[c][i][j]->SetMaximum(1.15);
              P_d[c][i][j]->GetXaxis()->SetLabelSize(0.06);
              P_d[c][i][j]->GetYaxis()->SetLabelSize(0.06);
              P_d[c][i][j]->SetTitle("");
              if(j==5) gPad->SetBottomMargin(.15);
              if(i==0) gPad->SetLeftMargin(.22);
              if(i==8 && j==5)
              {
                P_d[c][i][j]->GetXaxis()->SetTitle("#font[ 12]{z}");
                P_d[c][i][j]->GetXaxis()->SetTitleSize(0.08);
                P_d[c][i][j]->GetXaxis()->SetTitleOffset(.8);
              }
              P_d[c][i][j]->GetXaxis()->SetNdivisions(304,kTRUE);
              P_d[c][i][j]->GetYaxis()->SetNdivisions(304,kTRUE);
              if(i==2 && j==0)
              {
                P_d[c][i][j]->GetYaxis()->SetTitle("#font[12]{acceptance}^{#font[ 12]{#pi}}");
                P_d[c][i][j]->GetYaxis()->SetTitleSize(0.08);
              }
              P_d[c][i][j]->Draw("SAMEP");
              P_d[c][i][j]->GetXaxis()->SetLimits(0.1,0.9);
              P_d[c][i][j]->SetMinimum(0.6);
              P_d[c][i][j]->SetMaximum(1.15);
              c1.Range(0.,0.,1.,1.2);
            }
            else
            {
              P_d[c][i][j]->Draw("SAMEP");
              P_d[c][i][j]->GetXaxis()->SetLimits(0.1,0.9);
              P_d[c][i][j]->SetMinimum(0.6);
              P_d[c][i][j]->SetMaximum(1.15);
            }
          }
          c1.Update();
        }

        if(!k_empty)
        {
          c2.cd(i+1+9*j);
          if(K_d[c][i][j])
          {
            if(!c)
            {
              K_d[c][i][j]->Draw("SAMEPA");
              K_d[c][i][j]->GetXaxis()->SetLimits(0.1,0.9);
              K_d[c][i][j]->SetMinimum(0.7);
              K_d[c][i][j]->SetMaximum(1.15);
              K_d[c][i][j]->GetXaxis()->SetLabelSize(0.06);
              K_d[c][i][j]->GetYaxis()->SetLabelSize(0.06);
              K_d[c][i][j]->SetTitle("");
              if(j==5) gPad->SetBottomMargin(.15);
              if(i==0) gPad->SetLeftMargin(.22);
              if(i==8 && j==5)
              {
                K_d[c][i][j]->GetXaxis()->SetTitle("#font[ 12]{z}");
                K_d[c][i][j]->GetXaxis()->SetTitleSize(0.08);
                K_d[c][i][j]->GetXaxis()->SetTitleOffset(.8);
              }
              K_d[c][i][j]->GetXaxis()->SetNdivisions(304,kTRUE);
              K_d[c][i][j]->GetYaxis()->SetNdivisions(304,kTRUE);
              if(i==2 && j==0)
              {
                K_d[c][i][j]->GetYaxis()->SetTitle("#font[12]{acceptance}^{#font[ 12]{K}}");
                K_d[c][i][j]->GetYaxis()->SetTitleSize(0.08);
              }
              K_d[c][i][j]->Draw("SAMEP");
              K_d[c][i][j]->GetXaxis()->SetLimits(0.1,0.9);
              K_d[c][i][j]->SetMinimum(0.7);
              K_d[c][i][j]->SetMaximum(1.15);
              c2.Range(0.,0.,1.,1.2);
            }
            else
            {
              K_d[c][i][j]->Draw("SAMEP");
              K_d[c][i][j]->GetXaxis()->SetLimits(0.1,0.9);
              K_d[c][i][j]->SetMinimum(0.7);
              K_d[c][i][j]->SetMaximum(1.15);
            }
          }
          c2.Update();
        }
      }
    }
  }

  c1.Print("DVM_Pi.pdf");
  c2.Print("DVM_K.pdf");

}

void PlotAccComp()
{
  float dummy;
  double fAcc1[2][9][6][12][2];
  double fAcc1Err[2][9][6][12][2];
  double fAcc2[2][9][6][12][2];
  double fAcc2Err[2][9][6][12][2];
  double fComp[2][9][6][12][2];
  double fCompErr[2][9][6][12][2];

  Int_t fMarkerColor[6] = {2,95,209,226,4,221};
  Int_t fMarkerStyle[6][2] = {{24,20},{26,22},{25,21},{27,33},{28,34},{30,29}};
  Int_t fMarkerStyleb[2][2] = {{24,20},{27,33}};

  TCanvas c1("Comp","Comp",3200,1600);

  double z_range[12] = {.225,.275,.325,.375,.425,.475,.525,.575,.625,.675,.725,.8};

  c1.SetFillColor(0);

  c1.Divide(9,5,0,0);

  TGraphErrors* H_acc[2][9][6][2];

  ifstream Acc1("acceptance_P07.txt");
  for(int c=0; c<2; c++)
  {
    for(int i=0; i<9; i++)
    {
      for(int j=0; j<6; j++)
      {
        for(int k=0; k<12; k++)
        {
          for(int l=0; l<10; l++)
          {
            Acc1 >> dummy;
          }
          Acc1 >> fAcc1[c][i][j][k][1] >> fAcc1Err[c][i][j][k][1];
          for(int l=0; l<6; l++)
          {
            Acc1 >> dummy;
          }
          Acc1 >> fAcc1[c][i][j][k][0] >> fAcc1Err[c][i][j][k][0];
        }
      }
    }
  }
  Acc1.close();

  ifstream Acc2("acceptance_P09.txt");
  for(int c=0; c<2; c++)
  {
    for(int i=0; i<9; i++)
    {
      for(int j=0; j<6; j++)
      {
        for(int k=0; k<12; k++)
        {
          for(int l=0; l<10; l++)
          {
            Acc2 >> dummy;
          }
          Acc2 >> fAcc2[c][i][j][k][1] >> fAcc2Err[c][i][j][k][1];
          for(int l=0; l<6; l++)
          {
            Acc2 >> dummy;
          }
          Acc2 >> fAcc2[c][i][j][k][0] >> fAcc2Err[c][i][j][k][0];
          fComp[c][i][j][k][1] = fAcc1[c][i][j][k][1]/fAcc2[c][i][j][k][1];
          fComp[c][i][j][k][0] = fAcc1[c][i][j][k][0]/fAcc2[c][i][j][k][0];
          fCompErr[c][i][j][k][1] = fAcc1Err[c][i][j][k][1]+fAcc2Err[c][i][j][k][1];
          fCompErr[c][i][j][k][0] = fAcc1Err[c][i][j][k][0]+fAcc2Err[c][i][j][k][0];
        }
      }
    }
  }
  Acc2.close();

  TLine l1(0.1,1.05,0.9,1.05);
  TLine l2(0.1,0.95,0.9,0.95);

  l1.SetLineStyle(3); l2.SetLineStyle(3);

  for(int c=0; c<2; c++)
  {
    for(int i=0; i<9; i++)
    {
      for(int m=0; m<2; m++)
      {
        for(int j=0; j<5; j++)
        {
          std::vector<double> h_a;
          std::vector<double> h_err;
          std::vector<double> z_range_h;

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
              fComp[c][i][j][k][0]=0;
              fComp[c][i][j][k][1]=0;
            }
            h_a.push_back(fComp[c][i][j][k][m]);
            h_err.push_back(sqrt(fCompErr[c][i][j][k][m]));

            z_range_h.push_back(z_range[k]);
          }

          for(int k=12; k>0; k--)
          {
            if(!h_a[k-1]) {h_a.erase(h_a.begin()+k-1); h_err.erase(h_err.begin()+k-1); z_range_h.erase(z_range_h.begin()+k-1);}
          }

          bool h_a_empty = 0;

          if(!(int(h_a.size()))) h_a_empty = 1;

          H_acc[c][i][j][m] = new TGraphErrors(int(h_a.size()),&(z_range_h[0]),&(h_a[0]),0,&(h_err[0]));

          if(!c)
          {
            H_acc[c][i][j][m]->SetMarkerColor(fMarkerColor[4]);
          }
          else
          {
            H_acc[c][i][j][m]->SetMarkerColor(fMarkerColor[0]);
          }

          H_acc[c][i][j][m]->SetMarkerSize(1);

          if(!m)
          {
            H_acc[c][i][j][m]->SetMarkerStyle(fMarkerStyleb[1][c]);
          }
          else
          {
            H_acc[c][i][j][m]->SetMarkerStyle(fMarkerStyleb[0][c]);
          }

          H_acc[c][i][j][m]->GetYaxis()->SetTitle("");

          H_acc[c][i][j][m]->GetXaxis()->SetTitle("");

          H_acc[c][i][j][m]->SetTitle("");

          if(!h_a_empty)
          {
            c1.cd(i+1+9*j);
            gPad->SetFillStyle(4000);
            if(H_acc[c][i][j][m])
            {
              if(!c && !m)
              {
                H_acc[c][i][j][m]->Draw("SAMEPA");
                H_acc[c][i][j][m]->GetXaxis()->SetLimits(0.1,0.9);
                H_acc[c][i][j][m]->SetMinimum(0.75);
                H_acc[c][i][j][m]->SetMaximum(1.25);
                H_acc[c][i][j][m]->GetXaxis()->SetLabelSize(0.06);
                H_acc[c][i][j][m]->GetYaxis()->SetLabelSize(0.06);
                H_acc[c][i][j][m]->SetTitle("");
                if(j==5) gPad->SetBottomMargin(.15);
                if(i==0) gPad->SetLeftMargin(.22);
                if(i==8 && j==5)
                {
                  H_acc[c][i][j][m]->GetXaxis()->SetTitle("#font[ 12]{z}");
                  H_acc[c][i][j][m]->GetXaxis()->SetTitleSize(0.08);
                  H_acc[c][i][j][m]->GetXaxis()->SetTitleOffset(.8);
                }
                H_acc[c][i][j][m]->GetXaxis()->SetNdivisions(304,kTRUE);
                H_acc[c][i][j][m]->GetYaxis()->SetNdivisions(304,kTRUE);
                if(i==1 && j==0)
                {
                  H_acc[c][i][j][m]->GetYaxis()->SetTitle("#font[12]{acceptance}^{#font[ 12]{h}}");
                  H_acc[c][i][j][m]->GetYaxis()->SetTitleSize(0.08);
                }
                H_acc[c][i][j][m]->Draw("SAMEP");
                l1.Draw("SAME"); l2.Draw("SAME");
                H_acc[c][i][j][m]->GetXaxis()->SetLimits(0.1,0.9);
                H_acc[c][i][j][m]->SetMinimum(0.75);
                H_acc[c][i][j][m]->SetMaximum(1.25);
                c1.Range(0.1,0.,0.9,1.2);
              }
              else
              {
                H_acc[c][i][j][m]->Draw("SAMEP");
                l1.Draw("SAME"); l2.Draw("SAME");
                H_acc[c][i][j][m]->GetXaxis()->SetLimits(0.1,0.9);
                H_acc[c][i][j][m]->SetMinimum(0.75);
                H_acc[c][i][j][m]->SetMaximum(1.25);
              }
            }
            c1.Update();
          }
        }
      }
    }
  }

  c1.Print("Acceptance_comparison.pdf");
}

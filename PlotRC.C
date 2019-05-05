void PlotRC()
{
  Double_t fMult[2][9][6][12];
  Double_t fMult_err[2][9][6][12];

  int fMarkerColor[2] = {4,2};
  int fMarkerStyle[2][2] = {{24,20},{26,22}};
  int fMarkerColorAlt[6] = {2,95,209,226,4,221};
  int fMarkerStyleAlt[6][2] = {{24,20},{26,22},{25,21},{27,33},{28,34},{30,29}};

  TGraphErrors* H_d[2][9][6];

  string sdum;

  ifstream proton("data/newRC.txt");

  for(int c=0; c<2; c++)
  {
    for(int i=0; i<9; i++)
    {
      for(int j=0; j<6; j++)
      {
        for(int l=0; l<9; l++)
        {
          proton >> sdum;
        }
        for(int k=0; k<12; k++)
        {
          for(int l=0; l<7; l++)
          {
            proton >> sdum;
          }
          proton >> fMult[c][i][j][k] >> fMult_err[c][i][j][k];
          cout << fMult[c][i][j][k] << " " << fMult_err[c][i][j][k] << endl;
        }
        for(int l=0; l<9; l++)
        {
          proton >> sdum;
        }
      }
    }
  }
  proton.close();

  TCanvas c1("Hadron_Correction","Hadron_Correction",3200,1600);

  c1.SetFillColor(0);
  c1.Divide(5,2,0,0);

  TLine *l1[5];
  TLine *l2[5];

  for(int i=0; i<5; i++)
  {
    l1[i] = new TLine(0.175,1+i*0.1,0.85,1+i*0.1);
    l1[i]->SetLineStyle(4);
    l1[i]->SetLineColor(fMarkerColorAlt[i]);
    l2[i] = new TLine(0.,1+i*0.1,6.,1+i*0.1);
    l2[i]->SetLineStyle(4);
    l2[i]->SetLineColor(fMarkerColorAlt[i]);
  }

  double z_range[12] = {.225,.275,.325,.375,.425,.475,.525,.575,.625,.675,.725,.8};

  for(int c=0; c<2; c++)
  {
    for(int i=0; i<9; i++)
    {
      for(int j=0; j<5; j++)
      {
        std::vector<double> h_a;
        std::vector<double> h_err;
        std::vector<double> z_range_h;
        for(int l=0; l<12; l++)
        {
          z_range_h.push_back(z_range[l]);
        }
        for(int k=0; k<12; k++)
        {
          if((i==0 && j==4 && k>4)
            ||(i==0 && j==3 && k>9)
            ||(i==1 && j==4 && k>4)
            ||(i==1 && j==2 && k<2)
            ||(i==1 && j==1 && k<4)
            ||(i==2 && j==4 && k>4)
            ||(i==2 && j==2 && k<2)
            ||(i==2 && j==1 && k<4)
            ||(i==2 && j==0 && k<6)
            ||(i==3 && j==4 && k>4)
            ||(i==3 && j==2 && k<2)
            ||(i==3 && j==1 && k<4)
            ||(i==3 && j==0 && k<6)
            ||(i==4 && j==4 && k>4)
            ||(i==4 && j==2 && k<2)
            ||(i==4 && j==1 && k<4)
            ||(i==4 && j==0 && k<6)
            ||(i==5 && j==4 && k>4)
            ||(i==5 && j==2 && k<2)
            ||(i==5 && j==1 && k<4)
            ||(i==5 && j==0 && k<6)
            ||(i==6 && j==4 && k>4)
            ||(i==6 && j==3 && k>10)
            ||(i==6 && j==2 && k<2)
            ||(i==6 && j==1 && k<4)
            ||(i==6 && j==0 && k<6)
            ||(i==7 && j==4)
            ||(i==7 && j==3 && k>10)
            ||(i==7 && j==2 && k<2)
            ||(i==7 && j==1 && k<4)
            ||(i==7 && j==0 && k<6)
            ||(i==8 && j==4)
            ||(i==8 && j==3 && k>8)
            ||(i==8 && j==2 && k<2)
            ||(i==8 && j==1 && k<4)
            ||(i==8 && j==0)
          )
          {
            fMult[c][i][j][k]=0;
            fMult_err[c][i][j][k]=0;
          }
          h_a.push_back(fMult[c][i][j][k] ? fMult[c][i][j][k]+0.1*j : 0);
          h_err.push_back(fMult_err[c][i][j][k]);
        }

        for(int k=12; k>0; k--)
        {
          if(!h_a[k-1]) {h_a.erase(h_a.begin()+k-1); h_err.erase(h_err.begin()+k-1); z_range_h.erase(z_range_h.begin()+k-1);}
        }

        bool h_a_empty = 0;
        if(!(int(h_a.size()))) h_a_empty = 1;

        H_d[c][i][j] = new TGraphErrors(int(h_a.size()),&(z_range_h[0]),&(h_a[0]),0,&(h_err[0]));
        H_d[c][i][j]->SetMarkerColor(fMarkerColorAlt[j]);
        H_d[c][i][j]->SetMarkerSize(2);
        H_d[c][i][j]->SetMarkerStyle(fMarkerStyleAlt[j][c]);
        H_d[c][i][j]->GetYaxis()->SetTitle("");
        H_d[c][i][j]->GetXaxis()->SetTitle("");
        H_d[c][i][j]->SetTitle("");

        if(!h_a_empty)
        {
          c1.cd(i+1);
          gPad->SetFillStyle(4000);
          if(H_d[c][i][j])
          {
            if(c && j==3)
            {
              H_d[c][i][j]->Draw("SAMEPA");
              H_d[c][i][j]->GetXaxis()->SetLimits(0.175,0.85);
              H_d[c][i][j]->SetMinimum(0.98);
              H_d[c][i][j]->SetMaximum(1.65);
              H_d[c][i][j]->GetXaxis()->SetLabelSize(0.06);
              H_d[c][i][j]->GetYaxis()->SetLabelSize(0.06);
              H_d[c][i][j]->SetTitle("");
              if(i>4) gPad->SetBottomMargin(.15);
              if(i==0 || i==5) gPad->SetLeftMargin(.22);
              if(i==8)
              {
                H_d[c][i][j]->GetXaxis()->SetTitle("#font[ 12]{z}");
                H_d[c][i][j]->GetXaxis()->SetTitleSize(0.08);
                H_d[c][i][j]->GetXaxis()->SetTitleOffset(.8);
              }
              H_d[c][i][j]->GetXaxis()->SetNdivisions(304,kTRUE);
              H_d[c][i][j]->GetYaxis()->SetNdivisions(304,kTRUE);
              if(i==0)
              {
                H_d[c][i][j]->GetYaxis()->SetTitle("#font[12]{#eta}^{#font[ 12]{h}}+ #font[ 12]{#delta}");
                H_d[c][i][j]->GetYaxis()->SetTitleSize(0.08);
              }
              H_d[c][i][0]->Draw("SAMEP");
              H_d[c][i][0]->GetXaxis()->SetLimits(0.175,0.85);
              H_d[c][i][0]->SetMinimum(0.98);
              H_d[c][i][0]->SetMaximum(1.65);
              H_d[c][i][1]->Draw("SAMEP");
              H_d[c][i][1]->GetXaxis()->SetLimits(0.175,0.85);
              H_d[c][i][1]->SetMinimum(0.98);
              H_d[c][i][1]->SetMaximum(1.65);
              H_d[c][i][2]->Draw("SAMEP");
              H_d[c][i][2]->GetXaxis()->SetLimits(0.175,0.85);
              H_d[c][i][2]->SetMinimum(0.98);
              H_d[c][i][2]->SetMaximum(1.65);
              c1.Range(0.,.85,1.,1.65);
              l1[0]->Draw();
              l1[1]->Draw();
              l1[2]->Draw();
              l1[3]->Draw();
              l1[4]->Draw();
            }
            else
            {
              H_d[c][i][j]->Draw("SAMEP");
              H_d[c][i][j]->GetXaxis()->SetLimits(0.175,0.85);
              H_d[c][i][j]->SetMinimum(0.98);
              H_d[c][i][j]->SetMaximum(1.65);
            }
          }
          c1.Update();

        }
      }
    }
  }

  TLatex fTitle;

  c1.cd(1);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.6,"0.004#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.01");

  c1.cd(2);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.6,"0.01#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.02");

  c1.cd(3);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.6,"0.02#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.03");

  c1.cd(4);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.6,"0.03#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.04");

  c1.cd(5);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.6,"0.04#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.06");

  c1.cd(6);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.6,"0.06#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.1");

  c1.cd(7);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.6,"0.1#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.14");

  c1.cd(8);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.6,"0.14#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.18");

  c1.cd(9);
  fTitle.SetTextSize(0.078);
  fTitle.SetTextAlign(21);
  fTitle.DrawLatex(0.5, 1.6,"0.18#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.4");

  c1.cd(10);
  fTitle.SetTextSize(0.095);
  fTitle.SetTextAlign(11);
  fTitle.DrawLatex(0.05, 0.64,"#color[4]{0.50#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.70, #delta = 0.4}");
  fTitle.DrawLatex(0.05, 0.56,"#color[226]{0.30#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.50, #delta = 0.3}");
  fTitle.DrawLatex(0.05, 0.48,"#color[209]{0.20#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.30, #delta = 0.2}");
  fTitle.DrawLatex(0.05, 0.40,"#color[95]{0.15#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.20, #delta = 0.1}");
  fTitle.DrawLatex(0.05, 0.32,"#color[2]{0.10#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{y}#scale[0.5]{ }<#scale[0.5]{ }0.15, #delta = 0}");

  c1.Update();
  c1.Print("hadron_correction.pdf");


}

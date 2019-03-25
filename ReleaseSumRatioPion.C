void plot_sum_average_ratio()
{
  double zwidth[12]= {0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.1};
  double zwidthH[4]= {0.1,0.1,0.2,0.2};
  ////////////read in txt files
  double m_u[2][9][6][12];double me_u[2][9][6][12];
  double sys_a[2][9][6][12];
  double z[2][9][6][12];
  double x[2][9][6][12];
  for(int jj=0; jj<9; jj++)
    for(int kk=0; kk<6; kk++)
      for(int ll=0; ll<12; ll++)
	for(int ii=0; ii<2; ii++)
	  {
	    m_u[ii][jj][kk][ll]=0 ;
	    me_u[ii][jj][kk][ll]=0;
	    sys_a[ii][jj][kk][ll]=0.;
	    z[ii][jj][kk][ll]=0.;
	    x[ii][jj][kk][ll]=0.;
	  }//jj
  ///////////////////read in multiplicity
  int ii=-1; int jj=-1; int kk=-1;
  double jjj=0; double kkk=0; double iii=0;
  double xa_totalp=0; double ya_totalp=0; double q2a_totalp=0;double z_kaonp=0;
  double mp=0;double mep=0;double sysp=0;double leptoflagp=0;
  double xa_totaln=0;double ya_totaln=0;double q2a_totaln=0;double z_kaonn=0;
  double mn=0;double men=0;double sysn=0;double leptoflagn=0;

  ifstream IN1;
  IN1.open("./data/MultiplicityKaon_release_MCmeans_newnewRC_halfrad.txt", ifstream::in);
  for(int j=0; j<40000; j++)
    {
      if (! IN1.eof() )
	{
	  IN1>>iii>>jjj>>kkk>>xa_totalp>>ya_totalp>>q2a_totalp>>z_kaonp>>mp>>mep>>sysp>>leptoflagp>>xa_totaln>>ya_totaln>>q2a_totaln>>z_kaonn>>mn>>men>>sysn>>leptoflagn;
	  //
	  if(iii>0)ii=0;if(iii>0.0045)ii=1;if(iii>0.015)ii=2;if(iii>0.025)ii=3;if(iii>0.035)ii=4;if(iii>0.045)ii=5;if(iii>0.065)ii=6;if(iii>0.13)ii=7;if(iii>0.16)ii=8;
	  if(jjj>0)jj=0;if(jjj>0.105)jj=1;if(jjj>0.155)jj=2;if(jjj>0.205)jj=3;if(jjj>0.305)jj=4;if(jjj>0.505)jj=5;
	  if(kkk>0)kk=0;if(kkk>0.21)kk=1;if(kkk>0.251)kk=2;if(kkk>0.31)kk=3;if(kkk>0.351)kk=4;if(kkk>0.41)kk=5;if(kkk>0.451)kk=6;if(kkk>0.51)kk=7;if(kkk>0.551)kk=8;if(kkk>0.61)kk=9;if(kkk>0.651)kk=10;if(kkk>0.71)kk=11;
	  if(ii>-1&&jj>-1&&kk>-1){
	    m_u[1][ii][jj][kk]=mp;
	    me_u[1][ii][jj][kk]=mep;
	    sys_a[1][ii][jj][kk]=sysp;
	    z[1][ii][jj][kk]=z_kaonp;
	    x[1][ii][jj][kk]=xa_totalp;
	    //
	    m_u[0][ii][jj][kk]=mn;
	    me_u[0][ii][jj][kk]=men;
	    sys_a[0][ii][jj][kk]=sysn;
	    z[0][ii][jj][kk]=z_kaonn;
	    x[0][ii][jj][kk]=xa_totaln;
	    if((jj==0&&ii==8)||(jj==4&&(ii==7||ii==8))){
	      m_u[1][ii][jj][kk]=0;
	      me_u[1][ii][jj][kk]=9999999999999;
	      sys_a[1][ii][jj][kk]=0;
	      //
	      m_u[0][ii][jj][kk]=0;
	      me_u[0][ii][jj][kk]=9999999999999;
	      sys_a[0][ii][jj][kk]=0;
	    }
	    //
	    if(m_u[0][ii][jj][kk]==0)z[0][ii][jj][kk]=-10;
	    if(m_u[0][ii][jj][kk]==0)m_u[0][ii][jj][kk]=-10;
	    if(m_u[0][ii][jj][kk]==0)me_u[0][ii][jj][kk]=0;

	    if(m_u[1][ii][jj][kk]==0)z[1][ii][jj][kk]=-10;
	    if(m_u[1][ii][jj][kk]==0)m_u[1][ii][jj][kk]=-10;
	    if(m_u[1][ii][jj][kk]==0)me_u[1][ii][jj][kk]=-0;
	    jjj=kkk=iii=xa_totalp= ya_totalp= q2a_totalp= z_kaonp = mp= mep= sysp= leptoflagp=xa_totaln= ya_totaln= q2a_totaln= z_kaonn=mn= men= sysn= leptoflagn=0;
	  }
	}
    }
  IN1.close();
  //
  cout<<"done reading COMPASS kaon"<<endl;




  /////////////////read in q2 xBs///////////for a better average (has counts and x-values)

  double x_total[4][9][6][12];
  double x_counts[4][9][6][12];
  double q_total[4][9][6][12];
  double q_counts[4][9][6][12];
  double y_total[4][9][6][12];
  double y_counts[4][9][6][12];

  double qq;  double nq;  double xx;  double nx;  double yy;  double ny; double dummy;
  for(int hh=0; hh<4; hh++)
    for(int jj=0; jj<9; jj++)
      for(int kk=0; kk<6; kk++)
	for(int ll=0; ll<12; ll++)
	  {
	    x_total[hh][jj][kk][ll]=0;
	    x_counts[hh][jj][kk][ll]=0;
	    q_total[hh][jj][kk][ll]=0;
	    q_counts[hh][jj][kk][ll]=0;
	    y_total[hh][jj][kk][ll]=0;
	    y_counts[hh][jj][kk][ll]=0;
	  }
  hh=ii=jj=kk=-1;    qq=nq =xx =nx=yy =ny=0;
  IN1.open("./inputs/Q2_xB_y_for_allbins_MC.txt", ifstream::in);
  for(int j=0; j<40000; j++)
    {
      if (! IN1.eof() )
	{
	  IN1>>hh>>ii>>jj>>ll>>qq>>nq>>xx>>nx>>yy>>ny>>dummy>>dummy>>dummy>>dummy;
	  if(hh>-1&&ii>-1&&jj>-1&&ii<9&&ll>-1){
	    if(xx>0)x_total[hh][ii][jj][ll]=nx/xx;
	    if(xx>0)x_counts[hh][ii][jj][ll]=xx;
	    if(qq>0)q_total[hh][ii][jj][ll]=nq/qq;
	    if(qq>0)q_counts[hh][ii][jj][ll]=qq;
	    if(yy>0)y_total[hh][ii][jj][ll]=ny/yy;
	    if(yy>0)y_counts[hh][ii][jj][ll]=yy;
	    cout<<"x : "<<x_total[hh][ii][jj][ll]<<endl;
	    hh=ii=jj=kk=-1;    qq=nq =xx =nx=yy =ny=0;
	  }
	}
    }
  IN1.close();

 double zerror[12] ={0.};
/////////////////    AVERAGE   (done with statistical weights) ///////////////////////////

 double avg_sys[2][2][9][12];//correlated = avg_sys[1], uncorelated = avg_sys[0]
 double avg_sys_plotb[2][2][9][12]; //this is to make the white filling on K- systematic band on the average plot
 double avg_totsys[2][2][9][12]; //this is to make the systematic bands on the average plot

 double avg_m[2][2][9][12];//average multiplicity
 double avg_em[2][2][9][12];//avergae multiplicity statistical error
 double avg_x[9]; double temp_count[9];//average x for sum and ratio plots
 double avg_q[9]; double temp_countq[9];//average q2
 double avg_z[2][2][9][12]; //average z for the average plot
 double temp_countz[2][2][9][12];
 for(int dd=0;dd<2;dd++)
   for(int pp=0;pp<2;pp++)
     for(int ii=0;ii<9;ii++)
       for(int jj=0;jj<12;jj++){
	 avg_z[pp][dd][ii][jj]=0;
	 temp_countz[pp][dd][ii][jj]=0;
	 avg_m[pp][dd][ii][jj] =0;
	 avg_em[pp][dd][ii][jj] =0;
	 avg_sys[pp][dd][ii][jj] =0;
       }

     TGraphErrors *dmult_avg[2][2][9];

     TGraphErrors *h_ratio_2d[9];
     TGraphAsymmErrors *R00[9];
     double avg_eratio[9][12]; double avg_ratio[9][12]; double avg_sratio[9][12];
     //for systematic bands:
     double errorx[12]= {0.05/2,0.05/2,0.05/2,0.05/2,0.05/2,0.05/2,0.05/2,0.05/2,0.05/2,0.05/2,0.05/2,0.1/2};
     double z_avg[12]= {0.2+0.05/2,0.25+0.05/2,0.3+0.05/2,0.35+0.05/2,0.4+0.05/2,0.45+0.05/2,0.5+0.05/2,0.55+0.05/2,0.6+0.05/2,0.65+0.05/2,0.7+0.05/2,0.75+0.1/2};
     double z_avgb[12]= {0.2+0.05/2+0.005,0.25+0.05/2,0.3+0.05/2,0.35+0.05/2,0.4+0.05/2,0.45+0.05/2,0.5+0.05/2,0.55+0.05/2,0.6+0.05/2,0.65+0.05/2,0.7+0.05/2,0.75+0.1/2};
     double ysys[12] = {0};
     double yoffset[4][12];
     double yoffsetb[2][12];
     for(int ll=0; ll<12; ll++)
       {yoffset[1][ll] = -0.2;
	 yoffset[0][ll] = -0.35;

	 yoffsetb[0][ll] = -0.34;
	 yoffset[2][ll] = 0.65;
	 yoffset[3][ll] = 0.1;

	 yoffset[1][ll] = -0.035;
	 yoffset[0][ll] = -0.061;
	 yoffsetb[0][ll] = -0.0595;
       }
     TGraphAsymmErrors *A00[2][2][9];
     TGraphAsymmErrors *A00b[2][2][9];

     for(int jj=0; jj<9; jj++)
       {
	 for(int ii=0; ii<2; ii++)
	   {
	     for(int ll=0; ll<12; ll++)
	       {
		 for(int kk=0; kk<5; kk++)
		   {
		     if(((jj==0||jj==4)&&kk==8)||(jj==4&&kk==7));else{ //remove badbins
		       if(m_u[ii][jj][kk][ll]>0&&me_u[ii][jj][kk][ll]>0){
			 avg_m[1][ii][jj][ll] +=  (m_u[ii][jj][kk][ll])/pow(me_u[ii][jj][kk][ll],2);
			 avg_z[1][ii][jj][ll] +=  (z[ii][jj][kk][ll])/pow(me_u[ii][jj][kk][ll],2);
			 if(avg_m[1][ii][jj][ll]>0)avg_sys[1][ii][jj][ll] += (pow(0.8*sys_a[ii][jj][kk][ll],2))/pow(me_u[ii][jj][kk][ll],4);//fully correlated
			 if(avg_m[1][ii][jj][ll]>0)avg_sys[0][ii][jj][ll] += 0.6*sys_a[ii][jj][kk][ll]/pow(me_u[ii][jj][kk][ll],2);//fully uncorrelated
			 avg_em[1][ii][jj][ll] += (pow(me_u[ii][jj][kk][ll],-2));
		       }
		     }//remove badbins
		   }//kk
		 if(avg_m[1][ii][jj][ll]>0){
		   avg_z[1][ii][jj][ll] = avg_z[1][ii][jj][ll]/avg_em[1][ii][jj][ll];
		   avg_m[1][ii][jj][ll] =avg_m[1][ii][jj][ll]/avg_em[1][ii][jj][ll];
		   avg_sys[1][ii][jj][ll] = sqrt(avg_sys[1][ii][jj][ll]/pow(avg_em[1][ii][jj][ll],2));//fully correlated
		   avg_sys[0][ii][jj][ll] = avg_sys[0][ii][jj][ll]/avg_em[1][ii][jj][ll];//fully uncorrelated
		   avg_totsys[1][ii][jj][ll] = sqrt(pow(avg_sys[0][ii][jj][ll],2)+pow(avg_sys[1][ii][jj][ll],2));
		   avg_sys_plotb[1][ii][jj][ll] = avg_totsys[1][ii][jj][ll]-0.0032;
		   avg_em[1][ii][jj][ll] = 1/sqrt(avg_em[1][ii][jj][ll]);}
		 if(avg_m[0][ii][jj][ll]>0){
		   avg_m[0][ii][jj][ll] =avg_m[0][ii][jj][ll]/avg_em[0][ii][jj][ll];
		   avg_em[0][ii][jj][ll] = 1/sqrt(avg_em[0][ii][jj][ll]);}

		 if(ii==1) {
		   avg_eratio[jj][ll]=avg_sratio[jj][ll]=avg_ratio[jj][ll]=0;
		   avg_ratio[jj][ll] = avg_m[1][1][jj][ll]/avg_m[1][0][jj][ll];
		   avg_eratio[jj][ll] = avg_m[1][1][jj][ll]/avg_m[1][0][jj][ll]*sqrt(pow(avg_em[1][0][jj][ll]/avg_m[1][0][jj][ll],2)+pow(avg_em[1][1][jj][ll]/avg_m[1][1][jj][ll],2));
		   avg_sratio[jj][ll] = sqrt( pow(((avg_sys[1][1][jj][ll]+avg_m[1][1][jj][ll])/(avg_sys[1][0][jj][ll]+avg_m[1][0][jj][ll])-(avg_sys[1][1][jj][ll]-avg_m[1][1][jj][ll])/(avg_sys[1][0][jj][ll]-avg_m[1][0][jj][ll]))/2,2)+(avg_ratio[jj][ll]*(pow(avg_sys[0][1][jj][ll]/avg_m[1][1][jj][ll],2)+pow(avg_sys[0][0][jj][ll]/avg_m[1][0][jj][ll],2)))); //pippo


		   cout<<"AVERAGE: "<<avg_ratio[jj][ll]<<" "<<avg_eratio[jj][ll]<<" "<<avg_z[1][0][jj][ll]<<endl;
		 }
	       }//ll



	     if(jj>0)dmult_avg[0][ii][jj] = new TGraphErrors(12,avg_z[0][ii][jj],avg_m[0][ii][jj],zerror,avg_em[0][ii][jj]);
	     else dmult_avg[0][ii][jj] = new TGraphErrors(12,avg_z[0][ii][jj],avg_m[0][ii][jj],zerror,avg_em[0][ii][jj]);
	     if(jj>0)dmult_avg[1][ii][jj] = new TGraphErrors(12,avg_z[1][ii][jj],avg_m[1][ii][jj],zerror,avg_em[1][ii][jj]);
	     else dmult_avg[1][ii][jj] = new TGraphErrors(12,avg_z[1][ii][jj],avg_m[1][ii][jj],zerror,avg_em[1][ii][jj]);
	     if(jj>0 ){
	       A00[0][ii][jj] = new TGraphAsymmErrors(12,z_avg , yoffset[ii], errorx, errorx, ysys,avg_totsys[0][ii][jj] );
	       A00[1][ii][jj] =  new TGraphAsymmErrors(12,z_avg , yoffset[ii], errorx, errorx, ysys,avg_totsys[1][ii][jj] );
	       if(jj<4){
		 A00b[0][ii][jj] = new TGraphAsymmErrors(8,z_avgb , yoffsetb[ii], errorx, errorx, ysys,avg_sys_plotb[0][ii][jj] );
		 A00b[1][ii][jj] =  new TGraphAsymmErrors(8,z_avgb , yoffsetb[ii], errorx, errorx, ysys,avg_sys_plotb[1][ii][jj] );}
	       else if(jj<6) {
		 A00b[0][ii][jj] = new TGraphAsymmErrors(6,z_avgb , yoffsetb[ii], errorx, errorx, ysys,avg_sys_plotb[0][ii][jj] );
		 A00b[1][ii][jj] =  new TGraphAsymmErrors(6,z_avgb , yoffsetb[ii], errorx, errorx, ysys,avg_sys_plotb[1][ii][jj] );}
	       else {
		 A00b[0][ii][jj] = new TGraphAsymmErrors(4,z_avgb , yoffsetb[ii], errorx, errorx, ysys,avg_sys_plotb[0][ii][jj] );
		 A00b[1][ii][jj] =  new TGraphAsymmErrors(4,z_avgb , yoffsetb[ii], errorx, errorx, ysys,avg_sys_plotb[1][ii][jj] );}
	     }
	     else {
	       A00[0][ii][jj] = new TGraphAsymmErrors(10,z_avg , yoffset[ii], errorx, errorx, ysys,avg_totsys[0][ii][jj] );
	       A00[1][ii][jj] =  new TGraphAsymmErrors(10,z_avg , yoffset[ii], errorx, errorx, ysys,avg_totsys[1][ii][jj] );
	       A00b[0][ii][jj] = new TGraphAsymmErrors(8,z_avgb , yoffsetb[ii], errorx, errorx, ysys,avg_sys_plotb[0][ii][jj] );
	       A00b[1][ii][jj] =  new TGraphAsymmErrors(8,z_avgb , yoffsetb[ii], errorx, errorx, ysys,avg_sys_plotb[1][ii][jj] );}

	     for(int mm=0;mm<2;mm++){
	       if(ii==1){
		 dmult_avg[mm][ii][jj]->SetMarkerColor(kRed);  dmult_avg[mm][ii][jj]->SetLineColor(kRed);
		 dmult_avg[mm][ii][jj]->SetMarkerStyle(kFullCircle); dmult_avg[mm][ii][jj]->SetMarkerSize(1.5);
		 A00[mm][ii][jj]->SetFillColor(kRed);
	       }
	       else if(ii==0){
		 dmult_avg[mm][ii][jj]->SetMarkerColor(kBlue); dmult_avg[mm][ii][jj]->SetLineColor(kBlue);
		 dmult_avg[mm][ii][jj]->SetMarkerStyle(24); dmult_avg[mm][ii][jj]->SetMarkerSize(1.2);
		 A00[mm][ii][jj]->SetFillColor(kBlue);
		 A00b[mm][ii][jj]->SetFillColor(kWhite);
	       }}//mm
	   }//ii

	 h_ratio_2d[jj]= new TGraphErrors(12,avg_z[1][0][jj],avg_ratio[jj],zerror,avg_eratio[jj]);
	 h_ratio_2d[jj]->SetMarkerColor(kRed);  h_ratio_2d[jj]->SetLineColor(kRed);
	 h_ratio_2d[jj]->SetMarkerStyle(kFullCircle); h_ratio_2d[jj]->SetMarkerSize(1.5);
	 if(jj==0) R00[jj] = new TGraphAsymmErrors(10,z_avg , yoffset[2], errorx, errorx, ysys,avg_sratio[jj]);
	 else if(jj<5)R00[jj] = new TGraphAsymmErrors(12,z_avg , yoffset[2], errorx, errorx, ysys,avg_sratio[jj]);
	 else  R00[jj] = new TGraphAsymmErrors(12,z_avg , yoffset[3], errorx, errorx, ysys,avg_sratio[jj]);
	 R00[jj]->SetFillColor(kRed-10);
       }//jj

     double temp1=0;
     double temp2=0;
     double tempq1=0;
     double tempq2=0;
     for(int jj=1; jj<9; jj++)
       {
	 avg_x[jj-1] =0;
	 temp_count[jj-1]=0;
	 avg_q[jj-1] =0;
	 temp_countq[jj-1]=0;
	 for(int kk=0;kk<6;kk++){
	   for(int ll=0;ll<12;ll++){
	     avg_x[jj-1] +=        (x_total[1][jj][kk][ll]*x_counts[1][jj][kk][ll]);
	     temp_count[jj-1]+=(x_counts[1][jj][kk][ll]);
	     avg_q[jj-1] +=        (q_total[1][jj][kk][ll]*q_counts[1][jj][kk][ll] );
	     temp_countq[jj-1]+=(q_counts[1][jj][kk][ll]);
	   }
	 }
	 avg_x[jj-1]=avg_x[jj-1]/temp_count[jj-1];
	 avg_q[jj-1]=avg_q[jj-1]/temp_countq[jj-1];
       }
     for(int kk=0;kk<6;kk++){
       for(int ll=0;ll<12;ll++){
	 temp1 +=        (x_total[1][0][kk][ll]*x_counts[1][0][kk][ll]);
	 temp2 +=(x_counts[1][0][kk][ll]);
	 tempq1 +=        (q_total[1][0][kk][ll]*q_counts[1][0][kk][ll] );
	 tempq2 +=(q_counts[1][0][kk][ll]);

       }
     }

     ////////////sum over z

     double sum[9];
     double sum_e[9];
     double ratio[9];
     double ratio_e[9];
     double sum_m[2][9];
     double sum_em[2][9];
     double int_sysm[2][9];
     double int_sysm_uncorr[2][9];
     double summed_sysm[9];
     double ratio_sysm[9];

     for(int jj=1; jj<9; jj++)
       {
	 sum[jj-1]=0;
	 sum_e[jj-1]=0;
	 ratio[jj-1]=0;
	 ratio_e[jj-1]=0;
	 summed_sysm[jj-1]=0;
	 ratio_sysm[jj-1]=0;
	 for(int ii=0; ii<2; ii++)
	   {
	     int_sysm[ii][jj-1]=0;
	     int_sysm_uncorr[ii][jj-1]=0;
	     sum_m[ii][jj-1]=0;
	     sum_em[ii][jj-1]=0;

	     for(int ll=0; ll<12; ll++)
	       {
		 int_sysm[ii][jj-1]+=avg_sys[1][ii][jj][ll]*zwidth[ll];
		 int_sysm_uncorr[ii][jj-1]+=pow(avg_sys[0][ii][jj][ll]*zwidth[ll],2);
		 if(avg_m[1][ii][jj][ll]>0){sum_m[ii][jj-1] += avg_m[1][ii][jj][ll]*zwidth[ll];
		   sum_em[ii][jj-1] += (pow(avg_em[1][ii][jj][ll]*zwidth[ll],2));}

	       }
	     sum_em[ii][jj-1]=sqrt(sum_em[ii][jj-1]);
	     int_sysm_uncorr[ii][jj-1] = sqrt(int_sysm_uncorr[ii][jj-1]);


	   }
	 sum[jj-1]=sum_m[0][jj-1]+sum_m[1][jj-1];
	 sum_e[jj-1]=sqrt(pow(sum_em[0][jj-1],2)+pow(sum_em[1][jj-1],2));
	 summed_sysm[jj-1]=sqrt(pow(int_sysm[0][jj-1] + int_sysm[1][jj-1],2)+ pow(int_sysm_uncorr[0][jj-1],2)+pow(int_sysm_uncorr[1][jj-1],2));
	 cout<<"sum: "<<avg_x[jj-1]<<" "<<sum[jj-1]<<" stat: "<<sum_e[jj-1]<<" sys "<<summed_sysm[jj-1]<<endl;
	 ratio[jj-1]=sum_m[1][jj-1]/sum_m[0][jj-1];
	 ratio_e[jj-1]=ratio[jj-1]*sqrt(pow(sum_em[0][jj-1]/sum_m[0][jj-1],2)+pow(sum_em[1][jj-1]/sum_m[1][jj-1],2));
	 ratio_sysm[jj-1]=/*(int_sysm[1][jj-1]-ratio[jj-1]*int_sysm[0][jj-1])/(sum_m[0][jj-1]+int_sysm[0][jj-1]);
			   */sqrt( pow(((int_sysm[1][jj-1]+sum_m[1][jj-1])/(int_sysm[0][jj-1]+sum_m[0][jj-1])-(int_sysm[1][jj-1]-sum_m[1][jj-1])/(int_sysm[0][jj-1]-sum_m[0][jj-1]))/2,2)+(ratio[jj-1]*(pow(int_sysm_uncorr[1][jj-1]/sum_m[1][jj-1],2)+pow(int_sysm_uncorr[0][jj-1]/sum_m[0][jj-1],2))));


       }
     TGraphErrors *dmult_sum = new TGraphErrors(8,avg_x, sum ,zerror,sum_e);
     TGraphErrors *dmult_ratio = new TGraphErrors(8,avg_x, ratio,zerror,ratio_e);

     double offset1[8]={0.085,0.085,0.085,0.085,0.085,0.085,0.085,0.085};
     double offset2[9]={0.093,0.093,0.093,0.093,0.093,0.093,0.093,0.093,0.093};
     double offset2r[9]={1.,1.,1.,1.,1.,1.,1.,1.,1.};
     double offset3[7]={0.7,0.7,0.7,0.7,0.7,0.7,0.7};

     TGraphAsymmErrors *dsys_sum = new TGraphAsymmErrors(8,avg_x ,offset1 , errorx, errorx, ysys,summed_sysm );
     TGraphAsymmErrors *dsys_ratio = new TGraphAsymmErrors(8,avg_x ,offset2r , errorx, errorx,  ysys,ratio_sysm );

     dsys_sum->SetFillColor(kRed-10);  dsys_ratio->SetFillColor(kRed-10);
     dmult_sum->SetMarkerColor(kRed);  dmult_sum->SetLineColor(kRed);
     dmult_sum->SetMarkerStyle(kFullCircle); dmult_sum->SetMarkerSize(1.63);
     dmult_ratio->SetMarkerColor(kRed);  dmult_ratio->SetLineColor(kRed);
     dmult_ratio->SetMarkerStyle(kFullCircle); dmult_ratio->SetMarkerSize(1.63);


     /////////////COMPASS multiplicities from 2016

     double psys_a[2][9][12];
     double pm_u[2][9][12];
     double pme_u[2][9][12];
     double pz[2][9][12];
     double px[2][9][12];

     ifstream IN2;
     IN2.open("./data/MultiplicityKaon_yavg.txt", ifstream::in);
     for(int i=0; i<9; i++)
     {
       for(int k=0; k<12; k++)
       {
         IN2>>iii>>kkk>>xa_totalp>>q2a_totalp>>z_kaonp>>mp>>mep>>sysp>>leptoflagp>>xa_totaln>>q2a_totaln>>z_kaonn>>mn>>men>>sysn>>leptoflagn;
         pm_u[1][i][k]=mp;
   	     pme_u[1][i][k]=mep;
   	     psys_a[1][i][k]=sysp;
   	     pz[1][i][k]=z_kaonp;
   	     px[1][i][k]=xa_totalp;
     	   //
     	   pm_u[0][i][k]=mn;
     	   pme_u[0][i][k]=men;
     	   psys_a[0][i][k]=sysn;
     	   pz[0][i][k]=z_kaonn;
     	   px[0][i][k]=xa_totaln;
         kkk=iii=xa_totalp= ya_totalp= q2a_totalp= z_kaonp = mp= mep= sysp= leptoflagp=xa_totaln= ya_totaln= q2a_totaln= z_kaonn=mn= men= sysn= leptoflagn=0;
       }
     }

     IN2.close();
     //
     cout<<"Done reading COMPASS 2016 kaon"<<endl;


     double psum[9];
     double psum_e[9];
     double pratio[9];
     double pratio_e[9];
     double psum_m[2][9];
     double psum_em[2][9];
     //double pint_sysm[2][9];
     double pint_sysm_uncorr[2][9];
     double psummed_sysm[9];
     double pratio_sysm[9];

     for(int jj=1; jj<9; jj++)
     {
       psum[jj-1]=0;
       psum_e[jj-1]=0;
       pratio[jj-1]=0;
       pratio_e[jj-1]=0;
       psummed_sysm[jj-1]=0;
       pratio_sysm[jj-1]=0;
       for(int ii=0; ii<2; ii++)
       {
         pint_sysm[ii][jj-1]=0;
         pint_sysm_uncorr[ii][jj-1]=0;
         psum_m[ii][jj-1]=0;
         psum_em[ii][jj-1]=0;

         for(int ll=0; ll<12; ll++)
         {
           //pint_sysm[ii][jj-1]+=avg_sys[1][ii][ll]*zwidth[ll];
           pint_sysm_uncorr[ii][jj-1]+=psys_a[ii][jj][ll]*pow(zwidth[ll],2);
           if(pm_u[ii][jj][ll]>0){psum_m[ii][jj-1] += pm_u[ii][jj][ll]*zwidth[ll];
             psum_em[ii][jj-1] += pme_u[ii][jj][ll]*pow(zwidth[ll],2);}

         }

         psum_em[ii][jj-1]=sqrt(psum_em[ii][jj-1]);
         pint_sysm_uncorr[ii][jj-1] = sqrt(int_sysm_uncorr[ii][jj-1]);


       }

       psum[jj-1]=psum_m[0][jj-1]+psum_m[1][jj-1];
       psum_e[jj-1]=sqrt(pow(psum_em[0][jj-1],2)+pow(psum_em[1][jj-1],2));
       psummed_sysm[jj-1]=sqrt(pow(pint_sysm_uncorr[0][jj-1],2)+pow(pint_sysm_uncorr[1][jj-1],2));
       pratio[jj-1]=psum_m[1][jj-1]/psum_m[0][jj-1];
       pratio_e[jj-1]=pratio[jj-1]*sqrt(pow(psum_em[0][jj-1]/psum_m[0][jj-1],2)+pow(psum_em[1][jj-1]/psum_m[1][jj-1],2));
       pratio_sysm[jj-1]=/*(int_sysm[1][jj-1]-ratio[jj-1]*int_sysm[0][jj-1])/(sum_m[0][jj-1]+int_sysm[0][jj-1]);
             */sqrt( pow(((pint_sysm[1][jj-1]+psum_m[1][jj-1])/(pint_sysm[0][jj-1]+psum_m[0][jj-1])-(pint_sysm[1][jj-1]-psum_m[1][jj-1])/(pint_sysm[0][jj-1]-psum_m[0][jj-1]))/2,2)+(pratio[jj-1]*(pow(pint_sysm_uncorr[1][jj-1]/psum_m[1][jj-1],2)+pow(pint_sysm_uncorr[0][jj-1]/psum_m[0][jj-1],2))));
     }

     TGraphErrors *pmult_sum = new TGraphErrors(8,px, psum ,pzerror,psum_e);
     TGraphErrors *pmult_ratio = new TGraphErrors(8,px, pratio,pzerror,pratio_e);

     double poffset1[8]={0.085,0.085,0.085,0.085,0.085,0.085,0.085,0.085};
     double poffset2[9]={0.093,0.093,0.093,0.093,0.093,0.093,0.093,0.093,0.093};
     double poffset2r[9]={1.,1.,1.,1.,1.,1.,1.,1.,1.};
     double poffset3[7]={0.7,0.7,0.7,0.7,0.7,0.7,0.7};

     TGraphAsymmErrors *psys_sum = new TGraphAsymmErrors(8,pavg_x ,poffset1 , errorx, errorx, pysys,psummed_sysm );
     TGraphAsymmErrors *psys_ratio = new TGraphAsymmErrors(8,pavg_x ,poffset2r , errorx, errorx,  pysys,pratio_sysm );

     psys_sum->SetFillColor(kBlue-10);  psys_ratio->SetFillColor(kBlue-10);
     pmult_sum->SetMarkerColor(kBlue);  pmult_sum->SetLineColor(kBlue);
     pmult_sum->SetMarkerStyle(kFullCircle); pmult_sum->SetMarkerSize(1.63);
     pmult_ratio->SetMarkerColor(kBlue);  pmult_ratio->SetLineColor(kBlue);
     pmult_ratio->SetMarkerStyle(kFullCircle); pmult_ratio->SetMarkerSize(1.63);



     /////////////HERMES multiplicities from 2013/2014 publications

     double hsys_a[2][9][4];
     double hm_u[2][9][4];
     double hme_u[2][9][4];
     double hz[2][9][4];
     double hx[2][9][4];

     ifstream IN3;
     IN3.open("./data/HERMESK+Prot.txt", ifstream::in);
     for(int k=0; k<4; k++)
     {
       for(int i=0; i<9; i++)
       {
         IN3>>xa_totalp>>xa_totaln>>mp>>mep>>men>>sysp>>sysn;
         hm_u[1][i][k]=mp;
   	     hme_u[1][i][k]=mep;
   	     hsys_a[1][i][k]=sysp;
   	     hx[1][i][k]=xa_totalp;
     	   //
     	   pm_u[0][i][k]=mn;
     	   pme_u[0][i][k]=men;
     	   psys_a[0][i][k]=sysn;
     	   pz[0][i][k]=z_kaonn;
     	   px[0][i][k]=xa_totaln;
         kkk=iii=xa_totalp= ya_totalp= q2a_totalp= z_kaonp = mp= mep= sysp= leptoflagp=xa_totaln= ya_totaln= q2a_totaln= z_kaonn=mn= men= sysn= leptoflagn=0;
       }
     }

     IN3.close();
     //
     cout<<"done reading HERMES K+"<<endl;

     ifstream IN4;
     IN4.open("./data/HERMESK+Prot.txt", ifstream::in);
     for(int k=0; k<4; k++)
     {
       for(int i=0; i<9; i++)
       {
         IN4>>xa_totaln>>xa_totalp>>mn>>men>>mep>>sysn>>sysp;
     	   hm_u[0][i][k]=mn;
     	   hme_u[0][i][k]=men;
     	   hsys_a[0][i][k]=sysn;
     	   hx[0][i][k]=xa_totaln;
         kkk=iii=xa_totalp= ya_totalp= q2a_totalp= z_kaonp = mp= mep= sysp= leptoflagp=xa_totaln= ya_totaln= q2a_totaln= z_kaonn=mn= men= sysn= leptoflagn=0;
       }
     }

     IN4.close();
     //
     cout<<"done reading HERMES K+"<<endl;

     double hsum[9];
     double hsum_e[9];
     double hratio[9];
     double hratio_e[9];
     double hsum_m[2][9];
     double hsum_em[2][9];
     // double hint_sysm[2][9];
     double hint_sysm_uncorr[2][9];
     double hsummed_sysm[9];
     double hratio_sysm[9];

     for(int jj=1; jj<9; jj++)
     {
       hsum[jj-1]=0;
       hsum_e[jj-1]=0;
       hratio[jj-1]=0;
       hratio_e[jj-1]=0;
       hsummed_sysm[jj-1]=0;
       hratio_sysm[jj-1]=0;
       for(int ii=0; ii<2; ii++)
       {
         hint_sysm[ii][jj-1]=0;
         hint_sysm_uncorr[ii][jj-1]=0;
         hsum_m[ii][jj-1]=0;
         hsum_em[ii][jj-1]=0;

         for(int ll=0; ll<4; ll++)
         {
           // hint_sysm[ii][jj-1]+=avg_sys[ii][jj][ll]*zwidth[ll];
           hint_sysm_uncorr[ii][jj-1]+=psys_a[ii][jj][ll]*pow(zwidthH[ll],2);
           if(hm_u[ii][jj][ll]>0){hsum_m[ii][jj-1] += hm_u[ii][jj][ll]*zwidthH[ll];
             hsum_em[ii][jj-1] += hme_u[ii][jj][ll]*pow(zwidthH[ll],2);}

         }

         hsum_em[ii][jj-1]=sqrt(hsum_em[ii][jj-1]);
         hint_sysm_uncorr[ii][jj-1] = sqrt(hint_sysm_uncorr[ii][jj-1]);


       }

       hsum[jj-1]=hsum_m[0][jj-1]+hsum_m[1][jj-1];
       hsum_e[jj-1]=sqrt(pow(hsum_em[0][jj-1],2)+pow(hsum_em[1][jj-1],2));
       hsummed_sysm[jj-1]=sqrt(pow(hint_sysm_uncorr[0][jj-1],2)+pow(hint_sysm_uncorr[1][jj-1],2));
       hratio[jj-1]=hsum_m[1][jj-1]/hsum_m[0][jj-1];
       hratio_e[jj-1]=hratio[jj-1]*sqrt(pow(hsum_em[0][jj-1]/hsum_m[0][jj-1],2)+pow(hsum_em[1][jj-1]/hsum_m[1][jj-1],2));
       hratio_sysm[jj-1]=/*(int_sysm[1][jj-1]-ratio[jj-1]*int_sysm[0][jj-1])/(sum_m[0][jj-1]+int_sysm[0][jj-1]);
             */sqrt( pow(((hint_sysm[1][jj-1]+hsum_m[1][jj-1])/(hint_sysm[0][jj-1]+hsum_m[0][jj-1])-(hint_sysm[1][jj-1]-hsum_m[1][jj-1])/(hint_sysm[0][jj-1]-hsum_m[0][jj-1]))/2,2)+(hratio[jj-1]*(pow(hint_sysm_uncorr[1][jj-1]/hsum_m[1][jj-1],2)+pow(hint_sysm_uncorr[0][jj-1]/hsum_m[0][jj-1],2))));
     }

     for(int jj=0; jj<9; jj++)
     {
	      hsummed_sysmb[jj]=0;
	      hxb[jj] = hx[jj];
     	  hsummed_sysmb[jj]=hsummed_sysm[jj]-0.001;
        hratio_sysmb[jj]=hratio_sysm[jj]-0.001;
     }
     hxb[0]=hxb[0]+0.0005;
     hxb[8]=hxb[8]-0.006;


     double hoffseth2b[9]={1.07+0.005,1.07+0.005,1.07+0.005,1.07+0.005,1.07+0.005,1.07+0.005,1.07+0.005,1.07+0.005,1.07+0.005};
     double hoffseth2[9]={1.07,1.07,1.07,1.07,1.07,1.07,1.07,1.07,1.07};
     double hoffset2h[12]={0.093,0.093,0.093,0.093,0.093,0.093,0.093,0.093,0.093, 0.093, 0.093, 0.093};
     double hoffseth1b[12]={0.093+0.0005,0.093+0.0005,0.093+0.0005,0.093+0.0005,0.093+0.0005,0.093+0.0005,0.093+0.0005,0.093+0.0005,0.093+0.0005,0.093+0.0005,0.093+0.0005,0.093+0.0005};

     TGraphErrors *hmult_sum = new TGraphErrors(9,hx, hsum ,hzerror,hsum_e);
     TGraphAsymmErrors *hsys_sum = new TGraphAsymmErrors(9,hx ,hoffset2h , errorx, errorx, errorx,hsummed_sysm );
     TGraphAsymmErrors *hsys_sumb = new TGraphAsymmErrors(9,hxb ,hoffseth1b , errorx, errorx, errorx,hsummed_sysmb );

     TGraphErrors *hmult_ratio = new TGraphErrors(9,hx, hratio ,hzerror,hratio_e);
     TGraphAsymmErrors *hsys_ratio = new TGraphAsymmErrors(9,hx ,hoffseth2 , errorx, errorx, ysys,hratio_sysm);
     TGraphAsymmErrors *hsys_ratiob = new TGraphAsymmErrors(9,hxb ,hoffseth2b , errorx, errorx, ysys,hratio_sysmb );

     hsys_sum->SetFillColor(kBlack);   hsys_ratio->SetFillColor(kBlack);
     hsys_sumb->SetFillColor(kWhite);   hsys_ratiob->SetFillColor(kWhite);

     hmult_sum->SetMarkerColor(kBlack);  hmult_sum->SetLineColor(kBlack);
     hmult_sum->SetMarkerStyle(24); hmult_sum->SetMarkerSize(1.33);
     hmult_ratio->SetMarkerColor(kBlack);  hmult_ratio->SetLineColor(kBlack);
     hmult_ratio->SetMarkerStyle(24); hmult_ratio->SetMarkerSize(1.32);


   /////////////////template for plots///////////
     gROOT->SetStyle("Plain");
     gStyle->SetPalette(1);
     gStyle->SetOptStat(0);
     gStyle->SetErrorX(0);
     gStyle->SetLabelFont(132, "xy");
     gStyle->SetTitleFont(132,"xy");
     gStyle->SetTitleAlign(33);
     gStyle->SetTitleSize(0.06);
     gStyle->SetPadColor(0);
     gStyle->SetNdivisions(705,  "X");
     gStyle->SetNdivisions(805,  "Y");
     gStyle->SetLabelOffset(0.04, "Y");
     gStyle->SetPadRightMargin(0.02);
     gStyle->SetPadTopMargin(0.025);//0.02
     gStyle->SetPadLeftMargin(0.13);//0.135
     gStyle->SetPadBottomMargin(0.1); //0.13


     TLatex title3;
     title3.SetTextFont(132);
     title3.SetTextSize(0.06);
     title3.SetTextAlign(11);
     title3.SetNDC();
     TLatex y_t_axis_bins;
     y_t_axis_bins.SetTextFont( 132 );
     y_t_axis_bins.SetTextAlign(23);
     y_t_axis_bins.SetTextSize( 0.38 );

     TH2F *mAxis = new TH2F("mAxis","",100,0.,0.8,100,0.08,0.2);//-0.2,3.7}//-0.2,3.7
     mAxis->GetYaxis()->SetLabelSize(0.06);
     mAxis->GetXaxis()->SetLabelSize(0.06);
     mAxis->GetYaxis()->SetTitleSize(0.065);
     mAxis->GetYaxis()->SetTitleOffset(1);
     mAxis->GetYaxis()->SetLabelOffset(0.01);
     mAxis->GetXaxis()->SetTitleSize(0.06);
     mAxis->GetXaxis()->SetLabelOffset(-0.007);
     TCanvas *c_candd   = new TCanvas("c_candd", "Pion sum (x) compass data and theory",970,650);//
     c_candd->SetFillStyle(4000);
     c_candd->SetFrameFillStyle(4000);
     TPad *pad1dd = new TPad("pad1","",0.0,0.01,1,1);
     pad1dd->SetFillStyle(4000); //will be transparent
     pad1dd->SetFrameFillStyle(4000);
     pad1dd->Draw();
     pad1dd->cd();

     ////tag
     gPad->SetLogx();
     dmult_sum->SetMarkerColor(kBlack);
     dmult_sum->SetLineColor(kBlack);
     dsys_sum->SetFillColor(kBlack);
     mAxis->Draw("axis"); dmult_sum->Draw("P"); dsys_sum->Draw("3");
     dsys_sum->SetFillColor(kRed);
     //////////////
     mAxis->Draw("axissame");
     //////////////
     c_candd->cd();
     c_candd->Print("./figures/Mult_sum_test.eps");
     dmult_sum->SetMarkerColor(kRed);
     dmult_sum->SetLineColor(kRed);
     dsys_sum->SetFillColor(kRed-10);
     //endref








     TCanvas *c_can   = new TCanvas("c_acc", "Pion sum (x) compass data and theory",970,650);//
     TPad *spad000 = new TPad("spad000","The first subpad",0.0,0.01,1,1);
     c_can->SetFillStyle(4000);
     c_can->SetFrameFillStyle(4000);
     spad000->SetFillStyle(4000); //will be transparent
     spad000->SetFrameFillStyle(4000);
     spad000->Draw();
     spad000->cd();

     y_t_axis_bins.SetTextAlign(12);
     y_t_axis_bins.SetTextSize( 0.06 );
     TMarker markerp; TMarker markerd; TMarker markerh;
     TH2F *mAxis2 = new TH2F("mAxis2","",100,0.,1,100,0.08,0.2);//-0.2,3.7
     mAxis2->GetYaxis()->SetLabelSize(0.06);
     mAxis2->GetXaxis()->SetLabelSize(0.06);
     mAxis2->GetYaxis()->SetTitleSize(0.065);
     mAxis2->GetYaxis()->SetTitleOffset(1);
     mAxis2->GetYaxis()->SetLabelOffset(0.01);
     mAxis2->GetXaxis()->SetTitleSize(0.06);
     mAxis2->GetXaxis()->SetLabelOffset(-0.007);

     TCanvas *c_can2   = new TCanvas("c_acc2", "Kaon sum (x) compass and hermes data",970,650);//
     TPad *spad00 = new TPad("spad00","The first subpad",0.0,0.01,1,1);
     c_can2->SetFillStyle(4000);
     c_can2->SetFrameFillStyle(4000);
     spad00->SetFillStyle(4000); //will be transparent
     spad00->SetFrameFillStyle(4000);
     spad00->Draw();
     spad00->cd();
     gPad->SetLogx();
     mAxis2->Draw("axis");dsys_sum->Draw("3");psys_sum->Draw("3");  hsys_sum->Draw("3");hsys_sumb->Draw("3");dmult_sum->Draw("P");pmult_sum->Draw("P"); hmult_sum->Draw("P");
     //////////////
     c_can2->cd();
     y_t_axis_bins.SetTextSize( 0.065 );
     y_t_axis_bins.DrawLatex( 0.93,  0.03, "#font[ 12]{x}");
     y_t_axis_bins.SetTextSize( 0.06 );
     markerh.SetMarkerStyle(24);
     markerh.SetMarkerColor(kBlack);
     markerh.SetMarkerSize(1.32); markern.DrawMarker(0.2,.71);
     markerd.SetMarkerStyle(20);
     markerd.SetMarkerColor(kRed);
     markerd.SetMarkerSize(1.63); markerd.DrawMarker(0.2,.80);
     markerp.SetMarkerStyle(20);
     markerp.SetMarkerColor(kBlue);
     markerp.SetMarkerSize(1.63); markerp.DrawMarker(0.2,.89);
     y_t_axis_bins.DrawLatex( 0.22,  0.80, "COMPASS deuteron");
     y_t_axis_bins.DrawLatex( 0.22,  0.89, "COMPASS 2016 preliminary");
     y_t_axis_bins.DrawLatex( 0.22,  0.71, "HERMES proton");
     TMathText text;
     text.SetTextAngle(90);
     text.SetTextFont( 132 );
     text.SetTextSize( 0.07 );
     text.DrawMathText(.05,.64,"\\mathscr{M}^{ K^{+}}+\\mathscr{M}^{ K^{-}}");
     c_can2->Print("./figures/Mult_sum_hermes.eps");
     c_can2->Print("./root_files/Mult_sum_hermes.root");


















     //ratio
     TH2F *rmAxis = new TH2F("rmAxis","",100,0.,1,100,0.9,2.8);//-0.2,3.7
     rmAxis->GetYaxis()->SetLabelSize(0.06);
     rmAxis->GetXaxis()->SetLabelSize(0.06);
     rmAxis->GetYaxis()->SetTitleSize(0.065);
     rmAxis->GetYaxis()->SetTitleOffset(0.9);
     rmAxis->GetYaxis()->SetLabelOffset(0.01);
     rmAxis->GetXaxis()->SetTitleSize(0.06);
     rmAxis->GetXaxis()->SetLabelOffset(-0.007);
     TCanvas *c_can5   = new TCanvas("c_acc5", "Pion ratio (x) compass, hermes data",970,650);//
     TPad *spad1 = new TPad("spad1","The first subpad",0.0,0.01,1,1);
     c_can5->SetFillStyle(4000);
     c_can5->SetFrameFillStyle(4000);
     spad1->SetFillStyle(4000); //will be transparent
     spad1->SetFrameFillStyle(4000);
     spad1->Draw();
     spad1->cd();
     gPad->SetLogx();
     rmAxis->Draw("axis");dsys_ratio->Draw("3"); psys_ratio->Draw("3"); hsys_ratio->Draw("3"); hsys_ratiob->Draw("3"); dmult_ratio->Draw("P");pmult_ratio->Draw("P");hmult_ratio->Draw("P");
     c_can5->cd();
     y_t_axis_bins.SetTextAlign(12);
     y_t_axis_bins.SetTextSize( 0.065 );
     y_t_axis_bins.DrawLatex( 0.93,  0.05, "#font[ 12]{x}");
     y_t_axis_bins.SetTextSize( 0.06 );
     markerh.DrawMarker(0.2,.71);
     markerp.DrawMarker(0.2,.89);
     markerd.DrawMarker(0.2,.80);
     y_t_axis_bins.DrawLatex( 0.22,  0.80, "COMPASS deuteron");
     y_t_axis_bins.DrawLatex( 0.22,  0.89, "COMPASS 2016 preliminary");
     y_t_axis_bins.DrawLatex( 0.22,  0.71, "HERMES proton");
     text.DrawMathText(.05,.65,"\\mathscr{M}^{ K^{+}}/\\mathscr{M}^{ K^{-}}");
     c_can5->Print("./figures/Mult_ratio.eps");
     c_can5->Print("./root_files/Mult_ratio.root");








     // AVG plots
     /////////////////template for plots///////////
     gROOT->SetStyle("Plain");
     gStyle->SetPalette(1);
     gStyle->SetOptStat(0);
     gStyle->SetErrorX(0);
     gStyle->SetLabelFont(132, "xy");
     gStyle->SetTitleFont(132,"xy");
     gStyle->SetTitleAlign(33);
     gStyle->SetTitleSize(0.06);
     gStyle->SetPadColor(0);
     gStyle->SetNdivisions(705,  "X");
     gStyle->SetNdivisions(805,  "Y");
     gStyle->SetLabelOffset(0.04, "Y");
     gStyle->SetPadRightMargin(0.1);
     gStyle->SetPadTopMargin(0.1);
     gStyle->SetPadLeftMargin(0.235);
     gStyle->SetPadBottomMargin(0.13);

     y_t_axis_bins.SetTextAlign(23);
     y_t_axis_bins.SetTextSize( 0.38 );
     ///axis
     TH2F *aAxis[5];
     for(int ij=0; ij<5; ij++){ aAxis[ij] = new TH2F(Form("aAxis_%d",ij),"",100,0.15,0.89,100,-0.1,0.53);}//-0.2,3.7
     aAxis[0]->GetYaxis()->SetLabelSize(0.088);//0.078
     aAxis[0]->GetXaxis()->SetLabelSize(0.0);
     aAxis[0]->GetYaxis()->SetTitle("#frac{d#font[ 12]{M}}{d#font[ 12]{z}}^{#font[ 12]{K}}");
     aAxis[0]->GetYaxis()->SetTitleSize(0.11);
     aAxis[0]->GetYaxis()->SetTitleOffset(0.95);//1
     aAxis[0]->GetYaxis()->SetLabelOffset(0.01);
     aAxis[1]->GetYaxis()->SetLabelSize(0.11);//0.1
     aAxis[1]->GetYaxis()->SetLabelOffset(0.01);
     aAxis[1]->GetXaxis()->SetLabelSize(0.0);
     aAxis[2]->GetYaxis()->SetLabelSize(0.088); //0.078
     aAxis[2]->GetYaxis()->SetLabelOffset(0.01);
     aAxis[2]->GetXaxis()->SetLabelSize(0.084);//0.074
     aAxis[3]->GetYaxis()->SetLabelSize(0.0);
     aAxis[3]->GetXaxis()->SetLabelSize(0.104); //094
     aAxis[3]->GetXaxis()->SetLabelOffset(-0.011);//-0.007
     aAxis[4]->GetYaxis()->SetLabelSize(0.0);
     aAxis[4]->GetXaxis()->SetLabelSize(0.3); //2
     aAxis[4]->GetXaxis()->SetLabelOffset(-0.011);//-0.007
     TLine *La = new TLine(0.16,0.,0.89,0.);
     La->SetLineStyle(2);
     TCanvas *a_can   = new TCanvas("a_can", "multiplicity (x,z) bins with fits",2000,1100);//
     TPad *spada = new TPad("spada","The first subpad",0.0,0.04,1,1);
     a_can->SetFillStyle(4000);
     a_can->SetFrameFillStyle(4000);
     spada->SetFillStyle(4000); //will be transparent
     spada->SetFrameFillStyle(4000);

     spada->Draw();
     spada->cd();
     spada->Divide(5,2,0.001,0.00);
     spada->Draw();
     ////tag

     spada->cd(1);aAxis[0]->Draw("axis");La->Draw();dmult_avg[1][1][0]->Draw("P");dmult_avg[1][0][0]->Draw("P");A00[1][0][0]->Draw("3");A00b[1][0][0]->Draw("3");A00[1][1][0]->Draw("3");
     title3.SetTextSize(0.085);
     title3.SetTextAlign(21);
     title3.DrawLatex(0.65, 0.93,"0.004#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.01");
     spada->cd(2);aAxis[1]->Draw("axis");La->Draw();dmult_avg[1][1][1]->Draw("P");dmult_avg[1][0][1]->Draw("P");A00[1][0][1]->Draw("3");A00b[1][0][1]->Draw("3");A00[1][1][1]->Draw("3");
     title3.SetTextSize(0.102);
     title3.DrawLatex(0.5, 0.93,"0.01#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.02");
     spada->cd(3);aAxis[1]->Draw("axis");La->Draw();dmult_avg[1][1][2]->Draw("P");dmult_avg[1][0][2]->Draw("P");A00[1][0][2]->Draw("3");A00b[1][0][2]->Draw("3");A00[1][1][2]->Draw("3");
     title3.SetTextSize(0.102);
     title3.DrawLatex(0.5, 0.93,"0.02#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.03");
     spada->cd(4);aAxis[1]->Draw("axis");La->Draw();dmult_avg[1][1][3]->Draw("P");dmult_avg[1][0][3]->Draw("P");A00[1][0][3]->Draw("3");A00b[1][0][3]->Draw("3");A00[1][1][3]->Draw("3");
     title3.SetTextSize(0.102);
     title3.DrawLatex(0.5, 0.93,"0.03#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.04");
     spada->cd(5);aAxis[4]->Draw("axis");La->Draw();dmult_avg[1][1][4]->Draw("P");dmult_avg[1][0][4]->Draw("P");A00[1][0][4]->Draw("3");A00b[1][0][4]->Draw("3");A00[1][1][4]->Draw("3");
     title3.SetTextSize(0.102);
     title3.DrawLatex(0.5, 0.93,"0.04#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.06");
     spada->cd(6);aAxis[2]->Draw("axis");La->Draw();dmult_avg[1][1][5]->Draw("P");dmult_avg[1][0][5]->Draw("P");A00[1][0][5]->Draw("3");A00b[1][0][5]->Draw("3");A00[1][1][5]->Draw("3");
     title3.SetTextSize(0.08);
     title3.DrawLatex(0.65, 0.93,"0.06#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.10");
     spada->cd(7);aAxis[3]->Draw("axis");La->Draw();dmult_avg[1][1][6]->Draw("P");dmult_avg[1][0][6]->Draw("P");A00[1][0][6]->Draw("3");A00b[1][0][6]->Draw("3");A00[1][1][6]->Draw("3");
     title3.SetTextSize(0.095);
     title3.DrawLatex(0.5, 0.93,"0.10#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.14");
     spada->cd(8);aAxis[3]->Draw("axis");La->Draw();dmult_avg[1][1][7]->Draw("P");dmult_avg[1][0][7]->Draw("P");A00[1][0][7]->Draw("3");A00b[1][0][7]->Draw("3");A00[1][1][7]->Draw("3");
     title3.SetTextSize(0.095);
     title3.DrawLatex(0.5, 0.93,"0.14#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.18");
     spada->cd(9);aAxis[3]->Draw("axis");La->Draw();dmult_avg[1][1][8]->Draw("P");dmult_avg[1][0][8]->Draw("P");A00[1][0][8]->Draw("3");A00b[1][0][8]->Draw("3");A00[1][1][8]->Draw("3");
     title3.SetTextSize(0.095);
     title3.DrawLatex(0.5, 0.93,"0.18#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.40");
     spada->cd(10);
     TLatex x_axis_bins;
     x_axis_bins.SetTextFont( 132 );
     x_axis_bins.SetTextAlign(33);
     x_axis_bins.SetTextSize( 0.104 );
     x_axis_bins.DrawLatex( 0.13 ,  0.98, "0.2");
     x_axis_bins.DrawLatex( 0.39 ,  0.98, "0.4");
     x_axis_bins.DrawLatex( 0.65 ,  0.98, "0.6");
     x_axis_bins.DrawLatex( 0.93 ,  0.98, "0.8");
     a_can->cd();
     TArrow* z_hline = new TArrow(0.79, 0.529, 0.981,0.529, 0.02, "");
     z_hline -> Draw();
     TArrow* zy_hline = new TArrow(0.795, 0.529, 0.795,0.105, 0.02, "");
     zy_hline -> Draw();
     y_t_axis_bins.SetTextAlign(12);
     y_t_axis_bins.SetTextSize( 0.045 );
     y_t_axis_bins.DrawLatex( 0.76,  0.04, "#font[ 12]{z}");
     markern.SetMarkerStyle(24);
     markern.SetMarkerColor(kBlue);
     markern.SetMarkerSize(1.2); markern.DrawMarker(0.19,.75);
     markerp.SetMarkerStyle(20);
     markerp.SetMarkerColor(kRed);
     markerp.SetMarkerSize(1.5); markerp.DrawMarker(0.19,.8);
     y_t_axis_bins.SetTextSize( 0.04 );
     y_t_axis_bins.DrawLatex( 0.2,  0.75, "K^{#minus}");
     y_t_axis_bins.DrawLatex( 0.2,  0.8, "K^{#plus}");

     a_can->Print("./figures/Mult_K_newnu.eps");
     a_can->Print("./root_files/Mult_K_newnu.root");




  ////////////////ratio 3D


     TH2F *ratio_Axis[5];
     //ratio_Axis[0] = new TH2F(Form("ratio_Axis_%d",0),"",100,0.15,0.89,100,0.5,3.5/*20.7*/);
     //ratio_Axis[1] = new TH2F(Form("ratio_Axis_%d",1),"",100,0.15,0.89,100,0.5,3.5/*20.7*/);
     ratio_Axis[0] = new TH2F(Form("ratio_Axis_%d",0),"",100,0.15,0.89,100,0.3,4.4/*20.7*/);
     ratio_Axis[1] = new TH2F(Form("ratio_Axis_%d",1),"",100,0.15,0.89,100,0.3,4.4/*20.7*/);
     ratio_Axis[2] = new TH2F(Form("ratio_Axis_%d",2),"",100,0.15,0.89,100,-0.7,15.5/*20.7*/);
     ratio_Axis[3] = new TH2F(Form("ratio_Axis_%d",3),"",100,0.15,0.89,100,-0.7,15.5/*20.7*/);
     ratio_Axis[4] = new TH2F(Form("ratio_Axis_%d",4),"",100,0.15,0.89,100,-0.5,15.5/*20.7*/);

     ratio_Axis[0]->GetYaxis()->SetLabelSize(0.088);//0.078
     ratio_Axis[0]->GetXaxis()->SetLabelSize(0.0);
     ratio_Axis[0]->GetYaxis()->SetTitle("#frac{d#font[ 12]{M}}{d#font[ 12]{z}}^{#font[ 12]{K^{#plus}}}/#frac{d#font[ 12]{M}}{d#font[ 12]{z}}^{#font[ 12]{K^{#minus}}}");
     ratio_Axis[0]->GetYaxis()->SetTitleSize(0.11);
     ratio_Axis[0]->GetYaxis()->SetTitleOffset(0.95);//1
     ratio_Axis[0]->GetYaxis()->SetLabelOffset(0.01);
     ratio_Axis[1]->GetYaxis()->SetLabelSize(0.11);//0.1
     ratio_Axis[1]->GetYaxis()->SetLabelOffset(0.01);
     ratio_Axis[1]->GetXaxis()->SetLabelSize(0.0);
     ratio_Axis[2]->GetYaxis()->SetLabelSize(0.088); //0.078
     ratio_Axis[2]->GetYaxis()->SetLabelOffset(0.01);
     ratio_Axis[2]->GetXaxis()->SetLabelSize(0.084);//0.074
     ratio_Axis[3]->GetYaxis()->SetLabelSize(0.0);
     ratio_Axis[3]->GetXaxis()->SetLabelSize(0.104); //094
     ratio_Axis[3]->GetXaxis()->SetLabelOffset(-0.011);//-0.007
     ratio_Axis[4]->GetYaxis()->SetLabelSize(0.0);
     ratio_Axis[4]->GetXaxis()->SetLabelSize(0.3); //2
     ratio_Axis[4]->GetXaxis()->SetLabelOffset(-0.011);//-0.007
     TLine *La = new TLine(0.16,0.,0.89,0.);
     La->SetLineStyle(2);
     TCanvas *ratio_can   = new TCanvas("ratio_can", "multiplicity (x,z) bins with fits",2000,1100);//
     TPad *spada = new TPad("spada","The first subpad",0.0,0.04,1,1);
     ratio_can->SetFillStyle(4000);
     ratio_can->SetFrameFillStyle(4000);
     spada->SetFillStyle(4000); //will be transparent
     spada->SetFrameFillStyle(4000);
     spada->Draw();
     spada->cd();
     spada->Divide(5,2,0.001,0.00);
     spada->Draw();
     ////tag

     spada->cd(1);ratio_Axis[0]->Draw("axis");
     R00[0]->Draw("3");h_ratio_2d[0]->Draw("P");
     title3.SetTextSize(0.09);
     title3.SetTextSize(0.08);
     title3.SetTextSize(0.085);
     title3.SetTextAlign(21);
     title3.DrawLatex(0.65, 0.93,"0.004#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.01");
     spada->cd(2);ratio_Axis[1]->Draw("axis");h_ratio_2d[1]->Draw("P");
     R00[1]->Draw("3");
     title3.SetTextSize(0.102);
     title3.DrawLatex(0.5, 0.93,"0.01#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.02");
     spada->cd(3);ratio_Axis[1]->Draw("axis");h_ratio_2d[2]->Draw("P");
     R00[2]->Draw("3");
     title3.SetTextSize(0.102);
     title3.DrawLatex(0.5, 0.93,"0.02#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.03");
     spada->cd(4);ratio_Axis[1]->Draw("axis");h_ratio_2d[3]->Draw("P");
     R00[3]->Draw("3");
     title3.SetTextSize(0.102);
     title3.DrawLatex(0.5, 0.93,"0.03#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.04");
     spada->cd(5);ratio_Axis[1]->Draw("axis");h_ratio_2d[4]->Draw("P");
     R00[4]->Draw("3");
     title3.SetTextSize(0.102);
     title3.DrawLatex(0.5, 0.93,"0.04#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.06");
     spada->cd(6);ratio_Axis[2]->Draw("axis");h_ratio_2d[5]->Draw("P");
     R00[5]->Draw("3");
     title3.SetTextSize(0.08);
     title3.DrawLatex(0.65, 0.93,"0.06#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.10");
     spada->cd(7);ratio_Axis[3]->Draw("axis");h_ratio_2d[6]->Draw("P");
     R00[6]->Draw("3");
     title3.SetTextSize(0.095);
     title3.DrawLatex(0.5, 0.93,"0.10#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.14");
     spada->cd(8);ratio_Axis[3]->Draw("axis");h_ratio_2d[7]->Draw("P");
     R00[7]->Draw("3");
     title3.SetTextSize(0.095);
     title3.DrawLatex(0.5, 0.93,"0.14#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.18");
     spada->cd(9);ratio_Axis[3]->Draw("axis");h_ratio_2d[8]->Draw("P");
     R00[8]->Draw("3");
     title3.SetTextSize(0.095);
     title3.DrawLatex(0.5, 0.93,"0.18#scale[0.5]{ }<#scale[0.5]{ }#font[ 12]{x}#scale[0.5]{ }<#scale[0.5]{ }0.40");
     spada->cd(10);
     TLatex x_axis_bins;
     x_axis_bins.SetTextFont( 132 );
     x_axis_bins.SetTextAlign(33);
     x_axis_bins.SetTextSize( 0.104 );
     x_axis_bins.DrawLatex( 0.13 ,  0.98, "0.2");
     x_axis_bins.DrawLatex( 0.39 ,  0.98, "0.4");
     x_axis_bins.DrawLatex( 0.65 ,  0.98, "0.6");
     x_axis_bins.DrawLatex( 0.93 ,  0.98, "0.8");
     ratio_can->cd();
     TArrow* z_hline = new TArrow(0.79, 0.529, 0.981,0.529, 0.02, "");
     z_hline -> Draw();
     TArrow* zy_hline = new TArrow(0.795, 0.529, 0.795,0.105, 0.02, "");
     zy_hline -> Draw();
     y_t_axis_bins.SetTextAlign(12);
     y_t_axis_bins.SetTextSize( 0.045 );
     y_t_axis_bins.DrawLatex( 0.76,  0.04, "#font[ 12]{z}");


     ratio_can->Print("./figures/Mult_K_ratio.eps");
     ratio_can->Print("./root_files/Mult_K_ratio.root");



       }

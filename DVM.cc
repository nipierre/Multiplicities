#include "DVM.h"

//Inputs
#define target_file_2012 "data/target-107924-109081.dat"
#define target_file_2016 "data/target-274508-274901.dat"

// Flags
#define Y2006 0
#define Y2012 0
#define Y2016 1
#define SIDIS_XS 227010
#define RHO_XS 25592.8
#define PHI_XS 5995.35

using namespace std;

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
    float z, x, y, r, dummy;
    fin >> z >> dummy >> dummy >> dummy >> dummy >> r >> dummy >> x >> y;
    fZv.push_back(z);
    fXv.push_back(x);
    fYv.push_back(y);
    fRv.push_back(r);
  }
  cout<<"INFO : Target cell description loaded"<<endl;
}

void CellCenter(Double_t z, Double_t& xc, Double_t& yc, Double_t& R)
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

    Double_t rc1 = fRv[i];
    Double_t rc2 = fRv[i+1];

    Double_t dxcdz = (xc2-xc1)/(z2-z1);
    Double_t dycdz = (yc2-yc1)/(z2-z1);
    Double_t drcdz = (rc2-rc1)/(z2-z1);

    Double_t dz = z-z1;
    xc = xc1 + dxcdz*dz;
    yc = yc1 + dycdz*dz;
    R = rc1 + drcdz*dz;

    break;
  }
}

bool InTarget(Double_t xvtx, Double_t yvtx, Double_t zvtx)
{
  Double_t xc, yc, R;
  CellCenter(zvtx, xc, yc, R);
  Double_t dx = xvtx-xc;
  Double_t dy = yvtx-yc;
  Double_t r = sqrt(dx*dx + dy*dy);

  return( r <= R );
}

void readKinCuts(string pFile)
{
  string dummy;
  ifstream list(pFile);
  list >> dummy >> fXmin;
  list >> dummy >> fXmax;
  list >> dummy >> fYmin;
  list >> dummy >> fYmax;
  list >> dummy >> fWmin;
  list >> dummy >> fWmax;
  list >> dummy >> fPmin;
  list >> dummy >> fPmax;
  list.close();
}

void Extraction(string pFilelist, int pType)
{

  //Kinematics
  Double_t Q2 = 0;
  Double_t xBj = 0;
  Double_t yBj = 0;
  Double_t zBj = 0;
  Double_t wBj = 0;
  Double_t nu = 0;

  // Target cells
  if(Y2012) InitTargetFile(target_file_2012);
  else if(Y2016) InitTargetFile(target_file_2016);

  // List of files

  ifstream list(pFilelist);
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
    TBranch *HL04x = (TBranch*) tree->FindBranch("HL04x");
    TBranch *HL04y = (TBranch*) tree->FindBranch("HL04y");
    TBranch *HL05x = (TBranch*) tree->FindBranch("HL05x");
    TBranch *HL05y = (TBranch*) tree->FindBranch("HL05y");
    TBranch *HO03x = (TBranch*) tree->FindBranch("HO03x");
    TBranch *HO03y = (TBranch*) tree->FindBranch("HO03y");
    TBranch *HO04x = (TBranch*) tree->FindBranch("HO04x");
    TBranch *HO04y = (TBranch*) tree->FindBranch("HO04y");
    TBranch *HG01x = (TBranch*) tree->FindBranch("HG01x");
    TBranch *HG01y = (TBranch*) tree->FindBranch("HG01y");
    TBranch *HG021x = (TBranch*) tree->FindBranch("HG021x");
    TBranch *HG021y = (TBranch*) tree->FindBranch("HG021y");
    TBranch *HG022x = (TBranch*) tree->FindBranch("HG022x");
    TBranch *HG022y = (TBranch*) tree->FindBranch("HG022y");
    TBranch *saved = (TBranch*) tree->FindBranch("saved");
    TBranch *BPV = (TBranch*) tree->FindBranch("BPV");
    TBranch *isMuPrim = (TBranch*) tree->FindBranch("isMuPrim");
    TBranch *MZfirst = (TBranch*) tree->FindBranch("MZfirst");
    TBranch *beam_chi2 = (TBranch*) tree->FindBranch("beam_chi2");
    TBranch *mu_prim_chi2 = (TBranch*) tree->FindBranch("mu_prim_chi2");
    TBranch *cellsCrossed = (TBranch*) tree->FindBranch("cellsCrossed");
    TBranch *BMS = (TBranch*) tree->FindBranch("BMS");

    //Hadrons
    TBranch *p = (TBranch*) tree->FindBranch("Hadrons.P");
    TBranch *pt = (TBranch*) tree->FindBranch("Hadrons.pt");
    TBranch *th = (TBranch*) tree->FindBranch("Hadrons.th");
    TBranch *ph = (TBranch*) tree->FindBranch("Hadrons.ph");
    TBranch *ph_pl = (TBranch*) tree->FindBranch("Hadrons.ph_pl");
    TBranch *hXX0 = (TBranch*) tree->FindBranch("Hadrons.XX0");
    TBranch *inHCALacc = (TBranch*) tree->FindBranch("Hadrons.inHCALacc");
    TBranch *HCAL = (TBranch*) tree->FindBranch("Hadrons.HCAL");
    TBranch *charge = (TBranch*) tree->FindBranch("Hadrons.charge");
    TBranch *thRICH = (TBranch*) tree->FindBranch("Hadrons.thRICH");
    //TBranch *LH = (TBranch*) tree->FindBranch("Hadrons.LH");
    TBranch *MCpid = (TBranch*) tree->FindBranch("Hadrons.MCpid");
    //TBranch *MM01x = (TBranch*) tree->FindBranch("Hadrons.MM01x");
    //TBranch *MM01y = (TBranch*) tree->FindBranch("Hadrons.MM01y");
    //TBranch *MM02x = (TBranch*) tree->FindBranch("Hadrons.MM02x");
    //TBranch *MM02y = (TBranch*) tree->FindBranch("Hadrons.MM02y");
    //TBranch *MM03x = (TBranch*) tree->FindBranch("Hadrons.MM03x");
    //TBranch *MM03y = (TBranch*) tree->FindBranch("Hadrons.MM03y");
    //TBranch *Z2Ax = (TBranch*) tree->FindBranch("Hadrons.Z2Ax");
    //TBranch *Z2Ay = (TBranch*) tree->FindBranch("Hadrons.Z2Ay");
    //TBranch *Z2Bx = (TBranch*) tree->FindBranch("Hadrons.Z2Bx");
    //TBranch *Z2By = (TBranch*) tree->FindBranch("Hadrons.Z2By");
    TBranch *RICHx = (TBranch*) tree->FindBranch("Hadrons.RICHx");
    TBranch *RICHy = (TBranch*) tree->FindBranch("Hadrons.RICHy");
    TBranch *chi2_hadron = (TBranch*) tree->FindBranch("Hadrons.chi2_hadron");
    TBranch *HZfirst = (TBranch*) tree->FindBranch("Hadrons.HZfirst");
    TBranch *HZlast = (TBranch*) tree->FindBranch("Hadrons.HZlast");

    //DISMCEvt
    TBranch *MC_vx = (TBranch*) tree->FindBranch("MC_vx");
    TBranch *MC_vy = (TBranch*) tree->FindBranch("MC_vy");
    TBranch *MC_vz = (TBranch*) tree->FindBranch("MC_vz");
    TBranch *MC_p0x = (TBranch*) tree->FindBranch("MC_p0x");
    TBranch *MC_p0y = (TBranch*) tree->FindBranch("MC_p0y");
    TBranch *MC_p0z = (TBranch*) tree->FindBranch("MC_p0z");
    TBranch *MC_p1x = (TBranch*) tree->FindBranch("MC_p1x");
    TBranch *MC_p1y = (TBranch*) tree->FindBranch("MC_p1y");
    TBranch *MC_p1z = (TBranch*) tree->FindBranch("MC_p1z");
    TBranch *irad = (TBranch*) tree->FindBranch("irad");
    TBranch *mcWeight = (TBranch*) tree->FindBranch("mcWeight");
    TBranch *MC_nuTr = (TBranch*) tree->FindBranch("MC_nuTr");
    TBranch *MC_Q2Tr = (TBranch*) tree->FindBranch("MC_Q2Tr");
    TBranch *MC_w = (TBranch*) tree->FindBranch("MC_w");
    TBranch *MC_HM04x = (TBranch*) tree->FindBranch("MC_HM04x");
    TBranch *MC_HM04y = (TBranch*) tree->FindBranch("MC_HM04y");
    TBranch *MC_HM05x = (TBranch*) tree->FindBranch("MC_HM05x");
    TBranch *MC_HM05y = (TBranch*) tree->FindBranch("MC_HM05y");
    TBranch *MC_HL04x = (TBranch*) tree->FindBranch("MC_HL04x");
    TBranch *MC_HL04y = (TBranch*) tree->FindBranch("MC_HL04y");
    TBranch *MC_HL05x = (TBranch*) tree->FindBranch("MC_HL05x");
    TBranch *MC_HL05y = (TBranch*) tree->FindBranch("MC_HL05y");
    TBranch *MC_HO03x = (TBranch*) tree->FindBranch("MC_HO03x");
    TBranch *MC_HO03y = (TBranch*) tree->FindBranch("MC_HO03y");
    TBranch *MC_HO04x = (TBranch*) tree->FindBranch("MC_HO04x");
    TBranch *MC_HO04y = (TBranch*) tree->FindBranch("MC_HO04y");
    TBranch *MC_HG01x = (TBranch*) tree->FindBranch("MC_HG01x");
    TBranch *MC_HG01y = (TBranch*) tree->FindBranch("MC_HG01y");
    TBranch *MC_HG021x = (TBranch*) tree->FindBranch("MC_HG021x");
    TBranch *MC_HG021y = (TBranch*) tree->FindBranch("MC_HG021y");
    TBranch *MC_HG022x = (TBranch*) tree->FindBranch("MC_HG022x");
    TBranch *MC_HG022y = (TBranch*) tree->FindBranch("MC_HG022y");
    TBranch *recons = (TBranch*) tree->FindBranch("recons");
    TBranch *MC_yTr = (TBranch*) tree->FindBranch("MC_yTr");
    TBranch *MC_xTr = (TBranch*) tree->FindBranch("MC_xTr");
    TBranch *MC_TCx = (TBranch*) tree->FindBranch("MC_TCx");
    TBranch *MC_TCy = (TBranch*) tree->FindBranch("MC_TCy");

    //MCHadrons
    TBranch *MC_p = (TBranch*) tree->FindBranch("MCHadrons.P");
    TBranch *MC_th = (TBranch*) tree->FindBranch("MCHadrons.th");
    TBranch *MC_ph = (TBranch*) tree->FindBranch("MCHadrons.ph");
    TBranch *MC_charge = (TBranch*) tree->FindBranch("MCHadrons.charge");
    TBranch *MC_pid = (TBranch*) tree->FindBranch("MCHadrons.pid");
    TBranch *MC_recons = (TBranch*) tree->FindBranch("MCHadrons.recons");
    TBranch *MC_recHadIdx = (TBranch*) tree->FindBranch("MCHadrons.recHadIdx");

    map<int,int> idMCrec[2][4];
    map<int,Double_t> pMCrec[2][4];
    map<int,Double_t> zMCrec[2][4];
    map<int,int> prevMCrec[2][4];
    map<int,int> hidMCrec[2][4];

    // Loopy loop over the events
    int N = (int) tree->GetEntries();

    for (int ip = 0; ip < N; ip++)
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
      HL04x->GetEntry(ip);
      HL04y->GetEntry(ip);
      HL05x->GetEntry(ip);
      HL05y->GetEntry(ip);
      HO03x->GetEntry(ip);
      HO03y->GetEntry(ip);
      HO04x->GetEntry(ip);
      HO04y->GetEntry(ip);
      HG01x->GetEntry(ip);
      HG01y->GetEntry(ip);
      HG021x->GetEntry(ip);
      HG021y->GetEntry(ip);
      HG022x->GetEntry(ip);
      HG022y->GetEntry(ip);
      saved->GetEntry(ip);
      BPV->GetEntry(ip);
      isMuPrim->GetEntry(ip);
      MZfirst->GetEntry(ip);
      beam_chi2->GetEntry(ip);
      mu_prim_chi2->GetEntry(ip);
      cellsCrossed->GetEntry(ip);
      BMS->GetEntry(ip);

      //Hadrons
      p->GetEntry(ip);
      pt->GetEntry(ip);
      th->GetEntry(ip);
      ph->GetEntry(ip);
      ph_pl->GetEntry(ip);
      hXX0->GetEntry(ip);
      inHCALacc->GetEntry(ip);
      HCAL->GetEntry(ip);
      charge->GetEntry(ip);
      thRICH->GetEntry(ip);
      //LH->GetEntry(ip);
      MCpid->GetEntry(ip);
      //MM01x->GetEntry(ip);
      //MM01y->GetEntry(ip);
      //MM02x->GetEntry(ip);
      //MM02y->GetEntry(ip);
      //MM03x->GetEntry(ip);
      //MM03y->GetEntry(ip);
      //Z2Ax->GetEntry(ip);
      //Z2Ay->GetEntry(ip);
      //Z2Bx->GetEntry(ip);
      //Z2By->GetEntry(ip);
      RICHx->GetEntry(ip);
      RICHy->GetEntry(ip);
      chi2_hadron->GetEntry(ip);
      HZfirst->GetEntry(ip);
      HZlast->GetEntry(ip);

      //DISMCEvt
      MC_vx->GetEntry(ip);
      MC_vy->GetEntry(ip);
      MC_vz->GetEntry(ip);
      MC_p0x->GetEntry(ip);
      MC_p0y->GetEntry(ip);
      MC_p0z->GetEntry(ip);
      MC_p1x->GetEntry(ip);
      MC_p1y->GetEntry(ip);
      MC_p1z->GetEntry(ip);
      irad->GetEntry(ip);
      mcWeight->GetEntry(ip);
      MC_nuTr->GetEntry(ip);
      MC_Q2Tr->GetEntry(ip);
      MC_w->GetEntry(ip);
      recons->GetEntry(ip);
      MC_HM04x->GetEntry(ip);
      MC_HM04y->GetEntry(ip);
      MC_HM05x->GetEntry(ip);
      MC_HM05y->GetEntry(ip);
      MC_HL04x->GetEntry(ip);
      MC_HL04y->GetEntry(ip);
      MC_HL05x->GetEntry(ip);
      MC_HL05y->GetEntry(ip);
      MC_HO03x->GetEntry(ip);
      MC_HO03y->GetEntry(ip);
      MC_HO04x->GetEntry(ip);
      MC_HO04y->GetEntry(ip);
      MC_HG01x->GetEntry(ip);
      MC_HG01y->GetEntry(ip);
      MC_HG021x->GetEntry(ip);
      MC_HG021y->GetEntry(ip);
      MC_HG022x->GetEntry(ip);
      MC_HG022y->GetEntry(ip);
      MC_yTr->GetEntry(ip);
      MC_xTr->GetEntry(ip);
      MC_TCx->GetEntry(ip);
      MC_TCy->GetEntry(ip);

      //MCHadrons
      MC_p->GetEntry(ip);
      MC_th->GetEntry(ip);
      MC_ph->GetEntry(ip);
      MC_charge->GetEntry(ip);
      MC_pid->GetEntry(ip);
      MC_recons->GetEntry(ip);
      MC_recHadIdx->GetEntry(ip);

      //--------------------------------------------------------------------------
      //--------- Kinematics -----------------------------------------------------
      //--------------------------------------------------------------------------

      //Data
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
      trig = (trig&2047);


      //2006 ---

      if(Y2006)
      {
        if ((trig&256) && HM05x->GetLeaf("HM05x")->GetValue()<(HM05y->GetLeaf("HM05y")->GetValue()>0 ? 14.55-0.15 : 22.02864-0.12864) )
        {
          trig -= 256;
        }
      }

      //2006 ---

      //--------------------------------------------------------------------------
      //--------- Target ---------------------------------------------------------
      //--------------------------------------------------------------------------

      //2006 ---

      //MC target position new
      static const Double_t dz = 2;

      static const Double_t mcxU = -0.085;
      static const Double_t mcyU = 0.33;
      static const Double_t mczU_1 = -65+dz+4;
      static const Double_t mcxD = -0.085;
      static const Double_t mcyD = 0.33;
      static const Double_t mczD_2 = 65+dz;

      Double_t mcR    = 1.4;

      //target position data 2006
      static const Double_t xU = -0.1;
      static const Double_t yU = 0.33;
      static const Double_t zU_1 = -65+dz+4;
      static const Double_t xD = -0.07;
      static const Double_t yD = 0.33;
      static const Double_t zD_2 =  65+dz;

      Double_t R    = 1.4;
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

      if(pType==0)
      {
        SIDIS_EVENTS++;
        SIDIS_WEIGHT++;
      }
      else if(pType==1)
      {
        RHO_EVENTS++;
        RHO_WEIGHT += mcWeight->GetLeaf("mcWeight")->GetValue();
      }
      else if(pType==2)
      {
        PHI_EVENTS++;
        PHI_WEIGHT += mcWeight->GetLeaf("mcWeight")->GetValue();
      }


      // -----------------------------------------------------------------------
      // --------- DIS Selection -----------------------------------------------
      // -----------------------------------------------------------------------


      // -----------------------------------------------------------------------
      //  Data -----------------------------------------------------------------
      // -----------------------------------------------------------------------

      fAllDISflag = 0;

      // Best Primary Vertex

      // Reconstructed muon
      if((0<E_beam->GetLeaf("E_beam")->GetValue()))
      {

        //BMS (reconstructed beam track)
        if(true) //not used in acceptance
        {

          // Energy of the muon beam
          if((140<E_beam->GetLeaf("E_beam")->GetValue() && E_beam->GetLeaf("E_beam")->GetValue()<180))
          {

            //2006 ---
            if(Y2006)
            {
              // Z coordinate within target regions
              if(((-56<z->GetLeaf("z")->GetValue() && z->GetLeaf("z")->GetValue()<-35)
                    ||(-20<z->GetLeaf("z")->GetValue() && z->GetLeaf("z")->GetValue()<31)
                    ||(43<z->GetLeaf("z")->GetValue() && z->GetLeaf("z")->GetValue()<66)))
              {
                if((mcr < mcR &&  (y->GetLeaf("y")->GetValue()-mcyC)<yCUT
                     && r < R
                     &&  (y->GetLeaf("y")->GetValue()-yC)<yCUT
                     && ((z->GetLeaf("z")->GetValue()>(-65+2+7) && z->GetLeaf("z")->GetValue()<(-35+2-2))
                          ||(z->GetLeaf("z")->GetValue() > (-30+2+8) && z->GetLeaf("z")->GetValue() < (30+2-1))
                          ||(z->GetLeaf("z")->GetValue() > (35+2+6) && z->GetLeaf("z")->GetValue() < (65+2-1)))))
                {

                  // Cells crossing
                  if(true)
                  {

                    // IM/O triggers
                    if((trig&8 || trig&256))
                    {

                      // Q2 cut
                      if((Q2>1))
                      {

                        // y cut
                        if((fYmin<yBj && yBj<fYmax))
                        {

                          // W cut
                          if((fWmin<sqrt(wBj) && sqrt(wBj)<fWmax))
                          {

                            // x cut
                            if((fXmin<xBj && xBj<fXmax))
                            {
                              fAllDISflag = 1;
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }

            }
            //2006 ---

            //2012 ---
            else if(Y2012)
            {
              if(InTarget(x->GetLeaf("x")->GetValue(),y->GetLeaf("y")->GetValue(),z->GetLeaf("z")->GetValue()))
              {
                // Cells crossing
                if(true)
                {
                  if((trig&2 || trig&4 || trig&8))
                  {
                    // Q2 cut
                    if((Q2>1))
                    {
                      // y cut
                      if((fYmin<yBj && yBj<fYmax))
                      {
                        // W cut
                        if((fWmin<sqrt(wBj) && sqrt(wBj)<fWmax))
                        {
                          if((fXmin<xBj && xBj<fXmax))
                          {
                            fAllDISflag = 1;
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
            //2012 ---

            //2016 ---
            else if(Y2016)
            {
              if(InTarget(x->GetLeaf("x")->GetValue(),y->GetLeaf("y")->GetValue(),z->GetLeaf("z")->GetValue())
                  && (-325<z->GetLeaf("z")->GetValue() && z->GetLeaf("z")->GetValue()<-71))
              {
                if((beam_chi2->GetLeaf("beam_chi2")->GetValue()<10))
                {
                  // Cells crossing
                  if((cellsCrossed->GetLeaf("cellsCrossed")->GetValue()))
                  {
                    if((mu_prim_chi2->GetLeaf("mu_prim_chi2")->GetValue()<10))
                    {
                      if((MZfirst->GetLeaf("MZfirst")->GetValue()<350))
                      {
                        if((trig&2 || trig&4 || trig&8 || trig&512))
                        {
                          // Q2 cut
                          if((Q2>1))
                          {
                            // y cut
                            if((fYmin<yBj && yBj<fYmax))
                            {
                              // W cut
                              if((fWmin<sqrt(wBj) && sqrt(wBj)<fWmax))
                              {
                                if((fXmin<xBj && xBj<fXmax))
                                {
                                  fAllDISflag = 1;
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
            //2016 ---
          }
        }
      }

      // -----------------------------------------------------------------------
      // -----------------------------------------------------------------------
      // --------- DIS event calculation ---------------------------------------
      // -----------------------------------------------------------------------
      // -----------------------------------------------------------------------

      // -----------------------------------------------------------------------
      // MC --------------------------------------------------------------------
      // -----------------------------------------------------------------------

      // x Binning

      if(0.<xBj && xBj<0.01) xbin = 0;
      else if(0.01<=xBj && xBj<0.02) xbin = 1;
      else if(0.02<=xBj && xBj<0.03) xbin = 2;
      else if(0.03<=xBj && xBj<0.04) xbin = 3;
      else if(0.04<=xBj && xBj<0.06) xbin = 4;
      else if(0.06<=xBj && xBj<0.1) xbin = 5;
      else if(0.1<=xBj && xBj<0.14) xbin = 6;
      else if(0.14<=xBj && xBj<0.18) xbin = 7;
      else xbin = 8;

      // y Binning

      if(0.<yBj && yBj<0.15) ybin = 0;
      else if(0.15<=yBj && yBj<0.2) ybin = 1;
      else if(0.2<=yBj && yBj<0.3) ybin = 2;
      else if(0.3<=yBj && yBj<0.5) ybin = 3;
      else if(0.5<=yBj && yBj<0.7) ybin = 4;
      else ybin = 5;

      // -----------------------------------------------------------------------
      // -----------------------------------------------------------------------
      // --------- Hadrons Selection -------------------------------------------
      // -----------------------------------------------------------------------
      // -----------------------------------------------------------------------

      // -----------------------------------------------------------------------
      //  Data -----------------------------------------------------------------
      // -----------------------------------------------------------------------

      if(fAllDISflag)
      {
        if(pType==0) fNDIS_evt_SIDIS[xbin][ybin]++;
        else if(pType==1) fNDIS_evt_rho[xbin][ybin]++;
        else if(pType==2) fNDIS_evt_phi[xbin][ybin]++;

        for(int i=0; i<p->GetLeaf("Hadrons.P")->GetLen(); i++)
        {
          // **********************************************************************

          // Hadron identification cuts ------------------------------------------

          if(MCpid->GetLeaf("Hadrons.MCpid")->GetValue(i) == 8)//pi+
          {
            fId = 0;
          }
          else if(MCpid->GetLeaf("Hadrons.MCpid")->GetValue(i) == 9)//pi-
          {
            fId = 1;
          }
          else if(MCpid->GetLeaf("Hadrons.MCpid")->GetValue(i) == 11)//K+
          {
            fId = 2;
          }
          else if(MCpid->GetLeaf("Hadrons.MCpid")->GetValue(i) == 12)//K-
          {
            fId = 3;
          }
          else if(MCpid->GetLeaf("Hadrons.MCpid")->GetValue(i) == 14)//p
          {
            fId = 4;
          }
          else if(MCpid->GetLeaf("Hadrons.MCpid")->GetValue(i) == 15)//pb
          {
            fId = 5;
          }
          else if(MCpid->GetLeaf("Hadrons.MCpid")->GetValue(i) == 3)//e-
          {
            fId = 8;
          }
          else if(MCpid->GetLeaf("Hadrons.MCpid")->GetValue(i) == 2)//e+
          {
            fId = 9;
          }
          else//Hadron
          {
            if(charge->GetLeaf("Hadrons.charge")->GetValue(i)==1 && MCpid->GetLeaf("Hadrons.MCpid")->GetValue(i)>7)
            {
              fId = 6;
            }
            else if(charge->GetLeaf("Hadrons.charge")->GetValue(i)==-1 && MCpid->GetLeaf("Hadrons.MCpid")->GetValue(i)>7)
            {
              fId = 7;
            }
            else
            {
              continue;
            }
          }

          // **********************************************************************

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

          // Chi2/ndf
          if(!(chi2_hadron->GetLeaf("Hadrons.chi2_hadron")->GetValue(i) < 10)) continue;

          // Zfirst
          if(!(HZfirst->GetLeaf("Hadrons.HZfirst")->GetValue(i)<350)) continue;

          // Zlast
          if(!(350<HZlast->GetLeaf("Hadrons.HZlast")->GetValue(i))) continue;

          // Momentum cut (12 GeV to 40 GeV, increasing to 3 GeV to 40 GeV)
          if(!(fPmin<p->GetLeaf("Hadrons.P")->GetValue(i) && p->GetLeaf("Hadrons.P")->GetValue(i)<fPmax)) continue;

          // Theta cut
          if(!(0.01<thRICH->GetLeaf("Hadrons.thRICH")->GetValue(i) && thRICH->GetLeaf("Hadrons.thRICH")->GetValue(i)<0.12)) continue;

          // RICH position cut
          if(!(pow(RICHx->GetLeaf("Hadrons.RICHx")->GetValue(i),2)+pow(RICHy->GetLeaf("Hadrons.RICHy")->GetValue(i),2)>25)) continue;

          // z cut
          if(!(0.2<zBj && zBj<0.85)) continue;

          if(0.2<zBj && zBj<0.25) zbin = 0;
          else if(0.25<=zBj && zBj<0.30) zbin = 1;
          else if(0.30<=zBj && zBj<0.35) zbin = 2;
          else if(0.35<=zBj && zBj<0.40) zbin = 3;
          else if(0.40<=zBj && zBj<0.45) zbin = 4;
          else if(0.45<=zBj && zBj<0.50) zbin = 5;
          else if(0.50<=zBj && zBj<0.55) zbin = 6;
          else if(0.55<=zBj && zBj<0.60) zbin = 7;
          else if(0.60<=zBj && zBj<0.65) zbin = 8;
          else if(0.65<=zBj && zBj<0.70) zbin = 9;
          else if(0.70<=zBj && zBj<0.75) zbin = 10;
          else zbin = 11;


          // **********************************************************************

          // Save of hadrons

          if(fId==0)
          {
            if(pType==0)
            {
              fSIDIS[xbin][ybin][zbin].tab[1][0][0] += 1;
              fSIDIS[xbin][ybin][zbin].tab[1][0][3] += 1;
            }
            else if(pType==1)
            {
              fRho[xbin][ybin][zbin].tab[1][0][0] += 1;
              fRho[xbin][ybin][zbin].tab[1][0][3] += 1;
            }
            else if(pType==2)
            {
              fPhi[xbin][ybin][zbin].tab[1][0][0] += 1;
              fPhi[xbin][ybin][zbin].tab[1][0][3] += 1;
            }
          }
          else if(fId==1)
          {
            if(pType==0)
            {
              fSIDIS[xbin][ybin][zbin].tab[0][0][0] += 1;
              fSIDIS[xbin][ybin][zbin].tab[0][0][3] += 1;
            }
            else if(pType==1)
            {
              fRho[xbin][ybin][zbin].tab[0][0][0] += 1;
              fRho[xbin][ybin][zbin].tab[0][0][3] += 1;
            }
            else if(pType==2)
            {
              fPhi[xbin][ybin][zbin].tab[0][0][0] += 1;
              fPhi[xbin][ybin][zbin].tab[0][0][3] += 1;
            }
          }
          else if(fId==2)
          {
            if(pType==0)
            {
              fSIDIS[xbin][ybin][zbin].tab[1][0][1] += 1;
              fSIDIS[xbin][ybin][zbin].tab[1][0][3] += 1;
            }
            else if(pType==1)
            {
              fRho[xbin][ybin][zbin].tab[1][0][1] += 1;
              fRho[xbin][ybin][zbin].tab[1][0][3] += 1;
            }
            else if(pType==2)
            {
              fPhi[xbin][ybin][zbin].tab[1][0][1] += 1;
              fPhi[xbin][ybin][zbin].tab[1][0][3] += 1;
            }
          }
          else if(fId==3)
          {
            if(pType==0)
            {
              fSIDIS[xbin][ybin][zbin].tab[0][0][1] += 1;
              fSIDIS[xbin][ybin][zbin].tab[0][0][3] += 1;
            }
            else if(pType==1)
            {
              fRho[xbin][ybin][zbin].tab[0][0][1] += 1;
              fRho[xbin][ybin][zbin].tab[0][0][3] += 1;
            }
            else if(pType==2)
            {
              fPhi[xbin][ybin][zbin].tab[0][0][1] += 1;
              fPhi[xbin][ybin][zbin].tab[0][0][3] += 1;
            }
          }
          else if(fId==4)
          {
            if(pType==0)
            {
              fSIDIS[xbin][ybin][zbin].tab[1][0][2] += 1;
              fSIDIS[xbin][ybin][zbin].tab[1][0][3] += 1;
            }
            else if(pType==1)
            {
              fRho[xbin][ybin][zbin].tab[1][0][2] += 1;
              fRho[xbin][ybin][zbin].tab[1][0][3] += 1;
            }
            else if(pType==2)
            {
              fPhi[xbin][ybin][zbin].tab[1][0][2] += 1;
              fPhi[xbin][ybin][zbin].tab[1][0][3] += 1;
            }
          }
          else if(fId==5)
          {
            if(pType==0)
            {
              fSIDIS[xbin][ybin][zbin].tab[0][0][2] += 1;
              fSIDIS[xbin][ybin][zbin].tab[0][0][3] += 1;
            }
            else if(pType==1)
            {
              fRho[xbin][ybin][zbin].tab[0][0][2] += 1;
              fRho[xbin][ybin][zbin].tab[0][0][3] += 1;
            }
            else if(pType==2)
            {
              fPhi[xbin][ybin][zbin].tab[0][0][2] += 1;
              fPhi[xbin][ybin][zbin].tab[0][0][3] += 1;
            }
          }
          else if(fId==6)
          {
            if(pType==0)
            {
              fSIDIS[xbin][ybin][zbin].tab[1][0][3] += 1;
            }
            else if(pType==1)
            {
              fRho[xbin][ybin][zbin].tab[1][0][3] += 1;
            }
            else if(pType==2)
            {
              fPhi[xbin][ybin][zbin].tab[1][0][3] += 1;
            }
          }
          else if(fId==7)
          {
            if(pType==0)
            {
              fSIDIS[xbin][ybin][zbin].tab[0][0][3] += 1;
            }
            else if(pType==1)
            {
              fRho[xbin][ybin][zbin].tab[0][0][3] += 1;
            }
            else if(pType==2)
            {
              fPhi[xbin][ybin][zbin].tab[0][0][3] += 1;
            }
          }
          else
          {
            continue;
          }
        }
      }
    }

    cout << "\n-> Finished processing file " << filename << " <-\n" << endl;

    delete f;
  }
}

void DVMCalc()
{
  for(int i=0; i<9; i++)
  {
    for(int j=0; j<6; j++)
    {
      fNDIS_evt_SIDIS[i][j] /= SIDIS_WEIGHT/SIDIS_XS;
      fNDIS_evt_rho[i][j] /= RHO_WEIGHT/RHO_XS;
      fNDIS_evt_phi[i][j] /= PHI_WEIGHT/PHI_XS;
      fDVM_DIS_pi[i][j] = fNDIS_evt_rho[i][j]/fNDIS_evt_SIDIS[i][j];
      fDVM_DIS_K[i][j] = fNDIS_evt_phi[i][j]/fNDIS_evt_SIDIS[i][j];
      fNDIS_SIDIS_tot += fNDIS_evt_SIDIS[i][j];
      fNDIS_rho_tot += fNDIS_evt_rho[i][j];
      fNDIS_phi_tot += fNDIS_evt_phi[i][j];
      for(int k=0; k<12; k++)
      {
        fSIDIS[i][j][k].tab[1][0][0] /= SIDIS_WEIGHT/SIDIS_XS;
        fSIDIS[i][j][k].tab[0][0][0] /= SIDIS_WEIGHT/SIDIS_XS;
        fSIDIS[i][j][k].tab[1][0][1] /= SIDIS_WEIGHT/SIDIS_XS;
        fSIDIS[i][j][k].tab[0][0][1] /= SIDIS_WEIGHT/SIDIS_XS;
        fSIDIS[i][j][k].tab[1][0][3] /= SIDIS_WEIGHT/SIDIS_XS;
        fSIDIS[i][j][k].tab[0][0][3] /= SIDIS_WEIGHT/SIDIS_XS;
        fSIDIS_tot[1][0] += fSIDIS[i][j][k].tab[1][0][0];
        fSIDIS_tot[0][0] += fSIDIS[i][j][k].tab[0][0][0];
        fSIDIS_tot[1][1] += fSIDIS[i][j][k].tab[1][0][1];
        fSIDIS_tot[0][1] += fSIDIS[i][j][k].tab[0][0][1];
        fSIDIS_tot[1][3] += fSIDIS[i][j][k].tab[1][0][3];
        fSIDIS_tot[0][3] += fSIDIS[i][j][k].tab[0][0][3];
        fRho[i][j][k].tab[1][0][0] /= RHO_WEIGHT/RHO_XS;
        fRho[i][j][k].tab[0][0][0] /= RHO_WEIGHT/RHO_XS;
        fRho[i][j][k].tab[1][0][3] /= RHO_WEIGHT/RHO_XS;
        fRho[i][j][k].tab[0][0][3] /= RHO_WEIGHT/RHO_XS;
        fRho_tot[1][0] += fRho[i][j][k].tab[1][0][0];
        fRho_tot[0][0] += fRho[i][j][k].tab[0][0][0];
        fRho_tot[1][3] += fRho[i][j][k].tab[1][0][3];
        fRho_tot[0][3] += fRho[i][j][k].tab[0][0][3];
        fPhi[i][j][k].tab[1][0][1] /= PHI_WEIGHT/PHI_XS;
        fPhi[i][j][k].tab[0][0][1] /= PHI_WEIGHT/PHI_XS;
        fPhi[i][j][k].tab[1][0][3] /= PHI_WEIGHT/PHI_XS;
        fPhi[i][j][k].tab[0][0][3] /= PHI_WEIGHT/PHI_XS;
        fPhi_tot[1][1] += fPhi[i][j][k].tab[1][0][1];
        fPhi_tot[0][1] += fPhi[i][j][k].tab[0][0][1];
        fPhi_tot[1][3] += fPhi[i][j][k].tab[1][0][3];
        fPhi_tot[0][3] += fPhi[i][j][k].tab[0][0][3];
        fDVM_h[i][j][k].tab[1][0][0] = fRho[i][j][k].tab[1][0][0]/(fRho[i][j][k].tab[1][0][0]+fSIDIS[i][j][k].tab[1][0][0]);
        fDVM_h[i][j][k].tab[0][0][0] = fRho[i][j][k].tab[0][0][0]/(fRho[i][j][k].tab[0][0][0]+fSIDIS[i][j][k].tab[0][0][0]);
        fDVM_h[i][j][k].tab[1][0][1] = fPhi[i][j][k].tab[1][0][1]/(fPhi[i][j][k].tab[1][0][1]+fSIDIS[i][j][k].tab[1][0][1]);
        fDVM_h[i][j][k].tab[0][0][1] = fPhi[i][j][k].tab[0][0][1]/(fPhi[i][j][k].tab[0][0][1]+fSIDIS[i][j][k].tab[0][0][1]);

      }
    }
  }
}

void DVMSaver()
{
  ofstream ofs_dvm("DVM.dat", std::ofstream::out | std::ofstream::trunc);
  for(int i=0; i<9; i++)
  {
    for(int j=0; j<6; j++)
    {
      for(int k=0; k<12; k++)
      {
        ofs_dvm << i+1 << " " << j+1 << " " << k+1 << " "
                << fDVM_h[i][j][k].tab[1][0][0] << " " << fDVM_DIS_pi[i][j]
                << fDVM_h[i][j][k].tab[0][0][0] << " " << fDVM_DIS_pi[i][j]
                << fDVM_h[i][j][k].tab[1][0][1] << " " << fDVM_DIS_K[i][j]
                << fDVM_h[i][j][k].tab[0][0][1] << " " << fDVM_DIS_K[i][j] << endl;
      }
    }
  }
  ofs_dvm.close();
}

int main(int argc, char **argv)
{

  if(argc < 4)
  {
    cout << "ERROR : Not enough arguments." << endl;
    cout << "Asked : 4 *** Received : " << argc-1 << endl;
    cout << "./compMCMC [SIDIS filelist] [rho filelist] [phi filelist] [Cutfile]" << endl;

    return 1;
  }

  SIDIS_WEIGHT = 0;
  SIDIS_EVENTS = 0;
  RHO_WEIGHT = 0;
  RHO_EVENTS = 0;
  PHI_WEIGHT = 0;
  PHI_EVENTS = 0;

  readKinCuts(argv[4]);
  Extraction(argv[1],0);
  Extraction(argv[2],1);
  Extraction(argv[3],2);
  DVMCalc();
  DVMSaver();

  cout << "\n\n";
  cout << "             ********* Luminosity, number of events and number of hadrons ********* " << endl;
  cout << "             ---------------------------------------------------------------------- " << endl;

  cout << '|' << setw(30) << "" << '|' << setw(15) << "DJANGOH" << '|' << setw(15) << "Rho^0" << '|' << setw(15) << "Phi" << '|' << endl;
  cout << '|' << setw(30) << "Generated Events" << '|' << setw(15) << SIDIS_EVENTS << '|' << setw(15) << RHO_EVENTS << '|' << setw(15) << PHI_EVENTS << '|' << endl;
  cout << '|' << setw(30) << "Weighted Gen. Events" << '|' << setw(15) << SIDIS_WEIGHT << '|' << setw(15) << RHO_WEIGHT << '|' << setw(15) << PHI_WEIGHT << '|' << endl;
  cout << '|' << setw(30) << "Integrated XS [pb]" << '|' << setw(15) << SIDIS_XS << '|' << setw(15) << RHO_XS << '|' << setw(15) << PHI_XS << '|' << endl;
  cout << '|' << setw(30) << "MC Luminosity [pb-1]" << '|' << setw(15) << SIDIS_WEIGHT/SIDIS_XS << '|' << setw(15) << RHO_WEIGHT/RHO_XS << '|' << setw(15) << PHI_WEIGHT/PHI_XS << '|' << endl;
  cout << "             ------------------------------------------------------------------------------------------------------------------------------ " << endl;
  cout << '|' << setw(30) << "DIS Events [pb]" << '|' << setw(15) << fNDIS_SIDIS_tot << '|' << setw(15) << fNDIS_rho_tot << '|' << setw(15) << fNDIS_phi_tot << '|' << endl;
  cout << '|' << setw(30) << "h+ [pb]" << '|' << setw(15) << fSIDIS_tot[1][3] << '|' << setw(15) << fRho_tot[1][3] << '|' << setw(15) << fPhi_tot[1][3] << '|' << endl;
  cout << '|' << setw(30) << "h- [pb]" << '|' << setw(15) << fSIDIS_tot[0][3] << '|' << setw(15) << fRho_tot[0][3] << '|' << setw(15) << fPhi_tot[0][3] << '|' << endl;
  cout << '|' << setw(30) << "pi+ [pb]" << '|' << setw(15) << fSIDIS_tot[1][0] << '|' << setw(15) << fRho_tot[1][0] << '|' << setw(15) << "-" << '|' << endl;
  cout << '|' << setw(30) << "pi- [pb]" << '|' << setw(15) << fSIDIS_tot[0][0] << '|' << setw(15) << fRho_tot[0][0] << '|' << setw(15) << "-" << '|' << endl;
  cout << '|' << setw(30) << "K+ [pb]" << '|' << setw(15) << fSIDIS_tot[1][1] << '|' << setw(15) << "-" << '|' << setw(15) << fPhi_tot[1][1] << '|' << endl;
  cout << '|' << setw(30) << "K- [pb]" << '|' << setw(15) << fSIDIS_tot[0][1] << '|' << setw(15) << "-" << '|' << setw(15) << fPhi_tot[0][1] << '|' << endl;

  return 0;
}
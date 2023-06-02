const unsigned long nevents = 40000000;
int commTriggerID[]={2,3,16,17,18,19,20};
int prodTriggerID[]={900501,900502,900503,900504,900505,900506,900507};

void readMudst_upcjet(
	       const char* mudstfile = "/star/data12/reco/./production_AuAu_2023/ReversedFullField/dev/2023/148/24148005/st_upcjet_24148005_raw_0000009.MuDst.root", 
         const char* outfile="Mudst.test.output.root")
{
  gROOT->Macro("loadMuDst.C");
  gROOT->Macro("LoadLogger.C");
  //gROOT->Macro("LoadJetMaker.C");
  ///*
  gSystem->Load("StMcEvent");
  gSystem->Load("StMcEventMaker");
  gSystem->Load("StEEmcUtil");
  gSystem->Load("StEmcRawMaker");
  gSystem->Load("StEmcADCtoEMaker");
  gSystem->Load("StDetectorDbMaker");
  gSystem->Load("StTpcDb");
  gSystem->Load("StDbUtilities");
  gSystem->Load("StDaqLib");
  gSystem->Load("StDbBroker");
  gSystem->Load("St_db_Maker");
  gSystem->Load("StTrsMaker");	
  gSystem->Load("StEpdUtil");
 
  //*/
  // Create chain
  StChain* chain = new StChain;
  //
  // GEANT reader
  StMuDstMaker* muDstMaker = new StMuDstMaker(0,0,"",mudstfile,"",1000);
  StMuDbReader *muDstDb = StMuDbReader::instance();
  St_db_Maker* starDb = new St_db_Maker("StarDb","MySQL:StarDb");

  chain->Init();

  TFile *fout = new TFile(outfile, "recreate");
  
  const int num_trgs = 7;
  const int direction = 2;
  
  TH1F* hntrk[num_trgs][direction];
  TH1F *hvtxz[num_trgs][direction];
  TH1F* htrketa[num_trgs][direction];
  TH1F *epdeta[num_trgs][direction];
  TH1F *zdcadc_east[num_trgs][direction];
  TH1F *zdcadc_west[num_trgs][direction];
  TH1F *eemc[num_trgs][direction];
  
  TH2F *zdcadc[num_trgs][direction];
  TH2F *bbcadc[num_trgs][direction];
  TH2F *vpdadc[num_trgs][direction];
  TH2F *epdadc[num_trgs][direction];

  TH2F *pt2D[num_trgs][direction];
  TH2F *eta2D[num_trgs][direction];
  TH2F *phi2D[num_trgs][direction];

  for(int i=0;i<num_trgs;i++){
    for(int j=0;j<direction;j++){
      hntrk[i][j] = new TH1F(Form("hntrk_%d_%d",i,j),";N_{trk};counts",200,0,200);
      //vertex z
      hvtxz[i][j] = new TH1F(Form("hvtxz_%d_%d",i,j),";z;counts",200,-200,200);
      //track
      htrketa[i][j] = new TH1F(Form("htrketa_%d_%d",i,j),"#eta",100,-2.5,2.5);
      //epd
      epdeta[i][j] = new TH1F(Form("epdeta_%d_%d",i,j),";#eta;counts",256,-6,6);
      //zdc
      zdcadc_east[i][j] = new TH1F(Form("zdcadc_east_%d_%d",i,j),";ZDCE ADC;counts",256,0,2000);
      zdcadc_west[i][j] = new TH1F(Form("zdcadc_west_%d_%d",i,j),";ZDCW ADC;counts",256,0,2000);
      //eemc
      eemc[i][j] = new TH1F(Form("eemc_%d_%d",i,j),"; EEMC ADC;counts",256,0,25000);
      //zdc ADC
      zdcadc[i][j] = new TH2F(Form("zdcadc_%d_%d",i,j),";ZDCE ADC; ZDCW ADC;counts",256,0,2000,256,0,2000);
      //bbc ADC
      bbcadc[i][j] = new TH2F(Form("bbcadc_%d_%d",i,j),";BBCE ADC; BBCW ADC;counts",256,0,25000,256,0,25000);
      //vpd adc
      vpdadc[i][j] = new TH2F(Form("vpdadc_%d_%d",i,j),";VPDE ADC; VPDW ADC;counts",256,0,4000,256,0,4000);
      //epd adc
      epdadc[i][j] = new TH2F(Form("epdadc_%d_%d",i,j),";EPDE ADC; EPDW ADC;counts",256,0,40000,256,0,40000);
      //pt2D
      pt2D[i][j] = new TH2F(Form("pt2D_%d_%d",i,j),";leading p_{T} (GeV/c); subleading p_{T} (GeV/c);counts",100,0,20,100,0,20);
      //eta2D
      eta2D[i][j] = new TH2F(Form("eta2D_%d_%d",i,j),";leading #eta; subleading #eta;counts",100,-3,3,100,-3,3);
      //phi2D
      phi2D[i][j] = new TH2F(Form("phi2D_%d_%d",i,j),";leading #phi; subleading #phi;counts",100,-3.14,3.14,100,-3.14,3.14);
    }
  }

  for (int iEvent = 1; iEvent <= nevents; ++iEvent) {
    chain->Clear();
    int status = chain->Make(iEvent);
    if (status % 10 == kStEOF || status % 10 == kStFatal) break;
    StTriggerData *trgdata = StMuDst::event()->triggerData();
    /**********BBC trigger bits**************/
    int dsm = trgdata->lastDSM(1);
    int BBCE = 0, BBCW = 0;
    BBCE = dsm & 0x2;
    BBCW = dsm & 0x4;

    /************ZDC ADC***************/ 
    /******east = 0, west =1**********/
    int mNZdcPmt = 3;
    double zdcadc_e = 0;
    double zdcadc_w = 0;
    for (Int_t ipmt=1; ipmt<=mNZdcPmt; ipmt++) {
        zdcadc_e = zdcadc_e + trgdata->zdcADC(0, ipmt); 
        zdcadc_w = zdcadc_w + trgdata->zdcADC(1, ipmt);
    }
    // ZDC reqirement 1nXn
    int ZDCW=0;
    int ZDCE=0;
    if(zdcadc_w >= 250){ZDCW=1;}
    if(zdcadc_e >= 250){ZDCE=1;}

    //veto ZDC and EEMC
    // StZdcTriggerDetector &muZdc = StMuDst::event()->zdcTriggerDetector();
    // StEmcTriggerDetector &muemc = StMuDst::event()->emcTriggerDetector();

    //check BBC, EPD, VPD
    StBbcTriggerDetector &muBbc = StMuDst::event()->bbcTriggerDetector();
    StMuPrimaryVertex *vertex = StMuDst::primaryVertex();
    StEpdGeom * mEpdGeom = new StEpdGeom();

    //vertex selections:
    if(!vertex) continue;
    TVector3 pRcVx(vertex->position().x(), vertex->position().y(), vertex->position().z());
    //Tof or BEMC match
    if(vertex->nBTOFMatch()<1 && vertex->nBEMCMatch()<1)continue;

    //trigger decisions
    mMuDst = muDstMaker->muDst();
    StMuEvent *mMuEvent = mMuDst->event();
    //loop over triggers
    
    for(int i=0;i<num_trgs;i++){
      int trig_index=-1;
      if( mMuEvent->triggerIdCollection().nominal().isTrigger(commTriggerID[i]) ) trig_index=i;
      if(trig_index<0) continue; //protection

      //1) now decide direction based on BBC first.
      int direction_index=-1;
      if(trig_index==0){//bar1 trigger both sides are valid
        if(BBCE && !BBCW) direction_index=0;
        if(!BBCE && BBCW) direction_index=1;
      }
      else{//rest only west or neither
        if(!BBCE && !BBCW) direction_index=0;
        if(!BBCE && BBCW) direction_index=1;
      }
      if(direction_index<0) continue; //protection

      //2) filling ZDC distributions based on BBC configurations
      zdcadc_east[trig_index][direction_index]->Fill(zdcadc_e);
      zdcadc_west[trig_index][direction_index]->Fill(zdcadc_w);
      zdcadc[trig_index][direction_index]->Fill(zdcadc_e,zdcadc_w);
      
      //3) now add ZDC to BBC requirements
      direction_index=-1;
      if(trig_index==0){//bar1 trigger both sides are valid
        if(BBCE && !BBCW && ZDCE && !ZDCW) direction_index=0;
        if(!BBCE && BBCW && !ZDCE && ZDCW) direction_index=1;
      }
      else{//rest only west or neither
        if(!BBCE && !BBCW && !ZDCE && !ZDCW) direction_index=0;
        if(!BBCE && BBCW && !ZDCE && ZDCW) direction_index=1;
      }
      if(direction_index<0) continue; //protection again
     
    /*********** EEMC ************/
      int eemc_adc = 0;
      //cout<<muemc->patchEndcap()<<endl;
      for (Int_t ipatch=0; ipatch<90; ipatch++){
  	     eemc_adc = eemc_adc + trgdata->eemcJetPatch(ipatch);
      }
      eemc[trig_index][direction_index]->Fill(eemc_adc);

    /*************BBC*************/
      double bbcadc_e = 0;
      double bbcadc_w = 0;
      for (Int_t ipmt=1; ipmt<=16; ipmt++) {
  	      bbcadc_e = bbcadc_e + trgdata->bbcADC(0, ipmt);	
          bbcadc_w = bbcadc_w + trgdata->bbcADC(1, ipmt);
      }
      bbcadc[trig_index][direction_index]->Fill(bbcadc_e, bbcadc_w);

    /*************VPD*************/
      double vpdadc_e = 0;
      double vpdadc_w = 0;
      for (Int_t ipmt=1; ipmt<=16; ipmt++) {
      	vpdadc_e = vpdadc_e + trgdata->vpdADC(0, ipmt);	
        vpdadc_w = vpdadc_w + trgdata->vpdADC(1, ipmt);
      }
      vpdadc[trig_index][direction_index]->Fill(vpdadc_e, vpdadc_w);

    /*********EPD ************/
      int nhit = StMuDst::numberOfEpdHit();
      //cout<<"nhit"<<nhit<<endl;
      double epdadc_e = 0;
      double epdadc_w = 0;
      for(int iEpdHit = 0; iEpdHit < nhit; iEpdHit++){
        const StMuEpdHit* muEpdHit = StMuDst::epdHit(iEpdHit);
        int adc = muEpdHit->adc();
        if( muEpdHit->side()<0) epdadc_e += adc;
        else epdadc_w += adc;

        TVector3 StraightLine = mEpdGeom->RandomPointOnTile(muEpdHit->id()) - pRcVx;
        double phi_epd = StraightLine.Phi();
        double eta_epd = StraightLine.Eta();
        //cout<<"phi_epd "<<phi_epd<<" eta_epd "<<eta_epd <<endl;
        epdeta[trig_index][direction_index]->Fill(eta_epd);
      }
      epdadc[trig_index][direction_index]->Fill(epdadc_e, epdadc_w);

      float ranking = vertex->ranking();
      int nutrks = vertex->nTracksUsed();
      int notrks = vertex->noTracks();
      int nprim = StMuDst::numberOfPrimaryTracks();
      int nglob = StMuDst::numberOfGlobalTracks();

      int numtrk = 0;
      double leading_pt=0.;
      double leading_eta=0.;
      double leading_phi=0.;
      double subleading_pt=0.;
      double subleading_eta=0.;
      double subleading_phi=0.;

      for(int iTrack = 0; iTrack < nprim; iTrack++){
        const StMuTrack* muTrack = StMuDst::primaryTracks(iTrack);
        int id = muTrack->id();
        int nhits = muTrack->nHits();
        int nhitsfit = muTrack->nHitsFit();
        int nhitsposs = muTrack->nHitsPoss();
        bool type = muTrack->type();
        bool flag = muTrack->flag();
        float pt = muTrack->pt();
        float phi= muTrack->phi();
        double eta = muTrack->eta();

        if(pt<0.2)continue;
        if(nhits<15)continue;

        if(pt>leading_pt){
          subleading_pt=leading_pt;
          subleading_eta=leading_eta;
          subleading_phi=leading_phi;

          leading_pt=pt;
          leading_eta=eta;
          leading_phi=phi;
        }
        else if(pt>subleading_pt){
          subleading_pt=pt;
          subleading_eta=eta;
          subleading_phi=phi;
        }
        htrketa[trig_index][direction_index]->Fill(eta);
        numtrk++;
      }
      hntrk[trig_index][direction_index]->Fill(numtrk);
      hvtxz[trig_index][direction_index]->Fill(vertex->position().z());
      pt2D[trig_index][direction_index]->Fill(leading_pt, subleading_pt);
      eta2D[trig_index][direction_index]->Fill(leading_eta, subleading_eta);
      phi2D[trig_index][direction_index]->Fill(leading_phi, subleading_phi);
    }
  }

  fout->Write();
  fout->Close();
  return 0;
}



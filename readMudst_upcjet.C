const unsigned long nevents = 400;
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
  TH1F* htrketa[num_trgs][direction];
  TH1F *epdeta[num_trgs][direction];
  TH1F *zdcadc_east[num_trgs][direction];
  TH1F *zdcadc_west[num_trgs][direction];
  TH1F *eemc[num_trgs][direction];
  
  TH2F *zdcadc[num_trgs][direction];
  TH2F *bbcadc[num_trgs][direction];
  TH2F *vpdadc[num_trgs][direction];
  TH2F *epdadc[num_trgs][direction];

  for(int i=0;i<num_trgs;i++){
    for(int j=0;j<direction;j++){
      hntrk[i][j] = new TH1F(Form("hntrk_%d_%d",i,j),";N_{trk};counts",200,0,200);
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
      TH2F *zdcadc[i][j] = new TH2F(Form("zdcadc_%d_%d",i,j),";ZDCE ADC; ZDCW ADC;counts",256,0,2000,256,0,2000);
      //bbc ADC
      TH2F *bbcadc[i][j] = new TH2F(Form("bbcadc_%d_%d",i,j),";BBCE ADC; BBCW ADC;counts",256,0,25000,256,0,25000);
      //vpd adc
      TH2F *vpdadc[i][j] = new TH2F(Form("vpdadc_%d_%d",i,j),";VPDE ADC; VPDW ADC;counts",256,0,4000,256,0,4000);
      //epd adc
      TH2F *epdadc[i][j] = new TH2F(Form("epdadc_%d_%d",i,j),";EPDE ADC; EPDW ADC;counts",256,0,40000,256,0,40000);
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

    //veto ZDC and EEMC
    StZdcTriggerDetector &muZdc = StMuDst::event()->zdcTriggerDetector();
    StEmcTriggerDetector &muemc = StMuDst::event()->emcTriggerDetector();

    //check BBC, EPD, VPD
    StBbcTriggerDetector &muBbc = StMuDst::event()->bbcTriggerDetector();
    StMuPrimaryVertex *vertex = StMuDst::primaryVertex();
    StEpdGeom * mEpdGeom = new StEpdGeom();

    //trigger decisions
    mMuDst = muDstMaker->muDst();
    StMuEvent *mMuEvent = mMuDst->event();
    //loop over triggers
    int trig_index=-1;
    for(int i=0;i<num_trgs;i++){
      if( mMuEvent->triggerIdCollection().nominal().isTrigger(commTriggerID[i]) ) trig_index=i;
    }
    if(trig_index<0) continue; //protection

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
    if(zdcadc_w >= 350){ZDCW=1;}
    if(zdcadc_e >= 350){ZDCE=1;}

    //1) now decide direction based on BBC first.
    int direction_index=-1;
    if(BBCE && !BBCW) direction_index=0;
    if(!BBCE && BBCW) direction_index=1;
    if(direction_index<0) continue; //protection

    //2) filling ZDC distributions based on BBC configurations
    zdcadc_east[trig_index][direction_index]->Fill(zdcadc_e);
    zdcadc_west[trig_index][direction_index]->Fill(zdcadc_w);
    zdcadc[trig_index][direction_index]->Fill(zdcadc_e,zdcadc_w);
    
    //3) now add ZDC to BBC requirements
    direction_index=-1;
    if(BBCE && !BBCW && ZDCE && !ZDCW) direction_index=0;
    if(!BBCE && BBCW && !ZDCE && ZDCW) direction_index=1;
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

    //vertex selections:
    if(!vertex) continue;
    TVector3 pRcVx(vertex->position().x(), vertex->position().y(), vertex->position().z());
    //Tof or BEMC match
    if(vertex->nBTOFMatch()<1 && vertex->nBEMCMatch()<1)continue;

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
    //Printf("ranking=%lg nutrks=%d nprim=%d notrks=%d nglob=%d\n", ranking, nutrks, nprim, notrks, nglob);

    int numtrk = 0;
    for(int iTrack = 0; iTrack < nprim; iTrack++){
      const StMuTrack* muTrack = StMuDst::primaryTracks(iTrack);
      int id = muTrack->id();
      int nhits = muTrack->nHits();
      int nhitsfit = muTrack->nHitsFit();
      int nhitsposs = muTrack->nHitsPoss();
      bool type = muTrack->type();
      bool flag = muTrack->flag();
      float pt = muTrack->pt();
      if(pt<0.2)continue;
      double eta = muTrack->eta();
      htrketa[trig_index][direction_index]->Fill(eta);
      numtrk++;
    }
    hntrk[trig_index][direction_index]->Fill(numtrk);

  }

  fout->Write();
  fout->Close();
  return 0;
}



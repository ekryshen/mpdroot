// E.Andronov 2020
// adaptation of miniDstProcExample code for fluctuations analysis needs
// this code allows us to analyze miniDST files
// and extract a tree with the list of observables related to fluctuations
// such as multiplicity, sum of total transverse momenta etc.


// C++ headers
#include <vector>
// ROOT headers
#include "Rtypes.h"
#include "TChain.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TDatabasePDG.h"
#include <TMath.h>
#include <TROOT.h>
#include <TTree.h>

// MpdRoot part
R__ADD_INCLUDE_PATH($VMCWORKDIR)
#include "macro/mpd/mpdloadlibs.C"

//
// inFineName can take either filename.MiniDst.root file or
// a list of such files with .list or .lis extention (for example,
// somelist.list or somelist.lis)
//

//_________________
void miniDstProcExample(const Char_t* inFileName,char * outname = " ") {

  //TFile out("/eos/nica/mpd/users/evandron/output_smash.root","recreate");
  TFile out(outname,"recreate");

  TH2F* gMomXHist = new TH2F("gMomXHist","gMomXHist;gMomX;pMomX;entries",400,-20,20,400,-20,20);
  TH2F* gMomYHist = new TH2F("gMomYHist","gMomYHist;gMomY;pMomY;entries",400,-20,20,400,-20,20);
  TH2F* gMomZHist = new TH2F("gMomZHist","gMomZHist;gMomZ;pMomZ;entries",400,-20,20,400,-20,20);

  TH1F* etaGHist = new TH1F("etaGHist","etaGHist;#eta_{global};entries",1000,-10,10);
  TH1F* etaNonPHist = new TH1F("etaNonPHist","etaNonPHist;#eta_{not_primary};entries",1000,-10,10);
  TH1F*	etaPHist = new TH1F("etaPHist","etaPHist;#eta_{primary};entries",1000,-10,10);
  TH1F* etaMCHist = new TH1F("etaMCHist","etaMCHist;#eta_{mc};entries",1000,-10,10);

  TH2F* dcaXYHist = new TH2F("dcaXYHist","dcaXYHist;dca_{x};dca_{y};enries",1000,-100,100,1000,-100,100);
  TH1F* dcaZHist = new TH1F("dcaZHist","dcaZHist;dca_{z};enries",3000,-100,100);

  TH2F* dcaXYNonPHist = new TH2F("dcaXYNonPHist","dcaXYNonPHist;dca_{x};dca_{y};enries",1000,-100,100,1000,-100,100);
  TH1F* dcaZNonPHist = new TH1F("dcaZNonPHist","dcaZNonPHist;dca_{z};enries",3000,-300,300);

  TH2F* vertexXYHist = new TH2F("vertexXYHist","vertexXYHist;vtx_{x};vtx_{y};enries",1000,-1,1,1000,-1,1);
  TH1F* vertexZHist = new TH1F("vertexZHist","vertexZHist;vtx_{z};enries",3000,-300,300);

  TH1F* vertexZFFDHist = new TH1F("vertexZFFDHist","vertexZFFDHist;vtx_{z};enries",3000,-300,300);


  TH2F* vertexXYMCHist = new TH2F("vertexXYMCHist","vertexXYMCHist;vtx_{x,MC};vtx_{y,MC};enries",1000,-1,1,1000,-1,1);
  TH1F* vertexZMCHist = new TH1F("vertexZMCHist","vertexZMCHist;vtx_{z,MC};enries",3000,-300,300);

  TH1F* impParHist = new TH1F("impParHist","impParHist;b [fm];enries",300,-0.5,29.5);

  TH1F* reactionPlaneHist = new TH1F("reactionPlaneHist","reactionPlaneHist;#Phi [rad];enries",600,-6.5,6.5);


  TH1F* gMomXResHist = new TH1F("gMomXResHist","gMomXResHist;#Delta gMomX;entries",400,-20,20);
  TH1F* pMomXResHist = new TH1F("pMomXResHist","pMomXResHist;#Delta pMomX;entries",400,-20,20);

  TH1F* gMomYResHist = new TH1F("gMomYResHist","gMomYResHist;#Delta gMomY;entries",400,-20,20);
  TH1F* pMomYResHist = new TH1F("pMomYResHist","pMomYResHist;#Delta pMomY;entries",400,-20,20);

  TH1F* gMomZResHist = new TH1F("gMomZResHist","gMomZResHist;#Delta gMomZ;entries",400,-20,20);
  TH1F* pMomZResHist = new TH1F("pMomZResHist","pMomZResHist;#Delta pMomZ;entries",400,-20,20);


  int _nPtBins = 500;
  double ptMin = 0;
  double ptMax = 5;
  TH1D* h_Pt_P = new TH1D("h_Pt", "h_Pt", _nPtBins, ptMin, ptMax);       // create your histogra

  TH1D* h_Phi_P = new TH1D("h_Phi", "h_Phi", 400, -2*TMath::TwoPi(), 2*TMath::TwoPi());


  TH1D* h_nNumberOfHitsTPC = new TH1D("h_nNumberOfHitsTPC", "h_nNumberOfHitsTPC", 201, -0.5, 200.5);

  int LogBinsX = 581;
  double StartValueX = 0.1;//
  double StopValueX = 20;//
  double LogWidthX[LogBinsX+1];

  int LogBinsY = 581;
  double StartValueY = 1000;//
  double StopValueY = 20000;//
  double LogWidthY[LogBinsY+1];

  //calculate bins
  for(int i = 0; i <= LogBinsX; i++) LogWidthX[i] = pow(10,log10(StartValueX)+(log10(StopValueX)-log10(StartValueX))/double(LogBinsX)*double(i));
  for(int i = 0; i <= LogBinsY; i++) LogWidthY[i] = pow(10,log10(StartValueY)+(log10(StopValueY)-log10(StartValueY))/double(LogBinsY)*double(i));


  TH1D*  h_dEdx = new TH1D("h_dEdx","h_dEdx", LogBinsY, LogWidthY );
  TH2D*  h_dEdx_vs_mom_all = new TH2D("h_dEdx_vs_mom","h_dEdx_vs_mom", LogBinsX, LogWidthX,LogBinsY, LogWidthY);


  Float_t ZDC_energy_mpd[90];
  TH1F* zdcHist = new TH1F("zdcHist","zdcHist;zdc;entries",400,-0.5,30);
  TH2F* zdcCorrHist = new TH2F("zdcCorrHist","zdcCorrHist;zdc_{1};zdc_{2};",400,-0.5,15,400,-0.5,15);
  TH2F* zdcImpHist = new TH2F("zdcImpHist","zdcImpHist;zdc;b;",400,-0.5,30,300,-0.5,29.5);


  TH2F* zdcMultHist = new TH2F("zdcMultHist","zdcMultHist;zdc;N;",400,-0.5,30,1000,-0.5,999.5);
  TH2F* impMultHist = new TH2F("impMultHist","impMultHist;b;N;",400,-0.5,29.5,1000,-0.5,999.5);




  // ##### EVENT TREE:
  static const Int_t NMaxTrack = 12000;

  Int_t fEventId; //!
  Float_t fImpactPar;   //!
  Float_t fReactionPl;   //!
  Float_t fVertexX,fVertexY,fVertexZ; //!
  Float_t fVertexXMC,fVertexYMC,fVertexZMC; //!
  Float_t f_zdc_total_energy; //!

  gROOT->cd();
  TTree       *fEventTree = new TTree("tree_sim_rec","tree_sim_rec"); //! output Tree
  fEventTree->Branch("fEventId",&fEventId,"fEventId/I");
  fEventTree->Branch("fImpactPar",&fImpactPar,"impPar/F");
  fEventTree->Branch("fReactionPl",&fReactionPl,"fReactionPl/F");

  fEventTree->Branch("vertexX",&fVertexX,"vertexX/F");
  fEventTree->Branch("vertexY",&fVertexY,"vertexY/F");
  fEventTree->Branch("vertexZ",&fVertexZ,"vertexZ/F");
  fEventTree->Branch("vertexXMC",&fVertexXMC,"vertexXMC/F");
  fEventTree->Branch("vertexYMC",&fVertexYMC,"vertexYMC/F");
  fEventTree->Branch("vertexZMC",&fVertexZMC,"vertexZMC/F");

  fEventTree->Branch("zdc_total_energy", &f_zdc_total_energy,"zdc_total_energy/F");

  // SIM:
  Int_t f_SIM_NumberOfTracks; //!
  Int_t f_SIM_TrackId[NMaxTrack];
  Float_t f_SIM_TrackPt[NMaxTrack],f_SIM_TrackEta[NMaxTrack],f_SIM_TrackPhi[NMaxTrack], f_SIM_TrackY[NMaxTrack]; //
  Int_t f_SIM_TrackCharge[NMaxTrack];
  Int_t f_SIM_pid[NMaxTrack];
  Int_t f_SIM_fromGen[NMaxTrack];

  fEventTree->Branch("SIM_NumberOfTracks",&f_SIM_NumberOfTracks,"SIM_NumberOfTracks/I");
  fEventTree->Branch("SIM_TrackId",f_SIM_TrackId,"SIM_TrackId[SIM_NumberOfTracks]/I");

  fEventTree->Branch("SIM_TrackPt",f_SIM_TrackPt,"SIM_TrackPt[SIM_NumberOfTracks]/F");
  fEventTree->Branch("SIM_TrackPhi",f_SIM_TrackPhi,"SIM_TrackPhi[SIM_NumberOfTracks]/F");
  fEventTree->Branch("SIM_TrackEta",f_SIM_TrackEta,"SIM_TrackEta[SIM_NumberOfTracks]/F");
  fEventTree->Branch("SIM_TrackY",f_SIM_TrackY,"SIM_TrackY[SIM_NumberOfTracks]/F");
  fEventTree->Branch("SIM_TrackCharge",f_SIM_TrackCharge,"SIM_TrackCharge[SIM_NumberOfTracks]/I");
  fEventTree->Branch("SIM_pid",f_SIM_pid,"SIM_pid[SIM_NumberOfTracks]/I");
  fEventTree->Branch("SIM_fromGen",f_SIM_fromGen,"SIM_fromGen[SIM_NumberOfTracks]/I");

  // REC:
  Int_t fNumberOfTracks; //!

  Float_t fTrackPt[NMaxTrack],fTrackEta[NMaxTrack],fTrackPhi[NMaxTrack];
  Float_t fTrackPtPrimary[NMaxTrack],fTrackEtaPrimary[NMaxTrack],fTrackPhiPrimary[NMaxTrack];
  Int_t fTrackCharge[NMaxTrack];
  Float_t fTrackMom[NMaxTrack];
  Float_t fTrackMomPrimary[NMaxTrack];
  Float_t fdEdx[NMaxTrack];
  Float_t fChi2[NMaxTrack];
  Float_t fDCA_x[NMaxTrack];
  Float_t fDCA_y[NMaxTrack];
  Float_t fDCA_z[NMaxTrack];
  Int_t fnTPChits[NMaxTrack];
  Int_t fHasPrimaryTrack[NMaxTrack];

  fEventTree->Branch("nTracks",&fNumberOfTracks,"nTracks/I");

  fEventTree->Branch("trackPt",fTrackPt,"trackPt[nTracks]/F");
  fEventTree->Branch("trackPhi",fTrackPhi,"trackPhi[nTracks]/F");
  fEventTree->Branch("trackEta",fTrackEta,"trackEta[nTracks]/F");
  fEventTree->Branch("trackPtPrimary",fTrackPtPrimary,"trackPtPrimary[nTracks]/F");
  fEventTree->Branch("trackPhiPrimary",fTrackPhiPrimary,"trackPhiPrimary[nTracks]/F");
  fEventTree->Branch("trackEtaPrimary",fTrackEtaPrimary,"trackEtaPrimary[nTracks]/F");

  fEventTree->Branch("trackCharge",fTrackCharge,"trackCharge[nTracks]/I");
  fEventTree->Branch("trackMom",fTrackMom,"trackMom[nTracks]/F");
  fEventTree->Branch("trackMomPrimary",fTrackMomPrimary,"trackMomPrimary[nTracks]/F");
  fEventTree->Branch("dEdx",fdEdx,"dEdx[nTracks]/F");
  fEventTree->Branch("chi2",fChi2,"chi2[nTracks]/F");

  fEventTree->Branch("DCA_x",fDCA_x,"DCA_x[nTracks]/F");
  fEventTree->Branch("DCA_y",fDCA_y,"DCA_y[nTracks]/F");
  fEventTree->Branch("DCA_z",fDCA_z,"DCA_z[nTracks]/F");

  fEventTree->Branch("nTPChits",fnTPChits,"fnTPChits[nTracks]/I");
  fEventTree->Branch("hasPrimaryTrack",fHasPrimaryTrack,"fHasPrimaryTrack[nTracks]/I");




  // Instantiate reader
  MpdMiniDstReader* miniDstReader = new MpdMiniDstReader(inFileName);

  // Reader initialization
  miniDstReader->Init();

  // One can specify branches to read
  miniDstReader->SetStatus("*", 0);               // Turn off all branches
  miniDstReader->SetStatus("Event*", 1);          // Turn on specific branch
  miniDstReader->SetStatus("Track*", 1);
  miniDstReader->SetStatus("BTofHit*", 1);
  miniDstReader->SetStatus("BTofPidTraits*", 1);
  miniDstReader->SetStatus("BECalCluster*", 1);
  miniDstReader->SetStatus("FHCalHit*", 1);
  miniDstReader->SetStatus("TrackCovMatrix*", 0); // Turn off specific branch
  miniDstReader->SetStatus("McEvent*", 1);
  miniDstReader->SetStatus("McTrack*", 1);

  // Retrieve events in tree
  Long64_t events2read = miniDstReader->chain()->GetEntries();

  // Loop over events
  for (Long64_t i = 0; i < events2read; i++) {


    if (i%10000==0) {
      cout << "event #" << i << endl;
    }

    // Read next event
    Bool_t  isOk = miniDstReader->readMiniEvent( i );

    // Retrieve current miniDst (from the given .MiniDst.root file)
    MpdMiniDst *dst = miniDstReader->miniDst();

    //
    // Event information
    //

    // Get MiniEvent
    MpdMiniEvent *event = dst->event();
    if (!event)
      continue;

    Int_t evID = event->eventId();


    // Get primary vertex z-position
    Float_t x = event->primaryVertex().X();
    Float_t y = event->primaryVertex().Y();
    Float_t z = event->primaryVertex().Z();
    //Float_t zFFD = event->vzFfd();

    if (event->numberOfGlobalTracks()<2) {
      continue;
    }

    if (TMath::Abs(z)>20) {
      continue;
    }

    if (TMath::Abs(z)<0.01) {
      continue;
    }


    fVertexX = x;
    fVertexY = y;
    fVertexZ = z;

    vertexXYHist->Fill(x,y);
    vertexZHist->Fill(z);
    //vertexZFFDHist->Fill(zFFD);
    Float_t total_ZDC = 0;
    Float_t total_ZDC_1 = 0;
    Float_t total_ZDC_2 = 0;
    for (Int_t j = 0; j < 90; j++)
    {
        ZDC_energy_mpd[j] = event->fhcalEnergyDeposit(j);
        total_ZDC+=ZDC_energy_mpd[j];
        if (j<45) {
            total_ZDC_1+=ZDC_energy_mpd[j];
        } else {
            total_ZDC_2+=ZDC_energy_mpd[j];
        }
    }

    f_zdc_total_energy = total_ZDC;
    zdcHist->Fill(total_ZDC);

    zdcCorrHist->Fill(total_ZDC_1,total_ZDC_2);
    //
    // Track information
    //

    // Retrieve number of reconstructed tracks
    Int_t nGTracks = dst->numberOfTracks();

    Float_t primaryTrackCounter = 0;


    // ### REC track loop
    fNumberOfTracks = 0;
    // Track loop
    for (Int_t j = 0; j < nGTracks; j++) {

      // Retrieve j-th mini track
      MpdMiniTrack *miniTrack = dst->track(j);

      // Global track full momentum
      Float_t ptot = miniTrack->gMom().Mag();
      Float_t px = miniTrack->gMom().X();
      Float_t py = miniTrack->gMom().Y();
      Float_t pt = TMath::Sqrt(px*px+py*py);
      Float_t pz = miniTrack->gMom().Z();
      Float_t phi = TMath::ATan2(py,px);
      Float_t eta_global = 0.5 * TMath::Log( ( ptot + pz )/( ptot - pz ) );

      Float_t dedx = miniTrack->dEdx();
      Float_t chi2 = miniTrack->chi2();
      Int_t charge = miniTrack->charge();
      h_nNumberOfHitsTPC->Fill(miniTrack->nHits());
      Float_t nHits = miniTrack->nHits();
      if (nHits<30.0) {
        continue;
      }
      Float_t chi2perCluster = chi2/nHits;
      if (chi2perCluster>5.0) {
        continue;
      }

      //if (miniTrack->nHits() < 20) continue;
      Float_t dcaX = miniTrack->gDCAx(x);
      Float_t dcaY = miniTrack->gDCAy(y);
      Float_t dcaZ = miniTrack->gDCAz(z);




      //dcaXYHist->Fill(dcaX,dcaY);
      //dcaZHist->Fill(dcaZ);

      //if (dcaX*dcaX+dcaY*dcaY>6.25) {
      //  continue;
      //}

      //if (dcaZ>50.0) {
      //  continue;
      //}




      Float_t ptot_p = miniTrack->pMom().Mag();
      Float_t px_p = miniTrack->pMom().X();
      Float_t py_p = miniTrack->pMom().Y();
      Float_t pt_p = TMath::Sqrt(px_p*px_p+py_p*py_p);
      Float_t pz_p = miniTrack->pMom().Z();
      Float_t phi_p = TMath::ATan2(py_p,px_p);
      Float_t eta_primary = 0.5 * TMath::Log( ( ptot_p + pz_p )/( ptot_p - pz_p ) );
      if (ptot_p!=0){
        //if ( TMath::Abs(eta_primary) > 1.6 ) continue;

        //if (pt_p>2 || pt_p<0.2) continue;


        if (TMath::Abs(eta_primary) < 1.0 && pt_p>0.15 && pt_p<2.0 && dcaX*dcaX+dcaY*dcaY<6.25) {
          primaryTrackCounter = primaryTrackCounter + 1;
          dcaXYHist->Fill(dcaX,dcaY);
          dcaZHist->Fill(dcaZ);
        }

        etaPHist->Fill(eta_primary);
        gMomXHist->Fill(miniTrack->gMom().X(), miniTrack->pMom().X());
        gMomYHist->Fill(miniTrack->gMom().Y(), miniTrack->pMom().Y());
        gMomZHist->Fill(miniTrack->gMom().Z(), miniTrack->pMom().Z());

        Float_t mom_to_charge = fabs(ptot_p/charge);

        h_dEdx->Fill(dedx);
        h_dEdx_vs_mom_all->Fill( mom_to_charge, dedx );
        h_Pt_P->Fill(pt_p);
        h_Phi_P->Fill(phi_p);



        etaGHist->Fill(eta_global);


        fTrackPt[fNumberOfTracks] = pt;
        fTrackPhi[fNumberOfTracks] = phi;
        fTrackEta[fNumberOfTracks] = eta_global;

        fTrackPtPrimary[fNumberOfTracks] = pt_p;
        fTrackPhiPrimary[fNumberOfTracks] = phi_p;
        fTrackEtaPrimary[fNumberOfTracks] = eta_primary;


        fTrackCharge[fNumberOfTracks] = charge;
        fTrackMom[fNumberOfTracks] = ptot;
        fTrackMomPrimary[fNumberOfTracks] = ptot_p;
        fdEdx[fNumberOfTracks] = dedx;
        fChi2[fNumberOfTracks] = chi2;

        fDCA_x[fNumberOfTracks] = dcaX;
        fDCA_y[fNumberOfTracks] = dcaY;
        fDCA_z[fNumberOfTracks] = dcaZ;

        fnTPChits[fNumberOfTracks] = miniTrack->nHits();


        fHasPrimaryTrack[fNumberOfTracks] = 1;

        fNumberOfTracks++;



      }else{
        etaNonPHist->Fill(eta_global);
        dcaXYNonPHist->Fill(dcaX,dcaY);
        dcaZNonPHist->Fill(dcaZ);


        fTrackPt[fNumberOfTracks] = pt;
        fTrackPhi[fNumberOfTracks] = phi;
        fTrackEta[fNumberOfTracks] = eta_global;

        fTrackPtPrimary[fNumberOfTracks] = 0;
        fTrackPhiPrimary[fNumberOfTracks] = 0;
        fTrackEtaPrimary[fNumberOfTracks] = 0;


        fTrackCharge[fNumberOfTracks] = charge;
        fTrackMom[fNumberOfTracks] = ptot;
        fTrackMomPrimary[fNumberOfTracks] = 0;
        fdEdx[fNumberOfTracks] = dedx;
        fChi2[fNumberOfTracks] = chi2;

        fDCA_x[fNumberOfTracks] = dcaX;
        fDCA_y[fNumberOfTracks] = dcaY;
        fDCA_z[fNumberOfTracks] = dcaZ;

        fnTPChits[fNumberOfTracks] = miniTrack->nHits();

        fHasPrimaryTrack[fNumberOfTracks] = 0;

        fNumberOfTracks++;
        //Bool_t testMCCandidate = miniTrack->hasMcTrack();
        //if (testMCCandidate == true) {
        //    Int_t nonPrimaryMCCandidate = miniTrack->mcTrackIndex();
        //    MpdMiniMcTrack *mcCandidateTrack = dst->mcTrack( nonPrimaryMCCandidate );
        //    Int_t pdgId_MCCandidate = mcCandidateTrack->pdgId();
        //    if (TDatabasePDG::Instance()->GetParticle(pdgId_MCCandidate)==NULL) {
        //      continue;
        //    }
        //    cout << "non primary pdg = " << pdgId_MCCandidate << endl;
        //    cout << "isFromGenerator() = " << mcCandidateTrack->isFromGenerator() << endl;
        //}

      }



      //gMomXHist->Fill(miniTrack->gMom().X(), miniTrack->pMom().X());
      //gMomYHist->Fill(miniTrack->gMom().Y(), miniTrack->pMom().Y());
      //gMomZHist->Fill(miniTrack->gMom().Z(), miniTrack->pMom().Z());

      //
      // TOF-matched tracks
      //

      // Check if track matched TOF
      if ( miniTrack->isBTofTrack() ) {

	// Retrieve TOF information for the matched track
	MpdMiniBTofPidTraits *tofTrait = dst->btofPidTraits( miniTrack->bTofPidTraitsIndex() );

	// Retrive beta
	Float_t beta = tofTrait->beta();
      }
    }

    //cout << "fNumberOfTracks = " << fNumberOfTracks << " with " << primaryTrackCounter << "primaries" << endl;



      zdcMultHist->Fill(total_ZDC,primaryTrackCounter);

    //
    // Hit information
    //

    // Loop over barrel TOF hits
    for (Int_t j = 0; j < dst->numberOfBTofHits(); j++ ) {

      // Retrieve j-th hit information
      MpdMiniBTofHit *btofHit = dst->btofHit( j );

      // x position of hit
      Float_t xPosition = btofHit->btofHitPosX();
    }

    // Loop over FHCal hits
    for (Int_t j = 0; j < dst->numberOfFHCalHits(); j++) {

      // Retrieve j-th hit information
      MpdMiniFHCalHit *hit = dst->fhcalHit( j );

      // Module number that corresponds to the hit
      Int_t modId = hit->module();
    }

    //
    // BECal cluster information
    //

    // Loop over clusters
    for (Int_t j = 0; j < dst->numberOfBECalClusters(); j++) {

      // Retrieve j-th cluster
      MpdMiniBECalCluster *cluster = dst->becalCluster( j );

      // Retrieve cluster energy
      Float_t energy = cluster->energy();

      // Retrieve mcTrack IDs that produced cluster
      std::vector<UShort_t> mcTrackIds = cluster->mcTrackIds();
    }

    //
    // Monte Carlo information
    //

    // Retrieve MC event
    MpdMiniMcEvent *mcEvent = dst->mcEvent();

    Float_t x_MC = mcEvent->primaryVertex().X();
    Float_t y_MC = mcEvent->primaryVertex().Y();
    Float_t z_MC = mcEvent->primaryVertex().Z();

    vertexXYMCHist->Fill(x_MC,y_MC);
    vertexZMCHist->Fill(z_MC);

    // Impact parameter
    Float_t b = mcEvent->primaryVertexZ();

    impParHist->Fill(mcEvent->impactParameter());
    reactionPlaneHist->Fill(mcEvent->reactionPlaneAngle());

    zdcImpHist->Fill(total_ZDC, mcEvent->impactParameter());

    impMultHist->Fill(mcEvent->impactParameter(), primaryTrackCounter);

    UInt_t evMCID = mcEvent->eventId();

    fEventId = (Int_t) evMCID;
    fImpactPar = mcEvent->impactParameter();
    fReactionPl = mcEvent->reactionPlaneAngle();
    fVertexXMC = x_MC;
    fVertexYMC = y_MC;
    fVertexZMC = z_MC;

    f_SIM_NumberOfTracks = 0;
    Int_t f_SIM_NumberOfTracks_fromGen = 0;
    // MC track loop
    for (Int_t j = 0; j < dst->numberOfMcTracks(); j++) {

      // Retrieve j-th MC track
      MpdMiniMcTrack *mcTrack = dst->mcTrack( j );

      UShort_t mcID = mcTrack->id();
      Bool_t fromGen = mcTrack->isFromGenerator();


      Float_t ptot_MC = mcTrack->p().Mag();
      Float_t px_MC = mcTrack->p().X();
      Float_t py_MC = mcTrack->p().Y();
      Float_t pt_MC = TMath::Sqrt(px_MC*px_MC+py_MC*py_MC);
      Float_t pz_MC = mcTrack->p().Z();
      Float_t e_MC = mcTrack->e();
      Float_t phi_MC = TMath::ATan2(py_MC,px_MC);
      Float_t eta_MC = 0.5 * TMath::Log( ( ptot_MC + pz_MC )/( ptot_MC - pz_MC ) );
      Float_t rapidity_MC = 0.5 * TMath::Log( ( e_MC + pz_MC )/( e_MC - pz_MC ) );

      Int_t pdgId_MC = mcTrack->pdgId();
      if (TDatabasePDG::Instance()->GetParticle(pdgId_MC)==NULL) {
        continue;
      }


      Float_t charge_MC = (mcTrack->charge())/3.0;
      if (charge_MC==0) {
        continue;
      }

      if (pt_MC<0.15 || pt_MC>2.0) {
        continue;
      }

      if (TMath::Abs(eta_MC) > 1.0) {
        continue;
      }
      //cout << fromGen << " pdg = " << pdgId_MC << endl;

      if (fromGen!=0) {
        f_SIM_NumberOfTracks_fromGen++;
      }

        //if (pt_MC>=0.2 && pt_MC<=2 && TMath::Abs(eta_MC) <= 1.6) {
            // ### FOR TREE:
            f_SIM_TrackId[f_SIM_NumberOfTracks] = (Int_t) mcID;
            f_SIM_TrackPt[f_SIM_NumberOfTracks] = pt_MC;
            f_SIM_TrackPhi[f_SIM_NumberOfTracks] = phi_MC;
            f_SIM_TrackEta[f_SIM_NumberOfTracks] = eta_MC;
            f_SIM_TrackY[f_SIM_NumberOfTracks] = rapidity_MC;
            f_SIM_TrackCharge[f_SIM_NumberOfTracks] = charge_MC;
            f_SIM_pid[f_SIM_NumberOfTracks] = pdgId_MC;
            f_SIM_fromGen[f_SIM_NumberOfTracks] = (Int_t) fromGen;
            f_SIM_NumberOfTracks++;
        //}




      // Print MC track information
      //mcTrack->Print();

      std::vector< UShort_t > itsRecTracks = mcTrack->recoTrackIds();
      Int_t sizeOfitsRecTracks = itsRecTracks.size();


      if (sizeOfitsRecTracks>0) {
        etaMCHist->Fill(eta_MC);
      }

      //if (sizeOfitsRecTracks>1) {
      //  cout << "double candidate" << endl;
      //  cout << "sizeOfitsRecTracks = " << sizeOfitsRecTracks << endl;
      //  cout << "charge_MC = " << charge_MC << endl;
      //}
      for (Int_t i=0; i<sizeOfitsRecTracks; i++) {
        MpdMiniTrack *miniItsTrack = dst->track(itsRecTracks.at(i));
        Int_t charge_itsTrack = miniItsTrack->charge();
        Int_t nHits_itsTrack = miniItsTrack->nHits();
        Float_t ptot_its = miniItsTrack->gMom().Mag();
        Float_t px_its = miniItsTrack->gMom().X();
        Float_t py_its = miniItsTrack->gMom().Y();
        Float_t pt_its = TMath::Sqrt(px_its*px_its+py_its*py_its);
        Float_t pz_its = miniItsTrack->gMom().Z();

        Float_t ptot_p_its = miniItsTrack->pMom().Mag();
        Float_t px_p_its = miniItsTrack->pMom().X();
        Float_t py_p_its = miniItsTrack->pMom().Y();
        Float_t pt_p_its = TMath::Sqrt(px_p_its*px_p_its+py_p_its*py_p_its);
        Float_t pz_p_its = miniItsTrack->pMom().Z();
        Float_t eta_p_its = 0.5 * TMath::Log( ( ptot_p_its + pz_p_its )/( ptot_p_its - pz_p_its ) );

        Float_t dcaX_its = miniItsTrack->gDCAx(x);
        Float_t dcaY_its = miniItsTrack->gDCAy(y);
        Float_t dcaZ_its = miniItsTrack->gDCAy(z);

        if (ptot_p_its!=0 && nHits_itsTrack>=20 && dcaZ_its<50.0 && dcaX_its*dcaX_its+dcaY_its*dcaY_its<6.25 && TMath::Abs(eta_p_its) <= 1.6 && pt_p_its>=0.2 && pt_p_its<2) {
          Float_t tmp_x = (px_MC-px_its)/px_MC;
          gMomXResHist->Fill(tmp_x);
          Float_t tmp_y = (py_MC-py_its)/py_MC;
          gMomYResHist->Fill(tmp_y);
          Float_t tmp_z = (pz_MC-pz_its)/pz_MC;
          gMomZResHist->Fill(tmp_z);

          Float_t tmp_p_x = (px_MC-px_p_its)/px_MC;
          pMomXResHist->Fill(tmp_p_x);
          Float_t tmp_p_y = (py_MC-py_p_its)/py_MC;
          pMomYResHist->Fill(tmp_p_y);
          Float_t tmp_p_z = (pz_MC-pz_p_its)/pz_MC;
          pMomZResHist->Fill(tmp_p_z);

          //if (charge_itsTrack*charge_MC<0) {
          //  cout << "wrong charge!" << endl;
          //}

          //if (sizeOfitsRecTracks>1) {
          //  cout << "i = " << i << endl;
          //  cout << "charge_itsTrack = " << charge_itsTrack << " ; " << nHits_itsTrack << endl;
          //  cout << tmp_x << " " << tmp_y << " " << tmp_z << " " << TMath::Sqrt(tmp_x*tmp_x + tmp_y*tmp_y + tmp_z*tmp_z) << endl;
          //  cout << tmp_p_x << " " << tmp_p_y << " " << tmp_p_z << " " << TMath::Sqrt(tmp_p_x*tmp_p_x + tmp_p_y*tmp_p_y + tmp_p_z*tmp_p_z) << endl;
          //}

        }
      }

    }
    //cout << "f_SIM_NumberOfTracks = " << f_SIM_NumberOfTracks << " ; among them " << f_SIM_NumberOfTracks_fromGen << " fromGen" << endl;
    fEventTree->Fill();

  }

  // Finalize miniDst reader
  miniDstReader->Finish();

  /*TCanvas* c = new TCanvas("c","c",1200,800);
  c->Divide(3,1);
  c->cd(1);
  gMomXHist->Draw("colz");
  c->cd(2);
  gMomYHist->Draw("colz");
  c->cd(3);
  gMomZHist->Draw("colz");
  c->SaveAs("gMomXYZHist_smash.pdf");

  TCanvas* c2 = new TCanvas("c2","c2",1200,800);
  c2->cd();
  c2->cd()->SetLogy();
  etaMCHist->SetLineColor(kGreen);
  etaMCHist->Draw("hist");
  etaGHist->Draw("same");
  etaPHist->SetLineColor(kRed);
  etaPHist->Draw("same");
  c2->SaveAs("eta_smash.pdf");

  cout << "etaGHist->GetEntries() = " << etaGHist->GetEntries() << endl;
  cout << "etaPHist->GetEntries() = " << etaPHist->GetEntries() << endl;
  cout << "etaMCHist->GetEntries() = " << etaMCHist->GetEntries() << endl;

  TCanvas* c3 = new TCanvas("c3","c3",1200,800);
  c3->Divide(2,1);
  c3->cd(1);
  dcaXYHist->Draw("colz");
  c3->cd(2);
  c3->cd(2)->SetLogy();
  dcaZHist->Draw("hist");
  c3->SaveAs("dca_smash.pdf");


  TCanvas* c4 = new TCanvas("c4","c4",1200,800);
  c4->Divide(2,1);
  c4->cd(1);
  vertexXYHist->Draw("colz");
  c4->cd(2);
  c4->cd(2)->SetLogy();
  vertexZHist->Draw("hist");
  //vertexZFFDHist->SetLineColor(kRed);
  //vertexZFFDHist->Draw("same");
  c4->SaveAs("vertex_smash.pdf");


  TCanvas* c5 = new TCanvas("c5","c5",1200,800);
  c5->Divide(2,1);
  c5->cd(1);
  vertexXYMCHist->Draw("colz");
  c5->cd(2);
  c5->cd(2)->SetLogy();
  vertexZMCHist->Draw("hist");
  c5->SaveAs("vertexMC_smash.pdf");

  TCanvas* c6 = new TCanvas("c6","c6",1200,800);
  c6->Divide(2,1);
  c6->cd(1);
  impParHist->Draw("hist");
  c6->cd(2);
  reactionPlaneHist->Draw("hist");
  c6->SaveAs("imppar_smash.pdf");

  TCanvas* c7 = new TCanvas("c7","c7",1200,800);
  c7->Divide(2,1);
  c7->cd(1);
  c7->cd(1)->SetLogy();
  gMomXResHist->Draw("hist");
  c7->cd(2);
  c7->cd(2)->SetLogy();
  pMomXResHist->Draw("hist");
  c7->SaveAs("resX_smash.pdf");

  TCanvas* c8 = new TCanvas("c8","c8",1200,800);
  c8->Divide(2,1);
  c8->cd(1);
  c8->cd(1)->SetLogy();
  gMomYResHist->Draw("hist");
  c8->cd(2);
  c8->cd(2)->SetLogy();
  pMomYResHist->Draw("hist");
  c8->SaveAs("resY_smash.pdf");

  TCanvas* c9 = new TCanvas("c9","c9",1200,800);
  c9->Divide(2,1);
  c9->cd(1);
  c9->cd(1)->SetLogy();
  gMomZResHist->Draw("hist");
  c9->cd(2);
  c9->cd(2)->SetLogy();
  pMomZResHist->Draw("hist");
  c9->SaveAs("resZ_smash.pdf");
  */
  out.cd();
  dcaXYHist->Write();
  dcaZHist->Write();
  /*impParHist->Write();
  reactionPlaneHist->Write();
  zdcHist->Write();
  zdcCorrHist->Write();
  zdcImpHist->Write();
  zdcMultHist->Write();
  impMultHist->Write();
  vertexXYHist->Write();
  vertexZHist->Write();
  vertexXYMCHist->Write();
  vertexZMCHist->Write();
  etaGHist->Write();
  etaNonPHist->Write();
  etaPHist->Write();
  dcaXYHist->Write();
  dcaZHist->Write();
  dcaXYNonPHist->Write();
  dcaZNonPHist->Write();
  h_Pt_P->Write();
  h_Phi_P->Write();
  h_nNumberOfHitsTPC->Write();
  h_dEdx->Write();
  h_dEdx_vs_mom_all->Write();

  gMomXResHist->Write();
  pMomXResHist->Write();
  gMomYResHist->Write();
  pMomYResHist->Write();
  gMomZResHist->Write();
  pMomZResHist->Write();
  */
  fEventTree->Write();
  out.Close();

}

#if !defined(__CINT__) || defined(__MAKECINT__)
// MPD includes
#include "TpcPoint.h"
#include "MpdEvent.h"
#include "MpdTrack.h"
#include "MpdHelix.h"
#include "MpdItsKalmanTrack.h"
#include "MpdTpcKalmanTrack.h"
#include "MpdKalmanHit.h"
#include "MpdKalmanFilter.h"
#include "MpdMotherFitterPart.h"
#include "MpdKfV0Fitter.h"
#include "MpdParticle.h"
#include "MpdTrackFinderIts.h"
#include "MpdVertex.h"

// CBM includes
#include "FairMCPoint.h"
#include "FairMCTrack.h"
#include "FairRunAna.h"

// ROOT includes
#include <TBranch.h>
#include <TChain.h>
#include <TClonesArray.h>
//#include <TDatabasePDG.h>
#include <TFile.h>
#include <TFolder.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TParticlePDG.h>
#include <TRandom.h>
#include <TROOT.h>
#include <TTree.h>
#include <TVector3.h>
#include <TMinuit.h>
#include <Riostream.h>

#include <set>
#include <map>
#include <vector>

#endif

Int_t pdgCodePos= 2212; // proton
Int_t pdgCodeNeg= -211; // pi-
Int_t pdgCodeL0 = 3122; // Lambda0 (1.11568)
Int_t pdgCodeK0 = 310; // K0s
Int_t pdgCodeOm = 3334; // Omega- (1.67245 GeV/c)
Int_t pdgCodeXi = 3312; // Xi- (1.3213 GeV/c)
Int_t pdgCodeKm = -321; // K-
//Int_t nITS = 0, idPlus = -1;
Int_t *lays = 0x0;
MpdHelix trC;
MpdVertex *mpdVert;
TVector3 vtxN, momN, primVert;
TClonesArray *itsTracks, *mcTracks;
TMinuit *gMinuit = new TMinuit(2);  //initialize TMinuit with a maximum of 2 params
FILE *lun = 0x0; //fopen("event.dat","w");
FILE *lun1 = 0x0; //fopen("ids.dat","w");

// AZ:

const Double_t gC2p      = 9.;           //chi2 of p to PV
const Double_t gC2pi     = 11.;          //chi2 of pion to PV
const Double_t gC2Lpv    = 0.;           //chi2 of L to PV
const Double_t gC2Xipv   = 0.;           //chi2 of Xi to PV

const Double_t gDCAp     = 0.;           //cm - DCA of p to PV
const Double_t gDCApi    = 0.;           //cm - DCA of pion to PV
const Double_t gDCAL     = 0.;           //cm - DCA of L to PV
const Double_t gDCAXi    = 0.;           //cm - DCA of Xi to PV

const Double_t gDistL    = 9999.;        //cm - DCA between pion & p in V0
const Double_t gDistXi   = 9999.;        //cm - DCA between pion & L in Xi
const Double_t gPathL    = 2.0;          //cm - path to Lambda decay
const Double_t gPathXi   = 0.;           //cm - path to Xi decay
const Double_t gC2L      = 16.;//9999.;  //chi2 between pion & p in V0
const Double_t gC2Xi     = 16.;//9999.;  //chi2 between pion & L in Xi

const Double_t gDcaL0 = 0.; //0.15; 
const Double_t gDcaK = 0.; //0.3;
const Double_t gChi2K = 0.; //100; //50;
const Double_t gDecayOm = 0.; //1.0;
const Double_t gDistLK = 9999.; //1.15; //0.2;
const Double_t gDcaOm = 9999.; //0.1; //0.15;

//#define ITS
#ifdef ITS
typedef MpdItsKalmanTrack AzTrack;
#else
typedef MpdTpcKalmanTrack AzTrack;
#endif
void BuildLambda(vector<Int_t> &vecP, vector<Int_t> &vecPi, vector<MpdParticle*> &vecL,
		 vector<Double_t> &vecL1);
void BuildCascade(vector<Int_t> &vecK, vector<Int_t> &vecPi, vector<MpdParticle*> &vecL,
		  vector<Double_t> &vecL1);
MpdHelix MakeHelix(const MpdKalmanTrack *tr);
MpdHelix MakeHelix(const MpdParticle *part);
Double_t DistHelLin(MpdKalmanTrack *helix, MpdParticle *neu);
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

Float_t massh, pth, etah, yh, chi2h, disth, path, c2pv, d2pv, masshL, chi2hL, disthL, pathL;
Float_t etas[2], pts[2], chi2s[2], dcas[2], probs[2], chi2sL[2];
Double_t *kProb, *piProb, *pProb, *eProb;
Int_t evNo, origs[2], qs[2], dstNo[2], layMx[2];
Int_t *id2dst;
vector<pair<Double_t,Double_t> > vecL2;

//__________________________________________________________________________
void AnalOm(Int_t n1 = 0, Int_t n2 = 0)
{
  // Analyze TPC (ITS) reco data - reconstruct hyperons

  // Add particles to PDG database
  //gROOT->ProcessLine(".x $VMCWORKDIR/macro/mpd/AddToPdg.C");
  
  // Load basic libraries
  gROOT->ProcessLine(".x ~/mpd/loadlibs.C");

  gROOT->ProcessLine(".x ./Chain.C(200,\"./mc_1.root\")");
  //gROOT->ProcessLine(".x ./Chain.C(10,\"../its/AuAu_9AGeV_central_tpc_1.root\")");
  //Chain(10,"./mc_1.root");
  TChain *simMC = (TChain*) gROOT->FindObject("mpdsim");
  simMC->SetName("mpdsim1");
  //cout << simMC->GetEntries() << endl;
  gROOT->ProcessLine(".x ./Chain.C(200,\"./reco_1.root\")");
  //gROOT->ProcessLine(".x ./Chain.C(10,\"../its/AuAu_9AGeV_central.reco_tpc_1.root\")");
  //Chain(10,"./reco_1.root");
  TChain *simITS = (TChain*) gROOT->FindObject("mpdsim");
  //simITS->ls();
  //cout << simITS->GetEntries() << endl;

  //TTree *simITS = (TTree*) fileITS.Get("mpdsim");
  //TTree *simMC = (TTree*) fileMC.Get("mpdsim");

  TFile fileITS(simITS->GetListOfFiles()->First()->GetTitle());
  TFile fileMC(simMC->GetListOfFiles()->First()->GetTitle());
  fileMC.Get("FairBaseParSet");
  TClonesArray *vtxs = (TClonesArray*) fileITS.FindObjectAny("Vertex");
  simITS->SetBranchAddress("Vertex",&vtxs);
  TBranch *vtxB = simITS->GetBranch("Vertex");
  itsTracks = (TClonesArray*) fileITS.FindObjectAny("ItsTrack");
  simITS->SetBranchAddress("ItsTrack",&itsTracks);
  TBranch *itsRecoB = simITS->GetBranch("ItsTrack");
  if (itsTracks == 0x0) {
    itsTracks = (TClonesArray*) fileITS.FindObjectAny("TpcKalmanTrack");
    simITS->SetBranchAddress("TpcKalmanTrack",&itsTracks);
    itsRecoB = simITS->GetBranch("TpcKalmanTrack");
  }
  MpdEvent *event = 0x0;
  simITS->SetBranchAddress("MPDEvent.", &event);
  TClonesArray *mpdTracks = 0x0;

  TClonesArray *tpcPoints = (TClonesArray*) fileMC.FindObjectAny("TpcPoint");
  simMC->SetBranchAddress("TpcPoint",&tpcPoints);
  TBranch *tpcSimB = simMC->GetBranch("TpcPoint");
  TClonesArray *itsPoints = (TClonesArray*) fileMC.FindObjectAny("StsPoint");
  simMC->SetBranchAddress("StsPoint",&itsPoints);
  TBranch *itsSimB = simMC->GetBranch("StsPoint");
  mcTracks = (TClonesArray*) fileMC.FindObjectAny("MCTrack");
  simMC->SetBranchAddress("MCTrack",&mcTracks);
  TBranch *mcBranch = simMC->GetBranch("MCTrack");

  TFile out("hyperons.histo_Xi.root","recreate");

  FairRunAna ana;
  MpdKalmanFilter::Instance("KF")->Init();
  MpdTrackFinderIts* recoIts = new MpdTrackFinderIts("ITS track finder");
  recoIts->FillGeoScheme();

  // Book histos

  TH1D *hLambFlag = new TH1D("hLambFlag","Flags for lambda",10,0,10);
  TH1D *hRecognitionFlag = new TH1D("hRecognitionFlag","Flags for Recognition",10,0,10);
  TH1D *hXiFlag = new TH1D("hXiFlag","Flags for Xi",10,0,10);

  TH1D *hMassL = new TH1D("hMassL","Lambda mass",50,1.070,1.170);
  TH1D *hMassLsig = new TH1D("hMassLsig","Lambda mass (signal)",50,1.070,1.170);
  TH1D *hMassLbkg = new TH1D("hMassLbkg","Lambda mass (bkg.)",50,1.070,1.170);

  TH1D *hMassXi = new TH1D("hMassXi","Xi- mass",50,1.260,1.360);
  TH1D *hMassXiSig = new TH1D("hMassXiSig","Xi- mass (signal)",50,1.260,1.360);
  TH1D *hMassXiBkg = new TH1D("hMassXiBkg","Xi- mass (bkg.)",50,1.260,1.360);

  TH1D *hMassOm = new TH1D("hMassOm","Omega- mass",50,1.630,1.730);
  TH1D *hMassOmSig = new TH1D("hMassOmSig","Omega- mass (signal)",50,1.630,1.730);
  TH1D *hMassOmBkg = new TH1D("hMassOmBkg","Omega- mass (bkg.)",50,1.630,1.730);
  TH1D *hPdg = new TH1D("hPdg","PdgCod if is not Pion",1000,-2500,2500);
  
  TH2D *hProbTrueP = new TH2D("hProbTrueP","Probability for true Protons",50,0,1,50,0,1.1);
  TH2D *hProbP = new TH2D("hProbfalseP","Probability for Pions and identification Protons",50,0,1.1,50,0,1.1);
  TH2D *hProbPi = new TH2D("hProbfalsePi","Probability for Protons and identification Pions",50,0,1.1,50,0,1.1);
  TH2D *hProbTruePi = new TH2D("hProbTruePi","Probability for true Pions",50,0,1.1,50,0,1.1);
  TH2D *hProbTrueK = new TH2D("hProbTrueK","Probability for true Kaons",50,0,1.1,50,0,1.1);
  
  Double_t pmom, eta1, dpp, rorig, ptt;
  Int_t prim, idtr, np, moth, pdg;
  /*
  TTree *tree = new TTree("tracks","Tracks");
  tree->Branch("mom",&pmom,"pmom/F");
  tree->Branch("ptt",&ptt,"ptt/F");
  tree->Branch("eta",&eta1,"eta1/F");
  tree->Branch("dpp",&dpp,"dpp/F");
  tree->Branch("rorig",&rorig,"rorig/F");
  tree->Branch("prim",&prim,"prim/I");
  tree->Branch("id",&idtr,"idtr/I");
  tree->Branch("evNo",&evNo,"evNo/I");
  tree->Branch("np",&np,"np/I");
  tree->Branch("moth",&moth,"moth/I");
  tree->Branch("pdg",&pdg,"pdg/I");
  */

  TTree *thyps = new TTree("hypers","Hyperons");
  thyps->Branch("massh",&massh,"massh/F");
  thyps->Branch("pth",&pth,"pth/F");
  thyps->Branch("etah",&etah,"etah/F");
  thyps->Branch("yh",&yh,"yh/F");
  thyps->Branch("chi2h",&chi2h,"chi2h/F");
  thyps->Branch("disth",&disth,"disth/F");
  thyps->Branch("path",&path,"path/F");
  thyps->Branch("etas[2]",&etas,"etas[2]/F");
  thyps->Branch("pts[2]",&pts,"pts[2]/F");
  thyps->Branch("chi2s[2]",&chi2s,"chi2s[2]/F");
  thyps->Branch("dcas[2]",&dcas,"dcas[2]/F");
  thyps->Branch("probs[2]",&probs,"probs[2]/F");
  thyps->Branch("origs[2]",&origs,"origs[2]/I2");
  thyps->Branch("qs[2]",&qs,"qs[2]/I2");
  thyps->Branch("dstNo[2]",&dstNo,"dstNo[2]/I");
  thyps->Branch("layMx[2]",&layMx,"layMx[2]/I2");
  thyps->Branch("evNo",&evNo,"evNo/I");
  
  TTree *txi = new TTree("Xi","xi");
  txi->Branch("massXi",&massh,"massh/F");
  txi->Branch("ptXi",&pth,"pth/F");
  txi->Branch("etaXi",&etah,"etah/F");
  txi->Branch("yXi",&yh,"yh/F");
  txi->Branch("chi2Xi",&chi2h,"chi2h/F");
  txi->Branch("distXi",&disth,"disth/F");
  txi->Branch("pathXi",&path,"path/F");
  txi->Branch("c2Xipv",&c2pv,"c2pv/F");
  txi->Branch("dcaXi",&d2pv,"d2pv/F");
  txi->Branch("etas[2]",&etas,"etas[2]/F");
  txi->Branch("pts[2]",&pts,"pts[2]/F");
  txi->Branch("chi2s[2]",&chi2s,"chi2s[2]/F");
  txi->Branch("dcas[2]",&dcas,"dcas[2]/F");
  txi->Branch("probs[2]",&probs,"probs[2]/F");
  txi->Branch("origs[2]",&origs,"origs[2]/I2");
  txi->Branch("qs[2]",&qs,"qs[2]/I2");
  txi->Branch("dstNo[2]",&dstNo,"dstNo[2]/I");
  txi->Branch("layMx[2]",&layMx,"layMx[2]/I2");
  txi->Branch("evNo",&evNo,"evNo/I");
  txi->Branch("massL",&masshL,"masshL/F");
  txi->Branch("chi2L",&chi2hL,"chi2hL/F");
  txi->Branch("distL",&disthL,"disthL/F");
  txi->Branch("pathL",&pathL,"pathL/F");
  txi->Branch("chi2sL[2]",&chi2sL,"chi2sL[2]/F");

  Int_t events = simITS->GetEntries();
  if (n2 != 0) events = TMath::Min (events, n2);
  cout << " Number of events = " << events << endl;

  for (Int_t i = 0; i < events; ++i) {
    if (i < n1) continue;
    simMC->GetEntry(i);
    simITS->GetEntry(i);
    evNo = i + 1;

    // For ITS points
    TVector3 mom; 
    set<Int_t> idxs;
    if (itsPoints) {
      Int_t nITSp = itsPoints->GetEntriesFast();
      for (Int_t j = 0; j < nITSp; ++j) {
	FairMCPoint *p1 = (FairMCPoint*) itsPoints->UncheckedAt(j);
	Int_t id = p1->GetTrackID();
	if (idxs.find(id) != idxs.end()) continue;
	//idxs.insert(id);
	FairMCTrack* mcTr = (FairMCTrack*) mcTracks->UncheckedAt(id);
	if (mcTr->GetMotherId() == -1) continue;
	mcTr->GetMomentum(mom);
	if (TMath::Abs(mom.Eta()) < 1.3) {
	  pdg = mcTr->GetPdgCode();
	  moth = ((FairMCTrack*) mcTracks->UncheckedAt(mcTr->GetMotherId()))->GetPdgCode();
	  //if (TMath::Abs(pdg) == 211 && moth == 3122) hsecR[0]->Fill(mom.Pt()); // pion from lambda
	  //else if (TMath::Abs(pdg) == 2212 && moth == 3122) hsecR[1]->Fill(mom.Pt()); // proton
	  idxs.insert(id);
	}
      }
    }
    // For UrQMD tracks
    Int_t nMC = mcTracks->GetEntriesFast();
    Int_t skip = 0;
    for (Int_t j = 0; j < nMC; ++j) {
      FairMCTrack* mcTr = (FairMCTrack*) mcTracks->UncheckedAt(j);
      //if (mcTr->GetPdgCode() == pdgCodeXi) { skip = 1; break; }
      if (mcTr->GetPdgCode() == pdgCodeXi) hXiFlag->Fill(0);
      if (mcTr->GetPdgCode() == pdgCodeL0) hLambFlag->Fill(0);
      if (mcTr->GetMotherId() == -1) continue;
      TVector3 pos;
      mcTr->GetStartVertex(pos);
      mcTr->GetMomentum(mom);
      if (TMath::Abs(mom.Eta()) < 1.3) {
	pdg = mcTr->GetPdgCode();
	moth = ((FairMCTrack*) mcTracks->UncheckedAt(mcTr->GetMotherId()))->GetPdgCode();
	if (moth == 3122) {
	  //if (TMath::Abs(pdg) == 211) hsecS[0]->Fill(mom.Pt()); // pion from lambda
	  //else if (TMath::Abs(pdg) == 2212) hsecS[1]->Fill(mom.Pt()); // proton
	  if (lun1 && idxs.find(j) == idxs.end()) fprintf(lun1,"%d %d %f %f\n",pdg,j,mom.Pt(),mom.Eta());
	}
      }
    }
    if (skip) continue;

    Int_t nMpdTr = 0;
    Int_t nITS = itsTracks->GetEntriesFast();
    Int_t nVert = vtxs->GetEntriesFast();
    if (event) mpdTracks = event->GetGlobalTracks();
    if (mpdTracks) nMpdTr = mpdTracks->GetEntriesFast();

    cout << " *** Event No: " << i << ", reco tracks in TPC (ITS), global: " << " " << nITS 
	 << " " << nMpdTr << ", vertices: " << nVert << endl;
    MpdVertex *vtx = (MpdVertex*) vtxs->First();
    mpdVert = vtx;
    vtx->Position(primVert);
    TArrayI *indxs = vtx->GetIndices();
    Int_t nPrim = indxs->GetSize();
    set<int> indxVert;
    for (Int_t k = 0; k < nPrim; ++k) indxVert.insert((*indxs)[k]);
    cout << " Number of primary (used for vertex reco) tracks: " << indxVert.size() << endl;
    
    // Find TPC track IDs 
    Int_t nPoints = 0, idMax = 0;
    /*
    if (tpcPoints) {
      nPoints = tpcPoints->GetEntriesFast();
      for (Int_t j = 0; j < nPoints; ++j) {
	FairMCPoint *p1 = (FairMCPoint*) tpcPoints->UncheckedAt(j);
	idMax = TMath::Max(idMax,p1->GetTrackID());
      }
    } else {
    */
    for (Int_t j = 0; j < nITS; ++j) {
      MpdKalmanTrack *tr = (MpdKalmanTrack*) itsTracks->UncheckedAt(j);
      idMax = TMath::Max(idMax,tr->GetTrackID());
    } 
    //}
    cout << " Max ID: " << idMax << endl;
    Int_t *ids = new Int_t [idMax+1];
    lays = new Int_t [idMax+1];
    Int_t *moths = new Int_t [idMax+1];
    Int_t *pdgs = new Int_t [idMax+1];
    Double_t *pt = new Double_t [idMax+1];
    Double_t *th = new Double_t [idMax+1];
    Double_t *rad = new Double_t [idMax+1];
    kProb = new Double_t [idMax+1];
    piProb = new Double_t [idMax+1];
    pProb = new Double_t [idMax+1];
    eProb = new Double_t [idMax+1];
    //pFalseProb = new Double_t [idMax+1];
    //piFalseProb = new Double_t [idMax+1];
    id2dst = new Int_t [idMax+1];
    Double_t *dZ = new Double_t [idMax+1];
    FairMCPoint **point = new FairMCPoint* [idMax+1];
    AzTrack **track = new AzTrack* [idMax+1];

    for (Int_t j = 0; j <= idMax; ++j) { 
      ids[j] = lays[j] = 0; 
      point[j] = 0x0;
      track[j] = 0x0;
      dZ[j] = 999999;
    }

    // Get max. reached layer No.
    for (Int_t j = 0; j < nITS; ++j) {
      AzTrack *tr = (AzTrack*) itsTracks->UncheckedAt(j);
      Int_t id = tr->GetTrackID();
      ids[id]++;
      if (ids[id] > 1) cout << " More than 1 reco track ID: " << id << " " << ids[id] << endl;
      MpdKalmanHit *hit = (MpdKalmanHit*) tr->GetTrHits()->First();
      lays[id] = TMath::Max (hit->GetLayer(), lays[id]);
    }
    //cout << " *** " << tpcPoints << endl;

    // Exclude "clones" (multiple loops)
    for (Int_t j = 0; j < nITS; ++j) {
      AzTrack *tr = (AzTrack*) itsTracks->UncheckedAt(j);
      //if (tr->GetNode() != tr->GetNodeNew()) { cout << j << " " << tr->GetNode() << " | " << tr->GetNodeNew() << endl; exit(0); }
      //if (tr->GetNode() != "") { tr->SetChi2(-9.); continue; }
      Int_t id = tr->GetTrackID();
      if (track[id] == 0x0) track[id] = tr;
      // Get ITS info
      TClonesArray *hits = tr->GetTrHits();
      Int_t nHits = hits->GetEntriesFast();
      FairMCPoint *p1 = 0x0;
      for (Int_t ih = nHits-1; ih >= 0; --ih) {
	MpdKalmanHit *hit = (MpdKalmanHit*) hits->UncheckedAt(ih);
	//if (hit->GetDist() < rMin) continue; // ITS hit
	if (hit->GetUniqueID()) continue; // ITS hit
	//if (tpcPoints) {
	if (0) {
	  p1 = (FairMCPoint*) tpcPoints->UncheckedAt(hit->GetIndex());
	  //cout << p1 << " " << hit->GetUniqueID() << " " << ids[id] << " " << point[id] << " " << p1->GetTrackID() << endl;
	  if (p1->GetTrackID() != id) continue;
	  if (ids[id] > 1 && point[id]) {
	    // More than 1 reco track with the same ID
	    //cout << " Time: " << id << " " << p1 << " " << point[id] << " " << p1->GetTime() << " " << point[id]->GetTime() << endl;
	    if (p1 == point[id]) {
	      // The same 1st hit - take "better" track
	      if (tr->GetNofTrHits() - tr->GetNofWrong() < track[id]->GetNofTrHits() - track[id]->GetNofWrong()) {
		tr->SetChi2(-9.); // exclude this track from further consideration
		break;
	      } else {
		// Exclude previous track from further consideration
		track[id]->SetChi2(-9.);
		track[id] = tr;
	      }
	    } else if (p1->GetTime() > point[id]->GetTime()) {
	      tr->SetChi2(-9.); // exclude this track from further consideration
	      break;
	    } else {
	      // Exclude previous track from further consideration
	      track[id]->SetChi2(-9.);
	      track[id] = tr;
	    }
	  } 
	  point[id] = p1;
	} else {
	  // No MC points
          if (ids[id] > 1 && point[id]) {
            // More than 1 reco track with the same ID - take the one
	    // closer to z = 0
	    //cout << " z: " << id << " " << tr->GetParam(1) << " " << track[id]->GetParam(1) << " " << tr->Charge() << " " << track[id]->Charge() << endl;
	    if (TMath::Abs(tr->GetParam(1)) < TMath::Abs(track[id]->GetParam(1))) {
	      // Exclude previous track from further consideration
 	      track[id]->SetChi2(-9.);
	      track[id] = tr;
	    } else {
	      tr->SetChi2(-9.); // exclude this track from further consideration
              break;
	    }
	  }
	  point[id] = (FairMCPoint*)0x1;
	}
	break;
      } // for (Int_t ih = nHits-1; ih >= 0;

      // MC track
      FairMCTrack* mcTr = (FairMCTrack*) mcTracks->UncheckedAt(id);
      mcTr->GetMomentum(mom);
      pt[id] = mom.Pt();
      th[id] = mom.Theta();
    } // for (Int_t j = 0; j < nITS;

    //*
    // Loop over DST tracks
    for (Int_t j = 0; j < nMpdTr; ++j) {
      MpdTrack *mpdTr = (MpdTrack*) mpdTracks->UncheckedAt(j);
      Int_t id = mpdTr->GetID();
      if (id > idMax || track[id] == 0x0) continue;
      if (ids[id] == 1) {
	kProb[id] = mpdTr->GetPidProbKaon();
	piProb[id] = mpdTr->GetPidProbPion();
	pProb[id] = mpdTr->GetPidProbProton();
	eProb[id] = mpdTr->GetPidProbElectron();
	id2dst[id] = j;
      } else {
	if (TMath::Abs(mpdTr->GetFirstPointZ()-track[id]->GetParam(1)) < dZ[id]) {
	  dZ[id] = TMath::Abs(mpdTr->GetFirstPointZ()-track[id]->GetParam(1));
	  kProb[id] = mpdTr->GetPidProbKaon();
	  piProb[id] = mpdTr->GetPidProbPion();
	  pProb[id] = mpdTr->GetPidProbProton();
	  eProb[id] = mpdTr->GetPidProbElectron();
	  id2dst[id] = j;
	}
      }
    }
    //*/

    // Lambda acceptance
    
    multimap<Int_t,Int_t> mapLamb;
    for (Int_t j = 0; j <= idMax; ++j) {
      FairMCTrack* mcTr = (FairMCTrack*) mcTracks->UncheckedAt(j);
      mcTr->GetMomentum(mom);
      Int_t mothID = mcTr->GetMotherId();
      if (mothID < 0 && lays[j] != 0) {
	lays[j] = -lays[j]; // flag primary tracks
	//href->Fill(mom.Eta());
      }
      TVector3 pos;
      mcTr->GetStartVertex(pos);
      rad[j] = pos.Pt();
      moths[j] = 0;
      pdgs[j] = mcTr->GetPdgCode();
      if (mothID >= 0) {
	moths[j] = ((FairMCTrack*) mcTracks->UncheckedAt(mothID))->GetPdgCode();
	if (moths[j] == pdgCodeL0 && (pdgs[j] == pdgCodePos || pdgs[j] == pdgCodeNeg)) mapLamb.insert(pair<Int_t,Int_t>(mothID,j));
      }
      //if ((pdgs[j] == pdgCodePos || pdgs[j] == pdgCodeNeg) && moths[j] == pdgCodeL0) 
      //hPtVsEtaS->Fill(TMath::Abs(mom.Eta()),mom.Pt());
    }

    multimap<int,int>::iterator mit, mit1;
    pair<multimap<int,int>::iterator,multimap<int,int>::iterator> ret;

    for (mit = mapLamb.begin(); mit != mapLamb.end(); ++mit) {
      Int_t mothID = mit->first;
      if (mapLamb.count(mothID) != 2) continue; // only one decay particle
      //cout << mapLamb.count(mothID) << endl;
      ret = mapLamb.equal_range(mothID);
      Int_t nppi[2] = {0}, nok = 0;
      Int_t nok1 = 0, nok2 = 0, nok3 = 0;
      for (mit1 = ret.first; mit1 != ret.second; ++mit1) {
	FairMCTrack* mcTr = (FairMCTrack*) mcTracks->UncheckedAt(mit1->second);
	if (mcTr->GetPdgCode() == pdgCodePos) nppi[0] = 1; 
	else if (mcTr->GetPdgCode() == pdgCodeNeg) nppi[1] = 1;
	//cout << mcTr->GetPdgCode() << endl;
	mcTr->GetMomentum(mom);
	if (TMath::Abs(mom.Eta()) < 1.3) ++nok;
	if ((TMath::Abs(mom.Eta())< 1.3) && mom.Pt()> 0.05) ++nok1;
	if ((TMath::Abs(mom.Eta())< 1.3) && mom.Pt()> 0.1) ++nok2;
	if ((TMath::Abs(mom.Eta())< 1.3) && mom.Pt()> 0.2) ++nok3;
      }
      if (nppi[0] != 1 || nppi[1] != 1) { cout << " Wrong decay mode !!! " << endl; continue; }   

      if (nppi[0] == 1 && nppi[1] == 1) hLambFlag->Fill(1,0.5);
      if (nok == 2) hLambFlag->Fill(2,0.5); 
      if (nok1 == 2) hLambFlag->Fill(4,0.5); 
      if (nok2 == 2) hLambFlag->Fill(6,0.5); 
      if (nok3 == 2) hLambFlag->Fill(8,0.5); 
    }
  
    // Track selection 
    for (Int_t j = 0; j < nITS; ++j) {
      AzTrack *tr = (AzTrack*) itsTracks->UncheckedAt(j);
      if (tr->GetChi2() < -8) continue;
      //if (tr->ClassName().Contains("Its") && tr->GetNofIts() > 0) continue;
      Int_t id = tr->GetTrackID();
      Double_t thRec = tr->Theta();
      Double_t etaRec = tr->Momentum3().Eta();
      //if (TMath::Abs(lays[id]) < 41 || TMath::Abs(etaRec) > 1.3) tr->SetChi2(-9.); // flag
      if (TMath::Abs(lays[id]) < -41 || TMath::Abs(etaRec) > 1.3) tr->SetChi2(-9.); // flag
      Int_t iQ = tr->Charge();
      if (iQ > 0) {
	//if (pdgs[id] != pdgCodePos || tr->GetChi2Vertex() < gC2p) tr->SetChi2(-9.);
	if (tr->GetChi2Vertex() < gC2p) tr->SetChi2(-9.);
      } else {
	//if (pdgs[id] != pdgCodeNeg || tr->GetChi2Vertex() < gC2pi) tr->SetChi2(-9.); 
	if (tr->GetChi2Vertex() < gC2pi) tr->SetChi2(-9.); 
      }
#ifdef ITS
      if (TString(tr->ClassName()).Contains("Its") && tr->GetNofHits() - tr->GetNofIts() < 10) tr->SetChi2(-9.);
#else
      if (tr->GetNofHits() < 10) tr->SetChi2(-9.);
#endif
      if (tr->GetChi2() < -8) continue;
      // Create MpdHelix
      MpdHelix helix = MakeHelix(tr);
      /*
      Double_t r = tr->GetPosNew();
      Double_t phi = tr->GetParam(0) / r;
      Double_t x = r * TMath::Cos(phi);
      Double_t y = r * TMath::Sin(phi);
      Double_t dip = tr->GetParam(3);
      Double_t cur = 0.3 * 0.01 * 5 / 10; // 5 kG
      cur *= TMath::Abs (tr->GetParam(4));
      TVector3 o(x, y, tr->GetParam(1));
      Int_t h = (Int_t) TMath::Sign(1.1,tr->GetParam(4));
      MpdHelix helix(cur, dip, tr->GetParam(2)-TMath::PiOver2()*h, o, h);
      */
      // Get 3-D DCA to primary vertex
      TVector3 pca;
      Double_t s = helix.pathLength(primVert);
      pca = helix.at(s);
      pca -= primVert;
      if (iQ < 0) {
	if (pdgs[id] != pdgCodeKm && pca.Mag() < gDCApi) tr->SetChi2(-9.);
      }
      else if (iQ > 0 && pca.Mag() < gDCAp) tr->SetChi2(-9.);
    }

    // Collect "good" pions, kaons and protons
    vector<Int_t> vecPi, vecK, vecP;
    for (Int_t j = 0; j < nITS; ++j) {
      MpdKalmanTrack *tr = (MpdKalmanTrack*) itsTracks->UncheckedAt(j);
      if (tr->GetChi2() < -8) continue;
      //if (TString(tr->ClassName()).Contains("Its") && tr->GetNofIts() > 0) continue;
      Int_t id = tr->GetTrackID();
      FairMCTrack* mcTr = (FairMCTrack*) mcTracks->UncheckedAt(id);
      if (mcTr->GetPdgCode() == pdgCodePos && tr->Charge() == 1) vecP.push_back(j);
      else if (mcTr->GetPdgCode() == pdgCodeNeg && tr->Charge() == -1) vecPi.push_back(j);
      else if (mcTr->GetPdgCode() == pdgCodeKm && tr->Charge() == -1) vecK.push_back(j);

      /*
      if (mcTr->GetPdgCode() == pdgCodePos && tr->Charge() == 1) hRecognitionFlag->Fill(1);
	else if (mcTr->GetPdgCode() == pdgCodeNeg && tr->Charge() == -1) hRecognitionFlag->Fill(5);

      if (tr->Charge() == 1 && pProb[id] > piProb[id] && pProb[id] > 0.25) {
	vecP.push_back(j);
	
	// Fill if Proton
        if (mcTr->GetPdgCode() == pdgCodePos){
	hRecognitionFlag->Fill(2); //true proton
	hProbTrueP->Fill(pProb[id],piProb[id]);
	}
	// Fill if not Proton
        if (mcTr->GetPdgCode() != pdgCodePos) hRecognitionFlag->Fill(3); //false proton
	hProbP->Fill(pProb[id],piProb[id]);
      }	
      else if (tr->Charge() == -1 && piProb[id] > pProb[id] && piProb[id] > kProb[id] && piProb[id] > eProb[id] 	&& piProb[id] > 0.25) {
	vecPi.push_back(j);
        hProbTrueK->Fill(kProb[id],piProb[id]);
	// Fill if Pion
	if (mcTr->GetPdgCode() == pdgCodeNeg){
	  hRecognitionFlag->Fill(6); // true pion
	  hProbTruePi->Fill(pProb[id],piProb[id]);
	 // if (piProb[id] > eProb[id]) hRecognitionFlag->Fill(8);
	}
	// Fill if not Pion
	if (mcTr->GetPdgCode() != pdgCodeNeg) {
	  hRecognitionFlag->Fill(7); // false pion
	  //if (piProb[id] > eProb[id]) hRecognitionFlag->Fill(9);
	  hPdg->Fill(mcTr->GetPdgCode());
	  hProbPi->Fill(pProb[id],piProb[id]);
	}
      }	
      */
    }


    vector<MpdParticle*> vecL;
    vector<Double_t> vecL1;
    vecL.clear();
    BuildLambda(vecP, vecPi, vecL, vecL1);
    if (vecL.size()) BuildCascade(vecK, vecPi, vecL, vecL1);

    Int_t nLamb = vecL.size();
    for (Int_t ipart = 0; ipart < nLamb; ++ipart) delete vecL[ipart];

    delete [] lays;
    delete [] ids;
    delete [] moths;
    delete [] pdgs;
    delete [] pt;
    delete [] th;
    delete [] point;
    delete [] rad;
    delete [] kProb;
    delete [] piProb;
    delete [] pProb;
    delete [] eProb;
    delete [] dZ;
    delete [] id2dst;
    delete [] track;
  } // for (Int_t i = 0; i < events;
  
  //TFile out("tpc.histo.root","recreate");
  /*
  TH1D *hEff1 = (TH1D*) hrec->Clone("hEff");
  hEff1->Sumw2();
  hEff1->Divide(href);
  hLayEff[0]->Divide(hLayEff[1]);
  hLayEff[2]->Divide(hLayEff[3]);
  TH2D *hEffV0 = (TH2D*) hPtVsEtaR->Clone("hEffV0");
  hEffV0->Sumw2();
  hEffV0->Divide(hPtVsEtaS);
  */

  out.Write();
  out.Close();
  if (lun) fclose(lun);
  if (lun1) fclose(lun1);
}

//__________________________________________________________________________
MpdHelix MakeHelix(const MpdKalmanTrack *tr) 
{
  Double_t r = tr->GetPosNew();
  Double_t phi = tr->GetParam(0) / r;
  Double_t x = r * TMath::Cos(phi);
  Double_t y = r * TMath::Sin(phi);
  Double_t dip = tr->GetParam(3);
  Double_t cur = 0.3 * 0.01 * 5 / 10; // 5 kG
  cur *= TMath::Abs (tr->GetParam(4));
  TVector3 o(x, y, tr->GetParam(1));
  Int_t h = (Int_t) TMath::Sign(1.1,tr->GetParam(4));
  MpdHelix helix(cur, dip, tr->GetParam(2)-TMath::PiOver2()*h, o, h);
  return helix;
}

//__________________________________________________________________________
MpdHelix MakeHelix(const MpdParticle *part) 
{
  Double_t dip = TMath::PiOver2() - part->Theta();
  Double_t cur = TMath::Abs (part->GetMeas(4));
  if (part->GetCharge() == 0) cur = numeric_limits<double>::epsilon();
  Int_t h = (Int_t) TMath::Sign(1.1,part->GetMeas(4));
  Double_t phase = part->GetMeas(2) - TMath::PiOver2() * h;
  Double_t x = part->GetXY(0);
  Double_t y = part->GetXY(1);
  TVector3 o(x, y, part->GetMeas(1));
  MpdHelix helix(cur, dip, phase, o, h);
  return helix;
}  

//__________________________________________________________________________
void BuildLambda(vector<Int_t> &vecP, vector<Int_t> &vecPi, vector<MpdParticle*> &vecL, 
		 vector<Double_t> &vecL1)
{
  // Make lambdas

  Int_t nPi = vecPi.size(), nP = vecP.size();
  vector<MpdParticle*> vPart;
  vecL2.clear();

  for (Int_t ip = 0; ip < nP; ++ip) {
    // Proton
    MpdKalmanTrack *trP = (MpdKalmanTrack*) itsTracks->UncheckedAt(vecP[ip]);
    FairMCTrack *mcTr = (FairMCTrack*) mcTracks->UncheckedAt(trP->GetTrackID());
    Int_t mothId = mcTr->GetMotherId();
    MpdParticle prot(*trP, vecP[ip]);
    prot.SetPdg(pdgCodePos);
    prot.SetMass();
    qs[1] = 1;
    etas[1] = trP->Momentum3().Eta();
    pts[1] = trP->Pt();
    //chi2s[1] = TMath::Min (trP->GetChi2Vertex(),999.);
    chi2s[1] = TMath::Min (prot.Chi2Vertex(mpdVert),999.);
    layMx[1] = TMath::Abs (lays[trP->GetTrackID()]);
    MpdHelix helix = MakeHelix(trP);
    // Get 3-D DCA to primary vertex
    TVector3 pca;
    Double_t s = helix.pathLength(primVert);
    pca = helix.at(s);
    pca -= primVert;
    dcas[1] = pca.Mag();
    origs[1] = 0;
    if (mothId >= 0 && ((FairMCTrack*) mcTracks->UncheckedAt(mothId))->GetPdgCode() == pdgCodeL0)
      origs[1] = -1; // from lambda

    for (Int_t jpi = 0; jpi < nPi; ++jpi) {
      // Pion
      MpdKalmanTrack *trPi = (MpdKalmanTrack*) itsTracks->UncheckedAt(vecPi[jpi]);
      FairMCTrack *mcTr1 = (FairMCTrack*) mcTracks->UncheckedAt(trPi->GetTrackID());
      Int_t mothId1 = mcTr1->GetMotherId();
      origs[0] = 0;
      if (mothId1 >= 0 && ((FairMCTrack*) mcTracks->UncheckedAt(mothId1))->GetPdgCode() == pdgCodeL0)
	origs[0] = -1; // from lambda
      MpdParticle *pion = new MpdParticle(*trPi, vecPi[jpi]);
      pion->SetPdg(pdgCodeNeg);
      pion->SetMass();

      vPart.clear();
      vPart.push_back(new MpdParticle(prot));
      //vPart.push_back(&prot);
      vPart.push_back(pion);

      MpdParticle lambPart;
      Double_t chi2 = TMath::Abs (lambPart.BuildMother(vPart));
      TVector3 v0(lambPart.Getx()(0,0), lambPart.Getx()(1,0), lambPart.Getx()(2,0));
      v0 -= primVert;
      Double_t decay = v0.Mag();

      if (chi2 < gC2L && decay > gPathL) {
	if (origs[1] == 1) origs[1] = -1;
	FairMCTrack *moth = NULL;
	((TH1D*)gROOT->FindObjectAny("hMassL"))->Fill(lambPart.GetMass());
	if (mothId != mothId1 || mothId < 0) {
	  ((TH1D*)gROOT->FindObjectAny("hMassLbkg"))->Fill(lambPart.GetMass());
	} else {
	  //if (moth->GetPdgCode() == pdgCodeL0) {
	  if (origs[0] == -1) {
	    ((TH1D*)gROOT->FindObjectAny("hMassLsig"))->Fill(lambPart.GetMass());
	    origs[0] = origs[1] = 1;
	    moth = (FairMCTrack*) mcTracks->UncheckedAt(mothId);
	  }
	  else ((TH1D*)gROOT->FindObjectAny("hMassLbkg"))->Fill(lambPart.GetMass());
	}

	// Fill tree
	qs[0] = -1; // pion
	etas[0] = trPi->Momentum3().Eta();
	pts[0] = trPi->Pt();
	//chi2s[0] = TMath::Min (trPi->GetChi2Vertex(),999.);
	chi2s[0] = TMath::Min (pion->Chi2Vertex(mpdVert),999.);
	layMx[0] = TMath::Abs (lays[trPi->GetTrackID()]);
	MpdHelix helix1 = MakeHelix(trPi);
	// Get 3-D DCA to primary vertex
	s = helix1.pathLength(primVert);
	pca = helix1.at(s);
	pca -= primVert;
	dcas[0] = pca.Mag();

	massh = lambPart.GetMass();
	chi2h = chi2;
	path = decay;
	pth = lambPart.Pt(); // reconstructed
        if (pth > 0.001) etah = lambPart.Momentum3().Eta(); 
        else etah = TMath::Sign(100.,lambPart.Momentum3().Z()); 
	pair<Double_t,Double_t> paths = helix.pathLengths(helix1);
	TVector3 p1 = helix.at(paths.first);
	TVector3 p2 = helix1.at(paths.second);
	p1 -= p2;
	disth = p1.Mag(); // closest distance between daughters

	// Mass cut
	//if (massLamb < 1.10713 || massLamb > 1.12423) return; // lambda mass +- 5*1.71 MeV
	//if (massLamb < 1.10695 || massLamb > 1.12525) return; // lambda mass +- 5*1.83 MeV after the selection cut sigma=1.83
        if (lambPart.GetMass() >= 1.10695 && lambPart.GetMass() <= 1.12525) {
	  vecL.push_back(new MpdParticle(lambPart));
	  vecL1.push_back(disth);
	  vecL2.push_back(pair<Double_t,Double_t>(chi2s[0],chi2s[1]));
	  //cout << chi2 << " " << lambPart.GetMass() << " " << vecL.size() << " " << trPi->Charge() << " " << trP->Charge() << " " << trPi->GetTrackID() << " " << trP->GetTrackID() << endl;
	}
	
	if (origs[0] == 1) {
	  // True lambda
	  pth = moth->GetPt();
	  yh = moth->GetRapidity();
	}
	((TTree*)gROOT->FindObjectAny("hypers"))->Fill();
      }

      //delete vPart.back(); // delete pion
      Int_t nPart = vPart.size();
      for (Int_t ipart = 0; ipart < nPart; ++ipart) delete vPart[ipart];
    } // for (Int_t jpi = 0; jpi < nPi;
    //delete vPart.front(); // delete proton
  } // for (Int_t ip = 0; ip < nP;
}

//__________________________________________________________________________
void BuildCascade(vector<Int_t> &vecK, vector<Int_t> &vecPi, vector<MpdParticle*> &vecL,
		  vector<Double_t> &vecL1)
{
  // Make cascades (Xi- and Omega-)

  Int_t nPi = vecPi.size(), nK = vecK.size(), nL = vecL.size();
  vector<MpdParticle*> vPart;

  //cout << " Cascade: " << nL << endl;
  for (Int_t iL = 0; iL < nL; ++iL) {
    // Lambda
    MpdParticle *lamb = vecL[iL];
    //cout << " Lambda No: " << iL << " " << nPi << " " << nK << " " << lamb->GetMass() << " " 
    // << lamb->GetCharge() << endl;
    // Check whether lambda is a signal or background
    Double_t lambTrue = kFALSE;
    FairMCTrack *moth0 = NULL;
    for (Int_t j = 0; j < 2; ++j) {
      MpdKalmanTrack *tr = (MpdKalmanTrack*) itsTracks->UncheckedAt(lamb->DaughterInds()[j]);
      FairMCTrack *mcTr = (FairMCTrack*) mcTracks->UncheckedAt(tr->GetTrackID());
      if (mcTr->GetMotherId() < 0) break;
      FairMCTrack *moth = (FairMCTrack*) mcTracks->UncheckedAt(mcTr->GetMotherId());
      if (moth->GetPdgCode() != pdgCodeL0) break;
      if (moth == moth0) lambTrue = kTRUE;
      moth0 = moth;
    }
    Int_t gMothId = -1;
    if (lambTrue) gMothId = moth0->GetMotherId();
    // Put lambda info here
    qs[1] = 0;
    etas[1] = lamb->Momentum3().Eta();
    pts[1] = lamb->Pt();
    Double_t xyzv[3];
    primVert.GetXYZ(xyzv);
    chi2s[1] = TMath::Min (lamb->Chi2Vertex(mpdVert),999.);
    if (lambTrue) {
      if (gMothId >= 0 && ((FairMCTrack*) mcTracks->UncheckedAt(gMothId))->GetPdgCode() == pdgCodeXi) 
	origs[1] = -1;
      else origs[1] = 0;
    } else origs[1] = -2;
    masshL = lamb->GetMass();
    disthL = vecL1[iL];
    chi2hL = lamb->Chi2();
    TVector3 v00(lamb->Getx()(0,0), lamb->Getx()(1,0), lamb->Getx()(2,0));
    v00 -= primVert;
    pathL = v00.Mag();
    chi2sL[0] = vecL2[iL].first;
    chi2sL[1] = vecL2[iL].second;
    // Get 3-D DCA to primary vertex
    MpdHelix helix = MakeHelix(lamb);
    Double_t s = helix.pathLength(primVert);
    TVector3 pca = helix.at(s);
    pca -= primVert;
    dcas[1] = pca.Mag();

    for (Int_t ik = 0; ik < nK; ++ik) {
      // Kaon
      MpdKalmanTrack *trK = (MpdKalmanTrack*) itsTracks->UncheckedAt(vecK[ik]);
      FairMCTrack *mcTr = (FairMCTrack*) mcTracks->UncheckedAt(trK->GetTrackID());
      MpdParticle *kaon = new MpdParticle(*trK, vecK[ik]);
      kaon->SetPdg(pdgCodeKm);
      kaon->SetMass();

      vPart.clear();
      vPart.push_back(new MpdParticle(*lamb));
      //vPart.push_back(lamb);
      vPart.push_back(kaon);

      MpdParticle omegaPart;
      Double_t chi2 = TMath::Abs (omegaPart.BuildMother(vPart));
      TVector3 v0(omegaPart.Getx()(0,0), omegaPart.Getx()(1,0), omegaPart.Getx()(2,0));
      v0 -= primVert;
      Double_t decay = v0.Mag();

      if (chi2 < gC2Xi && decay > gDecayOm) {
	((TH1D*)gROOT->FindObjectAny("hMassOm"))->Fill(omegaPart.GetMass());
	FairMCTrack *moth = (FairMCTrack*) mcTracks->UncheckedAt(mcTr->GetMotherId());
	if (lambTrue && gMothId >= 0 && mcTr->GetMotherId() == gMothId && moth->GetPdgCode() == pdgCodeOm) {
	  ((TH1D*)gROOT->FindObjectAny("hMassOmSig"))->Fill(omegaPart.GetMass());
	} else {
	  ((TH1D*)gROOT->FindObjectAny("hMassOmBkg"))->Fill(omegaPart.GetMass());
	}

      }

      delete vPart.front(); // delete lambda
      delete vPart.back(); // delete kaon
    } // for (Int_t ik = 0; ik < nK;

    for (Int_t jpi = 0; jpi < nPi; ++jpi) {
      // Pion
      if (vecPi[jpi] == lamb->DaughterInds()[1]) continue; // the same pion as was used for lambda

      MpdKalmanTrack *trPi = (MpdKalmanTrack*) itsTracks->UncheckedAt(vecPi[jpi]);
      MpdParticle *pion = new MpdParticle(*trPi);
      FairMCTrack *mcTr = (FairMCTrack*) mcTracks->UncheckedAt(trPi->GetTrackID());
      pion->SetPdg(pdgCodeNeg);
      pion->SetMass();

      vPart.clear();
      vPart.push_back(new MpdParticle(*lamb));
      //vPart.push_back(lamb);

      // !!! True lambda mass
      vPart.back()->SetPdg(pdgCodeL0);
      vPart.back()->SetMass();

      vPart.push_back(pion);

      MpdParticle xiPart;
      Double_t chi2 = TMath::Abs (xiPart.BuildMother(vPart));
      TVector3 v0(xiPart.Getx()(0,0), xiPart.Getx()(1,0), xiPart.Getx()(2,0));
      v0 -= primVert;
      Double_t decay = v0.Mag();

      if (chi2 < gC2Xi && decay > gDecayOm) {
	((TH1D*)gROOT->FindObjectAny("hMassXi"))->Fill(xiPart.GetMass());
	FairMCTrack *moth = NULL;
	Int_t mID = mcTr->GetMotherId();
	if (mID >= 0) moth = (FairMCTrack*) mcTracks->UncheckedAt(mID);
	if (lambTrue && gMothId >= 0 && mID == gMothId && moth->GetPdgCode() == pdgCodeXi) {
	  ((TH1D*)gROOT->FindObjectAny("hMassXiSig"))->Fill(xiPart.GetMass());
	  origs[0] = origs[1] = 1;
	} else {
	  ((TH1D*)gROOT->FindObjectAny("hMassXiBkg"))->Fill(xiPart.GetMass());
	  if (moth && moth->GetPdgCode() == pdgCodeXi) origs[0] = -1;
	  else origs[0] = 0;
	  if (origs[1] == 1) origs[1] = -1;
	}
	// Fill tree
	c2pv = TMath::Min (xiPart.Chi2Vertex(mpdVert),999.);
	qs[0] = -1; // pion
	etas[0] = trPi->Momentum3().Eta();
	pts[0] = trPi->Pt();
	//chi2s[0] = TMath::Min (trPi->GetChi2Vertex(),999.);
	chi2s[0] = TMath::Min (pion->Chi2Vertex(mpdVert),999.);
	layMx[0] = TMath::Abs (lays[trPi->GetTrackID()]);
      	// Get 3-D DCA to primary vertex
        MpdHelix helix1 = MakeHelix(trPi);
	Double_t s1 = helix1.pathLength(primVert);
	TVector3 pca1 = helix1.at(s1);
	pca1 -= primVert;
	dcas[0] = pca1.Mag();
       
	massh = xiPart.GetMass();
	chi2h = chi2;
	path = decay;
	pth = xiPart.Pt(); // reconstructed
        if (pth > 0.001) etah = xiPart.Momentum3().Eta(); 
        else etah = TMath::Sign(100.,xiPart.Momentum3().Z()); 

	disth = DistHelLin(trPi, lamb);
      	// Get 3-D DCA to primary vertex
        helix1 = MakeHelix(&xiPart);
	s1 = helix1.pathLength(primVert);
	pca1 = helix1.at(s1);
	pca1 -= primVert;
	d2pv = pca1.Mag();

	if (origs[0] == 1) {
	  // True Xi
	  pth = moth->GetPt();
	  yh = moth->GetRapidity();
	}
	((TTree*)gROOT->FindObjectAny("Xi"))->Fill();
      }

      delete vPart.front(); // delete lambda
      delete vPart.back(); // delete pion
    } // for (Int_t jpi = 0; jpi < nPi;

  } // for (Int_t iL = 0; iL < nL; 
}

//__________________________________________________________________________
Double_t DistHelLin(MpdKalmanTrack *helix, MpdParticle *neu)
{
  // Compute distance between helix and straight line

  // Create MpdHelix
  trC = MakeHelix(helix);
  vtxN.SetXYZ(neu->Getx()(0,0), neu->Getx()(1,0), neu->Getx()(2,0));
  momN = neu->Momentum3();
  momN *= (1. / momN.Mag());

  gMinuit->SetFCN(fcn);

  Double_t arglist[10];
  Int_t ierflg = 0;
  arglist[0] = 1;
  gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);
  arglist[0] = -1; //1; //-1;
  gMinuit->mnexcm("SET PRINT", arglist, 1, ierflg);

  Double_t vstart[2] = {-0.1,-0.1};
  static Double_t step[2] = {0.1, 0.1};
  gMinuit->mnparm(0, "lengN", vstart[0], step[0], 0,0,ierflg);
  gMinuit->mnparm(1, "lengC", vstart[1], step[1], 0,0,ierflg);
    
  // Now ready for minimization step
  arglist[0] = 500;
  arglist[1] = 1.;
  gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);
  
  // Get results
  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  return amin;
}

//__________________________________________________________________________
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  // Compute distance between straight line and helix

  TVector3 mom = momN;
  mom *= par[0];
  TVector3 posN = vtxN;
  posN += mom;
  
  TVector3 posC = trC.at(par[1]);
  posC -= posN;
  f = posC.Mag();
  //cout << par[0] << " " << par[1] << " " << f << endl;
}

//__________________________________________________________________________

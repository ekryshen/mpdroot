//_________________________________________________________________________________
//
//                     Macro for Lambda reconstruction @ NICA/MPD
//                                Version 16.02.2018
//                      With PID in TPC & TOF and refit after PID
//                              With MpdParticle formalism
//                        For PV "diamond" (PV smearing along Z)
//                   Lambda acceptance: +-50 cm from reconstructed PV
//                               Decay L0 -> p + pi- 
//_________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
// MPD includes
#include "TpcPoint.h"
#include "MpdEvent.h"
#include "MpdTrack.h"
#include "MpdHelix.h"
#include "MpdItsKalmanTrack.h"
#include "MpdTpcKalmanFilter.h"
#include "MpdTpcKalmanTrack.h"
#include "MpdKalmanHit.h"
#include "MpdKalmanFilter.h"
#include "MpdMotherFitterPart.h"
#include "MpdKfV0Fitter.h"
#include "MpdParticle.h"
#include "MpdPid.h"
#include "MpdVertex.h"

// CBM includes
#include "FairMCPoint.h"
#include "FairMCTrack.h"
#include "FairRunAna.h"

// ROOT includes
#include <TBranch.h>
#include <TChain.h>
#include <TClonesArray.h>
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

// PDG codes of particles 
Int_t pdgCodePr= 2212;      // Proton
Int_t pdgCodeAPr= -2212;    // antiProton
Int_t pdgCodeNeg= -211;     // pi-
Int_t pdgCodePos= 211;      // pi+
Int_t pdgCodeL0 = 3122;     // Lambda (1.11568)
Int_t pdgCodeAL0 = -3122;   // antiLambda (1.11568)
Int_t pdgCodeXi = 3312;     // Xi- (1.3213)

Int_t *lays = 0x0;
MpdHelix trC(0,0,0,TVector3(0,0,0));
MpdVertex *mpdVert;
TVector3 vtxN, momN, primVert;
TClonesArray *itsTracks, *mcTracks, *mpdTracks;

// Initialize TMinuit with a maximum of 2 params
TMinuit *gMinuit = new TMinuit(2);
FILE *lun = 0x0; 
FILE *lun1 = 0x0;

// Cuts for Lambda selection
const Double_t gC2p      = 2.;       //chi2 of p to PV
const Double_t gC2pi     = 2.;       //chi2 of pi to PV
const Double_t gDCAp     = 0.;       //cm - DCA of p to PV
const Double_t gDCApi    = 0.;       //cm - DCA of pi to PV
const Double_t gDistL    = 9999.;    //cm - DCA between pi & p in V0
const Double_t gPathL    = 0.0;      //cm - path to Lambda decay
const Double_t gC2L      = 25.;      //chi2 between pi & p in V0

#ifdef ITS
typedef MpdItsKalmanTrack AzTrack;
#else
typedef MpdTpcKalmanTrack AzTrack;
#endif
MpdTpcKalmanFilter* recoTpc = NULL;
MpdTpcKalmanFilter* recoIts = NULL;

void BuildLambda(vector<Int_t> &vecP, vector<Int_t> &vecPi, vector<MpdParticle*> &vecL);
MpdHelix MakeHelix(const MpdKalmanTrack *tr);
MpdHelix MakeHelix(const MpdParticle *part);
Double_t DistHelLin(MpdKalmanTrack *helix, MpdParticle *neu);
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
void ApplyPid(MpdPid *pid, vector<Int_t> &vecP, vector<Int_t> &vecPi);
void RecoEff(vector<Int_t> &vecP, vector<Int_t> &vecPi); 
TChain *Chain(Int_t nFiles, TString firstFile);

Float_t massh, pth, ph, etah, yh, chi2h, disth, path;     // hyperon parameters
Float_t etas[2], ps[2], pts[2], chi2s[2], dcas[2], angle; // daughters parameters: [0]-pi, [1]-p;
Double_t *kProb, *piProb, *pProb, *eProb;                 // probability for PID
Int_t evNo, origs[2], qs[2], dstNo[2], layMx[2];
Int_t *id2dst;
vector<pair<Double_t,Double_t> > vecL1, vecL2;

//________________________________________________________________________________
//  
//             Analyze TPC (ITS) reco data - reconstruct Lambda hyperon
//                     (with track refit to account for dE/dx)
//---------------------------------------------------------------------------------

void AnalL0(Int_t n1 = 0, Int_t n2 = 0, Int_t firstFile = 1)
{   
  //gROOT->ProcessLine(".x ~/mpd/Chain.C(1,\"./mc_1.root\")");
  Chain(1,"./mc_1.root");
  TChain *simMC = (TChain*) gROOT->FindObject("cbmsim");
  simMC->SetName("cbmsim1");
  TString fileName = "./MPD-DST/dst-2018-01-18-a4bd06e-LAQGSM-1G3Mlem/AuAuss11_5mb5_mer.r12.mpddst-";
  fileName += firstFile;
  fileName += ".root";
  Chain(10,fileName);
  TChain *simITS = (TChain*) gROOT->FindObject("cbmsim");
  
  TFile fileITS(simITS->GetListOfFiles()->First()->GetTitle());
  TFile fileMC(simMC->GetListOfFiles()->First()->GetTitle());
  fileMC.Get("FairGeoParSet");
  TClonesArray *vtxs = (TClonesArray*) fileITS.FindObjectAny("Vertex");
  simITS->SetBranchAddress("Vertex",&vtxs);
  TBranch *vtxB = simITS->GetBranch("Vertex");
  itsTracks = (TClonesArray*) fileITS.FindObjectAny("ItsTrack");
  simITS->SetBranchAddress("ItsTrack",&itsTracks);
  TBranch *itsRecoB = simITS->GetBranch("ItsTrack");
#ifndef ITS
  itsTracks = NULL;
#endif
  if (itsTracks == 0x0) {
    itsTracks = (TClonesArray*) fileITS.FindObjectAny("TpcKalmanTrack");
    simITS->SetBranchAddress("TpcKalmanTrack",&itsTracks);
    itsRecoB = simITS->GetBranch("TpcKalmanTrack");
  }
  MpdEvent *event = 0x0;
  simITS->SetBranchAddress("MPDEvent.", &event);
  FairMCEventHeader *mchead = 0x0;
  simITS->SetBranchAddress("MCEventHeader.", &mchead);
  
  TClonesArray *tpcPoints = (TClonesArray*) fileMC.FindObjectAny("TpcPoint");
  simMC->SetBranchAddress("TpcPoint",&tpcPoints);
  TBranch *tpcSimB = simMC->GetBranch("TpcPoint");
  TClonesArray *itsPoints = (TClonesArray*) fileMC.FindObjectAny("StsPoint");
  simMC->SetBranchAddress("StsPoint",&itsPoints);
  TBranch *itsSimB = simMC->GetBranch("StsPoint");  
  mcTracks = (TClonesArray*) fileITS.FindObjectAny("MCTrack");
  simITS->SetBranchAddress("MCTrack",&mcTracks);
  TBranch *mcBranch = simITS->GetBranch("MCTrack");
  
  TFile out("L0_10k_PVsmear.test.root","recreate");
  
  FairRunAna ana;
  MpdKalmanFilter::Instance("KF")->Init();  
#ifdef ITS
  recoIts = new MpdTrackFinderIts5spd();
  recoIts->FillGeoScheme();
#endif
  recoTpc = new MpdTpcKalmanFilter("TPC Kalman filter");
  recoTpc->SetSectorGeo(MpdTpcSectorGeo::Instance());
  recoTpc->FillGeoScheme();
  
  // n-sigma bands for PID selection
  Double_t sigM = 4.0, sigE = 4.0, energy = 11.0, coef = 1.0; 
  TString generator = "LAQGSM", tracking = "CF";
  MpdPid *newPid = new MpdPid(sigM, sigE, energy, coef, generator, tracking, "pikarp");
  
  // Book histograms for Lambda
  TH1D *hLambFlag = new TH1D("hLambFlag","Flags for lambda",12,0,12);
  TH1D *hRecognitionFlag = new TH1D("hRecognitionFlag","Flags for Recognition",10,0,10);
  TH1D *hMassL = new TH1D("hMassL","Lambda mass",50,1.070,1.170);
  TH1D *hMassLsig = new TH1D("hMassLsig","Lambda mass (signal)",50,1.070,1.170);
  TH1D *hMassLbkg = new TH1D("hMassLbkg","Lambda mass (bkg.)",50,1.070,1.170);
  TH1D *hLambPTH = new TH1D("hLambPTH","Flags for lambdaPTH",12,0,12);
  TH1D *hPdg = new TH1D("hPdg","PdgCod if is not Pion",1000,-2500,2500);

  // Book histograms for PID check
  TH2D *hProbTrueP = new TH2D("hProbTrueP","Probability for true Protons",50,0,1,50,0,1.1);
  TH2D *hProbP = new TH2D("hProbfalseP","Probability for Pions and identification Protons",50,0,1.1,50,0,1.1);
  TH2D *hProbPi = new TH2D("hProbfalsePi","Probability for Protons and identification Pions",50,0,1.1,50,0,1.1);
  TH2D *hProbTruePi = new TH2D("hProbTruePi","Probability for true Pions",50,0,1.1,50,0,1.1);
  TH2D *hProbTrueK = new TH2D("hProbTrueK","Probability for true Kaons",50,0,1.1,50,0,1.1);
  new TH1D("hPIDflag","PID flags",12,0,12);
  
  Int_t moth, pdg;
    
  // TTree for Lambda selection
  TTree *thyps = new TTree("hypers","Hyperons");
  thyps->Branch("massh",&massh,"massh/F");        // inv. mass of hyperon
  thyps->Branch("pth",&pth,"pth/F");              // pt of hyperon
  thyps->Branch("yh",&yh,"yh/F");                 // rapidity of hyperon
  thyps->Branch("chi2h",&chi2h,"chi2h/F");        // chi2 between pi & p in V0
  thyps->Branch("disth",&disth,"disth/F");        // dca between pi & p in V0
  thyps->Branch("path",&path,"path/F");           // path of Lambda
  thyps->Branch("angle",&angle,"angle/F");        // angle between L0 momentum & PV to V0 direction
  thyps->Branch("chi2s[2]",&chi2s,"chi2s[2]/F");  // chi2 to PV of pi [0] & p [1]
  thyps->Branch("dcas[2]",&dcas,"dcas[2]/F");     // dca to PV of pi [0] & p [1]
  thyps->Branch("origs[2]",&origs,"origs[2]/I2"); // original p & pi (from Lambda)
  thyps->Branch("qs[2]",&qs,"qs[2]/I2");          // charge of pi [0] & p [1]
  thyps->Branch("layMx[2]",&layMx,"layMx[2]/I2"); // maximum padrows reached
  
  Int_t events = simITS->GetEntries();
  if (n2 != 0) events = TMath::Min (events, n2);
  cout << " Number of events = " << events << endl;
  
  //*******************************************************************************
  //*****************************   Loop over events    ***************************
  //*******************************************************************************
  
  for (Int_t i = 0; i < events; ++i) {
    if (i < n1) continue;
    simITS->GetEntry(i);
    evNo = i + 1;
    
    TVector3 genVert;
    mchead->GetVertex(genVert);

    //********************  Loop over ITS points   *********************
    
    TVector3 mom; 
    set<Int_t> idxs;
    if (itsPoints) {
      Int_t nITSp = itsPoints->GetEntriesFast();
      for (Int_t j = 0; j < nITSp; ++j) {
	FairMCPoint *p1 = (FairMCPoint*) itsPoints->UncheckedAt(j);
	Int_t id = p1->GetTrackID();
	if (idxs.find(id) != idxs.end()) continue;
	FairMCTrack* mcTr = (FairMCTrack*) mcTracks->UncheckedAt(id);
	if (mcTr->GetMotherId() == -1) continue;
	mcTr->GetMomentum(mom);
	if (TMath::Abs(mom.Eta()) < 1.3) {
	  pdg = mcTr->GetPdgCode();
	  moth = ((FairMCTrack*) mcTracks->UncheckedAt(mcTr->GetMotherId()))->GetPdgCode();	  
	  idxs.insert(id);
	}
      }
    } //***************  End of Loop over ITS points   *****************

    //                             *****
    
    //********************  Loop over MC tracks   **********************
 
    Int_t nMC = mcTracks->GetEntriesFast();
    Int_t skip = 0;
    for (Int_t j = 0; j < nMC; ++j) {
      FairMCTrack* mcTr = (FairMCTrack*) mcTracks->UncheckedAt(j);
      mcTr->GetMomentum(mom);
      if (mcTr->GetPdgCode() == pdgCodeL0) {
	// Check production vertex
	TVector3 pos;
	Double_t r = 0.0;
	if (mcTr->GetMotherId() >= 0) mcTr->GetStartVertex(pos);
	pos -= genVert;
	r = pos.Mag();
	if (r < 50.0) {
	  // Production vertex constraint 50 cm
	  hLambFlag->Fill(0);
	  hLambPTH->Fill(0);
	  if (mom.Pt() < 0.5) hLambPTH->Fill(2);
	  if (mom.Pt() > 0.5 && mom.Pt() < 1.0) hLambPTH->Fill(4);
	  if (mom.Pt() > 1.0 && mom.Pt() < 1.5) hLambPTH->Fill(6);
	  if (mom.Pt() > 1.5 && mom.Pt() < 2.0) hLambPTH->Fill(8);
	  if (mom.Pt() > 2.0) hLambPTH->Fill(10);
	}
      }
      if (mcTr->GetMotherId() < 0) continue;
      TVector3 pos;
      mcTr->GetStartVertex(pos);
      if (mom.Pt() < 0.001) continue;
      if (TMath::Abs(mom.Eta()) < 1.3) {
	pdg = mcTr->GetPdgCode();
	moth = ((FairMCTrack*) mcTracks->UncheckedAt(mcTr->GetMotherId()))->GetPdgCode();
	if (moth == 3122) {	  
	  if (lun1 && idxs.find(j) == idxs.end()) fprintf(lun1,"%d %d %f %f\n",pdg,j,mom.Pt(),mom.Eta());
	}
      }
    } //****************  End of Loop over MC tracks   ******************
    
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
    for (Int_t j = 0; j < nITS; ++j) {
      MpdKalmanTrack *tr = (MpdKalmanTrack*) itsTracks->UncheckedAt(j);
      idMax = TMath::Max(idMax,tr->GetTrackID());
    } 
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

    //*********************   Loop over tracks   ***********************
    
    //                Exclude "clones" (multiple loops)
    
    for (Int_t j = 0; j < nITS; ++j) {
      AzTrack *tr = (AzTrack*) itsTracks->UncheckedAt(j);     
      Int_t id = tr->GetTrackID();
      if (track[id] == 0x0) track[id] = tr;
      // Get ITS info
      TClonesArray *hits = tr->GetTrHits();
      Int_t nHits = hits->GetEntriesFast();
      FairMCPoint *p1 = 0x0;
      for (Int_t ih = nHits-1; ih >= 0; --ih) {
	MpdKalmanHit *hit = (MpdKalmanHit*) hits->UncheckedAt(ih);
	if (hit->GetUniqueID()) continue; // ITS hit
	if (0) {
	  p1 = (FairMCPoint*) tpcPoints->UncheckedAt(hit->GetIndex());
	  if (p1->GetTrackID() != id) continue;
	  if (ids[id] > 1 && point[id]) {
	    // More than 1 reco track with the same ID
	    if (p1 == point[id]) {
	      // The same 1st hit - take "better" track
	      if (tr->GetNofTrHits() - tr->GetNofWrong() < track[id]->GetNofTrHits() - track[id]->GetNofWrong()) {
		tr->SetChi2(-9.); // Exclude this track from further consideration
		break;
	      } else {
		track[id]->SetChi2(-9.); // Exclude previous track from further consideration
		track[id] = tr;
	      }
	    } else if (p1->GetTime() > point[id]->GetTime()) {
	      tr->SetChi2(-9.); // Exclude this track from further consideration
	      break;
	    } else {
	      track[id]->SetChi2(-9.); // Exclude previous track from further consideration
	      track[id] = tr;
	    }
	  } 
	  point[id] = p1;
	} else {
	  // No MC points
          if (ids[id] > 1 && point[id]) {            
	    if (TMath::Abs(tr->GetParam(1)) < TMath::Abs(track[id]->GetParam(1))) {
 	      track[id]->SetChi2(-9.); // Exclude previous track from further consideration
	      track[id] = tr;
	    } else {
	      tr->SetChi2(-9.); // Exclude this track from further consideration
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
    } //****************  End of Loop over tracks   *****************
    
    //                             *****
    
    //******************  Loop over DST tracks  *********************
    
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
    } //***************  End of Loop over DST tracks  *****************

    //                             *****

    //**********************  Lambda acceptance  **********************
  
    multimap<Int_t,Int_t> mapLamb;
    for (Int_t j = 0; j <= idMax; ++j) {
      FairMCTrack* mcTr = (FairMCTrack*) mcTracks->UncheckedAt(j);
      mcTr->GetMomentum(mom);
      Int_t mothID = mcTr->GetMotherId();
      if (mothID == -1 && lays[j] != 0) {
	lays[j] = -lays[j]; // flag primary tracks
      }
      TVector3 pos;
      mcTr->GetStartVertex(pos);
      rad[j] = pos.Pt();
      moths[j] = 0;
      pdgs[j] = mcTr->GetPdgCode();
      if (mothID >= 0) {
	// Check Lambda production vertex ( < 50 cm)
	FairMCTrack* moth = (FairMCTrack*) mcTracks->UncheckedAt(mothID);
	moth->GetStartVertex(pos);
	pos -= genVert;
	if (pos.Mag() < 50.0) {
	  moths[j] = moth->GetPdgCode();
	  if (moths[j] == pdgCodeL0 && (pdgs[j] == pdgCodePr || pdgs[j] == pdgCodeNeg)) mapLamb.insert(pair<Int_t,Int_t>(mothID,j));
	}
      }      
    } //******************  End of Lambda acceptance  *******************
    
    //                             *****

    //*************************   Pt selection  *************************
    
    multimap<int,int>::iterator mit, mit1;
    pair<multimap<int,int>::iterator,multimap<int,int>::iterator> ret;
    
    for (mit = mapLamb.begin(); mit != mapLamb.end(); ++mit) {
      Int_t mothID = mit->first;
      if (mapLamb.count(mothID) != 2) continue; // only one decay particle
      ret = mapLamb.equal_range(mothID);
      Int_t nppi[2] = {0}, nok = 0;
      Int_t nok1 = 0, nok2 = 0, nok3 = 0;
      for (mit1 = ret.first; mit1 != ret.second; ++mit1) {
	FairMCTrack* mcTr = (FairMCTrack*) mcTracks->UncheckedAt(mit1->second);
	if (mcTr->GetPdgCode() == pdgCodePr) nppi[0] = 1; 
	else if (mcTr->GetPdgCode() == pdgCodeNeg) nppi[1] = 1;
	mcTr->GetMomentum(mom);
	if (mom.Pt() < 0.001) continue;
	if (TMath::Abs(mom.Eta()) < 1.3) ++nok;
	if ((TMath::Abs(mom.Eta())< 1.3) && mom.Pt()> 0.05) ++nok1;
	if ((TMath::Abs(mom.Eta())< 1.3) && mom.Pt()> 0.1) ++nok2;
	if ((TMath::Abs(mom.Eta())< 1.3) && mom.Pt()> 0.2) ++nok3;
      }
      if (nppi[0] != 1 || nppi[1] != 1) { cout << " Wrong decay mode !!! " << endl; continue; } // not p - pi- decay   
      
      if (nppi[0] == 1 && nppi[1] == 1) hLambFlag->Fill(1,0.5); // all generated Lambda (+-50 cm from PV)
      if (nok == 2) hLambFlag->Fill(2,0.5);                     // all reconstructed Lambda (... + |eta| < 1.3)
      if (nok1 == 2) hLambFlag->Fill(4,0.5);                    // ... + pT > 0.05
      if (nok2 == 2) hLambFlag->Fill(6,0.5);                    // ... + pT > 0.1
      if (nok3 == 2) hLambFlag->Fill(8,0.5);                    // ... + pT > 0.2
      
    } //*******************  End of Pt selection  ********************

    //                             *****

    //**********************   Track selection  **********************
    
    for (Int_t j = 0; j < nITS; ++j) {
      AzTrack *tr = (AzTrack*) itsTracks->UncheckedAt(j);
      if (tr->GetChi2() < -8) continue;
      Int_t id = tr->GetTrackID();
      Double_t thRec = tr->Theta();
      Double_t etaRec = tr->Momentum3().Eta();
      if (TMath::Abs(lays[id]) < -41 || TMath::Abs(etaRec) > 1.3) tr->SetChi2(-9.); // flag
      Int_t iQ = tr->Charge();
      
#ifdef ITS
      if (TString(tr->ClassName()).Contains("Its") && tr->GetNofHits() - tr->GetNofIts() < 10) tr->SetChi2(-9.);
#else
      if (tr->GetNofHits() < 10) tr->SetChi2(-9.);
#endif
      if (tr->GetChi2() < -8) continue;

      // Create MpdHelix
      MpdHelix helix = MakeHelix(tr);
      
      // Get 3-D DCA to primary vertex
      TVector3 pca;
      Double_t s = helix.pathLength(primVert);
      pca = helix.at(s);
      pca -= primVert;
      if (iQ < 0) {
	if (pca.Mag() < gDCApi) tr->SetChi2(-9.);
      }
      else if (iQ > 0 && pca.Mag() < gDCAp) tr->SetChi2(-9.);
      
    } //******************  End of Track selection  ******************

    //                             *****

    //*************  Collect "good" pions and protons  ***************
    
    vector<Int_t> vecPi, vecK, vecP;
    for (Int_t j = 0; j < nITS; ++j) {
      MpdKalmanTrack *tr = (MpdKalmanTrack*) itsTracks->UncheckedAt(j);
      if (tr->GetChi2() < -8) continue;
      Int_t id = tr->GetTrackID();
      FairMCTrack* mcTr = (FairMCTrack*) mcTracks->UncheckedAt(id);
      if (mcTr->GetMotherId() == 0 && 
	  ((FairMCTrack*)mcTracks->UncheckedAt(0))->GetPdgCode() == 1010010030) continue; // !!! decay product of artificial H3L
      //*MC ID
      if (mcTr->GetPdgCode() == pdgCodePr && tr->Charge() == 1) vecP.push_back(j);
      else if (mcTr->GetPdgCode() == pdgCodeNeg && tr->Charge() == -1) vecPi.push_back(j);
      if (mcTr->GetPdgCode() == pdgCodePos && tr->Charge() == 1) hRecognitionFlag->Fill(1);
      else if (mcTr->GetPdgCode() == pdgCodeNeg && tr->Charge() == -1) hRecognitionFlag->Fill(5);
      if (tr->Charge() == 1 && pProb[id] > piProb[id] && pProb[id] > 0.25) {       
       	// Fill if Proton
        if (mcTr->GetPdgCode() == pdgCodePos){
	  hRecognitionFlag->Fill(2); // true Proton
	  hProbTrueP->Fill(pProb[id],piProb[id]);
	}
	// Fill if not Proton
        if (mcTr->GetPdgCode() != pdgCodePos) hRecognitionFlag->Fill(3); // false Proton
	hProbP->Fill(pProb[id],piProb[id]);
      }	
      else if (tr->Charge() == -1 && piProb[id] > pProb[id] && piProb[id] > kProb[id] && piProb[id] > eProb[id] 	
	       && piProb[id] > 0.25) {
        hProbTrueK->Fill(kProb[id],piProb[id]);
	// Fill if Pion
	if (mcTr->GetPdgCode() == pdgCodeNeg){
	  hRecognitionFlag->Fill(6); // true Pion
	  hProbTruePi->Fill
	    (pProb[id],piProb[id]);
	}
	// Fill if not Pion
	if (mcTr->GetPdgCode() != pdgCodeNeg) {
	  hRecognitionFlag->Fill(7); // false Pion
	  hPdg->Fill(mcTr->GetPdgCode());
	  hProbPi->Fill(pProb[id],piProb[id]);
	}
      }	
    } //*******************  End of Collection  ***********************
    
    cout << " Number of protons, pi: " << vecP.size() << " " << vecPi.size() << endl;
    RecoEff(vecP, vecPi);
    cout << " Number of protons, pi: " << vecP.size() << " " << vecPi.size() << endl;
    
    //*************************  Apply PID  ***************************

    ApplyPid(newPid, vecP, vecPi);
    
    vector<MpdParticle*> vecL;
    vecL.clear();
    BuildLambda(vecP, vecPi, vecL);
    
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
  }

  //*******************************************************************************
  //*************************   End of  Loop over events    ***********************
  //*******************************************************************************
  
  out.Write();
  out.Close();
  if (lun) fclose(lun);
  if (lun1) fclose(lun1);
}  

//_________________________________________________________________________________
//
//                                Lambda reconstruction
//_________________________________________________________________________________

void BuildLambda(vector<Int_t> &vecP, vector<Int_t> &vecPi, vector<MpdParticle*> &vecL) 
{
  Int_t nPi = vecPi.size(), nP = vecP.size();
  vector<MpdParticle*> vPart;
  vecL1.clear();
  vecL2.clear();
  
  //**********************  Loop over Proton  *********************

  for (Int_t ip = 0; ip < nP; ++ip) {
    AzTrack *trP = (AzTrack*) itsTracks->UncheckedAt(vecP[ip]);
    FairMCTrack *mcTr = (FairMCTrack*) mcTracks->UncheckedAt(trP->GetTrackID());
    Int_t mothId = mcTr->GetMotherId();
    AzTrack trCor = *trP;
    trCor.SetDirection(MpdKalmanTrack::kInward);
    if (recoIts) recoIts->Refit((MpdItsKalmanTrack*)&trCor, 0.93827, 1); // refit
    else recoTpc->Refit(&trCor, 0.93827, 1); // refit
    MpdParticle prot(trCor, vecP[ip]);
    prot.SetPdg(pdgCodePr);
    prot.SetMass();
    qs[1] = 1;
    etas[1] = trP->Momentum3().Eta();
    ps[1] = trP->Momentum();
    pts[1] = trP->Pt();
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
      origs[1] = -1; // Proton from Lambda
    
    //*********************  Loop over Pion  ************************
    
    for (Int_t jpi = 0; jpi < nPi; ++jpi) {
      MpdKalmanTrack *trPi = (MpdKalmanTrack*) itsTracks->UncheckedAt(vecPi[jpi]);
      FairMCTrack *mcTr1 = (FairMCTrack*) mcTracks->UncheckedAt(trPi->GetTrackID());
      Int_t mothId1 = mcTr1->GetMotherId();
      origs[0] = 0;
      if (mothId1 >= 0 && ((FairMCTrack*) mcTracks->UncheckedAt(mothId1))->GetPdgCode() == pdgCodeL0)
	origs[0] = -1; // Pion from Lambda
      MpdParticle *pion = new MpdParticle(*trPi, vecPi[jpi]);
      pion->SetPdg(pdgCodeNeg);
      pion->SetMass();

      vPart.clear();
      vPart.push_back(new MpdParticle(prot));
      vPart.push_back(pion);

      MpdParticle lambPart;
      Double_t chi2 = TMath::Abs (lambPart.BuildMother(vPart));
      TVector3 v0(lambPart.Getx()(0,0), lambPart.Getx()(1,0), lambPart.Getx()(2,0));
      v0 -= primVert;
      Double_t decay = v0.Mag();
      path = TMath::Sign (decay, v0*lambPart.Momentum3());

      if (chi2 < gC2L && path > gPathL) {
	if (origs[1] > 0) origs[1] = -1;
	FairMCTrack *moth = NULL;
	((TH1D*)gROOT->FindObjectAny("hMassL"))->Fill(lambPart.GetMass());
	if (mothId != mothId1 || mothId < 0) {
	  ((TH1D*)gROOT->FindObjectAny("hMassLbkg"))->Fill(lambPart.GetMass());
	} else {
	  if (origs[0] == -1) {
	    ((TH1D*)gROOT->FindObjectAny("hMassLsig"))->Fill(lambPart.GetMass());
	    origs[0] = origs[1] = 1;  // True Lambda
	    moth = (FairMCTrack*) mcTracks->UncheckedAt(mothId);
	  }
	  else ((TH1D*)gROOT->FindObjectAny("hMassLbkg"))->Fill(lambPart.GetMass());
	}

	// Fill tree
	qs[0] = -1; // Charge of Pion
	etas[0] = trPi->Momentum3().Eta();
	ps[0] = trPi->Momentum();
	pts[0] = trPi->Pt();
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
	pth = lambPart.Pt(); // reconstructed
	ph = lambPart.Momentum(); // reconstructed
	yh = lambPart.Rapidity();
	angle = v0.Angle(lambPart.Momentum3());
        if (pth > 0.001) etah = lambPart.Momentum3().Eta(); 
        else etah = TMath::Sign(100.,lambPart.Momentum3().Z()); 
	pair<Double_t,Double_t> paths = helix.pathLengths(helix1);
	TVector3 p1 = helix.at(paths.first);
	TVector3 p2 = helix1.at(paths.second);
	p1 -= p2;
	disth = p1.Mag(); // Closest distance between daughters

	// Mass cut
        if (lambPart.GetMass() >= 1.10695 && lambPart.GetMass() <= 1.12525) {
	  vecL.push_back(new MpdParticle(lambPart));
	  vecL1.push_back(pair<Double_t,Double_t>(disth,angle));
	  vecL2.push_back(pair<Double_t,Double_t>(chi2s[0],chi2s[1]));	  
	}
	
	if (origs[0] == 1) {         // True Lambda
	  lambPart.SetMass(1.11568); // Set true mass
	  // Check mother of Lambda
	  Int_t gMothId = moth->GetMotherId();
	  if (gMothId >= 0) origs[0] = origs[1] = 2; // Secondary Lambda
	}
	((TTree*)gROOT->FindObjectAny("hypers"))->Fill();
      }

      Int_t nPart = vPart.size();
      for (Int_t ipart = 0; ipart < nPart; ++ipart) delete vPart[ipart];
      
    } //*******************  End of Loop over Pion  **********************
  } //*******************  End of Loop over Proton  **********************
}

//_________________________________________________________________________________
//
//                                   Make helix
//_________________________________________________________________________________

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

//_________________________________________________________________________________
//
//               Compute distance between helix and straight line
//_________________________________________________________________________________

Double_t DistHelLin(MpdKalmanTrack *helix, MpdParticle *neu)
{
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
  arglist[0] = -1; 
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

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  TVector3 mom = momN;
  mom *= par[0];
  TVector3 posN = vtxN;
  posN += mom;
  
  TVector3 posC = trC.at(par[1]);
  posC -= posN;
  f = posC.Mag();
}

//_________________________________________________________________________________
//
//                                   Apply PID
//_________________________________________________________________________________

void ApplyPid(MpdPid *pid, vector<Int_t> &vecP, vector<Int_t> &vecPi)
{
  vecP.clear();
  vecPi.clear();

  Int_t nITS = itsTracks->GetEntriesFast();
  
  for (Int_t j = 0; j < nITS; ++j) {
    AzTrack *tr = (AzTrack*) itsTracks->UncheckedAt(j);
    if (tr->GetChi2() < -8) continue;
    Int_t id = tr->GetTrackID();
    FairMCTrack* mcTr = (FairMCTrack*) mcTracks->UncheckedAt(id);
    Int_t mothId = mcTr->GetMotherId();   
    
    MpdTrack* mpdTrack = (MpdTrack*) mpdTracks->UncheckedAt(j);
    if (mpdTrack->GetID() != id) { cout << id << " " << mpdTrack->GetID() << endl; Fatal("ApplyPid"," Different ID"); }
    
    Int_t ret = 0, charge = tr->Charge(), tofFlag = mpdTrack->GetTofFlag();
    Double_t dedx = tr->GetPartID(), m2 = mpdTrack->GetTofMass2();
    
    // Matching TPC+TOF (dE/dx & m^2)
    if (tofFlag == 2 || tofFlag == 6) 
      {
	ret = pid->FillProbs(tr->Momentum(), dedx, m2, charge);
	// Mismatching, only dE/dx used for momentum < 0.8 GeV/c
	if ( ret == 0 && tr->Momentum() < 0.8 ) ret = pid->FillProbs(tr->Momentum(), dedx, charge); 
      }
    
    // No matching TPC+TOF (only dE/dx available)
    if (tofFlag == 0 || tofFlag == 4)
      {
        ret = pid->FillProbs(tr->Momentum(), dedx, charge);
      }   

    // No PID
    TH1D *hFlag = (TH1D*) gROOT->FindObjectAny("hPIDflag");
    if (ret == 0) {
      if (mcTr->GetPdgCode() == pdgCodePos) hFlag->Fill(2.1); // lost Pion
      if (mcTr->GetPdgCode() == pdgCodeAPr) hFlag->Fill(6.1); // lost Proton
      continue;
    }   
    
    Double_t piThr = -0.75, probThr = -0.60;
    
    if (tr->Charge() < 0) {
      Double_t prob = pid->GetProbPi();
      if (prob > piThr && prob > pid->GetProbKa() && prob > pid->GetProbEl() && prob > pid->GetProbPr() &&
	  prob > pid->GetProbMu()) {
	// "Pion"
	if (mcTr->GetPdgCode() == pdgCodeNeg) hFlag->Fill(0.1); // correct Pion
	else if (mcTr->GetPdgCode() != pdgCodeNeg) hFlag->Fill(1.1); // false Pion
        Double_t chi2 = TMath::Min (tr->GetChi2Vertex(),999.);
        if (chi2 < gC2pi) continue;
	vecPi.push_back(j);
      } else if (mcTr->GetPdgCode() == pdgCodeNeg) hFlag->Fill(2.1); // lost Pion
    } else {
      Double_t prob = pid->GetProbPr();
      if (prob > probThr && prob > pid->GetProbKa() && prob > pid->GetProbPi() && prob > pid->GetProbDe()) {
	// "Proton"
	if (mcTr->GetPdgCode() == pdgCodePr) hFlag->Fill(4.1); // correct Proton
	else if (mcTr->GetPdgCode() != pdgCodePr) hFlag->Fill(5.1); // false Proton
	AzTrack trCor = *tr;
	trCor.SetDirection(MpdKalmanTrack::kInward);
	if (recoIts) recoIts->Refit((MpdItsKalmanTrack*)&trCor, 0.93827, 1); // refit
	else recoTpc->Refit(&trCor, 0.93827, 1); // refit
	MpdParticle prot(trCor, 0);
	prot.SetPdg(pdgCodePr);
	prot.SetMass();
	Double_t chi2 = TMath::Min (prot.Chi2Vertex(mpdVert),999.);
	if (chi2 < gC2p) continue;
 	vecP.push_back(j);
      } else if (mcTr->GetPdgCode() == pdgCodePr) hFlag->Fill(6.1); // lost Proton
    }
  }    
  cout << " Number of p, pi: " << vecP.size() << " " << vecPi.size() << endl;
}

//_________________________________________________________________________________
//
//                               Check reco efficiency
//_________________________________________________________________________________

void RecoEff(vector<Int_t> &vecP, vector<Int_t> &vecPi)
{
  Int_t nPi = vecPi.size(), nP = vecP.size();
  
  for (Int_t ip = nP - 1; ip >= 0; --ip) {
    // AntiProton
    AzTrack *trP = (AzTrack*) itsTracks->UncheckedAt(vecP[ip]);
    FairMCTrack *mcTr = (FairMCTrack*) mcTracks->UncheckedAt(trP->GetTrackID());
    Int_t mothId = mcTr->GetMotherId();
    if (mothId < 0) continue;
    FairMCTrack *moth = (FairMCTrack*) mcTracks->UncheckedAt(mothId);  
    if (moth->GetPdgCode() == pdgCodeL0) {
      Int_t mp = mothId;
      // Proton from Lambda
      for (Int_t jpi = nPi - 1; jpi >= 0; --jpi) {
	// Pion
	AzTrack *trPi = (AzTrack*) itsTracks->UncheckedAt(vecPi[jpi]);
	FairMCTrack *mcTr = (FairMCTrack*) mcTracks->UncheckedAt(trPi->GetTrackID());
	Int_t mothId = mcTr->GetMotherId();
	if (mothId < 0) continue;
	FairMCTrack *moth = (FairMCTrack*) mcTracks->UncheckedAt(mothId);  
	if (moth->GetPdgCode() == pdgCodeL0 && mp == mothId) {
	  ((TH1D*)gROOT->FindObjectAny("hLambFlag"))->Fill(10);
	  Int_t gmId = moth->GetMotherId();
	  if (gmId >= 0) {
	    FairMCTrack *gmoth = (FairMCTrack*) mcTracks->UncheckedAt(gmId);
	    if (gmoth->GetPdgCode() == pdgCodeXi) {
	      for (Int_t kpi = nPi - 1; kpi >= 0; --kpi) {
		// Pion
		AzTrack *trK = (AzTrack*) itsTracks->UncheckedAt(vecPi[kpi]);
		FairMCTrack *mcTr = (FairMCTrack*) mcTracks->UncheckedAt(trK->GetTrackID());
		Int_t mothId = mcTr->GetMotherId();
		if (mothId < 0) continue;
		FairMCTrack *moth = (FairMCTrack*) mcTracks->UncheckedAt(mothId);  
		if (moth->GetPdgCode() == pdgCodeXi && gmId == mothId) {
		  break;
		}
	      }
	    } 
	  }
	  break;
	}
      }
    }
  }
  
  for (Int_t ip = nP - 1; ip >= 0; --ip) {
    // Proton
    AzTrack *trP = (AzTrack*) itsTracks->UncheckedAt(vecP[ip]);
    AzTrack trCor = *trP;
    trCor.SetDirection(MpdKalmanTrack::kInward);
    if (recoIts) recoIts->Refit((MpdItsKalmanTrack*)&trCor, 0.93827, 1); // refit
    else recoTpc->Refit(&trCor, 0.93827, 1); // refit
    MpdParticle prot(trCor, vecP[ip]);
    prot.SetPdg(pdgCodePr);
    prot.SetMass();
    Double_t chi2 = TMath::Min (prot.Chi2Vertex(mpdVert),999.);
    if (chi2 < gC2p) vecP.erase(vecP.begin()+ip);
  }
  
  if (nP) {
    for (Int_t jpi = nPi - 1; jpi >= 0; --jpi) {
      // Pion
      AzTrack *trPi = (AzTrack*) itsTracks->UncheckedAt(vecPi[jpi]);
      AzTrack trCor = *trPi;
      Double_t chi2 = TMath::Min (trPi->GetChi2Vertex(),999.);
      if (chi2 < gC2pi) vecPi.erase(vecPi.begin()+jpi);
    }
  }
}

//____________________________________________________________________________
//
// Adds to chain files with names like "file_1.root" - ... - "file_10.root"
// or "file-1.root" - ... - "file-10.root"
//
TChain* Chain(Int_t nFiles, TString firstFile)
{
    // Get first file number
  Int_t leng = firstFile.Length(), i1 = 0, i2 = 0;
  //cout << leng << endl;
  TString numb, prefix, suffix, symb, symb0;
  //cout << numb.Length() << endl;
  for (Int_t i = leng-1; i > -1; --i) {
    symb = TString(firstFile(i,1));
    if (symb == "_" || symb == "-") {
      prefix = firstFile(0,i+1);
      i1 = i + 1;
      break;
    } else if (symb == ".") {
      suffix = firstFile(i,leng-i);
      i2 = i - 1;
    }
  }
  numb = TString(firstFile(i1,i2-i1+1));

  Int_t numb0 = numb.Atoi();
  cout << numb << endl;
  cout << numb0 << endl;
  cout << prefix << endl;
  cout << suffix << endl;

  TChain *chain = new TChain("cbmsim");
  TString fileName;
  nFiles += numb0;
  for (Int_t i = numb0; i < nFiles; ++i) {
    fileName = prefix;
    fileName += i;
    fileName += suffix;
    if (!gSystem->FindFile("./",fileName)) break;
    chain->AddFile(fileName);
  }
  chain->ls();
  return chain;
}

//__________________________________________________________________________
//
//           (\__/)
//           (='.'=) 
//        [:]||||||||||[:]          Happy End
//           (")_(")                            
//__________________________________________________________________________

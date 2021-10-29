#if !defined(__CINT__) || defined(__MAKECINT__)
// MPD includes
#include "TpcPoint.h"
#include "MpdTpcKalmanTrack.h"
#include "MpdKalmanFilter.h"
#include "MpdKfV0Fitter.h"
//#include "MpdStrawendcapPoint.h"
//#include "MpdEctKalmanTrack.h"

// CBM includes
#include "FairMCPoint.h"
#include "FairMCTrack.h"
#include "FairRunAna.h"

// ROOT includes
#include <TBranch.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TDatabasePDG.h>
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
#include <Riostream.h>

#include <map>
#include <vector>

#endif

FILE *lun = 0x0; //fopen("ids.dat","w");
const Double_t emass = 0.00051;
Int_t pdgCodeProt = 2212; // proton
Int_t pdgCodeK = 321; // Kaon
Int_t pdgCodeE = 11; // electron

vector<MpdTpcKalmanTrack*> bkgVec;
TH1D *hBkgMult;

Float_t mass, massgen, pt, etaf, angle, etas[2], pts[2], zimp[2], dca[2];
Int_t iev, mothid, charge, orig, nps[2];

//__________________________________________________________________________
Int_t ReadBkg(Int_t iseed)
{
  // Read selected tracks from UrQMD

  TFile fileBkg("reco_100_ang_conv.histo.root");
  //TH1D *h = (TH1D*) fileBkg.Get("hElNcutAng");
  TH1D *h = (TH1D*) fileBkg.Get("hElSel");
  hBkgMult = (TH1D*) h->Clone("hBkgMult");
  hBkgMult->SetDirectory(0);
  TTree *simBkg = (TTree*) fileBkg.Get("cbmsim");
  TClonesArray *bkgTracks = (TClonesArray*) fileBkg.FindObjectAny("TpcKalmanTrack");
  simBkg->SetBranchAddress("TpcKalmanTrack",&bkgTracks);

  Int_t events = simBkg->GetEntries();
  cout << " Number of UrQMD events: " << events << endl;

  for (Int_t i = 0; i < events; ++i) {
    simBkg->GetEntry(i);
    Int_t nTracks = bkgTracks->GetEntriesFast();
    //cout << " Event: " << i << ", tracks: " << nTracks << endl;
    if (nTracks < 1) continue;
    for (Int_t j = 0; j < nTracks; ++j) {
      MpdTpcKalmanTrack *tr = (MpdTpcKalmanTrack*) bkgTracks->UncheckedAt(j);
      MpdTpcKalmanTrack *trCopy = new MpdTpcKalmanTrack(*tr);
      bkgVec.push_back(trCopy);
    }
  }
  fileBkg.Close();
  //gRandom->SetSeed(time(NULL));
  gRandom->SetSeed(iseed);
  return events;
}

//__________________________________________________________________________
void AddBkg(map<Int_t,Int_t> ids, TClonesArray *tpcTracks) 
{
  // Add tracks from UrQMD

  //Int_t nAdd = Int_t (hBkgMult->GetRandom() + 0.5);
  Int_t nAdd = Int_t (hBkgMult->GetRandom() + 0.49999);
  TH1D *h = (TH1D*) gROOT->FindObjectAny("hBkgAdd");
  h->Fill(nAdd);

  if (nAdd == 0) return; 

  Int_t nBkgTot = bkgVec.size();
  Int_t *inds = new Int_t [nAdd];
  for (Int_t i = 0; i < nAdd; ++i) {
    inds[i] = gRandom->Integer(nBkgTot);
    //cout << iev << " " << nAdd << " " << i << " " << inds[i] << endl;
  }

  TTree *tree = (TTree*) gROOT->FindObjectAny("pairs");

  TLorentzVector lorVec[2];
  TVector3 mom, momRec;
  // Combinations of bkg. tracks
  for (Int_t i = 0; i < nAdd; ++i) {
    Int_t indx = inds[i];
    MpdTpcKalmanTrack *urqmd = bkgVec[indx];
    Int_t iCh = urqmd->Charge(); 
    momRec = urqmd->Momentum3();
    lorVec[0].SetXYZM(momRec.X(), momRec.Y(), momRec.Z(), emass);
    etas[0] = momRec.Eta();
    pts[0] = momRec.Pt();
    nps[0] = urqmd->GetNofHits();
    zimp[0] = TMath::Abs(urqmd->GetParam(1));
    dca[0] = urqmd->GetPosNew();

    for (Int_t j = i + 1; j < nAdd; ++j) {
      Int_t indx1 = inds[j];
      MpdTpcKalmanTrack *urqmd1 = bkgVec[indx1];
      Int_t iCh1 = urqmd1->Charge(); 
      mom = urqmd1->Momentum3();
      lorVec[1].SetXYZM(mom.X(), mom.Y(), mom.Z(), emass);
      TLorentzVector lorInv = lorVec[0] + lorVec[1];
      Double_t invMass = lorInv.Mag();
      etas[1] = mom.Eta();
      pts[1] = mom.Pt();
      nps[1] = urqmd1->GetNofHits();
      zimp[1] = TMath::Abs(urqmd1->GetParam(1));
      dca[1] = urqmd1->GetPosNew();
      angle = TMath::RadToDeg() * lorVec[0].Angle(lorVec[1].Vect());
      charge = iCh + iCh1;
      orig = -2;
      mothid = 0;
      mass = invMass;
      etaf = lorInv.Eta();
      pt = lorInv.Pt();

      tree->Fill();
    }
  }


  map<Int_t,Int_t>::iterator it;
  for (Int_t i = 0; i < nAdd; ++i) {
    Int_t indx = inds[i];
    MpdTpcKalmanTrack *urqmd = bkgVec[indx];
    Int_t iCh = urqmd->Charge(); 
    momRec = urqmd->Momentum3();
    lorVec[0].SetXYZM(momRec.X(), momRec.Y(), momRec.Z(), emass);
    etas[0] = momRec.Eta();
    pts[0] = momRec.Pt();
    nps[0] = urqmd->GetNofHits();
    zimp[0] = TMath::Abs(urqmd->GetParam(1));
    dca[0] = urqmd->GetPosNew();

    for (it = ids.begin(); it != ids.end(); ++it) {
      if (it->second < 0) continue;
      MpdTpcKalmanTrack *tr = (MpdTpcKalmanTrack*) tpcTracks->UncheckedAt(it->second);
      Int_t iCh1 = tr->Charge(); 
      mom = tr->Momentum3();
      //if (tr->GetNofHits() < 20 || TMath::Abs(mom.Eta()) > 1.2 ||
      //  tr->GetPosNew() > 0.354*3 || TMath::Abs(tr->GetParam(1)) > 0.457*3) continue;
      lorVec[1].SetXYZM(mom.X(), mom.Y(), mom.Z(), emass);
      TLorentzVector lorInv = lorVec[0] + lorVec[1];
      Double_t invMass = lorInv.Mag();
      etas[1] = mom.Eta();
      pts[1] = mom.Pt();
      nps[1] = tr->GetNofHits();
      zimp[1] = TMath::Abs(tr->GetParam(1));
      dca[1] = tr->GetPosNew();
      angle = TMath::RadToDeg() * lorVec[0].Angle(lorVec[1].Vect());
      charge = iCh + iCh1;
      orig = -1;
      mothid = 0;
      mass = invMass;
      etaf = lorInv.Eta();
      pt = lorInv.Pt();

      tree->Fill();
    }
  }
  delete [] inds;
}

//__________________________________________________________________________
void AnalDilep(Int_t n1 = 0, Int_t n2 = 0, Int_t iDataSet = 1)
{
  // Analyze mixed UrQMD + Pluto events

  // Load basic libraries
  //gROOT->ProcessLine(".x /1/zinchenk/zinch/new/loadlibs.C");
  //gROOT->ProcessLine(".x ../loadlibs.C");
  
  // Read background tracks
  Int_t nBkg = 0; //ReadBkg(iDataSet);
  cout << " " << nBkg << " " << bkgVec.size() << endl; 

  TFile fileMC("/opt/exp_soft/mpd/zinchenk/reco/dstAuAu_25gev_ee_central_301.root");

  gROOT->ProcessLine(".L ./Chain1.C");

  if (iDataSet == 1) gROOT->ProcessLine("simTPC=Chain(100,\"/scr/u/zinchenk/dst/dstAuAu_25gev_ee_central_1.root\")");
  else if (iDataSet == 2) gROOT->ProcessLine("simTPC=Chain(100,\"/scr/u/zinchenk/dst/dstAuAu_25gev_ee_central_101.root\")");
  else if (iDataSet == 3) gROOT->ProcessLine("simTPC=Chain(100,\"/opt/exp_soft/mpd/zinchenk/reco/dstAuAu_25gev_ee_central_201.root\")");
  else gROOT->ProcessLine("simTPC=Chain(100,\"/opt/exp_soft/mpd/zinchenk/reco/dstAuAu_25gev_ee_central_301.root\")");
  TChain *simTPC = (TChain*) gROOT->FindObject("cbmsim");
  TChain *simMC = simTPC;

  //TClonesArray *tpcTracks = (TClonesArray*) fileTPC.FindObjectAny("TpcKalmanTrack");
  TClonesArray *tpcTracks = (TClonesArray*) fileMC.FindObjectAny("TpcKalmanTrack");
  simTPC->SetBranchAddress("TpcKalmanTrack",&tpcTracks);
  TBranch *tpcRecoB = simTPC->GetBranch("TpcKalmanTrack");
  //TClonesArray *tpcPoints = (TClonesArray*) fileMC.FindObjectAny("TpcPoint");
  //simMC->SetBranchAddress("TpcPoint",&tpcPoints);
  //TBranch *tpcSimB = simMC->GetBranch("TpcPoint");
  TClonesArray *mcTracks = (TClonesArray*) fileMC.FindObjectAny("MCTrack");
  simMC->SetBranchAddress("MCTrack",&mcTracks);
  TBranch *mcBranch = simMC->GetBranch("MCTrack");

  TString fout = "tpc";
  //fout += iDataSet;
  fout += ".histo.root";
  //TFile out("tpc.histo.root","recreate");
  TFile out(fout,"recreate");
  // Book histos
  TH2D *hDedx = new TH2D("hDedx","dE/dx vs P",600,0,3,100,0,10);
  TH2D *hDedxGen = new TH2D("hDedxGen","dE/dx vs Pgen",600,0,3,100,0,1.e-5);
  TH2D *hDedxPi = new TH2D("hDedxPi","dE/dx vs P for pions",600,0,3,100,0,10);
  TH1D *hInvMass = new TH1D("hInvMass","e+e- invariant mass",100,0.,1.5);
  TH1D *hInvMassEe = new TH1D("hInvMassEe","e+e- invariant mass",100,0.,1.5);
  TH1D *hInvMassEeTrue = new TH1D("hInvMassEeTrue","e+e- invariant mass",100,0.,1.5);
  TH1D *hDpt = new TH1D("hDpt","dPt/Pt",100,-1.,1.);
  TH2D *hDptVsPt = new TH2D("hDptVsPt","dPt/Pt vs Pt",200,0,2,100,-1.,1.);
  TH2D *hDptVsEta = new TH2D("hDptVsEta","dPt/Pt vs Eta",100,-3,3,100,-1.,1.);
  TH2D *hDptVsNh = new TH2D("hDptVsNh","dPt/Pt vs number of hits",120,0,60,100,-1.,1.);
  TH2D *hImpVsZp = new TH2D("hImpVsZp","DCA vs Z",100,-50,50,100,-10.,50.);
  TH2D *hImpVsZs = new TH2D("hImpVsZs","DCA vs Z",100,-50,50,100,-10.,50.);
  TH2D *hPtVsM = new TH2D("hPtVsM","Source Pt vs inv. mass",100,0,1.5,100,0.,1.5);
  TH2D *hAngVsM = new TH2D("hAngVsM","Angle between e+e- vs inv. mass",100,0,1.5,100,-10,190);
  TH2D *hEtaVsM = new TH2D("hEtaVsM","Source eta vs inv. mass",100,0,1.5,100,-3,3);
  TH2D *hPtVsAng = new TH2D("hPtVsAng","Source Pt vs e+e- angle",100,-10,190,100,0.,1.5);
  TH2D *hPtVsMtrue = new TH2D("hPtVsMtrue","Source Pt vs inv. mass",100,0,1.5,100,0.,1.5);
  TH2D *hAngVsMtrue = new TH2D("hAngVsMtrue","Angle between e+e- vs inv. mass",100,0,1.5,100,-10,190);
  TH2D *hEtaVsMtrue = new TH2D("hEtaVsMtrue","Source eta vs inv. mass",100,0,1.5,100,-3,3);
  TH2D *hPtVsAngTrue = new TH2D("hPtVsAngTrue","Source Pt vs e+e- angle",100,-10,190,100,0.,1.5);
  TH2D *hPtVsAngTrueM = new TH2D("hPtVsAngTrueM","Source Pt vs e+e- angle",100,-10,190,100,0.,1.5);
  TH1D *hBkgAdd = new TH1D("hBkgAdd","Added background",100,-9.5,90.5);
  TH2D *hPerpA0 = new TH2D("hPerpA0","Angle v.r.t. transverse plane",100,-5,45,100,-1,1);
  TH2D *hPerpA1 = new TH2D("hPerpA1","Angle v.r.t. transverse plane",100,-5,45,100,-1,1);
  TH1D *hTrue = new TH1D("hTrue","True pairs ",30,-4.5,25.5);
  TH1D *hFalse = new TH1D("hFalse","False pairs ",30,-4.5,25.5);
  TH1D *hSel = new TH1D("hSel","Selected electrons ",30,-4.5,25.5);

  //Float_t mass, pt, etaf, angle, etas[2], pts[2], zimp[2], dca[2];
  //Int_t charge, orig, mothid, nps[2];

  TTree *tree = new TTree("pairs","Pairs");
  tree->Branch("mass",&mass,"mass/F");
  tree->Branch("massgen",&massgen,"massgen/F");
  tree->Branch("pt",&pt,"pt/F");
  tree->Branch("eta",&etaf,"etaf/F");
  tree->Branch("angle",&angle,"angle/F");
  tree->Branch("etas[2]",&etas,"etas[2]/F");
  tree->Branch("pts[2]",&pts,"pts[2]/F");
  tree->Branch("zimp[2]",&zimp,"zimp[2]/F");
  tree->Branch("dca[2]",&dca,"dca[2]/F");

  tree->Branch("iev",&iev,"iev/I");
  tree->Branch("moth",&mothid,"mothid/I");
  tree->Branch("charge",&charge,"charge/I2");
  tree->Branch("orig",&orig,"orig/I2");
  tree->Branch("nps[2]",&nps,"nps[2]/I2");

  TString inFile = "mc.root";
  TString outFile = "aaa.root";
  FairRunAna *fRun= new FairRunAna();
  //fRun->SetInputFile(inFile);
  fRun->SetInputFile(&fileMC);
  fRun->SetOutputFile(outFile);
  MpdKalmanFilter *pKF = MpdKalmanFilter::Instance("KF");
  pKF->Init();

  Int_t events = simTPC->GetEntries();
  if (n2 != 0) events = TMath::Min (events, n2);
  cout << " Number of events = " << events << endl;

  Int_t evNo = 0;
  for (Int_t i = 0; i < events; ++i) {
    if (i < n1) continue;
    simMC->GetEntry(i);
    evNo = i + 1;
    iev = i;

    // MC tracks
    Int_t nMC = mcTracks->GetEntriesFast(), nEl = 0, nElPrim = 0;
    for (Int_t j = 0; j < nMC; ++j) {
      FairMCTrack* mcTr = (FairMCTrack*) mcTracks->UncheckedAt(j);
      Int_t pdg = mcTr->GetPdgCode();
      if (TMath::Abs(pdg) != 11) continue;
      ++nEl;
      if (mcTr->GetMotherId() < 0) ++nElPrim;
    }
    //cout << " Electrons: " << nEl << " " << nElPrim << endl;

    // Find TPC reco tracks 
    Int_t nTPC = tpcTracks->GetEntriesFast();
    if (i % 500 == 0) cout << " Event No: " << i << ", reco tracks in TPC: " << nTPC << endl;
    //cout << " Event No: " << i << ", reco tracks in TPC: " << nTPC << endl;
    
    // Check for multiple reco tracks with the same ID
    map<Int_t, Int_t> ids;

    for (Int_t j = 0; j < nTPC; ++j) {
      MpdTpcKalmanTrack *tr = (MpdTpcKalmanTrack*) tpcTracks->UncheckedAt(j);
      Int_t id = tr->GetTrackID();
      if (ids.find(id) == ids.end()) {
	ids.insert(pair<Int_t,Int_t>(id,j));
	continue;
      }
      // Multiple tracks
      MpdTpcKalmanTrack *tr1 = (MpdTpcKalmanTrack*) tpcTracks->UncheckedAt(ids[id]);
      /*
      cout << " Mult " << id << " " << ids[id] << " " << j << " " 
	   << tr1->GetParam(1) << " " << tr->GetParam(1) << " " 
	   << tr1->GetNofWrong() << " " << tr->GetNofWrong() << " " 
	   << tr1->Pt() << " " << tr->Pt() << " " 
	   << tr1->GetNofHits() << " " << tr->GetNofHits() << " " 
	   << tr1->Momentum3().Eta() << " " << tr->Momentum3().Eta() << endl;
      */
      Int_t nHits = tr->GetNofTrHits();
      FairMCPoint *p = 0x0;
      Double_t z = 0, z1 = 0;
      for (Int_t k = nHits-1; k > -1; --k) {
	MpdKalmanHit *hit = (MpdKalmanHit*) tr->GetTrHits()->UncheckedAt(k);
	//p = (FairMCPoint*) tpcPoints->UncheckedAt(hit->GetIndex());
	//if (p->GetTrackID() == id) break;
	z = TMath::Abs (hit->GetMeas(1));
	break;
      }
      //z = TMath::Abs (p->GetZ());

      nHits = tr1->GetNofTrHits();
      for (Int_t k = nHits-1; k > -1; --k) {
        MpdKalmanHit *hit = (MpdKalmanHit*) tr1->GetTrHits()->UncheckedAt(k);
        //p = (FairMCPoint*) tpcPoints->UncheckedAt(hit->GetIndex());
        //if (p->GetTrackID() == id) break;
	z1 = TMath::Abs (hit->GetMeas(1));
	break;
      }
      //z1 = TMath::Abs (p->GetZ());

      if (z1 > z) { 
	ids.erase(id);
        ids.insert(pair<Int_t,Int_t>(id,j));
      }
    }

    map<Int_t,Int_t>::iterator it, it1;
    TLorentzVector lorVec[2], lorVecG[2];
    TVector3 momGen[2], momRec[2], zAxis(0,0,1), vtx;
    Int_t nTrue = 0, nFalse = 0;
    for (it = ids.begin(); it != ids.end(); ++it) {
      if (it->second < 0) continue;
      MpdTpcKalmanTrack *tr = (MpdTpcKalmanTrack*) tpcTracks->UncheckedAt(it->second);
      Int_t iCh = tr->Charge(); 
      Int_t id = tr->GetTrackID();
      FairMCTrack* mcTr = (FairMCTrack*) mcTracks->UncheckedAt(id);
      Int_t moth = mcTr->GetMotherId();
      mcTr->GetMomentum(momGen[0]);
      //if (iCh < 0) continue;
      Int_t pdg = TMath::Abs (mcTr->GetPdgCode());
      //if (pdg == pdgCodeProt || pdg == pdgCodeK) continue; // exclude proton and kaon
      momRec[0] = tr->Momentum3();

      // Cuts
      //if (tr->GetNofHits() < 20 || TMath::Abs(momRec.Eta()) > 1.2 ||
      //  tr->GetPosNew() > 0.354*3 || TMath::Abs(tr->GetParam(1)) > 0.457*3) {
      if (tr->GetNofHits() < 5 || TMath::Abs(momRec[0].Eta()) > 1.5 ||
	  tr->GetPosNew() >= 0.354*9 || TMath::Abs(tr->GetParam(1)) >= 0.457*9) {
	it->second = -1;
	continue;
      }

      lorVec[0].SetXYZM(momRec[0].X(), momRec[0].Y(), momRec[0].Z(), emass);
      lorVecG[0].SetXYZM(momGen[0].X(), momGen[0].Y(), momGen[0].Z(), emass);
      pt = momGen[0].Pt();
      etas[0] = momRec[0].Eta();
      pts[0] = momRec[0].Pt();
      nps[0] = tr->GetNofHits();
      zimp[0] = TMath::Abs(tr->GetParam(1));
      dca[0] = tr->GetPosNew();
      hDpt->Fill((tr->Pt()-pt)/pt);
      hDptVsPt->Fill(pt,(tr->Pt()-pt)/pt);
      hDptVsEta->Fill(momGen[0].Eta(),(tr->Pt()-pt)/pt);
      hDptVsNh->Fill(tr->GetNofHits()+0.1,(tr->Pt()-pt)/pt);
      if (moth < 0) hImpVsZp->Fill(tr->GetParam(1),tr->GetPosNew()); // primary
      else hImpVsZs->Fill(tr->GetParam(1),tr->GetPosNew());
      //for (it1 = ids.begin(); it1 != ids.end(); ++it1) {
      for (it1 = it; it1 != ids.end(); ++it1) {
	if (it1 == it) continue;
	if (it1->second < 0) continue;
	MpdTpcKalmanTrack *tr1 = (MpdTpcKalmanTrack*) tpcTracks->UncheckedAt(it1->second);
	Int_t iCh1 = tr1->Charge();
	//if (iCh1 == iCh) continue;
	id = tr1->GetTrackID();
	mcTr = (FairMCTrack*) mcTracks->UncheckedAt(id);
	Int_t moth1 = mcTr->GetMotherId();
	mcTr->GetMomentum(momGen[1]);
	Int_t pdg1 = TMath::Abs (mcTr->GetPdgCode());
	//if (pdg1 == pdgCodeProt || pdg1 == pdgCodeK) continue; // exclude proton and kaon
	momRec[1] = tr1->Momentum3();

	// Cuts
	//if (tr1->GetNofHits() < 20 || TMath::Abs(mom.Eta()) > 1.2 ||
	//tr1->GetPosNew() > 0.354*3 || TMath::Abs(tr1->GetParam(1)) > 0.457*3) {
	if (tr1->GetNofHits() < 5 || TMath::Abs(momRec[1].Eta()) > 1.5 ||
	    tr1->GetPosNew() >= 0.354*9 || TMath::Abs(tr1->GetParam(1)) >= 0.457*9) {
	  it1->second = -1;
	  continue;
	}

	// Cut conversion pairs
	if (iCh != iCh1) {
	  MpdTpcKalmanTrack track = *tr;
	  MpdTpcKalmanTrack track1 = *tr1;
	  Double_t chi2 = MpdKfV0Fitter::Instance()->FindVertex(&track, &track1, vtx);
	  TVector3 cross = momRec[0].Cross(momRec[1]);
	  cross *= (1./cross.Mag());
	  Double_t orient = cross * zAxis * TMath::Sign(1,-iCh);
	  if (moth >=0 && moth1 == moth) hPerpA0->Fill(vtx.Pt(),orient);
	  else hPerpA1->Fill(vtx.Pt(),orient);
	  if (orient > 0.5 && vtx.Pt() > 4) {
	    // "Conversion" pair
	    it->second = -1;
	    it1->second = -1;
	    if (moth >=0 && moth1 == moth) ++nTrue;
	    else ++nFalse;
	    break;
	  }
	}

	// Cuts
	if (tr->GetNofHits() < 20 || TMath::Abs(momRec[0].Eta()) > 1.2 ||
	    tr->GetPosNew() > 0.354*3 || TMath::Abs(tr->GetParam(1)) > 0.457*3) continue;
	if (tr1->GetNofHits() < 20 || TMath::Abs(momRec[1].Eta()) > 1.2 ||
	    tr1->GetPosNew() > 0.354*3 || TMath::Abs(tr1->GetParam(1)) > 0.457*3) continue;

	lorVec[1].SetXYZM(momRec[1].X(), momRec[1].Y(), momRec[1].Z(), emass);
	lorVecG[1].SetXYZM(momGen[1].X(), momGen[1].Y(), momGen[1].Z(), emass);
        TLorentzVector lorInv = lorVec[0] + lorVec[1];
        TLorentzVector lorInvG = lorVecG[0] + lorVecG[1];
	Double_t invMass = lorInv.Mag();
        hInvMass->Fill(invMass);
	etas[1] = momRec[1].Eta();
	pts[1] = momRec[1].Pt();
	nps[1] = tr1->GetNofHits();
	zimp[1] = TMath::Abs(tr1->GetParam(1));
	dca[1] = tr1->GetPosNew();
	angle = TMath::RadToDeg() * lorVec[0].Angle(lorVec[1].Vect());
	charge = iCh + iCh1;
	orig = 0;
	mothid = 0;
	mass = invMass;
	massgen = lorInvG.Mag();
	etaf = lorInv.Eta();
	pt = lorInv.Pt();
	if (pdg == pdgCodeE && pdg1 == pdgCodeE) {
	  hInvMassEe->Fill(invMass);
	  hPtVsM->Fill(invMass, lorInv.Pt());
	  hAngVsM->Fill(invMass, TMath::RadToDeg()*lorVec[0].Angle(lorVec[1].Vect()));
	  hEtaVsM->Fill(invMass, lorInv.Eta());
	  hPtVsAng->Fill(TMath::RadToDeg()*lorVec[0].Angle(lorVec[1].Vect()), lorVec[1].Pt());
	  if (moth == moth1 && moth < 0 && charge == 0) {
	    // True combination
	    mothid = moth;
	    orig = 1;
	    hInvMassEeTrue->Fill(invMass);
	    hPtVsMtrue->Fill(invMass, lorInv.Pt());
	    hAngVsMtrue->Fill(invMass, TMath::RadToDeg()*lorVec[0].Angle(lorVec[1].Vect()));
	    hEtaVsMtrue->Fill(invMass, lorInv.Eta());
	    hPtVsAngTrue->Fill(TMath::RadToDeg()*lorVec[0].Angle(lorVec[1].Vect()), lorVec[1].Pt());
	    if (moth < -40)  
	      hPtVsAngTrueM->Fill(TMath::RadToDeg()*lorVec[0].Angle(lorVec[1].Vect()), lorVec[1].Pt());
	  }
	}
	tree->Fill();
      } // for (it1 = it; it1 != ids.end();
    } // for (it = ids.begin(); it != ids.end();
    if (nTrue) hTrue->Fill(nTrue);
    if (nFalse) hFalse->Fill(nFalse);
    
    // Track quality cuts
    Int_t nSel = 0;
    for (it = ids.begin(); it != ids.end(); ++it) {
      if (it->second < 0) continue;
      MpdTpcKalmanTrack *tr = (MpdTpcKalmanTrack*) tpcTracks->UncheckedAt(it->second);
      if (tr->GetNofHits() < 20 || TMath::Abs(tr->Momentum3().Eta()) > 1.2 ||
	  tr->GetPosNew() > 0.354*3 || TMath::Abs(tr->GetParam(1)) > 0.457*3) {
	it->second = -1;
	continue;
      }
      ++nSel;
    }
    hSel->Fill(nSel);

    //AddBkg(ids, tpcTracks);

  } // for (Int_t i = 0; i < events;
  if (lun) fclose(lun);
  
  out.Write();
  out.Close();
}

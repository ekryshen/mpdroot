/// \ingroup physics
/// \class MpdAnaXiMix
/// \brief Lambda and Xi reconstruction with event mixing
///
/// \author Alexander Zinchenko, LHEP JINR Dubna - 21-03-2023

#include "MpdAnaXiMix.h"
#include "TpcPoint.h"
#include "MpdLambda.h"
#include "MpdKalmanFilter.h"
#include "MpdMCTrack.h"
#include "MpdItsKalmanTrack.h"
#include "MpdTrackFinderIts5spd.h"
#include "MpdTpcKalmanFilter.h"
#include "MpdTpcKalmanTrack.h"
//#include "MpdEmcClusterKI.h"
#include "MpdParticle.h"
#include "MpdPid.h"
#include "MpdTofMatchingData.h"
#include "MpdVertex.h"
#include "MpdEvent.h"
//#include "MpdEmcGeoUtils.h"

#include "FairFileSource.h"
#include "FairMCPoint.h"
#include "FairParRootFileIo.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"

#include "TFile.h"
#include "TGeoManager.h"
#include <TMinuit.h>

#include <iostream>
#include <fstream> // std::ifstream
#include <set>     //
#include <vector>  //

TVector3  fVtxN, fMomN;
MpdHelix *fTrC;

ClassImp(MpdAnaXiMix);

//-----------------------------------------------------------------------------------------------

MpdAnaXiMix::MpdAnaXiMix() : MpdAnalysisTask(), fRecoIts(nullptr) {}

//-----------------------------------------------------------------------------------------------

MpdAnaXiMix::MpdAnaXiMix(const char *name, const char *outputName)
   : MpdAnalysisTask(name, outputName), fRecoIts(nullptr)
{
   // mParamConfig = outputName;
}

//-----------------------------------------------------------------------------------------------

void MpdAnaXiMix::UserInit()
{
   // Prepare histograms etc.

   fOutputList = new TList();
   fOutputList->SetOwner(kTRUE);
   TH1F *hp = nullptr;

   hp                      = new TH1F("hLambFlag", "Flags for lambda", 14, 0, 14);
   fHistosF[hp->GetName()] = hp;
   hp                      = new TH1F("hXiFlag", "Flags for Xi", 14, 0, 14);
   fHistosF[hp->GetName()] = hp;
   hp                      = new TH1F("hVertGen", "True event vertex distribution", 100, -200., 200.);
   fHistosF[hp->GetName()] = hp;
   hp                      = new TH1F("hPIDflag", "PID flags", 12, 0, 12);
   fHistosF[hp->GetName()] = hp;

   hp                      = new TH1F("hMassL", "Lambda mass", 50, 1.070, 1.170);
   fHistosF[hp->GetName()] = hp;
   hp                      = new TH1F("hMassLsig", "Lambda mass (signal)", 50, 1.070, 1.170);
   fHistosF[hp->GetName()] = hp;
   hp                      = new TH1F("hMassLbkg", "Lambda mass (bkg.)", 50, 1.070, 1.170);
   fHistosF[hp->GetName()] = hp;

   hp                      = new TH1F("hMassXi", "Xi- mass", 50, 1.260, 1.360);
   fHistosF[hp->GetName()] = hp;
   hp                      = new TH1F("hMassXiSig", "Xi- mass (signal)", 50, 1.260, 1.360);
   fHistosF[hp->GetName()] = hp;
   hp                      = new TH1F("hMassXiBkg", "Xi- mass (bkg.)", 50, 1.260, 1.360);
   fHistosF[hp->GetName()] = hp;

   hp                      = new TH1F("hMassOm", "Omega- mass", 50, 1.630, 1.730);
   fHistosF[hp->GetName()] = hp;
   hp                      = new TH1F("hMassOmSig", "Omega- mass (signal)", 50, 1.630, 1.730);
   fHistosF[hp->GetName()] = hp;
   hp                      = new TH1F("hMassOmBkg", "Omega- mass (bkg.)", 50, 1.630, 1.730);
   fHistosF[hp->GetName()] = hp;

   hp                      = new TH1F("hPtProt", "Proton Pt", 20, 0, 5);
   fHistosF[hp->GetName()] = hp;
   hp                      = new TH1F("hPtProtT", "True Proton Pt", 20, 0, 5);
   fHistosF[hp->GetName()] = hp;
   hp                      = new TH1F("hPtProtF", "False Proton Pt", 20, 0, 5);
   fHistosF[hp->GetName()] = hp;

   for (auto it = fHistosF.begin(); it != fHistosF.end(); ++it) fOutputList->Add(it->second);

   fvvvL    = &fvLambdas;
   fvvvXi   = &fvXis;
   fvvvLpt  = &fvLambMpdgPtEtaY;
   fvvvXipt = &fvXiMpdgPtEtaY;
   fEvNo    = -1;

   //*
   fTree = new TTree("event", "Event");
   fTree->Branch("evNo", &fEvNo, "evNo/I");
   fTree->Branch("b0", &fB0, "b0/F");
   fTree->Branch("centr", &fCentr, "centr/F");
   fTree->Branch("zv", &fZv, "zv/F");
   fTree->Branch("zvgen", &fZvgen, "zvgen/F");
   fTree->Branch("ntr13", &fNtr13, "ntr13/I");
   fTree->Branch("nv", &fNv, "nv/I");
   TBranch *br = fTree->Branch("l0", "std::vector<MpdLambda>", &fvvvL);
   // br->Print();
   fTree->Branch("xi", "std::vector<MpdXi>", &fvvvXi);
   fTree->Branch("m_pt_eta_y_l0", "std::vector<tuple<int,float,float,float> >", &fvvvLpt);
   fTree->Branch("m_pt_eta_y_xi", "std::vector<tuple<int,float,float,float> >", &fvvvXipt);
   //*/
   fOutputList->Add(fTree);

   Double_t sigM = 3.0, sigE = 3.0, energy = 9.2, coef = 1.0; // n-sigma bands for PID selection
   // Double_t sigM = 4.0, sigE = 4.0, energy = 11.0, coef = 1.0; // n-sigma bands for PID selection
   // TString generator = "PHSD", tracking = "CF";
   TString generator = "NSIG", tracking = "CFHM";
   fPid = new MpdPid(sigM, sigE, energy, coef, generator, tracking, "pikaprdetrhe3he4");

   fTrC    = new MpdHelix(0, 0, 0, TVector3(0, 0, 0), 0);
   fMinuit = new TMinuit(2);

   cout << "[MpdAnaXiMix] -- Successful start --" << endl << endl;
}

//-----------------------------------------------------------------------------------------------

void MpdAnaXiMix::ProcessEvent(MpdAnalysisEvent &event)
{

   if (!fisInitialized) {
      // mKF = MpdKalmanFilter::Instance();
      // mKHit.SetType(MpdKalmanHit::kFixedR);
      // MpdKalmanFilter::Instance()->Init();
      TFile f("sim_Geo.root"); // V
      f.Get("FairGeoParSet");  // V
      /*
      FairRunAna &ana = *FairRunAna::Instance();
      //AZ-020822 - for magnetic field
      FairFileSource *fileSource = new FairFileSource("mc_0.root");
      ana.SetSource(fileSource);
      FairRuntimeDb* rtdb = ana.GetRuntimeDb();
      FairParRootFileIo* parInput1 = new FairParRootFileIo();
      parInput1->open("mc_0.root");
      rtdb->setFirstInput(parInput1);
      ana.SetOutputFile("/dev/null"); //AZ-230323
      ana.Init();
      */
      cout << "[MpdAnaXiMix] GeoMan: " << gGeoManager << endl;
      // MpdKalmanFilter::Instance()->Init();
      BaseTpcSectorGeo *secGeo = new TpcSectorGeoAZ();
      fRecoTpc                 = new MpdTpcKalmanFilter(*secGeo, "Kalman filter");
      fRecoTpc->FillGeoScheme();
      fisInitialized = true;
   }

   /*
   if (!selectEvent(event)) { //(V)
      return;
   }

   mKalmanTracks = event.fTPCKalmanTrack;

   if (isMC) {
      mMCTracks = event.fMCTrack;

      for (int i = 0; i < mMCTracks->GetEntriesFast(); i++) {
         MpdMCTrack *pr = (static_cast<MpdMCTrack *>(mMCTracks->At(i)));
         if (pr->GetPdgCode() == 333) {
            if (pr->GetStartX() * pr->GetStartX() + pr->GetStartY() * pr->GetStartY() < 1.) {
               TVector3 momentum;
               pr->GetMomentum(momentum);
               mInvGen->Fill(momentum.Pt());
               mInvGenBin[anaBin]->Fill(momentum.Pt());
            }
         }
      }
   } // isMC

   selectPosTrack(event);

   selectNegTrack(event);

   processHistograms(event);
   */

   if (!SelectEvent(event)) {
      return;
   }

   fMcTracks   = event.fMCTrack;
   fItsTracks  = event.fTPCKalmanTrack;
   fTofMatches = event.fTOFMatching;

   fB0 = fCentr = fZv = 0.0;
   fNtr13 = fNv = 0;
   fB0          = event.fMCEventHeader->GetB();
   fCentr       = event.getCentrTPC();
   fvLambMpdgPtEtaY.clear();
   fvXiMpdgPtEtaY.clear();
   fvLambdas.clear();

   SelectTracks(event);

   // MpdLambda l0;
   // fvLambdas.push_back(l0);
   // gDirectory->Print();

   fTree->Fill();
}

//-----------------------------------------------------------------------------------------------

void MpdAnaXiMix::Finish()
{
   // Post-scan processing not needed
}

//-----------------------------------------------------------------------------------------------

bool MpdAnaXiMix::SelectEvent(MpdAnalysisEvent &event)
{
   ++fEvNo;

   // Protection
   TClonesArray *vtxs = event.fVertex;
   if (vtxs == nullptr) return kFALSE;
   MpdVertex *vtx = (MpdVertex *)vtxs->First();
   if (vtx == nullptr) return kFALSE;
   if (event.fTPCKalmanTrack->GetEntriesFast() == 0) return kFALSE;

   // Check if this is an empty event (for UrQMD)
   TClonesArray *mcTracks = event.fMCTrack;
   Int_t         nNucl[2] = {0}, nMC = mcTracks->GetEntriesFast();

   for (Int_t j = 0; j < nMC; ++j) {
      // FairMCTrack* mcTr = (FairMCTrack*) mcTracks->UncheckedAt(j);
      MpdMCTrack *mcTr = (MpdMCTrack *)mcTracks->UncheckedAt(j);
      if (mcTr->GetMotherId() == -1 && mcTr->GetPdgCode() == 2212)
         ++nNucl[0]; // primary proton
      else if (mcTr->GetMotherId() == -1 && mcTr->GetPdgCode() == 2112)
         ++nNucl[1]; // primary neutron
   }
   // cout << nNucl[0] << " xxx " << nNucl[1] << endl;
   if (nNucl[0] == 166 && nNucl[1] == 252) return kFALSE; //!!! no interaction Bi+Bi !!!

   return kTRUE;
}

//-----------------------------------------------------------------------------------------------

/*
bool MpdPairKK::selectEvent(MpdAnalysisEvent &event)
{

   mhEvents->Fill(0.5);

   if (!event.fVertex) { // if even vertex not filled, skip event
      return false;
   }

   mhEvents->Fill(1.5);

   MpdVertex *vertex = (MpdVertex *)event.fVertex->First();
   vertex->Position(mPrimaryVertex);

   if (mPrimaryVertex.Z() == 0) { // not reconstructed (==0)
      return false;
   }

   if (fabs(mPrimaryVertex.Z()) > mParams.mZvtxCut) { // beyond the limits
      return false;
   }

   mhEvents->Fill(2.5);

   float cen = event.getCentrTPC();

   if (cen < 0 || cen >= 100) { // TPC centrality not defined
      return false;
   }

   mhEvents->Fill(3.5);

   mZvtxBin = 0.5 * (mPrimaryVertex.Z() / mParams.mZvtxCut + 1) * nMixEventZ;
   if (mZvtxBin < 0) mZvtxBin = 0;
   if (mZvtxBin >= nMixEventZ) mZvtxBin = nMixEventZ - 1;

   mCenBin = (cen / 100.) * nMixEventCent; // very rough
   if (mCenBin < 0) mCenBin = 0;
   if (mCenBin >= nMixEventCent) mCenBin = nMixEventCent - 1;

   mixBin = mZvtxBin * nMixEventCent + mCenBin;

   mhVertex->Fill(mPrimaryVertex.Z());
   mhCentrality->Fill(mCenBin);

   anaBin = -1;
   if (cen >= 0 && cen < 10) anaBin = 0;
   if (cen >= 10 && cen < 20) anaBin = 1;
   if (cen >= 20 && cen < 30) anaBin = 2;
   if (cen >= 30 && cen < 40) anaBin = 3;
   if (cen >= 40 && cen < 50) anaBin = 4;
   if (cen >= 50 && cen < 60) anaBin = 5;
   if (cen >= 60 && cen < 100) anaBin = 6;

   return true;
}
*/

//-----------------------------------------------------------------------------------------------

void MpdAnaXiMix::SelectTracks(MpdAnalysisEvent &event)
{
   // Select tracks in the event

   TVector3 genVert;
   event.fMCEventHeader->GetVertex(genVert);
   fZvgen                 = genVert.Z();
   TClonesArray *mcTracks = event.fMCTrack;
   Int_t         nMC      = mcTracks->GetEntriesFast();

   TVector3 mom, primVert;

   for (Int_t j = 0; j < nMC; ++j) {
      MpdMCTrack *mcTr = (MpdMCTrack *)mcTracks->UncheckedAt(j);
      // if (mcTr->GetPdgCode() == pdgCodeXi) { skip = 1; break; }
      mcTr->GetMomentum(mom);
      TVector3 pos;
      Double_t r = 0.0;
      if (mcTr->GetPdgCode() == fParams.pdgCodeXi) {
         // Check production vertex
         Int_t mpdg = -1;
         if (mcTr->GetMotherId() >= 0) {
            MpdMCTrack *moth = (MpdMCTrack *)mcTracks->UncheckedAt(mcTr->GetMotherId());
            mpdg             = moth->GetPdgCode();
         }
         mcTr->GetStartVertex(pos);
         pos -= genVert;
         r = pos.Mag();
         if (r < 50.0) {
            fHistosF["hXiFlag"]->Fill(0);
            Double_t pt = mom.Pt();
            if (mcTr->GetMotherId() == -1) pt *= -1; // negative pT for primaries
            Double_t eta = (TMath::Abs(pt) > 0.001) ? mom.Eta() : TMath::Sign(100., mom.Z());
            fvXiMpdgPtEtaY.push_back(make_tuple(mpdg, pt, eta, mcTr->GetRapidity()));
         }
      } else if (mcTr->GetPdgCode() == fParams.pdgCodeL0) {
         // Check production vertex
         Int_t mpdg = -1;
         if (mcTr->GetMotherId() >= 0) {
            MpdMCTrack *moth = (MpdMCTrack *)mcTracks->UncheckedAt(mcTr->GetMotherId());
            mpdg             = moth->GetPdgCode();
         }
         mcTr->GetStartVertex(pos);
         pos -= genVert;
         r = pos.Mag();
         if (r < 50.0) {
            // Production vertex constraint 50 cm
            fHistosF["hLambFlag"]->Fill(0);
            Double_t pt = mom.Pt();
            if (mcTr->GetMotherId() < 0) pt *= -1; // negative pT for primaries
            Double_t eta = (TMath::Abs(pt) > 0.001) ? mom.Eta() : TMath::Sign(100., mom.Z());
            fvLambMpdgPtEtaY.push_back(make_tuple(mpdg, pt, eta, mcTr->GetRapidity()));
         }
      }
   } // for (Int_t j = 0; j < nMC;

   Int_t nMpdTr = 0;
   // Int_t nITS = itsTracks->GetEntriesFast();
   TClonesArray *itsTracks = event.fTPCKalmanTrack, *vtxs = event.fVertex, *mpdTracks = nullptr;
   Int_t         nITS  = itsTracks->GetEntriesFast();
   Int_t         nVert = vtxs->GetEntriesFast();
   if (event.fMPDEvent) fMpdTracks = event.fMPDEvent->GetGlobalTracks();
   if (fMpdTracks) nMpdTr = fMpdTracks->GetEntriesFast();

   MpdVertex *vtx = (MpdVertex *)vtxs->First();
   fMpdVert       = vtx;
   vtx->Position(primVert);
   fZv            = primVert.Z();
   fNv            = vtx->GetNTracks();
   TArrayI *indxs = vtx->GetIndices();
   Int_t    nPrim = indxs->GetSize();
   set<int> indxVert;
   for (Int_t k = 0; k < nPrim; ++k) indxVert.insert((*indxs)[k]);
   if (fEvNo % 100 == -1) {
      cout << " *** Event No: " << fEvNo << ", reco tracks in TPC (ITS), global: "
           << " " << nITS << " " << nMpdTr << ", vertices: " << nVert << endl;
      cout << " Number of primary (used for vertex reco) tracks: " << indxVert.size() << endl;
   }

   // For mixing
   if (fParams.nMix && fMapVertexEvent.size() > fParams.nMix) {
      // Remove first stored event
      int iev0 = fMapVertexEvent.begin()->first;
      fMapVertexEvent.erase(iev0);
      fMapPiEvent.erase(iev0);
   }
   fMapVertexEvent[fEvNo] = *vtx;

   // Collect tracks
   fLays.clear();
   map<int, int>           ids, moths, pdgs, &lays = fLays;
   map<int, Double_t>      pts, ths, rads;
   map<int, FairMCPoint *> points;
   map<int, AzTrack *>     tracks;

   // Get max. reached layer No.
   for (Int_t j = 0; j < nITS; ++j) {
      AzTrack *tr = (AzTrack *)itsTracks->UncheckedAt(j);
      Int_t    id = tr->GetTrackID();
      ids[id]++;
      // if (i % 100 == 0 && ids[id] > 1) cout << " More than 1 reco track ID: " << id << " " << ids[id] << endl;
      MpdKalmanHit *hit = (MpdKalmanHit *)tr->GetTrHits()->First();
      if (lays.find(id) == lays.end())
         lays[id] = hit->GetLayer();
      else
         lays[id] = TMath::Max(hit->GetLayer(), lays[id]);
   }

   // Exclude "clones" (multiple loops)
   for (Int_t j = 0; j < nITS; ++j) {
      AzTrack *tr = (AzTrack *)itsTracks->UncheckedAt(j);
      Int_t    id = tr->GetTrackID();
      if (tracks.find(id) == tracks.end()) tracks[id] = tr;
      // Get ITS info
      TClonesArray *hits  = tr->GetTrHits();
      Int_t         nHits = hits->GetEntriesFast();
      FairMCPoint  *p1    = 0x0;

      for (Int_t ih = nHits - 1; ih >= 0; --ih) {
         MpdKalmanHit *hit = (MpdKalmanHit *)hits->UncheckedAt(ih);
         // if (hit->GetDist() < rMin) continue; // ITS hit
         if (hit->GetUniqueID()) continue; // ITS hit
         // if (tpcPoints) {
         if (0) {
            /*
           p1 = (FairMCPoint*) tpcPoints->UncheckedAt(hit->GetIndex());
           //cout << p1 << " " << hit->GetUniqueID() << " " << ids[id] << " " << point[id] << " " << p1->GetTrackID() <<
           endl; if (p1->GetTrackID() != id) continue; if (ids[id] > 1 && point[id]) {
             // More than 1 reco track with the same ID
             //cout << " Time: " << id << " " << p1 << " " << point[id] << " " << p1->GetTime() << " " <<
           point[id]->GetTime() << endl; if (p1 == point[id]) {
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
           */
         } else {
            // No MC points
            if (ids[id] > 1 && points[id]) {
               // More than 1 reco track with the same ID - take the one
               // closer to z = 0
               // cout << " z: " << id << " " << tr->GetParam(1) << " " << track[id]->GetParam(1) << " " << tr->Charge()
               // << " " << track[id]->Charge() << endl;
               if (TMath::Abs(tr->GetParam(1)) < TMath::Abs(tracks[id]->GetParam(1))) {
                  // Exclude previous track from further consideration
                  tracks[id]->SetChi2(-9.);
                  tracks[id] = tr;
               } else {
                  tr->SetChi2(-9.); // exclude this track from further consideration
                  break;
               }
            }
            points[id] = (FairMCPoint *)0x1;
         }
         break;
      } // for (Int_t ih = nHits-1; ih >= 0;

      // MC track
      MpdMCTrack *mcTr = (MpdMCTrack *)mcTracks->UncheckedAt(id);
      mcTr->GetMomentum(mom);
      pts[id] = mom.Pt();
      ths[id] = mom.Theta();
   } // for (Int_t j = 0; j < nITS;

   // Lambda acceptance

   multimap<Int_t, Int_t> mapLamb, mapXi;
   int                    idMax = ids.rbegin()->first;

   for (Int_t j = 0; j <= idMax; ++j) {
      MpdMCTrack *mcTr = (MpdMCTrack *)mcTracks->UncheckedAt(j);
      mcTr->GetMomentum(mom);
      Int_t mothID = mcTr->GetMotherId();
      if (mothID == -1 && lays.find(j) != lays.end() && lays[j] != 0) {
         lays[j] = -lays[j]; // flag primary tracks
                             // href->Fill(mom.Eta());
      }
      TVector3 pos;
      mcTr->GetStartVertex(pos);
      rads[j]  = pos.Pt();
      moths[j] = -1;
      pdgs[j]  = mcTr->GetPdgCode();
      if (mothID >= 0) {
         // Check lambda production vertex ( < 50 cm)
         MpdMCTrack *moth = (MpdMCTrack *)mcTracks->UncheckedAt(mothID);
         moth->GetStartVertex(pos);
         pos -= genVert;
         if (pos.Mag() < 50.0) {
            moths[j] = moth->GetPdgCode();
            if (moths[j] == fParams.pdgCodeL0 && (pdgs[j] == fParams.pdgCodePr || pdgs[j] == fParams.pdgCodeNeg))
               mapLamb.insert(pair<Int_t, Int_t>(mothID, j));
         }
         if (moths[j] == fParams.pdgCodeXi && (pdgs[j] == fParams.pdgCodeL0 || pdgs[j] == fParams.pdgCodeNeg))
            mapXi.insert(pair<Int_t, Int_t>(mothID, j));
      }
      // if ((pdgs[j] == pdgCodePos || pdgs[j] == pdgCodeNeg) && moths[j] == pdgCodeL0)
      // hPtVsEtaS->Fill(TMath::Abs(mom.Eta()),mom.Pt());
   }

   multimap<int, int>::iterator                                     mit, mit1;
   pair<multimap<int, int>::iterator, multimap<int, int>::iterator> ret;

   mit = mapLamb.begin();
   while (mit != mapLamb.end()) {
      Int_t mothID = mit->first;
      if (mapLamb.count(mothID) != 2) {
         mit = mapLamb.upper_bound(mothID);
         continue;
      } // only one decay particle
      // cout << mapLamb.count(mothID) << endl;
      ret           = mapLamb.equal_range(mothID);
      Int_t nppi[2] = {0}, nok = 0;
      Int_t nok1 = 0, nok2 = 0, nok3 = 0;

      for (mit1 = ret.first; mit1 != ret.second; ++mit1) {
         MpdMCTrack *mcTr = (MpdMCTrack *)mcTracks->UncheckedAt(mit1->second);
         if (mcTr->GetPdgCode() == fParams.pdgCodePr)
            nppi[0] = 1;
         else if (mcTr->GetPdgCode() == fParams.pdgCodeNeg)
            nppi[1] = 1;
         // cout << mcTr->GetPdgCode() << endl;
         mcTr->GetMomentum(mom);
         if (mom.Pt() < 0.001) continue;
         if (TMath::Abs(mom.Eta()) < 1.3) ++nok;
         if ((TMath::Abs(mom.Eta()) < 1.3) && mom.Pt() > 0.05) ++nok1;
         if ((TMath::Abs(mom.Eta()) < 1.3) && mom.Pt() > 0.1) ++nok2;
         if ((TMath::Abs(mom.Eta()) < 1.3) && mom.Pt() > 0.2) ++nok3;
      }
      if (nppi[0] != 1 || nppi[1] != 1) {
         // not p - p- decay
         // V cout << " Wrong decay mode !!! " << endl;
         mit = mapLamb.upper_bound(mothID);
         continue;
      }

      if (nppi[0] == 1 && nppi[1] == 1) fHistosF["hLambFlag"]->Fill(1);
      if (nok == 2) fHistosF["hLambFlag"]->Fill(2);
      if (nok1 == 2) fHistosF["hLambFlag"]->Fill(4);
      if (nok2 == 2) fHistosF["hLambFlag"]->Fill(6);
      if (nok3 == 2) fHistosF["hLambFlag"]->Fill(8);

      // Check Xi-
      MpdMCTrack *mcTr = (MpdMCTrack *)mcTracks->UncheckedAt(mothID);
      Int_t       gmID = mcTr->GetMotherId();

      if (mapXi.find(gmID) != mapXi.end()) {
         ret = mapXi.equal_range(gmID);

         for (mit1 = ret.first; mit1 != ret.second; ++mit1) {
            MpdMCTrack *mcTr = (MpdMCTrack *)mcTracks->UncheckedAt(mit1->second);
            if (mcTr->GetPdgCode() != fParams.pdgCodeNeg) continue;
            fHistosF["hXiFlag"]->Fill(1);
            mcTr->GetMomentum(mom);
            if (mom.Pt() < 0.001) continue;
            if (TMath::Abs(mom.Eta()) < 1.3 && nok == 2) fHistosF["hXiFlag"]->Fill(2);
            if (TMath::Abs(mom.Eta()) < 1.3 && mom.Pt() > 0.05 && nok1 == 2) fHistosF["hXiFlag"]->Fill(4);
            if (TMath::Abs(mom.Eta()) < 1.3 && mom.Pt() > 0.1 && nok2 == 2) fHistosF["hXiFlag"]->Fill(6);
            if (TMath::Abs(mom.Eta()) < 1.3 && mom.Pt() > 0.2 && nok3 == 2) fHistosF["hXiFlag"]->Fill(8);
         }
      }
      mit = mapLamb.upper_bound(mothID);
   } // while (mit != mapLamb.end())

   // Track selection

   for (Int_t j = 0; j < nITS; ++j) {
      AzTrack *tr = (AzTrack *)itsTracks->UncheckedAt(j);
      if (tr->GetNofTrHits() > 15 && tr->Pt() > 0.15 && TMath::Abs(tr->Momentum3().Eta()) < 0.5)
         ++fNtr13; // AZ-271022 - for centrality
      if (tr->GetChi2() < -8) continue;
      // if (tr->ClassName().Contains("Its") && tr->GetNofIts() > 0) continue;
      Int_t    id     = tr->GetTrackID();
      Double_t thRec  = tr->Theta();
      Double_t etaRec = tr->Momentum3().Eta();
      // AZ-171222 if (TMath::Abs(lays[id]) < -41 || TMath::Abs(etaRec) > 1.3) tr->SetChi2(-9.); // flag
      if (TMath::Abs(lays[id]) < -41 || TMath::Abs(etaRec) > 13) tr->SetChi2(-9.); // AZ-171222
      Int_t iQ = tr->Charge();
#ifdef ITS
      if (TString(tr->ClassName()).Contains("Its") && tr->GetNofHits() - tr->GetNofIts() < 10) tr->SetChi2(-9.);
#else
      if (tr->GetNofHits() < 10) tr->SetChi2(-9.);
#endif
      if (tr->GetChi2() < -8) continue;
      // Create MpdHelix
      /*
      MpdHelix helix = MakeHelix(tr);
      // Get 3-D DCA to primary vertex
      TVector3 pca;
      Double_t s = helix.pathLength(primVert);
      pca = helix.at(s);
      pca -= primVert;
      if (iQ < 0) {
   if (pdgs[id] != pdgCodeKm && pca.Mag() < gDCApi) tr->SetChi2(-9.);
      }
      else if (iQ > 0 && pca.Mag() < gDCAp) tr->SetChi2(-9.);
      */
      //++ntr13;
   } // for (Int_t j = 0; j < nITS;

   // Collect "good" pions, kaons and protons
   vector<Int_t> vecPi, vecK, vecP;

   for (Int_t j = 0; j < nITS; ++j) {
      MpdKalmanTrack *tr = (MpdKalmanTrack *)itsTracks->UncheckedAt(j);
      if (tr->GetChi2() < -8) continue;
      // if (TString(tr->ClassName()).Contains("Its") && tr->GetNofIts() > 0) continue;
      Int_t       id   = tr->GetTrackID();
      MpdMCTrack *mcTr = (MpdMCTrack *)mcTracks->UncheckedAt(id);
      // !!!
      if (mcTr->GetMotherId() == 0 && ((MpdMCTrack *)mcTracks->UncheckedAt(0))->GetPdgCode() == 1010010030)
         continue; // !!! decay product of artificial H3L
      // !!!
      //*MC ID
      // if (mcTr->GetPdgCode() == pdgCodePr && tr->Charge() == 1) vecP.push_back(j);
      // else if (mcTr->GetPdgCode() == pdgCodeNeg && tr->Charge() == -1) vecPi.push_back(j);
      // else if (mcTr->GetPdgCode() == pdgCodeKm && tr->Charge() == -1) vecK.push_back(j);
      if (mcTr->GetPdgCode() == fParams.pdgCodePr &&
          tr->Charge() == TMath::Nint(TDatabasePDG::Instance()->GetParticle(fParams.pdgCodePr)->Charge() / 3))
         vecP.push_back(j);
      else if (mcTr->GetPdgCode() == fParams.pdgCodeNeg &&
               tr->Charge() == TMath::Nint(TDatabasePDG::Instance()->GetParticle(fParams.pdgCodeNeg)->Charge() / 3))
         vecPi.push_back(j);
   }

   // V if (fEvNo % 100 == 0) cout << " Number of protons, pi: " << vecP.size() << " " << vecPi.size() << endl;
   RecoEff(vecP, vecPi, 1);
   // RecoEff(vecP, vecPi, 0);
   // V if (fEvNo % 100 == 0) cout << " Number of protons, pi: " << vecP.size() << " " << vecPi.size() << endl;

   // Apply PID
   // ApplyPid(newPid, vecP, vecPi);
   ApplyPid(vecP, vecPi);

   vector<MpdParticle *> vecL;
   vecL.clear();
   fvLambdas.clear();
   BuildLambda(vecP, vecPi, vecL);
   fvXis.clear();
   // if (vecL.size()) BuildCascade(vecK, vecPi, vecL);
   // vLambdas.clear(); // do not store lambda tree

   Int_t nLamb = vecL.size();
   for (Int_t ipart = 0; ipart < nLamb; ++ipart) delete vecL[ipart];
}

//-----------------------------------------------------------------------------------------------

void MpdAnaXiMix::RecoEff(vector<Int_t> &vecP, vector<Int_t> &vecPi, Int_t pid)
{
   // Check reco efficiency

   Int_t         nPi = vecPi.size(), nP = vecP.size();
   TClonesArray *itsTracks = fItsTracks, *mcTracks = fMcTracks;

   for (Int_t ip = nP - 1; ip >= 0; --ip) {
      // AntiProton
      AzTrack    *trP    = (AzTrack *)itsTracks->UncheckedAt(vecP[ip]);
      MpdMCTrack *mcTr   = (MpdMCTrack *)mcTracks->UncheckedAt(trP->GetTrackID());
      Int_t       mothId = mcTr->GetMotherId();
      if (mothId < 0) continue;
      MpdMCTrack *moth = (MpdMCTrack *)mcTracks->UncheckedAt(mothId);
      if (moth->GetPdgCode() == fParams.pdgCodeL0) {
         Int_t mp = mothId;
         // Proton from Lambda

         for (Int_t jpi = nPi - 1; jpi >= 0; --jpi) {
            // Pion
            AzTrack    *trPi   = (AzTrack *)itsTracks->UncheckedAt(vecPi[jpi]);
            MpdMCTrack *mcTr   = (MpdMCTrack *)mcTracks->UncheckedAt(trPi->GetTrackID());
            Int_t       mothId = mcTr->GetMotherId();
            if (mothId < 0) continue;
            MpdMCTrack *moth = (MpdMCTrack *)mcTracks->UncheckedAt(mothId);
            if (moth->GetPdgCode() == fParams.pdgCodeL0 && mp == mothId) {
               fHistosF["hLambFlag"]->Fill(12);
               if (TMath::Abs(trP->Momentum3().Eta()) < 1.3 && TMath::Abs(trPi->Momentum3().Eta()) < 1.3)
                  fHistosF["hLambFlag"]->Fill(10);
               // AZ - flag decay tracks to check PID influence later
               trP->SetUniqueID(mothId + 1);
               trPi->SetUniqueID(mothId + 1);
               //
               Int_t gmId = moth->GetMotherId();
               if (gmId >= 0) {
                  MpdMCTrack *gmoth = (MpdMCTrack *)mcTracks->UncheckedAt(gmId);
                  if (gmoth->GetPdgCode() == fParams.pdgCodeXi) {
                     for (Int_t kpi = nPi - 1; kpi >= 0; --kpi) {
                        // Pion
                        AzTrack    *trK    = (AzTrack *)itsTracks->UncheckedAt(vecPi[kpi]);
                        MpdMCTrack *mcTr   = (MpdMCTrack *)mcTracks->UncheckedAt(trK->GetTrackID());
                        Int_t       mothId = mcTr->GetMotherId();
                        if (mothId < 0) continue;
                        MpdMCTrack *moth = (MpdMCTrack *)mcTracks->UncheckedAt(mothId);
                        if (moth->GetPdgCode() == fParams.pdgCodeXi && gmId == mothId) {
                           fHistosF["hXiFlag"]->Fill(10);
                           // AZ - flag decay tracks to check PID influence later
                           trK->GetVertex().SetUniqueID(mothId + 1);
                           trP->GetVertex().SetUniqueID(mothId + 1);
                           trPi->GetVertex().SetUniqueID(mothId + 1);
                           //
                           break;
                        }
                     }
                  } // if (gmoth->GetPdgCode() == pdgCodeAXi)
               }
               break;
            }
         }
      }
   }

   if (pid) return; // skip the rest if PID is used

   for (Int_t ip = nP - 1; ip >= 0; --ip) {
      // Proton
      AzTrack *trP   = (AzTrack *)itsTracks->UncheckedAt(vecP[ip]);
      AzTrack  trCor = *trP;
      //*
      trCor.SetDirection(MpdKalmanTrack::kInward);
      if (fRecoIts)
         fRecoIts->Refit((MpdItsKalmanTrack *)&trCor, 0.93827, 1); // refit
      else
         fRecoTpc->Refit(&trCor, 0.93827, 1); // refit
      MpdParticle prot(trCor, vecP[ip]);
      prot.SetPdg(fParams.pdgCodePr);
      prot.SetMass();
      //*/
      Double_t chi2 = TMath::Min(prot.Chi2Vertex(fMpdVert), 999.);
      if (chi2 < fParams.gC2p) vecP.erase(vecP.begin() + ip);
   }

   if (nP) {
      for (Int_t jpi = nPi - 1; jpi >= 0; --jpi) {
         // Pion
         AzTrack *trPi  = (AzTrack *)itsTracks->UncheckedAt(vecPi[jpi]);
         AzTrack  trCor = *trPi;

         /*
         recoTpc->Refit(&trCor, 0.13957, 1); // refit
         MpdParticle pion(trCor, vecPi[jpi]);
         pion.SetPdg(pdgCodeNeg);
         pion.SetMass();
         */
         Double_t chi2 = TMath::Min(trPi->GetChi2Vertex(), 999.);
         // Double_t chi2 = TMath::Min (pion.Chi2Vertex(mpdVert),999.);
         if (chi2 < fParams.gC2pi) vecPi.erase(vecPi.begin() + jpi);
      }
   }
}

//-----------------------------------------------------------------------------------------------

void MpdAnaXiMix::ApplyPid(vector<Int_t> &vecP, vector<Int_t> &vecPi)
{
   // Apply PID

   // AZ - get information on hyperon decay products
   map<Int_t, set<Int_t>> mapL, mapXi, mapL13, mapXi13;
   TClonesArray          *itsTracks = fItsTracks;

   Int_t nP = vecP.size(), nPi = vecPi.size();

   for (Int_t ip = 0; ip < nP; ++ip) {
      AzTrack *trP = (AzTrack *)itsTracks->UncheckedAt(vecP[ip]);
      if (trP->GetUniqueID() > 0) {
         // Lambda decay product
         Int_t mid = trP->GetUniqueID();
         if (mapL.find(mid) == mapL.end()) {
            set<Int_t> aaa;
            mapL[mid] = aaa;
         }
         mapL[mid].insert(vecP[ip]);
         // AZ-120223 - |eta|<1.3
         if (TMath::Abs(trP->Momentum3().Eta()) < 1.3) {
            if (mapL13.find(mid) == mapL13.end()) {
               set<Int_t> aaa;
               mapL13[mid] = aaa;
            }
            mapL13[mid].insert(vecP[ip]);
         }
         trP->SetUniqueID(0); // reset
      }
      if (trP->GetVertex().GetUniqueID() > 0) {
         // Xi- decay product
         Int_t mid = trP->GetVertex().GetUniqueID();
         if (mapXi.find(mid) == mapXi.end()) {
            set<Int_t> aaa;
            mapXi[mid] = aaa;
         }
         mapXi[mid].insert(vecP[ip]);
      }
   }

   for (Int_t ip = 0; ip < nPi; ++ip) {
      AzTrack *trP = (AzTrack *)itsTracks->UncheckedAt(vecPi[ip]);
      if (trP->GetUniqueID() > 0) {
         // Lambda decay product
         Int_t mid = trP->GetUniqueID();
         if (mapL.find(mid) == mapL.end()) {
            set<Int_t> aaa;
            mapL[mid] = aaa;
         }
         mapL[mid].insert(vecPi[ip]);
         // AZ-120223 - |eta|<1.3
         if (TMath::Abs(trP->Momentum3().Eta()) < 1.3) {
            if (mapL13.find(mid) == mapL13.end()) {
               set<Int_t> aaa;
               mapL13[mid] = aaa;
            }
            mapL13[mid].insert(vecPi[ip]);
         }
         trP->SetUniqueID(0); // reset
      }
      if (trP->GetVertex().GetUniqueID() > 0) {
         // Xi- decay product
         Int_t mid = trP->GetVertex().GetUniqueID();
         if (mapXi.find(mid) == mapXi.end()) {
            set<Int_t> aaa;
            mapXi[mid] = aaa;
         }
         mapXi[mid].insert(vecPi[ip]);
      }
   }

   for (map<Int_t, set<Int_t>>::iterator mit = mapL.begin(); mit != mapL.end(); ++mit)
      if (mit->second.size() != 2) mit->second.insert(-999); // not 2 decay products reconstructed - add fake indx

   for (map<Int_t, set<Int_t>>::iterator mit = mapL13.begin(); mit != mapL13.end(); ++mit)
      if (mit->second.size() != 2) mit->second.insert(-999); // not 2 decay products reconstructed - add fake indx

   // Get TOF matches                                                           \

   Int_t             nTofMatch = fTofMatches->GetEntriesFast();
   map<Int_t, Int_t> mapTof;

   for (Int_t itof = 0; itof < nTofMatch; ++itof) {
      MpdTofMatchingData *match        = (MpdTofMatchingData *)fTofMatches->UncheckedAt(itof);
      mapTof[match->GetKFTrackIndex()] = itof;
   }

   vecP.clear();
   vecPi.clear();

   Int_t nITS = itsTracks->GetEntriesFast();

   for (Int_t j = 0; j < nITS; ++j) {
      AzTrack *tr = (AzTrack *)itsTracks->UncheckedAt(j);
      if (tr->GetChi2() < -8) continue;
      // if (TString(tr->ClassName()).Contains("Its") && tr->GetNofIts() > 0) continue;
      Int_t       id     = tr->GetTrackID();
      MpdMCTrack *mcTr   = (MpdMCTrack *)fMcTracks->UncheckedAt(id);
      Int_t       mothId = mcTr->GetMotherId();
      Int_t       uid    = tr->GetVertex().GetUniqueID();
      /*
      if (mothId == 0) {
        FairMCTrack* moth = (FairMCTrack*) mcTracks->UncheckedAt(mothId);
        if (moth->GetPdgCode() == pdgCodeH3L) continue; // track from artificial H3L
      }
      */

      MpdTrack *mpdTrack = (MpdTrack *)fMpdTracks->UncheckedAt(j);
      if (mpdTrack->GetID() != id) {
         // V cout << id << " " << mpdTrack->GetID() << endl;
         Fatal("ApplyPid", " Different ID");
      }
      // if (tr->GetNofTrHits() < 20) continue; // !!!

      Int_t    ret = 0, eta13 = 0, charge = tr->Charge(), tofFlag = mpdTrack->GetTofFlag();
      Double_t trEta = TMath::Abs(tr->Momentum3().Eta());
      // AZ Double_t dedx = tr->GetPartID(), m2 = mpdTrack->GetTofMass2();
      Double_t dedx = tr->GetDedx(), m2 = -1;
      if (mapTof.count(j) > 0) {
         MpdTofMatchingData *match = (MpdTofMatchingData *)fTofMatches->UncheckedAt(mapTof[j]);
         m2                        = match->GetMass2();
         tofFlag                   = 6;
      }
      if (TMath::Abs(tr->Momentum3().Eta()) < 1.3) eta13 = 1;

      Int_t asym = 0;
      if (tofFlag == 2 || tofFlag == 6) { // dE/dx+TOF
         // Apply asymmetric window
         // Double_t up = pid->Getm2TopBound(MpdPidUtils::kProton);
         // Double_t down = pid->Getm2BottomBound(MpdPidUtils::kProton);
         /*
         if (charge < 0 && tr->Momentum() > 1.0) {
      pid->Setm2BottomBound(1.5, MpdPidUtils::kProton); // 1 sigma for antiprotons above 1 GeV/c
      asym = 1;
         }
         */
         ret = fPid->FillProbs(tr->Momentum(), dedx, m2, charge);
         // pid->Setm2TopBound(up, MpdPidUtils::kProton);
         // pid->Setm2BottomBound(down, MpdPidUtils::kProton);
      }
      // if (tofFlag == 0 || tofFlag == 4)   // only dE/dx available
      if (ret == 0 && asym == 0) ret = fPid->FillProbs(tr->Momentum(), dedx, charge);

      //!!! No PID for (anti)protons above 2.5 GeV/c !!!
      // if (tr->Momentum() > 2.5 && charge < 0) ret = 1;
      if (tr->Momentum() > 2.5 && charge * fParams.pdgCodePr > 0) ret = 1; // AZ-130423

      if (ret == 0) {
         // No PID
         if (eta13) {
            if (mcTr->GetPdgCode() == fParams.pdgCodeNeg) fHistosF["hPIDflag"]->Fill(2.1); // lost pion
            if (mcTr->GetPdgCode() == fParams.pdgCodePr) fHistosF["hPIDflag"]->Fill(6.1);  // lost proton
         }
         continue;
      }

      Double_t piThr = -0.75, probThr = -0.60;

      if (fParams.pdgCodeL0 * tr->Charge() < 0) {
         Double_t prob = fPid->GetProbPi();
         if (prob > piThr && prob > fPid->GetProbKa() && prob > fPid->GetProbPr() && prob > fPid->GetProbDe() &&
             prob > fPid->GetProbTr() && prob > fPid->GetProbHe3() && prob > fPid->GetProbHe4()) {
            // "pion"
            if (mcTr->GetPdgCode() == fParams.pdgCodeNeg && eta13)
               fHistosF["hPIDflag"]->Fill(0.1); // correct pion
            else if (mcTr->GetPdgCode() != fParams.pdgCodeNeg && eta13)
               fHistosF["hPIDflag"]->Fill(1.1); // false pion
            //
            if (mapL.find(mothId + 1) != mapL.end() && mapL[mothId + 1].find(j) != mapL[mothId + 1].end())
               mapL[mothId + 1].erase(j);
            if (mapL13.find(mothId + 1) != mapL13.end() && mapL13[mothId + 1].find(j) != mapL13[mothId + 1].end())
               mapL13[mothId + 1].erase(j);
            if (mapXi.find(uid) != mapXi.end() && mapXi[uid].find(j) != mapXi[uid].end()) mapXi[uid].erase(j);
            //
            Double_t chi2 = TMath::Min(tr->GetChi2Vertex(), 999.);
            if (chi2 < fParams.gC2pi) continue;
            vecPi.push_back(j);
         } else if (mcTr->GetPdgCode() == fParams.pdgCodeNeg && eta13)
            fHistosF["hPIDflag"]->Fill(2.1); // lost pion
      } else {
         if (trEta < 0.5 && mcTr->GetPdgCode() == fParams.pdgCodePr) fHistosF["hPtProt"]->Fill(tr->Pt());
         Double_t prob = fPid->GetProbPr();
         // AZ-130423 if (tr->Momentum() > 2.5 && charge < 0) prob = 9.9; //!!! force antiproton !!!
         if (tr->Momentum() > 2.5 && charge * fParams.pdgCodePr > 0) prob = 9.9; //!!! force (anti)proton !!! AZ-130423

         if (prob > probThr && prob > fPid->GetProbKa() && prob > fPid->GetProbPi() && prob > fPid->GetProbDe()) {
            // "proton"

            if (mcTr->GetPdgCode() == fParams.pdgCodePr) {
               if (eta13) fHistosF["hPIDflag"]->Fill(4.1); // correct proton
               if (trEta < 0.5) fHistosF["hPtProtT"]->Fill(tr->Pt());
            } else if (mcTr->GetPdgCode() != fParams.pdgCodePr) {
               if (eta13) fHistosF["hPIDflag"]->Fill(5.1); // false proton
               if (trEta < 0.5) fHistosF["hPtProtF"]->Fill(tr->Pt());
            }
            //
            if (mapL.find(mothId + 1) != mapL.end() && mapL[mothId + 1].find(j) != mapL[mothId + 1].end())
               mapL[mothId + 1].erase(j);
            if (mapL13.find(mothId + 1) != mapL13.end() && mapL13[mothId + 1].find(j) != mapL13[mothId + 1].end())
               mapL13[mothId + 1].erase(j);
            if (mapXi.find(uid) != mapXi.end() && mapXi[uid].find(j) != mapXi[uid].end()) mapXi[uid].erase(j);
            //
            AzTrack trCor = *tr;
            trCor.SetDirection(MpdKalmanTrack::kInward);
            int ok = 0;
            if (fRecoIts)
               ok = fRecoIts->Refit((MpdItsKalmanTrack *)&trCor, 0.93827, 1); // refit
            else
               ok = fRecoTpc->Refit(&trCor, 0.93827, 1); // refit
            if (!ok) continue;
            MpdParticle prot(trCor, 0);
            prot.SetPdg(fParams.pdgCodePr);
            prot.SetMass();
            // Double_t chi2 = TMath::Min (trP->GetChi2Vertex(),999.);
            Double_t chi2 = TMath::Min(prot.Chi2Vertex(fMpdVert), 999.);
            if (chi2 < fParams.gC2p) continue;
            vecP.push_back(j);
         } else if (mcTr->GetPdgCode() == fParams.pdgCodePr && eta13)
            fHistosF["hPIDflag"]->Fill(6.1); // lost proton
      }
      // GetProbKa(), GetProbEl(), GetProbPr(), GetProbDe(), GetProbHe3
   } // for (Int_t j = 0; j < nITS;
   // V if (fEvNo % 100 == 0) cout << " Number of p, pi: " << vecP.size() << " " << vecPi.size() << endl;

   //
   Int_t nLok = 0, nXiok = 0, nLok13 = 0, nXiok13 = 0;

   for (map<Int_t, set<Int_t>>::iterator mit = mapL.begin(); mit != mapL.end(); ++mit) {
      if (mit->second.size() == 0) ++nLok;
      // else if (mit->second < 0) { cout << " Very strange 3 !!! " << mit->second << endl; exit(0); }
   }

   for (map<Int_t, set<Int_t>>::iterator mit = mapL13.begin(); mit != mapL13.end(); ++mit) {
      if (mit->second.size() == 0) ++nLok13;
      // else if (mit->second < 0) { cout << " Very strange 3 !!! " << mit->second << endl; exit(0); }
   }

   for (map<Int_t, set<Int_t>>::iterator mit = mapXi.begin(); mit != mapXi.end(); ++mit) {
      if (mit->second.size() == 0) ++nXiok;
      // if (mit->second < 0) { cout << " Very strange 4 !!! " << mit->second << endl; exit(0); }
   }
   fHistosF["hLambFlag"]->Fill(13, nLok);
   fHistosF["hLambFlag"]->Fill(11, nLok13);
   fHistosF["hXiFlag"]->Fill(11, nXiok);
   //
}

//-----------------------------------------------------------------------------------------------

void MpdAnaXiMix::BuildLambda(vector<Int_t> &vecP, vector<Int_t> &vecPi, vector<MpdParticle *> &vecL)
{
   // Make antilambdas

   int                   qs[2] = {0}, origs[2] = {0}, layMx[2] = {0}, dstNo[2] = {0};
   Float_t               etas[2] = {0}, ps[2] = {0}, pts[2] = {0}, chi2s[2] = {0}, dcas[2] = {0}, c2s[2] = {0};
   Float_t               path = 0, massh = 0, chi2h = 0, angle = 0, pth = 0, ph = 0, etah = 0, disth = 0, yh = 0;
   TClonesArray         *itsTracks = fItsTracks, *mcTracks = fMcTracks;
   Int_t                 nPi = vecPi.size(), nP = vecP.size(), saveMix = 0, mpdg = 0;
   vector<MpdParticle *> vPart;
   fVecL1.clear();
   fVecL2.clear();
   TVector3 primVert;
   fMpdVert->Position(primVert);

   for (Int_t ip = 0; ip < nP; ++ip) {
      // Proton
      AzTrack    *trP    = (AzTrack *)itsTracks->UncheckedAt(vecP[ip]);
      MpdMCTrack *mcTr   = (MpdMCTrack *)mcTracks->UncheckedAt(trP->GetTrackID());
      Int_t       mothId = mcTr->GetMotherId();
      AzTrack     trCor  = *trP;
      trCor.SetDirection(MpdKalmanTrack::kInward);
      if (fRecoIts)
         fRecoIts->Refit((MpdItsKalmanTrack *)&trCor, 0.93827, 1); // refit
      else
         fRecoTpc->Refit(&trCor, 0.93827, 1); // refit
      // MpdParticle prot(*trP, vecP[ip]);
      MpdParticle prot(trCor, vecP[ip]);
      prot.SetPdg(fParams.pdgCodePr);
      prot.SetMass();
      qs[1]   = TMath::Nint(TDatabasePDG::Instance()->GetParticle(fParams.pdgCodePr)->Charge() / 3);
      etas[1] = trP->Momentum3().Eta();
      ps[1]   = trP->Momentum();
      pts[1]  = trP->Pt();
      // chi2s[1] = TMath::Min (trP->GetChi2Vertex(),999.);
      chi2s[1]       = TMath::Min(prot.Chi2Vertex(fMpdVert), 999.);
      c2s[1]         = trP->GetChi2() / (trP->GetNofTrHits() * 2 - 5);
      layMx[1]       = TMath::Abs(fLays[trP->GetTrackID()]);
      MpdHelix helix = MakeHelix(trP);
      // Get 3-D DCA to primary vertex
      TVector3 pca;
      Double_t s = helix.pathLength(primVert);
      pca        = helix.at(s);
      pca -= primVert;
      dcas[1]  = pca.Mag();
      origs[1] = 0;
      if (mothId >= 0 && ((MpdMCTrack *)mcTracks->UncheckedAt(mothId))->GetPdgCode() == fParams.pdgCodeL0)
         origs[1] = -1; // from lambda
      ++saveMix;

      for (Int_t jpi = 0; jpi < nPi; ++jpi) {
         // Pion
         // 071122 MpdKalmanTrack *trPi = (MpdKalmanTrack*) itsTracks->UncheckedAt(vecPi[jpi]);
         AzTrack    *trPi    = (AzTrack *)itsTracks->UncheckedAt(vecPi[jpi]);
         MpdMCTrack *mcTr1   = (MpdMCTrack *)mcTracks->UncheckedAt(trPi->GetTrackID());
         Int_t       mothId1 = mcTr1->GetMotherId();
         origs[0]            = 0;
         if (mothId1 >= 0 && ((MpdMCTrack *)mcTracks->UncheckedAt(mothId1))->GetPdgCode() == fParams.pdgCodeL0)
            origs[0] = -1; // from lambda
         MpdParticle *pion = new MpdParticle(*trPi, vecPi[jpi]);
         pion->SetPdg(fParams.pdgCodeNeg);
         pion->SetMass();
         if (fParams.nMix > 0 && saveMix == 1) fMapPiEvent.insert(pair<int, AzTrack>(fEvNo, *trPi));

         vPart.clear();
         vPart.push_back(new MpdParticle(prot));
         // vPart.push_back(&prot);
         vPart.push_back(pion);

         MpdParticle lambPart;
         Double_t    chi2 = lambPart.BuildMother(vPart);
         TVector3    v0(lambPart.Getx()(0, 0), lambPart.Getx()(1, 0), lambPart.Getx()(2, 0));
         v0 -= primVert;
         Double_t decay = v0.Mag();
         path           = TMath::Sign(decay, v0 * lambPart.Momentum3());
         massh          = lambPart.GetMass();

         // AZ-061122 if (chi2 >= 0 && chi2 < gC2L && path > gPathL) {
         if (chi2 >= 0 && chi2 < fParams.gC2L && path > fParams.gPathL && massh < 1.2) { // AZ-061122
            if (origs[1] > 0) origs[1] = -1;
            MpdMCTrack *moth = NULL;
            fHistosF["hMassL"]->Fill(lambPart.GetMass());
            if (mothId != mothId1 || mothId < 0) {
               fHistosF["hMassLbkg"]->Fill(lambPart.GetMass());
            } else {
               // if (moth->GetPdgCode() == pdgCodeL0) {
               if (origs[0] == -1) {
                  fHistosF["hMassLsig"]->Fill(lambPart.GetMass());
                  origs[0] = origs[1] = 1;
                  moth                = (MpdMCTrack *)mcTracks->UncheckedAt(mothId);
               } else
                  fHistosF["hMassLbkg"]->Fill(lambPart.GetMass());
            }

            // Fill tree
            qs[0]   = TMath::Nint(TDatabasePDG::Instance()->GetParticle(fParams.pdgCodeNeg)->Charge() / 3); // pion
            etas[0] = trPi->Momentum3().Eta();
            ps[0]   = trPi->Momentum();
            pts[0]  = trPi->Pt();
            // chi2s[0] = TMath::Min (trPi->GetChi2Vertex(),999.);
            chi2s[0]        = TMath::Min(pion->Chi2Vertex(fMpdVert), 999.);
            c2s[0]          = trPi->GetChi2() / (trPi->GetNofTrHits() * 2 - 5);
            layMx[0]        = TMath::Abs(fLays[trPi->GetTrackID()]);
            MpdHelix helix1 = MakeHelix(trPi);
            // Get 3-D DCA to primary vertex
            s   = helix1.pathLength(primVert);
            pca = helix1.at(s);
            pca -= primVert;
            dcas[0] = pca.Mag();

            chi2h = chi2;
            angle = v0.Angle(lambPart.Momentum3());
            pth   = lambPart.Pt();       // reconstructed
            ph    = lambPart.Momentum(); // reconstructed
            if (pth > 0.001)
               etah = lambPart.Momentum3().Eta();
            else
               etah = TMath::Sign(100., lambPart.Momentum3().Z());
            pair<Double_t, Double_t> paths = helix.pathLengths(helix1);
            TVector3                 p1    = helix.at(paths.first);
            TVector3                 p2    = helix1.at(paths.second);
            p1 -= p2;
            disth = p1.Mag(); // closest distance between daughters

            dstNo[0] = vecPi[jpi]; // pion index
            dstNo[1] = vecP[ip];   // proton index

            // Mass cut
            // if (massLamb < 1.10713 || massLamb > 1.12423) return; // lambda mass +- 5*1.71 MeV
            // if (massLamb < 1.10695 || massLamb > 1.12525) return; // lambda mass +- 5*1.83 MeV after the selection
            // cut sigma=1.83 if (lambPart.GetMass() >= 1.10695 && lambPart.GetMass() <= 1.12525) {
            if (lambPart.GetMass() >= 1.10518 && lambPart.GetMass() <= 1.12668) { // lambda mass +- 5*2.15 MeV
               vecL.push_back(new MpdParticle(lambPart));
               fVecL1.push_back(pair<Double_t, Double_t>(disth, angle));
               fVecL2.push_back(pair<Double_t, Double_t>(chi2s[0], chi2s[1]));
               // cout << chi2 << " " << lambPart.GetMass() << " " << vecL.size() << " " << trPi->Charge() << " " <<
               // trP->Charge() << " " << trPi->GetTrackID() << " " << trP->GetTrackID() << endl;
            }

            lambPart.SetMass(1.11568); // set true mass
            yh   = lambPart.Rapidity();
            mpdg = -1;
            if (origs[0] == 1) {
               // True lambda
               // pth = moth->GetPt();
               // yh = moth->GetRapidity();
               // Check mother of lambda
               Int_t gMothId = moth->GetMotherId();
               if (gMothId >= 0) {
                  origs[0] = origs[1] = 2; // secondary lambda
                  mpdg                = ((MpdMCTrack *)mcTracks->UncheckedAt(gMothId))->GetPdgCode();
               }
            }
            //((TTree*)gROOT->FindObjectAny("hypers"))->Fill();

            MpdLambda l0(massh, pth, ph, etah, yh, chi2h, disth, path, angle, etas, ps, pts, chi2s, dcas, c2s, origs,
                         qs, layMx, mpdg);
            fvLambdas.push_back(l0);
         } // if (chi2 >= 0 && chi2 < gC2L

         Int_t nPart = vPart.size();
         for (Int_t ipart = 0; ipart < nPart; ++ipart) delete vPart[ipart];
      } // for (Int_t jpi = 0; jpi < nPi;

      if (fParams.nMix == 0) continue;

      // Event mixing

      for (map<int, MpdVertex>::iterator mit = fMapVertexEvent.begin(); mit != fMapVertexEvent.end(); ++mit) {
         if (mit->first == fEvNo) break;
         Double_t dz = mit->second.GetZ() - fMapVertexEvent.rbegin()->second.GetZ();
         pair<multimap<Int_t, AzTrack>::iterator, multimap<Int_t, AzTrack>::iterator> piter =
            fMapPiEvent.equal_range(mit->first);

         for (multimap<Int_t, AzTrack>::iterator mit1 = piter.first; mit1 != piter.second; ++mit1) {
            AzTrack  piTr = mit1->second;
            Double_t z0tr = piTr.GetParam(1);
            piTr.SetParam(1, z0tr - dz);

            origs[0]          = -9;
            MpdParticle *pion = new MpdParticle(piTr);
            pion->SetPdg(fParams.pdgCodeNeg);
            pion->SetMass();

            vPart.clear();
            vPart.push_back(new MpdParticle(prot));
            vPart.push_back(pion);

            MpdParticle lambPart;
            Double_t    chi2 = lambPart.BuildMother(vPart);
            TVector3    v0(lambPart.Getx()(0, 0), lambPart.Getx()(1, 0), lambPart.Getx()(2, 0));
            v0 -= primVert;
            Double_t decay = v0.Mag();
            path           = TMath::Sign(decay, v0 * lambPart.Momentum3());
            massh          = lambPart.GetMass();

            if (chi2 >= 0 && chi2 < fParams.gC2L && path > fParams.gPathL && massh < 1.2) {
               if (origs[1] > 0) origs[1] = -1;
               /*
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
               */
               // Fill tree
               qs[0]   = TMath::Nint(TDatabasePDG::Instance()->GetParticle(fParams.pdgCodeNeg)->Charge() / 3); // pion
               etas[0] = piTr.Momentum3().Eta();
               ps[0]   = piTr.Momentum();
               pts[0]  = piTr.Pt();
               // chi2s[0] = TMath::Min (trPi->GetChi2Vertex(),999.);
               chi2s[0] = TMath::Min(pion->Chi2Vertex(fMpdVert), 999.);
               c2s[0]   = piTr.GetChi2() / (piTr.GetNofTrHits() * 2 - 5);
               // layMx[0] = TMath::Abs (lays[trPi->GetTrackID()]);
               MpdHelix helix1 = MakeHelix(&piTr);
               // Get 3-D DCA to primary vertex
               s   = helix1.pathLength(primVert);
               pca = helix1.at(s);
               pca -= primVert;
               dcas[0] = pca.Mag();

               chi2h = chi2;
               angle = v0.Angle(lambPart.Momentum3());
               pth   = lambPart.Pt();       // reconstructed
               ph    = lambPart.Momentum(); // reconstructed
               if (pth > 0.001)
                  etah = lambPart.Momentum3().Eta();
               else
                  etah = TMath::Sign(100., lambPart.Momentum3().Z());
               lambPart.SetMass(1.11568); // set true mass
               yh                             = lambPart.Rapidity();
               pair<Double_t, Double_t> paths = helix.pathLengths(helix1);
               TVector3                 p1    = helix.at(paths.first);
               TVector3                 p2    = helix1.at(paths.second);
               p1 -= p2;
               disth = p1.Mag(); // closest distance between daughters
               mpdg  = -9;       // mixing

               MpdLambda l0(massh, pth, ph, etah, yh, chi2h, disth, path, angle, etas, ps, pts, chi2s, dcas, c2s, origs,
                            qs, layMx, mpdg);
               fvLambdas.push_back(l0);
            } // if (chi2 >= 0 && chi2 < gC2L

            Int_t nPart = vPart.size();
            for (Int_t ipart = 0; ipart < nPart; ++ipart) delete vPart[ipart];
         } // for (multimap<Int_t,AzTrack>::iterator mit1 = piter.first;
      }
   } // for (Int_t ip = 0; ip < nP;
}

//-----------------------------------------------------------------------------------------------

MpdHelix MpdAnaXiMix::MakeHelix(const MpdKalmanTrack *tr)
{
   Double_t r   = tr->GetPosNew();
   Double_t phi = tr->GetParam(0) / r;
   Double_t x   = r * TMath::Cos(phi);
   Double_t y   = r * TMath::Sin(phi);
   Double_t dip = tr->GetParam(3);
   Double_t cur = 0.3 * 0.01 * 5 / 10; // 5 kG
   cur *= TMath::Abs(tr->GetParam(4));
   TVector3 o(x, y, tr->GetParam(1));
   Int_t    h = (Int_t)TMath::Sign(1.1, tr->GetParam(4));
   MpdHelix helix(cur, dip, tr->GetParam(2) - TMath::PiOver2() * h, o, h);
   return helix;
}

//-----------------------------------------------------------------------------------------------

MpdHelix MpdAnaXiMix::MakeHelix(const MpdParticle *part)
{
   Double_t dip = TMath::PiOver2() - part->Theta();
   Double_t cur = TMath::Abs(part->GetMeas(4));
   if (part->GetCharge() == 0) cur = numeric_limits<double>::epsilon();
   Int_t    h     = (Int_t)TMath::Sign(1.1, part->GetMeas(4));
   Double_t phase = part->GetMeas(2) - TMath::PiOver2() * h;
   Double_t x     = part->GetXY(0);
   Double_t y     = part->GetXY(1);
   TVector3 o(x, y, part->GetMeas(1));
   MpdHelix helix(cur, dip, phase, o, h);
   return helix;
}

//-----------------------------------------------------------------------------------------------
//*
void MpdAnaXiMix::BuildCascade(vector<Int_t> &vecK, vector<Int_t> &vecPi, vector<MpdParticle *> &vecL)
{
   // Make cascades (Xi- and Omega-)

   Int_t                 nPi = vecPi.size(), nK = vecK.size(), nL = vecL.size();
   vector<MpdParticle *> vPart, vPart1;
   TClonesArray         *itsTracks = fItsTracks, *mcTracks = fMcTracks;
   int                   qs[2] = {0}, origs[2] = {0}, layMx[2] = {0}, mpdg = 0;
   Float_t               etas[2] = {0}, ps[2] = {0}, pts[2] = {0}, chi2s[2] = {0}, dcas[2] = {0}, chi2sL[2] = {0};
   Float_t               path = 0, massh = 0, chi2h = 0, angle = 0, pth = 0, ph = 0, etah = 0, disth = 0, yh = 0;
   Float_t               masshL = 0, disthL = 0, angL = 0, chi2hL = 0, pathL = 0, c2pv = 0, d2pv = 0;
   Float_t               chi2hL1 = 0, masshL1 = 0, c2s[2] = {0};
   TVector3              primVert;
   fMpdVert->Position(primVert);

   // cout << " Cascade: " << nL << endl;
   for (Int_t iL = 0; iL < nL; ++iL) {
      // Lambda
      MpdParticle *lamb = vecL[iL];
      // cout << " Lambda No: " << iL << " " << nPi << " " << nK << " " << lamb->GetMass() << " "
      //  << lamb->GetCharge() << endl;
      //  Check whether lambda is a signal or background
      Double_t     lambTrue = kFALSE;
      MpdMCTrack  *moth0    = NULL;
      MpdParticle *proton   = NULL;

      for (Int_t j = 0; j < 2; ++j) {
         AzTrack    *tr   = (AzTrack *)itsTracks->UncheckedAt(lamb->DaughterInds()[j]);
         MpdMCTrack *mcTr = (MpdMCTrack *)mcTracks->UncheckedAt(tr->GetTrackID());
         if (j == 0) {
            // Proton
            AzTrack trCor = *tr;
            trCor.SetDirection(MpdKalmanTrack::kInward);
            if (fRecoIts)
               fRecoIts->Refit((MpdItsKalmanTrack *)&trCor, 0.93827, 1); // refit
            else
               fRecoTpc->Refit(&trCor, 0.93827, 1); // refit
            proton = new MpdParticle(trCor, 0);
            proton->SetPdg(fParams.pdgCodePr);
            proton->SetMass();
            c2s[2 - j] = tr->GetChi2() / (tr->GetNofTrHits() * 2 - 5);
         }
         if (mcTr->GetMotherId() < 0) break;
         MpdMCTrack *moth = (MpdMCTrack *)mcTracks->UncheckedAt(mcTr->GetMotherId());
         if (moth->GetPdgCode() != fParams.pdgCodeL0) break;
         if (moth == moth0) lambTrue = kTRUE;
         moth0 = moth;
      }
      Int_t gMothId = -1;
      if (lambTrue) gMothId = moth0->GetMotherId();
      // Put lambda info here
      qs[1]   = 0;
      etas[1] = lamb->Momentum3().Eta();
      ps[1]   = lamb->Momentum();
      pts[1]  = lamb->Pt();
      Double_t xyzv[3];
      primVert.GetXYZ(xyzv);
      chi2s[1] = TMath::Min(lamb->Chi2Vertex(fMpdVert), 999.);
      if (lambTrue) {
         if (gMothId >= 0 && ((MpdMCTrack *)mcTracks->UncheckedAt(gMothId))->GetPdgCode() == fParams.pdgCodeXi)
            origs[1] = -1;
         else
            origs[1] = 0;
      } else
         origs[1] = -2;
      masshL = lamb->GetMass();
      disthL = fVecL1[iL].first;
      angL   = fVecL1[iL].second;
      chi2hL = lamb->Chi2();
      TVector3 v00(lamb->Getx()(0, 0), lamb->Getx()(1, 0), lamb->Getx()(2, 0));
      v00 -= primVert;
      pathL     = TMath::Sign(v00.Mag(), v00 * lamb->Momentum3());
      chi2sL[0] = fVecL2[iL].first;
      chi2sL[1] = fVecL2[iL].second;
      // Get 3-D DCA to primary vertex
      MpdHelix helix = MakeHelix(lamb);
      Double_t s     = helix.pathLength(primVert);
      TVector3 pca   = helix.at(s);
      pca -= primVert;
      dcas[1] = pca.Mag();

      for (Int_t ik = 0; ik < nK; ++ik) {
         // Kaon
         AzTrack     *trK  = (AzTrack *)itsTracks->UncheckedAt(vecK[ik]);
         MpdMCTrack  *mcTr = (MpdMCTrack *)mcTracks->UncheckedAt(trK->GetTrackID());
         MpdParticle *kaon = new MpdParticle(*trK, vecK[ik]);
         kaon->SetPdg(fParams.pdgCodeKm);
         kaon->SetMass();

         vPart.clear();
         vPart.push_back(new MpdParticle(*lamb));
         // vPart.push_back(lamb);
         vPart.push_back(kaon);

         MpdParticle omegaPart;
         Double_t    chi2 = TMath::Abs(omegaPart.BuildMother(vPart));
         TVector3    v0(omegaPart.Getx()(0, 0), omegaPart.Getx()(1, 0), omegaPart.Getx()(2, 0));
         v0 -= primVert;
         Double_t decay = v0.Mag();

         if (chi2 < fParams.gC2Xi && decay > fParams.gDecayOm) {
            fHistosF["hMassOm"]->Fill(omegaPart.GetMass());
            MpdMCTrack *moth = (MpdMCTrack *)mcTracks->UncheckedAt(mcTr->GetMotherId());
            if (lambTrue && gMothId >= 0 && mcTr->GetMotherId() == gMothId && moth->GetPdgCode() == fParams.pdgCodeOm) {
               fHistosF["hMassOmSig"]->Fill(omegaPart.GetMass());
            } else {
               fHistosF["hMassOmBkg"]->Fill(omegaPart.GetMass());
            }
         }

         delete vPart.front(); // delete lambda
         delete vPart.back();  // delete kaon
      }                        // for (Int_t ik = 0; ik < nK;

      for (Int_t jpi = 0; jpi < nPi; ++jpi) {
         // Pion
         if (vecPi[jpi] == lamb->DaughterInds()[1]) continue; // the same pion as was used for lambda

         AzTrack     *trPi = (AzTrack *)itsTracks->UncheckedAt(vecPi[jpi]);
         MpdParticle *pion = new MpdParticle(*trPi);
         MpdMCTrack  *mcTr = (MpdMCTrack *)mcTracks->UncheckedAt(trPi->GetTrackID());
         pion->SetPdg(fParams.pdgCodePos);
         pion->SetMass();

         vPart.clear();
         vPart.push_back(new MpdParticle(*lamb));
         // vPart.push_back(lamb);

         // !!! True lambda mass - it works
         vPart.back()->SetPdg(fParams.pdgCodeL0);
         vPart.back()->SetMass();

         vPart.push_back(new MpdParticle(*pion));

         MpdParticle xiPart;
         Double_t    chi2 = xiPart.BuildMother(vPart);
         TVector3    v0(xiPart.Getx()(0, 0), xiPart.Getx()(1, 0), xiPart.Getx()(2, 0));
         v0 -= primVert;
         Double_t decay = v0.Mag();
         massh          = xiPart.GetMass();

         // AZ-061122 if (chi2 >= 0 && chi2 < gC2Xi && decay > gDecayOm) {
         if (chi2 >= 0 && chi2 < fParams.gC2Xi && decay > fParams.gDecayOm && massh < 1.45) { // AZ-061122
            fHistosF["hMassXi"]->Fill(xiPart.GetMass());
            MpdMCTrack *moth = NULL;
            Int_t       mID  = mcTr->GetMotherId();
            if (mID >= 0) moth = (MpdMCTrack *)mcTracks->UncheckedAt(mID);
            if (lambTrue && gMothId >= 0 && mID == gMothId && moth->GetPdgCode() == fParams.pdgCodeXi) {
               fHistosF["hMassXiSig"]->Fill(xiPart.GetMass());
               origs[0] = origs[1] = 1;
            } else {
               fHistosF["hMassXiBkg"]->Fill(xiPart.GetMass());
               if (moth && moth->GetPdgCode() == fParams.pdgCodeXi)
                  origs[0] = -1;
               else
                  origs[0] = 0;
               if (origs[1] == 1) origs[1] = -1;
            }
            // Fill tree
            c2pv    = TMath::Min(xiPart.Chi2Vertex(fMpdVert), 999.);
            qs[0]   = TMath::Nint(TDatabasePDG::Instance()->GetParticle(fParams.pdgCodeNeg)->Charge() / 3); // pion
            etas[0] = trPi->Momentum3().Eta();
            ps[0]   = trPi->Momentum();
            pts[0]  = trPi->Pt();
            // chi2s[0] = TMath::Min (trPi->GetChi2Vertex(),999.);
            chi2s[0] = TMath::Min(pion->Chi2Vertex(fMpdVert), 999.);
            layMx[0] = TMath::Abs(fLays[trPi->GetTrackID()]);
            c2s[0]   = trPi->GetChi2() / (trPi->GetNofTrHits() * 2 - 5);
            // Get 3-D DCA to primary vertex
            MpdHelix helix1 = MakeHelix(trPi);
            Double_t s1     = helix1.pathLength(primVert);
            TVector3 pca1   = helix1.at(s1);
            pca1 -= primVert;
            dcas[0] = pca1.Mag();

            chi2h = chi2;
            path  = TMath::Sign(decay, v0 * xiPart.Momentum3());
            angle = v0.Angle(xiPart.Momentum3());
            pth   = xiPart.Pt(); // reconstructed
            ph    = xiPart.Momentum();
            if (pth > 0.001)
               etah = xiPart.Momentum3().Eta();
            else
               etah = TMath::Sign(100., xiPart.Momentum3().Z());

            disth = DistHelLin(trPi, lamb);
            // Get 3-D DCA to primary vertex
            helix1 = MakeHelix(&xiPart);
            s1     = helix1.pathLength(primVert);
            pca1   = helix1.at(s1);
            pca1 -= primVert;
            d2pv = pca1.Mag();

            xiPart.SetMass(1.3213); // set true mass
            yh   = xiPart.Rapidity();
            mpdg = -1;
            if (origs[0] == 1) {
               // True Xi
               // pth = moth->GetPt();
               // yh = moth->GetRapidity();
               // Check mother of Xi
               Int_t gMothId = moth->GetMotherId();
               if (gMothId >= 0) {
                  // origs[0] = origs[1] = 2; // secondary Xi
                  mpdg = ((MpdMCTrack *)mcTracks->UncheckedAt(gMothId))->GetPdgCode();
               }
            }
            // "Virtual lambda"
            vPart1.clear();
            vPart1.push_back(new MpdParticle(*proton));
            vPart1.push_back(pion);
            MpdParticle lambPart;
            chi2hL1 = lambPart.BuildMother(vPart1);
            masshL1 = lambPart.GetMass();
            delete vPart1.front();
            delete vPart1.back();

            //((TTree*)gROOT->FindObjectAny("Xi"))->Fill();

            MpdXi xi(massh, pth, ph, etah, yh, chi2h, disth, path, angle, c2pv, d2pv, etas, ps, pts, chi2s, dcas, c2s,
                     masshL, chi2hL, disthL, pathL, angL, chi2sL, masshL1, origs, qs, layMx, mpdg);
            fvXis.push_back(xi);
         } else
            delete pion;

         delete vPart.front(); // delete lambda
         delete vPart.back();  // delete pion
      }                        // for (Int_t jpi = 0; jpi < nPi;
      if (proton) {
         delete proton;
         proton = NULL;
      }

   } // for (Int_t iL = 0; iL < nL;
}
//*/
//-----------------------------------------------------------------------------------------------

Double_t MpdAnaXiMix::DistHelLin(MpdKalmanTrack *helix, MpdParticle *neu)
{
   // Compute distance between helix and straight line

   // Create MpdHelix
   *fTrC = MakeHelix(helix);
   fVtxN.SetXYZ(neu->Getx()(0, 0), neu->Getx()(1, 0), neu->Getx()(2, 0));
   fMomN = neu->Momentum3();
   fMomN *= (1. / fMomN.Mag());

   // gMinuit->SetFCN(fcn);
   fMinuit->SetFCN(Fcn);

   Double_t arglist[10];
   Int_t    ierflg = 0;
   arglist[0]      = 1;
   fMinuit->mnexcm("SET ERR", arglist, 1, ierflg);
   arglist[0] = -1; // 1; //-1;
   fMinuit->mnexcm("SET PRINT", arglist, 1, ierflg);
   fMinuit->mnexcm("SET NOWarnings", 0, 0, ierflg);

   Double_t        vstart[2] = {-0.1, -0.1};
   static Double_t step[2]   = {0.1, 0.1};
   fMinuit->mnparm(0, "lengN", vstart[0], step[0], 0, 0, ierflg);
   fMinuit->mnparm(1, "lengC", vstart[1], step[1], 0, 0, ierflg);

   // Now ready for minimization step
   arglist[0] = 500;
   arglist[1] = 1.;
   fMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);

   // Get results
   Double_t amin, edm, errdef;
   Int_t    nvpar, nparx, icstat;
   fMinuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);
   return amin;
}

//-----------------------------------------------------------------------------------------------

void MpdAnaXiMix::Fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
   // Compute distance between straight line and helix

   TVector3 mom = fMomN;
   mom *= par[0];
   TVector3 posN = fVtxN;
   posN += mom;

   TVector3 posC = fTrC->at(par[1]);
   posC -= posN;
   f = posC.Mag();
   // cout << par[0] << " " << par[1] << " " << f << endl;
}

//-----------------------------------------------------------------------------------------------

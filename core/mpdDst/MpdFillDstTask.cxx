// Author: Oleg Rogachevsky
// Update: 2009-10-07 17:56:17+0400
// Copyright: 2009 (C) MPD coll.
//
// fill dst task

#include "MpdFillDstTask.h"
#include "MpdKalmanTrack.h"
#include "MpdTpcKalmanTrack.h"
#include "MpdEctKalmanTrack.h"
#include "MpdTofMatchingData.h"
#include "MpdVertex.h"
#include "MpdParticleIdentification.h"
#include "MpdHelix.h"

#include "MpdFieldCreator.h"
#include "MpdMapPar.h"
#include "MpdMCTrack.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"

#include <TMath.h>
#include <TMatrixD.h>
#include "TGeoManager.h"
#include "TGeoBBox.h"
#include "TGeoTube.h"

#include <assert.h>
#include <set>
#include <iostream>
#include <algorithm>
#include <vector>
#include <bitset>

using namespace std;

// using std::cout;
// using std::flush;
// using std::endl;

// -----   constructor with names ------------------------------------

MpdFillDstTask::MpdFillDstTask(BaseQA *qaObject, const char *name, const char *title)
   : FairTask(name), fEvent(nullptr), fKFTracks(nullptr), fKFEctTracks(nullptr), fMCTracks(nullptr),
     fGenTracks(nullptr), fTpcHits(nullptr), fMCEventHeader(nullptr), fTofMatching(nullptr), fEtofMatching(nullptr),
     fVertex(nullptr), fZdcSkeletonesSaved(kFALSE), fHistZdc1En(nullptr), fHistZdc2En(nullptr),
     fELossZdc1Histo(nullptr), fELossZdc2Histo(nullptr), fELossZdc1Value(nullptr), fELossZdc2Value(nullptr),
     fhTrackMotherId(nullptr), fhTrackPrimaryPDG(nullptr), fhTrackVertex(nullptr), fhTruthVertex(nullptr),
     fPID(nullptr), fSharedHitArraySize(0), fSharedHitArray(nullptr)
{
   localQAptr = qaObject;
   if (localQAptr) LOG(info) << "QA mode on: BaseQA";
}
// -----   Destructor -----------------------------------------------

MpdFillDstTask::~MpdFillDstTask()
{
   if (fEvent) delete fEvent;
   if (fSharedHitArraySize > 0) delete[] fSharedHitArray;
}
// -------------------------------------------------------------------

InitStatus MpdFillDstTask::Init()
{
   fEvent = new MpdEvent();

   FairRootManager *manager = FairRootManager::Instance();

   fKFTracks      = (TClonesArray *)manager->GetObject("TpcKalmanTrack");
   fKFEctTracks   = (TClonesArray *)manager->GetObject("EctTrack");
   fMCTracks      = (TClonesArray *)manager->GetObject("MCTrack");
   fGenTracks     = (TClonesArray *)manager->GetObject("GenTracks");
   fTpcHits       = (TClonesArray *)manager->GetObject("TpcRecPoint");
   fMCEventHeader = (FairMCEventHeader *)manager->GetObject("MCEventHeader.");
   fTofMatching   = (TClonesArray *)manager->GetObject("TOFMatching");
   fEtofMatching  = (TClonesArray *)manager->GetObject("ETOFMatching");
   fVertex        = (TClonesArray *)manager->GetObject("Vertex"); // AZ

   fELossZdc1Value     = (TClonesArray *)manager->GetObject("ELossZdc1Value"); // EL
   fELossZdc2Value     = (TClonesArray *)manager->GetObject("ELossZdc2Value"); // EL
   fELossZdc1Histo     = (TClonesArray *)manager->GetObject("ELossZdc1Histo"); // EL
   fELossZdc2Histo     = (TClonesArray *)manager->GetObject("ELossZdc2Histo"); // EL
   fHistZdc1En         = (TH2F *)manager->GetObject("HistZdc1En");             // EL
   fHistZdc2En         = (TH2F *)manager->GetObject("HistZdc2En");             // EL
   fZdcSkeletonesSaved = 0;
   if (fPID == nullptr) {
      fPID = new MpdPid(0, 0, 0, 0, "DEFAULT", "CF", "");
      // PID with cluster finder
      // 	fPID->Init("DEFAULT","CF","");
   }
   if (fTpcHits == nullptr) {
      cout << "WARNING: no TpcRec array, can't fill hit maps!" << endl;
   }
   fSharedHitArraySize = 1000;
   fSharedHitArray     = new Short_t[fSharedHitArraySize];

   FairRootManager::Instance()->Register("MCEventHeader.", "MC", fMCEventHeader, kTRUE);
   FairRootManager::Instance()->Register("MCTrack", "MC", fMCTracks, kTRUE);
   if (fGenTracks) FairRootManager::Instance()->Register("GenTracks", "MC", fGenTracks, kTRUE);
   // FairRootManager::Instance()->Register("MpdEvent","MpdEvents", fEvents, kTRUE);
   FairRootManager::Instance()->Register("MPDEvent.", "MpdEvent", fEvent, kTRUE);
   // FairRootManager::Instance()->Register("EZdc1","EZdc1", fELossZdc1Histo, kTRUE);
   // FairRootManager::Instance()->Register("EZdc2","EZdc2", fELossZdc2Histo, kTRUE);

   // fhTrackMotherId = new TH1F("TrackMotherId","",1000,-5,995);
   // fhTrackPrimaryPDG = new TH1F("TrackPrimaryPDG","",1000,-500,500);
   // fhTrackVertex = new TH2F("TrackVertex","",1000,-500,500,1000,-500,500);
   // fhTruthVertex = new TH2F("TruthVertex","",1000,-500,500,1000,-500,500);

   return kSUCCESS;
}
// -------------------------------------------------------------------

void MpdFillDstTask::Exec(Option_t *option)
{
   Reset(); // Clear previos event information

   if (!fZdcSkeletonesSaved) { // empty skeletones saved only once
      FairRootManager *ioman = FairRootManager::Instance();
      if (fHistZdc1En) {
         fHistZdc1En->SetDirectory((TFile *)ioman->GetOutFile());
         fHistZdc2En->SetDirectory((TFile *)ioman->GetOutFile());
         fHistZdc1En->Write();
         fHistZdc2En->Write();
         fZdcSkeletonesSaved = 1;
      }
   }

   Int_t nReco    = fKFTracks ? fKFTracks->GetEntriesFast() : 0;
   Int_t nEctReco = fKFEctTracks ? fKFEctTracks->GetEntriesFast() : 0;
   cout << "\n-I- [MpdFillDstTask::Exec] " << nReco + nEctReco << " reconstruced tracks to write" << endl;

   if (fTpcHits) CalculateSharedArrayMap();

   FairRunAna *fRun = FairRunAna::Instance();
   fEvent->SetRunInfoRunId(fRun->GetRunId());

   FairField *fPar = fRun->GetField();
   fEvent->SetRunInfoMagneticFieldZ(fPar->GetType());

   if (fVertex) {
      MpdVertex *vertex = (MpdVertex *)fVertex->UncheckedAt(0);
      fEvent->SetPrimaryVerticesX(vertex->GetX());
      fEvent->SetPrimaryVerticesY(vertex->GetY());
      fEvent->SetPrimaryVerticesZ(vertex->GetZ());
   } else {
      // Vertex info is not available
      fEvent->SetPrimaryVerticesX(0.);
      fEvent->SetPrimaryVerticesY(0.);
      fEvent->SetPrimaryVerticesZ(0.);
   }
   // helix works with meters, primary vertices are in cm
   TVector3 recoVertex(fEvent->GetPrimaryVerticesX(), fEvent->GetPrimaryVerticesY(), fEvent->GetPrimaryVerticesZ());
   TVector3 mcVertex;
   fMCEventHeader->GetVertex(mcVertex);
   // check clone track into ECT
   typedef std::set<Int_t> trackSet;
   trackSet                EctTrackSet;
   for (Int_t index = 0; index < nEctReco; index++) // cycle by ECT KF tracks
   {
      MpdEctKalmanTrack *track = (MpdEctKalmanTrack *)fKFEctTracks->UncheckedAt(index);
      if (track->IsFromTpc()) EctTrackSet.insert(trackSet::value_type(track->GetTpcIndex()));
   }

   Int_t               nMatching = fTofMatching ? fTofMatching->GetEntriesFast() : 0;
   MpdTofMatchingData *pMatchingData;
   bool                matchingDataExist;

   // AZ MpdParticleIdentification *identificator = new MpdParticleIdentification();
   MpdParticleIdentification identificator;

   for (Int_t i = 0; i < nReco; i++) {
      // check clone track into ECT
      if (nEctReco > 0) {
         trackSet::iterator it = EctTrackSet.find(i);
         if (it != EctTrackSet.end()) continue;
      }

      MpdKalmanTrack    *kftrack    = (MpdKalmanTrack *)fKFTracks->UncheckedAt(i);
      MpdTpcKalmanTrack *kfTPCtrack = (MpdTpcKalmanTrack *)fKFTracks->UncheckedAt(i);

      MpdTrack *track = fEvent->AddGlobalTrack();
      if (kfTPCtrack->GetRecoQuality()) track->SetEdgeCut(kTRUE);
      track->SetID(kftrack->GetTrackID());
      track->SetNofHits(kftrack->GetNofHits());
      track->SetdEdXTPC(kftrack->GetPartID());
      Float_t Ppi, Pk, Pe, Pp;

      if (!identificator.GetTpcProbs(kftrack->Momentum3().Mag(), kftrack->GetPartID(), kftrack->GetNofHits(), Ppi, Pk,
                                     Pp, Pe, 0)) { // 0 - equal bayesian coefficients
         track->SetTPCpidProb(Pe, Ppi, Pk, Pp, BIT(2));
      }
      matchingDataExist = false;
      for (Int_t tofIndex = 0; tofIndex < nMatching; tofIndex++) {
         pMatchingData = (MpdTofMatchingData *)fTofMatching->UncheckedAt(tofIndex);
         if (pMatchingData->GetKFTrackIndex() == i) {
            matchingDataExist = true;
            break;
         } // first matching
      }

      if (matchingDataExist) {
         track->SetTofBeta(pMatchingData->GetBeta());
         track->SetTofMass2(pMatchingData->GetMass2());
         track->SetTofHitIndex(pMatchingData->GetTofHitIndex());

         if (!identificator.GetTofProbs(pMatchingData->GetMomentum().Mag(), pMatchingData->GetBeta(), Ppi, Pk, Pp, Pe,
                                        0)) {
            track->SetTOFpidProb(Pe, Ppi, Pk, Pp, BIT(1));
         }
      }
      Float_t tpcProbs[4] = {track->GetTPCPidProbPion(), track->GetTPCPidProbKaon(), track->GetTPCPidProbProton(),
                             track->GetTPCPidProbElectron()};
      Float_t tofProbs[4] = {track->GetTOFPidProbPion(), track->GetTOFPidProbKaon(), track->GetTOFPidProbProton(),
                             track->GetTOFPidProbElectron()};
      Float_t combProbs[4]; // probabilities combined from TOF & TPC
      identificator.GetCombinedProbs(tofProbs, tpcProbs, combProbs, 4);
      Ppi = combProbs[0];
      Pk  = combProbs[1];
      Pp  = combProbs[2];
      Pe  = combProbs[3];
      track->SetCombPidProb(Pe, Ppi, Pk, Pp);

      if (kftrack->GetParam(4) == 0.)
         track->SetPt(TMath::Sqrt(-1)); /*NaN*/
      else
         track->SetPt(1. / kftrack->GetParam(4)); /*signed Pt*/

      track->SetTheta(TMath::PiOver2() - kftrack->GetParam(3)); // Theta: angle from beam line
      track->SetPhi(kftrack->GetParam(2));                      // Phi

      TMatrixD Cov = *kftrack->GetCovariance();     // Error matrix
      track->SetPtError(TMath::Sqrt(Cov(4, 4)));    // Pt error
      track->SetThetaError(TMath::Sqrt(Cov(3, 3))); // Theta error
      track->SetPhiError(TMath::Sqrt(Cov(2, 2)));   // Phi error

      track->SetChi2(kftrack->GetChi2());

      Double_t phi = kftrack->GetParam(0) / kftrack->GetPosNew();
      track->SetFirstPointX(kftrack->GetPosNew() * TMath::Cos(phi)); // closest to beam line
      track->SetFirstPointY(kftrack->GetPosNew() * TMath::Sin(phi));
      track->SetFirstPointZ(kftrack->GetParam(1));

      track->SetLastPointX(0.); // AZ - currently not available
      track->SetLastPointY(0.); // AZ - currently not available
      track->SetLastPointZ(0.); // AZ - currently not available
      FillTrackDCA(track, &recoVertex, &mcVertex);
      FillTrackPID(track);
      FillTrackTpcHits(i, track);
   }

   Int_t               nEMatching = fEtofMatching ? fEtofMatching->GetEntries() : 0;
   MpdTofMatchingData *pEMatchingData;

   for (Int_t i = 0; i < nEctReco; i++) {
      MpdKalmanTrack *kftrack = (MpdKalmanTrack *)fKFEctTracks->UncheckedAt(i);

      MpdTrack *track = fEvent->AddGlobalTrack();
      track->SetID(kftrack->GetTrackID());
      track->SetNofHits(kftrack->GetNofHits());
      matchingDataExist = false;
      for (Int_t etofIndex = 0; etofIndex < nEMatching; etofIndex++) {
         pEMatchingData = (MpdTofMatchingData *)fEtofMatching->UncheckedAt(etofIndex);
         if (pEMatchingData->GetKFTrackIndex() == i) {
            matchingDataExist = true;
            break;
         } // first matching
      }

      if (matchingDataExist) {
         track->SetTofMass2(pEMatchingData->GetMass2());
         track->SetTofHitIndex(pEMatchingData->GetTofHitIndex());
      }

      if (kftrack->GetParam(4) == 0.)
         track->SetPt(TMath::Sqrt(-1)); /*NaN*/
      else
         track->SetPt(1. / kftrack->GetParam(4)); /*signed Pt*/

      track->SetTheta(TMath::PiOver2() - kftrack->GetParam(3)); // Theta: angle from beam line
      track->SetPhi(kftrack->GetParam(2));                      // Phi

      TMatrixD Cov = *kftrack->GetCovariance();     // Error matrix
      track->SetPtError(TMath::Sqrt(Cov(4, 4)));    // Pt error
      track->SetThetaError(TMath::Sqrt(Cov(3, 3))); // Theta error
      track->SetPhiError(TMath::Sqrt(Cov(2, 2)));   // Phi error

      track->SetChi2(kftrack->GetChi2());

      Double_t phi = kftrack->GetParam(0) / kftrack->GetPosNew();
      track->SetFirstPointX(kftrack->GetPosNew() * TMath::Cos(phi)); // closest to beam line
      track->SetFirstPointY(kftrack->GetPosNew() * TMath::Sin(phi));
      track->SetFirstPointZ(kftrack->GetParam(1));

      track->SetLastPointX(0.); // AZ - currently not available
      track->SetLastPointY(0.); // AZ - currently not available
      track->SetLastPointZ(0.); // AZ - currently not available
      FillTrackDCA(track, &recoVertex, &mcVertex);
      FillTrackPID(track);
   }
}

// -------------------------------------------------------------------
// Delete MC tracks being outside the MPD

void MpdFillDstTask::CleanMC()
{
   cout << "-I- [MpdFillDstTask::Exec] Cleaning from outer decays..." << flush;

   Int_t       nMC = fMCTracks->GetEntriesFast(), motherId;
   MpdMCTrack *track;

   // reading outer radius of TOF
   FairRunAna::Instance()->GetRuntimeDb()->getContainer("FairBaseParSet");
   assert(gGeoManager);
   TGeoTube *geoTube = (TGeoTube *)gGeoManager->GetVolume("tof1")->GetShape();
   Double_t  rMax    = geoTube->GetRmax();
   // cout <<"rMax="<<rMax<<endl;

   vector<int> removedMother;
   Int_t       i; //, j;
   for (i = 0; i < nMC; i++) {
      track = (MpdMCTrack *)fMCTracks->UncheckedAt(i);

      // filling hostograms
      motherId = track->GetMotherId();
      // fhTrackMotherId->Fill(motherId);
      if (motherId == -1) {
         // fhTrackPrimaryPDG->Fill(track->GetPdgCode());
      } else {
         // fhTrackVertex->Fill(track->GetStartX(),track->GetStartY());
         // if decay was out of MPD
         if (sqrt(track->GetStartX() * track->GetStartX() + track->GetStartY() * track->GetStartY()) > rMax) {
            fMCTracks->Remove(track);
            removedMother.push_back(i);
         } else {
            // motherId < childId (i) REQUIRED
            /*for (j = 0; j < cntRemovedMother; j++){
                if (motherId == removedMother[j]){
                    fMCTracks->Remove(track);
                    removedMother[cntRemovedMother] = i;
                    cntRemovedMother++;
                    break;
                }
            }
            if (j == cntRemovedMother) fhTruthVertex->Fill(track->GetStartX(),track->GetStartY());
            else cnt++;*/
            /*Int_t left = 0, right = cntRemovedMother-1, middle;
            while (left <= right){
                middle = (left + right)/2;
                if (motherId > removedMother[middle])
                    left = middle + 1;
                else{
                    if (motherId < removedMother[middle])
                        right = middle - 1;
                    else
                        break;
                }
            }
            if (left > right) fhTruthVertex->Fill(track->GetStartX(),track->GetStartY());
            else{
                fMCTracks->Remove(track);
                removedMother[cntRemovedMother] = i;
                cntRemovedMother++;
            }*/
            // binary search
            if (binary_search(removedMother.begin(), removedMother.end(), motherId)) {
               fMCTracks->Remove(track);
               removedMother.push_back(i);
            }
            // else
            //     fhTruthVertex->Fill(track->GetStartX(),track->GetStartY());
         }
      }
   }

   fMCTracks->Compress();
   cout << endl;
}

// -------------------------------------------------------------------

void MpdFillDstTask::Reset()
{
   fEvent->Reset();
}
// -------------------------------------------------------------------

void MpdFillDstTask::Finish()
{
   // cout<<"\n-I- [MpdFillDstTask::Finish] "<< endl;
   // fEvents->Dump();
   // cout << "\n";

   /* !this code gives some problems! -- !add extra event and so segfault!
   FairRootManager *fManager = FairRootManager::Instance();
   fManager->Fill();
    */
   // fhTrackMotherId->Write();
   // fhTrackPrimaryPDG->Write();
   // fhTrackVertex->Write();
   // fhTruthVertex->Write();

   delete fEvent;
   fEvent = nullptr;
}

void MpdFillDstTask::CalculateSharedArrayMap()
{
   if (fTpcHits->GetEntriesFast() > fSharedHitArraySize) {
      delete[] fSharedHitArray;
      fSharedHitArraySize = fTpcHits->GetEntriesFast();
      fSharedHitArray     = new Short_t[fSharedHitArraySize];
   }
   for (int i = 0; i < fTpcHits->GetEntriesFast(); i++) {
      fSharedHitArray[i] = 0;
   }
   for (int i = 0; i < fKFTracks->GetEntriesFast(); i++) {
      MpdTpcKalmanTrack *kalman = (MpdTpcKalmanTrack *)fKFTracks->UncheckedAt(i);
      TObjArray         *khits  = kalman->GetTrHits();
      if (khits)
         for (int j = 0; j < khits->GetEntriesFast(); j++) {
            MpdKalmanHit *hit = (MpdKalmanHit *)khits->UncheckedAt(j);
            Int_t         id  = hit->GetIndex(); // index of MpdTpcHit
            if (fSharedHitArray[id] < 3) fSharedHitArray[id]++;
         }
   }
}

void MpdFillDstTask::FillTrackTpcHits(Int_t particle_index, MpdTrack *track)
{
   MpdTpcKalmanTrack *kalman    = (MpdTpcKalmanTrack *)fKFTracks->UncheckedAt(particle_index);
   ULong64_t          layerHit  = 0;
   ULong64_t          sharedHit = 0;
   if (fTpcHits) { // calculate both maps
      TObjArray *khits = kalman->GetTrHits();
      if (khits)
         for (int j = 0; j < khits->GetEntriesFast(); j++) {
            MpdKalmanHit *hit   = (MpdKalmanHit *)khits->UncheckedAt(j);
            Int_t         id    = hit->GetIndex(); // index of MpdTpcHit
            Int_t         layer = hit->GetLayer();
            if (layer >= 0) {
               if (fSharedHitArray[id] > 1) SETBIT(sharedHit, layer);
               SETBIT(layerHit, layer);
            }
         }
   } else { // calculate only hit layer map
      TObjArray *khits = kalman->GetTrHits();
      if (khits)
         for (int j = 0; j < khits->GetEntriesFast(); j++) {
            MpdKalmanHit *hit   = (MpdKalmanHit *)khits->UncheckedAt(j);
            Int_t         layer = hit->GetLayer();
            if (layer >= 0) {
               SETBIT(layerHit, layer);
            }
         }
   }
   track->SetLayerHitMap(layerHit);
   track->SetSharedHitMap(sharedHit);
}

// -------------------------------------------------------------------

/*MpdEvent *MpdFillDstTask::AddEvent(Option_t * option)
{
  TClonesArray &events = *fEvents;
  Int_t size = events.GetEntriesFast(); // It is really always only one record
  MpdEvent* event = new(events[size]) MpdEvent();
  return event;
}*/

MpdTrack *MpdFillDstTask::AddPrimaryTrack()
{
   return nullptr;
}

void MpdFillDstTask::FillTrackDCA(MpdTrack *track, TVector3 *recoVertex, TVector3 *mcVertex)
{
   MpdHelix helix = track->GetHelix();
   Double_t path_at_mcVertex;
   Double_t path_at_recoVertex;
   path_at_mcVertex   = helix.pathLength(*mcVertex);
   path_at_recoVertex = helix.pathLength(*recoVertex);
   TVector3 DCA_MC    = helix.at(path_at_mcVertex);
   TVector3 DCA_RECO  = helix.at(path_at_recoVertex);
   // set dca global as dca to MC vertex DW
   track->SetDCAGlobalX(DCA_MC.X() - mcVertex->X());
   track->SetDCAGlobalY(DCA_MC.Y() - mcVertex->Y());
   track->SetDCAGlobalZ(DCA_MC.Z() - mcVertex->Z());
   // set dca as dca to reconstructed vertex DW
   track->SetDCAX(DCA_RECO.X() - recoVertex->X());
   track->SetDCAY(DCA_RECO.Y() - recoVertex->Y());
   track->SetDCAZ(DCA_RECO.Z() - recoVertex->Z());
}

void MpdFillDstTask::FillTrackPID(MpdTrack *track)
{
   TVector3 mom(track->GetPx(), track->GetPy(), track->GetPz());
   Double_t p        = mom.Mag();
   Double_t dedx     = track->GetdEdXTPC();
   Double_t dedx_el  = fPID->GetDedxElParam(p);
   Double_t dedx_pi  = fPID->GetDedxPiParam(p);
   Double_t dedx_ka  = fPID->GetDedxKaParam(p);
   Double_t dedx_pr  = fPID->GetDedxPrParam(p);
   Double_t sigma_el = fPID->GetDedxWidthValue(p, MpdPidUtils::kElectron) * dedx_el;
   Double_t sigma_pi = fPID->GetDedxWidthValue(p, MpdPidUtils::kPion) * dedx_pi;
   Double_t sigma_ka = fPID->GetDedxWidthValue(p, MpdPidUtils::kKaon) * dedx_ka;
   Double_t sigma_pr = fPID->GetDedxWidthValue(p, MpdPidUtils::kProton) * dedx_pr;

   if (sigma_el > FLT_EPSILON) track->SetNSigmaElectron((dedx - dedx_el) / sigma_el);

   if (sigma_pi > FLT_EPSILON) track->SetNSigmaPion((dedx - dedx_pi) / sigma_pi);

   if (sigma_ka > FLT_EPSILON) track->SetNSigmaKaon((dedx - dedx_ka) / sigma_ka);

   if (sigma_pr > FLT_EPSILON) track->SetNSigmaProton((dedx - dedx_pr) / sigma_pr);
}

// -------------------------------------------------------------------
ClassImp(MpdFillDstTask);

#include <iostream>
#include <fstream> // std::ifstream

#include "MpdMCTrack.h"
#include "MpdHelix.h"
#include "MpdTpcKalmanTrack.h"
#include "MpdEmcClusterKI.h"
#include "MpdParticle.h"
#include "MpdVertex.h"
#include "MpdEvent.h"
#include "MpdEmcGeoUtils.h"

#include "MpdPairPiLambda.h"
#include "TFile.h"

ClassImp(MpdPairPiLambda);

MpdPairPiLambda::MpdPairPiLambda(const char *name, const char *outputName) : MpdAnalysisTask(name, outputName)
{
   mParamConfig = outputName;
}

void MpdPairPiLambda::UserInit()
{
   cout << "[MpdPairPiLambda]: Initialization ... " << endl;

   mParams.ReadFromFile(mParamConfig);
   mParams.Print();

   // Prepare histograms etc.
   fOutputList = new TList();
   fOutputList->SetOwner(kTRUE);

   TH1::AddDirectory(kFALSE); // sets a global switch disabling the reference to histos in gROOT and their overwriting

   // General QA
   mhEvents = new TH1F("hEvents", "Number of events", 10, 0., 10.);
   fOutputList->Add(mhEvents);
   mhVertex = new TH1F("hVertex", "Event vertex distribution", 100, -200., 200.);
   fOutputList->Add(mhVertex);
   mhCentrality = new TH1F("hCentrality", "Centrality distribution", 100, 0., 100.);
   fOutputList->Add(mhCentrality);
   mhEvPlane = new TH1F("hEvPlane", "EvPlane distribution", 100, -4, 4);
   fOutputList->Add(mhEvPlane);

   // Minv
   mDcaPi = new TH2F("mDcaPi", "mDcaPi", 300, 0, 15, 300, 0, 15);
   fOutputList->Add(mDcaPi);
   mDcaPi1 = new TH2F("mDcaPi1", "mDcaPi1", 300, 0, 15, 300, 0, 15);
   fOutputList->Add(mDcaPi1);
   mDcaP = new TH2F("mDcaP", "mDcaP", 300, 0, 15, 300, 0, 15);
   fOutputList->Add(mDcaP);
   mDcaP1 = new TH2F("mDcaP1", "mDcaP1", 300, 0, 15, 300, 0, 15);
   fOutputList->Add(mDcaP1);

   mInvLambda = new TH2F("mInvLambda", "mInvLambda", 100, 0., 10., 1200, 1.0, 1.5);
   fOutputList->Add(mInvLambda);
   mInvLambda1 = new TH2F("mInvLambda1", "mInvLambda1", 100, 0., 10., 1200, 1.0, 1.5);
   fOutputList->Add(mInvLambda1);
   mInvLambda2 = new TH2F("mInvLambda2", "mInvLambda2", 100, 0., 10., 1200, 1.0, 1.5);
   fOutputList->Add(mInvLambda2);

   mAccEffRecTwoPID = new TH2F("mAccEffRecTwoPID", "mAccEffRecTwoPID", 100, -4, 4, 100, 0., 10.);
   fOutputList->Add(mAccEffRecTwoPID);

   mInvTwoPIDPip = new TH2F("mInvTwoPIDPip", "mInvTwoPIDPip", 100, 0., 10., 400, 0.0, 4.0);
   fOutputList->Add(mInvTwoPIDPip);
   mInvTwoPIDPim = new TH2F("mInvTwoPIDPim", "mInvTwoPIDPim", 100, 0., 10., 400, 0.0, 4.0);
   fOutputList->Add(mInvTwoPIDPim);
   mInvMixTwoPIDPip = new TH2F("mInvMixTwoPIDPip", "mInvMixTwoPIDPip", 100, 0., 10., 400, 0.0, 4.0);
   fOutputList->Add(mInvMixTwoPIDPip);
   mInvMixTwoPIDPim = new TH2F("mInvMixTwoPIDPim", "mInvMixTwoPIDPim", 100, 0., 10., 400, 0.0, 4.0);
   fOutputList->Add(mInvMixTwoPIDPim);

   for (int bin = 0; bin < nCenBinsAna; bin++) {
      mInvTwoPIDPipBin[bin] =
         new TH2F(Form("mInvTwoPIDPipBin_%d", bin), Form("mInvTwoPIDPipBin_%d", bin), 100, 0., 10., 400, 0.0, 4.0);
      fOutputList->Add(mInvTwoPIDPipBin[bin]);
      mInvTwoPIDPimBin[bin] =
         new TH2F(Form("mInvTwoPIDPimBin_%d", bin), Form("mInvTwoPIDPimBin_%d", bin), 100, 0., 10., 400, 0.0, 4.0);
      fOutputList->Add(mInvTwoPIDPimBin[bin]);

      mInvMixTwoPIDPipBin[bin] = new TH2F(Form("mInvMixTwoPIDPipBin_%d", bin), Form("mInvMixTwoPIDPipBin_%d", bin), 100,
                                          0., 10., 400, 0.0, 4.0);
      fOutputList->Add(mInvMixTwoPIDPipBin[bin]);
      mInvMixTwoPIDPimBin[bin] = new TH2F(Form("mInvMixTwoPIDPimBin_%d", bin), Form("mInvMixTwoPIDPimBin_%d", bin), 100,
                                          0., 10., 400, 0.0, 4.0);
      fOutputList->Add(mInvMixTwoPIDPimBin[bin]);
   }

   // MC
   if (isMC) {

      // Generated signals
      mInvGenSigmap = new TH1F("mInvGenSigmap", "mInvGenSigmap", 100, 0., 10.);
      fOutputList->Add(mInvGenSigmap);
      mInvGenSigmam = new TH1F("mInvGenSigmam", "mInvGenSigmam", 100, 0., 10.);
      fOutputList->Add(mInvGenSigmam);

      mInvTrueTwoPIDPip = new TH2F("mInvTrueTwoPIDPip", "mInvTrueTwoPIDPip", 100, 0., 10., 400, 0.0, 4.0);
      fOutputList->Add(mInvTrueTwoPIDPip);
      mInvTrueTwoPIDPim = new TH2F("mInvTrueTwoPIDPim", "mInvTrueTwoPIDPim", 100, 0., 10., 400, 0.0, 4.0);
      fOutputList->Add(mInvTrueTwoPIDPim);
      mInvTrueTwoPIDSigmap = new TH2F("mInvTrueTwoPIDSigmap", "mInvTrueTwoPIDSigmap", 100, 0., 10., 400, 0.0, 4.0);
      fOutputList->Add(mInvTrueTwoPIDSigmap);
      mMassResTwoPIDSigmap = new TH2F("mMassResTwoPIDSigmap", "mMassResTwoPIDSigmap", 100, 0., 10., 1000, -0.5, 0.5);
      fOutputList->Add(mMassResTwoPIDSigmap);
      mInvTrueTwoPIDSigmam = new TH2F("mInvTrueTwoPIDSigmam", "mInvTrueTwoPIDSigmam", 100, 0., 10., 400, 0.0, 4.0);
      fOutputList->Add(mInvTrueTwoPIDSigmam);
      mMassResTwoPIDSigmam = new TH2F("mMassResTwoPIDSigmam", "mMassResTwoPIDSigmam", 100, 0., 10., 1000, -0.5, 0.5);
      fOutputList->Add(mMassResTwoPIDSigmam);
      mAccEffGen = new TH2F("mAccEffGen", "mAccEffGen", 100, -4, 4, 100, 0., 10.);
      fOutputList->Add(mAccEffGen);

      for (int bin = 0; bin < nCenBinsAna; bin++) {
         mInvGenSigmapBin[bin] =
            new TH1F(Form("mInvGenSigmapBin_%d", bin), Form("mInvGenSigmapBin_%d", bin), 100, 0., 10.);
         fOutputList->Add(mInvGenSigmapBin[bin]);
         mInvGenSigmamBin[bin] =
            new TH1F(Form("mInvGenSigmamBin_%d", bin), Form("mInvGenSigmamBin_%d", bin), 100, 0., 10.);
         fOutputList->Add(mInvGenSigmamBin[bin]);
         mInvTrueTwoPIDSigmapBin[bin] = new TH2F(Form("mInvTrueTwoPIDSigmapBin_%d", bin),
                                                 Form("mInvTrueTwoPIDSigmapBin_%d", bin), 100, 0., 10., 400, 0.0, 4.0);
         fOutputList->Add(mInvTrueTwoPIDSigmapBin[bin]);
         mInvTrueTwoPIDSigmamBin[bin] = new TH2F(Form("mInvTrueTwoPIDSigmamBin_%d", bin),
                                                 Form("mInvTrueTwoPIDSigmamBin_%d", bin), 100, 0., 10., 400, 0.0, 4.0);
         fOutputList->Add(mInvTrueTwoPIDSigmamBin[bin]);
      }
   }

   for (long int i = 0; i < nMixTot; i++) {
      mixedEvents[i] = new TList();
   }

   cout << "[MpdPairPiLambda]: Reading Geo from sim_Geo.root for track refit ... " << endl;

   // Read-out TPC Geo for track Refit
   inFileSim = new TFile("sim_Geo.root", "READ");
   inTreeSim = (TTree *)inFileSim->Get("mpdsim");
   inTreeSim->SetName("mpdsim1");

   inFileSim->Get("FairGeoParSet");

   tpcPoints = (TClonesArray *)inFileSim->FindObjectAny("TpcPoint");

   inTreeSim->SetBranchAddress("TpcPoint", &tpcPoints);
   TBranch *tpcSimB = inTreeSim->GetBranch("TpcPoint");

   secGeo  = new TpcSectorGeoAZ();
   recoTpc = new MpdTpcKalmanFilter(*secGeo, "TPC Kalman filter");
   recoTpc->SetSectorGeo(*secGeo);
   recoTpc->FillGeoScheme();

   cout << "[MpdPairPiLambda]: Initialization done " << endl << endl;
}
//--------------------------------------
void MpdPairPiLambda::ProcessEvent(MpdAnalysisEvent &event)
{

   if (!isInitialized) {
      // mKF = MpdKalmanFilter::Instance();
      // mKHit.SetType(MpdKalmanHit::kFixedR);
      isInitialized = true;
   }

   if (!selectEvent(event)) {
      return;
   }

   mMpdGlobalTracks = event.fMPDEvent->GetGlobalTracks();
   mKalmanTracks    = event.fTPCKalmanTrack;

   if (isMC) {
      mMCTracks = event.fMCTrack;

      for (int i = 0; i < mMCTracks->GetEntriesFast(); i++) {
         MpdMCTrack *pr = (static_cast<MpdMCTrack *>(mMCTracks->At(i)));
         if (pr->GetPdgCode() == 3224 || pr->GetPdgCode() == 3114) {
            TVector3 momentum;
            pr->GetMomentum(momentum);
            mAccEffGen->Fill(pr->GetRapidity(), momentum.Pt());
            if (fabs(pr->GetRapidity()) < mParams.mYCut) {
               if (pr->GetStartX() * pr->GetStartX() + pr->GetStartY() * pr->GetStartY() < 1.) {
                  if (pr->GetPdgCode() == 3224) {
                     mInvGenSigmap->Fill(momentum.Pt());
                     mInvGenSigmapBin[anaBin]->Fill(momentum.Pt());
                  } // Sigmap
                  if (pr->GetPdgCode() == 3114) {
                     mInvGenSigmam->Fill(momentum.Pt());
                     mInvGenSigmamBin[anaBin]->Fill(momentum.Pt());
                  } // Sigmam
               }    // radius
            }       // rapidity
         }          // ID
      }
   } // isMC

   selectPionTrack(event);

   selectLambdaTrack(event);

   processHistograms(event);
}

void MpdPairPiLambda::Finish()
{
   // Post-scan processing not needed
}

//--------------------------------------
bool MpdPairPiLambda::selectEvent(MpdAnalysisEvent &event)
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

   cen = event.getCentrTPC();

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

   mRPBin = 0.5 * (event.fMpdEP.GetPhiEP_FHCal_F_all() / 3.14159 + 1) * nMixEventRP;
   if (mRPBin < 0) mRPBin = 0;
   if (mRPBin >= nMixEventRP) mRPBin = nMixEventRP - 1;

   mixBin = mZvtxBin * nMixEventCent * nMixEventRP + mCenBin * nMixEventRP + mRPBin;

   mhVertex->Fill(mPrimaryVertex.Z());
   mhCentrality->Fill(cen);

   mhEvPlane->Fill(event.fMpdEP.GetPhiEP_FHCal_F_all());

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
//--------------------------------------
void MpdPairPiLambda::selectPionTrack(MpdAnalysisEvent &event)
{
   // h+/-
   mP1.clear();

   int ntr = mMpdGlobalTracks->GetEntriesFast();

   for (long int i = 0; i < ntr; i++) {
      MpdTrack          *mpdtrack = (MpdTrack *)mMpdGlobalTracks->UncheckedAt(i);
      MpdTpcKalmanTrack *tr       = (MpdTpcKalmanTrack *)mKalmanTracks->UncheckedAt(i);
      if (!selectTrack(mpdtrack)) {
         continue;
      }

      int charge = mpdtrack->GetCharge();
      if (charge == 0) continue;

      long int trId = -1;
      if (isMC) {
         trId = mpdtrack->GetID();
      }

      int   isTOF    = -1;
      int   pid      = -1;
      int   isPi_TPC = -1;
      int   isPi_TOF = -1;
      float pmom     = sqrt(mpdtrack->GetPt() * mpdtrack->GetPt() + mpdtrack->GetPz() * mpdtrack->GetPz());

      if (fabs(mpdtrack->GetTofDphiSigma()) < 3.0 && fabs(mpdtrack->GetTofDzSigma()) < 3.0) isTOF = 1;

      if (fabs(mpdtrack->GetTPCNSigma(kPi)) < mParams.mPIDsigTPC) isPi_TPC = 1;
      if (isTOF == 1 && fabs(mpdtrack->GetTOFNSigma(kPi)) < mParams.mPIDsigTOF) isPi_TOF = 1;

      if (isPi_TPC == 1) pid = 1;
      if (isPi_TOF == 1) pid = 2;
      if (isPi_TPC == 1 || isPi_TOF == 1) pid = 3;
      if (isPi_TPC == 1 && isPi_TOF == 1) pid = 4;
      if ((isTOF != 1 && isPi_TPC == 1) || (isTOF == 1 && isPi_TPC == 1 && isPi_TOF == 1)) pid = 5;

      if (pid != 5) continue;

      float mPi = 0.13957;

      mP1.emplace_back(mpdtrack->GetPx(), mpdtrack->GetPy(), mpdtrack->GetPz(), sqrt(pmom * pmom + mPi * mPi));

      mP1.back().setPrimary1(trId);
      mP1.back().setPrimary2(-999);
      mP1.back().setTr1(i);
      mP1.back().setTr2(-999);
      mP1.back().setPID1(pid);
      mP1.back().setPID2(-999);
      mP1.back().setCh1(charge);
      mP1.back().setCh2(-999);
   } // i
}
//--------------------------------------
void MpdPairPiLambda::selectLambdaTrack(MpdAnalysisEvent &event)
{
   // Lambda ->p+pi-

   if (!eventM)
      eventM = new TClonesArray("MpdPairPiLambdaTrack", mMpdGlobalTracks->GetEntriesFast());
   else
      eventM->Clear();
   int iv = 0;

   mP2.clear();

   MpdVertex *vertex = (MpdVertex *)event.fVertex->First();

   int ntr = mMpdGlobalTracks->GetEntriesFast();

   for (long int i = 0; i < ntr; i++) {
      MpdTrack          *mpdtrack1 = (MpdTrack *)mMpdGlobalTracks->UncheckedAt(i);
      MpdTpcKalmanTrack *tr1       = (MpdTpcKalmanTrack *)mKalmanTracks->UncheckedAt(i);
      if (!selectTrackLambda(mpdtrack1)) {
         continue;
      }

      int charge1 = mpdtrack1->GetCharge();
      if (charge1 < 0) continue;

      long int trId1 = -1;
      if (isMC) {
         trId1 = mpdtrack1->GetID();
      }

      int isTOF1   = -1;
      int pid1     = -1;
      int isP_TPC1 = -1;
      int isP_TOF1 = -1;

      if (fabs(mpdtrack1->GetTofDphiSigma()) < 3.0 && fabs(mpdtrack1->GetTofDzSigma()) < 3.0) isTOF1 = 1;

      if (fabs(mpdtrack1->GetTPCNSigma(kP)) < mParams.mPIDsigTPC) isP_TPC1 = 1;
      if (isTOF1 == 1 && fabs(mpdtrack1->GetTOFNSigma(kP)) < mParams.mPIDsigTOF) isP_TOF1 = 1;

      if (isP_TPC1 == 1) pid1 = 1;
      if (isP_TOF1 == 1) pid1 = 2;
      if (isP_TPC1 == 1 || isP_TOF1 == 1) pid1 = 3;
      if (isP_TPC1 == 1 && isP_TOF1 == 1) pid1 = 4;
      if ((isTOF1 != 1 && isP_TPC1 == 1) || (isTOF1 == 1 && isP_TPC1 == 1 && isP_TOF1 == 1)) pid1 = 5;

      if (pid1 != 5) continue;

      float mP = 0.938272;

      MpdTpcKalmanTrack trCorP = *tr1;
      trCorP.SetDirection(MpdKalmanTrack::kInward);
      int ok1 = recoTpc->Refit(&trCorP, mP, 1); // refit
      if (!ok1) continue;

      MpdParticle proton(trCorP, i);
      proton.SetPdg(2212);
      proton.SetMass();

      double chi_p_pv = TMath::Min(proton.Chi2Vertex(vertex), 999.);
      mDcaP->Fill(chi_p_pv, sqrt(pow(mpdtrack1->GetNSigmaDCAx(), 2) + pow(mpdtrack1->GetNSigmaDCAy(), 2) +
                                 pow(mpdtrack1->GetNSigmaDCAz(), 2)));
      mDcaP1->Fill(chi_p_pv,
                   sqrt(pow(mpdtrack1->GetDCAX(), 2) + pow(mpdtrack1->GetDCAY(), 2) + pow(mpdtrack1->GetDCAZ(), 2)));

      if (chi_p_pv < mParams.mChi2ProtLam) continue;

      for (long int j = 0; j < ntr; j++) {
         MpdTrack          *mpdtrack2 = (MpdTrack *)mMpdGlobalTracks->UncheckedAt(j);
         MpdTpcKalmanTrack *tr2       = (MpdTpcKalmanTrack *)mKalmanTracks->UncheckedAt(j);
         if (!selectTrackLambda(mpdtrack2)) {
            continue;
         }

         int charge2 = mpdtrack2->GetCharge();
         if (charge2 > 0) continue;

         long int trId2 = -1;
         if (isMC) {
            trId2 = mpdtrack2->GetID();
         }

         int isTOF2    = -1;
         int pid2      = -1;
         int isPi_TPC2 = -1;
         int isPi_TOF2 = -1;

         if (fabs(mpdtrack2->GetTofDphiSigma()) < 3.0 && fabs(mpdtrack2->GetTofDzSigma()) < 3.0) isTOF2 = 1;

         if (fabs(mpdtrack2->GetTPCNSigma(kPi)) < mParams.mPIDsigTPC) isPi_TPC2 = 1;
         if (isTOF2 == 1 && fabs(mpdtrack2->GetTOFNSigma(kPi)) < mParams.mPIDsigTOF) isPi_TOF2 = 1;

         if (isPi_TPC2 == 1) pid2 = 1;
         if (isPi_TOF2 == 1) pid2 = 2;
         if (isPi_TPC2 == 1 || isPi_TOF2 == 1) pid2 = 3;
         if (isPi_TPC2 == 1 && isPi_TOF2 == 1) pid2 = 4;
         if ((isTOF2 != 1 && isPi_TPC2 == 1) || (isTOF2 == 1 && isPi_TPC2 == 1 && isPi_TOF2 == 1)) pid2 = 5;

         if (pid2 != 5) continue;

         MpdParticle pion(*tr2, j);
         pion.SetPdg(-211);
         pion.SetMass();

         double chi_pi_pv = TMath::Min(pion.Chi2Vertex(vertex), 999.);
         mDcaPi->Fill(chi_pi_pv, sqrt(pow(mpdtrack2->GetNSigmaDCAx(), 2) + pow(mpdtrack2->GetNSigmaDCAy(), 2) +
                                      pow(mpdtrack2->GetNSigmaDCAz(), 2)));
         mDcaPi1->Fill(chi_pi_pv, sqrt(pow(mpdtrack2->GetDCAX(), 2) + pow(mpdtrack2->GetDCAY(), 2) +
                                       pow(mpdtrack2->GetDCAZ(), 2)));
         if (chi_pi_pv < mParams.mChi2PionLam) continue;

         // Build Lambda
         vPartK.clear();
         vPartK.emplace_back(&proton);
         vPartK.emplace_back(&pion);

         MpdParticle LamPPi;
         double      chi2 = TMath::Abs(LamPPi.BuildMother(vPartK));

         TVector3 v0(LamPPi.Getx()(0, 0), LamPPi.Getx()(1, 0), LamPPi.Getx()(2, 0));
         v0 -= mPrimaryVertex;
         double decay = v0.Mag();

         double disth, angle;
         angle = v0.Angle(LamPPi.Momentum3());

         double pz_tmp = TMath::Sqrt(LamPPi.Momentum() * LamPPi.Momentum() - LamPPi.Pt() * LamPPi.Pt());

         MpdHelix             helix1 = MakeHelix(&trCorP);
         MpdHelix             helix2 = MakeHelix(tr2);
         pair<double, double> paths  = helix1.pathLengths(helix2);
         TVector3             p1     = helix1.at(paths.first);
         TVector3             p2     = helix2.at(paths.second);
         p1 -= p2;
         disth = p1.Mag();

         if (fabs(LamPPi.Rapidity()) > mParams.mLamEtaCut) continue;
         if (chi2 > mParams.mChi2Lam) continue;
         mInvLambda1->Fill(LamPPi.Pt(), LamPPi.GetMass());
         if (angle > mParams.mPALam) continue;
         if (decay < mParams.mDecayLam) continue;
         if (disth > mParams.mDistLam) continue;
         mInvLambda->Fill(LamPPi.Pt(), LamPPi.GetMass());

         float mLambda = 1.115683;
         if (fabs(LamPPi.GetMass() - mLambda) > mParams.mNSigmaLam * mParams.mWidthLam) continue;
         mInvLambda2->Fill(LamPPi.Pt(), LamPPi.GetMass());

         float px_lam = (LamPPi.Pt()) * TMath::Cos(LamPPi.Phi());
         float py_lam = (LamPPi.Pt()) * TMath::Sin(LamPPi.Phi());
         float pz_lam = TMath::Sign(pz_tmp, TMath::Cos(LamPPi.Theta()));
         // float E_lam = sqrt(px_lam*px_lam + py_lam*py_lam + pz_lam*pz_lam + LamPPi.GetMass()*LamPPi.GetMass());
         float E_lam = sqrt(px_lam * px_lam + py_lam * py_lam + pz_lam * pz_lam + mLambda * mLambda);

         mP2.emplace_back(px_lam, py_lam, pz_lam, E_lam);
         mP2.back().setPrimary1(trId1);
         mP2.back().setPrimary2(trId2);
         mP2.back().setTr1(i);
         mP2.back().setTr2(j);
         mP2.back().setPID1(pid1);
         mP2.back().setPID2(pid2);
         mP2.back().setCh1(1);
         mP2.back().setCh2(-1);

         new ((*eventM)[iv++]) MpdPairPiLambdaTrack(mP2.back());

      } // j
   }    // i
}

void MpdPairPiLambda::processHistograms(MpdAnalysisEvent &event)
{
   // Fill Real, Mixed distributions and update mixed array

   // Real
   int nPions = mP1.size();
   int nLam   = mP2.size();

   for (int i = 0; i < nPions; i++) {
      for (int j = 0; j < nLam; j++) {

         if (mP1[i].getTr1() == mP2[j].getTr1() || mP1[i].getTr1() == mP2[j].getTr2()) continue;

         TLorentzVector sum = mP1[i] + mP2[j];

         if (fabs(sum.Rapidity()) < mParams.mYCut) {

            if (mP1[i].getPID1() == 5 && mP2[j].getPID1() == 5 && mP2[j].getPID2() == 5) {
               if (mP1[i].getCh1() > 0) {
                  mInvTwoPIDPip->Fill(sum.Pt(), sum.M());
                  mInvTwoPIDPipBin[anaBin]->Fill(sum.Pt(), sum.M());
               }
               if (mP1[i].getCh1() < 0) {
                  mInvTwoPIDPim->Fill(sum.Pt(), sum.M());
                  mInvTwoPIDPimBin[anaBin]->Fill(sum.Pt(), sum.M());
               }
            }
         } // rapidity

         long int ip_lam = IsSameParent(mP2[j].primary1(), mP2[j].primary2());

         if (ip_lam >= 0) {
            long int ip = IsSameParent(mP1[i].primary1(), ip_lam);

            if (ip >= 0) {
               if (fabs(sum.Rapidity()) < mParams.mYCut) {
                  if (mP1[i].getPID1() == 5 && mP2[j].getPID1() == 5 && mP2[j].getPID2() == 5) {
                     if (mP1[i].getCh1() > 0) {
                        mInvTrueTwoPIDPip->Fill(sum.Pt(), sum.M());
                     }
                     if (mP1[i].getCh1() < 0) {
                        mInvTrueTwoPIDPim->Fill(sum.Pt(), sum.M());
                     }
                  } // pid
               }    // rapidity

               // Sigma+ -> Lambda + pi+
               if (static_cast<MpdMCTrack *>(mMCTracks->At(ip))->GetPdgCode() == 3224) {

                  MpdMCTrack *pr1 = (static_cast<MpdMCTrack *>(mMCTracks->At(mP1[i].primary1())));
                  MpdMCTrack *pr2 = (static_cast<MpdMCTrack *>(mMCTracks->At(ip_lam)));

                  if (pr1->GetPdgCode() == 211 && pr2->GetPdgCode() == 3122) {

                     TLorentzVector p1, p2, sump1p2;
                     p1.SetPxPyPzE(pr1->GetPx(), pr1->GetPy(), pr1->GetPz(), pr1->GetEnergy());
                     p2.SetPxPyPzE(pr2->GetPx(), pr2->GetPy(), pr2->GetPz(), pr2->GetEnergy());

                     sump1p2 = p1 + p2;

                     if (mP1[i].getPID1() == 5 && mP2[j].getPID1() == 5 && mP2[j].getPID2() == 5) {
                        mAccEffRecTwoPID->Fill(sum.Rapidity(), sum.Pt());
                     }

                     if (fabs(sum.Rapidity()) < mParams.mYCut) {
                        if (mP1[i].getPID1() == 5 && mP2[j].getPID1() == 5 && mP2[j].getPID2() == 5) {
                           mInvTrueTwoPIDSigmap->Fill(sum.Pt(), sum.M());
                           mInvTrueTwoPIDSigmapBin[anaBin]->Fill(sum.Pt(), sum.M());
                           mMassResTwoPIDSigmap->Fill(sum.Pt(), sum.M() - sump1p2.M());
                        } // pid
                     }    // rapidity
                  }       // pi+ + Lambda
               }          // Sigma+

               // Sigma- -> Lambda + pi-
               if (static_cast<MpdMCTrack *>(mMCTracks->At(ip))->GetPdgCode() == 3114) {

                  MpdMCTrack *pr1 = (static_cast<MpdMCTrack *>(mMCTracks->At(mP1[i].primary1())));
                  MpdMCTrack *pr2 = (static_cast<MpdMCTrack *>(mMCTracks->At(ip_lam)));

                  if (pr1->GetPdgCode() == -211 && pr2->GetPdgCode() == 3122) {

                     TLorentzVector p1, p2, sump1p2;
                     p1.SetPxPyPzE(pr1->GetPx(), pr1->GetPy(), pr1->GetPz(), pr1->GetEnergy());
                     p2.SetPxPyPzE(pr2->GetPx(), pr2->GetPy(), pr2->GetPz(), pr2->GetEnergy());

                     sump1p2 = p1 + p2;

                     if (mP1[i].getPID1() == 5 && mP2[j].getPID1() == 5 && mP2[j].getPID2() == 5) {
                        mAccEffRecTwoPID->Fill(sum.Rapidity(), sum.Pt());
                     }

                     if (fabs(sum.Rapidity()) < mParams.mYCut) {
                        if (mP1[i].getPID1() == 5 && mP2[j].getPID1() == 5 && mP2[j].getPID2() == 5) {
                           mInvTrueTwoPIDSigmam->Fill(sum.Pt(), sum.M());
                           mInvTrueTwoPIDSigmamBin[anaBin]->Fill(sum.Pt(), sum.M());
                           mMassResTwoPIDSigmam->Fill(sum.Pt(), sum.M() - sump1p2.M());
                        } // pid
                     }    // rapidity
                  }       // pi- + Lambda
               }          // Sigma-

            } // same parent for pi+ and Lambda
         }    // same parent for Lambda

      } // nLam

      // Mixing
      TIter         next(mixedEvents[mixBin]);
      TClonesArray *ma;
      while ((ma = (TClonesArray *)next())) {
         for (Int_t k = 0; k < ma->GetEntriesFast(); k++) {
            MpdPairPiLambdaTrack *p2 = (MpdPairPiLambdaTrack *)ma->At(k);
            double                m  = (mP1[i] + *p2).M();
            double                pt = (mP1[i] + *p2).Pt();

            if (fabs((mP1[i] + *p2).Rapidity()) < mParams.mYCut) {

               if (mP1[i].getPID1() == 5 && p2->getPID1() == 5 && p2->getPID2() == 5) {
                  if (mP1[i].getCh1() > 0) {
                     mInvMixTwoPIDPip->Fill(pt, m);
                     mInvMixTwoPIDPipBin[anaBin]->Fill(pt, m);
                  }
                  if (mP1[i].getCh1() < 0) {
                     mInvMixTwoPIDPim->Fill(pt, m);
                     mInvMixTwoPIDPimBin[anaBin]->Fill(pt, m);
                  }
               }
            } // rapidity
         }    // k
      }       // ma
   }          // nPions

   // Remove redandant events
   if (mixedEvents[mixBin]->GetSize() > nMixed) {
      TClonesArray *tmp = (TClonesArray *)mixedEvents[mixBin]->Last();
      mixedEvents[mixBin]->Remove(tmp);
      delete tmp;
   }
   mixedEvents[mixBin]->AddFirst(eventM);
   eventM = nullptr;

} // histos

bool MpdPairPiLambda::selectTrack(MpdTrack *mpdtrack)
{
   if (mpdtrack->GetNofHits() < mParams.mNofHitsCut) return false;

   if (fabs(mpdtrack->GetEta()) > mParams.mEtaCut) return false;

   float pt = fabs(mpdtrack->GetPt());
   if (pt < mParams.mPtminCut) return false;

   if (fabs(mpdtrack->GetNSigmaDCAx()) > mParams.mDCACut) return false; // |DCAx| < 2*sigmaxy(pT, eta, centrality).
   if (fabs(mpdtrack->GetNSigmaDCAy()) > mParams.mDCACut) return false; // |DCAy| < 2*sigmaxy(pT, eta, centrality).
   if (fabs(mpdtrack->GetNSigmaDCAz()) > mParams.mDCACut) return false; // |DCAz| < 2*sigmaz(pT, eta, centrality).

   return true;
}

bool MpdPairPiLambda::selectTrackLambda(MpdTrack *mpdtrack)
{
   if (mpdtrack->GetNofHits() < mParams.mNofHitsCut) return false;

   if (fabs(mpdtrack->GetEta()) > mParams.mEtaCut) return false;

   float pt = fabs(mpdtrack->GetPt());
   if (pt < mParams.mPtminCut) return false;

   return true;
}

long int MpdPairPiLambda::IsSameParent(long int prim1, long int prim2) const
{
   // Looks through parents and finds if there was commont pi0 among ancestors

   if (!isMC) return -1; // can not say anything

   while (prim1 != -1) {
      long int pr2 = prim2;

      while (pr2 != -1) {
         if (prim1 == pr2) {
            return prim1;
         }
         pr2 = (static_cast<MpdMCTrack *>(mMCTracks->At(pr2)))->GetMotherId();
      }
      prim1 = (static_cast<MpdMCTrack *>(mMCTracks->At(prim1)))->GetMotherId();
   }
   return -1;
}

MpdHelix MpdPairPiLambda::MakeHelix(const MpdTpcKalmanTrack *tr)
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

#include <iostream>
#include <fstream> // std::ifstream

#include "MpdMCTrack.h"
#include "MpdTpcKalmanTrack.h"
#include "MpdEmcClusterKI.h"
#include "MpdParticle.h"
#include "MpdVertex.h"
#include "MpdEvent.h"
#include "MpdEmcGeoUtils.h"

#include "MpdPairPiK.h"
#include "TFile.h"

ClassImp(MpdPairPiK);

MpdPairPiK::MpdPairPiK(const char *name, const char *outputName) : MpdAnalysisTask(name, outputName)
{
   mParamConfig = outputName;
}

void MpdPairPiK::UserInit()
{
   cout << "[MpdPairPiK]: Initialization ... " << endl;

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
   mInvNoPID = new TH2F("mInvNoPID", "mInvNoPID", 100, 0., 10., 250, 0.5, 4.0);
   fOutputList->Add(mInvNoPID);
   mInvOnePID = new TH2F("mInvOnePID", "mInvOnePID", 100, 0., 10., 250, 0.5, 4.0);
   fOutputList->Add(mInvOnePID);
   mInvTwoPID = new TH2F("mInvTwoPID", "mInvTwoPID", 100, 0., 10., 250, 0.5, 4.0);
   fOutputList->Add(mInvTwoPID);
   mAccEffRecNoPID = new TH2F("mAccEffRecNoPID", "mAccEffRecNoPID", 100, -4, 4, 100, 0., 10.);
   fOutputList->Add(mAccEffRecNoPID);
   mAccEffRecOnePID = new TH2F("mAccEffRecOnePID", "mAccEffRecOnePID", 100, -4, 4, 100, 0., 10.);
   fOutputList->Add(mAccEffRecOnePID);
   mAccEffRecTwoPID = new TH2F("mAccEffRecTwoPID", "mAccEffRecTwoPID", 100, -4, 4, 100, 0., 10.);
   fOutputList->Add(mAccEffRecTwoPID);

   // MinvMix
   mInvMixNoPID = new TH2F("mInvMixNoPID", "mInvMixNoPID", 100, 0., 10., 250, 0.5, 4.0);
   fOutputList->Add(mInvMixNoPID);
   mInvMixOnePID = new TH2F("mInvMixOnePID", "mInvMixOnePID", 100, 0., 10., 250, 0.5, 4.0);
   fOutputList->Add(mInvMixOnePID);
   mInvMixTwoPID = new TH2F("mInvMixTwoPID", "mInvMixTwoPID", 100, 0., 10., 250, 0.5, 4.0);
   fOutputList->Add(mInvMixTwoPID);

   for (int bin = 0; bin < nCenBinsAna; bin++) {
      mInvNoPIDBin[bin] =
         new TH2F(Form("mInvNoPIDBin_%d", bin), Form("mInvNoPIDBin_%d", bin), 100, 0., 10., 250, 0.5, 4.0);
      fOutputList->Add(mInvNoPIDBin[bin]);
      mInvOnePIDBin[bin] =
         new TH2F(Form("mInvOnePIDBin_%d", bin), Form("mInvOnePIDBin_%d", bin), 100, 0., 10., 250, 0.5, 4.0);
      fOutputList->Add(mInvOnePIDBin[bin]);
      mInvTwoPIDBin[bin] =
         new TH2F(Form("mInvTwoPIDBin_%d", bin), Form("mInvTwoPIDBin_%d", bin), 100, 0., 10., 250, 0.5, 4.0);
      fOutputList->Add(mInvTwoPIDBin[bin]);

      mInvMixNoPIDBin[bin] =
         new TH2F(Form("mInvMixNoPIDBin_%d", bin), Form("mInvMixNoPIDBin_%d", bin), 100, 0., 10., 250, 0.5, 4.0);
      fOutputList->Add(mInvMixNoPIDBin[bin]);
      mInvMixOnePIDBin[bin] =
         new TH2F(Form("mInvMixOnePIDBin_%d", bin), Form("mInvMixOnePIDBin_%d", bin), 100, 0., 10., 250, 0.5, 4.0);
      fOutputList->Add(mInvMixOnePIDBin[bin]);
      mInvMixTwoPIDBin[bin] =
         new TH2F(Form("mInvMixTwoPIDBin_%d", bin), Form("mInvMixTwoPIDBin_%d", bin), 100, 0., 10., 250, 0.5, 4.0);
      fOutputList->Add(mInvMixTwoPIDBin[bin]);
   }

   // MC
   if (isMC) {

      // Generated signals
      mInvGen = new TH1F("mInvGen", "mInvGen", 100, 0., 10.);
      fOutputList->Add(mInvGen);

      mInvTrueNoPID = new TH2F("mInvTrueNoPID", "mInvTrueNoPID", 100, 0., 10., 250, 0.5, 4.0);
      fOutputList->Add(mInvTrueNoPID);
      mInvTrueNoPIDKst892 = new TH2F("mInvTrueNoPIDKst892", "mInvTrueNoPIDKst892", 100, 0., 10., 250, 0.5, 4.0);
      fOutputList->Add(mInvTrueNoPIDKst892);
      mMassResNoPIDKst892 = new TH2F("mMassResNoPIDKst892", "mMassResNoPIDKst892", 100, 0., 10., 1000, -0.5, 0.5);
      fOutputList->Add(mMassResNoPIDKst892);
      mInvTrueOnePID = new TH2F("mInvTrueOnePID", "mInvTrueOnePID", 100, 0., 10., 250, 0.5, 4.0);
      fOutputList->Add(mInvTrueOnePID);
      mInvTrueOnePIDKst892 = new TH2F("mInvTrueOnePIDKst892", "mInvTrueOnePIDKst892", 100, 0., 10., 250, 0.5, 4.0);
      fOutputList->Add(mInvTrueOnePIDKst892);
      mMassResOnePIDKst892 = new TH2F("mMassResOnePIDKst892", "mMassResOnePIDKst892", 100, 0., 10., 1000, -0.5, 0.5);
      fOutputList->Add(mMassResOnePIDKst892);
      mInvTrueTwoPID = new TH2F("mInvTrueTwoPID", "mInvTrueTwoPID", 100, 0., 10., 250, 0.5, 4.0);
      fOutputList->Add(mInvTrueTwoPID);
      mInvTrueTwoPIDKst892 = new TH2F("mInvTrueTwoPIDKst892", "mInvTrueTwoPIDKst892", 100, 0., 10., 250, 0.5, 4.0);
      fOutputList->Add(mInvTrueTwoPIDKst892);
      mMassResTwoPIDKst892 = new TH2F("mMassResTwoPIDKst892", "mMassResTwoPIDKst892", 100, 0., 10., 1000, -0.5, 0.5);
      fOutputList->Add(mMassResTwoPIDKst892);
      mAccEffGen = new TH2F("mAccEffGen", "mAccEffGen", 100, -4, 4, 100, 0., 10.);
      fOutputList->Add(mAccEffGen);

      for (int bin = 0; bin < nCenBinsAna; bin++) {
         mInvGenBin[bin] = new TH1F(Form("mInvGenBin_%d", bin), Form("mInvGenBin_%d", bin), 100, 0., 10.);
         fOutputList->Add(mInvGenBin[bin]);
         mInvTrueNoPIDKst892Bin[bin] = new TH2F(Form("mInvTrueNoPIDKst892Bin_%d", bin),
                                                Form("mInvTrueNoPIDKst892Bin_%d", bin), 100, 0., 10., 250, 0.5, 4.0);
         fOutputList->Add(mInvTrueNoPIDKst892Bin[bin]);
         mInvTrueOnePIDKst892Bin[bin] = new TH2F(Form("mInvTrueOnePIDKst892Bin_%d", bin),
                                                 Form("mInvTrueOnePIDKst892Bin_%d", bin), 100, 0., 10., 250, 0.5, 4.0);
         fOutputList->Add(mInvTrueOnePIDKst892Bin[bin]);
         mInvTrueTwoPIDKst892Bin[bin] = new TH2F(Form("mInvTrueTwoPIDKst892Bin_%d", bin),
                                                 Form("mInvTrueTwoPIDKst892Bin_%d", bin), 100, 0., 10., 250, 0.5, 4.0);
         fOutputList->Add(mInvTrueTwoPIDKst892Bin[bin]);
      }
   }

   for (long int i = 0; i < nMixTot; i++) {
      mixedEvents[i] = new TList();
   }

   cout << "[MpdPairPiK]: Reading Geo from sim_1.root for track refit ... " << endl;

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

   cout << "[MpdPairPiK]: Initialization done " << endl << endl;
}
//--------------------------------------
void MpdPairPiK::ProcessEvent(MpdAnalysisEvent &event)
{

   if (!isInitialized) {
      // mKF = MpdKalmanFilter::Instance();
      // mKHit.SetType(MpdKalmanHit::kFixedR);
      isInitialized = true;
   }

   if (!selectEvent(event)) {
      return;
   }

   mKalmanTracks  = event.fTPCKalmanTrack;
   mpdTofMatching = event.fTOFMatching;

   if (isMC) {
      mMCTracks = event.fMCTrack;

      for (int i = 0; i < mMCTracks->GetEntriesFast(); i++) {
         MpdMCTrack *pr = (static_cast<MpdMCTrack *>(mMCTracks->At(i)));
         if (pr->GetPdgCode() == 313) {
            TVector3 momentum;
            pr->GetMomentum(momentum);
            mAccEffGen->Fill(pr->GetRapidity(), momentum.Pt());
            if (fabs(pr->GetRapidity()) < mParams.mYCut) {
               if (pr->GetStartX() * pr->GetStartX() + pr->GetStartY() * pr->GetStartY() < 1.) {
                  mInvGen->Fill(momentum.Pt());
                  mInvGenBin[anaBin]->Fill(momentum.Pt());
               } // radius
            }    // rapidity
         }       // ID
      }
   } // isMC

   selectPosTrack(event);

   selectNegTrack(event);

   processHistograms(event);
}

void MpdPairPiK::Finish()
{
   // Post-scan processing not needed
}

//--------------------------------------
bool MpdPairPiK::selectEvent(MpdAnalysisEvent &event)
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
void MpdPairPiK::selectPosTrack(MpdAnalysisEvent &event)
{
   // h+
   mP1.clear();

   mMpdGlobalTracks = event.fMPDEvent->GetGlobalTracks();
   int ntr          = mMpdGlobalTracks->GetEntriesFast();

   for (long int i = 0; i < ntr; i++) {
      MpdTrack          *mpdtrack = (MpdTrack *)mMpdGlobalTracks->UncheckedAt(i);
      MpdTpcKalmanTrack *tr       = (MpdTpcKalmanTrack *)mKalmanTracks->UncheckedAt(i);
      if (!selectTrack(mpdtrack)) {
         continue;
      }

      int charge = mpdtrack->GetCharge();
      if (charge < 0) continue;

      long int trId = -1;
      if (isMC) {
         trId = mpdtrack->GetID();
      }

      int   isTOF   = -1;
      int   pid     = -1;
      int   isK_TPC = -1;
      int   isK_TOF = -1;
      float pmom    = sqrt(mpdtrack->GetPt() * mpdtrack->GetPt() + mpdtrack->GetPz() * mpdtrack->GetPz());

      if (fabs(mpdtrack->GetTofDphiSigma()) < 3.0 && fabs(mpdtrack->GetTofDzSigma()) < 3.0) isTOF = 1;

      if (fabs(mpdtrack->GetTPCNSigma(kK)) < mParams.mPIDsigTPC) isK_TPC = 1;
      if (isTOF == 1 && fabs(mpdtrack->GetTOFNSigma(kK)) < mParams.mPIDsigTOF) isK_TOF = 1;

      if (isK_TPC == 1) pid = 1;
      if (isK_TOF == 1) pid = 2;
      if (isK_TPC == 1 || isK_TOF == 1) pid = 3;
      if (isK_TPC == 1 && isK_TOF == 1) pid = 4;
      if ((isTOF != 1 && isK_TPC == 1) || (isTOF == 1 && isK_TPC == 1 && isK_TOF == 1)) pid = 5;

      if (pid < 0) continue;

      float mK = 0.493677;

      MpdTpcKalmanTrack trCorK = *tr;
      trCorK.SetDirection(MpdKalmanTrack::kInward);
      int ok = recoTpc->Refit(&trCorK, mK, 1); // refit
      if (ok) {
         MpdParticle kaon(trCorK, i);
         kaon.SetPdg(321);
         kaon.SetMass();

         pmom = sqrt(kaon.Momentum3().Pt() * kaon.Momentum3().Pt() + kaon.Momentum3().Pz() * kaon.Momentum3().Pz());

         mP1.emplace_back(kaon.Momentum3().Px(), kaon.Momentum3().Py(), kaon.Momentum3().Pz(),
                          sqrt(pmom * pmom + mK * mK));
      } else {
         mP1.emplace_back(mpdtrack->GetPx(), mpdtrack->GetPy(), mpdtrack->GetPz(), sqrt(pmom * pmom + mK * mK));
      }

      mP1.back().setPrimary(trId);
      mP1.back().setTr1(i);
      mP1.back().setPID(pid);
   } // i
}
//--------------------------------------
void MpdPairPiK::selectNegTrack(MpdAnalysisEvent &event)
{
   // h-

   eventM = new TClonesArray("MpdPairPiKTrack", mMpdGlobalTracks->GetEntriesFast());
   int iv = 0;

   mP2.clear();

   mMpdGlobalTracks = event.fMPDEvent->GetGlobalTracks();
   int ntr          = mMpdGlobalTracks->GetEntriesFast();

   for (long int i = 0; i < ntr; i++) {
      MpdTrack          *mpdtrack = (MpdTrack *)mMpdGlobalTracks->UncheckedAt(i);
      MpdTpcKalmanTrack *tr       = (MpdTpcKalmanTrack *)mKalmanTracks->UncheckedAt(i);
      if (!selectTrack(mpdtrack)) {
         continue;
      }

      int charge = mpdtrack->GetCharge();
      if (charge > 0) continue;

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

      if (pid < 0) continue;

      float mPi = 0.13957;

      mP2.emplace_back(mpdtrack->GetPx(), mpdtrack->GetPy(), mpdtrack->GetPz(), sqrt(pmom * pmom + mPi * mPi));

      mP2.back().setPrimary(trId);
      mP2.back().setTr1(i);
      mP2.back().setPID(pid);

      new ((*eventM)[iv++]) MpdPairPiKTrack(mP2[mP2.size() - 1]);
   } // i
}

void MpdPairPiK::processHistograms(MpdAnalysisEvent &event)
{
   // Fill Real, Mixed distributions and update mixed array

   // Real
   int nPos = mP1.size();
   int nNeg = mP2.size();

   for (int i = 0; i < nPos; i++) {
      for (int j = 0; j < nNeg; j++) {
         TLorentzVector sum = mP1[i] + mP2[j];

         if (fabs(sum.Rapidity()) < mParams.mYCut) {
            mInvNoPID->Fill(sum.Pt(), sum.M());
            mInvNoPIDBin[anaBin]->Fill(sum.Pt(), sum.M());
            if (mP1[i].getPID() == 5 || mP2[j].getPID() == 5) {
               mInvOnePID->Fill(sum.Pt(), sum.M());
               mInvOnePIDBin[anaBin]->Fill(sum.Pt(), sum.M());
            }
            if (mP1[i].getPID() == 5 && mP2[j].getPID() == 5) {
               mInvTwoPID->Fill(sum.Pt(), sum.M());
               mInvTwoPIDBin[anaBin]->Fill(sum.Pt(), sum.M());
            }
         } // rapidity

         long int ip = IsSameParent(mP1[i].primary(), mP2[j].primary());
         if (ip >= 0) {

            if (fabs(sum.Rapidity()) < mParams.mYCut) {
               mInvTrueNoPID->Fill(sum.Pt(), sum.M());
               if (mP1[i].getPID() == 5 || mP2[j].getPID() == 5) mInvTrueOnePID->Fill(sum.Pt(), sum.M());
               if (mP1[i].getPID() == 5 && mP2[j].getPID() == 5) mInvTrueTwoPID->Fill(sum.Pt(), sum.M());
            } // rapidity

            if (static_cast<MpdMCTrack *>(mMCTracks->At(ip))->GetPdgCode() == 313) {

               MpdMCTrack *pr1 = (static_cast<MpdMCTrack *>(mMCTracks->At(mP1[i].primary())));
               MpdMCTrack *pr2 = (static_cast<MpdMCTrack *>(mMCTracks->At(mP2[j].primary())));

               if (pr1->GetPdgCode() == 321 && pr2->GetPdgCode() == -211) {

                  TLorentzVector p1, p2, sump1p2;
                  p1.SetPxPyPzE(pr1->GetPx(), pr1->GetPy(), pr1->GetPz(), pr1->GetEnergy());
                  p2.SetPxPyPzE(pr2->GetPx(), pr2->GetPy(), pr2->GetPz(), pr2->GetEnergy());

                  sump1p2 = p1 + p2;

                  mAccEffRecNoPID->Fill(sum.Rapidity(), sum.Pt());
                  if (mP1[i].getPID() == 5 || mP2[j].getPID() == 5) mAccEffRecOnePID->Fill(sum.Rapidity(), sum.Pt());
                  if (mP1[i].getPID() == 5 && mP2[j].getPID() == 5) mAccEffRecTwoPID->Fill(sum.Rapidity(), sum.Pt());

                  if (fabs(sum.Rapidity()) < mParams.mYCut) {
                     mInvTrueNoPIDKst892->Fill(sum.Pt(), sum.M());
                     mMassResNoPIDKst892->Fill(sum.Pt(), sum.M() - sump1p2.M());
                     mInvTrueNoPIDKst892Bin[anaBin]->Fill(sum.Pt(), sum.M());

                     if (mP1[i].getPID() == 5 || mP2[j].getPID() == 5) {
                        mInvTrueOnePIDKst892->Fill(sum.Pt(), sum.M());
                        mMassResOnePIDKst892->Fill(sum.Pt(), sum.M() - sump1p2.M());
                        mInvTrueOnePIDKst892Bin[anaBin]->Fill(sum.Pt(), sum.M());
                     }

                     if (mP1[i].getPID() == 5 && mP2[j].getPID() == 5) {
                        mInvTrueTwoPIDKst892->Fill(sum.Pt(), sum.M());
                        mMassResTwoPIDKst892->Fill(sum.Pt(), sum.M() - sump1p2.M());
                        mInvTrueTwoPIDKst892Bin[anaBin]->Fill(sum.Pt(), sum.M());
                     }
                  } // rapidity
               }    // K+pi-
            }       // 313
         }          // ip

      } // Neg

      // Mixing
      TIter         next(mixedEvents[mixBin]);
      TClonesArray *ma;
      while ((ma = (TClonesArray *)next())) {
         for (Int_t k = 0; k < ma->GetEntriesFast(); k++) {
            MpdPairPiKTrack *p2 = (MpdPairPiKTrack *)ma->At(k);
            double           m  = (mP1[i] + *p2).M();
            double           pt = (mP1[i] + *p2).Pt();

            if (fabs((mP1[i] + *p2).Rapidity()) < mParams.mYCut) {
               mInvMixNoPID->Fill(pt, m);
               mInvMixNoPIDBin[anaBin]->Fill(pt, m);
               if (mP1[i].getPID() == 5 || p2->getPID() == 5) {
                  mInvMixOnePID->Fill(pt, m);
                  mInvMixOnePIDBin[anaBin]->Fill(pt, m);
               }
               if (mP1[i].getPID() == 5 && p2->getPID() == 5) {
                  mInvMixTwoPID->Fill(pt, m);
                  mInvMixTwoPIDBin[anaBin]->Fill(pt, m);
               }
            } // rapidity
         }
      } // ma
   }    // nPos

   // Remove redandant events
   if (mixedEvents[mixBin]->GetSize() > nMixed) {
      TClonesArray *tmp = (TClonesArray *)mixedEvents[mixBin]->Last();
      mixedEvents[mixBin]->Remove(tmp);
      delete tmp;
   }
   mixedEvents[mixBin]->AddFirst(eventM);

} // histos

bool MpdPairPiK::selectTrack(MpdTrack *mpdtrack)
{
   if (mpdtrack->GetNofHits() < mParams.mNofHitsCut) return false; // nhits > 10

   if (fabs(mpdtrack->GetEta()) > mParams.mEtaCut) return false; // |eta| < 1.0

   float pt = fabs(mpdtrack->GetPt());
   if (pt < mParams.mPtminCut) return false; // pT > 50 MeV/c

   // if (fabs(mpdtrack->GetDCAX()) > mParams.mDCACut) return false; // |DCAx| < 2.
   // if (fabs(mpdtrack->GetDCAY()) > mParams.mDCACut) return false; // |DCAy| < 2.
   // if (fabs(mpdtrack->GetDCAZ()) > mParams.mDCACut) return false; // |DCAz| < 2.

   if (fabs(mpdtrack->GetNSigmaDCAx()) > mParams.mDCACut) return false; // |DCAx| < 2*sigmaxy(pT, eta, centrality).
   if (fabs(mpdtrack->GetNSigmaDCAy()) > mParams.mDCACut) return false; // |DCAy| < 2*sigmaxy(pT, eta, centrality).
   if (fabs(mpdtrack->GetNSigmaDCAz()) > mParams.mDCACut) return false; // |DCAz| < 2*sigmaz(pT, eta, centrality).

   return true;
}

long int MpdPairPiK::IsSameParent(long int prim1, long int prim2) const
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

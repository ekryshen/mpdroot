#include <iostream>
#include <fstream> // std::ifstream

#include "MpdKalmanFilter.h"
#include "MpdMCTrack.h"
#include "MpdTpcKalmanTrack.h"
#include "MpdEmcClusterKI.h"
#include "MpdParticle.h"
#include "MpdVertex.h"
#include "MpdEvent.h"
#include "MpdEmcGeoUtils.h"

#include "MpdPairKK.h"
#include "TFile.h"

ClassImp(MpdPairKK);

MpdPairKK::MpdPairKK(const char *name, const char *outputName) : MpdAnalysisTask(name, outputName)
{
   mParamConfig = outputName;
}

void MpdPairKK::UserInit()
{

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
   mhMultiplicity = new TH1F("hMultiplicity", "Multiplicity distribution", 2000, -0.5, 1999.5);
   fOutputList->Add(mhMultiplicity);

   // Minv
   mInvNoPID = new TH2F("mInvNoPID", "mInvNoPID", 100, 0., 10., 200, 0.9, 2.0);
   fOutputList->Add(mInvNoPID);
   mInvOnePID = new TH2F("mInvOnePID", "mInvOnePID", 100, 0., 10., 200, 0.9, 2.0);
   fOutputList->Add(mInvOnePID);
   mInvTwoPID = new TH2F("mInvTwoPID", "mInvTwoPID", 100, 0., 10., 200, 0.9, 2.0);
   fOutputList->Add(mInvTwoPID);

   // MinvMix
   mInvMixNoPID = new TH2F("mInvMixNoPID", "mInvMixNoPID", 100, 0., 10., 200, 0.9, 2.0);
   fOutputList->Add(mInvMixNoPID);
   mInvMixOnePID = new TH2F("mInvMixOnePID", "mInvMixOnePID", 100, 0., 10., 200, 0.9, 2.0);
   fOutputList->Add(mInvMixOnePID);
   mInvMixTwoPID = new TH2F("mInvMixTwoPID", "mInvMixTwoPID", 100, 0., 10., 200, 0.9, 2.0);
   fOutputList->Add(mInvMixTwoPID);

   // MC
   if (isMC) {
      mInvTrueNoPID = new TH2F("mInvTrueNoPID", "mInvTrueNoPID", 100, 0., 10., 200, 0.9, 2.0);
      fOutputList->Add(mInvTrueNoPID);
      mInvTrueNoPIDPhi = new TH2F("mInvTrueNoPIDPhi", "mInvTrueNoPIDPhi", 100, 0., 10., 200, 0.9, 2.0);
      fOutputList->Add(mInvTrueNoPIDPhi);
      mInvTrueOnePID = new TH2F("mInvTrueOnePID", "mInvTrueOnePID", 100, 0., 10., 200, 0.9, 2.0);
      fOutputList->Add(mInvTrueOnePID);
      mInvTrueOnePIDPhi = new TH2F("mInvTrueOnePIDPhi", "mInvTrueOnePIDPhi", 100, 0., 10., 200, 0.9, 2.0);
      fOutputList->Add(mInvTrueOnePIDPhi);
      mInvTrueTwoPID = new TH2F("mInvTrueTwoPID", "mInvTrueTwoPID", 100, 0., 10., 200, 0.9, 2.0);
      fOutputList->Add(mInvTrueTwoPID);
      mInvTrueTwoPIDPhi = new TH2F("mInvTrueTwoPIDPhi", "mInvTrueTwoPIDPhi", 100, 0., 10., 200, 0.9, 2.0);
      fOutputList->Add(mInvTrueTwoPIDPhi);
   }

   for (long int i = 0; i < nMixTot; i++) {
      mixedEvents[i] = new TList();
   }
}
//--------------------------------------
void MpdPairKK::ProcessEvent(MpdAnalysisEvent &event)
{

   if (!isInitialized) {
      mKF = MpdKalmanFilter::Instance();
      mKHit.SetType(MpdKalmanHit::kFixedR);
      isInitialized = true;
   }

   if (!selectEvent(event)) { //(V)
      return;
   }

   // mEMCClusters  = event.fEMCCluster;
   mKalmanTracks = event.fTPCKalmanTrack;
   if (isMC) {
      mMCTracks = event.fMCTrack;
   }

   selectPosTrack(event);

   selectNegTrack(event);

   processHistograms(event);
}

void MpdPairKK::Finish()
{
   // Post-scan processing not needed
}

//--------------------------------------
bool MpdPairKK::selectEvent(MpdAnalysisEvent &event)
{

   mhEvents->Fill(0.5);
   // first test if event filled?
   if (!event.fVertex) { // if even vertex not filled, skip event
      return false;
   }
   // Vertex z coordinate
   MpdVertex *vertex = (MpdVertex *)event.fVertex->First();
   vertex->Position(mPrimaryVertex);
   mhVertex->Fill(mPrimaryVertex.Z());
   if (fabs(mPrimaryVertex.Z()) > mParams.mZvtxCut) {
      return false;
   }
   mZvtxBin = 0.5 * (mPrimaryVertex.Z() / mParams.mZvtxCut + 1) * nMixEventZ;
   if (mZvtxBin < 0) mZvtxBin = 0;
   if (mZvtxBin >= nMixEventZ) mZvtxBin = nMixEventZ - 1;
   mhEvents->Fill(1.5);

   mMpdGlobalTracks = event.fMPDEvent->GetGlobalTracks();
   int ntr          = mMpdGlobalTracks->GetEntriesFast();

   mCenBin = (float(ntr) / 1200.) * nMixEventCent; // very rough
   if (mCenBin < 0) mCenBin = 0;
   if (mCenBin >= nMixEventCent) mCenBin = nMixEventCent - 1;

   // Multiplicity
   mhMultiplicity->Fill(ntr);

   // Centrality
   mhCentrality->Fill(mCenBin);
   // mCenBin = 0;
   mhEvents->Fill(2.5);

   // ZCD vs TPC (pileup?)
   mhEvents->Fill(3.5);

   // Eventplane  TODO
   mRPBin = 0;
   mhEvents->Fill(4.5);

   mixBin = (mCenBin + 1) * (mZvtxBin + 1) * (mRPBin + 1);
   // cout<<"Mixing bin: "<<mixBin<<" = "<<mCenBin<<" "<<mZvtxBin<<" "<<mRPBin<<endl;

   return true;
}
//--------------------------------------
void MpdPairKK::selectPosTrack(MpdAnalysisEvent &event)
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

      int charge;
      if (mpdtrack->GetPt() < 0)
         charge = 1;
      else
         charge = -1;
      if (charge < 0) continue;

      long int trId = -1;
      if (isMC) {
         trId = mpdtrack->GetID();
      }

      int   pid     = -1;
      int   isK_TPC = -1;
      int   isK_TOF = -1;
      float pmom    = sqrt(mpdtrack->GetPt() * mpdtrack->GetPt() + mpdtrack->GetPz() * mpdtrack->GetPz());

      if (fabs(dEdx_sigma_K(tr->GetDedx(), pmom)) < mParams.mPIDsigTPC) isK_TPC = 1;
      if ((mpdtrack->GetTofFlag() == 2 || mpdtrack->GetTofFlag() == 6) &&
          fabs(Beta_sigma_K(mpdtrack->GetTofBeta(), pmom)) < mParams.mPIDsigTOF)
         isK_TOF = 1;

      if (isK_TPC == 1) pid = 1;
      if (isK_TOF == 1) pid = 2;
      if (isK_TPC == 1 || isK_TOF == 1) pid = 3;
      if (isK_TPC == 1 && isK_TOF == 1) pid = 4;

      // if ( trId > -1 && fabs(static_cast<MpdMCTrack *>(mMCTracks->At(trId))->GetPdgCode()) == 321)
      // 	{
      // 	  cout<<"Pos: "<<pmom<<" dedx: "<<dEdx_sigma_K(tr->GetDedx(), pmom)<<"  tof: "<<mpdtrack->GetTofFlag()<<"
      // beta: "<<Beta_sigma_K(mpdtrack->GetTofBeta(), pmom)<<" isis: "<<isK_TPC<<" -  "<<isK_TOF<<"   PIDs:
      // "<<pid<<endl;
      // 	}

      if (pid < 0) continue;

      float mK = 0.493677;
      mP1.emplace_back(mpdtrack->GetPx(), mpdtrack->GetPy(), mpdtrack->GetPz(), sqrt(pmom * pmom + mK * mK));
      mP1.back().setPrimary(trId);
      mP1.back().setTr1(i);
      mP1.back().setPID(pid);
   } // i
}
//--------------------------------------
void MpdPairKK::selectNegTrack(MpdAnalysisEvent &event)
{
   // h-

   eventM = new TClonesArray("MpdPairKKTrack", mMpdGlobalTracks->GetEntriesFast());
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

      int charge;
      if (mpdtrack->GetPt() < 0)
         charge = 1;
      else
         charge = -1;
      if (charge > 0) continue;

      long int trId = -1;
      if (isMC) {
         trId = mpdtrack->GetID();
      }

      int   pid     = -1;
      int   isK_TPC = -1;
      int   isK_TOF = -1;
      float pmom    = sqrt(mpdtrack->GetPt() * mpdtrack->GetPt() + mpdtrack->GetPz() * mpdtrack->GetPz());

      if (fabs(dEdx_sigma_K(tr->GetDedx(), pmom)) < mParams.mPIDsigTPC) isK_TPC = 1;
      if ((mpdtrack->GetTofFlag() == 2 || mpdtrack->GetTofFlag() == 6) &&
          fabs(Beta_sigma_K(mpdtrack->GetTofBeta(), pmom)) < mParams.mPIDsigTOF)
         isK_TOF = 1;

      if (isK_TPC == 1) pid = 1;
      if (isK_TOF == 1) pid = 2;
      if (isK_TPC == 1 || isK_TOF == 1) pid = 3;
      if (isK_TPC == 1 && isK_TOF == 1) pid = 4;

      // if ( trId > -1 && fabs(static_cast<MpdMCTrack *>(mMCTracks->At(trId))->GetPdgCode()) == 321)
      // 	{
      // 	  cout<<"NEG: "<<pmom<<" dedx: "<<dEdx_sigma_K(tr->GetDedx(), pmom)<<"  tof: "<<mpdtrack->GetTofFlag()<<"
      // beta: "<<Beta_sigma_K(mpdtrack->GetTofBeta(), pmom)<<" isis: "<<isK_TPC<<" -  "<<isK_TOF<<"   PIDs:
      // "<<pid<<endl;
      // 	}

      if (pid < 0) continue;

      float mK = 0.493677;
      mP2.emplace_back(mpdtrack->GetPx(), mpdtrack->GetPy(), mpdtrack->GetPz(), sqrt(pmom * pmom + mK * mK));
      mP2.back().setPrimary(trId);
      mP2.back().setTr1(i);
      mP2.back().setPID(pid);

      new ((*eventM)[iv++]) MpdPairKKTrack(mP2[mP2.size() - 1]);
   } // i
}

void MpdPairKK::processHistograms(MpdAnalysisEvent &event)
{
   // Fill Real, Mixed distributions and update mixed array

   // Real
   int nPos = mP1.size();
   int nNeg = mP2.size();

   for (int i = 0; i < nPos; i++) {
      for (int j = 0; j < nNeg; j++) {
         TLorentzVector sum = mP1[i] + mP2[j];
         mInvNoPID->Fill(sum.Pt(), sum.M());
         if (mP1[i].getPID() == 4 || mP2[j].getPID() == 4) mInvOnePID->Fill(sum.Pt(), sum.M());
         if (mP1[i].getPID() == 4 && mP2[j].getPID() == 4) mInvTwoPID->Fill(sum.Pt(), sum.M());

         long int ip = IsSameParent(mP1[i].primary(), mP2[j].primary());
         if (ip >= 0) {

            mInvTrueNoPID->Fill(sum.Pt(), sum.M());
            if (mP1[i].getPID() == 4 || mP2[j].getPID() == 4) mInvTrueOnePID->Fill(sum.Pt(), sum.M());
            if (mP1[i].getPID() == 4 && mP2[j].getPID() == 4) mInvTrueTwoPID->Fill(sum.Pt(), sum.M());

            if (static_cast<MpdMCTrack *>(mMCTracks->At(ip))->GetPdgCode() == 333) {
               mInvTrueNoPIDPhi->Fill(sum.Pt(), sum.M());
               if (mP1[i].getPID() == 4 || mP2[j].getPID() == 4) mInvTrueOnePIDPhi->Fill(sum.Pt(), sum.M());
               if (mP1[i].getPID() == 4 && mP2[j].getPID() == 4) mInvTrueTwoPIDPhi->Fill(sum.Pt(), sum.M());
            }

         } // ip
      }    // nNeg

      // Mixing
      TIter         next(mixedEvents[mixBin]);
      TClonesArray *ma;
      while ((ma = (TClonesArray *)next())) {
         for (Int_t k = 0; k < ma->GetEntriesFast(); k++) {
            MpdPairKKTrack *p2 = (MpdPairKKTrack *)ma->At(k);
            double          m  = (mP1[i] + *p2).M();
            double          pt = (mP1[i] + *p2).Pt();
            mInvMixNoPID->Fill(pt, m);
            if (mP1[i].getPID() == 4 || p2->getPID() == 4) mInvMixOnePID->Fill(pt, m);
            if (mP1[i].getPID() == 4 && p2->getPID() == 4) mInvMixTwoPID->Fill(pt, m);
         }
      } // ma

   } // nPos

   // Remove redandant events
   if (mixedEvents[mixBin]->GetSize() > nMixed) {
      TClonesArray *tmp = (TClonesArray *)mixedEvents[mixBin]->Last();
      mixedEvents[mixBin]->Remove(tmp);
      delete tmp;
   }
   mixedEvents[mixBin]->AddFirst(eventM);

} // histos

bool MpdPairKK::selectTrack(MpdTrack *mpdtrack)
{
   if (mpdtrack->GetNofHits() < mParams.mNofHitsCut) return false; // nhits > 10

   if (fabs(mpdtrack->GetEta()) > mParams.mEtaCut) return false; //|eta| < 1.0

   float pt = TMath::Abs(mpdtrack->GetPt());
   if (pt < mParams.mPtminCut) return false; // pT > 50 MeV/c

   return true;
}

long int MpdPairKK::IsSameParent(long int prim1, long int prim2) const
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

//====================
float MpdPairKK::dEdx_sigma_K(float dEdx, float mom) const
{
   if (dEdx < 0) return -999;

   dEdx = log10(dEdx);

   if (mom < 0.075) mom = 0.075;
   if (mom > 2.5) mom = 2.5;

   float mean[7]  = {1.112458e-001, 1.366928e-002,  9.475245e+000, -2.230105e+000,
                    1.430526e-002, -7.711245e-001, 8.563292e-001};
   float width[7] = {2.883180e-003,  -2.542899e-002, 4.469939e+000, 7.136557e-001,
                     -1.194883e-001, 9.797266e-001,  5.707477e-002};

   float mean_exp, width_exp;

   mean_exp =
      mean[0] / mom / mom *
         (mean[1] * log(mom * mom) - mean[2] * mom * mom - mean[3] * mom - mean[4] - mean[5] * mom * mom * mom) +
      mean[6];
   width_exp =
      width[0] / mom / mom *
         (width[1] * log(mom * mom) - width[2] * mom * mom - width[3] * mom - width[4] - width[5] * mom * mom * mom) +
      width[6];

   return (dEdx - mean_exp) / width_exp;
}

float MpdPairKK::Beta_sigma_K(float beta, float mom) const
{
   if (mom < 0.2) return -999;

   float mean[7]  = {3.150000e+003,  1.086959e-007, 1.490652e+000, 9.061339e-005,
                    -7.965485e-006, 1.603893e-005, 4.696765e+003};
   float width[5] = {5.616303e-003, 5.634408e-003, -2.888325e-002, 1.108509e-006, -3.503470e-002};

   float mean_exp, width_exp;

   mean_exp =
      mean[0] / mom / mom *
         (mean[1] * log(mom * mom) - mean[2] * mom * mom - mean[3] * mom - mean[4] - mean[5] * mom * mom * mom) +
      mean[6] - 0.001;
   width_exp = width[0] + width[1] * pow(mom, width[2]) + width[3] / pow(mom - width[4], 6);

   return (beta - mean_exp) / width_exp;
}

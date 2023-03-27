#include <iostream>
#include <fstream> // std::ifstream

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
  cout << "[MpdPairKK]: Initialization ... " << endl;

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

   // Minv
   mInvNoPID = new TH2F("mInvNoPID", "mInvNoPID", 100, 0., 10., 200, 0.9, 2.0);
   fOutputList->Add(mInvNoPID);
   mInvOnePID = new TH2F("mInvOnePID", "mInvOnePID", 100, 0., 10., 200, 0.9, 2.0);
   fOutputList->Add(mInvOnePID);
   mInvTwoPID = new TH2F("mInvTwoPID", "mInvTwoPID", 100, 0., 10., 200, 0.9, 2.0);
   fOutputList->Add(mInvTwoPID);
   mAccEffRecNoPID = new TH2F("mAccEffRecNoPID", "mAccEffRecNoPID", 100, -4, 4, 100, 0., 10.);
   fOutputList->Add(mAccEffRecNoPID);
   mAccEffRecOnePID = new TH2F("mAccEffRecOnePID", "mAccEffRecOnePID", 100, -4, 4, 100, 0., 10.);
   fOutputList->Add(mAccEffRecOnePID);
   mAccEffRecTwoPID = new TH2F("mAccEffRecTwoPID", "mAccEffRecTwoPID", 100, -4, 4, 100, 0., 10.);
   fOutputList->Add(mAccEffRecTwoPID);

   // MinvMix
   mInvMixNoPID = new TH2F("mInvMixNoPID", "mInvMixNoPID", 100, 0., 10., 200, 0.9, 2.0);
   fOutputList->Add(mInvMixNoPID);
   mInvMixOnePID = new TH2F("mInvMixOnePID", "mInvMixOnePID", 100, 0., 10., 200, 0.9, 2.0);
   fOutputList->Add(mInvMixOnePID);
   mInvMixTwoPID = new TH2F("mInvMixTwoPID", "mInvMixTwoPID", 100, 0., 10., 200, 0.9, 2.0);
   fOutputList->Add(mInvMixTwoPID);

    for (int bin = 0; bin < nCenBinsAna; bin++) {
      mInvNoPIDBin[bin] = 
         new TH2F(Form("mInvNoPIDBin_%d", bin), Form("mInvNoPIDBin_%d", bin), 100, 0., 10., 200, 0.9, 2.0);
      fOutputList->Add(mInvNoPIDBin[bin]);
      mInvOnePIDBin[bin] = 
         new TH2F(Form("mInvOnePIDBin_%d", bin), Form("mInvOnePIDBin_%d", bin), 100, 0., 10., 200, 0.9, 2.0);
      fOutputList->Add(mInvOnePIDBin[bin]);
      mInvTwoPIDBin[bin] = 
         new TH2F(Form("mInvTwoPIDBin_%d", bin), Form("mInvTwoPIDBin_%d", bin), 100, 0., 10., 200, 0.9, 2.0);
      fOutputList->Add(mInvTwoPIDBin[bin]);

      mInvMixNoPIDBin[bin] = 
         new TH2F(Form("mInvMixNoPIDBin_%d", bin), Form("mInvMixNoPIDBin_%d", bin), 100, 0., 10., 200, 0.9, 2.0);
      fOutputList->Add(mInvMixNoPIDBin[bin]);
      mInvMixOnePIDBin[bin] = 
         new TH2F(Form("mInvMixOnePIDBin_%d", bin), Form("mInvMixOnePIDBin_%d", bin), 100, 0., 10., 200, 0.9, 2.0);
      fOutputList->Add(mInvMixOnePIDBin[bin]);
      mInvMixTwoPIDBin[bin] = 
         new TH2F(Form("mInvMixTwoPIDBin_%d", bin), Form("mInvMixTwoPIDBin_%d", bin), 100, 0., 10., 200, 0.9, 2.0);
      fOutputList->Add(mInvMixTwoPIDBin[bin]);
    }

   // MC
   if (isMC) {

     // Generated signals
     mInvGen = new TH1F("mInvGen", "mInvGen", 100, 0., 10.);
     fOutputList->Add(mInvGen);
     
     mInvTrueNoPID = new TH2F("mInvTrueNoPID", "mInvTrueNoPID", 100, 0., 10., 200, 0.9, 2.0);
     fOutputList->Add(mInvTrueNoPID);
     mInvTrueNoPIDPhi = new TH2F("mInvTrueNoPIDPhi", "mInvTrueNoPIDPhi", 100, 0., 10., 200, 0.9, 2.0);
     fOutputList->Add(mInvTrueNoPIDPhi);
     mMassResNoPIDPhi = new TH2F("mMassResNoPIDPhi", "mMassResNoPIDPhi", 100, 0., 10., 1000, -0.5, 0.5);
     fOutputList->Add(mMassResNoPIDPhi);
     mInvTrueOnePID = new TH2F("mInvTrueOnePID", "mInvTrueOnePID", 100, 0., 10., 200, 0.9, 2.0);
     fOutputList->Add(mInvTrueOnePID);
     mInvTrueOnePIDPhi = new TH2F("mInvTrueOnePIDPhi", "mInvTrueOnePIDPhi", 100, 0., 10., 200, 0.9, 2.0);
     fOutputList->Add(mInvTrueOnePIDPhi);
     mMassResOnePIDPhi = new TH2F("mMassResOnePIDPhi", "mMassResOnePIDPhi", 100, 0., 10., 1000, -0.5, 0.5);
     fOutputList->Add(mMassResOnePIDPhi);
     mInvTrueTwoPID = new TH2F("mInvTrueTwoPID", "mInvTrueTwoPID", 100, 0., 10., 200, 0.9, 2.0);
     fOutputList->Add(mInvTrueTwoPID);
     mInvTrueTwoPIDPhi = new TH2F("mInvTrueTwoPIDPhi", "mInvTrueTwoPIDPhi", 100, 0., 10., 200, 0.9, 2.0);
     fOutputList->Add(mInvTrueTwoPIDPhi);
     mMassResTwoPIDPhi = new TH2F("mMassResTwoPIDPhi", "mMassResTwoPIDPhi", 100, 0., 10., 1000, -0.5, 0.5);
     fOutputList->Add(mMassResTwoPIDPhi);
     mAccEffGen = new TH2F("mAccEffGen", "mAccEffGen", 100, -4, 4, 100, 0., 10.);
     fOutputList->Add(mAccEffGen);

     for (int bin = 0; bin < nCenBinsAna; bin++) {
       mInvGenBin[bin] = new TH1F(Form("mInvGenBin_%d", bin), Form("mInvGenBin_%d", bin), 100, 0., 10.);
       fOutputList->Add(mInvGenBin[bin]);
       mInvTrueNoPIDPhiBin[bin] = 
          new TH2F(Form("mInvTrueNoPIDPhiBin_%d", bin), Form("mInvTrueNoPIDPhiBin_%d", bin), 100, 0., 10., 200, 0.9, 2.0);
       fOutputList->Add(mInvTrueNoPIDPhiBin[bin]);
       mInvTrueOnePIDPhiBin[bin] = 
          new TH2F(Form("mInvTrueOnePIDPhiBin_%d", bin), Form("mInvTrueOnePIDPhiBin_%d", bin), 100, 0., 10., 200, 0.9, 2.0);
       fOutputList->Add(mInvTrueOnePIDPhiBin[bin]);
       mInvTrueTwoPIDPhiBin[bin] = 
          new TH2F(Form("mInvTrueTwoPIDPhiBin_%d", bin), Form("mInvTrueTwoPIDPhiBin_%d", bin), 100, 0., 10., 200, 0.9, 2.0);
       fOutputList->Add(mInvTrueTwoPIDPhiBin[bin]);
     }

   }

   for (long int i = 0; i < nMixTot; i++) {
      mixedEvents[i] = new TList();
   }

  cout << "[MpdPairKK]: Reading Geo for track refit ... " << endl;

   // Read-out TPC Geo for track Refit
  inFileSim = new TFile("sim_Geo.root","READ");
  inTreeSim = (TTree*) inFileSim->Get("mpdsim");
  inTreeSim->SetName("mpdsim1");

  inFileSim->Get("FairGeoParSet");

  tpcPoints = (TClonesArray*) inFileSim->FindObjectAny("TpcPoint");

  inTreeSim->SetBranchAddress("TpcPoint",&tpcPoints);
  TBranch *tpcSimB = inTreeSim->GetBranch("TpcPoint");

  secGeo = new TpcSectorGeoAZ();
  recoTpc = new MpdTpcKalmanFilter(*secGeo,"TPC Kalman filter");
  recoTpc->SetSectorGeo(*secGeo);
  recoTpc->FillGeoScheme();

  cout << "[MpdPairKK]: Reading DCA parameterizations ... " << endl;

  // Read-out DCA parameterization
  dcaFile = new TFile("DCAs.root","READ");

  for (Int_t etab=0; etab < neta_bins; etab++)
    {
      for (Int_t centb=0; centb < ncent_bins; centb++)
	{
	  f_dca_xy[etab][centb] = (TF1*) dcaFile->Get(Form("dcaxy_fitf_eta_%d_cent_%d",etab,centb));
	  f_dca_z[etab][centb]  = (TF1*) dcaFile->Get(Form("dcaz_fitf_eta_%d_cent_%d",etab,centb));
	}
    }        

  cout << "[MpdPairKK]: Initialization done " << endl << endl;
}
//--------------------------------------
void MpdPairKK::ProcessEvent(MpdAnalysisEvent &event)
{

   if (!isInitialized) {
     // mKF = MpdKalmanFilter::Instance();
     // mKHit.SetType(MpdKalmanHit::kFixedR);
     isInitialized = true;
   }

   if (!selectEvent(event)) { //(V)
      return;
   }

   mKalmanTracks = event.fTPCKalmanTrack;
   mpdTofMatching = event.fTOFMatching;

   if (isMC) {
      mMCTracks = event.fMCTrack;

      for (int i = 0; i < mMCTracks->GetEntriesFast(); i++) {
         MpdMCTrack *pr = (static_cast<MpdMCTrack *>(mMCTracks->At(i)));
         if (pr->GetPdgCode() == 333) {
	   TVector3 momentum;
	   pr->GetMomentum(momentum);
	   mAccEffGen->Fill(pr->GetRapidity(), momentum.Pt());
	   if (fabs(pr->GetRapidity()) < mParams.mYCut) {
	     if (pr->GetStartX() * pr->GetStartX() + pr->GetStartY() * pr->GetStartY() < 1.) {
               mInvGen->Fill(momentum.Pt());
               mInvGenBin[anaBin]->Fill(momentum.Pt());
	     } // radius
	   } // rapidity
         } // ID
      }
   } // isMC

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

   mCenBin   = (cen / 100.) * nMixEventCent; // very rough
   if (mCenBin < 0) mCenBin = 0;
   if (mCenBin >= nMixEventCent) mCenBin = nMixEventCent - 1;

   mixBin = mZvtxBin * nMixEventCent + mCenBin;

   mhVertex->Fill(mPrimaryVertex.Z());
   mhCentrality->Fill(mCenBin);

   anaBin = -1;
   if (cen >= 0  && cen < 10) anaBin = 0;
   if (cen >= 10 && cen < 20) anaBin = 1;
   if (cen >= 20 && cen < 30) anaBin = 2;
   if (cen >= 30 && cen < 40) anaBin = 3;
   if (cen >= 40 && cen < 50) anaBin = 4;
   if (cen >= 50 && cen < 60) anaBin = 5;
   if (cen >= 60 && cen < 100) anaBin = 6;

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

      int   isTOF    = -1;
      int   pid      = -1;
      int   isK_TPC  = -1;
      int   isK_TOF  = -1;
      float pmom     = sqrt(mpdtrack->GetPt() * mpdtrack->GetPt() + mpdtrack->GetPz() * mpdtrack->GetPz());
      
      if ( mpdtrack->GetTofFlag() == 2 || mpdtrack->GetTofFlag() == 6 ) {
	int matchingIndex = -1;
	      
	for (Int_t l=0; l < mpdTofMatching->GetEntries(); l++)
	  {
	    MpdTofMatchingData* Matching = (MpdTofMatchingData*) mpdTofMatching->At(l);
	    if (Matching->GetKFTrackIndex() == i) 
	      {
		matchingIndex=l; 
		break;
	      }
	  }
	      
	if (matchingIndex > 0)
	  {
	    MpdTofMatchingData* Matching = (MpdTofMatchingData*) mpdTofMatching->At(matchingIndex);
	    if(TestTofMatch(charge, fabs(mpdtrack->GetPt()), Matching->GetdPhi(), Matching->GetdZed()) == 1) isTOF = 1;
	  }	      
      }
   
      if (fabs(dEdx_sigma_K(tr->GetDedx(), pmom)) < mParams.mPIDsigTPC) isK_TPC = 1;
      if (isTOF == 1 && fabs(Beta_sigma_K(mpdtrack->GetTofBeta(), pmom)) < mParams.mPIDsigTOF) isK_TOF = 1;
      
      if (isK_TPC == 1) pid = 1;
      if (isK_TOF == 1) pid = 2;
      if (isK_TPC == 1 || isK_TOF == 1) pid = 3;
      if (isK_TPC == 1 && isK_TOF == 1) pid = 4;
      if ( (isTOF != 1 && isK_TPC == 1) || (isTOF == 1 && isK_TPC == 1 && isK_TOF == 1) ) pid = 5;
      
      if (pid < 0) continue;

      float mK = 0.493677;

      MpdTpcKalmanTrack trCorK = *tr;
      trCorK.SetDirection(MpdKalmanTrack::kInward);
      int ok = recoTpc->Refit(&trCorK, mK,  1); // refit
      if (ok) {
	MpdParticle kaon(trCorK, i);
	kaon.SetPdg(321);
	kaon.SetMass();
	
	pmom    = sqrt( kaon.Momentum3().Pt() * kaon.Momentum3().Pt() + kaon.Momentum3().Pz() * kaon.Momentum3().Pz());
	
	mP1.emplace_back(kaon.Momentum3().Px(), kaon.Momentum3().Py(), kaon.Momentum3().Pz(), sqrt(pmom * pmom + mK * mK));
      }
      else {
	mP1.emplace_back(mpdtrack->GetPx(), mpdtrack->GetPy(), mpdtrack->GetPz(), sqrt(pmom * pmom + mK * mK));
      }

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

      int   isTOF    = -1;
      int   pid      = -1;
      int   isK_TPC  = -1;
      int   isK_TOF  = -1;
      float pmom     = sqrt(mpdtrack->GetPt() * mpdtrack->GetPt() + mpdtrack->GetPz() * mpdtrack->GetPz());

      if ( mpdtrack->GetTofFlag() == 2 || mpdtrack->GetTofFlag() == 6 ) {
	int matchingIndex = -1;
	      
	for (Int_t l=0; l < mpdTofMatching->GetEntries(); l++)
	  {
	    MpdTofMatchingData* Matching = (MpdTofMatchingData*) mpdTofMatching->At(l);
	    if (Matching->GetKFTrackIndex() == i) 
	      {
		matchingIndex=l; 
		break;
	      }
	  }
	      
	if (matchingIndex > 0)
	  {
	    MpdTofMatchingData* Matching = (MpdTofMatchingData*) mpdTofMatching->At(matchingIndex);
	    if(TestTofMatch(charge, fabs(mpdtrack->GetPt()), Matching->GetdPhi(), Matching->GetdZed()) == 1) isTOF = 1;
	  }	      
      }

      if (fabs(dEdx_sigma_K(tr->GetDedx(), pmom)) < mParams.mPIDsigTPC) isK_TPC = 1;
      if (isTOF == 1 && fabs(Beta_sigma_K(mpdtrack->GetTofBeta(), pmom)) < mParams.mPIDsigTOF) isK_TOF = 1;

      if (isK_TPC == 1) pid = 1;
      if (isK_TOF == 1) pid = 2;
      if (isK_TPC == 1 || isK_TOF == 1) pid = 3;
      if (isK_TPC == 1 && isK_TOF == 1) pid = 4;
      if ( (isTOF != 1 && isK_TPC == 1) || (isTOF == 1 && isK_TPC == 1 && isK_TOF == 1) ) pid = 5;

      if (pid < 0) continue;

      float mK = 0.493677;

      MpdTpcKalmanTrack trCorK = *tr;
      trCorK.SetDirection(MpdKalmanTrack::kInward);
      int ok = recoTpc->Refit(&trCorK, mK,  1); // refit, sign seems to be not important ???
      if (ok) {
	MpdParticle kaon(trCorK, i);
	kaon.SetPdg(-321);
	kaon.SetMass();
	
	pmom    = sqrt( kaon.Momentum3().Pt() * kaon.Momentum3().Pt() + kaon.Momentum3().Pz() * kaon.Momentum3().Pz());
	
	mP2.emplace_back(kaon.Momentum3().Px(), kaon.Momentum3().Py(), kaon.Momentum3().Pz(), sqrt(pmom * pmom + mK * mK));
      }
      else {
	mP2.emplace_back(mpdtrack->GetPx(), mpdtrack->GetPy(), mpdtrack->GetPz(), sqrt(pmom * pmom + mK * mK));
      }

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
	 } //rapidity

	 long int ip = IsSameParent(mP1[i].primary(), mP2[j].primary());
	 if (ip >= 0) {
	   
	   if (fabs(sum.Rapidity()) < mParams.mYCut) {
	     mInvTrueNoPID->Fill(sum.Pt(), sum.M());
	     if (mP1[i].getPID() == 5 || mP2[j].getPID() == 5) mInvTrueOnePID->Fill(sum.Pt(), sum.M());
	     if (mP1[i].getPID() == 5 && mP2[j].getPID() == 5) mInvTrueTwoPID->Fill(sum.Pt(), sum.M());
	   } // rapidity

	   if (static_cast<MpdMCTrack *>(mMCTracks->At(ip))->GetPdgCode() == 333) {
	       
	     MpdMCTrack *pr1 = (static_cast<MpdMCTrack *>(mMCTracks->At(mP1[i].primary())));
	     MpdMCTrack *pr2 = (static_cast<MpdMCTrack *>(mMCTracks->At(mP2[j].primary())));
	     
	     TLorentzVector p1, p2, sump1p2;
	     p1.SetPxPyPzE(pr1->GetPx(), pr1->GetPy(), pr1->GetPz(), pr1->GetEnergy());
	     p2.SetPxPyPzE(pr2->GetPx(), pr2->GetPy(), pr2->GetPz(), pr2->GetEnergy());
	     
	     sump1p2 = p1 + p2;

	     mAccEffRecNoPID->Fill(sum.Rapidity(), sum.Pt());
	     if (mP1[i].getPID() == 5 || mP2[j].getPID() == 5) mAccEffRecOnePID->Fill(sum.Rapidity(), sum.Pt());
	     if (mP1[i].getPID() == 5 && mP2[j].getPID() == 5) mAccEffRecTwoPID->Fill(sum.Rapidity(), sum.Pt());

	     if (fabs(sum.Rapidity()) < mParams.mYCut) {
	       mInvTrueNoPIDPhi->Fill(sum.Pt(), sum.M());
	       mMassResNoPIDPhi->Fill(sum.Pt(), sum.M() - sump1p2.M());
	       mInvTrueNoPIDPhiBin[anaBin]->Fill(sum.Pt(), sum.M());
	       
	       if (mP1[i].getPID() == 5 || mP2[j].getPID() == 5) {
		 mInvTrueOnePIDPhi->Fill(sum.Pt(), sum.M());
		 mMassResOnePIDPhi->Fill(sum.Pt(), sum.M() - sump1p2.M());
		 mInvTrueOnePIDPhiBin[anaBin]->Fill(sum.Pt(), sum.M());
	       }	    
	       
	       if (mP1[i].getPID() == 5 && mP2[j].getPID() == 5) {
		 mInvTrueTwoPIDPhi->Fill(sum.Pt(), sum.M());
		 mMassResTwoPIDPhi->Fill(sum.Pt(), sum.M() - sump1p2.M());
		 mInvTrueTwoPIDPhiBin[anaBin]->Fill(sum.Pt(), sum.M());
	       }
	     } // rapidity
	   } // ID
	 } // ip
	 
      } // Neg

      // Mixing
      TIter         next(mixedEvents[mixBin]);
      TClonesArray *ma;
      while ((ma = (TClonesArray *)next())) {
         for (Int_t k = 0; k < ma->GetEntriesFast(); k++) {
            MpdPairKKTrack *p2 = (MpdPairKKTrack *)ma->At(k);
            double          m  = (mP1[i] + *p2).M();
            double          pt = (mP1[i] + *p2).Pt();

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

   if (fabs(mpdtrack->GetEta()) > mParams.mEtaCut) return false; // |eta| < 1.0

   float pt = fabs(mpdtrack->GetPt());
   if (pt < mParams.mPtminCut) return false; // pT > 50 MeV/c

   // if (fabs(mpdtrack->GetDCAX()) > mParams.mDCACut) return false; // |DCAx| < 2.
   // if (fabs(mpdtrack->GetDCAY()) > mParams.mDCACut) return false; // |DCAy| < 2.
   // if (fabs(mpdtrack->GetDCAZ()) > mParams.mDCACut) return false; // |DCAz| < 2.

   int cent_bin = -1;
   cent_bin = int(cen /(cent_max - cent_min) * float(ncent_bins));
   if (cent_bin == 9) cent_bin = 8; // 90-91% -> use parameterization for 80-90%
   if (cent_bin < 0 || cent_bin > ncent_bins-1) return false;

   int eta_bin = -1;	  
   eta_bin = int((mpdtrack->GetEta() + 1.5)/(eta_max - eta_min) * float(neta_bins));
   if (eta_bin < 0 || eta_bin > neta_bins-1) return false;
   
   TF1 sigma_fit_XY = *f_dca_xy[eta_bin][cent_bin];
   TF1 sigma_fit_Z = *f_dca_z[eta_bin][cent_bin];

   // use lower/upper limits outside of parameterization ranges
   float pt_dca = pt;
   if (pt_dca < 0.05) pt_dca = 0.05;
   double xfirst, xlast;
   sigma_fit_Z.GetRange(xfirst, xlast);
   if (pt_dca > xlast) pt_dca = xlast + 0.05;

   float sigma_exp_xy = sigma_fit_XY(pt_dca);
   float sigma_exp_z = sigma_fit_Z(pt_dca);

   float dcax_sig, dcay_sig, dcaz_sig;

   if (sigma_exp_xy != 0) { 
     dcax_sig = mpdtrack->GetDCAX()/sigma_exp_xy;
     dcay_sig = mpdtrack->GetDCAY()/sigma_exp_xy; 
   } else {dcax_sig = -999; dcay_sig = -999;}

   if (sigma_exp_z != 0) {
     dcaz_sig = mpdtrack->GetDCAZ()/sigma_exp_z;
   } else {dcaz_sig = -999;}

   if (fabs(dcax_sig) > mParams.mDCACut) return false; // |DCAx| < 2*sigmaxy(pT, eta, centrality).
   if (fabs(dcay_sig) > mParams.mDCACut) return false; // |DCAy| < 2*sigmaxy(pT, eta, centrality).
   if (fabs(dcaz_sig) > mParams.mDCACut) return false; // |DCAz| < 2*sigmaz(pT, eta, centrality).

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

// dE/dx parameterizations for el/Pi/K/p
float MpdPairKK::dEdx_sigma_El(float dEdx, float mom) const
{
  if (mom < 0.04) return -999;
  if (dEdx <= 0) return -999;
  
  dEdx = log(dEdx);

   if (mom < 0.06) mom = 0.06;
   if (mom > 2.0)  mom = 2.0;

   float mean[7] = { 6.869878e-003, -1.573146e-001, 1.847371e+000, -9.253678e-001, 1.073907e+000, -3.394239e+000, 6.451762e-001 };
   float width[7] = { -4.359368e+006, -8.508504e-012, -3.958364e-009, 1.526816e-009, 1.353776e-011, 3.426352e-009, 6.591542e-002 };

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

float MpdPairKK::dEdx_sigma_Pi(float dEdx, float mom) const
{
  if (mom < 0.06) return -999;
  if (dEdx <= 0) return -999;

   dEdx = log(dEdx);

   if (mom < 0.07) mom = 0.07;
   if (mom > 3.5)  mom = 3.5;

   float mean[7] = { -3.554770e-002, -2.666590e-001, 5.701508e+000, -5.220967e+000, 1.987420e+000, 1.087186e+000, 1.037435e-001 };
   float width[7] = { -3.495341e+007, 2.472950e-011, -1.422238e-010, 3.234257e-010, -1.242149e-010, -5.320861e-011, 7.771059e-002 };

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

float MpdPairKK::dEdx_sigma_K(float dEdx, float mom) const
{
  if (mom < 0.08) return -999;
  if (dEdx <= 0) return -999;

   dEdx = log(dEdx);

   if (mom < 0.09) mom = 0.09;
   if (mom > 3.5)  mom = 3.5;

   float mean[7]  = { -4.571263e-004, -1.255796e+001, -7.145174e+002, 1.023127e+003, 3.497842e+001, 2.890434e+002, -3.740436e-003 };
   float width[7] = { -9.820512e+006, -1.117400e-010, 1.188118e-009, -1.754282e-009, 7.875521e-010, -4.171753e-010, 8.100840e-002 };

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

float MpdPairKK::dEdx_sigma_P(float dEdx, float mom) const
{
  if (mom < 0.08) return -999;
  if (dEdx <= 0) return -999;

   dEdx = log(dEdx);

   if (mom < 0.09) mom = 0.09;
   if (mom > 3.5)  mom = 3.5;

   float mean[7] = { 1.183672e-001, -3.394823e-001, 2.418100e+001, -1.377881e+001, 2.534160e+000, -1.563054e+000, 2.010767e+000 };
   float width[7] = { -1.488536e+007, -2.075514e-010, 2.418077e-009, -2.996053e-009, 1.397806e-009, -4.161277e-010, 7.174217e-002 };

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

// Beta parameterizations for el/Pi/K/p
int MpdPairKK::TestTofMatch(int charge, float pt, float dphi, float dz) const
{
  if (pt < 0.1) return 0;

  float meanphiCoeff[7] = { 2.745900e-001, -1.437619e+000, 2.824655e+000, -2.644983e+000, 1.269623e+000, -3.012448e-001, 2.816063e-002 };
  float widthphiCoeff[7] = { 3.211516e+000, -2.325625e+001, 8.024492e+001, -1.407389e+002, 1.320427e+002, -6.304419e+001, 1.202106e+001 };
  float widthzCoeff[7] = { 1.646438e+000, -9.744003e+000, 3.054315e+001, -4.799731e+001, 3.994061e+001, -1.680153e+001, 2.810887e+000 };

  float meanphi, widthphi, meanz, widthz;

  if (pt < 3.0)
    {
      meanphi = meanphiCoeff[0] + 
	        meanphiCoeff[1] * pt +
	        meanphiCoeff[2] * pt*pt +
	        meanphiCoeff[3] * pt*pt*pt +
	        meanphiCoeff[4] * pt*pt*pt*pt +
	        meanphiCoeff[5] * pt*pt*pt*pt*pt +
	        meanphiCoeff[6] * pt*pt*pt*pt*pt*pt;
    }
  else
    {
      meanphi = 0.135;
    }

  if (charge < 0) meanphi = -meanphi;

  if (pt < 1.325)
    {
      widthphi = widthphiCoeff[0] + 
	         widthphiCoeff[1] * pt +
	         widthphiCoeff[2] * pt*pt +
	         widthphiCoeff[3] * pt*pt*pt +
	         widthphiCoeff[4] * pt*pt*pt*pt +
	         widthphiCoeff[5] * pt*pt*pt*pt*pt +
	         widthphiCoeff[6] * pt*pt*pt*pt*pt*pt;
    }
  else
    {
      widthphi = 0.454 + 0.00815163*pt;
    }

  meanz = 0.0;

  if (pt < 1.5)
    {
      widthz = widthzCoeff[0] + 
	       widthzCoeff[1] * pt +
	       widthzCoeff[2] * pt*pt +
	       widthzCoeff[3] * pt*pt*pt +
	       widthzCoeff[4] * pt*pt*pt*pt +
	       widthzCoeff[5] * pt*pt*pt*pt*pt +
	       widthzCoeff[6] * pt*pt*pt*pt*pt*pt;
    }
  else
    {
      widthz = 0.39844 -0.00660409*pt;
    }

  if ( fabs(dphi-meanphi)/widthphi < 3.0 && fabs(dz-meanz)/widthz < 3.0 ) return 1;

  return 0;
}

float MpdPairKK::Beta_sigma_El(float beta, float mom) const
{
  if (mom < 0.1) return -999;
  if (mom > 1.5) mom = 1.5;

  float mean[7] = { 3.150000e+003, -5.949618e-009, 5.973763e-002, -4.258060e-007, 7.841692e-008, -6.563439e-007, 1.891766e+002 };
  float width[7] = { -2.501047e+006, -1.253288e-011, 7.294240e-011, -1.362490e-010, 7.314968e-011, -5.672957e-010, 1.322390e-002 };

  float mean_exp, width_exp;

  mean_exp = mean[0]/mom/mom * (mean[1]*log(mom*mom) - mean[2]*mom*mom - mean[3]*mom - mean[4] - mean[5]*mom*mom*mom) + mean[6]-0.001;
  width_exp = width[0]/mom/mom * (width[1]*log(mom*mom) - width[2]*mom*mom - width[3]*mom - width[4] - width[5]*mom*mom*mom) + width[6]-0.001;

  return (beta - mean_exp)/width_exp;       
}

float MpdPairKK::Beta_sigma_Pi(float beta, float mom) const
{
  if (mom < 0.12) return -999;
  if (mom > 4.0) mom = 4.0;

  float mean[7] = { 3.150000e+003, -1.248103e-006, 1.108943e+000, -7.829288e-006, 7.832307e-006, -4.780232e-007, 3.494167e+003 };
  float width[7] = { -3.587374e+007, 1.391129e-012, 3.658606e-011, -2.025494e-011, -7.243473e-013, -3.150429e-011, 1.273604e-002 };

  float mean_exp, width_exp;

  mean_exp = mean[0]/mom/mom * (mean[1]*log(mom*mom) - mean[2]*mom*mom - mean[3]*mom - mean[4] - mean[5]*mom*mom*mom) + mean[6]-0.001;
  width_exp = width[0]/mom/mom * (width[1]*log(mom*mom) - width[2]*mom*mom - width[3]*mom - width[4] - width[5]*mom*mom*mom) + width[6]-0.001;

  return (beta - mean_exp)/width_exp;       
}

float MpdPairKK::Beta_sigma_K(float beta, float mom) const
{
  if (mom < 0.2) return -999;
  if (mom > 3.25) mom = 3.25;

  float mean[7] = { 3.150000e+003, -3.229103e-006, 1.500830e+000, 5.001743e-005, 9.494654e-006, 5.751109e-006, 4.728718e+003 };
  float width[7] = { -7.187794e+006, -9.693571e-011, 6.368089e-010, -1.548974e-009, 6.073912e-010, -3.236421e-010, 1.551237e-002 };

  float mean_exp, width_exp;

  mean_exp = mean[0]/mom/mom * (mean[1]*log(mom*mom) - mean[2]*mom*mom - mean[3]*mom - mean[4] - mean[5]*mom*mom*mom) + mean[6]-0.001;
  width_exp = width[0]/mom/mom * (width[1]*log(mom*mom) - width[2]*mom*mom - width[3]*mom - width[4] - width[5]*mom*mom*mom) + width[6]-0.001;

  return (beta - mean_exp)/width_exp;     
}

float MpdPairKK::Beta_sigma_P(float beta, float mom) const
{
  if (mom < 0.3) return -999;
  if (mom > 4.5) mom = 4.5;

  float mean[7] = { 3.150000e+003, 1.983733e-006, 1.499777e+000, 1.733988e-004, -3.201559e-005, 7.804674e-006, 4.725502e+003 };
  float width[7] = { -2.663514e+007, 1.388684e-011, 4.657636e-011, -2.603998e-011, -2.325289e-011, -1.370724e-011, 1.105921e-002 };

  float mean_exp, width_exp;

  mean_exp = mean[0]/mom/mom * (mean[1]*log(mom*mom) - mean[2]*mom*mom - mean[3]*mom - mean[4] - mean[5]*mom*mom*mom) + mean[6]-0.001;
  width_exp = width[0]/mom/mom * (width[1]*log(mom*mom) - width[2]*mom*mom - width[3]*mom - width[4] - width[5]*mom*mom*mom) + width[6]-0.001;

  return (beta - mean_exp)/width_exp;       
}

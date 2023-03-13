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
#include "TRandom.h"

#include "MpdCentralityAll.h"
#include "TFile.h"

ClassImp(MpdCentralityAll);

MpdCentralityAll::MpdCentralityAll(const char *name, const char *outputName) : MpdAnalysisTask(name, outputName)
{
   mParamConfig = outputName;
}

void MpdCentralityAll::UserInit()
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
   mhVertexAcc = new TH1F("hVertexAcc", "Accepted event vertex distribution", 100, -200., 200.);
   fOutputList->Add(mhVertexAcc);
   mhHits = new TH1F("hHits", "Number of TPC hits", 100, -0.5, 99.5);
   fOutputList->Add(mhHits);
   mhEta = new TH1F("hEta", "Eta", 100, -2., 2.);
   fOutputList->Add(mhEta);
   mhPt = new TH1F("hPt", "Pt", 300, 0., 3.);
   fOutputList->Add(mhPt);
   mhDca = new TH1F("hDca", "DCA", 100, -5., 5.);
   fOutputList->Add(mhDca);
   mhMultiplicity = new TH1F("hMultiplicity", "Multiplicity distribution", 2000, -0.5, 1999.5);
   fOutputList->Add(mhMultiplicity);
   mhMultiplicityEff = new TH1F("hMultiplicityEff", "Weighted multiplicity distribution", 2000, -0.5, 1999.5);
   fOutputList->Add(mhMultiplicityEff);
   mhCentrality = new TH1F("hCentrality", "Centrality distribution", 103, -2., 101.);
   fOutputList->Add(mhCentrality);
   mhCentConvert = new TH1F("hCentConvert", "nTr-Centrality converter", 2000, 0.5, 2000.5);
   fOutputList->Add(mhCentConvert);
   mhTrEff = new TH2F("hTrEff", "Track Efficiency", 100, -200, 200, 100, -4, 4);
   fOutputList->Add(mhTrEff);

   TDatime tim;
   RND.SetSeed(tim.GetTime());

   mhCentConvert->Reset();

   if (mParams.mProdGenerator != "ANY" && mParams.mInFileConvert != "ANY") {
      cout << "evCentrality: Reading out nTrack-to-Centrality conversion table from file " << mParams.mInFileConvert
           << endl;
      TFile *inr = new TFile((mParams.mInFileConvert).c_str());
      mhCentConvert->Add(((TH1F *)inr->Get("hCentr")));
      inr->Close();
   } else {
      for (int i = 0; i < 2000; i++) mhCentConvert->SetBinContent(i + 1, -1.);
   }

   mhTrEff->Reset();

   if (mParams.mInFileTrEff != "ANY") {
      cout << "evCentrality: Reading out track reconstruction values from file " << mParams.mInFileTrEff << endl
           << endl;
      TFile *inr = new TFile((mParams.mInFileTrEff).c_str());
      mhTrEff->Add(((TH2F *)inr->Get("RecEffHad_DCMSMM92")));
      inr->Close();
   } else {
      for (int i = 0; i < 100; i++) {
         for (int j = 0; j < 100; j++) mhTrEff->SetBinContent(i + 1, j + 1, 1.);
      }
   }

   // MC
   if (isMC) {
   }
}
//--------------------------------------
void MpdCentralityAll::ProcessEvent(MpdAnalysisEvent &event)
{
   event.setCentrTPC(-1);

   if (!isInitialized) {
      mKF = MpdKalmanFilter::Instance();
      mKHit.SetType(MpdKalmanHit::kFixedR);
      isInitialized = true;
   }

   if (!selectEvent(event)) {
      return;
   }

   int   nTrTPC    = 0;
   float nTrTPCEff = 0;

   mMpdGlobalTracks = event.fMPDEvent->GetGlobalTracks();

   for (long int i = 0; i < mMpdGlobalTracks->GetEntriesFast(); i++) {
      MpdTrack          *mpdtrack = (MpdTrack *)mMpdGlobalTracks->UncheckedAt(i);
      MpdTpcKalmanTrack *tr       = (MpdTpcKalmanTrack *)mKalmanTracks->UncheckedAt(i);
      if (!selectTrack(mpdtrack)) {
         continue;
      }

      nTrTPC++;

      int   binn   = mhTrEff->FindBin(mPrimaryVertex.Z(), mpdtrack->GetEta());
      float weight = mhTrEff->GetBinContent(binn);
      if (weight != 0) weight = 1. / weight;

      nTrTPCEff += weight;
   }

   // Emulation of trigger efficiecny by TrigEffMult(N_TPC)
   if (nTrTPC >= 15 || (nTrTPC < 15 && RND.Rndm() < TrigEffMult[nTrTPC])) {
      mhEvents->Fill(3.5);
      mhVertexAcc->Fill(mPrimaryVertex.Z());

      mhMultiplicity->Fill(nTrTPC);
      mhMultiplicityEff->Fill(nTrTPCEff);

      float cen = mhCentConvert->GetBinContent(int(nTrTPCEff));

      event.setCentrTPC(cen);

      // Centrality to be passed in the remaining wagons
      mhCentrality->Fill(event.getCentrTPC());
   }
}

void MpdCentralityAll::Finish()
{
   // Post-scan processing not needed
}

//--------------------------------------
bool MpdCentralityAll::selectEvent(MpdAnalysisEvent &event)
{
   mhEvents->Fill(0.5);

   // Reject empty events (UrQMD, PHSD)
   mMCTracks = event.fMCTrack;

   int nTrMc = 0;
   for (int i = 0; i < mMCTracks->GetEntriesFast(); i++) {
      MpdMCTrack *pr = (static_cast<MpdMCTrack *>(mMCTracks->At(i)));
      if (pr->GetMotherId() == -1) {
         nTrMc++;
      }
   } // i

   if (mParams.mProdGenerator == "Req25-UrQMD" && nTrMc <= 2 * 209) { // Just nucleons of Bi+Bi --> Request25-UrQMD
      return false;
   }

   if (mParams.mProdGenerator == "Req30-PHSD" && nTrMc <= 2 * 209) { // Just nucleons of Bi+Bi --> Request30-PHSD
      return false;
   }

   mhEvents->Fill(1.5);

   // Reject bad vertex
   if (!event.fVertex) {
      return false;
   }

   MpdVertex *vertex = (MpdVertex *)event.fVertex->First();
   vertex->Position(mPrimaryVertex);
   mhVertex->Fill(mPrimaryVertex.Z());

   if (mPrimaryVertex.Z() == 0) { // not reconstructed (==0)
      return false;
   }

   if (fabs(mPrimaryVertex.Z()) > mParams.mZvtxCut) { // beyond the limits
      return false;
   }

   mhEvents->Fill(2.5);

   return true;
}

bool MpdCentralityAll::selectTrack(MpdTrack *mpdtrack)
{
   if (mpdtrack->GetNofHits() < mParams.mNofHitsCut) return false; // nhits > 10

   if (fabs(mpdtrack->GetEta()) > mParams.mEtaCut) return false; // |eta| < 0.5

   float pt = TMath::Abs(mpdtrack->GetPt());
   if (pt < mParams.mPtminCut) return false; // pT > 1000 MeV/c

   if (fabs(mpdtrack->GetDCAX()) > mParams.mDcaCut) return false; // |DCAx| < 2.
   if (fabs(mpdtrack->GetDCAY()) > mParams.mDcaCut) return false; // |DCAy| < 2.
   if (fabs(mpdtrack->GetDCAZ()) > mParams.mDcaCut) return false; // |DCAz| < 2.

   mhHits->Fill(mpdtrack->GetNofHits());
   mhPt->Fill(pt);
   mhEta->Fill(mpdtrack->GetEta());
   mhDca->Fill(mpdtrack->GetDCAX());
   mhDca->Fill(mpdtrack->GetDCAY());
   mhDca->Fill(mpdtrack->GetDCAZ());

   return true;
}

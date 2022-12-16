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
   mhHits = new TH1F("hHits", "Number of TPC hits", 100, -0.5, 99.5);
   fOutputList->Add(mhHits);
   mhEta = new TH1F("hEta", "Eta", 100, -2., 2.);
   fOutputList->Add(mhEta);
   mhPt = new TH1F("hPt", "Pt", 300, 0., 3.);
   fOutputList->Add(mhPt);
   mhMultiplicity = new TH1F("hMultiplicity", "Multiplicity distribution", 2000, -0.5, 1999.5);
   fOutputList->Add(mhMultiplicity);
   mhCentrality = new TH1F("hCentrality", "Centrality distribution", 102, -2., 100.);
   fOutputList->Add(mhCentrality);

   // MC
   if (isMC) {
   }
}
//--------------------------------------
void MpdCentralityAll::ProcessEvent(MpdAnalysisEvent &event)
{

   if (!isInitialized) {
      mKF = MpdKalmanFilter::Instance();
      mKHit.SetType(MpdKalmanHit::kFixedR);
      isInitialized = true;
   }

   if (!selectEvent(event)) { //(V)
      return;
   }

   // TPC centrality
   int ntr          = 0;
   mMpdGlobalTracks = event.fMPDEvent->GetGlobalTracks();

   for (long int i = 0; i < mMpdGlobalTracks->GetEntriesFast(); i++) {
      MpdTrack          *mpdtrack = (MpdTrack *)mMpdGlobalTracks->UncheckedAt(i);
      MpdTpcKalmanTrack *tr       = (MpdTpcKalmanTrack *)mKalmanTracks->UncheckedAt(i);
      if (!selectTrack(mpdtrack)) {
         continue;
      }
      ntr++;
   }

   // Multiplicity
   mhMultiplicity->Fill(ntr);

   float cen = -1;
   if (ntr > 0) cen = 100. - 100. * float(ntr) / 500.; //(%) very rough, just a placeholder
   if (cen > 100.) cen = 100.;

   // Centrality
   mhCentrality->Fill(cen);

   event.setCentrTPC(cen);
}

void MpdCentralityAll::Finish()
{
   // Post-scan processing not needed
}

//--------------------------------------
bool MpdCentralityAll::selectEvent(MpdAnalysisEvent &event)
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

   return true;
}

bool MpdCentralityAll::selectTrack(MpdTrack *mpdtrack)
{
   if (mpdtrack->GetNofHits() < mParams.mNofHitsCut) return false; // nhits > 10

   if (fabs(mpdtrack->GetEta()) > mParams.mEtaCut) return false; //|eta| < 1.0

   float pt = TMath::Abs(mpdtrack->GetPt());
   if (pt < mParams.mPtminCut) return false; // pT > 50 MeV/c

   mhHits->Fill(mpdtrack->GetNofHits());
   mhPt->Fill(pt);
   mhEta->Fill(mpdtrack->GetEta());

   return true;
}

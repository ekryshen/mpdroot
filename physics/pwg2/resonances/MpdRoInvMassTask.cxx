// skelet

#ifndef ROOT_MpdRoInvMassTask
#include "MpdRoInvMassTask.h"
#endif

#include <TLorentzVector.h>
#include "FairRootManager.h"
#include "MpdMCTrack.h"
#include "MpdTrack.h"

#include <iostream>
using std::cout;
using std::endl;

// -----   Default constructor ---------------------------------------
MpdRoInvMassTask::MpdRoInvMassTask():
  FairTask(),
  fEventCounter(0)
{}

// -----   constructor with names ------------------------------------
MpdRoInvMassTask::MpdRoInvMassTask(const char *name, const char *title):
  FairTask(name),
  fEventCounter(0)
{}

// -----   Destructor -----------------------------------------------
MpdRoInvMassTask::~MpdRoInvMassTask()
{}

// -------------------------------------------------------------------
InitStatus MpdRoInvMassTask::Init()
{
  cout<<"\n-I- [MpdRoInvMassTask::Init] " <<endl;  
  
  FairRootManager *manager = FairRootManager::Instance();

  fMCTracks = (TClonesArray *) manager->GetObject("MCTrack");
  fDstEvent = (MpdEvent *) manager->GetObject("MPDEvent.");
  Register();
  
  fRoInvMass = new TH1F("RoInvMass","Invariant Mass: (#pi_{+}#pi_{-}) 90%CL",200,0.,2.);
  fRoInvMass->GetXaxis()->SetTitle("M_{inv}");
  fRoInvMassMC = new TH1F("RoInvMassMC","Invariant Mass (MC): (#pi_{+}#pi_{-}) 90%CL",200,0.,2.);
  fRoInvMassMC->GetXaxis()->SetTitle("M_{inv}");
  fRoInvMassMC->SetLineColor(kRed);


  fPDG = TDatabasePDG::Instance();

  return kSUCCESS;
}

// -------------------------------------------------------------------
void MpdRoInvMassTask::Exec(Option_t * option)
{
  fEventCounter++;
  cout<<"-I- [MpdRoInvMassTask::Exec] " << "{" << fEventCounter << "}" <<endl;
    
  //MpdEvent *mpdEvent = fDstEvent; // one event by one event
  //event->Dump(); 
  TClonesArray *mpdTracks = fDstEvent->GetGlobalTracks();
  Int_t nTracks = mpdTracks->GetEntriesFast();
  
  cout << "N of Reconstructed tracks = " << nTracks << endl;
    
  //mpdTracks->Dump();
  
  /* 
   * Let us calculate invariant mass of rho0 meson (rho0-> pi+ + pi-)
   */
  
  TLorentzVector lPos, lNeg;
  
  /* Events loop */
  for (Int_t i = 0; i < nTracks; i++)
  {
	MpdTrack *track1 = (MpdTrack *) mpdTracks->UncheckedAt(i);
	if ( !track1->GetTofFlag() ) continue;  // no tof identification
	
	Float_t pion1 = track1->GetPidProbPion();
    if ( !(pion1 > 0.9) ) continue; // accept only pi+ 90% CL, reject other
    //track1->Dump();
    lPos.SetXYZM( track1->GetPx(), track1->GetPy(), track1->GetPz(), fPDG->GetParticle(211)->Mass() );
    
    for (Int_t j = 0; j<i; j++)
    {
	  MpdTrack *track2 = (MpdTrack *) mpdTracks->UncheckedAt(j);
	  if ( !track2->GetTofFlag() ) continue;  // no tof identification
	  
	  Float_t pion2 = track2->GetPidProbPion();
      if( !(pion2 < 0.9) ) continue; // accept only pi- 90% CL, reject other
      //track2->Dump();
      lNeg.SetXYZM( track2->GetPx(), track2->GetPy(), track2->GetPz(), fPDG->GetParticle(-211)->Mass() );
      
      Float_t minv = (lPos + lNeg).Mag();
      
      fRoInvMass->Fill(minv);
    }
  }
  /* End of events loop */
  
  /* Events loop. The same but for MC */
  Int_t nMCTracks = fMCTracks->GetEntriesFast();  
  
  cout << "N of MC tracks = " << nMCTracks << endl;
  
  for (Int_t i = 0; i < nMCTracks; i++)
  {
    MpdMCTrack *MCtrack1 = (MpdMCTrack*) fMCTracks->UncheckedAt(i);
    //if( MCtrack1->GetMotherId() >0 ) continue; // start from initial vertex
    Int_t pdgCode1 = MCtrack1->GetPdgCode();
  
    if ( pdgCode1!=211 ) continue;  // pi+
  
    MCtrack1->Get4Momentum(lPos);
    
    for (Int_t j = 0; j<i; j++)
	{
	  MpdMCTrack *MCtrack2 = (MpdMCTrack*) fMCTracks->UncheckedAt(j);
	  //if( MCtrack2->GetMotherId() >0 ) continue; // start from initial vertex
	  Int_t pdgCode2 = MCtrack2->GetPdgCode();
      
   	  if ( pdgCode2!=-211 ) continue;  // pi-
      
	  MCtrack2->Get4Momentum(lNeg);
      
	  Float_t minv = (lPos + lNeg).Mag();
      
	  fRoInvMassMC->Fill(minv); // with probability weight
	}
    }
  /* End of MC events loop */

}

// -------------------------------------------------------------------
void MpdRoInvMassTask::Reset()
{}

// -------------------------------------------------------------------
void MpdRoInvMassTask::Finish()
{
  //cout<<"\n-I- [MpdRoInvMassTask::Finish] "<< endl;
  
  //TFile fileOut("MinvRho0.root","recreate");
  Double_t scale = 1./ (Double_t) fEventCounter;
  
  fRoInvMass->Scale(scale); // scale to total number of events
  fRoInvMassMC->Scale(scale); // scale to total number of events

  fRoInvMass->Write("");
  fRoInvMassMC->Write("");
}

// -------------------------------------------------------------------
void MpdRoInvMassTask::Register()
{
  //FairRootManager::Instance()->Register("MCTrack", "MC", fMCTracks, kTRUE);
}


// -------------------------------------------------------------------
ClassImp(MpdRoInvMassTask);

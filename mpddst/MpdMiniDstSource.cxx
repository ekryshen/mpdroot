/*
 * MpdMiniDstSource.cxx
 *
 *  Created on: 17 kwi 2020
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */

#include "MpdMiniDstSource.h"

#include <FairLogger.h>
#include <FairRootManager.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TTree.h>
#include <fstream>

MpdMiniDstSource::MpdMiniDstSource() : MpdMiniDstSource("data.root") {}

MpdMiniDstSource::MpdMiniDstSource(TString inFile)
   : fChain(nullptr), fEvent(nullptr), fTracks(nullptr), fTofInfo(nullptr), fEmcInfo(nullptr), fMcEvent(nullptr),
     fMcTracks(nullptr), fCovMatrix(nullptr), fMaxEventsNo(0)
{
   if (inFile.Length() == 0) return;
   if (inFile.EndsWith("root")) {
      fFileName.push_back(inFile);
   } else if (inFile.EndsWith("list")) {
      std::ifstream list;
      list.open(inFile);
      while (!list.eof()) {
         TString temp;
         list >> temp;
         fFileName.push_back(temp);
      }
      list.close();
   }
}

Bool_t MpdMiniDstSource::Init()
{
   FairRootManager *mngr = FairRootManager::Instance();
   fChain                = new TChain("MiniDst");
   for (unsigned int j = 0; j < fFileName.size(); j++) {
      LOG(DEBUG3) << "MpdMiniDstSource: opening root file: " << fFileName[j];
      fChain->Add(fFileName[j], 0);
   }

   fChain->SetBranchStatus("Event", 1);
   fChain->SetBranchStatus("Track", 1);
   fChain->SetBranchStatus("BTofPidTraits", 1);
   fChain->SetBranchStatus("TrackCovMatrix", 1);
   fChain->SetBranchStatus("BECalCluster", 1);
   fChain->SetBranchStatus("McEvent", 1);
   fChain->SetBranchStatus("McTrack", 1);
   fEvent     = new TClonesArray("MpdMiniEvent");
   fMcEvent   = new TClonesArray("MpdMiniMcEvent");
   fTracks    = new TClonesArray("MpdMiniTrack");
   fTofInfo   = new TClonesArray("MpdMiniBTofPidTraits");
   fEmcInfo   = new TClonesArray("MpdMiniBECalCluster");
   fMcTracks  = new TClonesArray("MpdMiniMcTrack");
   fCovMatrix = new TClonesArray("MpdMiniTrackCovMatrix");

   fChain->SetBranchAddress("Event", &fEvent);
   fChain->SetBranchAddress("Track", &fTracks);
   fChain->SetBranchAddress("BTofPidTraits", &fTofInfo);
   fChain->SetBranchAddress("McEvent", &fMcEvent);
   fChain->SetBranchAddress("McTrack", &fMcTracks);
   fChain->SetBranchAddress("TrackCovMatrix", &fCovMatrix);
   fChain->SetBranchAddress("BECalCluster", &fEmcInfo);
   fMaxEventsNo = fChain->GetEntries();

   mngr->SetInChain(fChain, -1);
   mngr->RegisterInputObject("Event", fEvent);
   mngr->RegisterInputObject("Track", fTracks);
   mngr->RegisterInputObject("BTofPidTraits", fTofInfo);
   mngr->RegisterInputObject("McEvent", fMcEvent);
   mngr->RegisterInputObject("McTrack", fMcTracks);
   mngr->RegisterInputObject("TrackCovMatrix", fCovMatrix);
   mngr->RegisterInputObject("BECalCluster", fEmcInfo);
   return kTRUE;
}

Int_t MpdMiniDstSource::ReadEvent(UInt_t unsignedInt)
{
   fChain->GetEntry(unsignedInt);
   return 0;
}

void MpdMiniDstSource::Close() {}

Int_t MpdMiniDstSource::CheckMaxEventNo(Int_t int1)
{
   return fMaxEventsNo;
}

MpdMiniDstSource::MpdMiniDstSource(const MpdMiniDstSource &other)
{
   fFileName    = other.fFileName;
   fMaxEventsNo = other.fMaxEventsNo;
   if (!other.fChain) {
      Init();
   }
}

MpdMiniDstSource::~MpdMiniDstSource() {}

void MpdMiniDstSource::AddFile(TString file)
{
   fFileName.push_back(file);
}

/*
 * NicaMpdV0Matcher.cxx
 *
 *  Created on: 27 lut 2021
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "NicaMpdV0Matcher.h"
#include "MpdMiniMcTrack.h"
#include "MpdMiniTrack.h"
#include "NicaLink.h"
#include "NicaTrackClones.h"
#include "NicaV0Track.h"
#include <FairLogger.h>

#include <TList.h>
#include <iostream>
#include <vector>

NicaMpdV0Matcher::NicaMpdV0Matcher() :
  fV0Tracks(nullptr), fV0Links(nullptr), fMcTracks(nullptr), fReTracks(nullptr), fWrite(kFALSE), fLastPrimary(-1) {}

InitStatus NicaMpdV0Matcher::Init() {
  FairRootManager* mngr = FairRootManager::Instance();
  mngr->UpdateBranches();
  if (mngr->GetObject("NicaV0Matches") != nullptr) {
    LOG(warning) << "NicaV0Matches present in tree";
    return kERROR;
  }
  TClonesArray* a = new TClonesArray("TLorentzVector");
  mngr->Register("test", "T", a, kFALSE);
  if (mngr->GetObject("NicaV0Tracks") == nullptr) {
    TList* l = mngr->GetBranchNameList();
    for (int i = 0; i < l->GetEntries(); i++) {
      TObjString* ob = (TObjString*) l->At(i);
      std::cout << "BR " << ob->GetString() << mngr->CheckBranch(ob->GetString()) << " " << mngr->GetObject(ob->GetString())
                << std::endl;
    }
    LOG(warning) << "Lack of V0Tracks needed for NicaMpdV0Matcher" << mngr->CheckBranch("NicaV0Tracks");
    return kERROR;
  }
  fV0Tracks = new NicaTrackClones("NicaV0Track", "NicaV0Tracks", "V0");
  fV0Links  = new NicaTrackClones("NicaLink", "NicaV0Matches", "V0");
  fMcTracks = new NicaTrackClones("MpdMiniMcTrack", "McTrack", "mc");
  fReTracks = new NicaTrackClones("MpdMiniTrack", "Track", "track");
  fV0Links->Register(fWrite);
  fReTracks->GetFromTree();
  fMcTracks->GetFromTree();
  fV0Tracks->GetFromTree();
  return kSUCCESS;
}

void NicaMpdV0Matcher::Exec(Option_t* opt) {
  fV0Links->GetArray()->Clear();
  fLastPrimary = -1;
  std::cout << "V0s " << fV0Tracks->GetEntriesFast() << std::endl;
  TClonesArray* matches = fV0Links->GetArray();
  matches->ExpandCreateFast(fV0Tracks->GetEntriesFast());
  for (int i = 0; i < fV0Tracks->GetEntriesFast(); i++) {
    NicaV0Track* v0       = (NicaV0Track*) fV0Tracks->UncheckedAt(i);
    Int_t posId           = v0->GetPosId();
    Int_t negId           = v0->GetNegId();
    Double_t mPos         = v0->GetMomPos().M();
    Double_t mNeg         = v0->GetMomNeg().M();
    MpdMiniTrack* recoPos = (MpdMiniTrack*) fReTracks->UncheckedAt(posId);
    MpdMiniTrack* recoNeg = (MpdMiniTrack*) fReTracks->UncheckedAt(negId);
    std::cout << posId << " " << negId << ">" << std::endl;
    posId = recoPos->mcTrackIndex();
    negId = recoNeg->mcTrackIndex();
    std::cout << posId << " " << negId << ">" << std::endl;
    std::cout << fMcTracks->GetEntriesFast() << std::endl;

    MpdMiniMcTrack* daughterPos = nullptr;
    MpdMiniMcTrack* daughterNeg = nullptr;


    if (posId >= 0 && posId < fMcTracks->GetEntriesFast()) daughterPos = (MpdMiniMcTrack*) fMcTracks->UncheckedAt(posId);
    if (negId >= 0 && negId < fMcTracks->GetEntriesFast()) daughterNeg = (MpdMiniMcTrack*) fMcTracks->UncheckedAt(negId);
    if (daughterPos) {
      if (daughterPos->isFromGenerator()) daughterPos = nullptr;
    }
    if (daughterNeg) {
      if (daughterNeg->isFromGenerator()) daughterNeg = nullptr;
    }
    MpdMiniMcTrack* daughter = nullptr;
    if (mPos > mNeg && daughterPos != nullptr)
      daughter = daughterPos;
    else
      daughter = daughterNeg;

    NicaLink* link = (NicaLink*) matches->UncheckedAt(i);
    if (daughter == nullptr) {
      link->SetLink(0, -1);
      continue;
    }
    if (daughter->isFromGenerator()) {
      link->SetLink(0, -1);
      continue;
    }
    Int_t motherID     = FindMotherIndex(daughter, v0);
    MpdMiniMcTrack* mc = (MpdMiniMcTrack*) fMcTracks->UncheckedAt(motherID);
    if (mc == nullptr) {
      link->SetLink(0, -1);
      continue;
    }
    link->SetLink(0, motherID);
  }
}

Int_t NicaMpdV0Matcher::FindMotherIndex(MpdMiniMcTrack* mcTrack, NicaV0Track* v0) {
  if (fLastPrimary == -1) {  // find last primary particle
    for (int i = 0; i < fMcTracks->GetEntriesFast(); i++) {
      MpdMiniMcTrack* mc = (MpdMiniMcTrack*) fMcTracks->UncheckedAt(i);
      if (!mc->isFromGenerator()) {
        fLastPrimary = i - 1;
        break;
      }
    }
  }
  Int_t index     = -1;
  Double_t minMag = 1E+9;
  for (int i = fLastPrimary; i < fMcTracks->GetEntriesFast(); i++) {
    MpdMiniMcTrack* mc = (MpdMiniMcTrack*) fMcTracks->UncheckedAt(i);
    if (mc->pdgId() == v0->GetPdg()) {  // potential parent
      Double_t dpx = mc->px() - v0->GetMom().Px();
      Double_t dpy = mc->py() - v0->GetMom().Py();
      Double_t dpz = mc->pz() - v0->GetMom().Pz();
      Double_t dp  = dpx * dpx + dpy * dpy + dpz * dpz;
      if (dp < minMag) {
        minMag = dp;
        index  = i;
      }
    }
  }
  return index;
}

NicaMpdV0Matcher::~NicaMpdV0Matcher() {
  // TODO
}

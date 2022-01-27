/*
 * MpdV0FinderHelix.cxx
 *
 *  Created on: 16 gru 2021
 *      Author: Daniel Wielanek
 *      E-mail: daniel.wielanek at pw.edu.pl
 *      Warsaw University of Technology, Faculty of Physics
 *      JINR,  Laboratory of High Energy Physics
 */

#include "MpdV0FinderHelix.h"

#include <FairTask.h>
#include <RtypesCore.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TObjArray.h>
#include <TVector3.h>
#include <iostream>
#include <utility>

#include "fairlogger/Logger.h"
#include "MpdEvent.h"
#include "MpdTrack.h"
#include "MpdMiniEvent.h"
#include "MpdMiniTrack.h"
#include "MpdV0CandidateCut.h"
#include "MpdV0DaughterCut.h"
#include "MpdV0Particle.h"
#include "MpdV0Track.h"

InitStatus MpdV0FinderHelix::Init()
{
   InitStatus stat = MpdV0FinderBasic::Init();
   if (stat != kSUCCESS) return stat;
   return kSUCCESS;
}

void MpdV0FinderHelix::ExecMiniDst(Option_t *option)
{
   MpdMiniEvent *event  = static_cast<MpdMiniEvent *>(fMiniEvents->UncheckedAt(0));
   TVector3      vertex = event->primaryVertex();
   Double_t      field  = event->bField();
   fFirstHelix.clear();
   fSecondHelix.clear();
   fFirstDaughterCut->SetMiniDstEventInfo(*event);
   fSecondDaughterCut->SetMiniDstEventInfo(*event);
   std::pair<TObject *, Int_t> data;

   for (int i = 0; i < fMiniTracks->GetEntriesFast(); i++) {
      MpdMiniTrack *track = (MpdMiniTrack *)fMiniTracks->UncheckedAt(i);
      data.first          = track;
      data.second         = i;
      if (fFirstDaughterCut->PassMiniDstTrack(*track)) {
         fFirstDaughters.push_back(data);
         fFirstHelix.push_back(NestedHelix(track->gMom(), track->gDCA(vertex), track->charge(), field * 0.1));
      } else if (fSecondDaughterCut->PassMiniDstTrack(*track)) {
         fSecondDaughters.push_back(data);
         fSecondHelix.push_back(NestedHelix(track->gMom(), track->gDCA(vertex), track->charge(), field * 0.1));
      }
   }

   MpdV0Track candidate;
   for (int i = 0; i < fFirstDaughters.size(); i++) {
      std::pair<TObject *, Int_t> data1;
      data1                = fFirstDaughters[i];
      MpdMiniTrack *track1 = (MpdMiniTrack *)data1.first;
      Double_t      p1Tot  = track1->gMom().Mag();
      candidate.SetFirstDaughterIndex(data1.second);
      NestedHelix &h1 = fFirstHelix[i];

      for (int j = 0; j < fSecondDaughters.size(); j++) {
         std::pair<TObject *, Int_t> data2;
         if (data1.second == data2.second) continue;
         data2                = fSecondDaughters[j];
         MpdMiniTrack *track2 = (MpdMiniTrack *)data2.first;
         Double_t      p2Tot  = track2->gMom().Mag();
         candidate.SetSecondDaughterIndex(data2.second);

         NestedHelix &             h2 = fSecondHelix[j];
         std::pair<double, double> s  = h1.pathLengths(h2);
         if (s.first < 0) s.first += h1.period();
         if (s.second < 0) s.second += h2.period();
         TVector3 pos1 = h1.at(s.first);
         TVector3 mom1 = h1.cat(s.first) * p1Tot;
         TVector3 pos2 = h2.at(s.second);
         TVector3 mom2 = h2.cat(s.second) * p2Tot;
         candidate.SetMomFirstDaughter(mom1);
         candidate.SetMomSecondDaughter(mom2);
         candidate.SetDecayPoint((pos1 + pos2) * 0.5);
         candidate.SetDau1to2((pos1 - pos2).Mag());
         candidate.SetFirstDaughterS(s.first);
         candidate.SetSecondDaughterS(s.second);
         candidate.Recalculate(vertex);

         if (fCandicateCut->Pass(candidate)) {
            MpdV0Track *v0 = (MpdV0Track *)fV0s->ConstructedAt(fV0s->GetEntriesFast());
            *v0            = candidate;
         }
      }
   }
   LOG(info) << "Found daughters: " << fFirstDaughters.size() << " " << fSecondDaughters.size();
}

void MpdV0FinderHelix::ExecDst(Option_t *option)
{
   TClonesArray *tracks = fMpdEvent->GetGlobalTracks();
   fFirstHelix.clear();
   fSecondHelix.clear();
   std::pair<TObject *, Int_t> data;
   for (int i = 0; i < fMiniTracks->GetEntriesFast(); i++) {
      MpdTrack *track = (MpdTrack *)fMiniTracks->UncheckedAt(i);
      data.first      = track;
      data.second     = i;
      if (fFirstDaughterCut->PassDstTrack(*track)) {
         fFirstDaughters.push_back(data);
         fFirstHelix.push_back(NestedHelix(track->GetHelix()));
      } else if (fSecondDaughterCut->PassDstTrack(*track)) {
         fSecondDaughters.push_back(data);
         fSecondHelix.push_back(NestedHelix(track->GetHelix()));
      }
   }

   MpdV0Track candidate;
   TVector3   vertex(fMpdEvent->GetPrimaryVerticesX(), fMpdEvent->GetPrimaryVerticesY(),
                   fMpdEvent->GetPrimaryVerticesZ());
   for (int i = 0; i < fFirstDaughters.size(); i++) {
      std::pair<TObject *, Int_t> data1;
      data1            = fFirstDaughters[i];
      MpdTrack *track1 = (MpdTrack *)data1.first;
      Double_t  pt1    = track1->GetPt();
      Double_t  pz1    = track1->GetPz();
      Double_t  p1Tot  = TMath::Sqrt(pt1 * pt1 + pz1 * pz1);
      candidate.SetFirstDaughterIndex(data1.second);
      NestedHelix &h1 = fFirstHelix[i];

      for (int j = 0; j < fSecondDaughters.size(); j++) {
         std::pair<TObject *, Int_t> data2;
         if (data1.second == data2.second) continue;
         data2            = fSecondDaughters[j];
         MpdTrack *track2 = (MpdTrack *)data2.first;
         Double_t  pt2    = track2->GetPt();
         Double_t  pz2    = track2->GetPz();
         Double_t  p2Tot  = TMath::Sqrt(pt2 * pt2 + pz2 * pz2);
         candidate.SetSecondDaughterIndex(data2.second);

         NestedHelix &             h2 = fSecondHelix[j];
         std::pair<double, double> s  = h1.pathLengths(h2);
         if (s.first < 0) s.first += h1.period();
         if (s.second < 0) s.second += h2.period();
         TVector3 pos1 = h1.at(s.first);
         TVector3 mom1 = h1.cat(s.first) * p1Tot;
         TVector3 pos2 = h2.at(s.second);
         TVector3 mom2 = h2.cat(s.second) * p2Tot;
         candidate.SetMomFirstDaughter(mom1);
         candidate.SetMomSecondDaughter(mom2);
         candidate.SetDecayPoint((pos1 + pos2) * 0.5);
         candidate.SetDau1to2((pos1 - pos2).Mag());
         candidate.SetFirstDaughterS(s.first);
         candidate.SetSecondDaughterS(s.second);
         candidate.Recalculate(vertex);

         if (fCandicateCut->Pass(candidate)) {
            MpdV0Track *v0 = (MpdV0Track *)fV0s->ConstructedAt(fV0s->GetEntriesFast());
            *v0            = candidate;
         }
      }
   }
}

MpdV0FinderHelix::MpdV0FinderHelix(Int_t pidMom, Int_t pidFirstDau, Int_t pidSecDau)
   : MpdV0FinderBasic(pidMom, pidFirstDau, pidSecDau), fFirstHelix(), fSecondHelix()
{
}

/*
 * MpdV0FinderBasic.cxx
 *
 *  Created on: 16 gru 2021
 *      Author: Daniel Wielanek
 *      E-mail: daniel.wielanek at pw.edu.pl
 *      Warsaw University of Technology, Faculty of Physics
 *      JINR,  Laboratory of High Energy Physics
 */

#include "MpdV0Finder.h"

#include <FairRootManager.h>
#include <TClonesArray.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TNamed.h>

#include "FairLogger.h"
#include "MpdEvent.h"
#include "MpdMiniEvent.h"
#include "MpdV0CandidateCut.h"
#include "MpdV0DaughterCut.h"
#include "MpdV0CandidateMonitor.h"
#include "MpdV0DaughterMonitor.h"

MpdV0Finder::MpdV0Finder(TString name, Int_t pidMom, Int_t pidFirstDau, Int_t pidSecDau)
   : FairTask(name), fInit(kFALSE), fWrite(kFALSE), fFirstV0(kFALSE), fPidDauPos(pidSecDau), fPidDauNeg(pidFirstDau),
     fPidV0(pidMom), fFormat(EFormat::kDst), fMpdEvent(nullptr), fMiniEvents(nullptr), fMiniTracks(nullptr),
     fMiniTofData(nullptr), fV0s(nullptr), fPositiveDaughterCut(nullptr), fNegativeDaughterCut(nullptr),
     fCandicateCut(nullptr), fDauMon1(nullptr), fDauMon2(nullptr), fV0Mon(nullptr)
{
}

MpdV0Finder::~MpdV0Finder()
{
   if (fPositiveDaughterCut) delete fPositiveDaughterCut;
   if (fNegativeDaughterCut) delete fNegativeDaughterCut;
   if (fCandicateCut) delete fCandicateCut;
   if (fDauMon1) delete fDauMon1;
   if (fDauMon2) delete fDauMon2;
   if (fV0Mon) delete fV0Mon;
}

MpdV0Finder::MpdV0Finder(const MpdV0Finder &other)
   : FairTask(), fInit(kFALSE), fWrite(other.fWrite), fFirstV0(other.fFirstV0), fPidDauPos(other.fPidDauPos),
     fPidDauNeg(other.fPidDauNeg), fPidV0(other.fPidV0), fFormat(other.fFormat), fMpdEvent(nullptr),
     fMiniEvents(nullptr), fMiniTracks(nullptr), fMiniTofData(nullptr), fV0s(nullptr), fPositiveDaughterCut(nullptr),
     fNegativeDaughterCut(nullptr), fCandicateCut(nullptr), fDauMon1(nullptr), fDauMon2(nullptr), fV0Mon(nullptr)
{
   if (other.fPositiveDaughterCut) {
      fPositiveDaughterCut = (MpdV0DaughterCut *)other.fPositiveDaughterCut->Clone();
   }
   if (other.fNegativeDaughterCut) {
      fNegativeDaughterCut = (MpdV0DaughterCut *)other.fNegativeDaughterCut->Clone();
   }
   if (other.fCandicateCut) {
      fCandicateCut = (MpdV0CandidateCut *)other.fCandicateCut->Clone();
   }
   if (other.fDauMon1) {
      fDauMon1 = other.fDauMon1->MakeCopy();
   }
   if (other.fDauMon2) {
      fDauMon2 = other.fDauMon1->MakeCopy();
   }
   if (other.fV0Mon) {
      fV0Mon = other.fV0Mon->MakeCopy();
   }
}

void MpdV0Finder::Exec(Option_t *option)
{
   if (fFirstV0) fV0s->Clear();
   fPositiveDaughters.clear();
   fNegativeDaughters.clear();
   switch (fFormat) {
   case EFormat::kDst: {

      SetDstData();
      ExecDst(option);
      break;
   }
   case EFormat::kMiniDst: {
      SetMiniDstData();
      ExecMiniDst(option);
   }
   }
}

InitStatus MpdV0Finder::Init()
{
   if (fInit) {
      LOG(ERROR) << "Task already initalized";
      return kSUCCESS;
   }
   FairRootManager *mngr = FairRootManager::Instance();
   fMpdEvent             = (MpdEvent *)mngr->GetObject("MPDEvent.");
   fMiniEvents           = (TClonesArray *)mngr->GetObject("Event");
   fMiniTracks           = (TClonesArray *)mngr->GetObject("Track");
   fMiniTofData          = (TClonesArray *)mngr->GetObject("BTofPidTraits");
   if (mngr->GetObject("MpdV0")) {
      fFirstV0 = kFALSE;
      fV0s     = (TClonesArray *)mngr->GetObject("MpdV0");
   } else {
      fFirstV0 = kTRUE;
      fV0s     = new TClonesArray("MpdV0Track");
   }

   mngr->Register("MpdV0", "V0", fV0s, fWrite);
   if (!fPositiveDaughterCut) {
      LOG(ERROR) << "Lack cut for first daughter";
      return kFATAL;
   }
   if (!fNegativeDaughterCut) {
      LOG(ERROR) << "Lack cut for second daughter";
      return kFATAL;
   }
   if (!fCandicateCut) {
      LOG(ERROR) << "Lack candicate cut";
      return kFATAL;
   }
   if (fMpdEvent) {
      fFormat = EFormat::kDst;
   } else if (fMiniEvents) {
      fFormat = EFormat::kMiniDst;
   } else {
      LOG(ERROR) << "No MpdEvent/MiniDstEvent in tree";
      return kFATAL;
   }
   if (fDauMon1) {
      fDauMon1->Init();
      fDauMon2->Init();
   }
   if (fV0Mon) fV0Mon->Init();
   fInit = kTRUE;
   return kSUCCESS;
}

void MpdV0Finder::SetPositiveDaughterCut(const MpdV0DaughterCut &cut)
{
   if (fPositiveDaughterCut) delete fPositiveDaughterCut;
   fPositiveDaughterCut = (MpdV0DaughterCut *)cut.Clone();
}

void MpdV0Finder::SetNegativeDaughterCut(const MpdV0DaughterCut &cut)
{
   if (fNegativeDaughterCut) delete fNegativeDaughterCut;
   fNegativeDaughterCut = (MpdV0DaughterCut *)cut.Clone();
}

void MpdV0Finder::SetCandicateCut(const MpdV0CandidateCut &cut)
{
   if (fCandicateCut) delete fCandicateCut;
   fCandicateCut = (MpdV0CandidateCut *)cut.Clone();
}

MpdV0Finder &MpdV0Finder::operator=(const MpdV0Finder &other)
{
   if (this == &other) return *this;
   fFormat      = other.fFormat;
   fMpdEvent    = nullptr;
   fMiniEvents  = nullptr;
   fMiniTracks  = nullptr;
   fMiniTofData = nullptr;
   if (fV0s) {
      delete fV0s;
      fV0s = nullptr;
   }
   fInit                = kFALSE;
   fWrite               = other.fInit;
   fPositiveDaughterCut = nullptr;
   fNegativeDaughterCut = nullptr;
   fCandicateCut        = nullptr;
   if (other.fPositiveDaughterCut) {
      fPositiveDaughterCut = (MpdV0DaughterCut *)other.fPositiveDaughterCut->Clone();
   }
   if (other.fNegativeDaughterCut) {
      fNegativeDaughterCut = (MpdV0DaughterCut *)other.fNegativeDaughterCut->Clone();
   }
   if (other.fCandicateCut) {
      fCandicateCut = (MpdV0CandidateCut *)other.fCandicateCut->Clone();
   }

   if (other.fDauMon1) {
      fDauMon1 = other.fDauMon1->MakeCopy();
   }
   if (other.fDauMon2) {
      fDauMon2 = other.fDauMon1->MakeCopy();
   }
   if (other.fV0Mon) {
      fV0Mon = other.fV0Mon->MakeCopy();
   }
   return *this;
}

void MpdV0Finder::SetDaugherMonitor(const MpdV0DaughterMonitor &mon)
{
   if (fDauMon1) {
      delete fDauMon1;
      delete fDauMon2;
   }
   fDauMon1 = mon.MakeCopy();
   fDauMon2 = mon.MakeCopy();
}

void MpdV0Finder::SetV0Monitor(const MpdV0CandidateMonitor &mon)
{
   if (fV0Mon) delete fV0Mon;
   fV0Mon = mon.MakeCopy();
}

void MpdV0Finder::SetDstData()
{
   if (fDauMon1) {
      fDauMon1->SetEventData(fMpdEvent);
      fDauMon2->SetEventData(fMpdEvent);
   }
   if (fV0Mon) fV0Mon->SetEventData(fMpdEvent);
   fNegativeDaughterCut->SetEventData(fMpdEvent);
   fPositiveDaughterCut->SetEventData(fMpdEvent);
   fCandicateCut->SetEventData(fMpdEvent);
}

void MpdV0Finder::SetMiniDstData()
{
   MpdMiniEvent *event = static_cast<MpdMiniEvent *>(fMiniEvents->UncheckedAt(0));
   if (fDauMon1) {
      fDauMon1->SetMiniEventData(event, fMiniTracks, fMiniTofData);
      fDauMon2->SetMiniEventData(event, fMiniTracks, fMiniTofData);
   }
   if (fV0Mon) {
      fV0Mon->SetMiniEventData(event, fMiniTracks, fMiniTofData);
   }

   fPositiveDaughterCut->SetMiniEventData(event, fMiniTracks, fMiniTofData);
   fNegativeDaughterCut->SetMiniEventData(event, fMiniTracks, fMiniTofData);
}

void MpdV0Finder::FinishTask()
{
   if (fDauMon1 != nullptr || fV0Mon != nullptr) {
      TDirectory *dir = (TDirectory *)gFile;
      dir->mkdir(GetName(), "V0Monitors");
      dir->cd(GetName());
      if (fDauMon1) {
         std::vector<TH1 *> vec1 = fDauMon1->GetHistograms("PosDau");
         std::vector<TH1 *> vec2 = fDauMon2->GetHistograms("NegDau");
         for (auto i : vec1) i->Write();
         for (auto i : vec2) i->Write();
      }
      if (fV0Mon) {
         std::vector<TH1 *> vec = fV0Mon->GetHistograms("V0Mon");
         for (auto i : vec) i->Write();
      }
      dir->cd();
   }
}

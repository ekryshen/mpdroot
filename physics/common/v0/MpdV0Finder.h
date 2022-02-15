/*
 * MpdV0FinderBasic.h
 *
 *  Created on: 16 gru 2021
 *      Author: Daniel Wielanek
 *      E-mail: daniel.wielanek at pw.edu.pl
 *      Warsaw University of Technology, Faculty of Physics
 *      JINR,  Laboratory of High Energy Physics
 */

#ifndef MPDROOT_PHYSICS_COMMON_V0_MPDV0FINDER_H_
#define MPDROOT_PHYSICS_COMMON_V0_MPDV0FINDER_H_

#include <FairTask.h>
#include <Rtypes.h>
#include <RtypesCore.h>
#include <TString.h>
#include <utility>
#include <vector>

class MpdV0CandidateMonitor;
class MpdV0DaughterMonitor;

class MpdV0CandidateCut;
class MpdV0DaughterCut;

/**
 * abstract class for finding vo particles
 */

class MpdEvent;

class MpdV0Finder : public FairTask {
protected:
   enum class EFormat { kDst, kMiniDst };
   Bool_t                                   fInit;
   Bool_t                                   fWrite;
   Bool_t                                   fFirstV0;
   Int_t                                    fPidDauPos, fPidDauNeg, fPidV0;
   EFormat                                  fFormat;
   MpdEvent *                               fMpdEvent;
   TClonesArray *                           fMiniEvents;
   TClonesArray *                           fMiniTracks;
   TClonesArray *                           fMiniTofData;
   TClonesArray *                           fV0s;
   std::vector<std::pair<TObject *, Int_t>> fPositiveDaughters;
   std::vector<std::pair<TObject *, Int_t>> fNegativeDaughters;
   virtual void                             ExecDst(Option_t *option)     = 0;
   virtual void                             ExecMiniDst(Option_t *option) = 0;
   virtual InitStatus                       Init();
   MpdV0DaughterCut *                       fPositiveDaughterCut;
   MpdV0DaughterCut *                       fNegativeDaughterCut;
   MpdV0CandidateCut *                      fCandicateCut;
   MpdV0DaughterMonitor *                   fDauMon1;
   MpdV0DaughterMonitor *                   fDauMon2;
   MpdV0CandidateMonitor *                  fV0Mon;
   void                                     SetDstData();
   void                                     SetMiniDstData();

public:
   MpdV0Finder(TString name = "LambdaFinder", Int_t pidMom = 3122, Int_t pidFirstDau = 211, Int_t pidSecDau = 2212);
   void         SaveV0s(Bool_t write) { fWrite = write; };
   void         SetDaugherMonitor(const MpdV0DaughterMonitor &mon);
   void         SetV0Monitor(const MpdV0CandidateMonitor &mon);
   virtual void Exec(Option_t *option);
   virtual ~MpdV0Finder();
   void SetPositiveDaughterCut(const MpdV0DaughterCut &cut);
   void SetNegativeDaughterCut(const MpdV0DaughterCut &cut);
   void SetCandicateCut(const MpdV0CandidateCut &cut);
   MpdV0Finder(const MpdV0Finder &other);
   MpdV0Finder &operator=(const MpdV0Finder &other);
   virtual void FinishTask();
   ClassDef(MpdV0Finder, 1)
};

#endif /* MPDROOT_PHYSICS_COMMON_V0_MPDV0FINDER_H_ */

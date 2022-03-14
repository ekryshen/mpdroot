#ifndef LHE_HITS_MAKER_H
#define LHE_HITS_MAKER_H

// *************************************************************************
//  Author: Oleg Rogachevsky             e-mail: rogach@sunhe.jinr.ru
//
//  collect all hits from tracking detectors in form used for
//  finding and fitting in LHE program
//
// Created: 1-07-07
// Modified:
//
// *************************************************************************

#include "TClonesArray.h"

#include "FairTask.h"
#include "FairDetector.h"
#include "FairHit.h"

#include "TpcLheHit.h"
#include "TpcLheTrack.h"
#include "TpcLheTrackCuts.h"

class TpcLheHitsMaker : public FairTask {

protected:
   Int_t fNHit;
   Int_t fNTrack;
   Int_t fNofEvents;

   TClonesArray *fTpcPoints; //! pointers to TPC points
   TClonesArray *fTstPoints; //! pointers to TST points

   TClonesArray *fLheHits; //! pointers to LHE hits

   TClonesArray *fListMCtracks; //! pointers to MC tracks
   TClonesArray *fGeantTracks;  //! pointers to selected Geant tracks

   TpcLheTrackCuts *fCuts;

   TString fOption; //  options to choose particles

public:
   virtual InitStatus Init(); //
   virtual InitStatus ReInit();
   virtual void       Exec(Option_t *option); //
   virtual void       Finish();               //

   void Reset();
   void Register();

   TpcLheHitsMaker(const char *name, const char *title = "CBM Task");
   TpcLheHitsMaker();
   virtual ~TpcLheHitsMaker();

   TpcLheHit   *AddHit();
   TpcLheTrack *AddTrack(Int_t mT);

   void GetTstHits();
   void GetTpcHits();
   void SetOption(Option_t *option = " ")
   {
      fOption = option;
      fOption.ToLower();
   }
   void SetTrack(TpcLheHit *hit);
   void CheckTracks();
   void PrintTracks(Int_t n);
   // void GetStripPoints();

   ClassDef(TpcLheHitsMaker, 1) // TpcLheHitsMaker
};

#endif

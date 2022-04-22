#ifndef LHE_TRACK_FINDER_H
#define LHE_TRACK_FINDER_H

// *************************************************************************
//  Author: Oleg Rogachevsky             e-mail: rogach@sunhe.jinr.ru
//
//  track finder with conformal mapping
//
// Created: 1-07-07
// Modified:
//
// *************************************************************************

#include "TObject.h"
#include "TObjArray.h"
#include "TArrayI.h"
#include "TBenchmark.h"
#include "Riostream.h"

#include "FairTask.h"

#include "TpcLheHit.h"
#include "TpcLheCMPoint.h"
#include "TpcLheCMTrack.h"
#include "TpcLheTrackCuts.h"
#include "TpcLheSegments.h"
#include "TpcLheTrackFitter.h"
#include "lhe.h"

class TpcLheTrackFinder : public FairTask {

protected:
   enum Momentum { HIGH, INTERMED, LOW };

   TpcLheTrackFitter *fFitter;
   TpcLhePoint       *fVertex; // pointer to the vertex

   TClonesArray *fCMHits;   // Array of hits transformed with conformal mapping
   TObjArray    *fCMTracks; // Array of tracks

   TpcLheSegments *fSegment;  //
   TObjArray      *fSegments; // segments with hits

   TClonesArray *fLheHits;     // Array of event's hits
   TClonesArray *fGeantTracks; // Array of geant tracks
   TClonesArray *fFoundTracks; // Array of found tracks

   TString fOption; //  options for operation

   TBenchmark *fBench; // benchmark object (just for run-time measurements)
   Float_t     mTime;  // total time consumption

   // Cuts
   TpcLheTrackCuts *fTrackCuts; // cuts for tracks

public:
   TpcLheTrackFinder(); //
   TpcLheTrackFinder(const char *name, const char *title = "CBM Task");
   virtual ~TpcLheTrackFinder(); //

   virtual void       Exec(Option_t *option);
   virtual InitStatus Init();   //
   virtual void       Finish(); //

   void Reset();    //
   void Register(); //
   void SetOption(Option_t *option = " ")
   {
      fOption = option;
      fOption.ToLower();
   }

   // getters

   void GetPhiRange(TArrayI *c, TpcLheCMTrack *cmtrack);
   void GetThetaRange(TArrayI *c, TpcLheCMTrack *cmtrack);
   void CopyClones(TClonesArray *cl1, TClonesArray *cl2);

   TpcLhePoint *GetVertex() { return fVertex; }
   Int_t        GetNumberOfHits() { return fCMHits->GetEntriesFast(); }
   TObjArray   *GetTracks() { return fCMTracks; }

   // setters:

   // Tracking procedures:
   void   DoTracking();                    // tracking of main vertex tracks
   void   LoopOverHits();                  // loops over hits
   void   CreateTrack(TpcLheCMPoint *hit); // create track with start at hit
   void   RemoveTrack(TpcLheCMTrack *track);
   void   GetClosestHit(TpcLheCMTrack *tr, Int_t st, Bool_t back);
   Bool_t TrackExtension(TpcLheCMTrack *track); //
   void   AddTrackForFit(TpcLheCMTrack *track); //
   void   TrackingInfo(char *info);             // information about the tracking
   void   CheckClones();

   ClassDef(TpcLheTrackFinder, 1); //
};

//________________________________________________________________
inline void TpcLheTrackFinder::RemoveTrack(TpcLheCMTrack *track)
{
   //--- Removes track from TObjArray and release the points.

   //  track->ClearPoints();        // release points
   fCMTracks->Remove(track); // remove from TObjArray
   delete track;             // delete track
   track = NULL;
}

#endif

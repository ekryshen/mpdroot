// Author: Oleg Rogachevsky
// Update: 2009-10-07 17:54:37+0400
// Copyright: 2009 (C) MPD coll.

#ifndef ROOT_MpdFillDstTask
#define ROOT_MpdFillDstTask

#include "MpdEvent.h"
#include "MpdPid.h"

#include "FairMCEventHeader.h"
#include "FairTask.h"

#include "TH1F.h"
#include "TH2F.h"

class MpdFillDstTask : public FairTask {

private:
   MpdEvent          *fEvent;
   TClonesArray      *fKFTracks, *fKFEctTracks;     // array of kalman filter tracks
   TClonesArray      *fMCTracks;                    // array of MC events to write
   TClonesArray      *fGenTracks;                   // array of MC events to write
   TClonesArray      *fTpcHits;                     // map with TPC hits
   FairMCEventHeader *fMCEventHeader;               // MC event header
   TClonesArray      *fTofMatching, *fEtofMatching; // tof information
   TClonesArray      *fVertex;                      // AZ

   Bool_t        fZdcSkeletonesSaved; // EL  ( flag of saved skeletones)
   TH2F         *fHistZdc1En;         // EL  ( empty skeletones of ZDC X-Y energy maps)
   TH2F         *fHistZdc2En;         // EL
   TClonesArray *fELossZdc1Histo;     // EL  ( arrays to fill empty skeletones)
   TClonesArray *fELossZdc2Histo;     // EL
   TClonesArray *fELossZdc1Value;     // EL  ( sum of signals at ZDC)
   TClonesArray *fELossZdc2Value;     // EL

   // Histograms of MC tracks
   TH1F    *fhTrackMotherId;
   TH1F    *fhTrackPrimaryPDG;
   TH2F    *fhTrackVertex;
   TH2F    *fhTruthVertex;
   MpdPid  *fPID;
   Int_t    fSharedHitArraySize;
   Short_t *fSharedHitArray; //[fTempArraySize]
   void     FillTrackDCA(MpdTrack *track, TVector3 *recoVertex, TVector3 *mcVertex);
   void     FillTrackPID(MpdTrack *track);
   void     FillTrackTpcHits(Int_t index, MpdTrack *track);
   void     CalculateSharedArrayMap();

public:
   MpdFillDstTask(const char *name = "MpdFillDstTask", const char *title = "MPD Task");
   virtual ~MpdFillDstTask(); // Destructor
   virtual void       SetPIDAlgorithm(MpdPid *pid) { fPID = pid; };
   virtual void       Exec(Option_t *option);
   virtual InitStatus Init();   //
   virtual void       Finish(); //

   void Reset(); //
   void SetOption(Option_t *option = " ")
   {
      fOption = option;
      fOption.ToLower();
   }
   // KG Delete MC tracks being outside the MPD - not implemented correct!!!
   void CleanMC();

   // MpdEvent *AddEvent(Option_t * option);
   MpdTrack *AddPrimaryTrack();

   ClassDef(MpdFillDstTask, 0); // fill MpdDst branch
};

#endif

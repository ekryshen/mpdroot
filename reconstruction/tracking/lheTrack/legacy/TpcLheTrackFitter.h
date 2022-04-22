#ifndef LHE_TRACK_FITTER_H
#define LHE_TRACK_FITTER_H

// *************************************************************************
//  Author: Oleg Rogachevsky             e-mail: rogach@sunhe.jinr.ru
//
//  track fitter
//
// Created: 1-07-07
// Modified:
//
// *************************************************************************

#include "FairTask.h"
#include "FairField.h"
//#include "FairTrackParH.h"

#include "TpcLheTrack.h"
#include "TpcLheTrackCuts.h"
#include "TpcLheCMTrack.h"
#include "TpcLheCMPoint.h"

class TH2F;
class TGeoTrack;

class TpcLheTrackFitter : public FairTask {

protected:
   FairField    *fMagField;  //
   TClonesArray *fTpcTracks; //
   TClonesArray *fTpcHits;   //!

   TClonesArray *fTpcTrCand;
   TClonesArray *fTpcTrFit;

   TH2F *fXYG, *fXYF, *fYZG, *fYZF, *fXZG, *fXZF, *fPXG, *fPXF, *fPYG, *fPYF, *fPZG, *fPZF;

   TClonesArray *fTpcPoints; //!
   TGeoTrack    *fTrCan, *fTrFit;

   TList  *fHistList;
   TString fOption; //  options to choose branches

   TpcLheTrackCuts *fTrackCuts; // RO cuts for tracks

   static TpcLheTrackFitter *ftInstance;

public:
   virtual void       Exec(Option_t *option);
   virtual InitStatus Init(); //
   virtual void       Finish();

   TpcLheTrackFitter(const char *name, const char *title = "Fair Task");
   TpcLheTrackFitter();
   virtual ~TpcLheTrackFitter();

   Int_t DeepFit(TpcLheTrack *tr);
   Int_t CircleFit(TpcLheTrack *tr);
   Int_t HelixFit(TpcLheTrack *tr);
   //  void Yvone(TpcLheTrack *track, Double_t psib);
   void Info4Fit(TpcLheTrack *track);
   void SetOption(Option_t *option = " ")
   {
      fOption = option;
      fOption.ToLower();
   }
   // void SaveHistograms();
   static TpcLheTrackFitter *Instance();

   ClassDef(TpcLheTrackFitter, 1); // TpcLheTrackFitter
};

#endif


#ifndef MPD_TPCTRACKFOLLOW2SFT_H
#define MPD_TPCTRACKFOLLOW2SFT_H

//#include "MpdKalmanHitZ.h"

#include "FairTask.h"
//#include "TLorentzVector.h"
#include "TH1.h"

class MpdKalmanHitZ;
class MpdKalmanHit;
class MpdKalmanTrack;
// class MpdEtofPoint;
// class TpcLheHit;
class MpdEctKalmanTrack;
class TpcLheKalmanTrack;
class TClonesArray;
class TVector2;
class TArrayI;

class MpdTpcTrackFollow2Sft : public FairTask {
public:
   /** Constructor **/
   MpdTpcTrackFollow2Sft(const char *name = "MpdTpcTrackFollow2Sft", Int_t iVerbose = 1);

   /** Destructor **/
   virtual ~MpdTpcTrackFollow2Sft();

   /// * FairTask methods

   /** Intialisation at begin of run. To be implemented in the derived class.
    *@value  Success   If not kSUCCESS, task will be set inactive.
    **/
   InitStatus Init();

   /** Reinitialisation.
    *@value  Success   If not kSUCCESS, task will be set inactive.
    **/
   InitStatus ReInit();

   /** Intialise parameter containers.
    **/
   void SetParContainers();

   void Exec(Option_t *option);

   /** Action after each event. **/
   void Finish();
   void Reset();

private:
   Int_t GetNofHitsInLayer(Int_t lay) { return (Int_t)fhLays->GetCellContent(lay + 1, 0); }
   Int_t GetHitsInLayer(Int_t lay) { return fLayPointers[lay]; }                ///< first index of hits in layer
   void  GetTrackSeeds(Int_t iPass);                                            // build track seeds
   void  Convert2ECT(const TpcLheKalmanTrack *trTPC, MpdEctKalmanTrack *trECT); // convert TPC to ECT track
   // void Convert2ECT(const MpdKalmanTrack *trTPC, MpdKalmanTrack *trECT); // convert TPC to ECT track

   void DoTracking(Int_t iPass); // run tracking

   void     FollowInTPC();                                                // follow tracks in TPC
   void     FollowInSFT();                                                // follow tracks in SFT
   void     MakeKalmanHitsSft();                                          // create Kalman hits for SFT
   Double_t EvalPt(const MpdKalmanHitZ *hit1, const MpdKalmanHitZ *hit2); // evaluate Pt
   void     EvalCovar(const MpdKalmanHitZ *hitOut, const MpdKalmanHitZ *hitIn,
                      MpdEctKalmanTrack *track);                                // eval. covariance
   Int_t    RunKalmanFilter(MpdEctKalmanTrack *track, Int_t layBeg);            ///< run Kalman filter
   Int_t    RunKalmanFilterTPC(MpdEctKalmanTrack *track, Int_t layBeg);         ///< run Kalman filter in TPC
   void     RemoveDoubles();                                                    ///< remove double tracks
   Int_t    NofCommonHits(MpdEctKalmanTrack *track, MpdEctKalmanTrack *track1); ///< get number of common hits
   void     AddHits();                                                          // add hit objects to tracks
   void     SetTrackID(MpdEctKalmanTrack *track);                               // set track ID from IDs of its hits
   void     StoreTracks();   // transfer tracks from fTrackCand to fTracks
   void     ExcludeHits(){}; // exclude used hits
   void     Write();
   void     Writedir2current(TObject *obj);
   Double_t Proxim(MpdKalmanTrack *track, const MpdKalmanHit *hit);    ///< adjust R-Phi coord. for continuity
   TVector2 GetDistance(MpdEctKalmanTrack *track, MpdKalmanHitZ *hit); ///< hit-track distance
   void     GoToVertex();                                              // propagate tracks to vertex
   Bool_t   PassInnerShell(MpdKalmanTrack *track);                     ///< propagate thru TPC inner shell

   Int_t       fNPass;
   TDirectory *fHistoDir;
   TH1F       *fhNBranches;

   TClonesArray *fSftPoints;   //! SFT points
   TClonesArray *fTpcPoints;   //! TPC points
   TClonesArray *fTpcHits;     //! TPC hits
   TClonesArray *fKHits;       //!< array of Kalman hits
   TClonesArray *fMCTracks;    //!< array of MC tracks
   TClonesArray *fEctTracks;   //!< array of ECT-TPC tracks
   TClonesArray *fTpcTracks;   //!< array of TPC tracks
   TClonesArray *fTracks;      //! SFT tracks
   Int_t        *fLayPointers; //!< locations of hits from different layers
   TH1F         *fhLays;       //!< histo with layer hit multiplicities

   // FairVertex *fPrimVtx;

private:
   // Some constants
   static const Double_t fgkChi2Cut; ///< max accepted Chi2 of hit for track

   ClassDef(MpdTpcTrackFollow2Sft, 1);
};

#endif

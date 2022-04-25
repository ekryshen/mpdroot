
#ifndef MPD_TRACKFINDERITS5SPD_H
#define MPD_TRACKFINDERITS5SPD_H

#include "MpdItsHit5spd.h"
#include "MpdItsKalmanTrack.h"
#include "FairTask.h"
//#include "TLorentzVector.h"
#include "TH1.h"
#include <map>

// class MpdKalmanStripHit;
class MpdKalmanHit;
class MpdKalmanTrack;
class MpdTpcKalmanFilter;
class TClonesArray;
// class TVector2;

class MpdTrackFinderIts5spd : public FairTask {
public:
   /** Constructor **/
   MpdTrackFinderIts5spd(Bool_t useVector = kFALSE, const char *name = "MpdTrackFinderIts5spd", Int_t iVerbose = 1);

   /** Destructor **/
   virtual ~MpdTrackFinderIts5spd();

   /// * FairTask methods

   /** Intialisation at begin of run. To be implemented in the derived class.
    *@return  Success   If not kSUCCESS, task will be set inactive.
    **/
   InitStatus Init();

   /** Reinitialisation.
    *@return  Success   If not kSUCCESS, task will be set inactive.
    **/
   InitStatus ReInit();

   /** Intialise parameter containers.
    **/
   void SetParContainers();

   void Exec(Option_t *option);

   /** Action after each event. **/
   void   Finish();
   void   Reset();
   void   FillGeoScheme();                                              // fill Kalman filter geometry manager
   Bool_t Refit(MpdItsKalmanTrack *track, Double_t mass, Int_t charge); // refit track

private:
   // Some constants
   static const Int_t fgkNlays = 5, fgkNlays2 = 10;

   Int_t GetNofHitsInLayer(Int_t lay) { return (Int_t)fhLays->GetBinContent(lay + 1, 0); }
   Int_t GetHitsInLayer(Int_t lay) { return fLayPointers[lay]; }     ///< first index of hits in layer
   void  GetTrackSeeds(Int_t iPass);                                 // build track seeds
   void  DoTracking(Int_t iPass);                                    // run tracking
   void  MakeKalmanHits();                                           // create Kalman hits
   Int_t RunKalmanFilterMod(MpdItsKalmanTrack *track, Int_t layBeg); ///< run Kalman filter - modular
   // Bool_t NavigateToLayer(Int_t lay, MpdItsKalmanTrack *curTr, MpdItsKalmanTrack *trackBr,
   //			 std::map<Int_t,Int_t>& trackBrM); ///< navigate track to layer
   Bool_t NavigateToLayer(Int_t lay, MpdItsKalmanTrack *curTr, std::map<Int_t, MpdItsKalmanTrack> &trackBr,
                          std::multimap<Int_t, Int_t> &hitsInWin);
   Int_t  TrackID(MpdKalmanHit *hit); ///< return track ID of the hit
   void   AddHits();                  // add hit objects to tracks
   void   StoreTracks();              // transfer tracks from fTrackCand to fTracks
   void   ExcludeHits(){};            // exclude used hits
   void   Write();
   void   Writedir2current(TObject *obj);
   // TVector2 GetDistance(TpcLheKalmanTrack *track, MpdKalmanStripHit *hit); ///< hit-track distance
   TVector2 GetDistance(MpdKalmanTrack *track, MpdKalmanHit *hit); ///< hit-track distance

   Int_t       fExact;
   Int_t       fNPass;
   TDirectory *fHistoDir;
   TH1F       *fhNBranches;

   TClonesArray *fItsPoints; //! ITS MC points
   TClonesArray *fItsHits;   //! ITS Hits
   TClonesArray *fKHits;     //!< array of Kalman hits
   TClonesArray *fTpcTracks; //!< array of TPC tracks
   TClonesArray *fMCTracks;  //!< array of MC tracks
   // TClonesArray *fSTSTrackMatch;
   TClonesArray                  *fTracks;              //! ITS tracks
   Int_t                         *fLayPointers;         //!< locations of hits from different layers
   TH1F                          *fhLays;               //!< histo with layer hit multiplicities
   Double_t                       fStereoA[2];          //!< stereo angles
   Double_t                       fDz[fgkNlays];        //!< layer half-lengths in Z
   Double_t                       fZmod[fgkNlays];      //!< module lengths in Z
   Double_t                       fRad[fgkNlays];       //!< layer radii (detector layers * 2 sides)
   Double_t                       fPipeR;               //!< beam pipe outer radius
   Int_t                          fNladders[fgkNlays2]; //!< numbers of ladders in layers
   Int_t                          fNsectors[fgkNlays2]; //!< numbers of sectors in layers
   Int_t                          fGeo;                 //!< geometry version
   Bool_t                         fUseVector;           //!< use Vector Finder Kalman Hits
   std::map<Double_t, Double_t>   fCables[fgkNlays];    //!< cable length (z1, z2) for layers 1-6
   MpdItsHit5spd                  fHitSts;              //!< STS hit - helper for detector numeration
   std::map<Int_t, Int_t>         fId2Id[fgkNlays2];    //!< STS hit ID to Kalman hit ID for layers 1-6
   std::multimap<Double_t, Int_t> fHitMapRphi[fgkNlays], fHitMapZ[fgkNlays];
   MpdTpcKalmanFilter            *fTpcKF; //!< TPC kalman filter

   // FairVertex *fPrimVtx;
   // FairStsKFTrackFitter fStsFitter;

private:
   // Some constants
   static const Double_t fgkChi2Cut; ///< max accepted Chi2 of hit for track

   ClassDef(MpdTrackFinderIts5spd, 0);
};

#endif

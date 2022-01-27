
#ifndef MPD_TRACKFINDERITS_H
#define MPD_TRACKFINDERITS_H

#include "MpdStsHit.h"
#include "FairTask.h"
//#include "TLorentzVector.h"
#include "TH1.h"
#include <map>

//class MpdKalmanStripHit;
class MpdKalmanHit;
class MpdKalmanTrack;
class MpdItsKalmanTrack;
class TClonesArray;
class TVector2;

class MpdTrackFinderIts :public FairTask
{
 public:

  /** Constructor **/
  MpdTrackFinderIts(const char *name="MpdTrackFinderIts", Int_t iVerbose = 1 );
  
  /** Destructor **/
  virtual ~MpdTrackFinderIts();
  
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

  void Exec(Option_t * option);

  /** Action after each event. **/
  void Finish();
  void Reset();
  void FillGeoScheme(); // fill Kalman filter geometry manager   

 private:

  Int_t GetNofHitsInLayer(Int_t lay) { return (Int_t)fhLays->GetCellContent(lay+1,0); }
  Int_t GetHitsInLayer(Int_t lay) { return fLayPointers[lay]; } ///< first index of hits in layer
  void GetTrackSeeds(Int_t iPass); // build track seeds                       
  void DoTracking(Int_t iPass); // run tracking               
  void MakeKalmanHits(); // create Kalman hits   
  Int_t RunKalmanFilterCyl(MpdItsKalmanTrack *track, Int_t layBeg); ///< run Kalman filter - cylindrical
  Int_t RunKalmanFilterMod(MpdItsKalmanTrack *track, Int_t layBeg); ///< run Kalman filter - modular
  Bool_t NavigateToLayer(Int_t lay, MpdItsKalmanTrack *curTr, MpdItsKalmanTrack *trackBr, 
			 std::map<Int_t,Int_t>& trackBrM); ///< navigate track to layer
  Int_t TrackID(MpdKalmanHit *hit); ///< return track ID of the hit
  void AddHits(); // add hit objects to tracks                               
  void StoreTracks(); // transfer tracks from fTrackCand to fTracks          
  void ExcludeHits() { }; // exclude used hits
  void Write();
  void Writedir2current( TObject *obj );
  //TVector2 GetDistance(TpcLheKalmanTrack *track, MpdKalmanStripHit *hit); ///< hit-track distance
  TVector2 GetDistance(MpdKalmanTrack *track, MpdKalmanHit *hit); ///< hit-track distance

  Int_t fNPass; 
  TDirectory *fHistoDir;
  TH1F *fhNBranches;

  TClonesArray *fItsPoints; //! ITS MC points
  TClonesArray *fItsHits; //! ITS Hits
  TClonesArray *fKHits; //!< array of Kalman hits
  TClonesArray *fTpcTracks; //!< array of TPC tracks
  TClonesArray *fEctTracks; //!< array of ECT tracks
  TClonesArray *fMCTracks; //!< array of MC tracks
  //TClonesArray *fSTSTrackMatch;
  TClonesArray *fTracks; //! ITS tracks
  Int_t *fLayPointers; //!< locations of hits from different layers             
  TH1F *fhLays; //!< histo with layer hit multiplicities                        
  Double_t fStereoA[2]; //!< stereo angles
  Double_t fDz[5]; //!< layer half-lengths in Z
  Double_t fZmod[5]; //!< module lengths in Z
  Double_t fRad[10]; //!< layer radii (detector layers * 2 sides)
  Double_t fPipeR; //!< beam pipe outer radius
  Int_t fNladders[10]; //!< numbers of ladders in layers
  Int_t fNsectors[10]; //!< numbers of sectors in layers
  Int_t fGeo; //!< geometry version 
  std::map<Double_t,Double_t> fCables[5]; //!< cable length (z1, z2) for layers 1-4
  MpdStsHit fHitSts; //!< STS hit - helper for detector numeration 
  std::map<Int_t,Int_t> fId2Id[10]; //!< STS hit ID to Kalman hit ID for layers 1-4

  //FairVertex *fPrimVtx;
  //FairStsKFTrackFitter fStsFitter;

 private:
  // Some constants                                                             
  static const Double_t fgkChi2Cut; ///< max accepted Chi2 of hit for track     

  ClassDef(MpdTrackFinderIts,1);

};

#endif

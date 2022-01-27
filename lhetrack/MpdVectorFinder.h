#ifndef MPD_VECTORFINDER_H 
#define MPD_VECTORFINDER_H 

/// \ingroup rec
/// \class MpdVectorFinder
/// \brief Track finder in MPD Inner Tracking System (ITS) using Vector Finder approach
///
/// \author Dmitri Zinchenko, Alexander Zinchenko, LHEP JINR Dubna

#include "MpdItsHit5spd.h"
#include "MpdVector.h"
#include "FairTask.h"
#include "TH1.h"
#include "TObject.h" 
#include "TClonesArray.h" 
#include "TVector2.h" 
#include <map>

class MpdKalmanHit;
class MpdKalmanTrack;
class MpdItsKalmanTrack;
class MpdTpcKalmanFilter;

class MpdVectorFinder :public FairTask 
{
 public:

  /** Constructor **/
  MpdVectorFinder(const char *name="MpdVectorFinder",  Bool_t sa = 0, Int_t iVerbose = 1);

  /** Destructor **/
  virtual ~MpdVectorFinder(); 
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
  Bool_t GoToBeamLine(MpdItsKalmanTrack *track); // go to beam line
  Bool_t Refit(MpdItsKalmanTrack *track, Double_t mass, Int_t charge);

 private:

  Int_t GetNofHitsInLayer(Int_t lay) { return (Int_t)fhLays->GetBinContent(lay+1,0); } /// was Cell
  Int_t GetHitsInLayer(Int_t lay) { return fLayPointers[lay]; } ///< first index of hits in layer
  void GetTrackSeeds(Int_t iPass); // build track seeds                       
  void DoTracking(Int_t iPass); // run tracking               
  ///void MakeKalmanHits(); // create Kalman hits //ms // deprecated? added to Build2DHits
  //  Double_t EvalPt(const MpdKalmanHit *hit1, const MpdKalmanHit *hit2); // evaluate Pt
  
  //AZ Double_t EvalPt(const MpdVector *track1, const MpdVector *track2);
  Double_t EvalPt(const MpdVector *track1, const MpdVector *track2, const MpdVector *track0 = NULL, Double_t *circle = NULL);
  ///
  Double_t EvalPhiPt(const MpdVector *track1, const MpdVector *track2);
  ///

  //void EvalCovar(const MpdCellTrack *hitOut, const MpdCellTrack *hitIn, Double_t *parOut, Double_t *parIn, MpdItsKalmanTrack *track, const MpdCellTrack *track1); ///< evaluate covar. matrix OLD Version
  //void EvalCovar(Double_t *parOut, Double_t *parIn, MpdItsKalmanTrack *track, const MpdVector *track4); ///< evaluate covar. matrix new version track4-> 1 layer OLD
  void EvalCovar(Double_t *parOut, Double_t *parIn, MpdItsKalmanTrack *track, const MpdVector *track1); ///< evaluate covar. matrix work version coment ot 25.12
 // void MakeKalmanHits2D();//create 2D hits //ms
  void Build2DHits(); // create 2D Hits // 06.05
  void MakeTrackCandidates(Int_t iPass);
  void ExtendCellTracks(Int_t iPass);//ms 03.03
  void ExtendSecondaryTracks(Int_t iPass); /// 26.2.20
  Bool_t PropagateHelix(Double_t rlay, TVector3 pos1, Double_t *circle, Double_t &phi); // propagate track as helix
  Int_t RunKalmanFilterCell(MpdItsKalmanTrack *track, Int_t iPass); ///< run Kalman filter (fitter) for cell track

  //Int_t TrackID(MpdKalmanHit *hit, Int_t indx = 0); ///< return track ID of the hit //ms
  Int_t TrackID(MpdKalmanHit *hit, Int_t indx = 0); ///< return track ID of the hit
  void AddHits(); // add hit objects to tracks
  void SetTrackID(MpdItsKalmanTrack* track); ///< set track ID as ID of majority of its hits
  Int_t GetHitID(MpdKalmanHit *hit); ///< get hit ID from MCPoint ID
  void StoreTracks(); // transfer tracks from fTrackCand to fTracks          
  void ExcludeHits();
  void Write();
  void Writedir2current( TObject *obj );
  void RemoveDoubles(); // remove double tracks 10.02
  TVector2 GetDistance(MpdKalmanTrack *track, MpdKalmanHit *hit); ///< hit-track distance //ms
  // Double_t Interp(Double_t angt, Int_t choice = 0, Int_t lay = 4); ///< parabolic interpolation
  Double_t Interp(Double_t angt, Int_t choice = 0, Int_t lay = 0); ///< parabolic interpolation
  //Bool_t AreTracksDoubles(MpdCellTrack *tr1, MpdCellTrack *tr2);//12.02
  Bool_t AreTracksDoubles(MpdItsKalmanTrack *tr1, MpdItsKalmanTrack *tr2);
  void GetShortTracks(); 
    
 private:
  Int_t fNPass; 
  Int_t fExact; //!< exact ID match if != 0
  TDirectory *fHistoDir;
  TH1F *fhNBranches;
  
  Int_t fNTracks;             // number of found tracks 10.02
  TClonesArray *fItsPoints; //! ITS MC points
  TClonesArray *fItsHits; //! ITS Hits
  TClonesArray *fKHits1; //!< array of Kalman hits
  TClonesArray *fKHits[5]; //!< array of Cell Tracks
  TClonesArray *f2DHits[5]; //!< array of Cell Tracks (2D Hits)
  TClonesArray *fTpcTracks; //!< array of TPC tracks
  TClonesArray *fMCTracks; //!< array of MC tracks
  TClonesArray *fTrackCand; //!< array of track candidates
  TClonesArray *fTrackExact; //!< array of track candidates
  //TClonesArray *fSTSTrackMatch;
  TClonesArray *fTracks; //! ITS tracks
  Int_t *fLayPointers; //!< locations of hits from different layers             
  TH1F *fhLays; //!< histo with layer hit multiplicities                        
  Double_t fStereoA[2]; //!< stereo angles
  Double_t fDz[5]; //!< layer half-lengths in Z
  Double_t fZmod[5]; //!< module lengths in Z
  Double_t fRad[5]; //!< layer radii (detector layers)
  Double_t fPipeR; //!< beam pipe outer radius
  Int_t fNladders[5]; //!< numbers of ladders in layers
  Int_t fNsectors[5]; //!< numbers of sectors in layers
  Int_t fGeo; //!< geometry version 

  ///
  Bool_t fSa; //!< standalone its tracking (without its-tpc matching)

  std::map<Double_t,Double_t> fCables[10]; //!< cable length (z1, z2) for layers 1-4 /// was 8
  MpdItsHit5spd fHitSts; //!< STS hit - helper for detector numeration 
  std::map<Int_t,Int_t> fId2Id[5]; //!< STS hit ID to Kalman hit ID for layers 1-5
  std::map<TString,Int_t> fCellMap;
  //FairVertex *fPrimVtx;
  std::vector<std::map<Int_t,MpdVector*> > fVectorCands; //!< vector combinations built
  MpdTpcKalmanFilter *fTpcKF; //!< TPC kalman filter
  
 private:
  // Some constants
  static const Double_t fgkChi2Cut; ///< max accepted Chi2 of hit for track

  ClassDef(MpdVectorFinder,0); 
};

#endif

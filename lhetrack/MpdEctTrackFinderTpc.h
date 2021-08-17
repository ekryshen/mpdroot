
#ifndef MPD_ECTTRACKFINDERTPC_H
#define MPD_ECTTRACKFINDERTPC_H

#include "TH1.h"
#include "FairTask.h"
//#include "FairStsKFTrackFitter.h"
//#include "TLorentzVector.h"

class MpdKalmanHitZ;
class MpdKalmanHit;
class MpdKalmanTrack;
//class MpdStrawendcapPoint;
class MpdEctKalmanTrack;
class FairGeoTransform;
class TClonesArray;
class TVector2;

class MpdEctTrackFinderTpc :public FairTask
{
 public:

  /** Constructor **/
  MpdEctTrackFinderTpc(const char *name="MpdEctTrackFinderTpc", Int_t iVerbose = 1 );
  
  /** Destructor **/
  virtual ~MpdEctTrackFinderTpc();
  
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

  Bool_t Mirror() { return fMirror; } ///< mirror flag (if TRUE - use ECT mirror hits)

 private:
  // Some constants
  static const Double_t fgkChi2Cut; ///< max accepted Chi2 of hit for track
  static const Int_t fgkNlays = 60 * 2, fgkNtube = 700, fgkNlays2 = 60; // nLaysMax, nTubesMax

 private:

  void InitGeo(); ///< initialize detector geometry                                  
  Int_t GetNofHitsInLayer(Int_t lay) { return (Int_t)fhLays->GetCellContent(lay+1,0); }
  Int_t GetHitsInLayer(Int_t lay) { return fLayPointers[lay]; } ///< first index of hits in layer
  void GetTrackSeeds(Int_t iPass); // build track seeds 
  void DoTracking(Int_t iPass); // run tracking 
  void MakeKalmanHits(); // create Kalman hits   
  void PassWall(MpdEctKalmanTrack *track, Double_t thick); ///< pass TPC end-plate
  Int_t RunKalmanFilter(MpdEctKalmanTrack *track, Int_t layBeg); ///< run Kalman filter       
  void SelectTracks(Int_t iPass); ///< do track selection and compute shared hit multiplicities
  void RemoveDoubles(); ///< remove double tracks
  Int_t NofCommonHits(MpdEctKalmanTrack* track, MpdEctKalmanTrack* track1, Int_t dir = 1); ///< get number of common hits
  void AddHits(); // add hit objects to tracks 
  Int_t HitMotherId(MpdKalmanHit* hit, Int_t idM, Int_t &id1); // check if hit has the same mother ID as idM
  void StoreTracks(Int_t iPass); // transfer tracks from fTrackCand to fTracks 
  void ExcludeHits(); // exclude used hits
  void Write();
  void Writedir2current( TObject *obj );
  Double_t Proxim(Double_t phi0, Double_t phi); // adjust angle for continuity
  Double_t Proxim(MpdKalmanTrack *track, const MpdKalmanHit *hit); ///< adjust R-Phi coord. for continuity
  TVector2 GetDistance(MpdEctKalmanTrack *track, MpdKalmanHit *hit); ///< hit-track distance

  Int_t fDetType; // detector type: 1 = ECT, 0 = TGEM, 2 = DCH
  Int_t fNPass; 
  Bool_t fMirror; // kTRUE if use ECT mirror hits 
  Bool_t fExact; // exact match flag (if 1 exact match of IDs)
  //Int_t fNofLays; //!< actual number of detector layers
  Int_t fNofTubes; //!< actual number of tubes in a layer
  Double_t frMinMax[3]; //!< min / max radii of ECT and tube radius
  Double_t fZplanes[fgkNlays];
  TDirectory *fHistoDir;
  TH1F *fhNBranches;
  Double_t fZtpc; // Z-edge of TPC

  TClonesArray *fEctPoints; //! ECT MC points
  TClonesArray *fEctHits; //! ECT Hits
  TClonesArray *fKHits; //!< array of Kalman hits
  TClonesArray *fTpcTracks; //!< array of TPC tracks
  TClonesArray *fMCTracks; //!< array of MC tracks
  //TClonesArray *fSTSTrackMatch;
  TClonesArray *fTracks; //! ECT tracks
  TClonesArray *fTrackCand; //! ECT track candidates
  Int_t *fLayPointers; //!< locations of hits from different layers
  TH1F *fhLays; //!< histo with layer hit multiplicities 
  //std::set<Int_t> fTpcUsed; //!< used TPC tracks

  FairGeoTransform *fTransf[fgkNlays][fgkNtube];
  //FairVertex *fPrimVtx;
  //FairStsKFTrackFitter fStsFitter;

  ClassDef(MpdEctTrackFinderTpc,1);

};

#endif

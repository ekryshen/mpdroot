#ifndef MPDTPCKALMANFILTER_H
#define MPDTPCKALMANFILTER_H

/// \ingroup rec
/// \class MpdTpcKalmanFilter
/// \brief Kalman filter track reconstructor in MPD central detector
///
/// \author Alexander Zinchenko, LHEP JINR Dubna

//#include "TpcPadPlane.h"
#include "TpcPoint.h"
#include "MpdTpcHit.h"
#include "TpcSectorGeoAZ.h"
#include "MpdTpcKalmanTrack.h"
#include "MpdKalmanHit.h"

#include "FairTask.h"
#include "FairField.h"

#include <TH1F.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include "TClonesArray.h"

#include <map>
//#include <tuple>
#include <vector>

class MpdTpcKalmanFilter : public FairTask {

public:
   struct matrix4 {
      TMatrixD xextr, xfilt, jacob, wfilt;
      matrix4() { resize(); }
      matrix4(TMatrixD a, TMatrixD b, TMatrixD c, TMatrixD d)
      {
         resize();
         xextr = a;
         xfilt = b;
         jacob = c;
         wfilt = d;
      }
      void resize()
      {
         xextr.ResizeTo(5, 1);
         xfilt.ResizeTo(5, 1);
         jacob.ResizeTo(5, 5);
         wfilt.ResizeTo(5, 5);
      }
   };

   MpdTpcKalmanFilter(); ///< Default ctor
   MpdTpcKalmanFilter(BaseTpcSectorGeo &secGeo, const char *name, const char *title = "TPC Kalman filter"); ///< Ctor
   virtual ~MpdTpcKalmanFilter(); ///< Destructor

   virtual InitStatus Init();
   virtual InitStatus ReInit();
   virtual void       Exec(Option_t *option);
   virtual void       Finish();

   void Reset();
   void Register();

   // virtual InitStatus ReInit();
   Int_t         GetNofHitsInLayer(Int_t lay) { return (Int_t)fhLays->GetCellContent(lay + 1, 0); }
   Int_t         GetHitsInLayer(Int_t lay) { return fLayPointers[lay]; } ///< first index of hits in layer
   TClonesArray *GetHits() { return fHits; }                             ///< get array of all hits
   TClonesArray *GetTracks() { return fTracks; }                         ///< get array of tracks
   Int_t         RunKalmanFilter(MpdTpcKalmanTrack *track);              ///< run Kalman filter
   Int_t         GetParticleId(Int_t id);                                ///< particle ID for track id
   void          SetModular(Int_t modular) { fModular = modular; } ///< set != 0 if modular geom. of r/out chambers
   void          SetSectorGeo(BaseTpcSectorGeo &secGeo)
   {
      fSecGeo = dynamic_cast<TpcSectorGeoAZ *>(&secGeo);
      if (!fSecGeo) Fatal("MpdTpcKalmanFilter::MpdTpcKalmanFilter", " !!! Wrong geometry type !!! ");
   } ///< set sector geometry
   // Bool_t Refit(MpdKalmanTrack *track, Double_t mass = 0.13957, Int_t charge = 1, Bool_t skip = kFALSE, Int_t iDir =
   // 1, Bool_t exclude = kFALSE, std::map<Double_t,std::tuple<TMatrixD,TMatrixD,TMatrixD,TMatrixD> > *cache = NULL);
   // ///< refit track using its points for given particle mass and charge
   Bool_t Refit(
      MpdKalmanTrack *track, Double_t mass = 0.13957, Int_t charge = 1, Bool_t skip = kFALSE, Int_t iDir = 1,
      Bool_t                       exclude = kFALSE,
      std::map<Double_t, matrix4> *cache   = NULL); ///< refit track using its points for given particle mass and charge
   Int_t  GetHitID(MpdKalmanHit *hit);            // get hit ID from MC point ID
   void   UseTpcHit(Bool_t useMC = kTRUE) { fUseMCHit = useMC; }; // to use TpcHit branch instead of MpdTpcHit
   void   FillGeoScheme(); // fill Kalman filter geometry manager info (for modular geometry of r/out chambers)
   Bool_t IsTpcHit() const { return fUseMCHit; } // tracking with hits or clusters ?
   void   Smooth(MpdTpcKalmanTrack *track, std::vector<std::pair<Double_t, TMatrixD>> &vecSmooth);

private:
   // Some constants
   static const Int_t fgkLays = 150; // number of padrows (with some margin)
   static const Int_t fgkSecs = 24;  // number of sectors (with some margin)

private:
   virtual void SetParContainers();
   void         GetTrackSeeds(Int_t iPass); // build track seeds
   void         GetTrackSeedsEndCaps();     // build track seeds from end-cap regions
   void         DoTracking(Int_t iPass);    // run tracking
   Double_t     EvalPt(const MpdKalmanHit *hit1, const MpdKalmanHit *hit2, Double_t *params); // evaluate Pt
   // this depends on class TpcCluster which is not built void Cluster2KalmanHits(); // create Kalman hits from clusters
   void   MakeKalmanHits();                    // create Kalman hits
   void   MakeKalmanHitsModul();               // create Kalman hits for modular geom. of r/out chambers
   void   GoOutward(MpdTpcKalmanTrack *track); // propagate track outward
   Bool_t BackTrace(MpdTpcKalmanTrack *track, TMatrixDSym &weight, Int_t iDir = 1,
                    Bool_t corDedx = kFALSE);                         // propagate track thru found hits
   void   GoOut();                                                    // backpropagate tracks outward
   void   RemoveDoubles(TClonesArray *tracks);                        // remove double tracks
   Int_t  GetNofCommonHits(MpdKalmanTrack *tr1, MpdKalmanTrack *tr2); // get common hits
   // Are compared tracks similar each other (if common hits in 2 tracks is greater limit)
   Bool_t   AreTracksDoubles(MpdKalmanTrack *tr1, MpdKalmanTrack *tr2);
   Double_t CorrectForLoss(Double_t pt, Double_t the, Int_t id); // correct for dE loss in pipe
   // Bool_t SameOrigin(TpcLheHit *hit, Int_t idKF, Int_t *mcTracks); // check hit origin
   Bool_t SameOrigin(TpcPoint *hit, Int_t idKF, Int_t *mcTracks); // check hit origin
   void   StoreTracks();                                          // transfer tracks from fTrackCand to fTracks
   void   ExcludeHits();                                          // exclude used hits
   void   RemoveShorts();                                         // remove short tracks (Nhits < 4)
   void   GoToBeamLine();                                         // propagate tracks to the beam line
   void   AddHits(Int_t indx0 = -1);                              // add hit objects to tracks
   void   SetTrackID(MpdTpcKalmanTrack *track);                   // set track ID from IDs of its hits
   // TpcPoint* GetPoint(MpdKalmanHit *hit); // get MCPoint pointer for the hit
   MpdTpcHit *GetTpcHit(MpdKalmanHit *hit); // get TpcHit pointer for the Kalman hit
   Int_t      SectorNo(const char *cpath);  // extract sector number from TGeo path
   Double_t   CorrectForLoss(Double_t pt, Double_t the, Double_t mass, Int_t charge = 1); ///< energy loss correction
   Double_t   CorrectForLossFluct(Double_t pt, Double_t the, Double_t mass,
                                  Int_t charge = 1); ///< energy loss fluct. correction
   void       MergeTracks(Int_t ipass);              // merge tracks
   Double_t Interp2d(const Double_t *moms, const Double_t *thes, std::vector<std::vector<Double_t>> &sigmas, Double_t p,
                     Double_t dip); // 2-d linear interpolation

   Int_t         fNofEvents;                   // number of events processed
   Int_t         fNTracks;                     // number of found tracks
   Int_t         fNPass;                       // number of reco passes to run
   TClonesArray *fHits;                        // array of hits
   TClonesArray *fKHits;                       // array of Kalman hits
   TClonesArray *fTracks;                      // array of tracks
   TClonesArray *fTrackCand;                   // array of track candidates
   TClonesArray *fLHEtracks;                   // array of "LHE" MC tracks
   TClonesArray *fMCtracks;                    // array of MC tracks
   TClonesArray *fTpcPoints;                   // array of TPC points
   Int_t        *fLayPointers;                 //! locations of hits from different layers
   Int_t fLaySecBegPointers[fgkLays][fgkSecs]; //! locations of hits from different layers and sectors (first hit)
   Int_t fLaySecEndPointers[fgkLays][fgkSecs]; //! locations of hits from different layers and sectors (last hit)
   TH1F *fhLays;                               // histo with layer hit multiplicities
   // std::multiset<Int_t> fLayset;      // layer set of hits
   Double_t fVertZ;   // primary vertex position estimate
   Int_t    fZflag;   // primary vertex position estimate quality flag
   Int_t    fModular; // not equal 0 if modular geometry of readout chambers
   // const TpcPadPlane *fPadPlane;        //! pointer to pad plane
   TpcSectorGeoAZ              *fSecGeo;   //! pointer to sector geometry
   Bool_t                       fUseMCHit; // to use TpcHit branch (hit producer) instead of MpdTpc (clusters)
   std::map<Double_t, matrix4> *fCache;    // cached track parameters for smoother

private:
   // Some constants
   static const Double_t fgkChi2Cut; // max accepted Chi2 of hit for track

   ClassDef(MpdTpcKalmanFilter, 0);
};
#endif

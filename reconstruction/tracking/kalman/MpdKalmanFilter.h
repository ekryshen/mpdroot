#ifndef MPDKALMANFILTER_H
#define MPDKALMANFILTER_H

/// \ingroup rec
/// \class MpdKalmanFilter
/// \brief Kalman filter track reconstructor in the MPD detector
///
/// \author Alexander Zinchenko, LHEP JINR Dubna

#include "MpdKalmanHit.h"
#include "MpdKalmanTrack.h"
#include "MpdKalmanGeoScheme.h"

#include "FairTask.h"
#include "FairField.h"

#include <TMatrixD.h>
#include <TMatrixDSym.h>
#include "TClonesArray.h"

class MpdKalmanFilter : public FairTask {
public:
   // KG public cpnstructor for TStreamer reading
   MpdKalmanFilter(const char *name = "MpdKalmanFilter", const char *title = "MPD Task"); ///< Ctor

   static MpdKalmanFilter *Instance();                                                 ///< get singleton instance
   static MpdKalmanFilter *Instance(const char *name, const char *title = "MPD Task"); ///< get singleton instance

   virtual void Exec(Option_t *option);

   void Reset();
   void Register();

   MpdKalmanGeoScheme *GetGeo() { return fGeoScheme; } ///< pointer to geo manager
   Bool_t PropagateToHit(MpdKalmanTrack *track, const MpdKalmanHit *hit, Bool_t calcLeng = kTRUE, Bool_t local = kFALSE,
                         Double_t stepBack = -1.0); ///< propagate track
   Bool_t PropagateParamToHit(MpdKalmanTrack *track, const MpdKalmanHit *hit, Bool_t calcLeng = kTRUE,
                              Bool_t local = kFALSE, Double_t stepBack = -1.0); ///< propagate track parameters
   Double_t  FilterHit(MpdKalmanTrack *track, const MpdKalmanHit *hit, TMatrixDSym &pointWeight,
                       TMatrixD &paramTmp); ///< compute Chi2 from hit addition (filter hit)
   Double_t  FilterHitR(MpdKalmanTrack *track, const MpdKalmanHit *hit, TMatrixDSym &pointWeight,
                        TMatrixD &paramTmp); ///< compute Chi2 from hit addition for "barrel" tracks and hits
   Double_t  FilterHitZ(MpdKalmanTrack *track, const MpdKalmanHit *hit, TMatrixDSym &pointWeight,
                        TMatrixD &parFilt); ///< compute Chi2 from hit addition for "end-cap" tracks and hits
   Double_t  FilterStrip(MpdKalmanTrack *track, const MpdKalmanHit *hit, TMatrixDSym &pointWeight,
                         TMatrixD &parFilt); ///< Kalman filtering for "strip" hits
   Double_t  FilterStripLocal(MpdKalmanTrack *track, const MpdKalmanHit *hit, TMatrixDSym &pointWeight,
                              TMatrixD &parFilt,
                              Double_t &posNew); ///< Kalman filtering for "strip" hits in local coordinates
   void      MnvertLocal(Double_t *a, Int_t l, Int_t, Int_t n, Int_t &ifail);
   Double_t  Proxim(Double_t phi0, Double_t phi, Double_t scale = 1.0); // adjust angle for continuity
   Double_t  Proxim(MpdKalmanTrack *track, const MpdKalmanHit *hit,
                    Double_t scale = 1.0); ///< adjust R-Phi coord. for continuity
   Bool_t    PropagateParamR(MpdKalmanTrack *track, const MpdKalmanHit *hit, Bool_t calcLeng,
                             Double_t *vert = (Double_t *)0x0); ///< propagate params
   Bool_t    PropagateParamP(MpdKalmanTrack *track, const MpdKalmanHit *hit, Bool_t calcLeng, Bool_t local = kFALSE,
                             Double_t stepBack = -1.0); ///< propagate params
   Bool_t    PropagateParamP(MpdKalmanTrack *track, const Double_t *plane, Bool_t calcLeng, Bool_t local = kFALSE,
                             Double_t stepBack = -1.0); ///< propagate params to plane
   Bool_t    PropagateParamZ(MpdKalmanTrack *track, const MpdKalmanHit *hit, Bool_t calcLeng); // propagate params
   Bool_t    FindPca(MpdKalmanTrack *track, Double_t *vert);                 ///< find PCA w.r.t. given point
   void      Convert(const MpdKalmanTrack *track1, MpdKalmanTrack *track2);  ///< convert track of diff. types
   Bool_t    Refit(MpdKalmanTrack *track, Int_t iDir, Bool_t leng = kFALSE); ///< refit track using its points
   Double_t  Interp(Int_t nDim, const Double_t *mom, const Double_t *dedxm, Double_t p); ///< energy loss interpolation
   Double_t  Scattering(MpdKalmanTrack *track, Double_t dx, TString mass2 = "0.0194797849",
                        Int_t charge = 1); // compute mult. scat. angle
   Double_t  Scattering(MpdKalmanTrack *track, Double_t x0, Double_t step, TString mass2 = "0.0194797849",
                        Int_t charge = 1);                                      // compute mult. scat. angle
   void      SetGeantParamB(MpdKalmanTrack *track, Double_t *v3, Double_t dir); // fill Geant vector
   void      GetField(Double_t *position, Double_t *field);                     // get mag. field vector
   Int_t     IsNumer() const { return fNumer; }                                 // if != 0 - numerical propagator
   void      SetNumer(Int_t numer) { fNumer = numer; }                          // if != 0 - numerical propagator
   Double_t *ExtrapOneStep(MpdKalmanTrack *track, Double_t step,
                           Int_t flag = 0); // propagate thru one step (use cache if flag != 0)
   virtual InitStatus Init();
   TMatrixD          &GetJacob() { return fJacob; } // return Jacobian

protected:
   virtual InitStatus ReInit();
   virtual void       Finish();
   virtual ~MpdKalmanFilter(); ///< Destructor

private:
   // automatic deleting when application exit
   static void DestroyInstance()
   {
      if (fgKF) delete fgKF;
   }

   static MpdKalmanFilter *fgKF; //! pointer to Singleton instance, may be ! required to exclude recursion

   Bool_t AnalParamZ(MpdKalmanTrack *track, const MpdKalmanHit *hit, Bool_t calcLeng); // propagate params analytically
   Bool_t AnalParamX(MpdKalmanTrack *track, const MpdKalmanHit *hit, Bool_t calcLeng,
                     Bool_t local = kFALSE); // propagate params analytically
   Bool_t PropagateWeight(MpdKalmanTrack *track, const MpdKalmanHit *hit, Double_t sign = 1.,
                          Double_t stepBack = -1.0); // propagate weight
   Bool_t AnalyticJacob(MpdKalmanTrack *track, const MpdKalmanHit *hit,
                        TMatrixD &jacob); // analytical Jacobian computation
   Bool_t AnalyticJacobX(MpdKalmanTrack *track, const MpdKalmanHit *hit,
                         TMatrixD &jacob);                                   // analytical Jacobian computation
   void   SetGeantParamE(MpdKalmanTrack *track, Double_t *v3, Double_t dir); // fill Geant vector
   Double_t Proxim(const MpdKalmanHit *hit0, const MpdKalmanHit *hit,
                   Double_t scale = 1.0);                                          // adjust R-Phi coord. for continuity
   void     GoOutward(MpdKalmanTrack *track);                                      // propagate track outward
   void     BackTrace(MpdKalmanTrack *track, TMatrixDSym &weight);                 // propagate track thru found hits
   void     GetFromGeantParamB(MpdKalmanTrack *track, Double_t *v3, Double_t dir); // fill params
   void     GetFromGeantParamE(MpdKalmanTrack *track, Double_t *v3, Double_t dir); // fill params
   void     ExtrapOneStepRungekutta(Double_t charge, Double_t step, Double_t *vect,
                                    Double_t *vout); // RungeKutta track propagation
   void     ExtrapOneStepHelix(Double_t charge, Double_t step, Double_t *vect,
                               Double_t *vout); // helix track propagation
   void     ExtrapOneStepHelix3(Double_t charge, Double_t step, Double_t *vect,
                                Double_t *vout) const; // helix track propagation
   void     SetField(FairField *field);                // set mag. field

   // KG exlude to copy by TStreamer because of get in Init()
   FairField *fMagField; //! magnetic field

private:
   MpdKalmanGeoScheme *fGeoScheme;  // pointer to geometry manager
   Int_t               fNumer;      // if != 0 - numerical propagator
   Double_t            fVectorG[8]; // cache for propagator
   TMatrixD            fJacob;
   // Some constants
   static const Int_t    fgkTriesMax; // max number of attempts to find exact position
   static const Double_t fgkEpsilon;  // tracking precision (cm)

   ClassDef(MpdKalmanFilter, 0);
};
#endif

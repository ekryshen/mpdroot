#ifndef MPDKFV0FITTER_H
#define MPDKFV0FITTER_H

/// \ingroup rec
/// \class MpdKfV0Fitter
/// \brief Kalman filter V0-fitter for the MPD detector
///
/// \author Alexander Zinchenko, LHEP JINR Dubna

class MpdKalmanTrack;
class MpdTpcKalmanTrack;
class TVector3;
#include "FairTask.h"
#include <TMatrixD.h>

class TClonesArray;

class MpdKfV0Fitter : public FairTask {

public:
   static MpdKfV0Fitter *Instance();                                                  ///< get singleton instance
   static MpdKfV0Fitter *Instance(const char *name, const char *title = "FAIR Task"); ///< get singleton instance

   virtual void Exec(Option_t *option);

   void Reset();
   void Register();

   void     EvalVertex(MpdKalmanTrack *tr[2]);                                         // evaluate V0 for 2 tracks
   Double_t FindVertex(MpdKalmanTrack *track1, MpdKalmanTrack *track2, TVector3 &vtx); // find V0 for 2 tracks
   void     Smooth(MpdKalmanTrack *tr[2]); // smooth (update track momenta and track lengths)
   void ComputeAandB(TMatrixD &xk0, const MpdKalmanTrack *track, const MpdKalmanTrack &trackM, TMatrixD &a, TMatrixD &b,
                     TMatrixD &ck0);                                 // compute matrices of derivatives
   void Proxim(const MpdKalmanTrack &track0, MpdKalmanTrack &track); // adjust track params
   TMatrixD GetCovariance() { return fCovar; }                       // access to covariance matrix

protected:
   virtual InitStatus Init();
   virtual InitStatus ReInit();
   virtual void       Finish();
   virtual ~MpdKfV0Fitter(); ///< Destructor

private:
   MpdKfV0Fitter();                                                  ///< Default ctor
   MpdKfV0Fitter(const char *name, const char *title = "FAIR Task"); ///< Ctor
   // automatic deleting
   static void DestroyInstance()
   {
      if (fgV0) delete fgV0;
   }

private:
   static MpdKfV0Fitter *fgV0;
   Double_t              fieldConst; ///< mag. field constant
   Double_t              fVtx[3];    ///< V0 position
   TMatrixD              fCovar;     ///< covariance matrix

   ClassDef(MpdKfV0Fitter, 1);
};
#endif

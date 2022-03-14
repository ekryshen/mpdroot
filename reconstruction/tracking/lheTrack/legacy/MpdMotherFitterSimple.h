#ifndef MPDMOTHERFITTERSIMPLE_H
#define MPDMOTHERFITTERSIMPLE_H

/// \ingroup rec
/// \class MpdMotherFitterSimple
/// \brief Simple mother fitter for the MPD detector
/// \bried (using helices and straight lines)
///
/// \author Alexander Zinchenko, LHEP JINR Dubna

class MpdKalmanTrack;
#include "MpdHelix.h"
#include "FairTask.h"
#include <TLorentzVector.h>
#include <TMatrixD.h>
#include <TVector3.h>
#include <vector>
using std::vector;

class TClonesArray;

void FcnDist(Int_t &npar, Double_t *gin, Double_t &rf, Double_t *par,
             Int_t iflag); // compute distance between straight line and helix

class MpdMotherFitterSimple : public FairTask {

public:
   static MpdMotherFitterSimple *Instance(); ///< get singleton instance
   static MpdMotherFitterSimple *Instance(const char *name,
                                          const char *title = "FAIR Task"); ///< get singleton instance

   virtual void Exec(Option_t *option);

   void Reset();
   void Register();

   // virtual InitStatus ReInit();
   TLorentzVector  BuildMother(vector<Double_t> &mass, vector<MpdKalmanTrack *> &vDaughtTr,
                               vector<TLorentzVector *> &vDaughLv); // build mother from tracks and/or Lorentz vectors
   void            FindVertex(vector<MpdKalmanTrack *> vDaught);    // find vertex from daughter tracks
   const Double_t *GetVertex() const { return fVtx; }               // return vertex position
   const TVector3 &GetMomN() const { return fMomN; }                // return neutral particle momentum
   MpdHelix       *GetHelix(Int_t i) const { return fTrC[i]; }      // return track helix
   Double_t        GetDca() const { return fDist; }                 // return distance between particles
   Double_t        Phi(MpdHelix *helix) const;                      // return track Phi
   Double_t        Pt(MpdHelix *helix) const { return TMath::Abs(fieldConst / helix->curvature()); } // track Pt
   Double_t        FieldConst() const { return fieldConst; } // return field constant

protected:
   virtual InitStatus Init();
   virtual InitStatus ReInit();
   virtual void       Finish();
   virtual ~MpdMotherFitterSimple(); ///< Destructor

private:
   MpdMotherFitterSimple();                                                  ///< Default ctor
   MpdMotherFitterSimple(const char *name, const char *title = "FAIR Task"); ///< Ctor
   static void DestroyInstance()
   {
      if (fgMF) delete fgMF;
   } // automatic deleting
     // void EvalVertex(vector<MpdKalmanTrack*> vDaught); // evaluate vertex from daughter tracks
   // void Smooth(vector<MpdKalmanTrack*> vDaught); // smooth (update track momenta and track lengths)
   // void ComputeAandB(TMatrixD &xk0, const MpdKalmanTrack *track, const MpdKalmanTrack &trackM,
   //                 TMatrixD &a, TMatrixD &b, TMatrixD &ck0); // compute matrices of derivatives
   // void Proxim(const MpdKalmanTrack &track0, MpdKalmanTrack &track); // adjust track params

private:
   static MpdMotherFitterSimple *fgMF;
   Double_t                      fieldConst; //! mag. field constant
   Double_t                      fVtx[3];    //! decay vertex position
   // TMatrixD fCovar; //! covariance matrix
   MpdHelix *fTrC[3]; //! charged track helices
   TVector3  fMomN;   //! neutral track momentum
   Double_t  fDist;   //! distance between particles

   ClassDef(MpdMotherFitterSimple, 1);
};
#endif

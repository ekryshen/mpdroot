#ifndef MPDMOTHERFITTERPART_H
#define MPDMOTHERFITTERPART_H

/// \ingroup rec
/// \class MpdMotherFitterPart
/// \brief Kalman filter mother particle fitter for the MPD detector
/// \brief (using MpdParticle)
///
/// \author Alexander Zinchenko, LHEP JINR Dubna

class MpdKalmanTrack;
class MpdVertex;
class TVector3;
#include "MpdParticle.h"
#include "FairTask.h"
#include <TMatrixD.h>
#include <vector>
using std::vector;

class TClonesArray;

class MpdMotherFitterPart : public FairTask {

public:
   static MpdMotherFitterPart *Instance();                                                  ///< get singleton instance
   static MpdMotherFitterPart *Instance(const char *name, const char *title = "FAIR Task"); ///< get singleton instance

   virtual void Exec(Option_t *option);

   void Reset();
   void Register();

   // virtual InitStatus ReInit();
   // void WeightAtDca(MpdParticle *part, MpdKalmanTrack &track, Double_t *vert); // obtain MpdParticle weight at DCA
   Double_t  BuildMother(MpdParticle *mother, vector<MpdParticle *> &vDaught); // build mother from daughters
   Double_t  BuildMother(MpdParticle *mother, vector<MpdKalmanTrack *> &vTracks, vector<MpdParticle *> &vDaught);
   Double_t  FindVertex(vector<MpdParticle *> vDaught, TVector3 &vtx); // find vertex from daughter tracks
   void      ComputeAandB(const TMatrixD &xk0, const MpdParticle &part, TMatrixD &a, TMatrixD &b, TMatrixD &ck0,
                          Bool_t flag = kTRUE); // compute matrices of derivatives
   Double_t  FieldConst() const { return fieldConst; }
   TMatrixD  ComputeQmatr(vector<MpdParticle *> vDaught);
   void      ParamsAtDca(MpdParticle *part);
   TMatrixD &GetCovariance() { return fCovar; }                   // access to covariance matrix
   Double_t  Chi2Vertex(MpdParticle *part, const MpdVertex *vtx); ///< compute Chi2 w.r.t. vertex

protected:
   virtual InitStatus Init();
   virtual InitStatus ReInit();
   virtual void       Finish();
   virtual ~MpdMotherFitterPart(); ///< Destructor

private:
   MpdMotherFitterPart();                                                  ///< Default ctor
   MpdMotherFitterPart(const char *name, const char *title = "FAIR Task"); ///< Ctor
   static void DestroyInstance()
   {
      if (fgMF) delete fgMF;
   }                                               // automatic deleting
   void EvalVertex(vector<MpdParticle *> vDaught); // evaluate vertex from daughter tracks
   void Smooth(vector<MpdParticle *> vDaught);     // smooth (update track momenta and track lengths)
   void ComputeAB0(const TMatrixD &xk0, const MpdParticle &part, TMatrixD &a, TMatrixD &b, TMatrixD &ck0,
                   Bool_t flag);                                     // compute matrices of derivatives for neutrals
   void Proxim(const MpdKalmanTrack &track0, MpdKalmanTrack &track); // adjust track params

private:
   static MpdMotherFitterPart *fgMF;
   Double_t                    fieldConst; ///< mag. field constant
   Double_t                    fVtx[3];    ///< decay vertex position
   TMatrixD                    fCovar;     ///< covariance matrix

   ClassDef(MpdMotherFitterPart, 1);
};
#endif

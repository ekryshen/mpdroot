/*
 * MpdV0FinderKFPackage.h
 *
 *  Created on: 19 sty 2022
 *      Author: Daniel Wielanek
 *      E-mail: daniel.wielanek at pw.edu.pl
 *      Warsaw University of Technology, Faculty of Physics
 *      JINR,  Laboratory of High Energy Physics
 */

#ifndef MPDROOT_PHYSICS_COMMON_V0_MPDV0FINDERKFPACKAGE_H_
#define MPDROOT_PHYSICS_COMMON_V0_MPDV0FINDERKFPACKAGE_H_

#include <FairTask.h>
#include <RtypesCore.h>
#include <TMatrixDfwd.h>
#include <TMatrixDSymfwd.h>
#include <vector>

#include "MpdV0FinderBasic.h"
#include "KFPTrackVector.h"

class MpdTpcKalmanTrack;

class MpdMiniTrack;

class KFParticleTopoReconstructor;

class MpdV0FinderKFPackage : public MpdV0FinderBasic {
protected:
   KFParticleTopoReconstructor *fKFFinder;
   TClonesArray *               fMiniCovMatrix;
   TVector3                     fEventVertex;
   KFPTrackVector               fInputTracks;

   virtual void       ExecDst(Option_t *option);
   virtual void       ExecMiniDst(Option_t *option);
   virtual InitStatus Init();

   std::vector<float> GetCovMatrixMini(Int_t index, Double_t *newParams);
   /**
    * derived from P.Batyuk MpdKfParticleFinder
    */
   void KalmanToInnerTpc(MpdTpcKalmanTrack &in, MpdTpcKalmanTrack &out);
   /**
    * derived from P.Batyuk MpdKfParticleFinder
    */
   std::vector<float> DoErrorPropagationToXYZPxPyPz(TMatrixDSym *cov, TMatrixD *params, MpdTpcKalmanTrack *track);
   void               WriteCandidates();

public:
   MpdV0FinderKFPackage(Int_t pidMom = 3122, Int_t pidFirstDau = 211, Int_t pidSecDau = 2212);
   virtual ~MpdV0FinderKFPackage();
   ClassDef(MpdV0FinderKFPackage, 1)
};

#endif /* MPDROOT_PHYSICS_COMMON_V0_MPDV0FINDERKFPACKAGE_H_ */

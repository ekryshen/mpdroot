//---------------------------------------------------------------------------
#ifndef __MPD_TPC_HIT_PRODUCER_H
#define __MPD_TPC_HIT_PRODUCER_H

/// \ingroup tpc
/// \class MpdTpcHitProducer
/// \brief Hit producer in MPD TPC
///
/// \author Alexander Zinchenko, LHEP JINR Dubna

#include "MpdTpcHit.h"

#include "FairTask.h"

#include <TVector3.h>
#include "TClonesArray.h"
#include "TpcSectorGeoAZ.h"

//---------------------------------------------------------------------------
class MpdTpcHitProducer : public FairTask {

public:
   MpdTpcHitProducer(BaseTpcSectorGeo &secGeo);
   ~MpdTpcHitProducer();

   virtual InitStatus Init();
   virtual void       Exec(Option_t *opt);

   void SetModular(Int_t modular) { fModular = modular; } ///< use simple or modular r/out chamber geometry
   void SetReal(Int_t real) { fReal = real; }             ///< if !=0 apply realistic effects
   void SetPersistance(Bool_t choice = kTRUE) { fPersistance = choice; } ///< set persistance flag
   void SetInput(TClonesArray *mcPoints) { fPointArray = mcPoints; }     // to run on the fly

private:
   virtual void SetParContainers();
   void         ExecModular(); ///< emulate geometry of readout chambers
   void         ExecReal();    ///< emulate realistic effects (resolution and 2-hit resolution)
   Bool_t       TwoHitRes(MpdTpcHit *hit0, MpdTpcHit *hit);          ///< apply 2-hit resolution
   void         ResVsAngle();                                        ///< apply resolution vs angle
   Double_t     Proxim(const MpdTpcHit *hit0, const MpdTpcHit *hit); // adjust R-Phi coord. for continuity
   MpdTpcHit *AddRawHit(Int_t hitIndx, Int_t detUID, const TVector3 &posHit, const TVector3 &posHitErr, Int_t pointIndx,
                        Int_t trackIndx);
   Bool_t     Interpolate(Int_t np, Int_t &ibeg, Double_t *yp, Double_t *xp, Double_t *zp, Double_t y0, Double_t &dir,
                          Double_t &xhit, Double_t &zhit);
   void       ExecNew(); ///< sectored and layered sensitive volume

   TClonesArray *fPointArray; // Input array of TpcPoints

   TClonesArray *fHitArray; // Output array of TpcHits

   TpcSectorGeoAZ *fSecGeo;
   Int_t           fModular;     // not equal 0 if modular structure of r/out chambers
   Int_t           fReal;        // not equal 0 if realistic effects (resolution and 2-hit resolution)
   Double_t        fZtpc;        // TPC half-length
   Bool_t          fPersistance; // to store or not hits

   ClassDef(MpdTpcHitProducer, 1);
};
//---------------------------------------------------------------------------
#endif // #ifndef __MPD_TPC_HIT_PRODUCER_H

// -------------------------------------------------------------------------
// -----                 MpdEmcHitCreator header file                 -----
// -------------------------------------------------------------------------

#ifndef MPDEMCHITCREATION_H
#define MPDEMCHITCREATION_H 1

#include <iostream>
#include "FairTask.h"
#include "MpdEmcHit.h"
#include "MpdEmcGeoParams.h"

class TClonesArray;

class MpdEmcHitCreation : public FairTask {
public:
   /** Default constructor **/
   MpdEmcHitCreation();

   /** Destructor **/
   ~MpdEmcHitCreation();

   /** Virtual method Init **/
   virtual InitStatus Init();

   /** Virtual method Exec **/
   virtual void Exec(Option_t *opt);
   void virtual Finish();

private:
   /** Input array of MpdEmcPoints **/
   TClonesArray *fPointArray;
   /** Input array of MCTracks **/
   TClonesArray *fMcTrackArray;

   /** Output array of MpdEmcHit **/
   TClonesArray    *fDigiArray;
   MpdEmcHit       *SearchHit(UInt_t detID);
   MpdEmcGeoParams *fGeoPar;

   ClassDef(MpdEmcHitCreation, 2);
};

#endif

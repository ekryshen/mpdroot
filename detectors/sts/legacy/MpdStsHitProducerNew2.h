// -------------------------------------------------------------------------
// -----                MpdStsHitProducerNew2 header file              -----
// ----                           18.04.2012                           -----
// -------------------------------------------------------------------------

#ifndef CBMHYPHITPRODUCERNEW2_H
#define CBMHYPHITPRODUCERNEW2_H 1

#include <map>
#include <iostream>
#include "FairTask.h"
#include "TH1F.h"
#include "TFile.h"
#include "MpdStsHit.h"
#include "TVector3.h"

class TClonesArray;
class TObjectArray;

class MpdStsHitProducerNew2 : public FairTask {

public:
   /** Default constructor **/
   MpdStsHitProducerNew2(const char *fileGeo);

   /** Destructor **/
   ~MpdStsHitProducerNew2();

   /** Virtual method Init **/
   virtual InitStatus Init();

   void SetParContainers(); // get par. containers

   /** Virtual method Exec **/
   virtual void Exec(Option_t *opt);

   MpdStsHit *AddHit(Int_t trackID, Int_t detID, Int_t side);
   void       CreateStructure();

private:
   /** Input array of MpdStsPoints **/
   TClonesArray *fPointArray;

   /** Output array of MpdStsHit **/
   TClonesArray *fDigiArray;

   TObjArray *fVolumeArray;

   /** Geo file to use **/
   TString fFileGeo;
   Float_t eneThr;

   typedef std::map<Int_t, Float_t> mapper;

   mapper emcX, emcY, emcZ, emcTheta, emcPhi, emcTau;

   ClassDef(MpdStsHitProducerNew2, 1);
};

#endif

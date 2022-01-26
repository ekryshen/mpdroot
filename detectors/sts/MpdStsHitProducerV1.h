// -------------------------------------------------------------------------
// -----                 MpdStsHitproducerV1 header file                 -----
// -------------------------------------------------------------------------

#ifndef CBMHYPHITPRODUCERV1_H
#define CBMHYPHITPRODUCERV1_H

#include <map>
#include <iostream>
#include "FairTask.h"
#include "TH1F.h"
#include "TFile.h"
#include "MpdStsHit.h"
#include "TVector3.h"

class TClonesArray;
class TObjectArray;

class MpdStsHitProducerV1 : public FairTask {

public:
   /** Default constructor **/
   MpdStsHitProducerV1(const char *fileGeo);

   /** Destructor **/
   ~MpdStsHitProducerV1();

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
   /*   map<Int_t, Float_t> emcX; */
   /*   map<Int_t, Float_t> emcY; */
   /*   map<Int_t, Float_t> emcZ; */
   /*   map<Int_t, Float_t> emcTheta; */
   /*   map<Int_t, Float_t> emcPhi; */
   /*   map<Int_t, Float_t> emcTau; */

   ClassDef(MpdStsHitProducerV1, 1);
};

#endif // #ifndef CBMHYPHITPRODUCERV1_H

// -------------------------------------------------------------------------
// -----                MpdItsHitProducer5spd header file              -----
// ----                           13.12.2016                           -----
// -------------------------------------------------------------------------

#ifndef MPDITSHITProducer5spd_H
#define MPDITSHITProducer5spd_H 1

#include "MpdItsHit5spd.h"
#include "FairTask.h"
#include "TH1F.h"
#include "TFile.h"
#include <TString.h>
#include "TVector3.h"
#include <iostream>
#include <map>

class TClonesArray;

class MpdItsHitProducer5spd : public FairTask {

public:
   /** Default constructor **/
   MpdItsHitProducer5spd();
   // MpdItsHitProducer5spd(const char* fileGeo);

   /** Destructor **/
   virtual ~MpdItsHitProducer5spd();

   /** Virtual method Init **/
   virtual InitStatus Init();

   void SetParContainers(); // get par. containers

   /** Virtual method Exec **/
   virtual void Exec(Option_t *opt);
   virtual void Finish();

   MpdItsHit5spd *AddHit(Int_t trackID, Int_t detID);
   void           CreateStructure();

private:
   /** Input array of MpdStsPoints **/
   TClonesArray *fPointArray;

   /** Output array of MpdStsHit **/
   TClonesArray *fDigiArray;

   /** Geo file to use **/
   TString fFileGeo;
   Float_t eneThr;

   // typedef std::map<Int_t, Float_t> mapper;

   // mapper emcX, emcY, emcZ, emcTheta, emcPhi, emcTau;

   ClassDef(MpdItsHitProducer5spd, 0);
};

#endif

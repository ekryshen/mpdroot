// -------------------------------------------------------------------------
// -----                 MpdBbcHitproducer header file                 -----
// -------------------------------------------------------------------------

#ifndef CBMHYPHITPRODUCER_H
#define CBMHYPHITPRODUCER_H 1

#include <map>
#include <iostream>
#include "FairTask.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TList.h"

#include "TFile.h"
#include "MpdBbcHit.h"
#include "TVector3.h"

class TClonesArray;
class TObjectArray;

class MpdBbcHitProducer : public FairTask {

public:
   /** Default constructor **/
   MpdBbcHitProducer(const char *fileGeo);

   /** Destructor **/
   ~MpdBbcHitProducer();

   /** Virtual method Init **/
   virtual InitStatus Init();

   /** Virtual method Exec **/
   virtual void Exec(Option_t *opt);

   MpdBbcHit *AddHit(Int_t trackID, Int_t detID, Float_t energy);
   void       CreateStructure();

   void virtual FinishTask();
   void virtual Finish();
   void MakeHists();

private:
   /** Input array of MpdBbcPoints **/
   TClonesArray *fPointArray;

   /** Output array of MpdBbcHit **/
   TClonesArray *fDigiArray;

   TObjArray *fVolumeArray;

   /** Geo file to use **/
   TString fFileGeo;
   Float_t eneThr;

   //_____ Histograms_____________
   TList *hlist;

   TH1F *fZ;
   TH1F *fR;

   TH2F *fXY;

   TH2F *fXY_small;
   TH2F *fXY_big;

   TH2F *fRphi;

   typedef std::map<Int_t, Float_t> mapper;

   mapper emcX, emcY, emcZ, emcTheta, emcPhi, emcTau;
   /*   map<Int_t, Float_t> emcX; */
   /*   map<Int_t, Float_t> emcY; */
   /*   map<Int_t, Float_t> emcZ; */
   /*   map<Int_t, Float_t> emcTheta; */
   /*   map<Int_t, Float_t> emcPhi; */
   /*   map<Int_t, Float_t> emcTau; */

   ClassDef(MpdBbcHitProducer, 1);
};

#endif

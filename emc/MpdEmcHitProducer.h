// -------------------------------------------------------------------------
// -----                 MpdEmcHitproducer header file                 -----
// -------------------------------------------------------------------------

#ifndef CBMHYPHITPRODUCER_H
#define CBMHYPHITPRODUCER_H 1


#include <map>
#include <iostream>
#include "FairTask.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TList.h"
//*
#include "TNtuple.h"
#include "TTree.h"
//*
#include "TFile.h"
#include "MpdEmcHit.h"
#include "TVector3.h"

class TClonesArray;
class TObjectArray;

class MpdEmcHitProducer : public FairTask
{

 public:

  /** Default constructor **/  
  MpdEmcHitProducer(const char* fileGeo);


  /** Destructor **/
  ~MpdEmcHitProducer();


  /** Virtual method Init **/
  virtual InitStatus Init();


  /** Virtual method Exec **/
  virtual void Exec(Option_t* opt);

  MpdEmcHit* AddHit(Int_t trackID, Int_t detID, Float_t energy);
  void CreateStructure();


  void virtual FinishTask();
  void virtual Finish();
  void MakeHists();
  

 private: 
   
  virtual void SetParContainers();

  /** Input array of MpdEmcPoints **/
  TClonesArray* fPointArray;

  /** Output array of MpdEmcHit **/
  TClonesArray* fDigiArray;  

  TObjArray *fVolumeArray;
 
  /** Geo file to use **/
  TString fFileGeo; 
  Float_t eneThr;
//_______Ntuple______________
  //TNtuple* nt;
 // TNtuple* ntxyz;
  //TNtuple* nl;
  

 
  TNtuple *nt;  
  //TNtuple *ntxyz;


  //_____ Histograms_____________ 
  TList *hlist; 
  TH1F *ffELoss;
  TH1F *fZ;
  TH1F *fR;
  TH2F *fZYp;
  TH2F *fXYp;
  TH2F *fZYm;
  TH2F *fXYm;
 


  TH2F *fXZ;

  TH2F *fRphi;
  Float_t fArgs[4];
  Float_t fxyz[9];
  
 
    
//  Float_t fELossXYZ[4];
  Float_t nnPoints;
  typedef std::map<Int_t, Float_t> mapper;
  
  mapper emcX, emcY, emcZ, emcTheta, emcPhi, emcTau;
/*   map<Int_t, Float_t> emcX; */
/*   map<Int_t, Float_t> emcY; */
/*   map<Int_t, Float_t> emcZ; */
/*   map<Int_t, Float_t> emcTheta; */
/*   map<Int_t, Float_t> emcPhi; */
/*   map<Int_t, Float_t> emcTau; */
  
  ClassDef(MpdEmcHitProducer,1);
  
};

#endif

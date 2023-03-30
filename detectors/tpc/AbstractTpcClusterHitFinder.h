//-----------------------------------------------------------
// Description:
//      AbstractTpcClusterHitFinder class is the interface
//      for the implementation of the specific class that
//      reads TPC digits, finds clusters and reconstructs hits
//
// Environment:
//      Software for the MPD Detector at NICA.
//
// Author:
//      Slavomir Hnatic LIT, JINR, Dubna - 5.2022
//
//-----------------------------------------------------------

#ifndef ABSTRACTTPCCLUSTERHITFINDER_HH
#define ABSTRACTTPCCLUSTERHITFINDER_HH

// FairRoot Class Headers ----------------
#include "FairRootManager.h"
#include "FairTask.h"

// MpdRoot Class Headers -----------------
#include "BaseTpcSectorGeo.h"

// ROOT Class Headers & Declarations ---------------
#include "TString.h"
#include "TClonesArray.h"

class AbstractTpcClusterHitFinder : public FairTask {
public:
   // Constructors/Destructors ---------
   // default constructor set to private, the subclasses must access/implement the custom one below
   AbstractTpcClusterHitFinder(BaseTpcSectorGeo &tpcGeo, const char *name, Bool_t val);
   virtual ~AbstractTpcClusterHitFinder();

   // FairRun methods left virtual for development purposes
   // specific implementation of this interface should not implement them
   virtual InitStatus Init();
   virtual void       Exec(Option_t *opt);
   void               Finish();

   // Methods to be implemented by specific algorithm
   virtual TString ModuleNameSuffix()   = 0;
   virtual void    TransformInputData() = 0;
   virtual void    FindClusters()       = 0;
   virtual void    FindHits()           = 0;

protected:
   FairRootManager  *ioman;
   BaseTpcSectorGeo *secGeo;

   // no getters/setters implemented, to avoid boilerplate code
   TClonesArray *digiArray;
   TClonesArray *clusArray;
   TClonesArray *hitArray;

   Bool_t persistence;

private:
   // disable default constructor
   AbstractTpcClusterHitFinder();

   InitStatus ReadInputRegisterOutput(); // read input digiArray, register outputs clusArray & hitArray
   void       ClearClustersHits();       // clears clusArray, hitArray

   ClassDef(AbstractTpcClusterHitFinder, 1);
};

#endif

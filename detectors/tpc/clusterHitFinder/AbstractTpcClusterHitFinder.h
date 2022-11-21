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

// ROOT Class Declarations ---------------
class TClonesArray;

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

   // Methods to be implemented by specific algorithm
   virtual void TransformInputData() = 0;
   virtual void FindClusters()       = 0;
   virtual void FindHits()           = 0;

protected:
   FairRootManager  *ioman;
   BaseTpcSectorGeo *secGeo;

   // there is no complex class hierarchy when deriving from AbstractClusterHitFinder class
   // no getters/setters are implemented, to avoid boilerplate code
   // it is up to developer to behave responsibly when accessing data
   TClonesArray *digiArray;
   TClonesArray *clusArray;
   TClonesArray *hitArray;

   // Geometry parameters
   // do not modify them in your implementation (const correctness not implemented)

   Bool_t persistence;

private:
   // disable default constructor
   AbstractTpcClusterHitFinder();

   InitStatus ReadGeometryParameters();  // read geometry parameters
   InitStatus ReadInputRegisterOutput(); // read input digiArray, register outputs clusArray & hitArray
   void       ClearClustersHits();       // clears clusArray, hitArray

   ClassDef(AbstractTpcClusterHitFinder, 1);
};

#endif

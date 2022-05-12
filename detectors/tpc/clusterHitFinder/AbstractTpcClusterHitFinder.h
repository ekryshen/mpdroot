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

// ROOT Class Declarations --
class TClonesArray;

class AbstractTpcClusterHitFinder : public FairTask {
public:
   // Constructors/Destructors ---------
   // default constructor set to private, the subclasses must access/implement the custom one below
   AbstractTpcClusterHitFinder(const char *name, Bool_t val);
   virtual ~AbstractTpcClusterHitFinder();

   // FairRun methods left virtual for development purposes
   // it is suggested to developer to not override them and to adapt to existing implementation
   virtual InitStatus Init();
   virtual void       Exec(Option_t *opt);
   virtual void       FinishTask(); 

   // Methods to be implemented by specific algorithm
   virtual void TransformInputData() = 0;
   virtual void FindClusters()       = 0;
   virtual void FindHits()           = 0;

protected:
   FairRootManager *ioman;

   // there is no complex class hierarchy when deriving from AbstractClusterHitFinder class
   // no getters/setters are implemented, to avoid boilerplate code
   // it is up to developer to behave responsibly when accessing data
   TClonesArray *digiArray;
   TClonesArray *clusArray;
   TClonesArray *hitArray;

   // Geometry parameters
   // do not modify them in your implementation (const correctness not implemented)
   Int_t  nSectors;  // number of sectors
   Int_t  nRows;     // number of rows
   Int_t *nPads;     // number of pads in a row
   Int_t  nTimeBins; // maximum number of time bins

private:
   // disable default constructor
   AbstractTpcClusterHitFinder();

   InitStatus ReadGeometryParameters();  // read geometry parameters
   InitStatus ReadInputRegisterOutput(); // read input digiArray, register outputs clusArray & hitArray
   void       ClearClustersHits();       // clears clusArray, hitArray

   Bool_t persistence;

   ClassDef(AbstractTpcClusterHitFinder, 1);
};

#endif

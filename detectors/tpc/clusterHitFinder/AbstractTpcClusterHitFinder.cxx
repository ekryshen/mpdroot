//-----------------------------------------------------------
// Description:
//      AbstractTpcClusterHitFinder interface source file
//
// Environment:
//      Software for the MPD Detector at NICA.
//
// Author:
//      Slavomir Hnatic LIT, JINR, Dubna - 5.2022
//
//-----------------------------------------------------------

#include "AbstractTpcClusterHitFinder.h"

// ROOT Headers --------
#include "TClonesArray.h"

#include <iostream>

//__________________________________________________________________________

AbstractTpcClusterHitFinder::AbstractTpcClusterHitFinder(const char *name, Bool_t val)
   : FairTask(name), persistence(val)
{
}

//__________________________________________________________________________

AbstractTpcClusterHitFinder::~AbstractTpcClusterHitFinder() {}

//__________________________________________________________________________

InitStatus AbstractTpcClusterHitFinder::Init()
{
   if (ReadGeometryParameters() == kERROR)
      return kERROR;
   else
      return ReadInputRegisterOutput();
}

//__________________________________________________________________________

void AbstractTpcClusterHitFinder::Exec(Option_t *opt)
{
   ClearClustersHits();
   TransformInputData();
   FindClusters();
   FindHits();
}

//__________________________________________________________________________

InitStatus AbstractTpcClusterHitFinder::ReadGeometryParameters()
{
   // get MpdTpc sector geometry
   secGeo = MpdTpcSectorGeo::Instance();
   if (!secGeo) {
      Error("AbstractTpcClusterHitFinder::readGeometryParameters", "MpdTpcSectorGeo not instantiated!");
      return kERROR;
   }

   nSectors  = secGeo->NofSectors() * 2;
   nRows     = secGeo->NofRows();
   nTimeBins = secGeo->GetNTimeBins();
   nPads     = secGeo->NPadsInRows();

   return kSUCCESS;
}

//__________________________________________________________________________

InitStatus AbstractTpcClusterHitFinder::ReadInputRegisterOutput()
{
   // get FairRoot Manager
   ioman = FairRootManager::Instance();
   if (!ioman) {
      Error("AbstractTpcClusterHitFinder::readInputRegisterOutput", "FairRootManager not instantiated!");
      return kERROR;
   }

   // Read input digi collection
   digiArray = (TClonesArray *)ioman->GetObject("MpdTpcDigit");
   if (!digiArray) {
      Error("AbstractTpcClusterHitFinder::readInputRegisterOutput", "Array of digits not found!");
      return kERROR;
   }

   // Get event header (for vertex Z-position)
   /////////////////////////////////////////////////////////////////////////
   // This is needed ONLY for AZ's MLEM algorithm,                        //
   // placed here to be able to get rid of Init() in his implementation   //
   // should be moved away from interface to his implementation of it     //
   ioman->GetObject("MCEventHeader");
   /////////////////////////////////////////////////////////////////////////

   // Create and register output arrays
   clusArray = new TClonesArray("MpdTpc2dCluster");
   ioman->Register("TpcCluster", "Tpc", clusArray, persistence);
   hitArray = new TClonesArray("MpdTpcHit");
   ioman->Register("TpcRecPoint", "Tpc", hitArray, persistence);

   return kSUCCESS;
}

//__________________________________________________________________________

void AbstractTpcClusterHitFinder::ClearClustersHits()
{
   clusArray->Delete();
   hitArray->Delete();
}

//__________________________________________________________________________

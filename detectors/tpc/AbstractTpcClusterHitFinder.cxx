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
#include "QA_TpcClusterHitFinder.h"
#include "FairLogger.h"

#include "AbstractTpc2dCluster.h"
#include "AbstractTpcHit.h"

// ROOT Headers --------
#include "TClonesArray.h"

#include <iostream>

//__________________________________________________________________________

AbstractTpcClusterHitFinder::AbstractTpcClusterHitFinder(BaseTpcSectorGeo &tpcGeo, BaseQA *qaObject, const char *name,
                                                         Bool_t val)
   : FairTask(name), persistence(val)
{
   secGeo = &tpcGeo;

   localQAptr = dynamic_cast<QA_TpcClusterHitFinder *>(qaObject);
   if (localQAptr) LOG(info) << "QA mode on: TpcClusterHitFinder";
}

//__________________________________________________________________________

AbstractTpcClusterHitFinder::~AbstractTpcClusterHitFinder() {}

//__________________________________________________________________________

InitStatus AbstractTpcClusterHitFinder::Init()
{
   if (localQAptr) localQAptr->moduleNameSuffix = ModuleNameSuffix();

   return ReadInputRegisterOutput();
}

//__________________________________________________________________________

void AbstractTpcClusterHitFinder::Exec(Option_t *opt)
{
   ClearClustersHits();
   TransformInputData();
   FindClusters();
   FindHits();

   if (localQAptr) GatherQAData();
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
   clusArray = (ModuleNameSuffix().EqualTo(TString("Mlem"))) ? new TClonesArray("MpdTpc2dCluster")
                                                             : new TClonesArray("Tpc2dClusterFast");
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

void AbstractTpcClusterHitFinder::GatherQAData()
{
   int currentEvent = ioman->GetEntryNr();
   localQAptr->eventNumber.push_back(currentEvent);
   localQAptr->eventClusArray.push_back(new TClonesArray());
   localQAptr->eventClusArray.back() = (TClonesArray *)clusArray->Clone();
   localQAptr->eventHitArray.push_back(new TClonesArray());
   localQAptr->eventHitArray.back() = (TClonesArray *)hitArray->Clone();
}

//__________________________________________________________________________

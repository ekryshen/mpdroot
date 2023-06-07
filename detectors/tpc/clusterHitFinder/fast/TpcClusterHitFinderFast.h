// Joint Institute for Nuclear Research (JINR), Dubna, Russia, Dec-2021
// Adapter for MPD TPC Hit Finder library TpcClustering.h
// author: A.V.Krylov JINR/LNP
// email: avkrylov@jinr.ru
//
// Interface port: Slavomir Hnatic MLIT JINR, June 2023
#ifndef TpcClusterHitFinderFast_HH
#define TpcClusterHitFinderFast_HH

// Interface Header
#include "AbstractTpcClusterHitFinder.h"

// Base Class Headers ----------------
#include "BaseTpcSectorGeo.h"

// Fast ClusterHitFinder library
#include "TpcClustering.h"

// Collaborating Class Declarations --
class TClonesArray;

class TpcClusterHitFinderFast : public AbstractTpcClusterHitFinder {
public:
   // Constructors/Destructors ---------
   TpcClusterHitFinderFast(BaseTpcSectorGeo &secGeo)
      : AbstractTpcClusterHitFinder(secGeo, "TPC Cluster finder Fast", kFALSE), _pSecGeo(&secGeo)
   {
   }
   ~TpcClusterHitFinderFast() {}

   // Interface Implementation
   TString ModuleNameSuffix() { return TString("Fast"); }
   void    TransformInputData() {}
   void    FindClusters() {}
   void    FindHits();

   // Modifiers -----------------------
   void SetPersistence(Bool_t opt = kTRUE) { fPersistence = opt; }

private:
   Bool_t            fPersistence;
   BaseTpcSectorGeo *_pSecGeo;

   void calcSector(const tpcClustering::EventClusters *pEventClusters);

   ClassDef(TpcClusterHitFinderFast, 1);
};

#endif

// Joint Institute for Nuclear Research (JINR), Dubna, Russia, Dec-2021/June 2023
// MPD TPC Hit Finder for NICA project
// Adapter of TpcClustering.h library to the abstract interface
// authors: A.V.Krylov JINR/LNP, S.Hnatic JINR/MLIT
// email: avkrylov@jinr.ru
//
#ifndef TpcClusterHitFinderFast_HH
#define TpcClusterHitFinderFast_HH

// Interface Header
#include "AbstractTpcClusterHitFinder.h"

// Base Class Headers ----------------
#include "FairTask.h"
#include <vector>

#include "TpcClustering.h"
#include "TpcSectorGeoAZ.h"

// Collaborating Class Declarations --
class TpcSectorGeoAZ;
class TClonesArray;

class TpcClusterHitFinderFast : public AbstractTpcClusterHitFinder {
public:
   // Constructors/Destructors ---------
   TpcClusterHitFinderFast(BaseTpcSectorGeo &secGeo);
   ~TpcClusterHitFinderFast();

   // Interface Implementation
   TString ModuleNameSuffix() { return TString("Fast"); }
   void    TransformInputData() {}
   void    FindClusters() {}
   void    FindHits();

   // Modifiers -----------------------
   void SetPersistence(Bool_t opt = kTRUE) { fPersistence = opt; }

   // Not needed yet methods
   inline void SetEvent(uint nEvent) { _nEvent = nEvent; }

   inline uint GetEvent() { return _nEvent; }
   inline uint GetSector() { return _nSector; }
   inline uint GetRow() { return _nRow; }
   inline uint GetTime() { return _nTime; }

private:
   int _nEvent;
   int _nSector;
   int _nRow;
   int _nTime;

   Bool_t          fPersistence;
   TpcSectorGeoAZ *_pSecGeo;

   void calcSector(const tpcClustering::EventClusters *pEventClusters);

   ClassDef(TpcClusterHitFinderFast, 1);
};

#endif

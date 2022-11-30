//
// Joint Institute for Nuclear Research (JINR), Dubna, Russia, Dec-2021
// MPD TPC Hit Finder for NICA project
// author: A.V.Krylov JINR/LNP
// email: avkrylov@jinr.ru
//
#ifndef TPCCLUSTERHITFINDERFAST_HH
#define TPCCLUSTERHITFINDERFAST_HH

// Base Class Headers ----------------
#include "FairTask.h"
#include <vector>

#include "TpcClustering.h"
#include "TpcSectorGeoAZ.h"

// Collaborating Class Declarations --
class TpcSectorGeoAZ;
class TClonesArray;

class TpcClusterHitFinderFast : public FairTask {
public:
   // Constructors/Destructors ---------
   TpcClusterHitFinderFast(BaseTpcSectorGeo &secGeo);
   ~TpcClusterHitFinderFast();

   // Modifiers -----------------------
   void SetPersistence(Bool_t opt = kTRUE) { fPersistence = opt; }

   virtual InitStatus Init();
   void               FinishTask();
   virtual void       Exec(Option_t *opt);

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

   TClonesArray *_pTpcBinDataAZlt; // input data
   TClonesArray *_pTpcHitFinder;   // output data

   // TH1 *ExecTime;
   // float fTime[1000];

   const char *_strFilePath;
   const char *_strInputBranchName;
   const char *_strOutputBranchName;

   Bool_t          fPersistence;
   TpcSectorGeoAZ *_pSecGeo;

   void calcSector(const tpcClustering::EventClusters *pEventClusters);

   ClassDef(TpcClusterHitFinderFast, 1);
};

#endif

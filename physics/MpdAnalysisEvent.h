#ifndef MPDANALYSISEVENT_H
#define MPDANALYSISEVENT_H

#include "FairEventHeader.h"
#include "FairMCEventHeader.h"
#include "MpdEvent.h"
#include "TClonesArray.h"
#include "MpdAnalysisEventPlane.h"

class MpdAnalysisEvent {

public:
   MpdAnalysisEvent() = default;
   virtual ~MpdAnalysisEvent(); // Destructor
   void Clear();

public:
   /// Branch name
   FairEventHeader   *fEventHeader   = nullptr; ///< EventHeader
   TClonesArray      *fVertex        = nullptr; ///< Vertex
   MpdEvent          *fMPDEvent      = nullptr; ///< MPDEvent
   FairMCEventHeader *fMCEventHeader = nullptr; ///< MCEventHeader
   TClonesArray      *fMCTrack       = nullptr; ///< MCTrack

   TObjArray    *fEMCCluster     = nullptr; ///< List of EMC clusters
   TClonesArray *fTPCKalmanTrack = nullptr; ///< TPCKalmanTrack
   TClonesArray *fTOFHit         = nullptr; ///< TOFHit
   TClonesArray *fTOFMatching    = nullptr; ///< TOFMatching
   TClonesArray *fZDCDigit       = nullptr; ///< ZDCDigi
   TClonesArray *fZDCEloss1Value = nullptr; ///< ElossZDC1Value
   TClonesArray *fZDCEloss2Value = nullptr; ///< ElossZDC2Value
   TClonesArray *fZDCEloss1Histo = nullptr; ///< ElossZDC1Histo
   TClonesArray *fZDCEloss2Histo = nullptr; ///< ElossZDC2Histo

   MpdAnalysisEventPlane fMpdEP; ///< Event plane info

   void  setCentrTPC(float cent) { centralityTPC = cent; }
   float getCentrTPC() { return centralityTPC; }

private:
   float centralityTPC = 0; // V Centrality

   ClassDef(MpdAnalysisEvent, 0);
};

#endif

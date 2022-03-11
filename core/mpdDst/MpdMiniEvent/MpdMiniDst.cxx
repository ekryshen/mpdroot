//
// MpdMiniDst holds pointers to TClonesArrays with all data
//

// MiniDst headers
#include "MpdMiniMessMgr.h"
#include "MpdMiniEvent.h"
#include "MpdMiniTrack.h"
#include "MpdMiniBTofHit.h"
#include "MpdMiniBTofPidTraits.h"
#include "MpdMiniBECalCluster.h"
#include "MpdMiniTrackCovMatrix.h"
#include "MpdMiniFHCalHit.h"
#include "MpdMiniMcEvent.h"
#include "MpdMiniMcTrack.h"
#include "MpdMiniDst.h" //MUST be the last one

TClonesArray **MpdMiniDst::miniArrays = 0;

//_________________
void MpdMiniDst::unset()
{
   miniArrays = 0;
}

//_________________
void MpdMiniDst::set(TClonesArray **theMiniArrays)
{
   miniArrays = theMiniArrays;
}

//_________________
void MpdMiniDst::printEvent() const
{
   LOG_INFO << "\n=========== Event header =============\n\n";
   event()->Print();
   LOG_INFO << "=====================================\n\n";
}

//_________________
void MpdMiniDst::printTracks()
{
   if (numberOfTracks() == 0) {
      LOG_INFO << "No tracks found!" << endm;
      return;
   }

   LOG_INFO << "\n+++++++++ track list ( " << numberOfTracks() << " entries )\n\n";
   for (UInt_t iTrk = 0; iTrk < numberOfTracks(); iTrk++) {
      LOG_INFO << "+++ track " << iTrk << "\n";
      track(iTrk)->Print();
      LOG_INFO << "\n";
   }

   LOG_INFO << endm;
}

//_________________
void MpdMiniDst::printBECalClusters()
{

   if (numberOfBECalClusters() == 0) {
      LOG_INFO << "No ECalClusters found!" << endm;
      return;
   }

   LOG_INFO << "\n+++++++++ ECalCluster list ( " << numberOfBECalClusters() << " entries )\n\n";
   for (UInt_t iEntry = 0; iEntry < numberOfBECalClusters(); iEntry++) {
      LOG_INFO << "+++ becalCluster " << iEntry << "\n";
      becalCluster(iEntry)->Print();
      LOG_INFO << "\n";
   }

   LOG_INFO << endm;
}

//_________________
void MpdMiniDst::printBTofHits()
{

   if (numberOfBTofHits() == 0) {
      LOG_INFO << "No BTofHit found!" << endm;
      return;
   }

   LOG_INFO << "\n+++++++++ BTof list ( " << numberOfBTofHits() << " entries )\n\n";
   for (UInt_t iEntry = 0; iEntry < numberOfBTofHits(); iEntry++) {
      LOG_INFO << "+++ btofHit " << iEntry << "\n";
      btofHit(iEntry)->Print();
      LOG_INFO << "\n";
   }

   LOG_INFO << endm;
}

//_________________
void MpdMiniDst::printBTofPidTraits()
{

   if (numberOfBTofPidTraits() == 0) {
      LOG_INFO << "No BTof pidTraits found!" << endm;
      return;
   }

   LOG_INFO << "\n+++++++++ BTof pidTraits list ( " << numberOfBTofPidTraits() << " entries )\n\n";
   for (UInt_t iEntry = 0; iEntry < numberOfBTofPidTraits(); iEntry++) {
      LOG_INFO << "+++ EmcPidTraits " << iEntry << "\n";
      btofPidTraits(iEntry)->Print();
      LOG_INFO << "\n";
   }

   LOG_INFO << endm;
}

//_________________
void MpdMiniDst::printTrackCovMatrices()
{

   if (numberOfTrackCovMatrices() == 0) {
      LOG_INFO << "No TrackCovMatrix found!" << endm;
      return;
   }

   LOG_INFO << "\n+++++++++ trackCovMatrix list ( " << numberOfTrackCovMatrices() << " entries )\n\n";
   for (UInt_t iEntry = 0; iEntry < numberOfTrackCovMatrices(); iEntry++) {
      LOG_INFO << "+++ trackCovMatrix " << iEntry << "\n";
      trackCovMatrix(iEntry)->Print();
      LOG_INFO << "\n";
   }

   LOG_INFO << endm;
}

//_________________
void MpdMiniDst::printFHCalHits()
{
   if (numberOfFHCalHits() == 0) {
      LOG_INFO << "No FHCal hits found!" << endm;
      return;
   }

   LOG_INFO << "\n+++++++++ fhcalHit list ( " << numberOfFHCalHits() << " entries )\n\n";
   for (UInt_t iEntry = 0; iEntry < numberOfFHCalHits(); iEntry++) {
      LOG_INFO << "+++ FHCal hit " << iEntry << "\n";
      fhcalHit(iEntry)->Print();
      LOG_INFO << "\n";
   }

   LOG_INFO << endm;
}

//_________________
void MpdMiniDst::printMcTracks()
{
   if (numberOfMcTracks() == 0) {
      LOG_INFO << "No MC tracks found!" << endm;
      return;
   }

   LOG_INFO << "\n+++++++++ MC track list ( " << numberOfMcTracks() << " entries )\n\n";
   for (UInt_t iTrk = 0; iTrk < numberOfMcTracks(); iTrk++) {
      LOG_INFO << "+++ MC track " << iTrk << "\n";
      mcTrack(iTrk)->Print();
      LOG_INFO << "\n";
   }

   LOG_INFO << endm;
}

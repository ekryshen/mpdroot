#include "MpdAnalysisEvent.h"
ClassImp(MpdAnalysisEvent);

MpdAnalysisEvent::~MpdAnalysisEvent()
{
   if (fEventHeader) {
      delete fEventHeader;
      fEventHeader = nullptr;
   }
   if (fTPCKalmanTrack) {
      delete fTPCKalmanTrack;
      fTPCKalmanTrack = nullptr;
   }
   if (fVertex) {
      delete fVertex;
      fVertex = nullptr;
   }
   if (fFfdHit) {
      delete fFfdHit;
      fFfdHit = nullptr;
   }
   if (fTOFHit) {
      delete fTOFHit;
      fTOFHit = nullptr;
   }
   if (fTOFMatching) {
      delete fTOFMatching;
      fTOFMatching = nullptr;
   }
   if (fEmcDigit) {
      delete fEmcDigit;
      fEmcDigit = nullptr;
   }
   if (fEMCCluster) {
      delete fEMCCluster;
      fEMCCluster = nullptr;
   }
   if (fZDCDigit) {
      delete fZDCDigit;
      fZDCDigit = nullptr;
   }
   if (fMCEventHeader) {
      delete fMCEventHeader;
      fMCEventHeader = nullptr;
   }
   if (fMCTrack) {
      delete fMCTrack;
      fMCTrack = nullptr;
   }
   if (fMPDEvent) {
      delete fMPDEvent;
      fMPDEvent = nullptr;
   }
   if (fV0) {
      delete fV0;
      fV0 = nullptr;
   }
}
void MpdAnalysisEvent::Clear()
{

   if (fEventHeader) fEventHeader->Clear("C");
   if (fTPCKalmanTrack) fTPCKalmanTrack->Delete();
   if (fVertex) fVertex->Clear("C");
   if (fFfdHit) fFfdHit->Clear("C");
   if (fTOFHit) fTOFHit->Clear("C");
   if (fTOFMatching) fTOFMatching->Clear("C");
   if (fEmcDigit) fEmcDigit->Clear("C");
   if (fEMCCluster) { // otherwise memry leaks
      fEMCCluster->Delete();
   }
   if (fZDCDigit) fZDCDigit->Delete();
   if (fMCEventHeader) fMCEventHeader->Clear("C");
   if (fMCTrack) fMCTrack->Clear("C");
   if (fMPDEvent) fMPDEvent->Clear("C");
   if (fV0) fV0->Clear("C");
}

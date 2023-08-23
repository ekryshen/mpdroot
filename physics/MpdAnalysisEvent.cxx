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

   if (fEventHeader) fEventHeader->Clear();
   if (fTPCKalmanTrack) fTPCKalmanTrack->Clear();
   if (fVertex) fVertex->Clear();
   if (fFfdHit) fFfdHit->Clear();
   if (fTOFHit) fTOFHit->Clear();
   if (fTOFMatching) fTOFMatching->Clear();
   if (fEmcDigit) fEmcDigit->Clear();
   if (fEMCCluster) fEMCCluster->Clear();
   if (fZDCDigit) fZDCDigit->Clear();
   if (fMCEventHeader) fMCEventHeader->Clear();
   if (fMCTrack) fMCTrack->Clear();
   if (fMPDEvent) fMPDEvent->Clear();
   if (fV0) fV0->Clear();
}

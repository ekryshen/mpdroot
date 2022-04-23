#include "MpdAnalysisEvent.h"
ClassImp(MpdAnalysisEvent);

MpdAnalysisEvent::~MpdAnalysisEvent()
{
   if (fEventHeader) {
      delete fEventHeader;
      fEventHeader = nullptr;
   }
   if (fVertex) {
      delete fVertex;
      fVertex = nullptr;
   }
   if (fMPDEvent) {
      delete fMPDEvent;
      fMPDEvent = nullptr;
   }
   if (fMCEventHeader) {
      delete fMCEventHeader;
      fMCEventHeader = nullptr;
   }
   if (fMCTrack) {
      delete fMCTrack;
      fMCTrack = nullptr;
   }
   if (fTPCKalmanTrack) {
      delete fTPCKalmanTrack;
      fTPCKalmanTrack = nullptr;
   }
   if (fTOFHit) {
      delete fTOFHit;
      fTOFHit = nullptr;
   }
   if (fTOFMatching) {
      delete fTOFMatching;
      fTOFMatching = nullptr;
   }
   if (fZDCDigit) {
      delete fZDCDigit;
      fZDCDigit = nullptr;
   }
   if (fZDCEloss1Value) {
      delete fZDCEloss1Value;
      fZDCEloss1Value = nullptr;
   }
   if (fZDCEloss2Value) {
      delete fZDCEloss2Value;
      fZDCEloss2Value = nullptr;
   }
   if (fZDCEloss1Histo) {
      delete fZDCEloss1Histo;
      fZDCEloss1Histo = nullptr;
   }
   if (fZDCEloss2Histo) {
      delete fZDCEloss2Histo;
      fZDCEloss2Histo = nullptr;
   }
}
void MpdAnalysisEvent::Clear()
{

   if (fEventHeader) fEventHeader->Clear();
   if (fVertex) fVertex->Clear();
   if (fMPDEvent) fMPDEvent->Clear();
   if (fMCEventHeader) fMCEventHeader->Clear();
   if (fMCTrack) fMCTrack->Clear();
   if (fTPCKalmanTrack) fTPCKalmanTrack->Clear();
   if (fTOFHit) fTOFHit->Clear();
   if (fTOFMatching) fTOFMatching->Clear();
   if (fZDCDigit) fZDCDigit->Clear();
   if (fZDCEloss1Value) fZDCEloss1Value->Clear();
   if (fZDCEloss2Value) fZDCEloss2Value->Clear();
   if (fZDCEloss1Histo) fZDCEloss1Histo->Clear();
   if (fZDCEloss2Histo) fZDCEloss2Histo->Clear();
}
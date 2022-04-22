//
// MpdMiniEvent: holds the information about the event
//

// C++ headers
#include <algorithm>
#include <limits>

// ROOT headers
#include "TMath.h"

// MiniDst headers
#include "MpdMiniMessMgr.h"
#include "MpdMiniEvent.h"

ClassImp(MpdMiniEvent);

   //_________________
   MpdMiniEvent::MpdMiniEvent()
   : TObject(), fRunId(0), fEventId(0), fFillId(0), fBunchCrossId(0), fBField(5.), fTime(0), fPrimaryVertexX(-999.),
     fPrimaryVertexY(-999.), fPrimaryVertexZ(-999.), fPrimaryVertexErrorX(0), fPrimaryVertexErrorY(0),
     fPrimaryVertexErrorZ(0), fNBECalMatch(0), fNBTOFMatch(0), fTriggerIds(), fRefMultNeg(0), fRefMultPos(0),
     fRefMult2NegEast(0), fRefMult2PosEast(0), fRefMult2NegWest(0), fRefMult2PosWest(0), fRefMultHalfNegEast(0),
     fRefMultHalfPosEast(0), fRefMultHalfNegWest(0), fRefMultHalfPosWest(0), fGRefMult(0), fNumberOfGlobalTracks(0),
     fBTofTrayMultiplicity(0), fNFfdHitsEast(0), fNFfdHitsWest(0), fVzFfd(-999.), fNTofT0(0), fFHCalX(0),
     fBackgroundRate(0), fFHCalEastRate(0), fFHCalWestRate(0), fFHCalEnergyDeposit{}
{
   // Default constructor
   if (!fTriggerIds.empty()) {
      fTriggerIds.clear();
   }
}

//_________________
MpdMiniEvent::MpdMiniEvent(const MpdMiniEvent &event) : TObject()
{

   // Copy constructor

   // Global info
   fRunId        = event.fRunId;
   fEventId      = event.fEventId;
   fFillId       = event.fFillId;
   fBunchCrossId = event.fBunchCrossId;
   fBField       = event.fBField;
   fTime         = event.fTime;

   // Primary vertex info
   fPrimaryVertexX      = event.fPrimaryVertexX;
   fPrimaryVertexY      = event.fPrimaryVertexY;
   fPrimaryVertexZ      = event.fPrimaryVertexZ;
   fPrimaryVertexErrorX = event.fPrimaryVertexErrorX;
   fPrimaryVertexErrorY = event.fPrimaryVertexErrorY;
   fPrimaryVertexErrorZ = event.fPrimaryVertexErrorZ;

   // Number of matched tracks to fast detectors
   fNBECalMatch = event.fNBECalMatch;
   fNBTOFMatch  = event.fNBTOFMatch;

   // Trigger ID collection
   fTriggerIds = event.fTriggerIds;

   // Reference multiplicities
   fRefMultNeg         = event.fRefMultNeg;
   fRefMultPos         = event.fRefMultPos;
   fRefMult2NegEast    = event.fRefMult2NegEast;
   fRefMult2PosEast    = event.fRefMult2PosEast;
   fRefMult2NegWest    = event.fRefMult2NegWest;
   fRefMult2PosWest    = event.fRefMult2PosWest;
   fRefMultHalfNegEast = event.fRefMultHalfNegEast;
   fRefMultHalfPosEast = event.fRefMultHalfPosEast;
   fRefMultHalfNegWest = event.fRefMultHalfNegWest;
   fRefMultHalfPosWest = event.fRefMultHalfPosWest;

   fGRefMult             = event.fGRefMult;
   fNumberOfGlobalTracks = event.fNumberOfGlobalTracks;

   // Hit mulitplicities
   fBTofTrayMultiplicity = event.fBTofTrayMultiplicity;

   // FFD info
   fNFfdHitsEast = event.fNFfdHitsEast;
   fNFfdHitsWest = event.fNFfdHitsWest;
   fNTofT0       = event.fNTofT0;
   fVzFfd        = event.fVzFfd;

   // Forward detector info
   fFHCalX         = event.fFHCalX;
   fBackgroundRate = event.fBackgroundRate;
   fFHCalEastRate  = event.fFHCalEastRate;
   fFHCalWestRate  = event.fFHCalWestRate;

   for (int iIter = 0; iIter < 90; iIter++) {
      fFHCalEnergyDeposit[iIter] = event.fFHCalEnergyDeposit[iIter];
   }
}

//_________________
MpdMiniEvent::~MpdMiniEvent()
{
   // Destructor
   /* empty */
}

//_________________
void MpdMiniEvent::Print(const Char_t *option __attribute__((unused))) const
{
   LOG_INFO << " year = " << year() << " day = " << day() << "\n"
            << " fill/run/event Id = " << fillId() << "/" << runId() << "/" << eventId() << "\n"

            << " vertex x = " << primaryVertex().X() << " vertex y = " << primaryVertex().Y()
            << " vertex z = " << primaryVertex().Z() << "\n"
            << " refMult = " << refMult() << " grefMult = " << grefMult() << "\n"
            << " nTofT0 = " << nTofT0() << " ffdVz = " << vzFfd() << endm;
}

//_________________
Int_t MpdMiniEvent::year() const
{
   return fRunId / 1000000 - 1 + 2000;
}

//_________________
Int_t MpdMiniEvent::day() const
{
   return (fRunId % 1000000) / 1000;
}

//_________________
Bool_t MpdMiniEvent::isTrigger(unsigned int id) const
{
   return std::find(fTriggerIds.begin(), fTriggerIds.end(), id) != fTriggerIds.end();
}

//_________________
void MpdMiniEvent::setTriggerId(UInt_t id)
{

   // If trigger list is not empty then loop over it
   // and check if the new trigger already in.
   if (!fTriggerIds.empty()) {

      // Assume that the new trigger is not in the list
      Bool_t isUsed = false;

      // Loop over the trigger list
      for (UInt_t iIter = 0; iIter < fTriggerIds.size(); iIter++) {

         // Compare triggers
         if (fTriggerIds.at(iIter) == id) {
            isUsed = true;
         }
      } //(unsigned int iIter=0; iIter<fTriggerIds.size(); iIter++)

      // If the trigger not in the list then add it
      if (!isUsed) {
         fTriggerIds.push_back(id);
      }
   } // if( !fTriggerIds.empty() )
   else {
      fTriggerIds.push_back(id);
   }
}

//_________________
void MpdMiniEvent::setTriggerIds(std::vector<unsigned int> newIds)
{

   // Protection: work only if input vector has entries
   if (!newIds.empty()) {

      // If trigger list is not empty then loop over it
      // and check if the new trigger already in.
      if (!fTriggerIds.empty()) {

         // For each entry in the input vector
         for (UInt_t iIter1 = 0; iIter1 < newIds.size(); iIter1++) {

            // Assume that the new trigger is not in the list
            Bool_t isUsed = false;

            // Loop over existing trigger list
            for (UInt_t iIter2 = 0; iIter2 < fTriggerIds.size(); iIter2++) {

               // Compare triggers
               if (fTriggerIds.at(iIter2) == newIds.at(iIter1)) {
                  isUsed = true;
               }
            } // for (unsigned int iIter2=0; iIter2<fTriggerIds.size(); iIter2++)

            // The entry is unique then add it to the list
            if (!isUsed) {
               fTriggerIds.push_back(newIds.at(iIter1));
            }

         } // for(unsigned int iIter1= 0; iIter1<newIds.size(); iIter1++)
      }    // if( !fTriggerIds.empty() )
      else {
         fTriggerIds = newIds;
      }
   } // if( !newIds.empty() )
}

//________________
void MpdMiniEvent::setFHCalEnergyDepositInModule(Int_t iModule, Float_t energyDeposit)
{
   if (iModule >= 0 && iModule <= 89) {
      fFHCalEnergyDeposit[iModule] = energyDeposit;
   } else {
      // Probably some information about this incident may be printed
   }
}

//________________
void MpdMiniEvent::setFHCalEastEnergyDepositInModule(Int_t iModule, Float_t energyDeposit)
{
   if (iModule >= 0 && iModule <= 44) {
      fFHCalEnergyDeposit[iModule] = energyDeposit;
   } else {
      // Probably some information about this incident may be printed
   }
}

//________________
void MpdMiniEvent::setFHCalWestEnergyDepositInModule(Int_t iModule, Float_t energyDeposit)
{
   if (iModule >= 0 && iModule <= 44) {
      fFHCalEnergyDeposit[45 + iModule] = energyDeposit;
   } else {
      // Probably some information about this incident may be printed
   }
}

//_________________
void MpdMiniEvent::setBunchId(Int_t id)
{
   if (id < 0) {
      LOG_INFO << "MpdMiniEvent::setBunchID() - negative bunch ID = " << id << endm;
   } else {
      fBunchCrossId =
         ((id > std::numeric_limits<unsigned short>::max()) ? std::numeric_limits<unsigned short>::max() : (UChar_t)id);
   }
}

//________________
Double_t MpdMiniEvent::fhcalModuleAziAngle(const Int_t iModule) const
{

   // Angle of the module
   Float_t phiAngle = -999.;

   // Direction of the x axis:
   //  1 - for negative eta
   // -1 - for positive eta
   Short_t  xAxisDirection = (iModule < 45) ? 1 : -1;
   Double_t x              = 0.;
   Double_t y              = 0.;
   Short_t  module         = 0;

   // Choose module number (0-44 = negative eta; 45-89 = positive eta)
   if (iModule >= 0 && iModule <= 44) {
      module = iModule;
   } else if (iModule >= 45 && iModule <= 89) {
      module = iModule - 45;
   } else {
      module = -1;
   }

   // Estimate angle based on the x and y BHCal module position and size
   if ((module >= 0) && (module <= 4)) {
      y        = 45.;
      x        = (module - 2) * 15.;
      phiAngle = TMath::ATan2(y, xAxisDirection * x);
   } else if ((module >= 5) && (module <= 39)) {
      y        = (3 - (module + 2) / 7) * 15;
      x        = (3 - (module + 2) % 7) * 15;
      phiAngle = TMath::ATan2(y, xAxisDirection * x);
   } else if ((module >= 40) && (module <= 44)) {
      y        = -45.;
      x        = (module - 42) * 15.;
      phiAngle = TMath::ATan2(y, xAxisDirection * x);
   } else {
      phiAngle = -999.;
   }

   return phiAngle;
}

//________________
Double_t MpdMiniEvent::fhcalEastModuleAziAngle(const Int_t iModule) const
{

   // Angle of the module
   Float_t phiAngle = -999.;

   // Direction of the x axis:
   //  1 - for negative eta
   // -1 - for positive eta
   Short_t  xAxisDirection = 1;
   Double_t x              = 0.;
   Double_t y              = 0.;
   Short_t  module         = 0;

   // Choose module number
   if (iModule >= 0 && iModule <= 44) {
      module = iModule;
   } else {
      module = -1;
   }

   // Estimate angle based on the x and y BHCal module position and size
   if ((module >= 0) && (module <= 4)) {
      y        = 45.;
      x        = (module - 2) * 15.;
      phiAngle = TMath::ATan2(y, xAxisDirection * x);
   } else if ((module >= 5) && (module <= 39)) {
      y        = (3 - (module + 2) / 7) * 15;
      x        = (3 - (module + 2) % 7) * 15;
      phiAngle = TMath::ATan2(y, xAxisDirection * x);
   } else if ((module >= 40) && (module <= 44)) {
      y        = -45.;
      x        = (module - 42) * 15.;
      phiAngle = TMath::ATan2(y, xAxisDirection * x);
   } else {
      phiAngle = -999.;
   }

   return phiAngle;
}

//________________
Double_t MpdMiniEvent::fhcalWestModuleAziAngle(const Int_t iModule) const
{

   // Angle of the module
   Float_t phiAngle = -999.;

   // Direction of the x axis:
   //  1 - for negative eta
   // -1 - for positive eta
   Short_t  xAxisDirection = -1;
   Double_t x              = 0.;
   Double_t y              = 0.;
   Short_t  module         = 0;

   // Choose module number
   if (iModule >= 0 && iModule <= 44) {
      module = iModule;
   } else {
      module = -1;
   }

   // Estimate angle based on the x and y BHCal module position and size
   if ((module >= 0) && (module <= 4)) {
      y        = 45.;
      x        = (module - 2) * 15.;
      phiAngle = TMath::ATan2(y, xAxisDirection * x);
   } else if ((module >= 5) && (module <= 39)) {
      y        = (3 - (module + 2) / 7) * 15;
      x        = (3 - (module + 2) % 7) * 15;
      phiAngle = TMath::ATan2(y, xAxisDirection * x);
   } else if ((module >= 40) && (module <= 44)) {
      y        = -45.;
      x        = (module - 42) * 15.;
      phiAngle = TMath::ATan2(y, xAxisDirection * x);
   } else {
      phiAngle = -999.;
   }

   return phiAngle;
}

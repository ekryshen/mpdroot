//--------------------------------------------------------------------------
//----                     MpdStsHit                                    ----
//--------------------------------------------------------------------------

#include <iostream>
using namespace std;

#include "MpdStsHit.h"
#include <TMath.h>
#include "TObject.h"

// -----   Default constructor   -------------------------------------------

MpdStsHit::MpdStsHit()
{
   fFlag = 1;
}

// -----   Standard constructor   ------------------------------------------

MpdStsHit::MpdStsHit(Int_t detID, TVector3 pos, TVector3 dpos, Double_t signal, Int_t index, Int_t flag)
   : FairHit(detID, pos, dpos, index)
{
   fFlag   = flag;
   fSignal = signal;
}

// -----   Constructor without flag  ------------------------------------------

MpdStsHit::MpdStsHit(Int_t detID, TVector3 pos, TVector3 dpos, Double_t signal, Int_t index)
   : FairHit(detID, pos, dpos, index)
{
   fSignal = signal;
}
//__________________________________________________________________________
/*
MpdStsHit::MpdStsHit (const MpdStsHit& hit)
  : FairHit(hit),
    fDetectorID(hit.fDetectorID),
    fHitType(hit.fHitType),
    fDist(hit.fDist)
{
}
*/
// -----   Destructor   ----------------------------------------------------
MpdStsHit::~MpdStsHit() {}

// -----   Compare  --------------------------------------------------------

Int_t MpdStsHit::Compare(const TObject *hit) const
{
   /// "Compare" function to sort in ascending order of the layer number

   MpdStsHit *kHit = (MpdStsHit *)hit;

   // Check layers
   if (Layer() > kHit->Layer())
      return 1;
   else if (Layer() < kHit->Layer())
      return -1;
   return 0;
}

// -----  Print  -----------------------------------------------------------
void MpdStsHit::Print(const Option_t *opt) const
{
   cout << "-I- MpdStsHit" << endl;
   cout << "    DetectorID: " << fDetectorID << endl;
   cout << "    Position: (" << fX << ", " << fY << ", " << fZ << ") cm" << endl;
   cout << "    Position error: (" << fDx << ", " << fDy << ", " << fDz << ") cm" << endl;
   cout << "    Flag: " << fFlag << endl;
}

// -------------------------------------------------------------------------
Int_t MpdStsHit::SetDetId(Int_t sectorType, Int_t layer, Int_t ladder, Int_t det, Int_t side)
{
   // Set detector ID - helper for the tracking

   if (sectorType < 0) {
      // Layout with sensors
      fDetectorID = side;
      fDetectorID |= ((layer - 1) << 11);
      fDetectorID |= (ladder << 6);
      fDetectorID |= (det << 1);
   } else {
      // Layout with sectors
      fDetectorID = side;                       // # of side      on bit 0
      fDetectorID |= (det << kSensorS);         // # of sector    on bit 1 - 5
      fDetectorID |= (ladder << kLadderS);      // # of ladder    on bit 6 - 10
      fDetectorID |= (layer << kLayerS);        // # of layer     on bit 11 - 13
      fDetectorID |= (sectorType << kSecTypeS); // type of sector (1,2 or 3) on bit 14-16
   }
   return fDetectorID;
}

ClassImp(MpdStsHit)

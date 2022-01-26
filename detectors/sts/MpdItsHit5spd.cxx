//--------------------------------------------------------------------------
//----                     MpdItsHit5spd                                    ----
//--------------------------------------------------------------------------
//---VK 13.12/2016

#include <iostream>
using namespace std;

#include "MpdItsHit5spd.h"
#include <TMath.h>
#include "TObject.h"

// -----   Default constructor   -------------------------------------------

MpdItsHit5spd::MpdItsHit5spd()
{
   fFlag = 1;
}

// -----   Standard constructor   ------------------------------------------

MpdItsHit5spd::MpdItsHit5spd(Int_t detID, TVector3 pos, TVector3 dpos, Double_t signal, Int_t index, Int_t flag)
   : FairHit(detID, pos, dpos, index)
{
   fFlag   = flag;
   fSignal = signal;
}

// -----   Constructor without flag  ------------------------------------------

MpdItsHit5spd::MpdItsHit5spd(Int_t detID, TVector3 pos, TVector3 dpos, Double_t signal, Int_t index)
   : FairHit(detID, pos, dpos, index)
{
   fSignal = signal;
}
//__________________________________________________________________________
/*
MpdItsHit_5spd::MpdItsHit_5spd (const MpdItsHit_5spd& hit)
  : FairHit(hit),
    fDetectorID(hit.fDetectorID),
    fHitType(hit.fHitType),
    fDist(hit.fDist)
{
}
*/
// -----   Destructor   ----------------------------------------------------
MpdItsHit5spd::~MpdItsHit5spd() {}

// -----   Compare  --------------------------------------------------------

Int_t MpdItsHit5spd::Compare(const TObject *hit) const
{
   /// "Compare" function to sort in ascending order of the layer number

   MpdItsHit5spd *kHit = (MpdItsHit5spd *)hit;

   // Check layers
   if (Layer() > kHit->Layer())
      return 1;
   else if (Layer() < kHit->Layer())
      return -1;
   return 0;
}

// -----  Print  -----------------------------------------------------------
void MpdItsHit5spd::Print(const Option_t *opt) const
{
   cout << "-I- MpdItsHit5spd" << endl;
   cout << "    DetectorID: " << fDetectorID << endl;
   cout << "    Position: (" << fX << ", " << fY << ", " << fZ << ") cm" << endl;
   cout << "    Position error: (" << fDx << ", " << fDy << ", " << fDz << ") cm" << endl;
   cout << "    Flag: " << fFlag << endl;
}

// -------------------------------------------------------------------------
Int_t MpdItsHit5spd::SetDetId(Int_t layer, Int_t ladder, Int_t sector)
{
   // Set detector ID - helper for the tracking

   // Layout with sectors
   fDetectorID = 0;
   fDetectorID |= (sector << kSensorS); // # of sector    on bit 0 - 6
   fDetectorID |= (ladder << kLadderS); // # of ladder    on bit 7 - 12
   fDetectorID |= (layer << kLayerS);   // # of layer     on bit 13 - 15

   return fDetectorID;
}

ClassImp(MpdItsHit5spd)

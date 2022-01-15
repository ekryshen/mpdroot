//--------------------------------------------------------------------------
//----                     MpdBbcHit                                    ----
//--------------------------------------------------------------------------

#include <iostream>
using namespace std;

#include "MpdBbcHit.h"

// -----   Default constructor   -------------------------------------------

MpdBbcHit::MpdBbcHit()
{
   fFlag = 1;
}

// -----   Standard constructor   ------------------------------------------

MpdBbcHit::MpdBbcHit(Int_t detID, TVector3 pos, TVector3 dpos, Int_t index, Int_t flag)
   : FairHit(detID, pos, dpos, index)
{
   fFlag = flag;
}

// -----   Constructor without flag  ------------------------------------------

MpdBbcHit::MpdBbcHit(Int_t detID, TVector3 pos, TVector3 dpos, Int_t index) : FairHit(detID, pos, dpos, index) {}

// -----   Destructor   ----------------------------------------------------
MpdBbcHit::~MpdBbcHit() {}

// -----  Print  -----------------------------------------------------------
void MpdBbcHit::Print(const Option_t *opt) const
{
   cout << "-I- MpdBbcHit" << endl;
   cout << "    DetectorID: " << fDetectorID << endl;
   cout << "    Position: (" << fX << ", " << fY << ", " << fZ << ") cm" << endl;
   cout << "    Position error: (" << fDx << ", " << fDy << ", " << fDz << ") cm" << endl;
   cout << "    Flag: " << fFlag << endl;
}
// -------------------------------------------------------------------------

ClassImp(MpdBbcHit)

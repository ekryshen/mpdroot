#include "Riostream.h"
#include <iomanip>

#include "TpcLheHit.h"

using namespace std;

ClassImp(TpcLheHit);

   //________________________________________________________________
   TpcLheHit::TpcLheHit()
   : TObject()
{
   //---

   SetUsage(kFALSE);
   SetHitNumber(-1);
   SetTrackID(-1);

   SetX(0);
   SetY(0);
   SetZ(0);

   SetXerr(0);
   SetYerr(0);
   SetZerr(0);

   SetFlag(1);
   SetEdep(0);
   SetStep(0);
}

//________________________________________________________________
TpcLheHit::TpcLheHit(TpcPoint *point) : TObject()
{

   SetUsage(kFALSE);
   SetHitNumber(-1);
   SetTrackID(-1);

   Float_t error = 0.0;
   SetXerr(error);
   SetYerr(error);
   SetZerr(error);

   SetFlag(1);
   SetEdep(0);
   SetStep(0);
}

//________________________________________________________________
TpcLheHit::TpcLheHit(Double_t *x, Int_t stn) : TObject()
{
   // Constructor which takes the coodrinates and the station number.

   SetUsage(kFALSE);
   SetHitNumber(-1);

   SetTrackID(-1);

   SetX(x[0]);
   SetY(x[1]);
   SetZ(x[2]);

   SetXerr(0.);
   SetYerr(0.);
   SetZerr(0.);

   SetFlag(1);
   SetEdep(0);
   SetStep(0);
}

//________________________________________________________________
TpcLheHit::~TpcLheHit()
{
   // ---
}

//________________________________________________________________
void TpcLheHit::Print()
{

   cout << flush;
   cout << " " << fTrackID;
   cout << " " << fHitNumber;
   cout << " " << Int_t(fUsed);

   cout.setf(ios::right, ios::floatfield);
   cout.precision(3);
   cout << " " << setw(7) << fCoord.X();
   cout << " " << setw(7) << fCoord.Y();
   cout << " " << setw(5) << fCoord.Z();
   cout.setf(ios::scientific);
   cout.precision(5);
   cout << " " << setw(10) << fError.X();
   cout << " " << setw(10) << fError.Y();
   cout << " " << setw(10) << fError.Z();
   cout << endl;
   cout.setf(ios::right, ios::floatfield);
}

//________________________________________________________________
Int_t TpcLheHit::Compare(const TObject *hit) const
{
   /// "Compare" function to sort in descending order in radius

   TpcLheHit *tpcHit = (TpcLheHit *)hit;
   if (fLayer < tpcHit->GetLayer())
      return 1;
   else if (fLayer > tpcHit->GetLayer())
      return -1;
   else {
      if (fR < tpcHit->GetR())
         return 1;
      else if (fR > tpcHit->GetR())
         return -1;
      else
         return 0;
   }
}

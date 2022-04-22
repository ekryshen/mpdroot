#include "TpcLhePoint.h"

ClassImp(TpcLhePoint);

   //________________________________________________________________________
   TpcLhePoint::TpcLhePoint()
{
   // Default constructor.
   SetX(0.);
   SetY(0.);
   SetZ(0.);

   SetXerr(0.);
   SetYerr(0.);
   SetZerr(0.);
}

//________________________________________________________________________
TpcLhePoint::TpcLhePoint(Double_t pos[3], Double_t err[3])
{
   SetX((Double_t)pos[0]);
   SetY((Double_t)pos[1]);
   SetZ((Double_t)pos[2]);

   SetXerr((Double_t)err[0]);
   SetYerr((Double_t)err[1]);
   SetZerr((Double_t)err[2]);
}

//________________________________________________________________________
TpcLhePoint::TpcLhePoint(Double_t x, Double_t y, Double_t z)
{
   SetX(x);
   SetY(y);
   SetZ(z);
}

//________________________________________________________________________
TpcLhePoint::~TpcLhePoint()
{
   //---
}

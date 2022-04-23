#include <Riostream.h>
#include "TMath.h"

#include "TpcLheCMPoint.h"

using namespace std;

ClassImp(TpcLheCMPoint);

//________________________________________________________________
TpcLheCMPoint::TpcLheCMPoint() : TpcLheHit()
{
   //--- clean all

   SetPhi(0.);
   SetTheta(0.);

   fXprime    = 0.;
   fYprime    = 0.;
   fXprimeerr = 0.;
   fYprimeerr = 0.;

   SetIntPoint(0., 0., 0., 0., 0., 0.);
   SetShiftedCoord();
   SetDist(0., 0.);
}

//________________________________________________________________
TpcLheCMPoint::TpcLheCMPoint(TpcLheHit *point) : TpcLheHit(*point)
{
   //---
}

//________________________________________________________________
TpcLheCMPoint::~TpcLheCMPoint()
{
   //---
}

//________________________________________________________________
void TpcLheCMPoint::SetIntPoint(const Double_t in_x, const Double_t in_y, const Double_t in_z, const Double_t in_x_err,
                                const Double_t in_y_err, const Double_t in_z_err)
{
   // Defines a new interaction point. This point is needed to
   // calculate the conformal coordinates.

   SetXt(in_x);
   SetYt(in_y);
   SetZt(in_z);

   SetXterr(in_x_err);
   SetYterr(in_y_err);
   SetZterr(in_z_err);
}

//________________________________________________________________
void TpcLheCMPoint::SetAllCoord(const TpcLheCMPoint *hit)
{

   // Sets the interaction point, the shifted coordinates,
   // and the conformal mapping coordinates.
   // These values are calculated from the interaction point of the given hit
   // which should be  already found hit on the same track.

   SetIntPoint(hit->GetX(), hit->GetY(), hit->GetZ(), hit->GetXerr(), hit->GetYerr(), hit->GetZerr());
   SetShiftedCoord();
   SetConfCoord();
}

//________________________________________________________________
void TpcLheCMPoint::SetShiftedCoord()
{
   // Sets the coordinates with respect to the given vertex point

   SetXv(GetX() - fXt);
   SetYv(GetY() - fYt);
   SetZv(GetZ() - fZt);

   SetXverr(TMath::Sqrt(GetXerr() * GetXerr() + fXterr * fXterr));
   SetYverr(TMath::Sqrt(GetYerr() * GetYerr() + fYterr * fYterr));
   SetZverr(TMath::Sqrt(GetZerr() * GetZerr() + fZterr * fZterr));
}

//________________________________________________________________
void TpcLheCMPoint::SetConfCoord()
{
   // Calculates the conformal coordinates of one hit.

   Double_t r2 = fXv * fXv + fYv * fYv;

   if (r2 != 0.) {
      fXprime = fXv / r2;
      fYprime = -fYv / r2;

      fXprimeerr =
         TMath::Sqrt(TMath::Power((-fXv * fXv + fYv * fYv) * fXverr, 2) + TMath::Power(2 * fXv * fYv * fYverr, 2)) /
         TMath::Power(fXv * fXv + fYv * fYv, 2);

      fYprimeerr = TMath::Sqrt(TMath::Power((-fXv * fXv - 3 * fYv * fYv) * fYverr, 2) +
                               TMath::Power(-2 * fXv * fYv * fXverr, 2)) /
                   TMath::Power(fXv * fXv + fYv * fYv, 2);
   }

   else {

      fXprime    = 0.;
      fYprime    = 0.;
      fXprimeerr = 0.;
      fYprimeerr = 0.;
   }
}

//________________________________________________________________
void TpcLheCMPoint::SetAngles()
{
   // Calculates the azimutal angle phi and the polar theta for each hit.

   TVector3 v(fXv, fYv, fZv);
   fPhi   = (v.Phi() > 0.) ? fPhi = v.Phi() : fPhi = v.Phi() + 2. * TMath::Pi();
   fTheta = v.Theta();
}

//________________________________________________________________
void TpcLheCMPoint::Setup(TpcLhePoint *vertex)
{
   // ---  right order is important

   SetIntPoint(vertex->GetX(), vertex->GetY(), vertex->GetZ(), vertex->GetXerr(), vertex->GetYerr(), vertex->GetZerr());

   // The angles are set properly if they are set after the interaction
   // point and the shifted coordinates

   SetShiftedCoord();
   SetConfCoord();

   SetAngles();
   SetDist(0., 0.);
}

//________________________________________________________________
void TpcLheCMPoint::Print()
{

   cout << flush;
   cout << " " << this->GetTrackID();
   cout << " " << this->GetHitNumber();
   cout << " " << this->GetUsage();

   cout.setf(ios::right, ios::floatfield);
   cout.precision(3);
   cout << " " << setw(7) << this->GetX();
   cout << " " << setw(7) << this->GetY();
   cout << " " << setw(5) << this->GetZ();
   cout.setf(ios::scientific);
   cout.precision(5);
   cout << " " << setw(10) << fPhi;
   cout << " " << setw(10) << fTheta;
   cout << endl;
   cout.setf(ios::right, ios::floatfield);
}

#include "TpcLheTrackCuts.h"

ClassImp(TpcLheTrackCuts);

TpcLheTrackCuts *TpcLheTrackCuts::fInstance = 0;

//________________________________________________________________
TpcLheTrackCuts *TpcLheTrackCuts::Instance()
{
   //--- Returns instance.

   if (fInstance == 0) {
      fInstance = new TpcLheTrackCuts();
   }
   return fInstance;
}

//________________________________________________________________
TpcLheTrackCuts::TpcLheTrackCuts()
{
   //---
}

//________________________________________________________________
TpcLheTrackCuts::~TpcLheTrackCuts()
{
   //---
   if (fInstance) delete fInstance;
}

//________________________________________________________________
void TpcLheTrackCuts::Reset() {}

//________________________________________________________________
Double_t const TpcLheTrackCuts::DistFromLine(Double_t Xh, Double_t Yh, Double_t *coeff)
{
   // Returns the distance of a point to a straight line.
   // straight line is given by its to coefficients: y = coeff[0]*x + coeff[1].

   Double_t x = (coeff[0] / (1 + coeff[0] * coeff[0])) * (1 / coeff[0] * Xh + Yh - coeff[1]);

   return TMath::Sqrt(TMath::Power(x - Xh, 2) + TMath::Power(coeff[0] * x + coeff[1] - Yh, 2));
}

//________________________________________________________________
void TpcLheTrackCuts::SetHighMomTrackCuts()
{
   //--- Sets cuts of tracks

   //  Range for prediction

   fTrackCuts[0].fDelX = .001; //

   fTrackCuts[1].fDelX = 3 * .0055;
   // 2 -- seed
   fTrackCuts[0].fMaxHitDist = 3.;
   fTrackCuts[1].fMaxHitDist = 5.;
}

//________________________________________________________________
void TpcLheTrackCuts::SetLowMomTrackCuts()
{
   // Sets cuts of tracks for the given vertex constraint.

   //  Range for prediction

   fTrackCuts[0].fDelX = 3 * .025; //
   fTrackCuts[1].fDelX = 3 * .034;
}

//________________________________________________________________
Bool_t TpcLheTrackCuts::VerifyTrack(TpcLheCMTrack *track, TpcLheCMPoint *hit, Bool_t back)
{
   // --------------------------------

   if (track->GetNumberOfPoints() < 3) {
      // tracklet
      return kTRUE;
   } else {
      // track

      Double_t coeff[3];
      coeff[0] = coeff[1] = coeff[2] = 0;

      CMLineFit(track, coeff);

      Double_t dist_cm = DistFromLine(hit->GetXprime(), hit->GetYprime(), coeff);

      if (dist_cm < fTrackCuts[0].fDelX) {
         return kTRUE;
      } else
         return kFALSE;
   }

} ///:~

//________________________________________________________________
Double_t TpcLheTrackCuts::GetChi2Bend(Int_t pl)
{
   //---
   return fTrackCuts[pl].fBendChi2;
}

//________________________________________________________________
Double_t TpcLheTrackCuts::GetChi2Deep(Int_t pl)
{
   //---
   return fTrackCuts[pl].fDeepChi2;
}

//________________________________________________________________
Double_t TpcLheTrackCuts::GetDelX(Int_t pl)
{
   //---
   return fTrackCuts[pl].fDelX;
}

//________________________________________________________________
void TpcLheTrackCuts::Circle3pnts(Double_t x[], Double_t y[], Double_t r[])
{
   // calc center and R of circle from 3 points

   Double_t m_a = (y[1] - y[0]) / (x[1] - x[0]);
   Double_t m_b = (y[2] - y[1]) / (x[2] - x[1]);

   if (m_a != m_b) {
      Double_t x_c = .5 * (m_a * m_b * (y[0] - y[2]) + m_b * (x[0] + x[1]) - m_a * (x[1] + x[2])) / (m_b - m_a);

      Double_t y_c = -(x_c - .5 * (x[0] + x[1])) / m_a + .5 * (y[0] + y[1]);

      Double_t dely = y_c - y[0];
      Double_t delx = x_c - x[0];
      Double_t rad  = TMath::Sqrt(delx * delx + dely * dely);

      r[0] = x_c;
      r[1] = y_c;
      r[2] = rad;
   }
}

//________________________________________________________________
int TpcLheTrackCuts::Circle_Circle_Intersection(double x0, double y0, double r0, double x1, double y1, double r1,
                                                double *xi, double *yi, double *xi_prime, double *yi_prime)
{
   //---

   double a, dx, dy, d, h, rx, ry;
   double x2, y2;

   /* dx and dy are the vertical and horizontal distances between
    * the circle centers.
    */
   dx = x1 - x0;
   dy = y1 - y0;

   /* Determine the straight-line distance between the centers. */
   d = sqrt((dy * dy) + (dx * dx));

   /* Check for solvability. */
   if (d > (r0 + r1)) {
      /* no solution. circles do not intersect. */
      return 0;
   }

   if (d < TMath::Abs(r0 - r1)) {
      /* no solution. one circle is contained in the other */
      return 0;
   }

   /* 'point 2' is the point where the line through the circle
    * intersection points crosses the line between the circle
    * centers.
    */

   /* Determine the distance from point 0 to point 2. */
   a = ((r0 * r0) - (r1 * r1) + (d * d)) / (2.0 * d);

   /* Determine the coordinates of point 2. */
   x2 = x0 + (dx * a / d);
   y2 = y0 + (dy * a / d);

   /* Determine the distance from point 2 to either of the intersection
    * points.
    */
   h = sqrt((r0 * r0) - (a * a));

   /* Now determine the offsets of the intersection points from point
    * 2.
    */
   rx = -dy * (h / d);
   ry = dx * (h / d);

   /* Determine the absolute intersection points. */
   *xi       = x2 + rx;
   *xi_prime = x2 - rx;
   *yi       = y2 + ry;
   *yi_prime = y2 - ry;

   return 1;
}

//________________________________________________________________
Double_t TpcLheTrackCuts::GetPhiPrediction(TpcLheCMTrack *track)
{
   //---

   Int_t      last = (track->GetRHits())->GetLast();
   TpcLheHit *hit0 = (TpcLheHit *)(track->GetRHits())->At(last);
   TpcLheHit *hit1 = (TpcLheHit *)(track->GetRHits())->At(last - 1);

   Double_t x[3], y[3], rez[3];

   x[0] = 0.; // hit0->GetX();
   x[1] = hit0->GetX();
   x[2] = hit1->GetX();

   y[0] = 0.; // hit0->GetY();
   y[1] = hit0->GetY();
   y[2] = hit1->GetY();

   Circle3pnts(x, y, rez);

   Double_t x_c = rez[0];
   Double_t y_c = rez[1];
   Double_t Rad = rez[2];
   track->SetRadius(Rad);
   track->SetVertex(x_c, y_c, 0.);

   TVector3 v(x_c, y_c, 0.);
   Double_t phi_c = (v.Phi() >= 0.) ? phi_c = v.Phi() : phi_c = v.Phi() + 2. * TMath::Pi();

   return phi_c;
}

//________________________________________________________________
Double_t TpcLheTrackCuts::GetThetPrediction(TpcLheCMTrack *track, Double_t zst, Bool_t back)
{
   //--- return alpha angle in the next plane

   TpcLheHit *hit1, *hit0;

   if (back) {
      hit1 = (TpcLheHit *)(track->GetRHits())->At(1);
      hit0 = (TpcLheHit *)(track->GetRHits())->At(0);
   } else {
      Int_t last = (track->GetRHits())->GetLast();
      hit0       = (TpcLheHit *)(track->GetRHits())->At(last);
      hit1       = (TpcLheHit *)(track->GetRHits())->At(last - 1);
   }

   Double_t a0 = (hit0->GetY() - hit1->GetY()) / (hit0->GetZ() - hit1->GetZ());
   Double_t a1 = (hit0->GetY() * hit1->GetZ() - hit1->GetY() * hit0->GetZ()) / (hit1->GetZ() - hit0->GetZ());

   Double_t y_cross = a0 * zst + a1;

   //  return  TMath::ATan2(y_cross, zst);
   return y_cross;
}

//________________________________________________________________
void TpcLheTrackCuts::CMLineFit(TpcLheCMTrack *track, Double_t *a)
{
   //---

   TObjArray *trackpoints = track->GetCMHits();
   Int_t      n           = trackpoints->GetEntriesFast();

   TpcLheCMPoint *trackpoint = NULL;

   TArrayD *x = new TArrayD();
   x->Set(n);
   TArrayD *y = new TArrayD();
   y->Set(n);
   TArrayD *delx = new TArrayD();
   delx->Set(n);
   TArrayD *dely = new TArrayD();
   dely->Set(n);

   Int_t ip = 0;
   for (Int_t is = 0; is < n; is++) {
      trackpoint = (TpcLheCMPoint *)trackpoints->At(is);
      //      trackpoint->Print();
      if (trackpoint->GetXprime() != 0.) {

         x->AddAt(trackpoint->GetXprime(), ip);
         y->AddAt(trackpoint->GetYprime(), ip);

         delx->AddAt(trackpoint->GetXprimeerr(), ip);
         dely->AddAt(trackpoint->GetYprimeerr(), ip);
         ip++;
      }
   }

   Double_t coeff[3];

   LineFit(x, delx, y, dely, coeff, ip);

   a[0] = coeff[0];
   a[1] = coeff[1];
   a[2] = coeff[2];

   if (coeff[1] != 0.) {
      Double_t D = TMath::Sqrt((coeff[0] * coeff[0] + 1.) / (4. * coeff[1] * coeff[1]));
      track->SetRadius(D);
   } else {
      track->SetRadius(0.);
   }

   delete x;
   delete y;
   delete delx;
   delete dely;
}

//________________________________________________________________
void TpcLheTrackCuts::LineFit(TArrayD *x, TArrayD *delx, TArrayD *y, TArrayD *dely, Double_t *coeff, Int_t np)
{

   // fit points to a  straight line y = ax + b
   //        NRC p. 661

   Double_t L11 = 0.;
   Double_t L12 = 0.;
   Double_t L22 = 0.;
   Double_t g1  = 0.;
   Double_t g2  = 0.;

   for (Int_t i = 0; i < np; i++) {

      L11 += 1.;                  // S
      L12 += x->At(i);            // Sx
      L22 += x->At(i) * x->At(i); // Sxx
      g1 += y->At(i);             // Sy
      g2 += x->At(i) * y->At(i);  // Sxy
   }

   Double_t D = L11 * L22 - L12 * L12; // S*Sxx - (Sx)^2

   coeff[0] = (g2 * L11 - g1 * L12) / D; //   a = (Sxy*S - Sy*Sx)/D
   coeff[1] = (g1 * L22 - g2 * L12) / D; //   b = (Sy*Sxx - Sxy*Sx)/D

   // Calculate chi squared

   Double_t chi2 = 0.;
   if (np > 2) {
      for (Int_t i = 0; i < np; i++) {

         chi2 +=
            TMath::Power(y->At(i) - coeff[0] * x->At(i) - coeff[1], 2.) / TMath::Abs(coeff[0] * x->At(i) + coeff[1]);
      }
      chi2 = chi2 / float(np - 2);
   }
   coeff[2] = chi2;
}

//________________________________________________________________
void TpcLheTrackCuts::DeepAngleFit(const TpcLheCMTrack *track, Double_t *a)
{
   //---

   TObjArray *trackpoints = track->GetCMHits();
   Int_t      n           = trackpoints->GetEntriesFast();

   if (n > 2) {
      TArrayD *zv = new TArrayD();
      zv->Set(n);
      TArrayD *yv = new TArrayD();
      yv->Set(n);
      TArrayD *delzv = new TArrayD();
      delzv->Set(n);
      TArrayD *delyv = new TArrayD();
      delyv->Set(n);

      for (Int_t is = 0; is < n; is++) {
         TpcLheCMPoint *trackpoint = (TpcLheCMPoint *)trackpoints->At(is);
         //      trackpoint->Print();

         zv->AddAt(trackpoint->GetZv(), is);
         yv->AddAt(trackpoint->GetYv(), is);

         delzv->AddAt(trackpoint->GetZverr(), is);
         delyv->AddAt(trackpoint->GetYverr(), is);
      }

      Double_t coeff[3];

      LineFit(zv, delzv, yv, delyv, coeff, n);

      a[0] = coeff[0];
      a[1] = coeff[1];
      a[2] = coeff[2];

      delete zv;
      delete yv;
      delete delzv;
      delete delyv;
   } else {
      TpcLheCMPoint *hit1 = (TpcLheCMPoint *)trackpoints->At(0);
      TpcLheCMPoint *hit2 = (TpcLheCMPoint *)trackpoints->At(1);
      a[0]                = (hit1->GetY() - hit2->GetY()) / (hit1->GetZ() - hit2->GetZ());
      a[1]                = (hit1->GetY() * hit2->GetZ() - hit2->GetY() * hit1->GetZ()) / (hit2->GetZ() - hit1->GetZ());
   }
}

//__________________________________________________________________
Bool_t TpcLheTrackCuts::IsGoodFoundTrack(TpcLheCMTrack *track)
{
   //---

   if (track->GetNumberOfPoints() < 3)
      return kFALSE;
   else
      return kTRUE;
}

//__________________________________________________________________
Bool_t TpcLheTrackCuts::IsGoodGeantTrack(TpcLheTrack *track)
{
   //--- Returns true if the given track fulfills all requirements to
   // be a "good" track.

   return kTRUE;
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

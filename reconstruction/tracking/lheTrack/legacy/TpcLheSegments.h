#ifndef LHE_SEGMENTS_H
#define LHE_SEGMENTS_H

// *************************************************************************
//  Author: Oleg Rogachevsky             e-mail: rogach@sunhe.jinr.ru
//
//  arrange all hits in theta-phi cells
//
// Created: 1-07-07
// Modified:
//
// *************************************************************************

#include "Riostream.h"
#include "TClonesArray.h"
#include "assert.h"

#include "TpcLheCMPoint.h"

using namespace std;

class TpcLheSegments {

protected:
   TObjArray *fSegments; // array of volume (theta, phi) elements
   TVector3  *fVertex;

   Int_t fNumThetaSegment; // number of theta segments
   Int_t fNumPhiSegment;   // number of phi segments
   Int_t fBounds;          // Number of cells (segments)

   // max. min. phi values
   Double_t fPhiMin;
   Double_t fPhiMax;

   // max. min. theta values
   Double_t fThetaMin;
   Double_t fThetaMax;

public:
   TpcLheSegments(); //
   TpcLheSegments(Int_t nTheta, Int_t nPhi);
   virtual ~TpcLheSegments(); //

   void FillSegments(TClonesArray *hits); //
   void PrintSegmentContents(Int_t n);
   void PrintSegments();
   void Init();
   void Clear();

   // getters

   Int_t GetNumPhiSegments() { return fNumPhiSegment; }     // returns # of phi
   Int_t GetNumThetaSegments() { return fNumThetaSegment; } // returns # of theta
   // Int_t  GetRadiusSegm(const TpcLheCMPoint *hit);  // returns # of station
   // Int_t  GetRadiusSegm(Int_t segm);        //  number of station segment
   Int_t GetThetaSegm(const TpcLheCMPoint *hit); // returns # of theta
   Int_t GetThetaSegm(Int_t segm);               //  number of theta segment
   Int_t GetThetaSegm(Double_t s);               //  number of theta segment
   Int_t GetPhiSegm(const TpcLheCMPoint *hit);   // returns # of phi
   Int_t GetPhiSegm(Int_t segm);                 //  number of phi
   Int_t GetPhiSegm(Double_t al);                //  number of phi segment
   // Int_t  GetStation(Int_t segm);          // returns the station

   Double_t GetTheta(Int_t segm);               // returns theta of a given theta segment
   Double_t GetPhi(Int_t segm);                 // returns phi of a given phi segment
   Double_t GetThetaMax() { return fThetaMax; } //
   Double_t GetPhiMax() { return fPhiMax; }     //

   TObjArray *GetSegments() { return fSegments; }  //
   Int_t      GetSegm(Int_t theta_s, Int_t phi_s); // returns number of segment

   ClassDef(TpcLheSegments, 1); //
};

//_________________________________________________________________
inline Int_t TpcLheSegments::GetThetaSegm(const TpcLheCMPoint *hit)
{
   // Returns number of theta segment of a specific hit.

   Double_t theta = hit->GetTheta();
   Int_t    theta_segm;
   if (theta < fThetaMin) {
      cout << " hit is smaller than fThetaMin\n";
      theta = fThetaMin;
   }

   if (theta > fThetaMax) {
      cout << " hit is larger than fThetaMax\n";
      theta = fThetaMax;
   }

   theta_segm = (Int_t)((theta - fThetaMin) / (fThetaMax - fThetaMin) * fNumThetaSegment);

   return theta_segm;
}

//_________________________________________________________________
inline Int_t TpcLheSegments::GetThetaSegm(Double_t theta)
{
   // Returns number of theta segment of a hit theta.

   if (theta > fThetaMax) return fNumThetaSegment - 1; // theta = fThetaMax;
   if (theta < fThetaMin) theta = fThetaMin;

   return (Int_t)((theta - fThetaMin) / (fThetaMax - fThetaMin) * fNumThetaSegment);
}

//_________________________________________________________________
inline Double_t TpcLheSegments::GetTheta(Int_t segm)
{
   // Returns the angle theta of a given segment.

   return segm * (fThetaMax - fThetaMin) / ((float)fNumThetaSegment) + fThetaMin;
}

//_________________________________________________________________
inline Int_t TpcLheSegments::GetPhiSegm(const TpcLheCMPoint *hit)
{
   // Returns number of phi segment of a hit.

   Double_t phi      = hit->GetPhi();
   Int_t    phi_segm = (Int_t)((phi - fPhiMin) * fNumPhiSegment / (fPhiMax - fPhiMin));
   return phi_segm;
}

//_________________________________________________________________
inline Int_t TpcLheSegments::GetPhiSegm(Double_t phi)
{
   // Returns number of phi segment of a specific phi.

   if (phi > fPhiMax) phi = fPhiMax;
   if (phi < fPhiMin) phi = fPhiMin;

   return (Int_t)((phi - fPhiMin) / (fPhiMax - fPhiMin) * fNumPhiSegment);
}

//_________________________________________________________________
inline Int_t TpcLheSegments::GetPhiSegm(Int_t segm)
{
   // Returns number of phi segment of a specifiv segment.

   return (segm % (fNumThetaSegment * fNumPhiSegment)) % (fNumPhiSegment);
}

//_________________________________________________________________
inline Double_t TpcLheSegments::GetPhi(Int_t segm)
{
   // Returns the  phi of the given segment.

   return segm * (fPhiMax - fPhiMin) / (fNumPhiSegment) + fPhiMin;
}

//_________________________________________________________________
inline Int_t TpcLheSegments::GetSegm(Int_t theta_segm, Int_t phi_segm)
{
   // Calculates the volume segment number from the segmented volumes
   //   (segm = segm(theta,phi)).

   Int_t segm = theta_segm * (fNumPhiSegment) + phi_segm;

   if (segm >= fBounds) {
      cout << " Error in segment: "
           << " segm " << segm << " theta_segm " << theta_segm << " phi_segm " << phi_segm << endl;
      assert(kFALSE);
      return fBounds - 1;
   } else
      return segm;
}

#if 0
//_________________________________________________________________
inline Int_t TpcLheSegments::GetRadiusSegm(Int_t segm) {
  // Returns number of station segment of a specifiv segment.

  return (segm - GetPhiSegm(segm) - GetThetaSegm(segm)) /
    (fNumThetaSegment * fNumPhiSegment);
}
#endif

//_________________________________________________________________
inline Int_t TpcLheSegments::GetThetaSegm(Int_t segm)
{
   // Returns number of theta segment of a specifiv segment.

   return (segm - GetPhiSegm(segm)) % (fNumThetaSegment * fNumPhiSegment) / (fNumPhiSegment);
}

#endif

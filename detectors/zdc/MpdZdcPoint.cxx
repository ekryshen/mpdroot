/////////////////////////////////////////////////////////////
//
//  MpdZdcPoint
//
//  Geant point for Hyp detector
//
//  Created 14/08/06  by S.Spataro
//
///////////////////////////////////////////////////////////////

#include <iostream>
#include "MpdZdcPoint.h"

// -----   Default constructor   -------------------------------------------
MpdZdcPoint::MpdZdcPoint() : FairMCPoint()
{
   //   fTrackID    = -1;
   //   fDetectorID = -1;
   //   //  fEventID    = -1;
   //   fX          = fY  = fZ =  0.;
   //   fPx         = fPy = fPz = 0.;
   //   fTime       =  0.;
   //   fLength     =  0.;
   //   fELoss      =  0.;
   //   fEventId

   nCopy       = -1;
   nCopyMother = -1; // module
   nCopyZdc    = -1; // zdc (left,right)

   // nCopy_h = -1;
   // nCopyMother_h=-1; //module with hole
}
// -------------------------------------------------------------------------

// -----   Standard constructor   ------------------------------------------
// MpdZdcPoint::MpdZdcPoint(Int_t trackID, Int_t detID, Int_t copyNo, Int_t copyNoMother,  Int_t copyNo_h, Int_t
// copyNoMother_h, TVector3 pos,TVector3 mom, Double_t tof, Double_t length,Double_t eLoss, UInt_t EventId)
// MpdZdcPoint::MpdZdcPoint(Int_t trackID, Int_t detID, Int_t copyNo, Int_t copyNoMother, TVector3 pos,TVector3 mom,
// Double_t tof, Double_t length,Double_t eLoss, UInt_t EventId)
MpdZdcPoint::MpdZdcPoint(Int_t trackID, Int_t detID, Int_t copyNo, Int_t copyNoMother, Int_t copyNoZdc, TVector3 pos,
                         TVector3 mom, Double_t tof, Double_t length, Double_t eLoss, UInt_t EventId)
   : FairMCPoint(trackID, detID, pos, mom, tof, length, eLoss, EventId)
{
   nCopy       = copyNo;
   nCopyMother = copyNoMother;
   nCopyZdc    = copyNoZdc;

   // nCopy_h = copyNo_h;
   // nCopyMother_h =  copyNoMother_h;
   // nCopy_h = copyNo;
   // nCopyMother_h =  copyNoMother;
   //  fTrackID    = trackID;
   //  fDetectorID = detID;

   //  fX          = pos.X();
   //  fY          = pos.Y();
   //  fZ          = pos.Z();
   //  fPx         = mom.Px();
   //  fPy         = mom.Py();
   //  fPz         = mom.Pz();
   //  fTime       = tof;
   //  fLength     = length;
   //  fELoss      = eLoss;
   //  fEventId
}

// -----   Destructor   ----------------------------------------------------
MpdZdcPoint::~MpdZdcPoint() {}
// -------------------------------------------------------------------------

// -----   Public method Print   -------------------------------------------
void MpdZdcPoint::Print(const Option_t *opt) const
{
   cout << "-I- MpdZdcPoint: MUO Point for track " << fTrackID << " in detector " << fDetectorID << " " << nCopyMother
        << " " << nCopy << " " << nCopyZdc << endl;
   cout << "    Position (" << fX << ", " << fY << ", " << fZ << ") cm" << endl;
   cout << "    Momentum (" << fPx << ", " << fPy << ", " << fPz << ") GeV" << endl;
   cout << "    Time " << fTime << " ns,  Length "
        << fLength
        //<< " cm,  Energy loss " << fELoss*1.0e06 << " keV" << endl;
        << " cm,  Energy loss " << fELoss << " GeV" << endl;
}
// -------------------------------------------------------------------------

void MpdZdcPoint::IncreaseTimesVR(Float_t e, Float_t time)
{
  // Update hit information
  if (fContribTimes.find(time) == fContribTimes.end()) fContribTimes[time] = e;
  else fContribTimes[time] += e;
}

// -------------------------------------------------------------------------

ClassImp(MpdZdcPoint);

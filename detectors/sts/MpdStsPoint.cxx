//------------------------------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------
// -----                      MpdStsPoint source file                  -----
// -------------------------------------------------------------------------

#include <iostream>
#include "MpdStsPoint.h"

//------------------------------------------------------------------------------------------------------------------------
MpdStsPoint::MpdStsPoint() : FairMCPoint() {}
//------------------------------------------------------------------------------------------------------------------------
MpdStsPoint::MpdStsPoint(Int_t trackID, Int_t detID, TVector3 posIn, TVector3 momIn, TVector3 posOut, Double_t tof,
                         Double_t length, Double_t eLoss)
   : FairMCPoint(trackID, detID, posIn, momIn, tof, length, eLoss)
{
   SetPositionOut(posOut);
}
//------------------------------------------------------------------------------------------------------------------------
MpdStsPoint::~MpdStsPoint() {}
//------------------------------------------------------------------------------------------------------------------------
void MpdStsPoint::Print(const Option_t *opt) const
{
   cout << "-I- MpdStsPoint: STS point for track " << fTrackID << " in detector " << fDetectorID << endl;
   cout << "    Position (" << fX << ", " << fY << ", " << fZ << ") cm" << endl;
   cout << "    Momentum (" << fPx << ", " << fPy << ", " << fPz << ") GeV" << endl;
   cout << "    Time " << fTime << " ns,  Length " << fLength << " cm,  Energy loss " << fELoss * 1.0e06 << " keV"
        << endl;
}
//------------------------------------------------------------------------------------------------------------------------
ClassImp(MpdStsPoint);

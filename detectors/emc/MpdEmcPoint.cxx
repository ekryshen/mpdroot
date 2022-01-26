#include "MpdEmcPoint.h"

#include <iostream>
using std::cout;
using std::endl;


// -----   Default constructor   -------------------------------------------
MpdEmcPoint::MpdEmcPoint() : FairMCPoint() { }
// -------------------------------------------------------------------------

// -----   Standard constructor   ------------------------------------------
MpdEmcPoint::MpdEmcPoint(Int_t trackID, Int_t detID, 
					 TVector3 pos, TVector3 mom, 
					 Double_t tof, Double_t length,
					 Double_t ELoss)
  : FairMCPoint(trackID, detID, pos, mom, tof, length, ELoss) { }
// -------------------------------------------------------------------------

// -----   Destructor   ----------------------------------------------------
MpdEmcPoint::~MpdEmcPoint() { }
// -------------------------------------------------------------------------

// -----   Public method Print   -------------------------------------------
void MpdEmcPoint::Print(const Option_t* opt) const {
  cout << "-I- MpdEmcPoint: TutorialDet point for track " << fTrackID
       << " in detector " << fDetectorID << endl;
  cout << "    Position (" << fX << ", " << fY << ", " << fZ
       << ") cm" << endl;
  cout << "    Momentum (" << fPx << ", " << fPy << ", " << fPz
       << ") GeV" << endl;
  cout << "    Time " << fTime << " ns,  Length " << fLength
       << " cm,  Energy loss " << fELoss*1.0e06 << " keV" << endl;
}
// -------------------------------------------------------------------------

ClassImp(MpdEmcPoint)
  

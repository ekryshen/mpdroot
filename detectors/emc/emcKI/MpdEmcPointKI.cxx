#include "MpdEmcPointKI.h"

#include <iostream>
using std::cout;
using std::endl;

// -----   Default constructor   -------------------------------------------
MpdEmcPointKI::MpdEmcPointKI() : FairMCPoint() {}
// -------------------------------------------------------------------------

// -----   Standard constructor   ------------------------------------------
MpdEmcPointKI::MpdEmcPointKI(Int_t trackID, Int_t detID, TVector3 pos, TVector3 mom, Double_t tof, Double_t length,
                         Double_t ELoss)
  : FairMCPoint(trackID, detID, pos, mom, tof, length, ELoss)
{
}
// -------------------------------------------------------------------------

// -----   Destructor   ----------------------------------------------------
MpdEmcPointKI::~MpdEmcPointKI() {}
// -------------------------------------------------------------------------

// -----   Public method Print   -------------------------------------------
void MpdEmcPointKI::Print(const Option_t* opt) const
{
  cout << "-I- MpdEmcPointKI: TutorialDet point for track " << fTrackID << " in detector " << fDetectorID << endl;
  cout << "    Position (" << fX << ", " << fY << ", " << fZ << ") cm" << endl;
  cout << "    Momentum (" << fPx << ", " << fPy << ", " << fPz << ") GeV" << endl;
  cout << "    Time " << fTime << " ns,  Length " << fLength << " cm,  Energy loss " << fELoss * 1.0e06 << " keV"
       << endl;
}
// -------------------------------------------------------------------------

Bool_t MpdEmcPointKI::operator<(const MpdEmcPointKI& rhs) const
{
  if (GetDetectorID() == rhs.GetDetectorID())
    return GetTrackID() < rhs.GetTrackID();
  return GetDetectorID() < rhs.GetDetectorID();
}

Bool_t MpdEmcPointKI::operator==(const MpdEmcPointKI& rhs) const
{
  return ((GetDetectorID() == rhs.GetDetectorID()) && (GetTrackID() == rhs.GetTrackID()));
}

MpdEmcPointKI& MpdEmcPointKI::operator+=(const MpdEmcPointKI& rhs)
{
  if (rhs.GetEnergyLoss() > GetEnergyLoss())
    SetTime(rhs.GetTime());
  SetEnergyLoss(GetEnergyLoss() + rhs.GetEnergyLoss());
  return *this;
}

MpdEmcPointKI MpdEmcPointKI::operator+(const MpdEmcPointKI& rhs) const
{
  MpdEmcPointKI result(*this);
  if (rhs.GetEnergyLoss() > result.GetEnergyLoss())
    result.SetTime(rhs.GetTime());
  result.SetEnergyLoss(result.GetEnergyLoss() + rhs.GetEnergyLoss());
  return *this;
}
Int_t MpdEmcPointKI::Compare(const TObject* obj) const
{
  // Compares two MpdEmcPointKI with respect to its Id
  // to sort according increasing Id

  const MpdEmcPointKI* rhs = dynamic_cast<const MpdEmcPointKI*>(obj);

  if (!rhs) {
    return 1;
  }
  if (GetDetectorID() == rhs->GetDetectorID()) {
    if (GetTrackID() < rhs->GetTrackID()) {
      return -1;
    }
    if (GetTrackID() == rhs->GetTrackID()) {
      return 0;
    } else {
      return 1;
    }
  }
  if (GetDetectorID() < rhs->GetDetectorID()) {
    return -1;
  } else {
    return 1;
  }

  return 0; // never happens
}

ClassImp(MpdEmcPointKI)

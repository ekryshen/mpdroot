////////////////////////////////////////////////////////////////
//                                                            //
//  Authors : D.Peresunko, KI                                 //
//                                                            //
////////////////////////////////////////////////////////////////

#include "MpdEmcDigitKI.h"
#include "MpdEmcGeoUtils.h"
#include "MpdEmcPointKI.h"

#include "TMath.h"

#include <iostream>

using namespace std;

// -----   Default constructor   -------------------------------------------

MpdEmcDigitKI::MpdEmcDigitKI() : fDetId(0), fE(0), fTime(0), fNprimary(0), fPrimary(nullptr), fPrimEdep(nullptr) {}

MpdEmcDigitKI::MpdEmcDigitKI(MpdEmcPointKI* point)
  : fDetId(point->GetDetectorID()),
    fE(point->GetEnergyLoss()),
    fTime(point->GetTime()),
    fNprimary(0),
    fPrimary(nullptr),
    fPrimEdep(nullptr)

{
  if (fE > 0) { // only primaries with non-zero energy deposition count
    fNprimary = 1;
    fPrimary = new Int_t[1];
    fPrimary[0] = point->GetTrackID();
    fPrimEdep = new Float_t[1];
    fPrimEdep[0] = point->GetEnergyLoss();
  }
}

MpdEmcDigitKI::MpdEmcDigitKI(int cellId, float energy, float time, int trackId)
  : fDetId(cellId), fE(energy), fTime(time), fNprimary(0), fPrimary(nullptr), fPrimEdep(nullptr)
{
  if (trackId >= 0 && energy > 0) {
    fNprimary = 1;
    fPrimary = new Int_t[1];
    fPrimary[0] = trackId;
    fPrimEdep = new Float_t[1];
    fPrimEdep[0] = energy;
  }
}
MpdEmcDigitKI::~MpdEmcDigitKI()
{
  if (fPrimary) {
    delete[] fPrimary;
    fPrimary = nullptr;
  }
  if (fPrimEdep) {
    delete[] fPrimEdep;
    fPrimEdep = nullptr;
  }
}
Bool_t MpdEmcDigitKI::CanAdd(MpdEmcPointKI* point) const
{
  // Check, if point can be added (hits from the same Tower)

  return fDetId == (UInt_t) point->GetDetectorID();
}

void MpdEmcDigitKI::AddPoint(MpdEmcPointKI* point)
{
  // Adds point (add energy, change time if necessary, add primary)
  if (fDetId != (UInt_t) point->GetDetectorID()) {
    return;
  }
  if (point->GetEnergyLoss() == 0.) { // do nothing
    return;
  }
  if (point->GetEnergyLoss() > fE) {
    fTime = point->GetTime();
  }
  // Check if track already exist
  Int_t iprim = point->GetTrackID();
  Bool_t found = false;
  for (Int_t i = 0; i < fNprimary; i++) {
    if (fPrimary[i] == iprim) {
      fPrimEdep[i] += point->GetEnergyLoss();
      found = true;
      break;
    }
  }
  if (!found) {
    Int_t* tmp = fPrimary;
    fPrimary = new Int_t[fNprimary + 1];
    for (int i = 0; i < fNprimary; i++) {
      fPrimary[i] = tmp[i];
    }
    fPrimary[fNprimary] = point->GetTrackID();
    if (tmp)
      delete[] tmp;
    Float_t* tmpE = fPrimEdep;
    fPrimEdep = new Float_t[fNprimary + 1];
    for (int i = 0; i < fNprimary; i++) {
      fPrimEdep[i] = tmpE[i];
    }
    fPrimEdep[fNprimary] = point->GetEnergyLoss();

    if (tmpE) {
      delete[] tmpE;
    }
    fNprimary++;
  }
  fE += point->GetEnergyLoss();
}

// -----  Print  -----------------------------------------------------------

void MpdEmcDigitKI::Print(const Option_t* opt) const
{
  cout << "MpdEmcDigitKI: " << endl;
  cout << "\tDetId: " << fDetId << "\tDeposited energy: " << fE << "\ttime: " << fTime << endl;
  cout << "\tNumber of tracks in tower: " << fNprimary << endl;
  cout << "Tracks: " << endl;
  for (int i = 0; i < fNprimary; i++) {
    cout << i << " TrackId " << fPrimary[i] << ", deposited E " << fPrimEdep[i] << endl;
  }
}

Int_t MpdEmcDigitKI::Compare(const TObject* obj) const
{
  const MpdEmcDigitKI* rhs = dynamic_cast<const MpdEmcDigitKI*>(obj);
  if (!rhs) {
    return 1;
  }
  if (fDetId < rhs->fDetId) {
    return -1;
  } else {
    if (fDetId == rhs->fDetId) {
      return 0;
    } else {
      return 1;
    }
  }
}

double MpdEmcDigitKI::GetZcenter()
{
  // return Z coordinate in global system
  MpdEmcGeoUtils* geom = MpdEmcGeoUtils::GetInstance();
  double x, y, z;
  geom->DetIdToGlobalPosition(fDetId, x, y, z);
  return z;
}
double MpdEmcDigitKI::GetPhiCenter()
{
  // return phi angle in global system
  MpdEmcGeoUtils* geom = MpdEmcGeoUtils::GetInstance();
  double x, y, z;
  geom->DetIdToGlobalPosition(fDetId, x, y, z);
  return TMath::ATan2(y, x);
}

// -------------------------------------------------------------------------

ClassImp(MpdEmcDigitKI)
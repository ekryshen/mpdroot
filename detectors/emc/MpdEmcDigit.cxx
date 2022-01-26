//-----------------------------------------------------------
//
// Description:
//      Implementation of class MpdEmcDigit
//      see MpdEmcDigit.h for details
//
// Author List:
//      Alexander Zinchenko LHEP, JINR, Dubna - 8-May-2016
//
//-----------------------------------------------------------

#include "MpdEmcDigit.h"

#include <iostream>

using namespace std;


// -----   Default constructor   -------------------------------------------

MpdEmcDigit::MpdEmcDigit() : FairHit(),
  fChanZId(-1),
  fChanPhiId(-1),
  fE(-1.0),
  fTrackID(-1),
  fFlag(0),
  fPDG(0),
  fNumTracks(0),
  fZCenter(0.0),
  fPhiCenter(0.0) {
}


// -----   Standard constructor   ------------------------------------------

MpdEmcDigit::MpdEmcDigit(Int_t detID, TVector3 pos, TVector3 dpos, Int_t index, Int_t flag)
: FairHit(detID, pos, dpos, index) 
{
  fFlag = flag;
}

// -----   Constructor without flag  ------------------------------------------

MpdEmcDigit::MpdEmcDigit(Int_t detID, TVector3 pos, TVector3 dpos, Int_t index)
: FairHit(detID, pos, dpos, index) 
{
}

// -----   Constructor without flag  ------------------------------------------

MpdEmcDigit::MpdEmcDigit(UInt_t side, UInt_t sector, UInt_t row, UInt_t channel) : FairHit(),
  fE(0),
  fTrackID(-1),
  fFlag(0),
  fPDG(0),
  fNumTracks(0) 
{
  SetDetId(side, sector, row, channel);
  //cout << " Decode: " << Side() << " " << Sector() << " " << Row() << " " << Channel() << endl;
}

// -----   Destructor   ----------------------------------------------------

MpdEmcDigit::~MpdEmcDigit() 
{
}

// -----  Print  -----------------------------------------------------------

void MpdEmcDigit::Print(const Option_t* opt) const {
    cout << "MpdEmcDigit: " << endl;
    //cout << "\tSec: " << fSecId << "   Row: " << fRowId << "   Super module: " << fSupModId << "   Module:" << fModId << endl;
    cout << "\tDeposited energy: " << fE << "   ZCenter: " << fZCenter << "   PhiCenter: " << fPhiCenter << endl;
    cout << "\tNumber of tracks in module: " << fNumTracks << endl;
    if (fNumTracks == 1) cout << "PDG code: " << fPDG << "   Track ID: " << fTrackID << endl;

}
// -------------------------------------------------------------------------

void MpdEmcDigit::IncreaseEnergy(Float_t e, Int_t trId)
{
  // Update hit information

  fE += e;
  if (fContrib.find(trId) == fContrib.end()) fContrib[trId] = e;
  else fContrib[trId] += e;
  fDy = TMath::Max(Double_t(e),fDy);
}
// -------------------------------------------------------------------------

UInt_t MpdEmcDigit::SetDetId(UInt_t side, UInt_t sector, UInt_t row, UInt_t channel)
{
  // Code DetectorID

  fDetectorID = side;                       // side # (Z<>0)  on bit 0
  fDetectorID |= (channel << kChannelS);    // channel #      on bit 1 - 6
  fDetectorID |= (row << kRowS);            // row #          on bit 7 - 12
  fDetectorID |= (sector << kSectorS);      // sector #       on bit 13 - 16
  return fDetectorID;
}

// -------------------------------------------------------------------------

ClassImp(MpdEmcDigit)

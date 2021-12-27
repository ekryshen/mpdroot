//--------------------------------------------------------------------------
//----                     MpdEmcHit                                    ----
//--------------------------------------------------------------------------

////////////////////////////////////////////////////////////////
//                                                            //
//  Updated : Martemianov M., ITEP, 2021                      //
//                                                            //
////////////////////////////////////////////////////////////////

#include "MpdEmcHit.h"

#include <iostream>

using namespace std;


// -----   Default constructor   -------------------------------------------

MpdEmcHit::MpdEmcHit() :
fE(-1.0),
fTime(-1.0),
fTrackID(-1),
fDetectorID(-1),
fFlag(0),
fPDG(0),
fNumTracks(0),
fPhiCenter(0.0), 
fThetaCenter(0.0)

{}


// -----   Standard constructor   ------------------------------------------

MpdEmcHit::MpdEmcHit(Int_t detID, TVector3 pos, TVector3 dpos, Int_t index, Int_t flag)
: FairHit(detID, pos, dpos, index) {
    fFlag = flag;
}

// -----   Constructor without flag  ------------------------------------------

MpdEmcHit::MpdEmcHit(Int_t detID, TVector3 pos, TVector3 dpos, Int_t index)
: FairHit(detID, pos, dpos, index) {
}

MpdEmcHit::MpdEmcHit(UInt_t detID, Float_t e, Float_t time) :
fDetectorID(detID), 
fE(e),
fTime(time),
fTrackID(-1),
fFlag(kFALSE),
fPDG(0),
fNumTracks(0) {
}

// -----   Destructor   ----------------------------------------------------

MpdEmcHit::~MpdEmcHit() {
}

// -----  Print  -----------------------------------------------------------

void MpdEmcHit::Print(const Option_t* opt) const {
    cout << "MpdEmcHit: " << endl;
    cout << "\tdetID: " << fDetectorID << endl;
    cout << "\tDeposited energy: " << fE << "\tMean time: " << fTime << 
    "   RhoCenter: " << sqrt(fX*fX+fY*fY) << "   ZCenter: " << fZ << "   PhiCenter: " << fPhiCenter << 
    "   ThetaCenter: " << fThetaCenter << endl;
    cout << "\tNumber of tracks in module: " << fNumTracks << endl;
    if (fNumTracks == 1) cout << "PDG code: " << fPDG << "   Track ID: " << fTrackID << endl;

}
// -------------------------------------------------------------------------


ClassImp(MpdEmcHit)

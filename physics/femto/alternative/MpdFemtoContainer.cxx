#include <MpdFemtoContainer.h>

//--------------------------------------------------------------------------

MpdFemtoContainer::MpdFemtoContainer() :
fEventNumber(0),
fPhi(0.),
fTheta(0.)
{
    fMo.SetPxPyPzE(.0, .0, .0, .0);
    fCo.SetPxPyPzE(.0, .0, .0, .0);
}

//--------------------------------------------------------------------------

MpdFemtoContainer::MpdFemtoContainer(Int_t ev, TLorentzVector mom, TLorentzVector coord, Float_t phi, Float_t theta, Int_t trackID) {
    fEventNumber = ev;
    fMo = mom;
    fCo = coord;
    fPhi = phi;
    fTheta = theta;
    fTrackID = trackID;
}

//--------------------------------------------------------------------------

MpdFemtoContainer::~MpdFemtoContainer() {

}
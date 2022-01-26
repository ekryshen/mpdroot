//------------------------------------------------------------------------------------------------------------------------
/// \class MpdFfdPoint
///
/// \brief
/// \author Sergei Lobastov (LHE, JINR, Dubna)
//------------------------------------------------------------------------------------------------------------------------

#include <iostream>

#include "MpdFfdPoint.h"

using namespace std;

ClassImp(MpdFfdPoint) MpdFfdPoint::FFDPointMode MpdFfdPoint::fCurrentMode = FFDPointMode::kPhotoElectron;
//------------------------------------------------------------------------------------------------------------------------
MpdFfdPoint::MpdFfdPoint() : FairMCPoint() {}
//------------------------------------------------------------------------------------------------------------------------
MpdFfdPoint::MpdFfdPoint(Int_t tid, Int_t suid) : FairMCPoint(tid, suid, TVector3(), TVector3(), 0., 0., 0.)
{
   fMode = fCurrentMode;
}
//------------------------------------------------------------------------------------------------------------------------
MpdFfdPoint::~MpdFfdPoint()
{
   fData.clear();
}
//------------------------------------------------------------------------------------------------------------------------
void MpdFfdPoint::SaveParentTrackParams(const TVector3 &posIn, const TVector3 &posOut, const TLorentzVector &P,
                                        Double_t time, Double_t length)
{
   SetPosition(posIn);
   fPositionOut = posOut;
   SetMomentum(P.Vect());
   SetTime(time);
   SetLength(length);
   fBeta = P.Beta();

   isClosed = true;
}
//------------------------------------------------------------------------------------------------------------------------
void MpdFfdPoint::Print(const Option_t *opt) const
{
   double step = (fPositionOut - TVector3(fX, fY, fZ)).Mag();
   if (opt) cout << opt;
   cout << "\n-I- MpdFfdPoint(" << isClosed << "): tid=" << fTrackID << ", suid=" << fDetectorID << ", mode=" << fMode
        << ", posIn=(" << fX << ", " << fY << ", " << fZ << ") cm,  posOut=(" << fPositionOut.X() << ", "
        << fPositionOut.Y() << ", " << fPositionOut.Z() << ") cm, time=" << fTime << " ns, trackLength=" << fLength
        << " cm, quartzStep=" << step << ", beta=" << fBeta << " " << this;
   if (step > 0.) cout << ", Nmb per 1cm =" << fData.size() / step;
}
//------------------------------------------------------------------------------------------------------------------------

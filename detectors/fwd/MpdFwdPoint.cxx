//------------------------------------------------------------------------------------------------------------------------
/// \class MpdFwdPoint
///
/// \brief
/// \author Evgeny Kryshen (PNPI, Gatchina)
//------------------------------------------------------------------------------------------------------------------------
#include "MpdFwdPoint.h"
#include "TLorentzVector.h"

ClassImp(MpdFwdPoint)

MpdFwdPoint::MpdFwdPoint() : FairMCPoint() {}

MpdFwdPoint::MpdFwdPoint(Int_t trackID, TLorentzVector posIn, TLorentzVector posOut, Double_t trackLength, Double_t eLoss): 
FairMCPoint(trackID, 0, posIn.Vect(),posOut.Vect(),posIn.T(),trackLength, eLoss,0){}

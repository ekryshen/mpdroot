//------------------------------------------------------------------------------------------------------------------------
/// \class MpdFwdPoint
///
/// \brief
/// \author Evgeny Kryshen (PNPI, Gatchina)
//------------------------------------------------------------------------------------------------------------------------

#ifndef __HH_MPDFWDPOINT_H
#define __HH_MPDFWDPOINT_H

#include "FairMCPoint.h"
#include "TLorentzVector.h"
//------------------------------------------------------------------------------------------------------------------------
class MpdFwdPoint : public FairMCPoint {
public:
   MpdFwdPoint();
   MpdFwdPoint(Int_t trackID, TLorentzVector posIn, TLorentzVector posOut, Double_t trackLength, Double_t eLoss);
   virtual ~MpdFwdPoint(){}
private:
   ClassDef(MpdFwdPoint, 1)
};
//------------------------------------------------------------------------------------------------------------------------
#endif // #ifndef __HH_MPDFWDPOINT_H

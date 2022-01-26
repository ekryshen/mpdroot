/*
 * MpdMcordPoint.h
 *
 *  Created on: 21 maj 2019
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef MCORD_MCORD_MPDMCORDPOINT_H_
#define MCORD_MCORD_MPDMCORDPOINT_H_
#include "TObject.h"
#include "TLorentzVector.h"
#include "FairMCPoint.h"
class MpdMcordPoint : public FairMCPoint {
private:
public:
   MpdMcordPoint();
   MpdMcordPoint(Int_t trackId, Int_t detId, TVector3 pos, TVector3 mom, Double_t time, Double_t length, Double_t eloss,
                 Int_t eventId = 0);
   virtual ~MpdMcordPoint();
   ClassDef(MpdMcordPoint, 1)
};

#endif /* MCORD_MCORD_MPDMCORDPOINT_H_ */

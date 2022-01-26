/*
 * MpdMcordPoint.cxx
 *
 *  Created on: 21 maj 2019
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "MpdMcordPoint.h"

MpdMcordPoint::MpdMcordPoint() : FairMCPoint() {}

MpdMcordPoint::~MpdMcordPoint() {}

MpdMcordPoint::MpdMcordPoint(Int_t trackId, Int_t detId, TVector3 pos, TVector3 mom, Double_t time, Double_t length,
                             Double_t eloss, Int_t eventId)
   : FairMCPoint(trackId, detId, pos, mom, time, length, eloss, eventId)
{
}

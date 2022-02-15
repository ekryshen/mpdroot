/*
 * MpdV0DaughterCut.h
 *
 *  Created on: 16 gru 2021
 *      Author: Daniel Wielanek
 *      E-mail: daniel.wielanek at pw.edu.pl
 *      Warsaw University of Technology, Faculty of Physics
 *      JINR,  Laboratory of High Energy Physics
 */

#ifndef MPDROOT_PHYSICS_COMMON_V0_MPDV0DAUGHTERCUT_H_
#define MPDROOT_PHYSICS_COMMON_V0_MPDV0DAUGHTERCUT_H_

#include <Rtypes.h>
#include <RtypesCore.h>
#include <TVector3.h>

#include "MpdV0FinderCut.h"

class MpdMiniEvent;

class MpdMiniTrack;
class MpdTrack;
class TClonesArray;

class MpdV0DaughterCut : public MpdV0FinderCut {

public:
   MpdV0DaughterCut(){};
   virtual ~MpdV0DaughterCut(){};
   virtual Bool_t PassDstTrack(MpdTrack &track) const         = 0;
   virtual Bool_t PassMiniDstTrack(MpdMiniTrack &track) const = 0;
   MpdV0DaughterCut(const MpdV0DaughterCut &other)            = default;
   MpdV0DaughterCut &operator=(const MpdV0DaughterCut &other) = default;
   ClassDef(MpdV0DaughterCut, 1)
};

#endif /* MPDROOT_PHYSICS_COMMON_V0_MPDV0DAUGHTERCUT_H_ */

/*
 * MpdDCAOnTheFly.h
 *
 *  Created on: 14 lut 2022
 *      Author: Daniel Wielanek
 *      E-mail: daniel.wielanek at pw.edu.pl
 *      Warsaw University of Technology, Faculty of Physics
 *      JINR,  Laboratory of High Energy Physics
 */

#ifndef MPDROOT_PHYSICS_COMMON_V0_MPDDCAONTHEFLY_H_
#define MPDROOT_PHYSICS_COMMON_V0_MPDDCAONTHEFLY_H_

#include <FairTask.h>
#include <Rtypes.h>
#include <RtypesCore.h>

class MpdEvent;

class MpdDCAOnTheFly : public FairTask {
   MpdEvent *fEvent;
   Double_t  fBz;

protected:
   virtual InitStatus Init();

public:
   MpdDCAOnTheFly();
   void         SetBz(Double_t bz) { fBz = bz; };
   virtual void Exec(Option_t *option);
   virtual ~MpdDCAOnTheFly();
   ClassDef(MpdDCAOnTheFly, 1)
};

#endif /* MPDROOT_PHYSICS_COMMON_V0_MPDDCAONTHEFLY_H_ */

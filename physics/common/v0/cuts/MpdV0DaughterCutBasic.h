/*
 * MpdV0DaughterCutBasic.h
 *
 *  Created on: 16 gru 2021
 *      Author: Daniel Wielanek
 *      E-mail: daniel.wielanek at pw.edu.pl
 *      Warsaw University of Technology, Faculty of Physics
 *      JINR,  Laboratory of High Energy Physics
 */

#ifndef MPDROOT_PHYSICS_COMMON_V0_MPDV0DAUGHTERCUTBASIC_H_
#define MPDROOT_PHYSICS_COMMON_V0_MPDV0DAUGHTERCUTBASIC_H_

#include <RtypesCore.h>

#include "MpdV0DaughterCut.h"
#include "MpdV0Namespace.h"

class MpdV0DaughterCutBasic : public MpdV0DaughterCut {
protected:
   Int_t                   fCharge;
   Double_t                fDcaXY2;
   Double_t                fDcaZ;
   Double_t                fSigmaLow;
   Double_t                fSigmaHigh;
   Double_t                fNHitsTpcLow;
   Double_t                fNHitsTpcHigh;
   MpdV0::ESigmaType fSigmaType;

public:
   MpdV0DaughterCutBasic();
   MpdV0DaughterCutBasic(const MpdV0DaughterCutBasic &other) = default;
   MpdV0DaughterCutBasic &operator=(const MpdV0DaughterCutBasic &other) = default;
   void                   SetChargeCut(Int_t charge) { fCharge = charge; };
   void                   SetNTpcHitsCut(Double_t low, Double_t high = 56);
   void                   SetDcaMinCut(Double_t dcaXY, Double_t dcaZ);
   void                   SetSigmaCut(Double_t low, Double_t high, MpdV0::ESigmaType type);
   virtual Bool_t         PassDstTrack(MpdTrack &track) const;
   virtual Bool_t         PassMiniDstTrack(MpdMiniTrack &track) const;
   virtual ~MpdV0DaughterCutBasic();
   ClassDef(MpdV0DaughterCutBasic, 1)
};

#endif /* MPDROOT_PHYSICS_COMMON_V0_MPDV0DAUGHTERCUTBASIC_H_ */

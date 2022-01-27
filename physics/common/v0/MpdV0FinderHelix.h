/*
 * MpdV0FinderHelix.h
 *
 *  Created on: 16 gru 2021
 *      Author: Daniel Wielanek
 *      E-mail: daniel.wielanek at pw.edu.pl
 *      Warsaw University of Technology, Faculty of Physics
 *      JINR,  Laboratory of High Energy Physics
 */

#ifndef MPDROOT_PHYSICS_COMMON_V0_MPDV0FINDERHELIX_H_
#define MPDROOT_PHYSICS_COMMON_V0_MPDV0FINDERHELIX_H_

#include "MpdV0FinderBasic.h"
#include "MpdHelix.h"
#include <vector>

class MpdV0FinderHelix : public MpdV0FinderBasic {
public:
   class NestedHelix : public MpdHelix {
   public:
      NestedHelix() : MpdHelix(){};
      NestedHelix(const MpdHelix &h) : MpdHelix(h){};
      NestedHelix(TVector3 mom, TVector3 o, Double_t charge, Double_t Bz = 0.5) : MpdHelix(mom, o, charge, Bz){};
      virtual ~NestedHelix(){};
   };

protected:
   std::vector<NestedHelix> fFirstHelix;
   std::vector<NestedHelix> fSecondHelix;
   virtual void             ExecMiniDst(Option_t *option);
   virtual void             ExecDst(Option_t *option);
   virtual InitStatus       Init();

public:
   MpdV0FinderHelix(Int_t pidMom = 3122, Int_t pidFirstDau = 211, Int_t pidSecDau = 2212);
   virtual ~MpdV0FinderHelix(){};
   MpdV0FinderHelix(const MpdV0FinderHelix &other) = default;
   MpdV0FinderHelix &operator=(const MpdV0FinderHelix &other) = default;
   ClassDef(MpdV0FinderHelix, 1)
};

#endif /* MPDROOT_PHYSICS_COMMON_V0_MPDV0FINDERHELIX_H_ */

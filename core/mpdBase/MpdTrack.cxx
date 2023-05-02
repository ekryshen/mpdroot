// Author: Oleg Rogachevsky
// Update: 2009-09-17 18:43:28+0400
// Update: 2023-04-24 Dmitri Peresunko
// Copyright: 2009 (C) MPD coll.
//
// Track container

#include "MpdTrack.h"
#include "FairRunAna.h"
#include "FairField.h"

ClassImp(MpdTrack);

MpdHelix MpdTrack::GetHelix() const
{
   TVector3 mom(GetPx(), GetPy(), GetPz());
   TVector3 pos(GetFirstPointX(), GetFirstPointY(), GetFirstPointZ());
   Double_t charge = GetCharge();
   Double_t Bz     = 0.5;
   if (FairRunAna::Instance()) {
      if (FairRunAna::Instance()->GetField()) Bz = FairRunAna::Instance()->GetField()->GetBz(0, 0, 0) * 0.1;
   }
   return MpdHelix(mom, pos, charge, Bz);
}

int MpdTrack::GetNSharedTpcHits() const
{
   // count number of non-zero bits in the map
   int       shared = 0;
   ULong64_t a{fSharedHitMap};
   while (a) {
      shared++;
      a &= a - 1;
   }
   return shared;
}

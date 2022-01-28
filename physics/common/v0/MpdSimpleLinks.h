/*
 * MpdV0Links.h
 *
 *  Created on: 28 sty 2022
 *      Author: Daniel Wielanek
 *      E-mail: daniel.wielanek at pw.edu.pl
 *      Warsaw University of Technology, Faculty of Physics
 *      JINR,  Laboratory of High Energy Physics
 */

#ifndef MPDROOT_PHYSICS_COMMON_V0_MPDSIMPLELINKS_H_
#define MPDROOT_PHYSICS_COMMON_V0_MPDSIMPLELINKS_H_

#include <TNamed.h>
#include <vector>

template <class T>
class MpdSimpleLinks : public TNamed {
   std::vector<T> fArray;

public:
   MpdSimpleLinks(){};
   void        Reserve(Int_t size) { fArray.reserve(size); };
   void        Clear() { fArray.clear(); }
   void        PushBack(T val) { fArray.push_back(val); }
   inline void SetLink(Int_t id, T link) { fArray[id] = link; }
   inline T    GetLink(Int_t id) const { return fArray[id]; }
   Int_t       GetSize() const { return fArray.size(); }
   virtual ~MpdSimpleLinks(){};
   ClassDef(MpdSimpleLinks, 1)
};

#endif /* MPDROOT_PHYSICS_COMMON_V0_MPDSIMPLELINKS_H_ */

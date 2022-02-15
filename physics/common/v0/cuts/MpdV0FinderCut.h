/*
 * MpdV0Cut.h
 *
 *  Created on: 15 lut 2022
 *      Author: Daniel Wielanek
 *      E-mail: daniel.wielanek at pw.edu.pl
 *      Warsaw University of Technology, Faculty of Physics
 *      JINR,  Laboratory of High Energy Physics
 */

#ifndef MPDROOT_PHYSICS_COMMON_V0_CUTS_MPDV0FINDERCUT_H_
#define MPDROOT_PHYSICS_COMMON_V0_CUTS_MPDV0FINDERCUT_H_

#include <Rtypes.h>
#include <RtypesCore.h>
#include <TClonesArray.h>
#include <TObjArray.h>
#include <TVector3.h>

#include "MpdMiniTrack.h"

class MpdEvent;
class MpdMiniEvent;
class MpdMiniBTofPidTraits;
class TClonesArray;

class MpdV0FinderCut : public TObject {
private:
   MpdEvent *    fMpdEvent;   //!
   MpdMiniEvent *fMiniEvent;  //!
   TClonesArray *fMiniTracks; //!
   TClonesArray *fTofTraits;  //!
   TVector3      fVertex;

protected:
   MpdEvent *            GetDstEvent() const { return fMpdEvent; };
   MpdMiniEvent *        GetMiniDstEvent() const { return fMiniEvent; }
   MpdMiniTrack *        GetMiniTrack(Int_t index) const { return (MpdMiniTrack *)fMiniTracks->UncheckedAt(index); }
   MpdMiniBTofPidTraits *GetTofPidTraits(Int_t index) const
   {
      return (MpdMiniBTofPidTraits *)fTofTraits->UncheckedAt(index);
   };
   inline const TVector3 &GetVertex() const { return fVertex; }

public:
   MpdV0FinderCut();
   void SetEventData(MpdEvent *event);
   void SetMiniEventData(MpdMiniEvent *event, TClonesArray *tracks, TClonesArray *tof);
   MpdV0FinderCut(const MpdV0FinderCut &other);
   MpdV0FinderCut &operator=(const MpdV0FinderCut &other);
   virtual ~MpdV0FinderCut();
   ClassDef(MpdV0FinderCut, 1)
};

#endif /* MPDROOT_PHYSICS_COMMON_V0_CUTS_MPDV0FINDERCUT_H_ */

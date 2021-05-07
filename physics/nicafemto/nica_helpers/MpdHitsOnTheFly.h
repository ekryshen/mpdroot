/*
 * MpdHitsOnTheFly.h
 *
 *  Created on: 25 paź 2018
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef MPDROOT_NICA_MPD_NICA_HELPERS_MPDHITSONTHEFLY_H_
#define MPDROOT_NICA_MPD_NICA_HELPERS_MPDHITSONTHEFLY_H_
#include "FairTask.h"
#include "MpdEvent.h"
#include "MpdPid.h"
#include "MpdTrack.h"
#include "TVector3.h"
class MpdHitsOnTheFly : public FairTask {
  MpdEvent *fEvent;
  TClonesArray *fKalmans;
  TClonesArray *fMpdTpcHits;
  Int_t fHitMapSize;
  Int_t *fHitMap;  //[fHitMapSize]
  void GetHitMaps(Int_t particle_index, ULong64_t &layerHit,
                  ULong64_t &sharedHit);

 public:
  MpdHitsOnTheFly();
  virtual InitStatus Init();
  virtual void Exec(Option_t *opt);
  virtual ~MpdHitsOnTheFly();
  ClassDef(MpdHitsOnTheFly, 1)
};

#endif /* MPDROOT_NICA_MPD_NICA_HELPERS_MPDHITSONTHEFLY_H_ */

/*
 * MpdPairSharedHitsCut.h
 *
 *  Created on: 26 paź 2018
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef MPDROOT_NICA_MPD_CUTS_PAIRCUTS_FULL_MPDPAIRSHAREDHITSCUT_H_
#define MPDROOT_NICA_MPD_CUTS_PAIRCUTS_FULL_MPDPAIRSHAREDHITSCUT_H_
#include "MpdTpcPadsPairCut.h"
/**
 * cut that check fraction of shared hits
 */
namespace MpdPadsFormat {
class MpdPairSharedHitsCut : public MpdTpcPadsPairCut {
 public:
  MpdPairSharedHitsCut();
  Bool_t Pass(NicaTwoTrack *pair);
  virtual ~MpdPairSharedHitsCut();
  ClassDef(MpdPairSharedHitsCut, 1)
};
}  // namespace MpdPadsFormat

#endif /* MPDROOT_NICA_MPD_CUTS_PAIRCUTS_FULL_MPDPAIRSHAREDHITSCUT_H_ */

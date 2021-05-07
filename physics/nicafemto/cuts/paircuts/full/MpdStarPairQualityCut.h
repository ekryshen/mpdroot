/*
 * MpdStarQualityPairCut.h
 *
 *  Created on: 30 paź 2018
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef MPDROOT_NICA_MPD_CUTS_PAIRCUTS_TPC_MPDSTARQUALITYPAIRCUT_H_
#define MPDROOT_NICA_MPD_CUTS_PAIRCUTS_TPC_MPDSTARQUALITYPAIRCUT_H_

#include "NicaTwoTrackCut.h"
/**
 * calculate pair quality like in STAR experiment
 */
namespace MpdPadsFormat {
class MpdStarPairQualityCut : public NicaTwoTrackCut {
 public:
  MpdStarPairQualityCut();
  Bool_t Pass(NicaTwoTrack *pair);
  virtual ~MpdStarPairQualityCut();
  ClassDef(MpdStarPairQualityCut, 1)
};
}  // namespace MpdPadsFormat
#endif /* MPDROOT_NICA_MPD_CUTS_PAIRCUTS_TPC_MPDSTARQUALITYPAIRCUT_H_ */

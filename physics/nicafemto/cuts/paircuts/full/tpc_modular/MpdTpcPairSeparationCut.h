/*
 * MpdTpcSeparationCut.h
 *
 *  Created on: 30 paź 2018
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef MPDROOT_NICA_MPD_CUTS_PAIRCUTS_MPDTPCPAIRSEPARATIONCUT_H_
#define MPDROOT_NICA_MPD_CUTS_PAIRCUTS_MPDTPCPAIRSEPARATIONCUT_H_

#include "MpdNominalTpcPairPadsDistanceCut.h"

/**
 * calculate separations between two helices in TPC active area
 */
namespace MpdPadsFormat {
class MpdTpcPairSeparationCut : public MpdNominalTpcPairPadsDistanceCut {
 public:
  MpdTpcPairSeparationCut();
  virtual Bool_t Pass(NicaTwoTrack *pair);
  static Int_t TpcEntrance() { return 0; };
  static Int_t TpcExit() { return 1; };
  static Int_t TpcAverage() { return 2; };
  static Int_t TpcMinimal() { return 3; };
  static Int_t TpcMaximal() { return 4; };
  virtual ~MpdTpcPairSeparationCut();
  ClassDef(MpdTpcPairSeparationCut, 1)
};
}  // namespace MpdPadsFormat
#endif /* MPDROOT_NICA_MPD_CUTS_PAIRCUTS_MPDTPCPAIRSEPARATIONCUT_H_ */

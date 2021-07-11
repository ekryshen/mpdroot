/*
 * MpdAlicePairQuality.h
 *
 *  Created on: 26 paź 2018
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef MPDROOT_NICA_MPD_CUTS_PAIRCUTS_FULL_MPDALICEPAIRQUALITY_H_
#define MPDROOT_NICA_MPD_CUTS_PAIRCUTS_FULL_MPDALICEPAIRQUALITY_H_
#include "MpdTpcPadsPairCut.h"
/**
 * calculate quality of pair like in ALICE experiment
 */
namespace MpdPadsFormat {
class MpdAlicePairQuality : public MpdTpcPadsPairCut {
 public:
  MpdAlicePairQuality();
  Bool_t Pass(NicaTwoTrack *pair);
  virtual ~MpdAlicePairQuality();
  ClassDef(MpdAlicePairQuality, 1)
};
}  // namespace MpdPadsFormat
#endif /* MPDROOT_NICA_MPD_CUTS_PAIRCUTS_FULL_MPDALICEPAIRQUALITY_H_ */

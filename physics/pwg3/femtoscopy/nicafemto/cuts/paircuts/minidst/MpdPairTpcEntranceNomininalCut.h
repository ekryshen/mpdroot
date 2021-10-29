/*
 * MpdPairTpcEntranceNomininalCut.h
 *
 *  Created on: 1 paź 2020
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef MPDROOT_NICA_MPD_CUTS_PAIRCUTS_MINIDST_MPDPAIRTPCENTRANCENOMININALCUT_H_
#define MPDROOT_NICA_MPD_CUTS_PAIRCUTS_MINIDST_MPDPAIRTPCENTRANCENOMININALCUT_H_

#include "MpdFemtoPairCut.h"
namespace MpdHbtDst {
class MpdPairTpcEntranceNomininalCut : public MpdFemtoPairCut {
 public:
  MpdPairTpcEntranceNomininalCut();
  static Int_t Z() { return 0; };
  static Int_t XY() { return 1; };
  static Int_t XYZ() { return 2; }
  virtual Bool_t Pass(NicaTwoTrack *pair);
  virtual ~MpdPairTpcEntranceNomininalCut();
  ClassDef(MpdPairTpcEntranceNomininalCut, 1)
};
}  // namespace MpdHbtDst

#endif /* MPDROOT_NICA_MPD_CUTS_PAIRCUTS_MINIDST_MPDPAIRTPCENTRANCENOMININALCUT_H_ \
        */

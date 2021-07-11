/*
 * MpdSharedPadsCut.h
 *
 *  Created on: 1 paź 2020
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef MPDROOT_NICA_MPD_CUTS_PAIRCUTS_MINIDST_MPDSHAREDPADSCUT_H_
#define MPDROOT_NICA_MPD_CUTS_PAIRCUTS_MINIDST_MPDSHAREDPADSCUT_H_

#include "MpdFemtoPairCut.h"
namespace MpdHbtDst {
class MpdSharedPadsCut : public MpdFemtoPairCut {
 public:
  MpdSharedPadsCut();
  static Int_t SharedPadsNo() { return 0; };
  static Int_t SharedPadsFraction() { return 1; }
  void SetSharedPadsNo(Int_t min, Int_t max) {
    SetMinMax(min, max, SharedPadsNo());
  }
  void SetSharedPadsFraction(Int_t min, Int_t max) {
    SetMinMax(min, max, SharedPadsFraction());
  }
  Bool_t Pass(NicaTwoTrack *pair);
  virtual ~MpdSharedPadsCut();
  ClassDef(MpdSharedPadsCut, 1)
};
}  // namespace MpdHbtDst
#endif /* MPDROOT_NICA_MPD_CUTS_PAIRCUTS_MINIDST_MPDSHAREDPADSCUT_H_ */

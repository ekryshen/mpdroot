/*
 * MpdSharedPadsCut.h
 *
 *  Created on: 30 paź 2018
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef MPDROOT_NICA_MPD_CUTS_PAIRCUTS_PAD_CUTS_MPDSHAREDPADSCUT_H_
#define MPDROOT_NICA_MPD_CUTS_PAIRCUTS_PAD_CUTS_MPDSHAREDPADSCUT_H_

#include "MpdModularTpcPairCut.h"
/**
 * calculate fraction of shared pads
 */
namespace MpdPadsFormat {
class MpdPairSharedPadsCut : public MpdModularTpcPairCut {
  Bool_t fShift;

 public:
  /**
   *
   * @param shift shift all helicities by minus primary vertex
   */
  MpdPairSharedPadsCut(Bool_t shift = kTRUE);
  virtual Bool_t Pass(NicaTwoTrack *pair);
  static Int_t OverlappedPads() { return 0; }
  static Int_t SharedPads() { return 1; };
  virtual ~MpdPairSharedPadsCut();
  ClassDef(MpdPairSharedPadsCut, 1)
};
}  // namespace MpdPadsFormat

#endif /* MPDROOT_NICA_MPD_CUTS_PAIRCUTS_PAD_CUTS_MPDSHAREDPADSCUT_H_ */

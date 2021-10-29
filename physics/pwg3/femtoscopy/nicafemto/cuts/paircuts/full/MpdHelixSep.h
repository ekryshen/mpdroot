/*
 * MpdHelixSep.h
 *
 *  Created on: 18 paź 2018
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef MPDROOT_NICA_MPD_CUTS_PAIRCUTS_FULL_MPDHELIXSEP_H_
#define MPDROOT_NICA_MPD_CUTS_PAIRCUTS_FULL_MPDHELIXSEP_H_
#include "MpdTpcPadsPairCut.h"
/**
 * like MpdTcpEntranceCut but shift particles by DCA
 */
namespace MpdPadsFormat {
class MpdHelixSep : public MpdTpcPadsPairCut {
 public:
  MpdHelixSep();
  virtual Bool_t Pass(NicaTwoTrack *pair);
  virtual ~MpdHelixSep();
  ClassDef(MpdHelixSep, 1)
};
}  // namespace MpdPadsFormat

#endif /* MPDROOT_NICA_MPD_CUTS_PAIRCUTS_FULL_MPDHELIXSEP_H_ */

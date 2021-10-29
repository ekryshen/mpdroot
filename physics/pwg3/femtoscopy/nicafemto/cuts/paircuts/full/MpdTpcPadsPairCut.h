/*
 * MpdTpcPadsPairCut.h
 *
 *  Created on: 1 paź 2020
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef MPDROOT_NICA_MPD_CUTS_PAIRCUTS_FULL_MPDTPCPADSPAIRCUT_H_
#define MPDROOT_NICA_MPD_CUTS_PAIRCUTS_FULL_MPDTPCPADSPAIRCUT_H_

#include "NicaTwoTrackCut.h"
namespace MpdPadsFormat {
class MpdTpcPadsPairCut : public NicaTwoTrackCut {
 public:
  MpdTpcPadsPairCut(Int_t size = 1);
  virtual Bool_t Init(Int_t task_id);
  virtual ~MpdTpcPadsPairCut();
  ClassDef(MpdTpcPadsPairCut, 1);
};
}  // namespace MpdPadsFormat
#endif /* MPDROOT_NICA_MPD_CUTS_PAIRCUTS_FULL_MPDTPCPADSPAIRCUT_H_ */

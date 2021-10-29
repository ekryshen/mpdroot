/*
 * MpdFemtoPairCut.h
 *
 *  Created on: 1 paź 2020
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef MPDROOT_NICA_MPD_CUTS_PAIRCUTS_MINIDST_MPDFEMTOPAIRCUT_H_
#define MPDROOT_NICA_MPD_CUTS_PAIRCUTS_MINIDST_MPDFEMTOPAIRCUT_H_

#include "NicaTwoTrackCut.h"
namespace MpdHbtDst {
class MpdFemtoPairCut : public NicaTwoTrackCut {
 public:
  MpdFemtoPairCut(Int_t size = 1);
  virtual Bool_t Init(Int_t task_id);
  virtual ~MpdFemtoPairCut();
  ClassDef(MpdFemtoPairCut, 1)
};
}  // namespace MpdFemtoMiniDst

#endif /* MPDROOT_NICA_MPD_CUTS_PAIRCUTS_MINIDST_MPDFEMTOPAIRCUT_H_ */

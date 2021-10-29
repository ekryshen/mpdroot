/*
 * MpdRejectSPlittedPairsCut.h
 *
 *  Created on: 28 lip 2018
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef MPDROOT_NICA_MPD_CUTS_PAIRCUTS_COMMON_MPDSPLITTEDPAIRSCUT_H_
#define MPDROOT_NICA_MPD_CUTS_PAIRCUTS_COMMON_MPDSPLITTEDPAIRSCUT_H_
#include "NicaComplexTrack.h"
#include "NicaTwoTrackCut.h"

/**
 * by default remove pairs with the same matched ID (parents that have the same
 * mathed MC parent)
 */
class MpdSplittedPairsCut : public NicaTwoTrackCut {
  Bool_t fReject;

 public:
  MpdSplittedPairsCut();
  virtual Bool_t Init(Int_t task_id);
  virtual Bool_t Pass(NicaTwoTrack *pair);
  void AcceptSplitted() { fReject = kFALSE; };
  virtual ~MpdSplittedPairsCut();
  ClassDef(MpdSplittedPairsCut, 1)
};

#endif /* MPDROOT_NICA_MPD_CUTS_PAIRCUTS_COMMON_MPDSPLITTEDPAIRSCUT_H_ */

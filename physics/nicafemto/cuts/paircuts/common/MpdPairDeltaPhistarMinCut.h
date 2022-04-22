/*
 * MpdPairDeltaPhistarMinCut.h
 *
 *  Created on: 22 lip 2021
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef MPDROOT_PHYSICS_NICAFEMTO_CUTS_PAIRCUTS_COMMON_MPDPAIRDELTAPHISTARMINCUT_H_
#define MPDROOT_PHYSICS_NICAFEMTO_CUTS_PAIRCUTS_COMMON_MPDPAIRDELTAPHISTARMINCUT_H_


#include "NicaTwoTrackCut.h"

class MpdPairDeltaPhistarMinCut : public NicaTwoTrackCut {
  Double_t fStep;

public:
  MpdPairDeltaPhistarMinCut();
  void SetNSteps(Int_t steps);
  virtual Bool_t Pass(NicaTwoTrack* pair);
  virtual NicaPackage* Report() const;
  virtual ~MpdPairDeltaPhistarMinCut();
  ClassDef(MpdPairDeltaPhistarMinCut, 1);
};


#endif /* MPDROOT_PHYSICS_NICAFEMTO_CUTS_PAIRCUTS_COMMON_MPDPAIRDELTAPHISTARMINCUT_H_ */

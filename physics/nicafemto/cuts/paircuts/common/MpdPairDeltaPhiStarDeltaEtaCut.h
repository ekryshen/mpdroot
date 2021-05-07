/*
 * MpdPairDeltaPhiStarDeltaEtaCut.h
 *
 *  Created on: 1 paź 2020
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef MPDPAIRDELTAPHISTARDELTAETACUT_H_
#define MPDPAIRDELTAPHISTARDELTAETACUT_H_

#include "NicaTwoTrackCut.h"

/**
 * calculate delta phi delta eta value by using semi-analytical formula,
 * base on AliFemtoPairCutRadialDistance
 */
class MpdPairDeltaPhiStarDeltaEtaCut : public NicaTwoTrackCut {
  Double_t fR;


public:
  MpdPairDeltaPhiStarDeltaEtaCut();
  MpdPairDeltaPhiStarDeltaEtaCut(const MpdPairDeltaPhiStarDeltaEtaCut& other);
  void SetDeltaEtaCut(Double_t min, Double_t max);
  void SetDeltaPhiStarCut(Double_t min, Double_t max);
  /**
   *
   * @param R radius in [cm]
   */
  void SetR(Double_t R) { fR = R * 0.01; };
  Bool_t Init(Int_t task_id);
  virtual Bool_t Pass(NicaTwoTrack* pair);
  NicaPackage* Report() const;
  virtual ~MpdPairDeltaPhiStarDeltaEtaCut();
  ClassDef(MpdPairDeltaPhiStarDeltaEtaCut, 1)
};

class MpdPairDeltaPhiStarDeltaEtaMinCut : public MpdPairDeltaPhiStarDeltaEtaCut {
public:
  MpdPairDeltaPhiStarDeltaEtaMinCut() {};
  virtual Bool_t Pass(NicaTwoTrack* pair);
  NicaPackage* Report() const;
  virtual ~MpdPairDeltaPhiStarDeltaEtaMinCut() {};
  ClassDef(MpdPairDeltaPhiStarDeltaEtaMinCut, 1)
};

#endif /* MPDROOT_NICA_MPD_CUTS_PAIRCUTS_COMMON_MPDPAIRDELTAPHISTARDELTAETACUT_H_                                                \
        */

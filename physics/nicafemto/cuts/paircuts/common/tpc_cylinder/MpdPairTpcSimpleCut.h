/*
 * MpdPairTpcSimpleCut.h
 *
 *  Created on: 22 lis 2018
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef MPDROOT_CYLINDERMPDPAIRTPCSIMPLECUT_H_
#define MPDROOT_CYLINDERMPDPAIRTPCSIMPLECUT_H_

#include "TVector3.h"

#include "MpdCylinderTpcPairCut.h"
#include "NicaExpTrack.h"

class MpdPairTpcSimpleCut : public MpdCylinderTpcPairCut {
 protected:
  TVector3 PosNominal(NicaExpTrackHelix *track, Double_t R) const;
  TVector3 PosReal(NicaExpTrackHelix *track, Double_t R) const;

 public:
  MpdPairTpcSimpleCut(Int_t size = 1);
  virtual ~MpdPairTpcSimpleCut();
  ClassDef(MpdPairTpcSimpleCut, 1)
};

#endif /* MPDROOT_CYLINDERMPDPAIRTPCSIMPLECUT_H_ */

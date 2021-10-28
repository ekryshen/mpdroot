/*
 * MpdFemtoDistAtLayerCut.h
 *
 *  Created on: 21 lip 2021
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef MPDROOT_PHYSICS_NICAFEMTO_CUTS_PAIRCUTS_MINIDST_MPDFEMTODISTATLAYERCUT_H_
#define MPDROOT_PHYSICS_NICAFEMTO_CUTS_PAIRCUTS_MINIDST_MPDFEMTODISTATLAYERCUT_H_


#include "MpdFemtoPairCut.h"

namespace MpdHbtDst {
  class MpdFemtoDistAtLayerCut : public MpdFemtoPairCut {
    Int_t fLayerNo;

  public:
    MpdFemtoDistAtLayerCut();
    void SetLayerNo(Int_t lay);
    void SetXYCut(Double_t min, Double_t max) { SetMinMax(min, max, 0); };
    void SetZCut(Double_t min, Double_t max) { SetMinMax(min, max, 1); };
    void SetXYZCut(Double_t min, Double_t max) { SetMinMax(min, max, 2); };
    Bool_t Pass(NicaTwoTrack* pair);
    virtual NicaPackage* Report() const;
    virtual ~MpdFemtoDistAtLayerCut();
  };

} /* namespace MpdHbtDst */

#endif /* MPDROOT_PHYSICS_NICAFEMTO_CUTS_PAIRCUTS_MINIDST_MPDFEMTODISTATLAYERCUT_H_ */

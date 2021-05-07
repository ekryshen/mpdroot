/*
 * NicaMpdConst.h
 *
 *  Created on: 1 paź 2020
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef MPDROOT_NICA_MPD_NICA_HELPERS_NICAMPDCONST_H_
#define MPDROOT_NICA_MPD_NICA_HELPERS_NICAMPDCONST_H_

#include "NicaDataFormat.h"
#include <Rtypes.h>

namespace NicaMpdConst {
  const Double_t TpcInnerRadius      = 27;
  const Double_t TpcOuterRadius      = 140;
  const Double_t TpcInnerDriftRadius = 24;
  const Double_t TPcOuterDriftRadius = 133;
  const Double_t TpcLenght           = 340;
  const Double_t TpcDriftLenght      = 326;
  const Int_t TpcLayers              = 53;
}  // namespace NicaMpdConst

namespace MpdDetectorID {
  const UInt_t kTPC  = NicaDetectorID::kTPC;
  const UInt_t kTOF  = NicaDetectorID::kTOF;
  const UInt_t kSTS  = NicaDetectorID::kSTS;
  const UInt_t kZDC  = 45840;
  const UInt_t kEMC  = 18948;
  const UInt_t kCPC  = 16464;
  const UInt_t kETOF = 691647;
  const UInt_t kFFD  = 19993;
};     // namespace MpdDetectorID
#endif /* MPDROOT_NICA_MPD_NICA_HELPERS_NICAMPDCONST_H_ */

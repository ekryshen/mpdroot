/*
 * MpdPairCutModularTpc.h
 *
 *  Created on: 22 lis 2018
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef MPDROOT_NICA_MPD_CUTS_PAIRCUTS_FULL_TPC_MODULAR_MPDMODULARTPCPAIRCUT_H_
#define MPDROOT_NICA_MPD_CUTS_PAIRCUTS_FULL_TPC_MODULAR_MPDMODULARTPCPAIRCUT_H_

#include "MpdTpcPadsPairCut.h"
namespace MpdPadsFormat {
class MpdModularTpcPairCut : public MpdTpcPadsPairCut {
 public:
  MpdModularTpcPairCut(Int_t size);
  MpdModularTpcPairCut& operator=(const MpdModularTpcPairCut& other);
  virtual ~MpdModularTpcPairCut();
  ClassDef(MpdModularTpcPairCut, 1)
};
}  // namespace MpdPadsFormat

#endif /* MPDROOT_NICA_MPD_CUTS_PAIRCUTS_FULL_TPC_MODULAR_MPDMODULARTPCPAIRCUT_H_ \
        */

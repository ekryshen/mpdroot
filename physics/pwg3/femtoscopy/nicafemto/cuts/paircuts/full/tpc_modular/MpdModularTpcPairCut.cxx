/*
 * MpdPairCutModularTpc.cxx
 *
 *  Created on: 22 lis 2018
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "MpdModularTpcPairCut.h"
#include "NicaTwoTrack.h"

namespace MpdPadsFormat {
MpdModularTpcPairCut::MpdModularTpcPairCut(Int_t size)
    : MpdTpcPadsPairCut(size) {}

MpdModularTpcPairCut& MpdModularTpcPairCut::operator=(
    const MpdModularTpcPairCut& other) {
  if (&other == this) return *this;
  MpdTpcPadsPairCut::operator=(other);
  return *this;
}

MpdModularTpcPairCut::~MpdModularTpcPairCut() {}
}  // namespace MpdPadsFormat

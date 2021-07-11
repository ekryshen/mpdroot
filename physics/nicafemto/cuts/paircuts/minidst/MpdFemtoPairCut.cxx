/*
 * MpdFemtoPairCut.cxx
 *
 *  Created on: 1 paź 2020
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "MpdFemtoPairCut.h"
#include "NicaDataFormatManager.h"
namespace MpdHbtDst {
MpdFemtoPairCut::MpdFemtoPairCut(Int_t size) : NicaTwoTrackCut(size) {}

MpdFemtoPairCut::~MpdFemtoPairCut() {}

Bool_t MpdFemtoPairCut::Init(Int_t task_id) {
  return FormatInhertis("NicaMpdHbtEvent", task_id,
                        ENicaFormatDepth::kBuffered);
}
}  // namespace MpdHbtDst

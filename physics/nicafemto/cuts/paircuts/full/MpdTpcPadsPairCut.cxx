/*
 * MpdTpcPadsPairCut.cxx
 *
 *  Created on: 1 paź 2020
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "MpdTpcPadsPairCut.h"
namespace MpdPadsFormat {
MpdTpcPadsPairCut::MpdTpcPadsPairCut(Int_t size) : NicaTwoTrackCut(size) {}

MpdTpcPadsPairCut::~MpdTpcPadsPairCut() {}

Bool_t MpdPadsFormat::MpdTpcPadsPairCut::Init(Int_t task_id) {
  return FormatInhertis("NicaMpdEventTpcPads", task_id,
                        ENicaFormatDepth::kBuffered);
}
}  // namespace MpdPadsFormat

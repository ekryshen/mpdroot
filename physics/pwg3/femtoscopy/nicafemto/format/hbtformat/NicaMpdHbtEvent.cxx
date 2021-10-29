/*
 * MpdHbtEvent.cxx
 *
 *  Created on: 1 paź 2020
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */

#include "NicaMpdHbtEvent.h"


#include "NicaCout.h"
#include "NicaStd.h"

NicaMpdHbtEvent::NicaMpdHbtEvent() : NicaExpEvent("NicaMpdHbtTrack") {}

Bool_t NicaMpdHbtEvent::IsCompatible(const NicaEvent* buffered) const {
  if (buffered->InheritsFrom("NicaMpdMiniDstEvent")) return kTRUE;
  return kFALSE;
}

NicaPackage* NicaMpdHbtEvent::Report() const {
  NicaPackage* rep = NicaEvent::Report();
  return rep;
}

NicaMpdHbtEvent::NicaMpdHbtEvent(const NicaMpdHbtEvent& other)
    : NicaExpEvent(other) {}

NicaMpdHbtEvent& NicaMpdHbtEvent::operator=(const NicaMpdHbtEvent& other) {
  if (this == &other) return *this;
  NicaExpEvent::operator=(other);
  return *this;
}

void NicaMpdHbtEvent::CreateSource() {
  NicaCout::PrintInfo("Cannot create source for NicaMpdHbtEvent",
                      ENicaInfo::kLessError);
}

NicaMpdHbtEvent::~NicaMpdHbtEvent() {}

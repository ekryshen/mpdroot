/*
 * MpdV0Cut.cxx
 *
 *  Created on: 15 lut 2022
 *      Author: Daniel Wielanek
 *      E-mail: daniel.wielanek at pw.edu.pl
 *      Warsaw University of Technology, Faculty of Physics
 *      JINR,  Laboratory of High Energy Physics
 */

#include "MpdV0FinderCut.h"

#include "MpdEvent.h"
#include "MpdMiniEvent.h"

#include <TClonesArray.h>

MpdV0FinderCut::MpdV0FinderCut() : fMpdEvent(nullptr), fMiniEvent(nullptr), fMiniTracks(nullptr), fTofTraits(nullptr) {}

MpdV0FinderCut::MpdV0FinderCut(const MpdV0FinderCut &other) : MpdV0FinderCut() {}

MpdV0FinderCut &MpdV0FinderCut::operator=(const MpdV0FinderCut &other)
{
   if (this == &other) return *this;
   fMpdEvent   = nullptr;
   fMiniEvent  = nullptr;
   fMiniTracks = nullptr;
   fTofTraits  = nullptr;
   fVertex.SetXYZ(0, 0, 0);
   return *this;
}

MpdV0FinderCut::~MpdV0FinderCut() {}

void MpdV0FinderCut::SetEventData(MpdEvent *event)
{
   fMpdEvent = event;
   fVertex.SetXYZ(fMpdEvent->GetPrimaryVerticesX(), fMpdEvent->GetPrimaryVerticesY(), fMpdEvent->GetPrimaryVerticesZ());
}

void MpdV0FinderCut::SetMiniEventData(MpdMiniEvent *event, TClonesArray *tracks, TClonesArray *tof)
{
   fMiniEvent  = event;
   fTofTraits  = tof;
   fMiniTracks = tracks;
   fVertex     = event->primaryVertex();
}

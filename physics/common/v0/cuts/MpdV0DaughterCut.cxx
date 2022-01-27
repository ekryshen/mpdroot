/*
 * MpdV0DaughterCut.cxx
 *
 *  Created on: 16 gru 2021
 *      Author: Daniel Wielanek
 *      E-mail: daniel.wielanek at pw.edu.pl
 *      Warsaw University of Technology, Faculty of Physics
 *      JINR,  Laboratory of High Energy Physics
 */

#include "MpdV0DaughterCut.h"

#include "MpdMiniEvent.h"

void MpdV0DaughterCut::SetMiniDstEventInfo(MpdMiniEvent &event)
{
   fVertex = event.primaryVertex();
}

/*
 * MpdV0Track.cxx
 *
 *  Created on: 16 gru 2021
 *      Author: Daniel Wielanek
 *      E-mail: daniel.wielanek at pw.edu.pl
 *      Warsaw University of Technology, Faculty of Physics
 *      JINR,  Laboratory of High Energy Physics
 */

#include "MpdV0Track.h"

MpdV0Track::MpdV0Track() : MpdV0Particle(), fPositiveDaughterS(0), fNegativeDaughterS(0), fChi2(-1) {}

MpdV0Track::~MpdV0Track() {}

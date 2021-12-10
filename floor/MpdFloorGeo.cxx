/*
 * MpdFloorGeo.cxx
 *
 *  Created on: 21 maj 2019
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "MpdFloorGeo.h"

MpdFloorGeo::MpdFloorGeo() : FairGeoSet() {
  fName      = "fl";
  maxSectors = 128;
  maxModules = 500;
}

MpdFloorGeo::~MpdFloorGeo() {
  // TODO Auto-generated destructor stub
}

const char* MpdFloorGeo::getModuleName(Int_t m) {
  sprintf(modName, "fl%02d", m + 1);
  return modName;
}

const char* MpdFloorGeo::getEleName(Int_t m) {
  sprintf(eleName, "fl%02d", m + 1);
  return eleName;
}

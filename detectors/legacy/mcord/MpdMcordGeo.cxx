/*
 * MpdMcordGeo.cxx
 *
 *  Created on: 21 maj 2019
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "MpdMcordGeo.h"

MpdMcordGeo::MpdMcordGeo() : FairGeoSet()
{
   fName      = "mcd";
   maxSectors = 128;
   maxModules = 500;
}

MpdMcordGeo::~MpdMcordGeo()
{
   // TODO Auto-generated destructor stub
}

const char *MpdMcordGeo::getModuleName(Int_t m)
{
   sprintf(modName, "md%02d", m + 1);
   return modName;
}

const char *MpdMcordGeo::getEleName(Int_t m)
{
   sprintf(eleName, "md%02d", m + 1);
   return eleName;
}

#include "FairGeoCradle.h"

#include "TString.h"
#include <string.h>

FairGeoCradle::FairGeoCradle() : FairGeoSet()
{

   fName = "cradle";
   strcpy(modName, "c");
   strcpy(eleName, "c");
   maxSectors = 0;
   maxModules = 1;
}

ClassImp(FairGeoCradle);
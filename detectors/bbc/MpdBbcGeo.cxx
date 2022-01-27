/////////////////////////////////////////////////////////////
// MpdBbcGeo
//
// Class for geometry of BBC
//
/////////////////////////////////////////////////////////////

#include "MpdBbcGeo.h"
#include "FairGeoNode.h"

ClassImp(MpdBbcGeo);

MpdBbcGeo::MpdBbcGeo()
{
   // Constructor
   fName      = "bbc";
   maxSectors = 0;
   maxModules = 4;
}

const char *MpdBbcGeo::getModuleName(Int_t m)
{
   // Returns the module name of bbc number m
   sprintf(modName, "bbc0%i", m + 1);
   return modName;
}

const char *MpdBbcGeo::getEleName(Int_t m)
{
   // Returns the element name of bbc number m
   sprintf(eleName, "bbc0%i", m + 1);
   return eleName;
}

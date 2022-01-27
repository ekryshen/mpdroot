/////////////////////////////////////////////////////////////
// MpdStsGeo
//
// Class for geometry of STS
//
/////////////////////////////////////////////////////////////

#include "MpdStsGeo.h"
#include "FairGeoNode.h"

ClassImp(MpdStsGeo);

MpdStsGeo::MpdStsGeo()
{
   // Constructor
   fName      = "sts";
   maxSectors = 0;
   maxModules = 4;
}

const char *MpdStsGeo::getModuleName(Int_t m)
{
   // Returns the module name of sts number m
   sprintf(modName, "sts0%i", m + 1);
   return modName;
}

const char *MpdStsGeo::getEleName(Int_t m)
{
   // Returns the element name of sts number m
   sprintf(eleName, "sts0%i", m + 1);
   return eleName;
}

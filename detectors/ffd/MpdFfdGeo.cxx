/////////////////////////////////////////////////////////////
// MpdFfdGeo
//
// Class for geometry of FFD
//
/////////////////////////////////////////////////////////////

#include "MpdFfdGeo.h"
#include "FairGeoNode.h"

ClassImp(MpdFfdGeo)

   MpdFfdGeo::MpdFfdGeo()
{
   // Constructor
   fName      = "ffd";
   maxSectors = 0;
   maxModules = 4;
}

const char *MpdFfdGeo::getModuleName(Int_t m)
{
   // Returns the module name of cpc number m
   sprintf(modName, "Ffd0%i", m + 1);
   return modName;
}

const char *MpdFfdGeo::getEleName(Int_t m)
{
   // Returns the element name of cpc number m
   sprintf(eleName, "Ffd0%i", m + 1);
   return eleName;
}

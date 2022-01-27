#include "FairGeoNode.h"

#include "MpdTofGeo.h"
//------------------------------------------------------------------------------------------------------------------------
MpdTofGeo::MpdTofGeo()
{
   fName      = "tof";
   maxSectors = 0;
   maxModules = 124;
}
//------------------------------------------------------------------------------------------------------------------------
const char *MpdTofGeo::getModuleName(Int_t m)
{
   // Returns the module name of tof number m
   sprintf(modName, "tof%i", m + 1);
   return modName;
}
//------------------------------------------------------------------------------------------------------------------------
const char *MpdTofGeo::getEleName(Int_t m)
{
   // Returns the element name of tof number m
   sprintf(eleName, "t%i", m + 1);
   return eleName;
}
//------------------------------------------------------------------------------------------------------------------------
ClassImp(MpdTofGeo);

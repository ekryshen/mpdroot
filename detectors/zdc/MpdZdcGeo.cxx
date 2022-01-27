/*************************************************************************************
 *
 *         Class MpdZdcGeo
 *
 *  Adopted for MPD by:   Elena Litvinenko
 *  e-mail:   litvin@nf.jinr.ru
 *  Version:  8-Apr-2008
 *
 ************************************************************************************/

#include <iostream>

#include "MpdZdcGeo.h"
#include "FairGeoNode.h"
#include <iostream>

using std::cout;
using std::endl;

ClassImp(MpdZdcGeo);

// -----   Default constructor   -------------------------------------------
MpdZdcGeo::MpdZdcGeo()
{
   // Constructor
   fName      = "zdc";
   maxSectors = 0;
   maxModules = 4;
}
// -------------------------------------------------------------------------

const char *MpdZdcGeo::getModuleName(Int_t m)
{
   // Returns the module name of muo number m

   sprintf(modName, "zdc0%i", m + 1);
   return modName;
   cout << "MODNAME: " << modName << endl;
}

const char *MpdZdcGeo::getEleName(Int_t m)
{
   // Returns the element name of muo number m

   sprintf(eleName, "s%i", m + 1);
   return eleName;
   cout << "ELENAME: " << eleName << endl;
}

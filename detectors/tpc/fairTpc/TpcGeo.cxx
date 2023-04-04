//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      Implementation of class TpcGeo
//      see TpcGeo.hh for details
//
// Environment:
//      Software developed for the PANDA Detector at FAIR.
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------

// Panda Headers ----------------------

// This Class' Header ------------------
#include "TpcGeo.h"

// C/C++ Headers ----------------------

// Collaborating Class Headers --------
#include "FairGeoNode.h"

// Class Member definitions -----------

ClassImp(TpcGeo);

// -----   Default constructor   -------------------------------------------
TpcGeo::TpcGeo()
{
   // Constructor
   fName      = "tpc";
   maxSectors = 9;
   maxModules = 9;
}

// -------------------------------------------------------------------------

const char *TpcGeo::getModuleName(Int_t m)
{

   // Returns the module name of Tpc number m

   // Setting tpc here means that all modules names in the ASCII file should start with tpc otherwise they will not be
   // constructed

   sprintf(modName, "tpcChamber%i", m + 1);
   return modName;
}

const char *TpcGeo::getEleName(Int_t m)
{
   // Returns the element name of Det number m
   sprintf(eleName, "tpc0%i", m + 1);
   return eleName;
}

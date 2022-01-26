//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      Implementation of class BmdGeo
//      see BmdGeo.hh for details
//
// Environment:
//      Software developed for the PANDA Detector at FAIR.
//
// Author List:
//      Pedro González    TUM            (original author)
//
//
//-----------------------------------------------------------
// Panda Headers --------------------------------------------

// This Class' Header ------------------
#include "BmdGeo.h"

// C/C++ Headers ----------------------


// Collaborating Class Headers --------
#include "FairGeoNode.h"

// Class Member definitions -----------

ClassImp(BmdGeo)

// -----   Default constructor   -------------------------------------------
BmdGeo::BmdGeo() {
  // Constructor
  fName="bmd";
  maxSectors=2;
  maxModules=160;
}

// -------------------------------------------------------------------------

const char* BmdGeo::getModuleName(Int_t m) {

  // Returns the module name of Bmd number m

  // Setting tpc here means that all modules names in the ASCII file should start with tpc otherwise they will not be constructed


  //sprintf(modName,"bmdChamber%i",m+1); 
  sprintf(modName,"BMD");
  return modName;
}

const char* BmdGeo::getEleName(Int_t m) {
    
    printf("Hola entro a BmdGeo");
  // Returns the element name of Det number m
  if( m == 0 )
  sprintf(eleName,"BMDPLUSA");
  else if ( m == 1 )
  sprintf(eleName,"BMDPLUSC");  
  return eleName;
}

 


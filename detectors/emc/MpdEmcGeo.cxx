#include "MpdEmcGeo.h"
#include "FairGeoNode.h"

ClassImp(MpdEmcGeo)

// -----   Default constructor   -------------------------------------------
MpdEmcGeo::MpdEmcGeo() {
  // Constructor
  fName="emc";
  maxSectors=0;
  maxModules=10;
}

// -------------------------------------------------------------------------

const char* MpdEmcGeo::getModuleName(Int_t m) {
  /** Returns the module name of TutorialDet number m
      Setting MyDet here means that all modules names in the 
      ASCII file should start with TutorialDet otherwise they will 
      not be constructed
  */
  sprintf(modName,"emc%i",m+1); 
  return modName;
}

const char* MpdEmcGeo::getEleName(Int_t m) {
  /** Returns the element name of Det number m */
  sprintf(eleName,"emc%i",m+1);
  return eleName;
}

//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      Bmd geometry module
//
//
// Environment:
//      Software developed for the PANDA Detector at FAIR.
//
// Author List:
//      Pedro Gonzalez    TUM            (original author)
//
//
//-----------------------------------------------------------

#ifndef BMDGEO_HH
#define BMDGEO_HH

// Base Class Headers ----------------
#include "FairGeoSet.h"

// Collaborating Class Headers -------

// Collaborating Class Declarations --



class BmdGeo : public FairGeoSet {
public:

  // Constructors/Destructors ---------
   BmdGeo();
  ~BmdGeo(){}

  // Operators
  

  // Accessors -----------------------
  const char* getModuleName(Int_t);
  const char* getEleName(Int_t);
  //inline Int_t getModNumInMod(const TString&);


  // Modifiers -----------------------


  // Operations ----------------------
protected:
  char modName[20];  // name of module
  char eleName[20];  // substring for elements in module
private:

  // Private Data Members ------------


  // Private Methods -----------------

public:
  ClassDef(BmdGeo,1)

};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------

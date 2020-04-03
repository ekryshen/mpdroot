//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      Monte Carlo Point in the MBB
//
//
// Environment:
//      Software developed for the PANDA Detector at FAIR.
//
// Author List:
//     Pedro González    TUM            (original author)
//     Luis Valenzuela-Cazares
//
//-----------------------------------------------------------

#ifndef MBBPOINT_HH
#define MBBPOINT_HH

// Base Class Headers ----------------
#include "FairMCPoint.h"

// Collaborating Class Headers -------
#include <ostream> // remove if you do not need streaming op

// Collaborating Class Declarations --



class MbbPoint : public FairMCPoint {
public:

  // Constructors/Destructors ---------
  MbbPoint();
  MbbPoint(Int_t trackID, Int_t detID, TVector3 pos, TVector3 mom,
             Double_t tof, Double_t length, Double_t eLoss, Int_t statusCode);
  virtual ~MbbPoint(){;}

  // Operators


  // Accessors -----------------------


  // Modifiers -----------------------


  // Operations ----------------------
  virtual void Print(const Option_t* opt=0) const ;
  Double_t GetStep() const { return fStep; } ///> return step length
  void SetStep(Double_t step) { fStep = step; } ///> set step length
  Int_t GetStatusCode(){return fStatusCode;}

private:

  // Private Data Members ------------
  Double32_t fStep; // step length
  Double_t  pT;
  Double_t  eta;
  Int_t     fStatusCode;

  //Int_t     pdgCode;


  // Private Methods -----------------

public:
  ClassDef(MbbPoint,1)

};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------

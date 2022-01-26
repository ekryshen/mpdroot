//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      Monte Carlo Point in the BMD
//
//
// Environment:
//      Software developed for the PANDA Detector at FAIR.
//
// Author List:
//     Pedro González    TUM            (original author)
//
//
//-----------------------------------------------------------

#ifndef BMDPOINT_HH
#define BMDPOINT_HH

// Base Class Headers ----------------
#include "FairMCPoint.h"

// Collaborating Class Headers -------
#include <ostream> // remove if you do not need streaming op

// Collaborating Class Declarations --



class BmdPoint : public FairMCPoint {
public:

  // Constructors/Destructors ---------
  BmdPoint();
  BmdPoint(Int_t trackID, Int_t detID, TVector3 pos, TVector3 mom,
             Double_t tof, Double_t length, Double_t eLoss);
  virtual ~BmdPoint(){;}

  // Operators
  

  // Accessors -----------------------


  // Modifiers -----------------------


  // Operations ----------------------
  virtual void Print(const Option_t* opt=0) const ;
  Double_t GetStep() const { return fStep; } ///> return step length
  void SetStep(Double_t step) { fStep = step; } ///> set step length

private:

  // Private Data Members ------------
  Double32_t fStep; // step length
  Double_t  pT;
  Double_t  eta;
  //Int_t     pdgCode;
  

  // Private Methods -----------------

public:
  ClassDef(BmdPoint,1)

};

#endif

//--------------------------------------------------------------
// $Log$
//--------------------------------------------------------------


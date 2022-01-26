//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      Implementation of class BmdPoint
//      see BmdPoint.hh for details
//
// Environment:
//      Software developed for the PANDA Detector at FAIR.
//
// Author List:
//      Pedro Gonzalez    TUM            (original author)
//
//
//-----------------------------------------------------------

// Panda Headers ----------------------

// This Class' Header ------------------
#include "BmdPoint.h"

// C/C++ Headers ----------------------
#include <iostream>

// Collaborating Class Headers --------


// Class Member definitions -----------

BmdPoint::BmdPoint()
  : FairMCPoint()
{}

BmdPoint::BmdPoint(Int_t trackID, Int_t detID, TVector3 pos, TVector3 mom,
             Double_t tof, Double_t length, Double_t eLoss)
  : FairMCPoint(trackID, detID, pos, mom, tof, length, eLoss)
{
  
  pT = mom.Pt();
  eta = mom.Eta();
  
  
  
}

void 
BmdPoint::Print(const Option_t*) const {
  std::cout<<"BmdPoint\n"
	   <<" Pos("<<fX<<","<<fY<<","<<fZ<<")\n"
	   <<" dE="<< fELoss << fTrackID
	   <<std::endl;	   
}



ClassImp(BmdPoint)

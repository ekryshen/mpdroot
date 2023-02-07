//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      Implementation of class TpcPoint
//      see TpcPoint.hh for details
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
#include "TpcPoint.h"

// C/C++ Headers ----------------------
#include <iostream>

// Collaborating Class Headers --------

// Class Member definitions -----------

TpcPoint::TpcPoint() : FairMCPoint() {}

TpcPoint::TpcPoint(Int_t trackID, Int_t detID, TVector3 pos, TVector3 mom, Double_t tof, Double_t length,
                   Double_t eLoss)
   : FairMCPoint(trackID, detID, pos, mom, tof, length, eLoss)
{
}

void TpcPoint::Print(const Option_t *) const
{
   std::cout << "TpcPoint\n"
             << " Pos(" << fX << "," << fY << "," << fZ << ")\n"
             << " dE=" << fELoss << fTrackID << std::endl;
}

ClassImp(TpcPoint);

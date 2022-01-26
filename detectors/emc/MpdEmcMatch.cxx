//-----------------------------------------------------------
//
// Description:
//      Implementation of class MpdEmcMatch
//      see MpdEmcMatch.h for details
//
// Author List:
//      Alexander Zinchenko LHEP, JINR, Dubna - 2-June-2016
//
//-----------------------------------------------------------

#include "MpdEmcMatch.h"

#include <iostream>

using namespace std;

// -----   Default constructor   -------------------------------------------

MpdEmcMatch::MpdEmcMatch() : TObject(),
  fEnergy(0.0), fChi2(999.0)

{
  fDist[0] = fDist[1] = 0.0;
  fInds[0] = fInds[1] = -1;
}


// -----   Standard constructor   ------------------------------------------

MpdEmcMatch::MpdEmcMatch(Int_t trackInd, Int_t recpInd, Double_t distT, Double_t distL, 
			 Double_t ener, Double_t trLeng) : TObject(),
  fEnergy(ener), fTrLeng(trLeng), fChi2(999.0), fChi2pi(999.0)
{
  fDist[0] = distT;
  fDist[1] = distL;
  fInds[0] = trackInd;
  fInds[1] = recpInd;
}

// -----   Destructor   ----------------------------------------------------

MpdEmcMatch::~MpdEmcMatch() 
{
}

// -------------------------------------------------------------------------

ClassImp(MpdEmcMatch)

// -------------------------------------------------------------------------
//                            MpdRegion source file                    -----
//        Created 23/07/13  by P. Batyuk (MPD) <batyuk@jinr.ru>        -----
//                          from PndRegion (PNDROOT)                   -----
// -------------------------------------------------------------------------

/// Last modified: 24.07.2013, P.B.

#include "MpdRegion.h"

// -------------   Default constructor  ----------------------------------
MpdRegion::MpdRegion( Double_t Zmin, Double_t Zmax)
  :TObject()
{
  fZmin=Zmin;
  fZmax=Zmax;
}

// ------------   Destructor   --------------------------------------------
MpdRegion::~MpdRegion() { }

// ------------   Check if inside this region-------------------------------
Bool_t MpdRegion::IsInside(Double_t Z){
  if( Z>=fZmin && Z<=fZmax ) return kTRUE;
  else return kFALSE;
}

ClassImp(MpdRegion)

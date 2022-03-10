// -------------------------------------------------------------------------
//                            MpdFieldMapData source file              -----
//        Created 23/07/13  by P. Batyuk (MPD) <batyuk@jinr.ru>        -----
//                          from MpdFieldMapData (PNDROOT)             -----
// -------------------------------------------------------------------------

/// Last modified: 24.07.2013, P.B.

#include "TArrayF.h"

#include "MpdFieldMap.h"
#include "MpdFieldMapData.h"

// -------------   Default constructor  ----------------------------------
MpdFieldMapData::MpdFieldMapData() {
  fType  = 1;
  fXmin  = fYmin  = fZmin  = 0.;
  fXmax  = fYmax  = fZmax  = 0.;
  fNx    = fNy    = fNz    = 0;
  fBx    = fBy    = fBz    = NULL;
}

// -------------   Standard constructor   ---------------------------------
MpdFieldMapData::MpdFieldMapData(const char* mapName)
  : TNamed(mapName, "PND Field Map Data") {
  fType  = 1;
  fXmin  = fYmin  = fZmin  = 0.;
  fXmax  = fYmax  = fZmax  = 0.;
  fNx    = fNy    = fNz    = 0;
  fBx    = fBy    = fBz    = NULL;
}

// -----   Constructor from MpdFieldMap   ------------------------------
MpdFieldMapData::MpdFieldMapData(const char* name,
				 const MpdFieldMap& map) 
: TNamed(name, "PND Field Map Data") {

  fType  = map.GetType();
  fXmin  = map.GetXmin();
  fYmin  = map.GetYmin();
  fZmin  = map.GetZmin();
  fXmax  = map.GetXmax();
  fYmax  = map.GetYmax();
  fZmax  = map.GetZmax();
  fNx    = map.GetNx();
  fNy    = map.GetNy();
  fNz    = map.GetNz();
  fBx = new TArrayF(*(map.GetBx()));
  fBy = new TArrayF(*(map.GetBy()));
  fBz = new TArrayF(*(map.GetBz()));
  fUnit = map.GetUnit();
  
  // Take out scaling factor and convert from kG to T
  Double_t factor = map.GetScale() * 10.; 
  Int_t index = 0;
  for (Int_t ix=0; ix<fNx; ix++) {
    for (Int_t iy=0; iy<fNy; iy++) {
      for (Int_t iz=0; iz<fNz; iz++) {
	index = ix*fNy*fNz + iy*fNz + iz;
	if ( fBx ) (*fBx)[index] = (*fBx)[index] / factor;
	if ( fBy ) (*fBy)[index] = (*fBy)[index] / factor;
	if ( fBz ) (*fBz)[index] = (*fBz)[index] / factor;
      }    // z loop
    }      // y loop
  }        // x loop
}

// ------------   Destructor   --------------------------------------------
MpdFieldMapData::~MpdFieldMapData() {
  if ( fBx ) delete fBx;
  if ( fBy ) delete fBy;
  if ( fBz ) delete fBz;
}

ClassImp(MpdFieldMapData)

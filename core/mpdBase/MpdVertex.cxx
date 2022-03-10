// -------------------------------------------------------------------------
// -----                       MpdVertex source file                   -----
// -----                  Created 10/07/09  by A. Zinchenko            -----
// -------------------------------------------------------------------------
#include "MpdVertex.h"

#include <iostream>

using std::cout;
using std::endl;


// -----   Default constructor   -------------------------------------------
MpdVertex::MpdVertex() : TNamed("Vertex", "Global") {
  fX = fY = fZ = fChi2 = 0.;
  fNDF = fNTracks = 0;
  for (Int_t i=0; i<6; i++) fCovMatrix[i] = 0;
  fCovMatrix[0] = fCovMatrix[3] = fCovMatrix[5] = 1.0;
  fTrInd = new TArrayI(1);
}
// -------------------------------------------------------------------------

// -----   Constructor with name and title   -------------------------------
MpdVertex::MpdVertex(const char* name, const char* title) 
  : TNamed(name, title) {
  fX = fY = fZ = fChi2 = 0.;
  fNDF = fNTracks = 0;
  for (Int_t i=0; i<6; i++) fCovMatrix[i] = 0;
  fCovMatrix[0] = fCovMatrix[3] = fCovMatrix[5] = 1.0;
  fTrInd = new TArrayI(1);
}
// -------------------------------------------------------------------------
  
// -----   Constructor with all parameters   -------------------------------
MpdVertex::MpdVertex(const char* name, const char* title,
		     Double_t x, Double_t y, Double_t z, Double_t chi2,
		     Int_t ndf, Int_t nTracks, 
		     const TMatrixFSym& covMat) 
  : TNamed(name, title) {
  fTitle   = title;
  fX       = x;
  fY       = y;
  fZ       = z;
  fChi2    = chi2;
  fNDF     = ndf;
  fNTracks = nTracks;
  Int_t index = 0;
  for (Int_t i=0; i<3; i++) {
    for (Int_t j=i; j<3; j++) fCovMatrix[index++] = covMat[i][j];
  }
  fTrInd = new TArrayI(1);
}
// -------------------------------------------------------------------------
MpdVertex::MpdVertex(const MpdVertex& vert)
  : TNamed(vert),
    fX(vert.fX), fY(vert.fY), fZ(vert.fZ), fChi2(vert.fChi2), 
    fNDF(vert.fNDF), fNTracks(vert.fNTracks), fTrInd(new TArrayI(1))
{
  ///copy constructor

  for (Int_t i = 0; i < 6; ++i) fCovMatrix[i] = vert.fCovMatrix[i];
  *fTrInd = *(vert.fTrInd);
}

//__________________________________________________________________________
MpdVertex & MpdVertex::operator= (const MpdVertex& vert) 
{
  /// Asignment operator

  // check assignement to self
  if (this == &vert) return *this;

  // base class assignement
  TNamed::operator=(vert);

  fX = vert.fX;
  fY = vert.fY;
  fZ = vert.fZ;
  fChi2 = vert.fChi2;
  fNDF = vert.fNDF;
  fNTracks = vert.fNTracks;
  fTrInd = new TArrayI(1);

  for (Int_t i = 0; i < 6; ++i) fCovMatrix[i] = vert.fCovMatrix[i];
  *fTrInd = *(vert.fTrInd);
  return *this;
}

// -----   Destructor   ----------------------------------------------------
MpdVertex::~MpdVertex() { delete fTrInd; }
// -------------------------------------------------------------------------

// -----   Public method Print   -------------------------------------------
void MpdVertex::Print() {
  Double_t chi2ndf;
  if (fNDF) chi2ndf = fChi2 / Double_t(fNDF);
  else chi2ndf = 0.;	       
  cout << "Vertex coord. (" << fX << "," << fY << "," << fZ << ") cm, "
       << "chi2/ndf = " << chi2ndf << ", " << fNTracks
       << " tracks used" << endl;
}
// -------------------------------------------------------------------------

// -----   Accessor to covariance matrix    --------------------------------
void MpdVertex::CovMatrix(TMatrixFSym& covMat) const {
  Int_t index = 0;
  for (int i=0; i<3; i++) {
    for (int j=i; j<3; j++) {
      covMat[i][j] = fCovMatrix[index];
      covMat[j][i] = fCovMatrix[index];
      index++;
    }
  }
}
// -------------------------------------------------------------------------

// -----   Accessor to covariance matrix elements   ------------------------
Double_t MpdVertex::GetCovariance(Int_t i, Int_t j) const {
  TMatrixFSym* mat = new TMatrixFSym(3);
  CovMatrix(*mat);
  Double_t element = (*mat)[i][j];
  delete mat;
  return element;
}
// -------------------------------------------------------------------------

// -----   Public method SetVertex   ---------------------------------------
void MpdVertex::SetVertex(Double_t x, Double_t y, Double_t z, Double_t chi2,
			  Int_t ndf, Int_t nTracks, 
			  const TMatrixFSym& covMat) {
  fX       = x;
  fY       = y;
  fZ       = z;
  fChi2    = chi2;
  fNDF     = ndf;
  fNTracks = nTracks;
  Int_t index = 0;
  for (Int_t i=0; i<3; i++) {
    for (Int_t j=i; j<3; j++) fCovMatrix[index++] = covMat[i][j];
  }
}
// -------------------------------------------------------------------------

// -----   Public method Reset   -------------------------------------------
void MpdVertex::Reset() {
  fX = fY = fZ = fChi2 = 0.;
  fNDF = fNTracks = 0;
  for(Int_t i=0; i<6; i++) fCovMatrix[i] = 0;
}  
// -------------------------------------------------------------------------

ClassImp(MpdVertex)

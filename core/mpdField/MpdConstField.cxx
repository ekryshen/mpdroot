// -------------------------------------------------------------------------
//                            MpdContField source file                 -----
//        Created 23/07/13  by P. Batyuk (MPD) <batyuk@jinr.ru>        -----
//                          from MpdContField (PNDROOT)                -----
// -------------------------------------------------------------------------

/// Last modified: 24.07.2013, P.B.

#include <iomanip>
#include <iostream>

#include "MpdConstField.h"
#include "MpdConstPar.h"

using namespace std;

// -----   Default constructor   -------------------------------------------
MpdConstField::MpdConstField() {
  fXmin = fXmax = fYmin = fYmax = fZmin = fZmax = fBx = fBy = fBz = 0.;
  fType = 0;
}

// -----   Standard constructor   ------------------------------------------
MpdConstField::MpdConstField(const char* name, Double_t xMin, 
			     Double_t xMax, Double_t yMin, 
			     Double_t yMax, Double_t zMin,
			     Double_t zMax, Double_t bX, 
			     Double_t bY, Double_t bZ) 
  : FairField(name) {
  fXmin = xMin;
  fXmax = xMax;
  fYmin = yMin;
  fYmax = yMax;
  fZmin = zMin;
  fZmax = zMax;
  fBx   = bX;
  fBy   = bY;
  fBz   = bZ;
  fType = 0;
}

// --------   Constructor from MpdFieldPar   -------------------------------
MpdConstField::MpdConstField(MpdConstPar* fieldPar) 
  :FairField()
{
  if ( ! fieldPar ) {
    cerr << "-W- MpdConstField::MpdConstField: empty parameter container!"
	 << endl;
    fXmin = fXmax = fYmin = fYmax = fZmin = fZmax = fBx = fBy = fBz = 0.;
    fType= -1;
  }
  else {
    fXmin = fieldPar->GetXmin();
    fXmax = fieldPar->GetXmax();
    fYmin = fieldPar->GetYmin();
    fYmax = fieldPar->GetYmax();
    fZmin = fieldPar->GetZmin();
    fZmax = fieldPar->GetZmax();
    fBx   = fieldPar->GetBx();
    fBy   = fieldPar->GetBy();
    fBz   = fieldPar->GetBz();
    fType = fieldPar->GetType();
  }
}

// -----   Destructor   ----------------------------------------------------
MpdConstField::~MpdConstField() { }

// -----   Set field region   ----------------------------------------------
void MpdConstField::SetFieldRegion(Double_t xMin, Double_t xMax, 
				   Double_t yMin, Double_t yMax, 
				   Double_t zMin, Double_t zMax) {
  fXmin = xMin;
  fXmax = xMax;
  fYmin = yMin;
  fYmax = yMax;
  fZmin = zMin;
  fZmax = zMax;
}

// -----   Set field values   ----------------------------------------------
void MpdConstField::SetField(Double_t bX, Double_t bY, Double_t bZ) {
  fBx   = bX;
  fBy   = bY;
  fBz   = bZ;
}

// -----   Get x component of field   --------------------------------------
Double_t MpdConstField::GetBx(Double_t x, Double_t y, Double_t z) {
  if ( x <= fXmin  ||  x >= fXmax  ||
       y <= fYmin  ||  y >= fYmax  ||
       z <= fZmin  ||  z >= fZmax ) return 0.;
  return fBx;
}

// -----   Get y component of field   --------------------------------------
Double_t MpdConstField::GetBy(Double_t x, Double_t y, Double_t z) {
  if ( x <= fXmin  ||  x >= fXmax  ||
       y <= fYmin  ||  y >= fYmax  ||
       z <= fZmin  ||  z >= fZmax ) return 0.;
  return fBy;
}

// -----   Get z component of field   --------------------------------------
Double_t MpdConstField::GetBz(Double_t x, Double_t y, Double_t z) {
  if ( x <= fXmin  ||  x >= fXmax  ||
       y <= fYmin  ||  y >= fYmax  ||
       z <= fZmin  ||  z >= fZmax ) return 0.;
  return fBz;
}

// -------------------------------------------------------------------------
void MpdConstField::GetFieldValue(const Double_t point[3], Double_t* bField)
{
  bField[0]=GetBx(point[0],point[1],point[2]);
  bField[1]=GetBy(point[0],point[1],point[2]);
  bField[2]=GetBz(point[0],point[1],point[2]);
}

// -------------------------------------------------------------------------
void MpdConstField::GetBxyz(const Double_t point[3], Double_t* bField)
{
  GetFieldValue(point,bField);
}

// -----   Screen output   -------------------------------------------------
void MpdConstField::Print() {
  cout << "======================================================" << endl;
  cout << "----  " << fTitle << " : " << fName << endl;
  cout << "----" << endl;
  cout << "----  Field type    : constant" << endl;
  cout << "----" << endl;
  cout << "----  Field regions : " << endl;
  cout << "----        x = " << setw(4) << fXmin << " to " << setw(4) 
       << fXmax << " cm" << endl;
  cout << "----        y = " << setw(4) << fYmin << " to " << setw(4) 
       << fYmax << " cm" << endl;
  cout << "----        z = " << setw(4) << fZmin << " to " << setw(4)
       << fZmax << " cm" << endl;
  cout.precision(4);
  cout << "----  B = ( " << fBx << ", " << fBy << ", " << fBz << " ) kG"
       << endl;
  cout << "======================================================" << endl;
}

ClassImp(MpdConstField)

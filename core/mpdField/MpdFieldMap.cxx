// -------------------------------------------------------------------------
//                            MpdFieldMap source file                  -----
//        Created 23/07/13  by P. Batyuk (MPD) <batyuk@jinr.ru>        -----
//                          from PndFieldMap (PNDROOT)                 -----
// -------------------------------------------------------------------------

/// Last modified: 24.07.2013, P.B.
/// Last modified: 30.07.2013, P.B.
/// Last modified: 31.07.2013, P.B.

#include <iomanip>
#include <iostream>
#include <fstream>
#include "stdlib.h"

#include "TArrayF.h"
#include "TFile.h"
#include "TMath.h"

#include "MpdFieldMap.h"
#include "MpdFieldMapData.h"
#include "MpdFieldPar.h"

using namespace std;

// -------------   Default constructor  ----------------------------------
MpdFieldMap::MpdFieldMap() {
  fPosX  = fPosY  = fPosZ  = 0.;
  fXmin  = fYmin  = fZmin  = 0.;
  fXmax  = fYmax  = fZmax  = 0.;
  fXstep = fYstep = fZstep = 0.;
  fNx    = fNy    = fNz    = 0;
  fScale = 1.;
  funit = 10.0;
  fBx    = fBy    = fBz    = NULL;
  fPosX = fPosY = fPosZ = 0.;
  fName     = "";
  fFileName = "";
  fType = 1;
}

// -------------   Standard constructor   ---------------------------------
MpdFieldMap::MpdFieldMap(const char* mapName, const char* fileType)
  : FairField(mapName) {
  fPosX  = fPosY  = fPosZ  = 0.;
  fXmin  = fYmin  = fZmin  = 0.;
  fXmax  = fYmax  = fZmax  = 0.;
  fXstep = fYstep = fZstep = 0.;
  fNx    = fNy    = fNz    = 0;
  fScale = 1.;
  funit = 10.0;
  fBx    = fBy    = fBz    = NULL;
  fName  = mapName;
  TString dir = getenv("VMCWORKDIR");
  fFileName = dir + "/input/" + mapName;
  if ( fileType[0] == 'R' ) fFileName += ".root";
  else                      fFileName += ".dat";
  fType = 1;
}

// ------------   Constructor from MpdFieldPar   --------------------------
MpdFieldMap::MpdFieldMap(MpdFieldPar* fieldPar) {
  fType = 1;
  fPosX  = fPosY  = fPosZ  = 0.;
  fXmin  = fYmin  = fZmin  = 0.;
  fXmax  = fYmax  = fZmax  = 0.;
  fXstep = fYstep = fZstep = 0.;
  fNx    = fNy    = fNz    = 0;
  fScale = 1.;
  funit = 10.0;
  fBx    = fBy    = fBz    = NULL;
  if ( ! fieldPar ) {
    cerr << "-W- MpdConstField::MpdConstField: empty parameter container!"
	 << endl;
    fName     = "";
    fFileName = "";
    fType     = -1;
  }
  else {
    fieldPar->MapName(fName);
    fPosX  = fieldPar->GetPositionX();
    fPosY  = fieldPar->GetPositionY();
    fPosZ  = fieldPar->GetPositionZ();
    fScale = fieldPar->GetScale();
    TString dir = getenv("VMCWORKDIR");
    fFileName = dir + "/input/" + fName + ".root";
    fType = fieldPar->GetType();
  }
}

// ------------   Destructor   --------------------------------------------
MpdFieldMap::~MpdFieldMap() {
  if ( fBx ) delete fBx;
  if ( fBy ) delete fBy;
  if ( fBz ) delete fBz;
}

// -----------   Intialization   ------------------------------------------
void MpdFieldMap::Init() {
  if      (fFileName.EndsWith(".root")) ReadRootFile(fFileName, fName);
  else if (fFileName.EndsWith(".dat"))  ReadAsciiFile(fFileName);
  else {
    cerr << "-E- MpdFieldMap::Init: No proper file name defined! ("
	 << fFileName << ")" << endl;
    Fatal("Init", "No proper file name");
  }
}

// -----------   Check whether a point is inside the map   ----------------
Bool_t MpdFieldMap::IsInside(Double_t x, Double_t y, Double_t z) {
  
  /// Define geometry sizes of TPC to interpolate mag. field....
  /// Extrapolation used to define mag. field components in vicinity of the TPC-volume 
  /// r1extrap <= R <= r2extrap
  /// r1 - 2cm <= R <= r2 + 2 cm

  Double_t delta_r(2.);

  Double_t r1(40.3), r2(120.3); /// Real internal and external radiuses of the TPC

  //Double_t fZmin(-170.), fZmax(170.);

  //Double_t fXmin_r2 = -r2, fXmax_r2 = r2;
  //Double_t fYmin_r2 = fXmin_r2, fYmax_r2 = fXmax_r2;

  //Double_t fXmin_r1 = -r1, fXmax_r1 = r1;
  //Double_t fYmin_r1 = fXmin_r1, fYmax_r1 = fXmax_r1;
  
  /// Checking only correct field map initialization
  /// In fact, the mag. field used at the point (0, 0, 0) is considered to be (0, 0, 5) kG 
  if (x == 0 && y == 0 && z== 0) return kTRUE;
 
  /// Checking a point to be inside the TPC (B-field Map) (r1 - delta_r = 40.3 cm - 2 cm, r2 + delta_r = 120.3 cm + 2 cm)
  else if ( TMath::Power(x, 2) + TMath::Power(y, 2) >= TMath::Power(r1 - delta_r, 2) && 
	TMath::Power(x, 2) + TMath::Power(y, 2) <= TMath::Power(r2 + delta_r, 2)) return kTRUE;
  
  else return kFALSE;
}

// ----------   Write the map to an ASCII file   --------------------------
void MpdFieldMap::WriteAsciiFile(const char* fileName) {
  
  // Open file
  cout << "-I- MpdFieldMap: Writing field map to ASCII file " 
       << fileName << endl;
  ofstream mapFile(fileName);
  if ( ! mapFile.is_open() ) {
    cerr << "-E- MpdFieldMap:ReadAsciiFile: Could not open file! " << endl;
    return;
  }
  
  // Write field map grid parameters
  mapFile.precision(4);
  mapFile << showpoint;
  if ( fType == 1 ) mapFile << "nosym" << endl;
  if ( funit == 10.0  ) mapFile << "T" << endl;
  else if ( funit == 0.001 ) mapFile << "G" << endl;
  else if ( funit == 1.0   ) mapFile << "kG" << endl;
  
  mapFile << fXmin << " " << fXmax << " " << fNx << endl;
  mapFile << fYmin << " " << fYmax << " " << fNy << endl;
  mapFile << fZmin << " " << fZmax << " " << fNz << endl;

  // Write field values
  Double_t factor = funit * fScale;  // Takes out scaling 
  cout << right;
  Int_t nTot = fNx * fNy * fNz;
  cout << "-I- MpdFieldMap: " << fNx*fNy*fNz << " entries to write... " 
       << setw(3) << 0 << " % ";
  Int_t index=0;
  div_t modul;
  Int_t iDiv = TMath::Nint(nTot/100.) + 1;
  for(Int_t ix=0; ix<fNx; ix++) {
    for(Int_t iy=0; iy<fNy; iy++) {
      for(Int_t iz=0; iz<fNz; iz++) {
	index =ix*fNy*fNz + iy*fNz + iz;
	modul = div(index,iDiv);
	if ( modul.rem == 0 ) {
	  Double_t perc = TMath::Nint(100.*index/nTot);
	  cout << "\b\b\b\b\b\b" << setw(3) << perc << " % " << flush;
	}
	mapFile << fBx->At(index)/factor << " " << fBy->At(index)/factor 
		<< " " << fBz->At(index)/factor << endl;
      } // z-Loop
    }   // y-Loop
  }     // x-Loop
  cout << "   " << index+1 << " written" << endl;
  mapFile.close();		
}	

// -------   Write field map to a ROOT file   -----------------------------
void MpdFieldMap::WriteRootFile(const char* fileName,
				const char* mapName) {
  
  MpdFieldMapData* data = new MpdFieldMapData(mapName, *this);
  TFile* oldFile = gFile;
  TFile* file = new TFile(fileName, "RECREATE");
  data->Write();
  file->Close();
  if(oldFile) oldFile->cd();
}

// -----  Set the position of the field centre in global coordinates  -----
void MpdFieldMap::SetPosition(Double_t x, Double_t y, Double_t z) {
  fPosX = x;
  fPosY = y;
  fPosZ = z;
}

// ---------   Screen output   --------------------------------------------
void MpdFieldMap::Print() {
  TString type = "Map";
  if (fType == 1 ) type = "Nosym Map";
  else if ( fType == 2 ) type = "Solenoid Map ";
  else if ( fType == 3 ) type = "Dipole Map ";
  else if ( fType == 4 ) type = "Trans Map ";
  
  cout << "======================================================" << endl;
  
  cout.precision(4);
  cout << showpoint;
  cout << "----  " << fTitle << " : " << fName << endl;
  cout << "----" << endl;
  cout << "----  Field type     : " << type << endl;
  cout << "----" << endl;
  cout << "----  Field map grid : " << endl;
  cout << "----  x = " << setw(4) << fXmin << " to " << setw(4) << fXmax 
       << " cm, " << fNx << " grid points, dx = " << fXstep << " cm" << endl;
  cout << "----  y = " << setw(4) << fYmin << " to " << setw(4) << fYmax 
       << " cm, " << fNy << " grid points, dy = " << fYstep << " cm" << endl;
  cout << "----  z = " << setw(4) << fZmin << " to " << setw(4) << fZmax 
       << " cm, " << fNz << " grid points, dz = " << fZstep << " cm" << endl;
  cout << endl;
  cout << "----  Field centre position: ( " << setw(6) << fPosX << ", "
       << setw(6) << fPosY << ", " << setw(6) << fPosZ << ") cm" << endl;
  cout << "----  Field scaling factor: " << fScale << endl;
  
  Double_t bx = GetBx(0.,0.,0.);
  Double_t by = GetBy(0.,0.,0.);
  Double_t bz = GetBz(0.,0.,0.);
  
  cout << "----" << endl;
  cout << "----  Field at origin is ( " << setw(6) << bx << ", " << setw(6)
       << by << ", " << setw(6) << bz << ") kG" << endl;
  cout << "======================================================" << endl;
}

// ---------    Reset parameters and data (private)  ----------------------
void MpdFieldMap::Reset() {
  fPosX = fPosY = fPosZ = 0.;
  fXmin = fYmin = fZmin = 0.;
  fXmax = fYmax = fZmax = 0.;
  fXstep = fYstep = fZstep = 0.;
  fNx = fNy = fNz = 0;
  fScale = 1.;
  funit = 10.0;
  if ( fBx ) { delete fBx; fBx = NULL; }
  if ( fBy ) { delete fBy; fBy = NULL; }
  if ( fBz ) { delete fBz; fBz = NULL; }
}

// -----   Read field map from ASCII file (private)   ---------------------
void MpdFieldMap::ReadAsciiFile(const char* fileName) {

  Double_t xx=0., yy=0., zz=0.;
  Double_t bx=0., by=0., bz=0.;

  // Open file
  cout << "-I- MpdFieldMap: Reading field map from ASCII file " 
       << fileName << endl;
  ifstream mapFile(fileName);
  if ( ! mapFile.is_open() ) {
    cerr << "-E- MpdFieldMap:ReadAsciiFile: Could not open file! " << endl;
    Fatal("ReadAsciiFile","Could not open file");
  }

  // Read map type
  TString type;
  mapFile >> type;

  Int_t iType = 0;
  if ( type == "nosym" ) iType = 1;
  else if ( type == "Solenoid") iType = 2;
  else if ( type == "Dipole"  ) iType = 3;
  else if ( type == "Trans"  ) iType = 4;
  if ( fType != iType ) {
    cout << "-E- MpdFieldMap::ReadAsciiFile: Incompatible map types!"
	 << endl;
    cout << "    Field map is of type " << fType 
	 << " but map on file is of type " << iType << endl;
    Fatal("ReadAsciiFile","Incompatible map types");
  }
  // Read Units
  TString unit;
  mapFile >> unit;
  if ( unit == "G" ) funit = 0.001;
  else if ( unit == "T"  ) funit = 10.0;
  else if ( unit == "kG"  ) funit=1.0;
  else {
    cout << "-E- FieldMap::ReadAsciiFile: No units!"
	 << endl;
    Fatal("ReadAsciiFile","No units defined");
  }
   
  // Read grid parameters
  
  mapFile >>fXmin >> fXmax >> fNx;
  mapFile >>fYmin >> fYmax >> fNy;
  mapFile >>fZmin >> fZmax >> fNz;
  fXstep = ( fXmax - fXmin ) / Double_t( fNx - 1 );
  fYstep = ( fYmax - fYmin ) / Double_t( fNy - 1 );
  fZstep = ( fZmax - fZmin ) / Double_t( fNz - 1 );
  
  // Create field arrays
  fBx = new TArrayF(fNx * fNy * fNz);
  fBy = new TArrayF(fNx * fNy * fNz);
  fBz = new TArrayF(fNx * fNy * fNz);
  
  // Read the field values
  Double_t factor = fScale * funit;   // Factor 1/1000 for G -> kG
  cout << right;
  Int_t nTot = fNx * fNy * fNz;
  cout << "-I- MpdFieldMap: " << nTot << " entries to read... " 
       << setw(3) << 0 << " % ";
  Int_t index = 0;
  div_t modul;
  Int_t iDiv = TMath::Nint(nTot/100.) + 1;
  for (Int_t ix=0; ix<fNx; ix++) {
    for (Int_t iy = 0; iy<fNy; iy++) {
      for (Int_t iz = 0; iz<fNz; iz++) {
	if (! mapFile.good()) cerr << "-E- MpdFieldMap::ReadAsciiFile: "
				   << "I/O Error at " << ix << " "
				   << iy << " " << iz << endl;
	index = ix*fNy*fNz + iy*fNz + iz;
	modul = div(index,iDiv);
	if ( modul.rem == 0 ) {
	  Double_t perc = TMath::Nint(100.*index/nTot);
	  cout << "\b\b\b\b\b\b" << setw(3) << perc << " % " << flush;
	}
	//mapFile >>  bx >> by >> bz;
	mapFile >> xx >> yy >> zz >> bx >> by >> bz ;
	//cout  << " x= " <<xx <<" y= " << yy<<" z= " << zz<<" bx= " <<  bx <<" by= " <<by <<" bz= " << bz<< endl;
	fBx->AddAt(factor*bx, index);
	fBy->AddAt(factor*by, index);
	fBz->AddAt(factor*bz, index);
	if ( mapFile.eof() ) {
	  cerr << endl << "-E- MpdFieldMap::ReadAsciiFile: EOF"
	       << " reached at " << ix << " " << iy << " " << iz << endl;
	  mapFile.close();
	  break;
	}
      }   // z-Loop
    }     // y-Loop0)
  }       // x-Loop
  
  cout << "   " << index+1 << " read" << endl;
  
  mapFile.close();
}

// -------------   Read field map from ROOT file (private)  ---------------
void MpdFieldMap::ReadRootFile(const char* fileName, 
			       const char* mapName) {

   cout<<"Field Map is read from ASCII file ..."<<endl;
}

// ------------   Set field parameters and data (private)  ----------------
void MpdFieldMap::SetField(const MpdFieldMapData* data) {

  // Check compatibility
  if ( data->GetType() != fType ) {
    cout << "-E- MpdFieldMap::SetField: Incompatible map types!"
	 << endl;
    cout << "    Field map is of type " << fType 
	 << " but map on file is of type " << data->GetType() << endl;
    Fatal("SetField","Incompatible map types");
  }
  
  fXmin = data->GetXmin();
  fYmin = data->GetYmin();
  fZmin = data->GetZmin();
  fXmax = data->GetXmax();
  fYmax = data->GetYmax();
  fZmax = data->GetZmax();
  fNx = data->GetNx();
  fNy = data->GetNy();
  fNz = data->GetNz();
  fXstep = ( fXmax - fXmin ) / Double_t( fNx - 1 );
  fYstep = ( fYmax - fYmin ) / Double_t( fNy - 1 );
  fZstep = ( fZmax - fZmin ) / Double_t( fNz - 1 );
  if ( fBx ) delete fBx;
  if ( fBy ) delete fBy;
  if ( fBz ) delete fBz;
  fBx = new TArrayF(*(data->GetBx()));
  fBy = new TArrayF(*(data->GetBy()));
  fBz = new TArrayF(*(data->GetBz()));
  
  // Scale and convert from G(or T) to kG
  Double_t factor = fScale * funit;
  Int_t index = 0;
  for (Int_t ix=0; ix<fNx; ix++) {
    for (Int_t iy=0; iy<fNy; iy++) {
      for (Int_t iz=0; iz<fNz; iz++) {
	index = ix*fNy*fNz + iy*fNz + iz;
	if ( fBx ) (*fBx)[index] = (*fBx)[index] * factor;
	if ( fBy ) (*fBy)[index] = (*fBy)[index] * factor;
	if ( fBz ) (*fBz)[index] = (*fBz)[index] * factor;
      }
    }
  }
}

// ------------------------------------------------------------------------  
Double_t MpdFieldMap::MPDCalcPhi(Double_t x, Double_t y) {
  
  /// We assume that Phi is located in the interval [0; 2pi), P.B.
  
  Double_t phi(0.);
  
  if (x > 0. && y>=0.) phi = TMath::ATan(y/x);
  else if (x > 0. && y < 0.) phi = TMath::ATan(y/x) + 2*TMath::Pi();
  else if (x < 0.) phi = TMath::ATan(y/x) + TMath::Pi();
  else if (x == 0. && y > 0.) phi = 0.5*TMath::Pi();
  else if (x == 0. && y < 0.) phi = 1.5*TMath::Pi();
  else if (x == 0. && y == 0.) phi = 0.;
  
  return phi;
}

//-----------------------------------------------------------------------------------
Double_t MpdFieldMap::MPDCalcBr(Double_t x, Double_t y) {
  
  /// Mag. Field in T!!!
  
  Double_t br_coeff[25] = {3.55683796461738e-005, -1.45190262790289e-005, -4.22589216944936e-005, -3.7825589770834e-006,
			   5.31016888137445e-006, -0.000357993220605977, 0.000104590154606646, 0.000421949428873677,
			   1.14589001127912e-005, 7.04773281719065e-007, 0.000289792998812542, -0.000269810303185643,
			   -0.000367922801279625, 8.60754473338966e-005, 0.000114899755175272, -0.000301080291666739, 
			   0.000447886344136466, 0.000250622861695885, -0.000333504572726095, -0.000111012107349364, 
			   7.82187500000153e-005, -0.00015617776565864, -8.1476324428263e-005, 0.000159232809881339, 
			   2.93802936497351e-005};
  
  Double_t sum(0.);
  
  for (int i(0); i<5; i++) sum += br_coeff[i]*pow(y,i);
  for (int i(5); i<10; i++) sum += br_coeff[i]*x*pow(y,i-5);
  for (int i(10); i<15; i++) sum += br_coeff[i]*pow(x,2)*pow(y,i-10);
  for (int i(15); i<20; i++) sum += br_coeff[i]*pow(x,3)*pow(y,i-15);
  for (int i(20); i<25; i++) sum += br_coeff[i]*pow(x,4)*pow(y,i-20);
  
  return sum;
}

//-----------------------------------------------------------------------------------
Double_t MpdFieldMap::MPDCalcBz(Double_t x, Double_t y) {
  
  /// Mag. Field in T!!!
  
  Double_t bz_coeff[25] = {0.500047838143456, 0.00044540880088384, -5.35189259780466e-005, -0.000247723183378179,
			   2.63504568969453e-005, 1.23804543665921e-005, 0.000177653223608631, 1.66011921364229e-005, 
			   6.08797488443891e-005, -8.64585638589155e-005, -3.187574503429e-005, -0.000208086128594243,
			   -9.1879573536513e-005, -1.4640341932548e-005, 0.000284560265815692, 7.60046874443087e-005,
			   0.000394401600524841, -0.000120094441661944, 0.000101870491589828, -0.000203451132764876, 
			   -1.8619791659889e-005, -0.000151059550349997, 8.41518325664925e-005, -3.90510111913933e-005,
			   2.77271406987634e-005};
    
  Double_t sum(0.);
  
  for (int i(0); i<5; i++) sum += bz_coeff[i]*pow(y,i);
  for (int i(5); i<10; i++) sum += bz_coeff[i]*x*pow(y,i-5);
  for (int i(10); i<15; i++) sum += bz_coeff[i]*pow(x,2)*pow(y,i-10);
  for (int i(15); i<20; i++) sum += bz_coeff[i]*pow(x,3)*pow(y,i-15);
  for (int i(20); i<25; i++) sum += bz_coeff[i]*pow(x,4)*pow(y,i-20);
  
  return  sum; 
}

/// MPDGetBx(y, z) return Mag. Field in kG!

//-------------------------------------------------------------------------------------
Double_t MpdFieldMap::GetBx(Double_t x, Double_t y, Double_t z) {
  
  if (IsInside(x, y, z) == kTRUE) {
    
    x /= 100; y /= 100; z /=100;
    
    Double_t phi = MPDCalcPhi(x, y);
    
    Double_t r = TMath::Sqrt(x*x + y*y);
    
    Double_t Br = MPDCalcBr(r,z);
    
    return Br*TMath::Cos(phi)*fScale*funit;
  }
  
  else return 0;
}

//--------------------------------------------------------------------------------------
Double_t MpdFieldMap::GetBy(Double_t x, Double_t y, Double_t z) {
  
  if (IsInside(x, y, z) == kTRUE) {
    
    x /= 100; y /= 100; z /=100;
    
    Double_t phi = MPDCalcPhi(x, y);
    
    Double_t r = TMath::Sqrt(x*x + y*y);
    
    Double_t Br = MPDCalcBr(r,z);
    
    return Br*TMath::Sin(phi)*fScale*funit;
  }
  
  else return 0;
}

//--------------------------------------------------------------------------------------
Double_t MpdFieldMap::GetBz(Double_t x, Double_t y, Double_t z) {
  
  if (IsInside(x, y, z) == kTRUE) {
    
    x /= 100; y /= 100; z /= 100;
    
    Double_t phi = MPDCalcPhi(x, y);
    
    Double_t r = TMath::Sqrt(x*x + y*y);
    
    return MPDCalcBz(r,z)*fScale*funit;
  }
  
  else return 5.;
}

//---------------------------------------------------------------------------------------
void MpdFieldMap::GetBxyz(const Double_t point[3], Double_t* bField) {
  
  bField[0] = GetBx(point[0], point[1], point[2]);
  bField[1] = GetBy(point[0], point[1], point[2]);
  bField[2] = GetBz(point[0], point[1], point[2]);
  
  //cout<<"x = " <<point[0] <<" y = " <<point[1] <<" z = " <<point[2];
  //cout<<" Bx = "<<bField[0]<<" By = "<<bField[1]<<" Bz = "<<bField[2]<<endl;
}

//---------------------------------------------------------------------------------------
void MpdFieldMap::GetFieldValue(const Double_t point[3], Double_t* bField) {
  
  bField[0] = GetBx(point[0], point[1], point[2]);
  bField[1] = GetBy(point[0], point[1], point[2]);
  bField[2] = GetBz(point[0], point[1], point[2]);
}

Double_t MpdFieldMap::Interpolate(Double_t dx, Double_t dy, Double_t dz)
{
    cout<<"MpdFieldMap::Interpolate( is not implemented!!!"<<endl;
    return 0;
}

ClassImp(MpdFieldMap)

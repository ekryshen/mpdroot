
//************************************************
//  Script to prepare geometric file for ITS
//************************************************ 
class FairGeoRotation;

#include <TVector3.h>

#include <iostream>
#include <fstream>
#include <string>

class Mpdshape: public TObject {

public:
  Mpdshape();
  Mpdshape(ofstream* file, std::string  volume_name,
	   std::string  mother_volume_name, std::string shape,
	   std::string media, std::string  points,
	   std::string position = "0. 0. 0.", 
	   std::string rotation = "1. 0. 0. 0. 1. 0. 0. 0. 1.", 
	   Int_t segment = 0, Int_t mother_segment = 0); 
  Mpdshape(ofstream* file, std::string  volume_name,
	   std::string  mother_volume_name, std::string shape,
	   std::string media); 
  ~Mpdshape();
  
  void SetFile(ofstream* file) { fFile = file; }
  void SetVolumeName(std::string  name) { fVolumeName = name; }
  void SetMotherVolumeName(std::string  name) { fMotherVolumeName = name; }
  void SetShape(std::string  shape) { fShape = shape; }
  void SetMedia(std::string  media) { fMedia = media; }
  void SetPoints(std::string  points) { fPoints = points; }

  void SetPosition(std::string  position) { fPosition = position; }
  void SetPosition(Double_t x, Double_t y, Double_t z, Int_t precision=6);
  void SetPosition(TVector3 vec);

  void SetRotation(std::string  rotation) { fRotation = rotation; }
  void SetRotation(Double_t p1, Double_t p2, Double_t p3,
		   Double_t p4, Double_t p5, Double_t p6,
		   Double_t p7, Double_t p8, Double_t p9);
  void SetRotation(Double_t z, Double_t y1, Double_t z1);

  void SetPositionRotation(std::string  position, std::string  rotation)
  { fPosition = position; fRotation = rotation; }
  void SetSegment(Int_t segment) { fSegment = segment; }
  void SetMotherSegment(Int_t segment) { fMotherSegment = segment; }

  void DumpToFile();
  void DumpWithIncrement() { DumpToFile(); fSegment++; }
  static std::string itoa(Double_t x);

  std::string GetVolumeName() const { return fVolumeName;}
  Int_t GetSegment() const { return fSegment;}
  void Clear();
  void Fill_TUBE(Double_t length_z, Double_t r_max, Double_t r_min);
  void Fill_TUBS(Double_t zmin, Double_t zmax, Double_t rmin,Double_t rmax, Double_t dr_st, Double_t dr_end);
  void Fill_TRAP(Double_t dx11, Double_t dx12, Double_t dy1, Double_t dx21, Double_t dx22, 
		 Double_t dy2, Double_t dz);
  void Fill_TRAP(Double_t x, Double_t X,Double_t x_small_f, Double_t x_large_f, Double_t yW, Double_t yW2,Double_t zW, Double_t zW2, Int_t ra);
  //void Fill_RECTRAP(Double_t x_small, Double_t x_large,Double_t xx_small, Double_t xx_large,
	//	 Double_t yWidth, Double_t zWidth);	 
  void Fill_SPHE( Double_t rzmin, Double_t rmax, Double_t thetamin, Double_t thetamax, Double_t phimin, Double_t phimax);
  void Fill_PGON( Double_t zmin, Double_t zmax, Double_t rmin, Double_t rmax, Double_t phi1);
  //void Fill_BOX(Double_t x, Double_t y, Double_t z);
  void Fill_BRIK(Double_t xWidth1, Double_t xWidth2, Double_t yWidth, Double_t zWidth1, Double_t zWidth2, Int_t ra);
  
  
protected:
  std::ofstream* fFile;
  std::string  fVolumeName;
  std::string  fMotherVolumeName;
  std::string  fShape;
  std::string  fMedia;
  std::string  fPoints;
  std::string  fPosition;
  std::string  fRotation;
  Int_t fSegment;
  Int_t fMotherSegment;
  FairGeoRotation* rot;
    
  ClassDef(Mpdshape,0);
};

ClassImp(Mpdshape);

#include "FairGeoRotation.h"

#include <TROOT.h>
#include <TSystem.h>

#include <sstream>
using namespace std;

Mpdshape::Mpdshape() {

  gROOT->LoadMacro("$VMCWORKDIR/gconfig/basiclibs.C");
  gSystem->Load("libGeoBase");
  gSystem->Load("libParBase");
  gSystem->Load("libBase");
  rot = new FairGeoRotation();
}

Mpdshape::Mpdshape(ofstream* file, std::string  volume_name,
		   std::string  mother_volume_name, std::string  shape,
		   std::string  media, std::string  points, std::string  position,
		   std::string  rotation, Int_t segment, Int_t mother_segment) 
  : fFile(file), fVolumeName(volume_name), fMotherVolumeName(mother_volume_name),
    fShape(shape), fMedia(media), fPoints(points), fPosition(position),
    fRotation(rotation), fSegment(segment), fMotherSegment(mother_segment) {

    gSystem->Load("libGeoBase");
    rot = new FairGeoRotation();
}

Mpdshape::Mpdshape(ofstream* file, std::string  volume_name,
		   std::string  mother_volume_name, std::string  shape,
		   std::string  media):
  fFile(file),  fVolumeName(volume_name), fMotherVolumeName(mother_volume_name),
  fShape(shape), fMedia(media) {

  gSystem->Load("libGeoBase");
  rot = new FairGeoRotation();
  fSegment = 0;
  fMotherSegment = 0;
}

Mpdshape::~Mpdshape() {}

void Mpdshape::DumpToFile() {
  *fFile << "//*************************************************" << endl;
  if (fSegment > 0) {
    *fFile << fVolumeName << "#" << fSegment << endl;;
  } else {
    *fFile << fVolumeName << endl;
  }
  if (fMotherSegment > 0) {
    *fFile << fMotherVolumeName << "#" << fMotherSegment << endl;
  } else {
    *fFile << fMotherVolumeName << endl;
  }

  if (fSegment <= 1) {
    *fFile << fShape << endl << fMedia << endl << fPoints << endl;
  }
  *fFile << fPosition << endl << fRotation << endl;
}

std::string Mpdshape::itoa(Double_t x) {
  std::ostringstream o;
  if (!(o << x)) return "ERROR";
  return o.str();
} 

void Mpdshape::SetPosition(TVector3 vec) {
  ostringstream o;
  o.precision(6);
  o.setf(ios::showpoint); 
  o.setf(ios::fixed);
  o << vec.x() << " " << vec.y() << " " << vec.z();
  fPosition.append(o.str().c_str());
}

void Mpdshape::SetPosition(Double_t x, Double_t y, Double_t z, Int_t precision) {
  ostringstream o;
  o.precision(precision);
  o.setf(ios::showpoint); 
  o.setf(ios::fixed);
  o << x << " " << y << " " << z;
  fPosition = o.str();
}

void Mpdshape::SetRotation(Double_t p1, Double_t p2, Double_t p3,
		      Double_t p4, Double_t p5, Double_t p6,
		      Double_t p7, Double_t p8, Double_t p9) {
	
  std::ostringstream o;
  o.setf(ios::showpoint); 
  o.setf(ios::fixed);
  o.precision(10);
  o << p1 << " " << p2 << " " << p3 << " "
    << p4 << " " << p5 << " " << p6 << " "
    << p7 << " " << p8 << " " << p9;
  fRotation = o.str();
}

void Mpdshape::SetRotation(Double_t z, Double_t y1, Double_t z1) {
  // z = phi
  // y1 = theta
  rot->setEulerAngles(z, y1, z1);
  std::ostringstream o;
  o.setf(ios::showpoint); 
  o.setf(ios::fixed);
  o.precision(10);
  for (Int_t i = 0; i < 9; i++) {
    o << rot->operator()(i) << " ";
  }
  fRotation = o.str();
}

std::string itoa(Double_t x) {
  std::ostringstream o;
  if (!(o << x)) return "ERROR";
  return o.str();
} 

 void its_geometry_bp38() {
 cout << "Start " << endl;

 const Int_t ICables = 0;          // Cables 0 = Off, 1 = On
  Double_t radiusPipe = 0.5*38.;    //mm// Radius of beam pipe
  Double_t outerRadius =0.5*496.;   //mm// Outer radius of barrel
  Double_t halfLength = 0.5*1554.;  //mm Half length of barrel

// spd elements for ITS layers (#1-5)

  Double_t spd_sensor_l = 30.0;      //mm// Length of sensor
  Double_t spd_sensor_h = 15.0;      //mm// Height of sensor
  Double_t spd_sensorThick_IB1 = 0.050;   //mm// Thickness of sensor for inner barrels
  Double_t spd_sensorThick_IB2 = 0.050;   //mm// Thickness of sensor for inner barrels
  Double_t spd_sensorThick_IB3 = 0.050;   //mm// Thickness of sensor for inner barrels
  Double_t spd_sensorThick_OB4 = 0.700;   //mm// Thickness of sensor for outer barrels
  Double_t spd_sensorThick_OB5 = 0.700;   //mm// Thickness of sensor for outer barrels

  Double_t spd_deltaX = 0.1;        //mm// Distance in X between sensors in ladder
  Double_t spd_deltaZ = 0.1;        //mm// Distance in Z between sensors in ladder
  Double_t spd_stock = 0.2;         //mm// Distance from sensors to the frame of mother volume 
  Double_t delta = 4.0;             //mm// Overlapping of 2 ladders in 3- ayers

// layers

  const Int_t Nlayers = 5; // Number of layers in barrel

  Double_t halfZ[Nlayers];   // Half-size of layer along z in mm
  halfZ[0] = 0.5*750.0;   // layer #1
  halfZ[1] = 0.5*750.0;   // layer #2
  halfZ[2] = 0.5*750.0;   // layer #3
  halfZ[3] = 0.5*1526.0;   // layer #4
  halfZ[4] = 0.5*1526.8;  // layer #5

 
  Double_t radiusMin[Nlayers];        // Layer minimum radius 
  radiusMin[0] = 22.4;   // layer #1
  radiusMin[1] = 40.7;   // layer #2
  radiusMin[2] = 59.8;   // layer #3
  radiusMin[3] = 144.5;  // layer #4
  radiusMin[4] = 194.4;  // layer #5
 

  Double_t radiusMax[Nlayers];        // Layer maximum radius 
  radiusMax[0] = 26.7;   // layer #1
  radiusMax[1] = 45.9;   // layer #2
  radiusMax[2] = 65.1;   // layer #3
  radiusMax[3] = 147.9;   // layer #4
  radiusMax[4] = 197.6;  // layer #5

  Double_t radiusCenter[Nlayers];     // layer center radius (for volume definition)

// ladders

  Double_t spd_ladderThick = 1.0;     //mm// Thickness of spd ladder

  Int_t Nladders[Nlayers];   // Number of ladders in each layer
  Nladders[0] = 12;         // layer #1
  Nladders[1] = 22;         // layer #2
  Nladders[2] = 32;         // layer #3
  Nladders[3] = 36;         // layer #4
  Nladders[4] = 48;         // layer #5

  Double_t roll[Nlayers]; // ladder roll in each layer in degrees
  roll[0] = 0.;          // layer #1
  roll[1] = 0.;          // layer #2
  roll[2] = 0.;          // layer #3
  roll[3] = 0.;          // layer #4
  roll[4] = 0.;          // layer #5

  
  Int_t Nelements[Nlayers]; // Number of sensors in one ladder in each layer
  Nelements[0] = 24;  // layer #1
  Nelements[1] = 24;  // layer #2
  Nelements[2] = 24;  // layer #3
  Nelements[3] = 98;  // layer #4
  Nelements[4] = 98;  // layer #5
 
  Int_t nSec1[Nlayers];     // Number of sectors 1 in one ladder in each layer
  Int_t nSec2[Nlayers];     // Number of sectors 2 in one ladder in each layer
  Int_t nSec3[Nlayers];     // Number of sectors 3 in one ladder in each layer
  nSec1[0] = 24; // layer #1
  nSec2[0] = 0; // layer #1
  nSec3[0] = 0; // layer #1
  nSec1[1] = 24; // layer #2
  nSec2[1] = 0; // layer #2
  nSec3[1] = 0; // layer #2
  nSec1[2] = 24; // layer #3
  nSec2[2] = 0;  // layer #3
  nSec3[2] = 0; // layer #3
  nSec1[3] = 98; // layer #4
  nSec2[3] = 0; // layer #4
  nSec3[3] = 0; // layer #4
  nSec1[4] = 98; // layer #5
  nSec2[4] = 0; // layer #5
  nSec3[4] = 0; // layer #5
 
  std::ostringstream cablename;

//====== end declaration =======================================

  // output file
  const char* filename = "its_5spd_ib50_ob700_bp38.geo";
  ofstream* f = new ofstream(filename, ios::out | ios::trunc);

  // helper streams
  std::ostringstream points, position, rotation;
  points.setf(ios::showpoint); points.setf(ios::fixed);
  points << setprecision(6);
  position.setf(ios::showpoint); position.setf(ios::fixed);
  position << setprecision(6);
  rotation.setf(ios::showpoint); rotation.setf(ios::fixed);
  rotation << setprecision(6);

  // its volume (tube)
  points.str(""); position.str(""); rotation.str("");
  points << 0.0 << " " << 0.0 << " " << halfLength + 5.0 << endl;
  points << radiusPipe << " " << outerRadius << endl;
  points << 0.0 << " " << 0.0 << " " << -halfLength - 5.0;

  // sts tube definition
  Mpdshape* tube = new Mpdshape(f, "sts01", "cave", "TUBE", "air",
				points.str());
  tube->SetSegment(0);
  tube->DumpWithIncrement();

  tube->SetSegment(0);
  tube->SetMotherVolumeName("sts01");
  tube->SetPosition(0.,0.,0.);
 
  Double_t x,y,z;
  
  // Radius pre-calculation
  
     for (Int_t i = 0; i < Nlayers; i++) {
	    radiusCenter[i] = 0.5*(radiusMin[i]+radiusMax[i]);
     }


// Layer 111111111111111111111111111111111111111111111111111

  // Ladder of 1-st layer (ladder1)
  Double_t xWidth = 0.5*spd_sensor_h + 0.5*spd_stock;
  Double_t yWidth = 0.5*(spd_ladderThick);
  Double_t zWidth = halfZ[0] + 0.5*spd_stock;

  points.str(""); position.str(""); rotation.str("");
  points << xWidth << " " << -yWidth << " " << -zWidth << endl;
  points << xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << -zWidth << endl;
  points << xWidth << " " << -yWidth << " " << zWidth << endl;
  points << xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << zWidth;


  // ladder1 (ladder definition for layer#1)
  Mpdshape* ladder1 = new Mpdshape(f, "sts01ladder1", "sts01", "BOX", "air",
				   points.str(), position.str());
  ladder1->SetSegment(1);

  // ladder1 positioning
  Double_t beta = 360./Nladders[0];
  Double_t beta_rad = beta*TMath::DegToRad();
  Double_t gamma_deg = 90.;
  Double_t gamma_rad = gamma_deg*TMath::DegToRad();
  Double_t alpha = 90.;

  radiusCenter[0]=TMath::Sqrt(0.5*(radiusMin[0]*radiusMin[0]+radiusMax[0]*radiusMax[0]-
0.5*spd_sensor_h*spd_sensor_h));
  Double_t theta=TMath::ASin((radiusMax[0]*radiusMax[0]-radiusCenter[0]*radiusCenter[0]-
0.25*spd_sensor_h*spd_sensor_h)/(spd_sensor_h*radiusCenter[0]));
  roll[0]=theta*TMath::RadToDeg();
  cout <<"Layer#1:  Radius= "<< radiusCenter[0] << " mm " << "Inclination angle = " << roll[0] << " deg " << endl;

  for (Int_t i = 0; i < Nladders[0]; i++){
    Double_t x = radiusCenter[0] * TMath::Cos(gamma_rad+beta_rad*i);
    Double_t y = radiusCenter[0] * TMath::Sin(gamma_rad+beta_rad*i);
    points.str(""); position.str(""); rotation.str("");
    points << xWidth << " " << -yWidth << " " << -zWidth << endl;
    points << xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << -zWidth << endl;
    points << xWidth << " " << -yWidth << " " << zWidth << endl;
    points << xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << zWidth;
      
    ladder1->SetPosition(x,y,0.);
    ladder1->SetRotation(roll[0]+alpha+(gamma_deg+beta*i),0.,180.);
    ladder1->SetPoints(points.str());
    ladder1->DumpWithIncrement();
  }

// one spd sensor 11 of ladder1 (sensor 15*30 mm)
  xWidth = 0.5*spd_sensor_h; 
  yWidth = 0.5*spd_sensorThick_IB1; 
  zWidth = 0.5*spd_sensor_l;
  
  points.str(""); position.str(""); rotation.str("");
  points << xWidth << " " << -yWidth << " " << -zWidth << endl;
  points << xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << -zWidth << endl;  
  points << xWidth << " " << -yWidth << " " << zWidth << endl;
  points << xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << zWidth;
  
  // sensor 11 definition
  Mpdshape* sector11 = new Mpdshape(f, "sts01sector11", "sts01ladder1", 
				     "BOX", "silicon", points.str(), position.str());
  sector11->SetSegment(1);
  sector11->SetMotherSegment(1);
  
  // sectors 11 positioning into ladder of 1-st layer  

 for (Int_t j = 0; j < nSec1[0]/2; j++) {    	      
    
      x = 0; 
      y = 0;
      z = 0.5*spd_sensor_l + j*(spd_sensor_l+spd_deltaZ);
 
    points.str(""); position.str(""); rotation.str("");
    points << xWidth << " " << -yWidth << " " << -zWidth << endl;
    points << xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << -zWidth << endl;	
    points << xWidth << " " << -yWidth << " " << zWidth << endl;
    points << xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << zWidth;
    
    sector11->SetPosition(x,y,z);
    sector11->SetPoints(points.str());
    sector11->DumpWithIncrement();        	 
  }  

 for (Int_t j = 0; j < nSec1[0]/2; j++) {    	      
    
      x = 0; 
      y = 0; 
      z = -spd_deltaZ - 0.5*spd_sensor_l - j*(spd_sensor_l+spd_deltaZ);

    points.str(""); position.str(""); rotation.str("");
    points << xWidth << " " << -yWidth << " " << -zWidth << endl;
    points << xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << -zWidth << endl;	
    points << xWidth << " " << -yWidth << " " << zWidth << endl;
    points << xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << zWidth;
    
    sector11->SetPosition(x,y,z);
    sector11->SetPoints(points.str());
    sector11->DumpWithIncrement();        	 
  }  

// Layer 222222222222222222222222222222222222222222222222222

  // Ladder of 2 layer (ladder2)
  xWidth = 0.5*spd_sensor_h + 0.5*spd_stock;
  yWidth = 0.5*(spd_ladderThick);
  zWidth = halfZ[1] + 0.5*spd_stock;

  points.str(""); position.str(""); rotation.str("");
  points << xWidth << " " << -yWidth << " " << -zWidth << endl;
  points << xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << -zWidth << endl;
  points << xWidth << " " << -yWidth << " " << zWidth << endl;
  points << xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << zWidth;


  // ladder2 (ladder definition for layer#2)
  Mpdshape* ladder2 = new Mpdshape(f, "sts01ladder2", "sts01", "BOX", "air",
				   points.str(), position.str());
  ladder2->SetSegment(1);

  // ladder2 positioning
  beta = 360./Nladders[1];
  beta_rad = beta*TMath::DegToRad();
  gamma_deg = 90.;
  gamma_rad = gamma_deg*TMath::DegToRad();
  alpha = 90.;
      radiusCenter[1]=TMath::Sqrt(0.5*(radiusMin[1]*radiusMin[1]+radiusMax[1]*radiusMax[1]-
0.5*spd_sensor_h*spd_sensor_h));
  theta=TMath::ASin((radiusMax[1]*radiusMax[1]-radiusCenter[1]*radiusCenter[1]-
0.25*spd_sensor_h*spd_sensor_h)/(spd_sensor_h*radiusCenter[1]));
  roll[1]=theta*TMath::RadToDeg();
  cout <<"Layer#2:  Radius= "<< radiusCenter[1] << " mm " << "Inclination angle = " << roll[1] << " deg " << endl;

  for (Int_t i = 0; i < Nladders[1]; i++){
    Double_t x = radiusCenter[1] * TMath::Cos(gamma_rad+beta_rad*i);
    Double_t y = radiusCenter[1] * TMath::Sin(gamma_rad+beta_rad*i);
    points.str(""); position.str(""); rotation.str("");
    points << xWidth << " " << -yWidth << " " << -zWidth << endl;
    points << xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << -zWidth << endl;
    points << xWidth << " " << -yWidth << " " << zWidth << endl;
    points << xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << zWidth;
      
    ladder2->SetPosition(x,y,0.);
    ladder2->SetRotation(roll[1]+alpha+(gamma_deg+beta*i),0.,180.);
    ladder2->SetPoints(points.str());
    ladder2->DumpWithIncrement();
  }

// one spd sensor 12 of ladder2 (sensor 15*30 mm)
  xWidth = 0.5*spd_sensor_h; 
  yWidth = 0.5*spd_sensorThick_IB2; 
  zWidth = 0.5*spd_sensor_l;
  
  points.str(""); position.str(""); rotation.str("");
  points << xWidth << " " << -yWidth << " " << -zWidth << endl;
  points << xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << -zWidth << endl;  
  points << xWidth << " " << -yWidth << " " << zWidth << endl;
  points << xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << zWidth;
  
  // sensor 12 definition
  Mpdshape* sector12 = new Mpdshape(f, "sts01sector12", "sts01ladder2", 
				     "BOX", "silicon", points.str(), position.str());
  sector12->SetSegment(1);
  sector12->SetMotherSegment(1);
  
  // sectors 12 positioning into ladder of 2 layer  

 for (Int_t j = 0; j < nSec1[1]/2; j++) {    	      
    
      x = 0; 
      y = 0;
      z = 0.5*spd_sensor_l + j*(spd_sensor_l+spd_deltaZ);
 
    points.str(""); position.str(""); rotation.str("");
    points << xWidth << " " << -yWidth << " " << -zWidth << endl;
    points << xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << -zWidth << endl;	
    points << xWidth << " " << -yWidth << " " << zWidth << endl;
    points << xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << zWidth;
    
    sector12->SetPosition(x,y,z);
    sector12->SetPoints(points.str());
    sector12->DumpWithIncrement();        	 
  }  

 for (Int_t j = 0; j < nSec1[1]/2; j++) {    	      
    
      x = 0; 
      y = 0; 
      z = -spd_deltaZ - 0.5*spd_sensor_l - j*(spd_sensor_l+spd_deltaZ);

    points.str(""); position.str(""); rotation.str("");
    points << xWidth << " " << -yWidth << " " << -zWidth << endl;
    points << xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << -zWidth << endl;	
    points << xWidth << " " << -yWidth << " " << zWidth << endl;
    points << xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << zWidth;
    
    sector12->SetPosition(x,y,z);
    sector12->SetPoints(points.str());
    sector12->DumpWithIncrement();        	 
  }  
 
// Layer 33333333333333333333333333333333333333333333333333333

  // Ladder of 3 layer (ladder3)
  xWidth = 0.5*spd_sensor_h + 0.5*spd_stock;
  yWidth = 0.5*(spd_ladderThick);
  zWidth = halfZ[2] + 0.5*spd_stock;

  points.str(""); position.str(""); rotation.str("");
  points << xWidth << " " << -yWidth << " " << -zWidth << endl;
  points << xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << -zWidth << endl;
  points << xWidth << " " << -yWidth << " " << zWidth << endl;
  points << xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << zWidth;


  // ladder3 (ladder definition for layer#3)
  Mpdshape* ladder3 = new Mpdshape(f, "sts01ladder3", "sts01", "BOX", "air",
				   points.str(), position.str());
  ladder3->SetSegment(1);

  // ladder3 positioning
  beta = 360./Nladders[2];
  beta_rad = beta*TMath::DegToRad();
  gamma_deg = 90.;
  gamma_rad = gamma_deg*TMath::DegToRad();
  alpha = 90.;

  radiusCenter[2]=TMath::Sqrt(0.5*(radiusMin[2]*radiusMin[2]+radiusMax[2]*radiusMax[2]-
0.5*spd_sensor_h*spd_sensor_h));
  theta=TMath::ASin((radiusMax[2]*radiusMax[2]-radiusCenter[2]*radiusCenter[2]-
0.25*spd_sensor_h*spd_sensor_h)/(spd_sensor_h*radiusCenter[2]));
  roll[2]=theta*TMath::RadToDeg();
  cout <<"Layer#3:  Radius= "<< radiusCenter[2] << " mm " << "Inclination angle = " << roll[2] << " deg " << endl;

  for (Int_t i = 0; i < Nladders[2]; i++){
    Double_t x = radiusCenter[2] * TMath::Cos(gamma_rad+beta_rad*i);
    Double_t y = radiusCenter[2] * TMath::Sin(gamma_rad+beta_rad*i);
    points.str(""); position.str(""); rotation.str("");
    points << xWidth << " " << -yWidth << " " << -zWidth << endl;
    points << xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << -zWidth << endl;
    points << xWidth << " " << -yWidth << " " << zWidth << endl;
    points << xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << zWidth;
      
    ladder3->SetPosition(x,y,0.);
    ladder3->SetRotation(roll[2]+alpha+(gamma_deg+beta*i),0.,180.);
    ladder3->SetPoints(points.str());
    ladder3->DumpWithIncrement();
  }

// one spd sensor 13 of ladder3 (sensor 15*30 mm)
  xWidth = 0.5*spd_sensor_h; 
  yWidth = 0.5*spd_sensorThick_IB3; 
  zWidth = 0.5*spd_sensor_l;
  
  points.str(""); position.str(""); rotation.str("");
  points << xWidth << " " << -yWidth << " " << -zWidth << endl;
  points << xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << -zWidth << endl;  
  points << xWidth << " " << -yWidth << " " << zWidth << endl;
  points << xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << zWidth;
  
  // sensor 13 definition
  Mpdshape* sector13 = new Mpdshape(f, "sts01sector13", "sts01ladder3", 
				     "BOX", "silicon", points.str(), position.str());
  sector13->SetSegment(1);
  sector13->SetMotherSegment(1);
  
  // sectors 13 positioning into ladder of 3 layer  

 for (Int_t j = 0; j < nSec1[2]/2; j++) {    	      
    
      x = 0; 
      y = 0;
      z = 0.5*spd_sensor_l + j*(spd_sensor_l+spd_deltaZ);
 
    points.str(""); position.str(""); rotation.str("");
    points << xWidth << " " << -yWidth << " " << -zWidth << endl;
    points << xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << -zWidth << endl;	
    points << xWidth << " " << -yWidth << " " << zWidth << endl;
    points << xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << zWidth;
    
    sector13->SetPosition(x,y,z);
    sector13->SetPoints(points.str());
    sector13->DumpWithIncrement();        	 
  }  

 for (Int_t j = 0; j < nSec1[2]/2; j++) {    	      
    
      x = 0; 
      y = 0; 
      z = -spd_deltaZ - 0.5*spd_sensor_l - j*(spd_sensor_l+spd_deltaZ);

    points.str(""); position.str(""); rotation.str("");
    points << xWidth << " " << -yWidth << " " << -zWidth << endl;
    points << xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << -zWidth << endl;	
    points << xWidth << " " << -yWidth << " " << zWidth << endl;
    points << xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << zWidth;
    
    sector13->SetPosition(x,y,z);
    sector13->SetPoints(points.str());
    sector13->DumpWithIncrement();        	 
  }  
 

// Layer 4444444444444444444444444444444444444444444444444444

  // Ladder of 4-th layer (ladder4)
  xWidth = spd_sensor_h + 0.5*spd_stock;
  yWidth = 0.5*(spd_ladderThick);
  zWidth = halfZ[3] + 0.5*spd_stock;

  points.str(""); position.str(""); rotation.str("");
  points << xWidth << " " << -yWidth << " " << -zWidth << endl;
  points << xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << -zWidth << endl;
  points << xWidth << " " << -yWidth << " " << zWidth << endl;
  points << xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << zWidth;

  // ladder4 (ladder definition for layer#4)
  Mpdshape* ladder4 = new Mpdshape(f, "sts01ladder4", "sts01", "BOX", "air",  points.str(), position.str());
  ladder4->SetSegment(1);

  // ladder4 positioning
  beta = 360./Nladders[3];
  beta_rad = beta*TMath::DegToRad();
  gamma_deg = 90.;
  gamma_rad = gamma_deg*TMath::DegToRad();
  alpha = 90.;
  Double_t shift1 = TMath::ATan((0.5*spd_sensor_l-delta)/radiusMin[3]);
  Double_t shift2 = TMath::ATan(0.5*spd_sensor_l/radiusMax[3]);

 Double_t Radius = TMath::Sqrt(radiusMin[3]*radiusMin[3]+(0.5*spd_sensor_l-delta)*(0.5*spd_sensor_l-delta));
    for (Int_t i = 0; i < Nladders[3]/2; i++) {
    Double_t x = Radius * TMath::Cos(gamma_rad+2.*beta_rad*i+shift1);
    Double_t y = Radius * TMath::Sin(gamma_rad+2.*beta_rad*i+shift1);
    points.str(""); position.str(""); rotation.str("");
    points << xWidth << " " << -yWidth << " " << -zWidth << endl;
    points << xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << -zWidth << endl;
    points << xWidth << " " << -yWidth << " " << zWidth << endl;
    points << xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << zWidth;

    ladder4->SetPosition(x,y,0.);
    ladder4->SetRotation(alpha+(gamma_deg+2.*beta*i),0.,180.);
    ladder4->SetPoints(points.str());
    ladder4->DumpWithIncrement();
  }

    Radius = TMath::Sqrt(radiusMax[3]*radiusMax[3]+(0.5*spd_sensor_l)*(0.5*spd_sensor_l));
    for (Int_t i = 0; i < Nladders[3]/2; i++) {
    Double_t x = Radius * TMath::Cos(gamma_rad+2.*beta_rad*i-shift2);
    Double_t y = Radius * TMath::Sin(gamma_rad+2.*beta_rad*i-shift2);
    points.str(""); position.str(""); rotation.str("");
    points << xWidth << " " << -yWidth << " " << -zWidth << endl;
    points << xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << -zWidth << endl;
    points << xWidth << " " << -yWidth << " " << zWidth << endl;
    points << xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << zWidth;

    ladder4->SetPosition(x,y,0.);
    ladder4->SetRotation(alpha+(gamma_deg+2.*beta*i),0.,180.);
    ladder4->SetPoints(points.str());
    ladder4->DumpWithIncrement();
  }

// one spd sensor 14 of ladder4 (sensor 15*30 mm)
  xWidth = 0.5*spd_sensor_h; 
  yWidth = 0.5*spd_sensorThick_OB4; 
  zWidth = 0.5*spd_sensor_l;
  
  points.str(""); position.str(""); rotation.str("");
  points << xWidth << " " << -yWidth << " " << -zWidth << endl;
  points << xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << -zWidth << endl;  
  points << xWidth << " " << -yWidth << " " << zWidth << endl;
  points << xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << zWidth;
  
  // sector 14 definition
  Mpdshape* sector14 = new Mpdshape(f, "sts01sector14", "sts01ladder4", "BOX", "silicon", points.str(), position.str());
  sector14->SetSegment(1);
  sector14->SetMotherSegment(1);
  
  // sectors 14 positioning into ladder of 4-th layer  

 for (Int_t j = 0; j <= nSec1[3]/4; j++) {    	      
    
      x = 0.5*spd_sensor_h; 
      y = 0;
      z = j*(spd_sensor_l+spd_deltaZ);
 
    points.str(""); position.str(""); rotation.str("");
    points << xWidth << " " << -yWidth << " " << -zWidth << endl;
    points << xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << -zWidth << endl;	
    points << xWidth << " " << -yWidth << " " << zWidth << endl;
    points << xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << zWidth;
    
    sector14->SetPosition(x,y,z);
    sector14->SetPoints(points.str());
    sector14->DumpWithIncrement();        	 
  }  

 for (Int_t j = 1; j <= nSec1[3]/4; j++) {    	      
    
      x = 0.5*spd_sensor_h; 
      y = 0; 
      z = -j*(spd_sensor_l+spd_deltaZ);

    points.str(""); position.str(""); rotation.str("");
    points << xWidth << " " << -yWidth << " " << -zWidth << endl;
    points << xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << -zWidth << endl;	
    points << xWidth << " " << -yWidth << " " << zWidth << endl;
    points << xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << zWidth;
    
    sector14->SetPosition(x,y,z);
    sector14->SetPoints(points.str());
    sector14->DumpWithIncrement();        	 
  }  

 for (Int_t j = 0; j <= nSec1[3]/4; j++) {    	      
    
      x = -0.5*spd_sensor_h - spd_deltaX; 
      y = 0; 
      z = j*(spd_sensor_l+spd_deltaZ);

    points.str(""); position.str(""); rotation.str("");
    points << xWidth << " " << -yWidth << " " << -zWidth << endl;
    points << xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << -zWidth << endl;	
    points << xWidth << " " << -yWidth << " " << zWidth << endl;
    points << xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << zWidth;
    
    sector14->SetPosition(x,y,z);
    sector14->SetPoints(points.str());
    sector14->DumpWithIncrement();        	 
  }  

for (Int_t j = 1; j <= nSec1[3]/4; j++) {    	      
    
      x = -0.5*spd_sensor_h - spd_deltaX; 
      y = 0;
      z = -j*(spd_sensor_l+spd_deltaZ);  
  
    points.str(""); position.str(""); rotation.str("");
    points << xWidth << " " << -yWidth << " " << -zWidth << endl;
    points << xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << -zWidth << endl;	
    points << xWidth << " " << -yWidth << " " << zWidth << endl;
    points << xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << zWidth;
    
    sector14->SetPosition(x,y,z);
    sector14->SetPoints(points.str());
    sector14->DumpWithIncrement();        	 
  }  

// Layer 5555555555555555555555555555555555555555555555555555555

  // Ladder of 5-th layer (ladder5)
  xWidth = spd_sensor_h + 0.5*spd_stock;
  yWidth = 0.5*(spd_ladderThick);
  zWidth = halfZ[4] + 0.5*spd_stock;

  points.str(""); position.str(""); rotation.str("");
  points << xWidth << " " << -yWidth << " " << -zWidth << endl;
  points << xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << -zWidth << endl;
  points << xWidth << " " << -yWidth << " " << zWidth << endl;
  points << xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << zWidth;

  // ladder5 (ladder definition for layer#5)
  Mpdshape* ladder5 = new Mpdshape(f, "sts01ladder5", "sts01", "BOX", "air",  points.str(), position.str());
  ladder5->SetSegment(1);

  // ladder5 positioning
  beta = 360./Nladders[4];
  beta_rad = beta*TMath::DegToRad();
  gamma_deg = 90.;
  gamma_rad = gamma_deg*TMath::DegToRad();
  alpha = 90.;
  shift1 = TMath::ATan((0.5*spd_sensor_l-delta)/radiusMin[4]);
  shift2 = TMath::ATan(0.5*spd_sensor_l/radiusMax[4]);

    Radius = TMath::Sqrt(radiusMin[4]*radiusMin[4]+(0.5*spd_sensor_l-delta)*(0.5*spd_sensor_l-delta));
    for (Int_t i = 0; i < Nladders[4]/2; i++) {
    Double_t x = Radius * TMath::Cos(gamma_rad+2.*beta_rad*i+shift1);
    Double_t y = Radius * TMath::Sin(gamma_rad+2.*beta_rad*i+shift1);
    points.str(""); position.str(""); rotation.str("");
    points << xWidth << " " << -yWidth << " " << -zWidth << endl;
    points << xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << -zWidth << endl;
    points << xWidth << " " << -yWidth << " " << zWidth << endl;
    points << xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << zWidth;

    ladder5->SetPosition(x,y,0.);
    ladder5->SetRotation(alpha+(gamma_deg+2.*beta*i),0.,180.);
    ladder5->SetPoints(points.str());
    ladder5->DumpWithIncrement();
  }

    Radius = TMath::Sqrt(radiusMax[4]*radiusMax[4]+(0.5*spd_sensor_l)*(0.5*spd_sensor_l));
    for (Int_t i = 0; i < Nladders[4]/2; i++) {
    Double_t x = Radius * TMath::Cos(gamma_rad+2.*beta_rad*i-shift2);
    Double_t y = Radius * TMath::Sin(gamma_rad+2.*beta_rad*i-shift2);
    points.str(""); position.str(""); rotation.str("");
    points << xWidth << " " << -yWidth << " " << -zWidth << endl;
    points << xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << -zWidth << endl;
    points << xWidth << " " << -yWidth << " " << zWidth << endl;
    points << xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << zWidth;

    ladder5->SetPosition(x,y,0.);
    ladder5->SetRotation(alpha+(gamma_deg+2.*beta*i),0.,180.);
    ladder5->SetPoints(points.str());
    ladder5->DumpWithIncrement();
  }

// one spd sensor 15 of ladder5 (sensor 15*30 mm)
  xWidth = 0.5*spd_sensor_h; 
  yWidth = 0.5*spd_sensorThick_OB5; 
  zWidth = 0.5*spd_sensor_l;
  
  points.str(""); position.str(""); rotation.str("");
  points << xWidth << " " << -yWidth << " " << -zWidth << endl;
  points << xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << -zWidth << endl;  
  points << xWidth << " " << -yWidth << " " << zWidth << endl;
  points << xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << zWidth;
  
  // sector 15 definition
  Mpdshape* sector15 = new Mpdshape(f, "sts01sector15", "sts01ladder5", "BOX", "silicon", points.str(), position.str());
  sector15->SetSegment(1);
  sector15->SetMotherSegment(1);
  
  // sectors 15 positioning into ladder of 5-th layer  

 for (Int_t j = 0; j <= nSec1[4]/4; j++) {    	      
    
      x = 0.5*spd_sensor_h; 
      y = 0;
      z = j*(spd_sensor_l+spd_deltaZ);
 
    points.str(""); position.str(""); rotation.str("");
    points << xWidth << " " << -yWidth << " " << -zWidth << endl;
    points << xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << -zWidth << endl;	
    points << xWidth << " " << -yWidth << " " << zWidth << endl;
    points << xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << zWidth;
    
    sector15->SetPosition(x,y,z);
    sector15->SetPoints(points.str());
    sector15->DumpWithIncrement();        	 
  }  

 for (Int_t j = 1; j <= nSec1[4]/4; j++) {    	      
    
      x = 0.5*spd_sensor_h; 
      y = 0; 
      z = -j*(spd_sensor_l+spd_deltaZ);

    points.str(""); position.str(""); rotation.str("");
    points << xWidth << " " << -yWidth << " " << -zWidth << endl;
    points << xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << -zWidth << endl;	
    points << xWidth << " " << -yWidth << " " << zWidth << endl;
    points << xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << zWidth;
    
    sector15->SetPosition(x,y,z);
    sector15->SetPoints(points.str());
    sector15->DumpWithIncrement();        	 
  }  

 for (Int_t j = 0; j <= nSec1[4]/4; j++) {    	      
    
      x = -0.5*spd_sensor_h - spd_deltaX; 
      y = 0; 
      z = j*(spd_sensor_l+spd_deltaZ);

    points.str(""); position.str(""); rotation.str("");
    points << xWidth << " " << -yWidth << " " << -zWidth << endl;
    points << xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << -zWidth << endl;	
    points << xWidth << " " << -yWidth << " " << zWidth << endl;
    points << xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << zWidth;
    
    sector15->SetPosition(x,y,z);
    sector15->SetPoints(points.str());
    sector15->DumpWithIncrement();        	 
  }  

for (Int_t j = 1; j <= nSec1[4]/4; j++) {    	      
    
      x = -0.5*spd_sensor_h - spd_deltaX; 
      y = 0;
      z = -j*(spd_sensor_l+spd_deltaZ);  
  
    points.str(""); position.str(""); rotation.str("");
    points << xWidth << " " << -yWidth << " " << -zWidth << endl;
    points << xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << -zWidth << endl;	
    points << xWidth << " " << -yWidth << " " << zWidth << endl;
    points << xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << zWidth;
    
    sector15->SetPosition(x,y,z);
    sector15->SetPoints(points.str());
    sector15->DumpWithIncrement();        	 
  }  

   f->close(); 
   cout << "Finish " << endl;
}

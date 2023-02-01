{
  // 25.03.2010

#include "mpdshape.class.C"
  const Int_t Nlayer = 4; // Number of layers in barrel
  Double_t radiusPipe = 0.5*100.; //mm// Radius of pipe 
  
  // silicon elements of barrel
  Double_t detector_l = 75.0; //mm// Length of detector
  Double_t detector_h = 42.0; //mm// Height of detector
  Double_t dead_frame = 1.;//mm// Framework from silicon_dead in detector
  Double_t stock_0_2 = 0.2;//mm// Stock in mother_volume (detector)
  Double_t stock_0_3 = 0.3;//mm// Stock in mother_volume (ladder)    
  Double_t sensor_l = detector_l - 2.*dead_frame; //mm// Length of sensitive area of detector
  Double_t sensor_h = detector_h - 2.*dead_frame; //mm// Height of sensitive area of detector
  Double_t sensorThick = 0.3; //mm// Thickness of sensor and detector

  Double_t deltaRadius2 = 4.; //mm// Distance between circles in 2-nd layer
  Double_t deltaRadius3 = 3.; //mm// Distance between circles in 3-d layer
  Double_t deltaRadius4 = 3.; //mm// Distance between circles in 4-d layer
  Double_t deltaSensor = 0.3; //mm// Distance between detectors (or sensors) in ladder
  Double_t delta = 3.; //mm// Overlaping of detectors in ladder along z

  Double_t radius[Nlayer]; // Radius of layer  
  Double_t radius[0] = 0.5*124.+ sensorThick + 0.5*deltaSensor; //mm// 1-st layer
  Double_t radius[1] = 0.5*212.+ sensorThick + 0.5*deltaSensor; //mm// 2-nd layer
  Double_t radius[2] = 0.5*314.+ sensorThick + 0.5*deltaSensor; //mm// 3-d layer
  Double_t radius[3] = 0.5*404.+ sensorThick + 0.5*deltaSensor; //mm// 4-th layer 

  Double_t halfZ[Nlayer]; // Half-size of layers of barrel along z
  Double_t halfZ[0] = 0.5*627.; //mm// in 1-st layer
  Double_t halfZ[1] = 0.5*744.; //mm// in 2-nd layer
  Double_t halfZ[2] = 0.5*822.; //mm// in 3-d layer
  Double_t halfZ[3] = 0.5*939.; //mm// in 4-th layer

  Int_t Nladders[Nlayer]; // Number of ladders in each layer
  Int_t Nladders[0] = 6;  // in 1-st layer 
  Int_t Nladders[1] = 10; // in 2-nd layer 
  Int_t Nladders[2] = 14; // in 3-d layer 
  Int_t Nladders[3] = 18; // in 4-th layer 

  Int_t Nelements[Nlayer]; // Number of detectors in one ladder in each layer
  Int_t Nelements[0] = 16; // in 1-st layer 
  Int_t Nelements[1] = 19; // in 2-nd layer 
  Int_t Nelements[2] = 21; // in 3-d layer 
  Int_t Nelements[3] = 24; // in 4-th layer 
   
  // silicon elements of disk
  Double_t initDist = 0.5*1280.; //mm// Distance from interacting point to 1-st disk 
  Double_t trapThick = 0.5*0.3; //mm// Half-thickness of silicon trapezoid
  Double_t gamma = 15.; //degrees// Angle size of elements on disk
  Double_t gamma_rad = 2.0 * TMath::Pi() * gamma/360.0; // The same in radian      
  Double_t distTrap= 0.3; //mm// Distance between trapezoids in disk
  Double_t distDisk = 50.; //mm// Distance between disks along z
  Double_t h1 = 0.5 * 56.; //mm// Half-height of 1-st trapezoid (trap1)
  Double_t h2 = 0.5 * 72.; //mm// Half-height of 2-nd trapezoid (trap2)
  Double_t h3 = 0.5 * 67.2; //mm// Half-height of 3-d trapezoid (trap3)
  Double_t a0 = radiusPipe + distTrap; //mm// Distance to 1-st base of trap1 
  Double_t a1 = a0 + 2.*h1 + distTrap; //mm// Distance to 1-st base of trap2 
  Double_t a2 = a1 + 2.*h2 + distTrap; //mm// Distance to 1-st base of trap3 
  Double_t a3 = a2 + 2.*h3; //mm// Distance to 2-nd base of trap3
  Double_t b1 = TMath::Tan(0.5*gamma_rad)*a0; //mm// 1-st base of trap1
  Double_t b11 = TMath::Tan(0.5*gamma_rad)*(a0+2.*h1); //mm// 2-nd base of trap1
  Double_t b2 = TMath::Tan(0.5*gamma_rad)*a1; //mm// 1-st base of trap2
  Double_t b22 = TMath::Tan(0.5*gamma_rad)*(a1+2.*h2); //mm// 2-nd base of trap2
  Double_t b3 = TMath::Tan(0.5*gamma_rad)*a2; //mm// 1-st base of trap3
  Double_t b33 = TMath::Tan(0.5*gamma_rad)*a3; //mm// 2-nd base of trap3
  Double_t innerRadius = radiusPipe + trapThick; //mm// Inner radius of disk
  Double_t outerRadius = 0.5*500.; // Outer radius of disk

  //====== end declaration =======================================

  // output file
  const char* filename = "sts_v3.geo";
  ofstream* f = new ofstream(filename, ios::out | ios::trunc);                                                                            
  // helper streams
  std::ostringstream points, position, rotation;
  points.setf(ios::showpoint); points.setf(ios::fixed);
  points << setprecision(6);
  position.setf(ios::showpoint); position.setf(ios::fixed);
  position << setprecision(6);
  rotation.setf(ios::showpoint); rotation.setf(ios::fixed);
  rotation << setprecision(6);

  //***************************************************************

  // sts volume tube
  points.str(""); position.str(""); rotation.str("");
  points << 0.0 << " " << 0.0 << " " << halfZ[3] << endl;
  points << innerRadius << " " << outerRadius << endl;
  points << 0.0 << " " << 0.0 << " " << -halfZ[3];
  
  // sts tube definition
  Mpdshape* tube = new Mpdshape(f, "sts01", "cave", "TUBE", "air",
				points.str());
  tube->SetSegment(0);
  tube->DumpWithIncrement();
  
  tube->SetSegment(0);
  tube->SetMotherVolumeName("sts01");
  tube->SetPosition(0.,0.,0.);
  
  //111111111111111111111111111111111111111111111111111111111111111
  
  // one ladder of 1-st layer (ladder1)
  Double_t xWidth = 0.5*(detector_l + stock_0_3); 
  Double_t yWidth = 0.5*(2.*sensorThick + deltaSensor + stock_0_3); 
  Double_t zWidth = halfZ[0] + 0.5*stock_0_3;
  
  points.str(""); position.str(""); rotation.str("");
  points << xWidth << " " << -yWidth << " " << -zWidth << endl;
  points << xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << -zWidth << endl;  
  points << xWidth << " " << -yWidth << " " << zWidth << endl;
  points << xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << zWidth;
  
  // ladder1 definition
  Mpdshape* ladder1 = new Mpdshape(f, "sts01ladder1", "sts01", "BOX", "air",   
				   points.str(), position.str());
  ladder1->SetSegment(1);
  
  // output of ladder1
  Double_t beta = 360./Nladders[0];             
  Double_t beta_rad = beta*TMath::DegToRad(); 
  Double_t alpha = 96.;
  Double_t Radius = radius[0]; 	 
  Int_t ii = 0, jj = 0; 
  for (Int_t i = 0; i < Nladders[0]; i++) {    	
    
    Double_t x = Radius * TMath::Cos(beta_rad*i);
    Double_t y = Radius * TMath::Sin(beta_rad*i);     
    
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
    ladder1->SetRotation(alpha + beta*i,0.,0.);
    ladder1->SetPoints(points.str());
    ladder1->DumpWithIncrement();        	 
  }
  
  // one detector of ladder1 (detector1)
  Double_t xWidth = 0.5*(detector_l + stock_0_2); 
  Double_t yWidth = 0.5*(sensorThick + stock_0_2); 
  Double_t zWidth = 0.5*(detector_h + stock_0_2);
  
  points.str(""); position.str(""); rotation.str("");
  points << xWidth << " " << -yWidth << " " << -zWidth << endl;
  points << xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << -zWidth << endl;  
  points << xWidth << " " << -yWidth << " " << zWidth << endl;
  points << xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << zWidth;
  
  // detector1 definition
  Mpdshape* detector1 = new Mpdshape(f, "sts01detector1", "sts01ladder1", 
				     "BOX", "air", points.str(), position.str());
  detector1->SetSegment(1);
  detector1->SetMotherSegment(1);
  
  // output detectors of 1-st laeyr  
  for (Int_t j = 0; j < Nelements[0]; j++) {    	      
    jj=jj+1;
    if(jj%2 == 0) y = -0.5*(sensorThick + deltaSensor);
    if(jj%2 != 0) y = 0.5*(sensorThick + deltaSensor);  	   
    Double_t z = (-halfZ[0]+0.5*detector_h) + j*(detector_h-delta);
    
    points.str(""); position.str(""); rotation.str("");
    points << xWidth << " " << -yWidth << " " << -zWidth << endl;
    points << xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << -zWidth << endl;	
    points << xWidth << " " << -yWidth << " " << zWidth << endl;
    points << xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << zWidth;
    
    
    detector1->SetPosition(0.,y,z);
    detector1->SetPoints(points.str());
    detector1->DumpWithIncrement();        	 
  }
  
  // sensitive area of detector1 (sensor1)
  Double_t xWidth = 0.5*sensor_l; 
  Double_t yWidth = 0.5*sensorThick; 
  Double_t zWidth = 0.5*sensor_h;
  
  points.str(""); position.str(""); rotation.str("");
  points << xWidth << " " << -yWidth << " " << -zWidth << endl;
  points << xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << -zWidth << endl;  
  points << xWidth << " " << -yWidth << " " << zWidth << endl;
  points << xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << zWidth;
  
  // sensor1 definition
  Mpdshape* sensor1 = new Mpdshape(f, "sts01sensor1", "sts01detector1", 
				   "BOX", "silicon", points.str(), position.str());
  sensor1->SetSegment(0);
  sensor1->SetMotherSegment(1);
  sensor1->SetPosition(0.,0.,0.);   
  sensor1->DumpToFile(); 
  
  // vertical frame of detector
  Double_t xWidth = 0.5*dead_frame; 
  Double_t yWidth = 0.5*sensorThick; 
  Double_t zWidth = 0.5*detector_h;
  
  points.str(""); position.str(""); rotation.str("");
  points << xWidth << " " << -yWidth << " " << -zWidth << endl;
  points << xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << -zWidth << endl;  
  points << xWidth << " " << -yWidth << " " << zWidth << endl;
  points << xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << zWidth;
  
  Mpdshape* frame_v = new Mpdshape(f, "sts01frame_v", "sts01detector1", 
				   "BOX", "silicon_dead", points.str(), position.str());
  
  frame_v->SetSegment(1);
  frame_v->SetMotherSegment(1);
  
  frame_v->SetPosition(-0.5*sensor_l - 0.5*dead_frame,0.,0.);
  frame_v->DumpWithIncrement(); 
  
  frame_v->SetPosition(0.5*sensor_l + 0.5*dead_frame,0.,0.);
  frame_v->DumpWithIncrement(); 
  
  // horizontal frame of detector  
  Double_t xWidth = 0.5*sensor_l; 
  Double_t yWidth = 0.5*sensorThick; 
  Double_t zWidth = 0.5*dead_frame;
  
  points.str(""); position.str(""); rotation.str("");
  points << xWidth << " " << -yWidth << " " << -zWidth << endl;
  points << xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << -zWidth << endl;  
  points << xWidth << " " << -yWidth << " " << zWidth << endl;
  points << xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << zWidth;
  
  Mpdshape* frame_h = new Mpdshape(f, "sts01frame_h", "sts01detector1", 
				   "BOX", "silicon_dead", points.str(), position.str());
  
  frame_h->SetSegment(1);
  frame_h->SetMotherSegment(1);
  
  frame_h->SetPosition(0.,0.,0.5*sensor_h+0.5*dead_frame);
  frame_h->DumpWithIncrement();
  
  frame_h->SetPosition(0.,0.,-0.5*sensor_h-0.5*dead_frame);
  frame_h->DumpWithIncrement();
  
  //22222222222222222222222222222222222222222222222222222222222222222222222222222222
    
  // one ladder of 2-st layer (ladder2)
  Double_t xWidth = 0.5*(detector_l + stock_0_3); 
  Double_t yWidth = 0.5*(2.*sensorThick + deltaSensor + stock_0_3);   
  Double_t zWidth = halfZ[1] + 0.5*stock_0_3;
  
  points.str(""); position.str(""); rotation.str("");
  points << xWidth << " " << -yWidth << " " << -zWidth << endl;
  points << xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << -zWidth << endl;  
  points << xWidth << " " << -yWidth << " " << zWidth << endl;
  points << xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << zWidth;
  
  // ladder2 definition
  Mpdshape* ladder2 = new Mpdshape(f, "sts01ladder2", "sts01", "BOX", "air",   
				   points.str(), position.str());
  ladder2->SetSegment(1);
  
  // output of ladder2
  Double_t beta = 360./Nladders[1];      
  Double_t beta_rad = beta*TMath::DegToRad();        
  Double_t alpha = 90.; 	 
  Int_t ii = 0, jj = 0;
  
  for (Int_t i = 0; i < Nladders[1]; i++) {    	
    ii = ii+1;
    if(ii%2 == 0) Radius = radius[1] + deltaRadius2;
    if(ii%2 != 0) Radius = radius[1];
    Double_t x = Radius * TMath::Cos(beta_rad*i);
    Double_t y = Radius * TMath::Sin(beta_rad*i);     
    
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
    ladder2->SetRotation(alpha + beta*i,0.,0.);
    ladder2->SetPoints(points.str());
    ladder2->DumpWithIncrement();        	 
  }
  
  // one detector of ladder2 (detector2)
  Double_t xWidth = 0.5*(detector_l + stock_0_2); 
  Double_t yWidth = 0.5*(sensorThick + stock_0_2); 
  Double_t zWidth = 0.5*(detector_h + stock_0_2);
  
  points.str(""); position.str(""); rotation.str("");
  points << xWidth << " " << -yWidth << " " << -zWidth << endl;
  points << xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << -zWidth << endl;  
  points << xWidth << " " << -yWidth << " " << zWidth << endl;
  points << xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << zWidth;
  
  // detector2 definition
  Mpdshape* detector2 = new Mpdshape(f, "sts01detector2", "sts01ladder2", 
				     "BOX", "air", points.str(), position.str());
  detector2->SetSegment(1);
  detector2->SetMotherSegment(1);
  
  // output detectors of 2-st laeyr  
  for (Int_t j = 0; j < Nelements[1]; j++) {    	      
    jj=jj+1;
    if(jj%2 == 0) y = -0.5*(sensorThick + deltaSensor);
    if(jj%2 != 0) y = 0.5*(sensorThick + deltaSensor);      	   
    Double_t z = (-halfZ[1]+0.5*detector_h) + j*(detector_h-delta);
    
    points.str(""); position.str(""); rotation.str("");
    points << xWidth << " " << -yWidth << " " << -zWidth << endl;
    points << xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << -zWidth << endl;	
    points << xWidth << " " << -yWidth << " " << zWidth << endl;
    points << xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << zWidth;
    
    detector2->SetPosition(0.,y,z);
    detector2->SetPoints(points.str());
    detector2->DumpWithIncrement();        	 
  }
  
  // sensitive area of detector2 (sensor2)
  Double_t xWidth = 0.5*sensor_l; 
  Double_t yWidth = 0.5*sensorThick; 
  Double_t zWidth = 0.5*sensor_h;
  
  points.str(""); position.str(""); rotation.str("");
  points << xWidth << " " << -yWidth << " " << -zWidth << endl;
  points << xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << -zWidth << endl;  
  points << xWidth << " " << -yWidth << " " << zWidth << endl;
  points << xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << zWidth;
  
  // sensor2 definition
  Mpdshape* sensor2 = new Mpdshape(f, "sts01sensor2", "sts01detector2", 
				   "BOX", "silicon", points.str(), position.str());
  sensor2->SetSegment(0);
  sensor2->SetMotherSegment(1);
  sensor2->SetPosition(0.,0.,0.);   
  sensor2->DumpToFile(); 
  
  // vertical frame of detector
  Double_t xWidth = 0.5*dead_frame; 
  Double_t yWidth = 0.5*sensorThick; 
  Double_t zWidth = 0.5*detector_h;
  
  points.str(""); position.str(""); rotation.str("");
  points << xWidth << " " << -yWidth << " " << -zWidth << endl;
  points << xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << -zWidth << endl;  
  points << xWidth << " " << -yWidth << " " << zWidth << endl;
  points << xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << zWidth;
  
  Mpdshape* frame_v = new Mpdshape(f, "sts01frame_v", "sts01detector2", 
				   "BOX", "silicon_dead", points.str(), position.str());
  frame_v->SetSegment(1);
  frame_v->SetMotherSegment(1);
  
  frame_v->SetPosition(-0.5*sensor_l - 0.5*dead_frame,0.,0.);
  frame_v->DumpWithIncrement(); 
  
  frame_v->SetPosition(0.5*sensor_l + 0.5*dead_frame,0.,0.);
  frame_v->DumpWithIncrement();
    
  // horizontal frame of detector
  Double_t xWidth = 0.5*sensor_l; 
  Double_t yWidth = 0.5*sensorThick; 
  Double_t zWidth = 0.5*dead_frame;
  
  points.str(""); position.str(""); rotation.str("");
  points << xWidth << " " << -yWidth << " " << -zWidth << endl;
  points << xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << -zWidth << endl;  
  points << xWidth << " " << -yWidth << " " << zWidth << endl;
  points << xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << zWidth;
  
  Mpdshape* frame_h = new Mpdshape(f, "sts01frame_h", "sts01detector2", 
				   "BOX", "silicon_dead", points.str(), position.str());
  frame_h->SetSegment(1);
  frame_h->SetMotherSegment(1);
  
  frame_h->SetPosition(0.,0.,0.5*sensor_h+0.5*dead_frame);
  frame_h->DumpWithIncrement();
  
  frame_h->SetPosition(0.,0.,-0.5*sensor_h-0.5*dead_frame);
  frame_h->DumpWithIncrement();
  
  
  //33333333333333333333333333333333333333333333333333333333333333333333333333333333
    
  // one ladder of 3-st layer (ladder3)
  Double_t xWidth = 0.5*(detector_l + stock_0_3); 
  Double_t yWidth = 0.5*(2.*sensorThick + deltaSensor + stock_0_3);  
  Double_t zWidth = halfZ[2] + 0.5*stock_0_3;
  
  points.str(""); position.str(""); rotation.str("");
  points << xWidth << " " << -yWidth << " " << -zWidth << endl;
  points << xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << -zWidth << endl;  
  points << xWidth << " " << -yWidth << " " << zWidth << endl;
  points << xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << zWidth;
  
  // ladder3 definition
  Mpdshape* ladder3 = new Mpdshape(f, "sts01ladder3", "sts01", "BOX", "air",   
				   points.str(), position.str());
  ladder3->SetSegment(1);
  
  // output of ladder3
  Double_t beta = 360./Nladders[2];            
  Double_t beta_rad = beta*TMath::DegToRad();  
  Double_t alpha = 90.; 	 
  Int_t ii = 0, jj = 0; 

  for (Int_t i = 0; i < Nladders[2]; i++) {    	
    ii = ii+1;
    if(ii%2 == 0) Radius = radius[2] + deltaRadius3;
    if(ii%2 != 0) Radius = radius[2];    
    Double_t x = Radius * TMath::Cos(beta_rad*i);
    Double_t y = Radius * TMath::Sin(beta_rad*i);     
    
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
    ladder3->SetRotation(alpha + beta*i,0.,0.);
    ladder3->SetPoints(points.str());
    ladder3->DumpWithIncrement();        	 
  }
  
  // one detector of ladder3 (detector3)
  Double_t xWidth = 0.5*(detector_l + stock_0_2); 
  Double_t yWidth = 0.5*(sensorThick + stock_0_2); 
  Double_t zWidth = 0.5*(detector_h + stock_0_2); 
  
  points.str(""); position.str(""); rotation.str("");
  points << xWidth << " " << -yWidth << " " << -zWidth << endl;
  points << xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << -zWidth << endl;  
  points << xWidth << " " << -yWidth << " " << zWidth << endl;
  points << xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << zWidth;
  
  // detector3 definition
  Mpdshape* detector3 = new Mpdshape(f, "sts01detector3", "sts01ladder3", 
				     "BOX", "air", points.str(), position.str());
  detector3->SetSegment(1);
  detector3->SetMotherSegment(1);
  
  // output detectors of 3-st laeyr   
  for (Int_t j = 0; j < Nelements[2]; j++) {    	      
    jj=jj+1;
    if(jj%2 == 0) y = -0.5*(sensorThick + deltaSensor);
    if(jj%2 != 0) y = 0.5*(sensorThick + deltaSensor);            	   
    Double_t z = (-halfZ[2]+0.5*detector_h) + j*(detector_h-delta);
    
    points.str(""); position.str(""); rotation.str("");
    points << xWidth << " " << -yWidth << " " << -zWidth << endl;
    points << xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << -zWidth << endl;	
    points << xWidth << " " << -yWidth << " " << zWidth << endl;
    points << xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << zWidth;
    
    detector3->SetPosition(0.,y,z);
    detector3->SetPoints(points.str());
    detector3->DumpWithIncrement();        	 
  }
  
  // sensitive area of detector3 (sensor3)
  Double_t xWidth = 0.5*sensor_l; 
  Double_t yWidth = 0.5*sensorThick; 
  Double_t zWidth = 0.5*sensor_h;
  
  points.str(""); position.str(""); rotation.str("");
  points << xWidth << " " << -yWidth << " " << -zWidth << endl;
  points << xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << -zWidth << endl;  
  points << xWidth << " " << -yWidth << " " << zWidth << endl;
  points << xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << zWidth;
  
  // sensor3 definition
  Mpdshape* sensor3 = new Mpdshape(f, "sts01sensor3", "sts01detector3", 
				   "BOX", "silicon", points.str(), position.str());
  sensor3->SetSegment(0);
  sensor3->SetMotherSegment(1);
  sensor3->SetPosition(0.,0.,0.);   
  sensor3->DumpToFile(); 
  
  // vertical frame of detector
  Double_t xWidth = 0.5*dead_frame; 
  Double_t yWidth = 0.5*sensorThick; 
  Double_t zWidth = 0.5*detector_h;
  
  points.str(""); position.str(""); rotation.str("");
  points << xWidth << " " << -yWidth << " " << -zWidth << endl;
  points << xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << -zWidth << endl;  
  points << xWidth << " " << -yWidth << " " << zWidth << endl;
  points << xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << zWidth;
  
  Mpdshape* frame_v = new Mpdshape(f, "sts01frame_v", "sts01detector3", 
				   "BOX", "silicon_dead", points.str(), position.str());
  frame_v->SetSegment(1);
  frame_v->SetMotherSegment(1);
  
  frame_v->SetPosition(-0.5*sensor_l - 0.5*dead_frame,0.,0.);
  frame_v->DumpWithIncrement(); 
  
  frame_v->SetPosition(0.5*sensor_l + 0.5*dead_frame,0.,0.);
  frame_v->DumpWithIncrement();  
  
  // horizontal frame of detector
  Double_t xWidth = 0.5*sensor_l; 
  Double_t yWidth = 0.5*sensorThick; 
  Double_t zWidth = 0.5*dead_frame;
  
  points.str(""); position.str(""); rotation.str("");
  points << xWidth << " " << -yWidth << " " << -zWidth << endl;
  points << xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << -zWidth << endl;  
  points << xWidth << " " << -yWidth << " " << zWidth << endl;
  points << xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << zWidth;
  
  Mpdshape* frame_h = new Mpdshape(f, "sts01frame_h", "sts01detector3", 
				   "BOX", "silicon_dead", points.str(), position.str());
  frame_h->SetSegment(1);
  frame_h->SetMotherSegment(1);
  
  frame_h->SetPosition(0.,0.,0.5*sensor_h+0.5*dead_frame);
  frame_h->DumpWithIncrement();
  
  frame_h->SetPosition(0.,0.,-0.5*sensor_h-0.5*dead_frame);
  frame_h->DumpWithIncrement();
  
  //44444444444444444444444444444444444444444444444444444444444444444444444444444444
    
  // one ladder of 4-st layer (ladder4)
  Double_t xWidth = 0.5*(detector_l + stock_0_3); 
  Double_t yWidth = 0.5*(2.*sensorThick + deltaSensor + stock_0_3);   
  Double_t zWidth = halfZ[3] + 0.5*stock_0_3;
  
  points.str(""); position.str(""); rotation.str("");
  points << xWidth << " " << -yWidth << " " << -zWidth << endl;
  points << xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << -zWidth << endl;  
  points << xWidth << " " << -yWidth << " " << zWidth << endl;
  points << xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << zWidth;
  
  // ladder4 definition
  Mpdshape* ladder4 = new Mpdshape(f, "sts01ladder4", "sts01", "BOX", "air",   
				   points.str(), position.str());
  ladder4->SetSegment(1);
  
  // output of ladder4
  Double_t beta = 360./Nladders[3];  
  Double_t beta_rad = beta*TMath::DegToRad();             
  Double_t alpha = 90.;  	 
  Int_t ii = 0, jj = 0; 

  for (Int_t i = 0; i < Nladders[3]; i++) {    	
    ii = ii+1;
    if(ii%2 == 0) Radius = radius[3] + deltaRadius4;
    if(ii%2 != 0) Radius = radius[3];     
    Double_t x = Radius * TMath::Cos(beta_rad*i);
    Double_t y = Radius * TMath::Sin(beta_rad*i);     
    
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
    ladder4->SetRotation(alpha + beta*i,0.,0.);
    ladder4->SetPoints(points.str());
    ladder4->DumpWithIncrement();        	 
  }
  
  // one detector of ladder4 (detector4)
  Double_t xWidth = 0.5*(detector_l + stock_0_2); 
  Double_t yWidth = 0.5*(sensorThick + stock_0_2); 
  Double_t zWidth = 0.5*(detector_h + stock_0_2);  
  
  points.str(""); position.str(""); rotation.str("");
  points << xWidth << " " << -yWidth << " " << -zWidth << endl;
  points << xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << -zWidth << endl;  
  points << xWidth << " " << -yWidth << " " << zWidth << endl;
  points << xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << zWidth;
  
  // detector4 definition
  Mpdshape* detector4 = new Mpdshape(f, "sts01detector4", "sts01ladder4", 
				     "BOX", "air", points.str(), position.str());
  detector4->SetSegment(1);
  detector4->SetMotherSegment(1);
  
  // output detectors of 4-st laeyr  
  for (Int_t j = 0; j < Nelements[3]; j++) {    	      
    jj=jj+1;
    if(jj%2 == 0) y = -0.5*(sensorThick + deltaSensor);
    if(jj%2 != 0) y = 0.5*(sensorThick + deltaSensor);   
    Double_t z = (-halfZ[3]+0.5*detector_h) + j*(detector_h-delta);
    
    points.str(""); position.str(""); rotation.str("");
    points << xWidth << " " << -yWidth << " " << -zWidth << endl;
    points << xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << yWidth << " " << -zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << -zWidth << endl;	
    points << xWidth << " " << -yWidth << " " << zWidth << endl;
    points << xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << yWidth << " " << zWidth << endl;
    points << -xWidth << " " << -yWidth << " " << zWidth;
    
    detector4->SetPosition(0.,y,z);
    detector4->SetPoints(points.str());
    detector4->DumpWithIncrement();        	 
  }
  
  // sensitive area of detector4 (sensor4)
  Double_t xWidth = 0.5*sensor_l; 
  Double_t yWidth = 0.5*sensorThick; 
  Double_t zWidth = 0.5*sensor_h;
  
  points.str(""); position.str(""); rotation.str("");
  points << xWidth << " " << -yWidth << " " << -zWidth << endl;
  points << xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << -zWidth << endl;  
  points << xWidth << " " << -yWidth << " " << zWidth << endl;
  points << xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << zWidth;
  
  // sensor4 definition
  Mpdshape* sensor4 = new Mpdshape(f, "sts01sensor4", "sts01detector4", 
				   "BOX", "silicon", points.str(), position.str());
  sensor4->SetSegment(0);
  sensor4->SetMotherSegment(1);
  sensor4->SetPosition(0.,0.,0.);   
  sensor4->DumpToFile(); 
  
  // vertical frame of detector
  Double_t xWidth = 0.5*dead_frame; 
  Double_t yWidth = 0.5*sensorThick; 
  Double_t zWidth = 0.5*detector_h;
  
  points.str(""); position.str(""); rotation.str("");
  points << xWidth << " " << -yWidth << " " << -zWidth << endl;
  points << xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << -zWidth << endl;  
  points << xWidth << " " << -yWidth << " " << zWidth << endl;
  points << xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << zWidth;
  
  Mpdshape* frame_v = new Mpdshape(f, "sts01frame_v", "sts01detector4", 
				   "BOX", "silicon_dead", points.str(), position.str());
  frame_v->SetSegment(1);
  frame_v->SetMotherSegment(1);
  
  frame_v->SetPosition(-0.5*sensor_l - 0.5*dead_frame,0.,0.);
  frame_v->DumpWithIncrement(); 
  
  frame_v->SetPosition(0.5*sensor_l + 0.5*dead_frame,0.,0.);
  frame_v->DumpWithIncrement();  

  // horizontal frame of detector
  Double_t xWidth = 0.5*sensor_l; 
  Double_t yWidth = 0.5*sensorThick; 
  Double_t zWidth = 0.5*dead_frame;
  
  points.str(""); position.str(""); rotation.str("");
  points << xWidth << " " << -yWidth << " " << -zWidth << endl;
  points << xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << yWidth << " " << -zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << -zWidth << endl;  
  points << xWidth << " " << -yWidth << " " << zWidth << endl;
  points << xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << yWidth << " " << zWidth << endl;
  points << -xWidth << " " << -yWidth << " " << zWidth;
  
  Mpdshape* frame_h = new Mpdshape(f, "sts01frame_h", "sts01detector4", 
				   "BOX", "silicon_dead", points.str(), position.str());
  frame_h->SetSegment(1);
  frame_h->SetMotherSegment(1);
  
  frame_h->SetPosition(0.,0.,0.5*sensor_h+0.5*dead_frame);
  frame_h->DumpWithIncrement();
  
  frame_h->SetPosition(0.,0.,-0.5*sensor_h-0.5*dead_frame);
  frame_h->DumpWithIncrement();
    
//******************************************************************
/*
  // silicon disk helper
  points.str(""); position.str(""); rotation.str("");
  points << 0.0 << " " << 0.0 << " " << -trapThick << endl;
  points << innerRadius << " " << outerRadius << endl;
  points << 0.0 << " " << 0.0 << " " << trapThick;
  position << 0.0 << " " << 0.0 << " " << initDist;
  
  // silicon disk definition
  Mpdshape* disk = new Mpdshape(f, "sts01disk", "cave", "TUBE",
				"silicon", points.str());      	
  disk->SetSegment(1);

  //1-st layer of trapezoids on disk
  Double_t x0Width = b1;
  Double_t x1Width = b11;
  Double_t yWidth = h1;
  Double_t zWidth = trapThick;

  points.str(""); position.str(""); rotation.str("");
  points << x1Width << " " << -yWidth << " " << -zWidth << endl;
  points << x0Width << " " << yWidth << " " << -zWidth << endl;
  points << -x0Width << " " << yWidth << " " << -zWidth << endl;
  points << -x1Width << " " << -yWidth << " " << -zWidth << endl;	    
  points << x1Width << " " << -yWidth << " " << zWidth << endl;
  points << x0Width << " " << yWidth << " " << zWidth << endl;
  points << -x0Width << " " << yWidth << " " << zWidth << endl;
  points << -x1Width << " " << -yWidth << " " << zWidth;
  
  // trapezoid definition
  Mpdshape* trap1 = new Mpdshape(f, "sts01trap1", "sts01disk", "TRAP",
				 "silicon", points.str());
  trap1->SetSegment(1);
  trap1->SetMotherSegment(1);

  //2-nd layer of trapezoid on disk
  x0Width = b2;
  x1Width = b22;
  yWidth = h2;
  zWidth = trapThick;

  points.str(""); position.str(""); rotation.str("");
  points << x1Width << " " << -yWidth << " " << -zWidth << endl;
  points << x0Width << " " << yWidth << " " << -zWidth << endl;
  points << -x0Width << " " << yWidth << " " << -zWidth << endl;
  points << -x1Width << " " << -yWidth << " " << -zWidth << endl;	    
  points << x1Width << " " << -yWidth << " " << zWidth << endl;
  points << x0Width << " " << yWidth << " " << zWidth << endl;
  points << -x0Width << " " << yWidth << " " << zWidth << endl;
  points << -x1Width << " " << -yWidth << " " << zWidth;
  
  // trapezoid definition
  Mpdshape* trap2 = new Mpdshape(f, "sts01trap2", "sts01disk", "TRAP",
				 "silicon", points.str());     
  trap2->SetSegment(1);
  trap2->SetMotherSegment(1);
  
  //3-d layer of trapezoid on disk
  x0Width = b3;
  x1Width = b33;
  yWidth = h3;
  zWidth = trapThick;

  points.str(""); position.str(""); rotation.str("");
  points << x1Width << " " << -yWidth << " " << -zWidth << endl;
  points << x0Width << " " << yWidth << " " << -zWidth << endl;
  points << -x0Width << " " << yWidth << " " << -zWidth << endl;
  points << -x1Width << " " << -yWidth << " " << -zWidth << endl;	    
  points << x1Width << " " << -yWidth << " " << zWidth << endl;
  points << x0Width << " " << yWidth << " " << zWidth << endl;
  points << -x0Width << " " << yWidth << " " << zWidth << endl;
  points << -x1Width << " " << -yWidth << " " << zWidth;
  
  // trapezoid definition 
  Mpdshape* trap3 = new Mpdshape(f, "sts01trap3", "sts01disk", "TRAP",
				 "silicon", points.str());     
  trap3->SetSegment(1);
  trap3->SetMotherSegment(1);
  
  // output silicon disks
  
  disk->SetPosition(0.,0.,initDist + distDisk);
  disk->DumpWithIncrement();
  disk->SetPosition(0.,0.,initDist + 2.*distDisk);
  disk->DumpWithIncrement();
  disk->SetPosition(0.,0.,initDist + 3.*distDisk);
  disk->DumpWithIncrement();
  disk->SetPosition(0.,0.,-initDist - distDisk);
  disk->DumpWithIncrement();
  disk->SetPosition(0.,0.,-initDist - 2.*distDisk);
  disk->DumpWithIncrement();
  disk->SetPosition(0.,0.,-initDist - 3.*distDisk);
  disk->DumpWithIncrement();
  
  // output 1-st layer of trapezoids
  
  for (Int_t k = 0; k < 24; k++) {              
    x = (a0+h1) * TMath::Cos(gamma_rad*k);
    y = (a0+h1) * TMath::Sin(gamma_rad*k);
    beta = 90.;        
    x0Width = b1;
    x1Width = b11;
    yWidth = h1;
    zWidth = trapThick;

    points.str(""); position.str(""); rotation.str("");
    points << x1Width << " " << -yWidth << " " << -zWidth << endl;
    points << x0Width << " " << yWidth << " " << -zWidth << endl;
    points << -x0Width << " " << yWidth << " " << -zWidth << endl;
    points << -x1Width << " " << -yWidth << " " << -zWidth << endl;	    
    points << x1Width << " " << -yWidth << " " << zWidth << endl;
    points << x0Width << " " << yWidth << " " << zWidth << endl;
    points << -x0Width << " " << yWidth << " " << zWidth << endl;
    points << -x1Width << " " << -yWidth << " " << zWidth;
    
    trap1->SetPosition(x,y,0.);
    trap1->SetRotation(beta+gamma*k,0.,0.);
    trap1->DumpWithIncrement();      
  }
  
  //output 2-nd layer of trapezoids
  
  for (Int_t k = 0; k < 24; k++) {       
    x = (a1+h2) * TMath::Cos(gamma_rad*k);
    y = (a1+h2) * TMath::Sin(gamma_rad*k);     
    beta = 90.; 
    x0Width = b2;
    x1Width = b22;
    yWidth = h2;
    zWidth = trapThick;

    points.str(""); position.str(""); rotation.str("");    
    points << x1Width << " " << -yWidth << " " << -zWidth << endl;
    points << x0Width << " " << yWidth << " " << -zWidth << endl;
    points << -x0Width << " " << yWidth << " " << -zWidth << endl;
    points << -x1Width << " " << -yWidth << " " << -zWidth << endl;	    
    points << x1Width << " " << -yWidth << " " << zWidth << endl;
    points << x0Width << " " << yWidth << " " << zWidth << endl;
    points << -x0Width << " " << yWidth << " " << zWidth << endl;
    points << -x1Width << " " << -yWidth << " " << zWidth;
    
    trap2->SetPoints(points.str());
    trap2->SetPosition(x,y,0.);
    trap2->SetRotation(beta+gamma*k,0.,0.);
    trap2->DumpWithIncrement();
  }
  
  //output 3-d layer of trapezoids
  
  for (Int_t k = 0; k < 24; k++) {       
    x = (a2+h3) * TMath::Cos(gamma_rad*k);
    y = (a2+h3) * TMath::Sin(gamma_rad*k);   
    beta = 90.;      
    x0Width = b3;
    x1Width = b33;
    yWidth = h3;
    zWidth = trapThick;

    points.str(""); position.str(""); rotation.str("");    
    points << x1Width << " " << -yWidth << " " << -zWidth << endl;
    points << x0Width << " " << yWidth << " " << -zWidth << endl;
    points << -x0Width << " " << yWidth << " " << -zWidth << endl;
    points << -x1Width << " " << -yWidth << " " << -zWidth << endl;	    
    points << x1Width << " " << -yWidth << " " << zWidth << endl;
    points << x0Width << " " << yWidth << " " << zWidth << endl;
    points << -x0Width << " " << yWidth << " " << zWidth << endl;
    points << -x1Width << " " << -yWidth << " " << zWidth;
    
    trap3->SetPoints(points.str());
    trap3->SetPosition(x,y,0.);
    trap3->SetRotation(beta+gamma*k,0.,0.);
    trap3->DumpWithIncrement();
  }
*/       
  // close geometry file                                                       
  f->close(); 
}

  


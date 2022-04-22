#include "EmcGeom.h"


ClassImp(EmcGeom);

EmcGeom::EmcGeom() {}

EmcGeom::EmcGeom(ofstream* f) {
  fGeoFile = f;
}

EmcGeom::~EmcGeom() {
  delete rot;
  
}

//_____________________________________________________________
void EmcGeom::BuildEMC() {
 Mpdshape* tube = new Mpdshape(fGeoFile, "emc1Chamber1#1", "cave", "TUBE","air");
 
 tube->Fill_TUBE(emc1Chamber_z_l/2, inr, outr);
 tube->SetSegment(0);
 tube->SetPosition(0., 0., 0.);
 tube->SetRotation(1., 0., 0., 0., 1., 0., 0., 0., 1.);

 tube->DumpToFile();

 BuildSensVolume();
  
  
 delete tube;
}

//_______________________________________________________________
void EmcGeom::BuildSensVolume() {
   
  Mpdshape* emc1Chamber2= new Mpdshape(fGeoFile, "emc1Chamber2", "emc1Chamber1#1", "TUBE", "air");

  emc1Chamber2->Fill_TUBE(emc1_box1_u/2, inr, outr);
  emc1Chamber2->SetSegment(1);
  emc1Chamber2->SetPosition(0. ,0. , -emc1Chamber_z_th);
  emc1Chamber2->SetRotation(1., 0., 0., 0., 1., 0., 0., 0., 1.);
  emc1Chamber2->DumpWithIncrement();
  
  Float_t zst;
  zst=0;
  for(int i=1; i<nbox1_inl; i++){
    zst= -1*emc1Chamber_z_th + emc1_box1_u*i;
    emc1Chamber2->SetSegment(i+1);
    emc1Chamber2->SetPosition(z, z, zst);
    emc1Chamber2->DumpToFile();
  }
  
//*************************************************
  Mpdshape* emc1_trap1 = new Mpdshape(fGeoFile, "emc1_trap1", "emc1Chamber2#1", "TUBS", "air");
  
  emc1_trap1->Fill_TUBS(-emc1_box1_u/2, emc1_box1_u/2, inr, outr, draw_st_pos-r/2., draw_st_pos+r/2.);
  emc1_trap1->SetSegment(1);
  emc1_trap1->SetPosition(z, z, z);
  emc1_trap1->SetRotation(1., 0., 0., 0., 1., 0., 0., 0., 1.);
  emc1_trap1->DumpWithIncrement();

  Float_t cs, si;
  for (int i=1;i<number_emc1Chamber;i++){
    si=TMath::Sin(TMath::Pi()*(r*i)/180);
    cs=TMath::Cos(TMath::Pi()*(r*i)/180);
    
   emc1_trap1->SetSegment(i+1);
   emc1_trap1->SetPosition(z, z, z);
   emc1_trap1->SetRotation(cs, -si, z, si, cs, z, z, z, 1);
   emc1_trap1->DumpToFile();    
}

//*****************************************************
  Mpdshape* emc1_tr = new Mpdshape(fGeoFile, "emc1_tr", "emc1_trap1#1", "TRAP", "air");
  
  emc1_tr->Fill_TRAP(emc1_box1_d/2., emc1_box1_u/2, emc1_box1_h/2,emc1_box1_u/2, 0);
  emc1_tr->SetSegment(1);
  emc1_tr->SetPosition(z, inr+emc1_box1_h/2, z);
  emc1_tr->SetRotation(1., 0., 0., 0., 1., 0., 0., 0., 1.);
  emc1_tr->DumpToFile();

//*****************************************************
  Mpdshape* emc1_b_box = new Mpdshape(fGeoFile, "emc1_b_box", "emc1_tr#1", "BOX", "air");
  
  emc1_b_box->Fill_BOX(base/2, emc1_box1_h/2, emc1_box1_u/2);
  emc1_b_box->SetSegment(1);
  emc1_b_box->SetPosition(z, z, z);
  emc1_b_box->SetRotation(1., 0., 0., 0., 1., 0., 0., 0., 1.);
  emc1_b_box->DumpToFile(); 

//*****************************************************  
  Mpdshape* emc1_box = new Mpdshape(fGeoFile, "emc1_box", "emc1_b_box#1", "BOX", "air");
  
  emc1_box->Fill_BOX(base/2, emc1_box1_h/2, base/2);
  emc1_box->SetSegment(1); 
  emc1_box->SetPosition(z, z, z);
  emc1_box->SetRotation(1., 0., 0., 0., 1., 0., 0., 0., 1.);
  emc1_box->DumpWithIncrement();
  
  Float_t bbn;
  for(int i=1;i<3;i++){   
    if(i==1){
      bbn=base*i;
      emc1_box->SetSegment(i+1);
      emc1_box->SetPosition(z, z, bbn);
      emc1_box->DumpToFile(); 
    }
    if(i==2){
      bbn=base;
      emc1_box->SetSegment(i+1);
      emc1_box->SetPosition(z,z,-bbn);
      emc1_box->DumpToFile(); 
    }
  } 
 
//*****************************************************
  Float_t sumps;
  sumps=0.;
  for(int i=1;i<(st+1);i++){
    
    Mpdshape* emc1_box_pb = new Mpdshape(fGeoFile, "emc1_box_pb", "emc1_box#1", "BOX", "lead"); 
    
    emc1_box_pb->Fill_BOX(base/2,  box_pb_h/2, base/2);  
    emc1_box_pb->SetSegment(i);
    emc1_box_pb->SetMedia("lead");
    emc1_box_pb->SetPosition(z, -emc1_box1_h/2+box_pb_h/2+sumps, z);
    emc1_box_pb->SetRotation(1., 0., 0., 0., 1., 0., 0., 0., 1.);
    emc1_box_pb->DumpToFile();

//*****************************************************    
    Mpdshape* emc1_box_sc = new Mpdshape(fGeoFile, "emc1_box_sc", "emc1_box#1", "BOX", "FscScint");
    
    emc1_box_sc->Fill_BOX(base/2, box_sc_h/2, base/2);
    emc1_box_sc->SetSegment(i);
    emc1_box_sc->SetMedia("FscScint");
    emc1_box_sc->SetPosition(z, -emc1_box1_h/2+sumps+box_pb_h+box_sc_h/2, z);
    emc1_box_sc->SetRotation(1., 0., 0., 0., 1., 0., 0., 0., 1.); 
    emc1_box_sc->DumpToFile();
   
    sumps=(box_sc_h+box_pb_h)*i;
  }    
  
 
//*****************************************************
  Mpdshape* emc1_t = new Mpdshape(fGeoFile, "emc1_t", "emc1_tr#1", "TRAP", "air"); 
 
  emc1_t->Fill_TRAP(box_s_d/2, base-box_s_d/2, emc1_box1_h/2, emc1_box1_u/2, 1 );
  emc1_t->SetSegment(1);
  emc1_t->SetPosition(base/2+box_s_d/2, z, z);
  emc1_t->SetRotation(1., 0., 0., 0., 1., 0., 0., 0., 1.); 
  emc1_t->DumpWithIncrement();
  

  emc1_t->SetSegment(2);
  emc1_t->SetPosition(-base/2-box_s_d/2,z,z);
  emc1_t->SetRotation(-1., 0., 0., 0., 1., 0., 0., 0., 1.);
  emc1_t->DumpToFile();   

//*****************************************************
  Mpdshape* emc1_tcl = new Mpdshape(fGeoFile, "emc1_tcl", "emc1_t#1", "TRAP", "air");
 
  emc1_tcl->Fill_TRAP(box_s_d/2, base-box_s_d/2, emc1_box1_h/2, base/2, 1);
  emc1_tcl->SetSegment(1);
  emc1_tcl->SetPosition(z, z, z);
  emc1_tcl->SetRotation(1., 0., 0., 0., 1., 0., 0., 0., 1.);
  emc1_tcl->DumpWithIncrement();
 
  Float_t bbn=base;
  for(int i=1;i<3;i++){  
     if(i==1){
       emc1_tcl->SetSegment(i+1);    
       emc1_tcl->SetPosition(z,z,bbn);
       emc1_tcl->DumpToFile();
     }
     if(i==2){
       emc1_tcl->SetSegment(i+1);    
       emc1_tcl->SetPosition(z,z,-bbn);
       emc1_tcl->DumpToFile(); 
     }
   }
 
//*****************************************************
  std::ostringstream name;
  name.str("");
  
  Double_t x, x1, y, p, xx,xxx, xx1,xx2, xr,xrr;
  x=0.;
  xx=0.;
  y=apb;
  p=asc; 
  xr=0.;
  xrr=0.;
  for(int i=1;i<(((base-box_s_d)/tan(q))/box_sum);i++){ 
    if(ad<(base-asc-apb)){ 

	name.str("");
	name << "emc1_cl_pb" << i;
	
	Mpdshape*  emc1_cl_pb = new Mpdshape(fGeoFile, (name.str()).c_str(), "emc1_tcl#1", "TRAP", "lead");
	
	emc1_cl_pb->Fill_TRAP((box_s_d+x)/2,((box_s_d+y+x)-(box_s_d+x)/2), box_pb_h/2., base/2., 1);
	emc1_cl_pb->SetSegment(1);
	emc1_cl_pb->SetMedia("lead");
	emc1_cl_pb->SetPosition((box_s_d+x)/2-box_s_d/2, -emc1_box1_h/2.+box_pb_h/2.+x1, z);
	emc1_cl_pb->SetRotation(1., 0., 0., 0., 1., 0., 0., 0., 1.);
	emc1_cl_pb->DumpToFile();
	
	name.str("");
	name << "emc1_cl_sc" << i;
	Mpdshape*  emc1_cl_sc = new Mpdshape(fGeoFile,(name.str()).c_str(), "emc1_tcl#1", "TRAP", "FscScint");
	
	emc1_cl_sc->Fill_TRAP((box_s_d+x+y)/2,((box_s_d+y+x+p)-(box_s_d+x+y)/2), box_sc_h/2., base/2.,1); 
	emc1_cl_sc->SetSegment(1);
	emc1_cl_sc->SetMedia("FscScint");
	emc1_cl_sc->SetPosition((box_s_d+x+y)/2-box_s_d/2, -emc1_box1_h/2.+box_pb_h+x1+box_sc_h/2., z);
	emc1_cl_sc->SetRotation(1., 0., 0., 0., 1., 0., 0., 0., 1.);
	emc1_cl_sc->DumpToFile(); 
	
	x=x+y+p;
	x1=(box_sc_h+box_pb_h)*i; 
	ad=box_s_d+x;
	
    }
  }
}

#if 0

#endif
















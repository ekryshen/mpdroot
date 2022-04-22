#include "EmcGeom.h"
#include <iostream>
#include <string>
#include <sstream>

ClassImp(EmcGeom);

using namespace emc;

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
 
//***************************_______All Barrel ECal____**************  
 tube->Fill_TUBE(emc1Chamber_z_l/2, inr, outro);
 tube->SetSegment(0);
 tube->SetPosition(0., 0., 0.);
 tube->SetRotation(1., 0., 0., 0., 1., 0., 0., 0., 1.);

 tube->DumpToFile();
 BuildSensVolume();
  
  
 delete tube;
}
//_____________________________________________________________

void EmcGeom::BuildSensVolume() { 

//***************************_______1/2 Barrel____**************  
  Mpdshape* emc1ChH= new Mpdshape(fGeoFile, "emc1ChH", "emc1Chamber1#1", "TUBE", "air");
  emc1ChH->Fill_TUBE(emc1Chamber_z/2, inr, outro);
  emc1ChH->SetSegment(1);
  emc1ChH->SetPosition(0. ,0. , -emc1Chamber_z_th);
  emc1ChH->SetRotation(1., 0., 0., 0., 1., 0., 0., 0., 1.);
  emc1ChH->DumpWithIncrement();

  emc1ChH->SetSegment(2);
  emc1ChH->SetPosition(0. ,0. , emc1Chamber_z_th);
  emc1ChH->SetRotation(-1., 0., 0., 0., 1., 0., 0., 0., 1.);
  emc1ChH->DumpWithIncrement();

//***************************_______Sector____**************
  Mpdshape* emc1Sector = new Mpdshape(fGeoFile, "emc1Sector", "emc1ChH#1", "TUBS", "air");
  emc1Sector->Fill_TUBS(-emc1Chamber_z/2, emc1Chamber_z/2, inr, outro, draw_st_pos-angle*4/2, draw_st_pos+angle*4/2);
  emc1Sector->SetSegment(1);
  emc1Sector->SetPosition(z, z, z);
  emc1Sector->SetRotation(1., 0., 0., 0., 1., 0., 0., 0., 1.);
  emc1Sector->DumpWithIncrement();
  
  for(int i=1; i<NSector; i++){
    si=TMath::Sin(TMath::Pi()*(angle*4*i)/180);
    cs=TMath::Cos(TMath::Pi()*(angle*4*i)/180);
    
    emc1Sector->SetSegment(i+1);
    emc1Sector->SetPosition(z, z, z);
    emc1Sector->SetRotation(cs, -si, z, si, cs, z, z, z, 1); 
    emc1Sector->DumpWithIncrement();
  }  
//***************************_______Parameter by the position in Tube____**************
  angleSteel = (180*(TMath::ATan((stThM/2)/(inr+box_pad)))/TMath::Pi());       //half angle in degree
  angleModule = (180*(TMath::ATan((emc1_box1_d/2)/(inr+box_pad)))/TMath::Pi());//half angle in degree
 
//***************************_______Tube of Steel____**************
 Mpdshape* emc1TubeSt = new Mpdshape(fGeoFile, "emc1TubeSt", "emc1Sector#1", "TUBS", "steel");
  emc1TubeSt->Fill_TUBS(-emc1Chamber_z/2, emc1Chamber_z/2, inr+box_pad, outro, draw_st_pos-angle*4/2, draw_st_pos-angle*4/2+angleSteel*2);
  emc1TubeSt->SetSegment(1);       
  emc1TubeSt->SetPosition(z, z, z);
  emc1TubeSt->SetRotation(1., 0., 0., 0., 1., 0., 0., 0., 1.);
  emc1TubeSt->DumpWithIncrement(); 
   
  for (int i=1;i<5;i++){
    si=TMath::Sin(TMath::Pi()*((angleModule*2+angleSteel*2)*i)/180);
    cs=TMath::Cos(TMath::Pi()*((angleModule*2+angleSteel*2)*i)/180);
    emc1TubeSt->SetSegment(i+1);
    emc1TubeSt->SetPosition(z, z, z);
    emc1TubeSt->SetRotation(cs, -si, z, si, cs, z, z, z, 1);
    emc1TubeSt->DumpWithIncrement();  
  }
  
//***************************_______Tube of Module____**************
 Mpdshape* emc1Tube = new Mpdshape(fGeoFile, "emc1Tube", "emc1Sector#1", "TUBS", "air");
  emc1Tube->Fill_TUBS(-emc1_box1_z/2, emc1_box1_z/2, inr, outro, draw_st_pos-angle*4/2+angleSteel*2, draw_st_pos-angle*4/2+angleSteel*2+angleModule*2);
  emc1Tube->SetSegment(1);       
  emc1Tube->SetPosition(z, z, -emc1Chamber_z/2 + emc1_box1_z/2);
  emc1Tube->SetRotation(1., 0., 0., 0., 1., 0., 0., 0., 1.);
  emc1Tube->DumpWithIncrement();
 
//***************************_______Fill Sector with Tube in Phi and Theta plane
  for (int i=0;i<4;i++){
    for(int j=0; j<NmodZ; j++){
        si=TMath::Sin(TMath::Pi()*((angleModule*2+angleSteel*2)*i)/180);
        cs=TMath::Cos(TMath::Pi()*((angleModule*2+angleSteel*2)*i)/180);
        sch++;
        if(sch>1){
            emc1Tube->SetSegment(sch+1);
            emc1Tube->SetPosition(z, z, -emc1Chamber_z/2 + emc1_box1_z/2+emc1_box1_z*j);
            emc1Tube->SetRotation(cs, -si, z, si, cs, z, z, z, 1);
            emc1Tube->DumpWithIncrement();
        }
    }
  }  

//***************************_______Plastic pad(base) to the module____**************
  Mpdshape* emc1ModulePad = new Mpdshape(fGeoFile, "emc1ModulePad", "emc1Tube#1", "TRAP", "polypropylene"); 
  emc1ModulePad->Fill_TRAP(emc1_box1_d/2, emc1_box1_d/2, emc1_box1_d/2, emc1_box1_d/2, box_pad,  emc1_box1_z/2, 0);

  angleModulePad = (180*(TMath::ATan((emc1_box1_d/2)/inr))/TMath::Pi());//half angle in degree
  dA = (4*angle) - angleModulePad*2*4;
  dA = dA/5; 

  dAlf = (4*angle)/2 - angleModulePad - dA; 
  Alf = TMath::Pi()*(dAlf)/180;
  dx =   inr*(TMath::Tan(Alf));
  dx1 = (inr+box_pad)*(TMath::Tan(Alf));  

  si=TMath::Sin(TMath::Pi()*(360-dAlf)/180);
  cs=TMath::Cos(TMath::Pi()*(360-dAlf)/180);
  cout<<angleModulePad<<" "<<dA<<" "<<dAlf<<" "<<Alf<<" "<<dx<<endl;

  emc1ModulePad->SetSegment(1);
  emc1ModulePad->SetPosition(dx, inr, z);
  emc1ModulePad->SetRotation(cs, -si, z, si, cs, z, z, z, 1); 
  emc1ModulePad->DumpWithIncrement();

//***************************_______Module____**************
  Mpdshape* emc1Module = new Mpdshape(fGeoFile, "emc1Module", "emc1Tube#1", "TRAP", "air"); 
  emc1Module->Fill_TRAP(emc1_box1_d/2, emc1_box1_z/2, 0, 0, lenghtM,  emc1_box1_z/2, 0);
  emc1Module->SetSegment(1);
  emc1Module->SetPosition(dx1, inr+box_pad, z);
  emc1Module->SetRotation(cs, -si, z, si, cs, z, z, z, 1); 
  emc1Module->DumpWithIncrement();

//***************************_______Box in Module____*********
 Mpdshape* emc1_box = new Mpdshape(fGeoFile, "emc1_box", "emc1Module#1", "TRAP", "air"); 
  emc1_box->Fill_TRAP(base/2, 0, 0, 0, lenghtM,  base/2, 3);
  emc1_box->SetSegment(1);
  emc1_box->SetPosition(z, z, -base);
  emc1_box->SetRotation(1., 0., 0., 0., 1., 0., 0., 0., 1.); 
  emc1_box->DumpWithIncrement();
    
  Float_t bbn;
  for(int i=1;i<3;i++){
      bbn=base*i;
      emc1_box->SetSegment(i+1);
      emc1_box->SetPosition(z,z,-base+base*i);
      emc1_box->SetRotation(1., 0., 0., 0., 1., 0., 0., 0., 1.); 
      emc1_box->DumpWithIncrement(); 
  } 

//***************************_______Fill Box____*********
 Mpdshape* emc1_box_pb = new Mpdshape(f, "emc1_box_pb", "emc1_box#1", "TRAP", "lead"); 
  emc1_box_pb->Fill_TRAP(base/2, 0, 0, 0, box_pb_h,  base/2, 3);
  emc1_box_pb->SetSegment(1);   
  emc1_box_pb->SetPosition(z,z,z);
  emc1_box_pb->SetRotation(1., 0., 0., 0., 1., 0., 0., 0., 1.); 
  emc1_box_pb->DumpToFile();
  
 Mpdshape* emc1_box_sc = new Mpdshape(f, "emc1_box_sc", "emc1_box#1", "TRAP", "FscScint");
  emc1_box_sc->Fill_TRAP(base/2, 0, 0, 0, box_sc_h,  base/2, 3);
  emc1_box_sc->SetSegment(1);
  emc1_box_sc->SetMedia("FscScint");
  emc1_box_sc->SetPosition(z,box_pb_h,z);
  emc1_box_sc->SetRotation(1., 0., 0., 0., 1., 0., 0., 0., 1.); 
  emc1_box_sc->DumpToFile();
      
  Float_t sumps = 0;
  for(int il=1;il<NScPb; il++){
    sumps=(box_sc_h+box_pb_h)*il;
    emc1_box_pb->SetSegment(il+1); 
    emc1_box_pb->SetPosition(z,sumps,z);
    emc1_box_pb->SetRotation(1., 0., 0., 0., 1., 0., 0., 0., 1.); 
    emc1_box_pb->DumpToFile();  
      
    emc1_box_sc->SetSegment(il+1);
    emc1_box_sc->SetPosition(z,box_pb_h+sumps,z);
    emc1_box_sc->SetRotation(1., 0., 0., 0., 1., 0., 0., 0., 1.); 
    emc1_box_sc->DumpToFile();   
  }

//***************************_______TRAP in Module____*********  
 Mpdshape* emc1_bt_box = new Mpdshape(fGeoFile, "emc1_bt_box", "emc1Module#1", "TRAP", "air"); 
  emc1_bt_box->Fill_TRAP((emc1_box1_d-base)/4, (base-(emc1_box1_d-base)/4), 0, 0, lenghtM,  base/2, 1);
  Int_t ssr = 0; 
  for(int it = 0; it<2; it++){
    for(int jt=0; jt<3; jt++){      //po "z"
      ssr+=1;
      if(it==0){
	    emc1_bt_box->SetSegment(ssr);
	    emc1_bt_box->SetPosition(base/2+(emc1_box1_d-base)/4,z,-base+base*jt);
	    emc1_bt_box->SetRotation(1., 0., 0., 0., 1., 0., 0., 0., 1.); 
	    emc1_bt_box->DumpWithIncrement();   
	  }
	  if(it==1){
	    emc1_bt_box->SetSegment(ssr);
	    emc1_bt_box->SetPosition(-(base/2+(emc1_box1_d-base)/4),z,-base+base*jt);
	    emc1_bt_box->SetRotation(-1., 0., 0., 0., 1., 0., 0., 0., -1.); 
	    emc1_bt_box->DumpWithIncrement();   
	  }	  
    }
  }

//***************************_______Fill TRAP in Module____*********    
  std::ostringstream name, nameP;
  Double_t x, y, p, x1, t;

  x=0.;
  x1 = 0;
  y=2*apb;
  p=2*asc; 
  cout<<angleTrap<<" "<<apb<<" "<<asc<<endl;

  for(int i=1;i<NScPb+1;i++){             //(((base-box_s_d)/tan(q))/box_sum);i++){ 
	nameP.str("");
	nameP << "emc1_cl_pb" << i;
	Mpdshape*  emc1_cl_pb = new Mpdshape(fGeoFile,(nameP.str().c_str()), "emc1_bt_box#1", "TRAP", "lead");	
	emc1_cl_pb->Fill_TRAP((emc1_box1_d-base)/4+x/2,(emc1_box1_d-base)/4+apb+x/2, 0, 0, box_pb_h, base/2.,1); 
	emc1_cl_pb->SetSegment(1);
	emc1_cl_pb->SetPosition((x)/2, x1, z);
	emc1_cl_pb->SetRotation(1., 0., 0., 0., 1., 0., 0., 0., 1.);
	emc1_cl_pb->DumpToFile(); 
 
	name.str("");
	name << "emc1_cl_sc" << i;
	Mpdshape*  emc1_cl_sc = new Mpdshape(fGeoFile,(name.str().c_str()), "emc1_bt_box#1", "TRAP", "FscScint");
	
	emc1_cl_sc->Fill_TRAP((emc1_box1_d-base)/4+(x+apb)/2,(emc1_box1_d-base)/4+(apb+x)/2+asc, 0, 0, box_sc_h, base/2.,1); 
	emc1_cl_sc->SetSegment(1);
	emc1_cl_sc->SetPosition((x+apb)/2, x1+box_pb_h, z);
	emc1_cl_sc->SetRotation(1., 0., 0., 0., 1., 0., 0., 0., 1.);
	emc1_cl_sc->DumpToFile(); 

	x=x+apb+asc;
	x1=(box_sc_h+box_pb_h)*i; 

  }
}

#if 0

#endif








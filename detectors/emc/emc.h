#include "TMath.h"

#ifndef EMC_H
#define EMC_H 1

namespace emc {
Double_t z = 0.0;
Double_t draw_st_pos = 90.;

UShort_t nModInSuperModByZ = 3;
UShort_t nModInSuperModByPhi = 3;
UShort_t nSuperModInRow = 4;

Double_t h = 407.8; //lenght of all module (Pb+Sc)*n + plastics
Double_t box_pb_h = 0.3; //thickness of Pb
Double_t box_sc_h = 1.5; //thickness of Sc
Double_t box_pad = 10; //thickness of base of the module
Double_t base = 40.; //width of Channel
Int_t NScPb = 221; //number of layers Sc+Pb
Int_t NmodZ = 23; //modules in one sector
Int_t NSector = 28; //number of sector in 1/2 Barrel
Int_t NmodPHY = 112; //number of modules in Phi plane
Double_t emc1_box1_d = 98.; //size of down side of Trap
Double_t stThM = 0.5; //thickness of steel ~between modules - 0.5mm; ~between sector 1mm


Double_t angleSector = 360. / NSector; //rotation angel in Phi plane for sector
Double_t angleSuperModule = angleSector / nSuperModInRow; //rotation angel in Phi plane for super-module
Double_t angleModule = angleSuperModule / nModInSuperModByPhi; //rotation angel in Phi plane for module
Double_t emc1_box1_z = base * nModInSuperModByZ; //120 mm
Double_t emc1Chamber_z = NmodZ*emc1_box1_z; //lenght "z" in 1/2 Barrel
Double_t emc1Chamber_z_l = emc1Chamber_z * 2; //lenght "z" in all Barrel
Double_t emc1Chamber_z_th = emc1Chamber_z / 2; //"z" start position of 1/2 Barrel 
Double_t lenghtM = NScPb*(box_pb_h + box_sc_h); //real lenght of module

Double_t inr = 1748; //Inner rarius
Double_t outr = inr + h; //Outr radius -> module + fron's plastic and end's plastic
Double_t outro = inr + 464.4; //modul+base+space for electronics
//Double_t angleModule; //rotation module in sector
Double_t angleSteel; //rotation steel in sector
Double_t angleModulePad; //cut angle in module
Double_t angleTrap = TMath::ATan(((emc1_box1_z - emc1_box1_d) / 2) / lenghtM);

Double_t apb = box_pb_h*TMath::Tan(angleTrap);
Double_t asc = box_sc_h*TMath::Tan(angleTrap);

Float_t cs, si, dx, dx1, dxM, dxSt, dAlf, Alf, dAlfM, AlfM, dAlfSt, AlfSt, dA;
Int_t sch;

}
#endif

//---------------------------
//M.Golubeva (INR RAS), February 2016
//---------------------------

#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoMedium.h"
#include "TFile.h"
#include "TMath.h"

R__ADD_INCLUDE_PATH($VMCWORKDIR)
//#include "macro/mpd/mpdloadlibs.C"
#include "../mpdloadlibs.C"
//#include "ZDC_geom_par.h"

//Detector's position
const Double_t ZDC_Xpos = 0.0;
const Double_t ZDC_Ypos = 0.0;
//const Double_t ZDC_Zpos = 319.5;//cm to ZDC surface
//const Double_t ZDC_center_Zpos = 381.9;//3.819 m (319.5+62.4);//cm 10 sect
const Double_t ZDC_center_Zpos = 363.54;//3.6354 m (319.5+44.04);//cm 7 sect

//Detector's construct parameters ===================

//Global dimensions for ZDC in cm
Double_t Zdc_X_size = 52.5; //105./2.
Double_t Zdc_Y_size = 52.5; //105./2.
//Double_t Zdc_Z_size = 62.4; //10 sect
Double_t Zdc_Z_size = 44.04; //7 sect

Double_t Zdc_mother_inner_radius = 7.5; //cm
Double_t Zdc_mother_outer_radius = 64.518; //cm //sqrt(52.5*52.5+37.5*37.5))

//Modules
Double_t Zdc_NColumns = 7; 
Double_t Zdc_NRows = 7; 
Double_t Zdc_module_X_size = 7.5; //15./2.
Double_t Zdc_module_Y_size = 7.5; //15./2.
//Double_t Zdc_module_Z_size = 62.4; //10 sect
Double_t Zdc_module_Z_size = 44.04; //7 sect

//Fe slices
Double_t Zdc_Fe_X_size = 7.4; //(15.-2.*0.1)/2.
Double_t Zdc_Fe_Y_size = 7.35; //(15.-2.*0.15)/2.
Double_t Zdc_Fe_Z_size = 1.; //2./2.

//Pb slices
//Double_t Zdc_NPb_in_module = 59;//10sect
Double_t Zdc_NPb_in_module = 41;//7sect
Double_t Zdc_Pb_X_size = 7.4; //(15.-2.*0.1)/2.
Double_t Zdc_Pb_Y_size = 7.35; //(15.-2.*0.15)/2.
Double_t Zdc_Pb_Z_size = 0.8; //1.6/2.

//Scin. slices
//Double_t Zdc_NSc_in_module = 1;//put into Tyvec
Double_t Zdc_Sc_X_size = 7.38; //(15.-2.*0.1-2.*0.02)/2.
Double_t Zdc_Sc_Y_size = 7.33; //(15.-2.*0.15-2.*0.02)/2.
Double_t Zdc_Sc_Z_size = 0.2; //0.4/2.

//Tyvec slices
//Double_t Zdc_NTyvec_in_module = 60;//10sect
Double_t Zdc_NTyvec_in_module = 42;//7sect
Double_t Zdc_Tyvec_X_size = 7.4; //(15.-2.*0.1)/2.
Double_t Zdc_Tyvec_Y_size = 7.35; //(15.-2.*0.15)/2.
Double_t Zdc_Tyvec_Z_size = 0.22; //0.44/2.

//Slot slices (divided by 2)
Double_t Zdc_Y_slot_position = 0.75; 

Double_t Zdc_VRFM_X_size = 0.1; //Fe
Double_t Zdc_VRFM_Y_size = 6.1; 
Double_t Zdc_VRFM_Z_size = 1.; 

Double_t Zdc_VRPM_X_size = 0.1; //Pb
Double_t Zdc_VRPM_Y_size = 6.1; 
Double_t Zdc_VRPM_Z_size = 0.8; 

Double_t Zdc_VRYM_X_size = 0.01; //Tyvec_Y
Double_t Zdc_VRYM_Y_size = 6.1; 
Double_t Zdc_VRYM_Z_size = 0.22; 

Double_t Zdc_VRAM_X_size = 0.09; //Tyvec_A 
Double_t Zdc_VRAM_Y_size = 6.1; 
Double_t Zdc_VRAM_Z_size = 0.01; 

Double_t Zdc_VRSM_X_size = 0.09; //Scint 
Double_t Zdc_VRSM_Y_size = 6.1; 
Double_t Zdc_VRSM_Z_size = 0.2; 

//Hole
Double_t Zdc_Hole_inner_radius = 0.; 
Double_t Zdc_Hole_outer_radius = 5.; 
Double_t Zdc_Hole_Tyvec_Z_size = 0.01; 


//media ============================================
TGeoMedium *pMedVacuum = 0;
TGeoMedium *pMedAir = 0;
TGeoMedium *pMedIron = 0;
TGeoMedium *pMedLead = 0;
TGeoMedium *pMedScint = 0;
TGeoMedium *pMedTyvec = 0;
TGeoMedium *pMedSlot = 0;

class FairGeoMedia;
class FairGeoBuilder;

void FillModule(TGeoVolume *mother_volume,TGeoVolume *mother_volume_hole);
void DefineRequiredMedia(FairGeoMedia* geoMedia, FairGeoBuilder* geoBuild) {
    //vacuum medium
    FairGeoMedium* mVacuum = geoMedia->getMedium("vacuum");
    if ( ! mVacuum ) Fatal("Main", "FairMedium vacuum not found");
    geoBuild->createMedium(mVacuum);
    pMedVacuum = gGeoManager->GetMedium("vacuum");
    if ( ! pMedVacuum ) Fatal("Main", "Medium vacuum not found");

    //air medium
    FairGeoMedium* mAir = geoMedia->getMedium("air");
    if ( ! mAir ) Fatal("Main", "FairMedium air not found");
    geoBuild->createMedium(mAir);
    pMedAir = gGeoManager->GetMedium("air");
    if ( ! pMedAir ) Fatal("Main", "Medium air not found");

    //iron medium
    FairGeoMedium* mIron = geoMedia->getMedium("iron");
    if ( ! mIron ) Fatal("Main", "FairMedium iron not found");
    geoBuild->createMedium(mIron);
    pMedIron = gGeoManager->GetMedium("iron");
    if ( ! pMedIron ) Fatal("Main", "Medium iron not found");

    //tyvec medium
    FairGeoMedium* mTyvec = geoMedia->getMedium("tyvec");
    if ( ! mTyvec ) Fatal("Main", "FairMedium tyvec not found");
    geoBuild->createMedium(mTyvec);
    pMedTyvec = gGeoManager->GetMedium("tyvec");
    if ( ! pMedTyvec ) Fatal("Main", "Medium tyvec not found");

    //lead medium
    FairGeoMedium* mLead = geoMedia->getMedium("lead");
    if ( ! mLead ) Fatal("Main", "FairMedium lead not found");
    geoBuild->createMedium(mLead);
    pMedLead = gGeoManager->GetMedium("lead");
    if ( ! pMedLead ) Fatal("Main", "Medium lead not found");
    /*
    // medium polystyrene (Scint.)
    FairGeoMedium* mPolystyrene = geoMedia->getMedium("polystyrene");
    if ( ! mPolystyrene ) Fatal("Main", "FairMedium polystyrene not found");
    geoBuild->createMedium(mPolystyrene);
    pMedScint = gGeoManager->GetMedium("polystyrene");
    if ( ! pMedScint ) Fatal("Main", "Medium polystyrene not found");
    */
    // medium FscScint (Scint.)
    FairGeoMedium* mPolystyrene = geoMedia->getMedium("FscScint");
    if ( ! mPolystyrene ) Fatal("Main", "FairMedium FscScint not found");
    geoBuild->createMedium(mPolystyrene);
    pMedScint = gGeoManager->GetMedium("FscScint");
    if ( ! pMedScint ) Fatal("Main", "Medium FscScint not found");

    // medium plastic (slot)
    FairGeoMedium* mPlastic = geoMedia->getMedium("plastic");
    if ( ! mPlastic ) Fatal("Main", "FairMedium plastic not found");
    geoBuild->createMedium(mPlastic);
    pMedSlot = gGeoManager->GetMedium("plastic");
    if ( ! pMedSlot ) Fatal("Main", "Medium plastic not found");
}
//DefineRequiredMedia

void create_rootgeom_zdc_oldnames_7sect_v1_no_overlaps_w_pipe_magnet() {

    // Load necessary libraries
    //gROOT->LoadMacro("$VMCWORKDIR/macro/mpd/mpdloadlibs.C");
    mpdloadlibs(); // load libraries

    // ----  set working directory  --------------------------------------------
    TString gPath = gSystem->Getenv("VMCWORKDIR");

    // -------   Geometry file name (output)   ---------------------------------
    const TString geoDetectorName = "zdc_oldnames_7sect";
    const TString geoDetectorVersion = "v1";
    TString geoFileName = gPath + "/geometry/" + geoDetectorName + "_"+ geoDetectorVersion + "_no_overlaps_w_pipe_magnet.root";
    // -------------------------------------------------------------------------

    // ----  global geometry parameters  ---------------------------------------
    FairGeoLoader*    geoLoad = new FairGeoLoader("TGeo","FairGeoLoader");
    FairGeoInterface* geoFace = geoLoad->getGeoInterface();

    // -------   Load media from media file   ----------------------------------
    TString medFile = gPath + "/geometry/media.geo";
    geoFace->setMediaFile(medFile);
    geoFace->readMedia();
    FairGeoMedia*   geoMedia = geoFace->getMedia();
    FairGeoBuilder* geoBuild = geoLoad->getGeoBuilder();
    //geoMedia->list();
    //geoMedia->print();
    DefineRequiredMedia(geoMedia, geoBuild);
    // -------------------------------------------------------------------------

    // --------------   Create geometry and global top volume  ------------------------
    gGeoManager = (TGeoManager*)gROOT->FindObject("FAIRGeom");
    gGeoManager->SetName(geoDetectorName + "_geom");
    TGeoVolume* top = new TGeoVolumeAssembly("TOP");
    top->SetMedium(pMedAir);
    gGeoManager->SetTopVolume(top);
    //gGeoMan->SetTopVisible(1);
    // --------------------------------------------------------------------------

    // Define ZDC Geometry
    TString zdc_common_name = "ZDC_common";
    TGeoVolume *ZDC_common = new TGeoVolumeAssembly(zdc_common_name);
    ZDC_common->SetMedium(pMedAir);
    //ZDC_common->SetMedium(pMedVacuum);

    TString zdc_name = "zdc01";
    TGeoTube *ZDCS = new TGeoTube(zdc_name, Zdc_mother_inner_radius, Zdc_mother_outer_radius, Zdc_Z_size);
    TGeoVolume *ZDCV = new TGeoVolume(zdc_name, ZDCS);
    ZDCV->SetMedium(pMedAir);

    TString zdcModule_name = "zdc01m";

    TGeoBBox *zdcModuleS = new TGeoBBox(zdcModule_name, Zdc_module_X_size, Zdc_module_Y_size, Zdc_module_Z_size);
    TGeoVolume *zdcModuleV = new TGeoVolume(zdcModule_name, zdcModuleS);

    zdcModuleV->SetMedium(pMedIron);

    //zdcModuleV->SetLineColor(TColor::GetColor("#fefe22"));
    ////zdcModuleV->SetVisibility(false);
    //zdcModuleV->SetTransparency(95);

    TString zdcModule_hole_name = "module_hole";

    TGeoBBox *zdcModuleSH = new TGeoBBox(zdcModule_hole_name, Zdc_module_X_size, Zdc_module_Y_size, Zdc_module_Z_size);
    TGeoVolume *zdcModuleVH = new TGeoVolume(zdcModule_hole_name, zdcModuleSH);

    zdcModuleVH->SetMedium(pMedIron);

    Double_t xxxxNICA[91], yyyyNICA[91];
    //Fill ZDC with modules
    Float_t xCur = Zdc_X_size + Zdc_module_X_size;//60 cm
    Float_t yCur = Zdc_Y_size + Zdc_module_Y_size;//60 cm

    Int_t iMod=-1, modNb=-1;
    for(Int_t iy=0; iy<Zdc_NRows; iy++) {
      yCur=yCur - 15.;
      for(Int_t ix=0; ix<Zdc_NColumns; ix++) {
	xCur=xCur - 15.;
	
	TGeoTranslation *zdcModuleV_position = new TGeoTranslation(xCur,yCur,0.);
	iMod=iy*7 + ix;
	
	if((iMod+1)>=2 && (iMod+1)<=6) {
	  modNb=iMod+1-1;  //mods 1-5
	  ZDCV->AddNode(zdcModuleV,modNb,zdcModuleV_position);
	  xxxxNICA[modNb-1] = xCur; 
	  yyyyNICA[modNb-1] = yCur; 
	  cout <<"1-5 iMod,modNb,xxxx,yyyy " <<iMod+1 <<" " <<modNb <<" " <<xCur <<" " <<yCur <<endl;
	}
	if((iMod+1)>=8 && (iMod+1)<=24) {
	  modNb=iMod+1-2;  //mods 6-22
	  ZDCV->AddNode(zdcModuleV,modNb,zdcModuleV_position);
	  xxxxNICA[modNb-1] = xCur; 
	  yyyyNICA[modNb-1] = yCur; 
	  cout <<"6-22 iMod,modNb,xxxx,yyyy " <<iMod+1 <<" " <<modNb <<" " <<xCur <<" " <<yCur <<endl;
	}

	//hole in place of mods 23 and 68

	if((iMod+1)>=26 && (iMod+1)<=42) {
	  modNb=iMod+1-2;  //mods 24-40
	  ZDCV->AddNode(zdcModuleV,modNb,zdcModuleV_position);
	  xxxxNICA[modNb-1] = xCur; 
	  yyyyNICA[modNb-1] = yCur; 
	  cout <<"24-40 iMod,modNb,xxxx,yyyy " <<iMod+1 <<" " <<modNb <<" " <<xCur <<" " <<yCur <<endl;
	}
	if((iMod+1)>=44 && (iMod+1)<=48) {
	  modNb=iMod+1-3;  //mods 41-45
	  ZDCV->AddNode(zdcModuleV,modNb,zdcModuleV_position);
	  xxxxNICA[modNb-1] = xCur; 
	  yyyyNICA[modNb-1] = yCur; 
	  cout <<"41-45 iMod,modNb,xxxx,yyyy " <<iMod+1 <<" " <<modNb <<" " <<xCur <<" " <<yCur <<endl;
	}
	
	if(ix==6) xCur=60.;
	
      }//for(Int_t ix==0; ix<7; ix++)
    }//for(Int_t iy==0; iy<7; iy++)

    Double_t xxxx[91], yyyy[91];
    for(Int_t im=0; im<90; im++) {
      if(im<45) {
	xxxx[im] = xxxxNICA[im];
	yyyy[im] = yyyyNICA[im];
      }
      else {
	xxxx[im] = -xxxxNICA[im-45];
	yyyy[im] = yyyyNICA[im-45];
      }
    }
    xxxx[22] = 0; yyyy[22] = 0; 
    xxxx[67] = 0; yyyy[67] = 0; 
    for(Int_t im=0; im<90; im++) { 
      if(im<9) cout <<im+1 << setw(11) <<xxxx[im] << setw(10) <<yyyy[im] <<endl;
      else cout <<im+1 << setw(10) <<xxxx[im] << setw(10) <<yyyy[im] <<endl;
    }
    //return 1;

    //fill module and module with hole with Fe-Sc-Lead slices
    FillModule(zdcModuleV,zdcModuleVH);

    //Adding the zdc mother volume to the global TOP Volume
    TGeoTranslation *zdcZDC = new TGeoTranslation("zdcZDC",ZDC_Xpos, ZDC_Ypos, ZDC_center_Zpos);
    ZDC_common->AddNode(ZDCV, 1, zdcZDC);//ZDC right

    TGeoRotation *rotZDC = new TGeoRotation();
    rotZDC->RotateY(180.);
    TGeoCombiTrans *zdcZDC2 = new TGeoCombiTrans(ZDC_Xpos, ZDC_Ypos, -ZDC_center_Zpos, rotZDC);
    ZDC_common->AddNode(ZDCV, 2, zdcZDC2);//ZDC left

    top->AddNode(ZDC_common, 0);

    top->SetVisContainers(kTRUE);

    // ---------------   Finish   ----------------------------------------------
    gGeoManager->CloseGeometry();
    gGeoManager->CheckOverlaps(0.001);
    gGeoManager->PrintOverlaps();

    gGeoManager->Test();

    TFile* geoFile = new TFile(geoFileName, "RECREATE");
    top->Write();
    geoFile->Close();

    //top->Draw("ogl");
    //TGLViewer *v = (TGLViewer*)gPad->GetViewer3D();
    //v->SetStyle(TGLRnrCtx::kOutline);

}//create_rootgeom_zdc_v1

void FillModule(TGeoVolume *mother_volume,TGeoVolume *mother_volume_hole) {

    TString Fe_name = "Fe";
    TString slot_Fe_name = "Fe_slot";
    TString Pb_name = "zdc01a";
    TString slot_Pb_name = "Pb_slot";
    TString SensitiveVolume_name = "zdc01s";
    TString slot_Sc_name = "Sc_slot";
    TString Tyvec_name = "Tyvec";
    TString slot_Tyvec_Y_name = "Tyvec_slot_Y";
    TString slot_Tyvec_A_name = "Tyvec_slot_A";

    //Shapes
    TGeoShape *FeS = new TGeoBBox(Fe_name, Zdc_Fe_X_size, Zdc_Fe_Y_size, Zdc_Fe_Z_size);
    TGeoVolume *FeV = new TGeoVolume(Fe_name, FeS);
    TGeoMedium *Fe_medium = pMedIron; //set medium
    if(Fe_medium) FeV->SetMedium(Fe_medium);
    else Fatal("Main", "Invalid medium for Fe!");

    TGeoShape *PbS = new TGeoBBox(Pb_name, Zdc_Pb_X_size, Zdc_Pb_Y_size, Zdc_Pb_Z_size);
    TGeoVolume *PbV = new TGeoVolume(Pb_name, PbS);
    TGeoMedium *Pb_medium = pMedLead; //set medium
    if(Pb_medium) PbV->SetMedium(Pb_medium);
    else Fatal("Main", "Invalid medium for Pb!");

    TGeoShape *zdcSensitiveVolumeS = new TGeoBBox(SensitiveVolume_name, Zdc_Sc_X_size, Zdc_Sc_Y_size, Zdc_Sc_Z_size);
    TGeoVolume *zdcSensitiveVolumeV = new TGeoVolume(SensitiveVolume_name, zdcSensitiveVolumeS);
    TGeoMedium *Sc_medium = pMedScint; //set medium
    if(Sc_medium) zdcSensitiveVolumeV->SetMedium(Sc_medium);
    else Fatal("Main", "Invalid medium for Sc!");


    TGeoShape *TyvecS = new TGeoBBox(Tyvec_name, Zdc_Tyvec_X_size, Zdc_Tyvec_Y_size, Zdc_Tyvec_Z_size);
    TGeoVolume *TyvecV = new TGeoVolume(Tyvec_name, TyvecS);
    TGeoMedium *Tyvec_medium = pMedTyvec; //set medium
    if(Tyvec_medium) TyvecV->SetMedium(Tyvec_medium);
    else Fatal("Main", "Invalid medium for Tyvec!");

    //Slot's shapes
    TGeoShape *slot_FeS = new TGeoBBox(slot_Fe_name, Zdc_VRFM_X_size, Zdc_VRFM_Y_size, Zdc_VRFM_Z_size);
    TGeoVolume *slot_FeV = new TGeoVolume(slot_Fe_name, slot_FeS);
    TGeoMedium *slot_Fe_medium = pMedSlot; //set medium
    if(slot_Fe_medium) slot_FeV->SetMedium(slot_Fe_medium);
    else Fatal("Main", "Invalid medium for slot_Fe!");

    TGeoShape *slot_PbS = new TGeoBBox(slot_Pb_name, Zdc_VRPM_X_size, Zdc_VRPM_Y_size, Zdc_VRPM_Z_size);
    TGeoVolume *slot_PbV = new TGeoVolume(slot_Pb_name, slot_PbS);
    TGeoMedium *slot_Pb_medium = pMedSlot; //set medium
    if(slot_Pb_medium) slot_PbV->SetMedium(slot_Pb_medium);
    else Fatal("Main", "Invalid medium for slot_Pb!");

    TGeoShape *slot_ScS = new TGeoBBox(slot_Sc_name, Zdc_VRSM_X_size, Zdc_VRSM_Y_size, Zdc_VRSM_Z_size);
    TGeoVolume *slot_ScV = new TGeoVolume(slot_Sc_name, slot_ScS);
    TGeoMedium *slot_Sc_medium = pMedSlot; //set medium
    if(slot_Sc_medium) slot_ScV->SetMedium(slot_Sc_medium);
    else Fatal("Main", "Invalid medium for slot_Sc!");

    TGeoShape *slot_Tyvec_YS = new TGeoBBox(slot_Tyvec_Y_name, Zdc_VRYM_X_size, Zdc_VRYM_Y_size, Zdc_VRYM_Z_size);
    TGeoVolume *slot_Tyvec_YV = new TGeoVolume(slot_Tyvec_Y_name, slot_Tyvec_YS);
    TGeoMedium *slot_Tyvec_Y_medium = pMedSlot; //set medium
    if(slot_Tyvec_Y_medium) slot_Tyvec_YV->SetMedium(slot_Tyvec_Y_medium);
    else Fatal("Main", "Invalid medium for slot_Tyvec_Y!");

    TGeoShape *slot_Tyvec_AS = new TGeoBBox(slot_Tyvec_A_name, Zdc_VRAM_X_size, Zdc_VRAM_Y_size, Zdc_VRAM_Z_size);
    TGeoVolume *slot_Tyvec_AV = new TGeoVolume(slot_Tyvec_A_name, slot_Tyvec_AS);
    TGeoMedium *slot_Tyvec_A_medium = pMedSlot; //set medium
    if(slot_Tyvec_A_medium) slot_Tyvec_AV->SetMedium(slot_Tyvec_A_medium);
    else Fatal("Main", "Invalid medium for slot_Tyvec_A!");
    
    //insert Fe
    mother_volume->AddNode(FeV,1,new TGeoTranslation(0, 0, -(Zdc_module_Z_size-Zdc_Fe_Z_size)));
    FeV->AddNode(slot_FeV,1,new TGeoTranslation(Zdc_Fe_X_size-Zdc_VRFM_X_size, Zdc_Y_slot_position, 0.));
    mother_volume->AddNode(FeV,2,new TGeoTranslation(0, 0, (Zdc_module_Z_size-Zdc_Fe_Z_size)));
    
    //insert Tyvec
    mother_volume->AddNode(TyvecV, 1, new TGeoTranslation(0,0, -(Zdc_module_Z_size-2.*Zdc_Fe_Z_size-Zdc_Tyvec_Z_size)));
    cout <<"tyvec 1 " <<-(Zdc_module_Z_size-2.*Zdc_Fe_Z_size-Zdc_Tyvec_Z_size) <<endl;
    TyvecV->AddNode(slot_Tyvec_YV,1,new TGeoTranslation((Zdc_Tyvec_X_size-Zdc_VRYM_X_size), Zdc_Y_slot_position, 0));
    TyvecV->AddNode(slot_Tyvec_AV, 1, new TGeoTranslation((Zdc_Tyvec_X_size-Zdc_VRAM_X_size-2.*Zdc_VRYM_X_size), Zdc_Y_slot_position, -(Zdc_Tyvec_Z_size-Zdc_VRAM_Z_size)));
    TyvecV->AddNode(slot_Tyvec_AV, 2, new TGeoTranslation((Zdc_Tyvec_X_size-Zdc_VRAM_X_size-2.*Zdc_VRYM_X_size), Zdc_Y_slot_position,(Zdc_Tyvec_Z_size-Zdc_VRAM_Z_size)));
    
    //insert Scint into Tyvec
    TyvecV->AddNode(zdcSensitiveVolumeV, 1, new TGeoTranslation(0., 0., 0.));
    zdcSensitiveVolumeV->AddNode(slot_ScV, 1, new TGeoTranslation((Zdc_Tyvec_X_size-Zdc_VRAM_X_size-2.*Zdc_VRYM_X_size), Zdc_Y_slot_position, 0));
    
    //insert Pb
    mother_volume->AddNode(PbV, 2, new TGeoTranslation(0, 0, -(Zdc_module_Z_size-2.*Zdc_Fe_Z_size-2.*Zdc_Tyvec_Z_size-Zdc_Pb_Z_size)));
    cout <<"Pb 2 " <<-(Zdc_module_Z_size-2.*Zdc_Fe_Z_size-2.*Zdc_Tyvec_Z_size-Zdc_Pb_Z_size) <<endl;
    PbV->AddNode(slot_PbV,1,new TGeoTranslation((Zdc_Pb_X_size-Zdc_VRPM_X_size), Zdc_Y_slot_position, 0));
    slot_PbV->SetLineColor(kBlue);
    
    //insert Tyvec slices
    Double_t itInit = -(Zdc_module_Z_size-2.*Zdc_Fe_Z_size-Zdc_Tyvec_Z_size); 
    Double_t itCur=0, itDelta=2.*(Zdc_Pb_Z_size+Zdc_Tyvec_Z_size);
    for(Int_t it=1+1; it<Zdc_NTyvec_in_module+1; it++) {//60 slices (1-60), the first slice is already inserted
      itCur =itInit + (it-1)*itDelta;
      mother_volume->AddNode(TyvecV,it,new TGeoTranslation(0,0,itCur));
      cout <<"it,itCur Ty " <<it <<" " <<itCur <<endl;
    }

    //insert Pb slices    
    itInit=-(Zdc_module_Z_size-2.*Zdc_Fe_Z_size-2.*Zdc_Tyvec_Z_size-Zdc_Pb_Z_size); 
    itCur=0; itDelta=2.*(Zdc_Pb_Z_size+Zdc_Tyvec_Z_size);
    for(Int_t it=1+2; it<Zdc_NPb_in_module+2; it++) {//59 slices (2-60), the first slice is already inserted
      itCur =itInit + (it-2)*itDelta;
      mother_volume->AddNode(PbV,it,new TGeoTranslation(0,0,itCur));
      cout <<"it,itCur Pb " <<it <<" " <<itCur <<endl;
    }
    

    //Fill module with hole
    TString FeH_name = "FeH";
    TString PbH_name = "PbH";
    TString SensitiveVolumeH_name = "ScH";
    TString TyvecH_name = "TyvecH";

    TString FeHTube_name = "FeHTube";
    TString PbHTube_name = "PbHTube";
    TString ScHTube_name = "ScHTube";
    TString TyvecHTube_name = "TyvecHTube";

    //Volumes where hole will be inserted
    TGeoShape *FeHS = new TGeoBBox(FeH_name, Zdc_Fe_X_size, Zdc_Fe_Y_size, Zdc_Fe_Z_size);
    TGeoVolume *FeHV = new TGeoVolume(FeH_name, FeHS);
    TGeoMedium *FeH_medium = pMedIron; //set medium
    if(FeH_medium) FeHV->SetMedium(FeH_medium);
    else Fatal("Main", "Invalid medium for FeH!");

    //FeHV->SetLineColor(TColor::GetColor("#fab700"));
    ////FeHV->SetVisibility(false);
    ////FeHV->SetTransparency(0);

    TGeoShape *PbHS = new TGeoBBox(PbH_name, Zdc_Pb_X_size, Zdc_Pb_Y_size, Zdc_Pb_Z_size);
    TGeoVolume *PbHV = new TGeoVolume(PbH_name, PbHS);
    TGeoMedium *PbH_medium = pMedLead; //set medium
    if(PbH_medium) PbHV->SetMedium(PbH_medium);
    else Fatal("Main", "Invalid medium for PbH!");

    TGeoShape *zdcSensitiveVolumeHS = new TGeoBBox(SensitiveVolumeH_name, Zdc_Sc_X_size, Zdc_Sc_Y_size, Zdc_Sc_Z_size);
    TGeoVolume *zdcSensitiveVolumeHV = new TGeoVolume(SensitiveVolumeH_name, zdcSensitiveVolumeHS);
    TGeoMedium *ScH_medium = pMedScint; //set medium
    if(ScH_medium) zdcSensitiveVolumeHV->SetMedium(Sc_medium);
    else Fatal("Main", "Invalid medium for ScH!");

    TGeoShape *TyvecHS = new TGeoBBox(TyvecH_name, Zdc_Tyvec_X_size, Zdc_Tyvec_Y_size, Zdc_Tyvec_Z_size);
    TGeoVolume *TyvecHV = new TGeoVolume(TyvecH_name, TyvecHS);
    TGeoMedium *TyvecH_medium = pMedTyvec; //set medium
    if(TyvecH_medium) TyvecHV->SetMedium(TyvecH_medium);
    else Fatal("Main", "Invalid medium for TyvecH!");

    //Hole elements
    TGeoTube *FeHTubeS = new TGeoTube(FeHTube_name, Zdc_Hole_inner_radius, Zdc_Hole_outer_radius, Zdc_Fe_Z_size);
    TGeoVolume *FeHTubeV = new TGeoVolume(FeHTube_name, FeHTubeS);
    FeHTubeV->SetMedium(pMedVacuum);

    TGeoTube *PbHTubeS = new TGeoTube(PbHTube_name, Zdc_Hole_inner_radius, Zdc_Hole_outer_radius, Zdc_Pb_Z_size);
    TGeoVolume *PbHTubeV = new TGeoVolume(PbHTube_name, PbHTubeS);
    PbHTubeV->SetMedium(pMedVacuum);

    TGeoTube *ScHTubeS = new TGeoTube(ScHTube_name, Zdc_Hole_inner_radius, Zdc_Hole_outer_radius, Zdc_Sc_Z_size);
    TGeoVolume *ScHTubeV = new TGeoVolume(ScHTube_name, ScHTubeS);
    ScHTubeV->SetMedium(pMedVacuum);

    TGeoTube *TyvecHTubeS = new TGeoTube(TyvecHTube_name, Zdc_Hole_inner_radius, Zdc_Hole_outer_radius, Zdc_Hole_Tyvec_Z_size);
    TGeoVolume *TyvecHTubeV = new TGeoVolume(TyvecHTube_name, TyvecHTubeS);
    TyvecHTubeV->SetMedium(pMedVacuum);
    
    //insert Fe into module with hole
    mother_volume_hole->AddNode(FeHV,1,new TGeoTranslation(0, 0, -(Zdc_module_Z_size-Zdc_Fe_Z_size)));
    FeHV->AddNode(slot_FeV,1,new TGeoTranslation(Zdc_Fe_X_size-Zdc_VRFM_X_size, Zdc_Y_slot_position, 0));
    FeHV->AddNode(FeHTubeV,1,new TGeoTranslation(0., 0., 0.));
    mother_volume_hole->AddNode(FeHV,2,new TGeoTranslation(0, 0, (Zdc_module_Z_size-Zdc_Fe_Z_size)));

            
    //insert Pb into module with hole
    mother_volume_hole->AddNode(PbHV, 2, new TGeoTranslation(0, 0, -(Zdc_module_Z_size-2.*Zdc_Fe_Z_size-2.*Zdc_Tyvec_Z_size-Zdc_Pb_Z_size)));
    cout <<"Pb 2 hole " <<-(Zdc_module_Z_size-2.*Zdc_Fe_Z_size-2.*Zdc_Tyvec_Z_size-Zdc_Pb_Z_size) <<endl;
    PbHV->AddNode(slot_PbV,1,new TGeoTranslation((Zdc_Pb_X_size-Zdc_VRPM_X_size), Zdc_Y_slot_position, 0));
    PbHV->AddNode(PbHTubeV,1,new TGeoTranslation(0., 0., 0.));


    //insert Tyvec into module with hole
    mother_volume_hole->AddNode(TyvecHV, 1, new TGeoTranslation(0,0, -(Zdc_module_Z_size-2.*Zdc_Fe_Z_size-Zdc_Tyvec_Z_size)));
    cout <<"tyvec 1 hole " <<-(Zdc_module_Z_size-2.*Zdc_Fe_Z_size-Zdc_Tyvec_Z_size) <<endl;
    TyvecHV->AddNode(slot_Tyvec_YV,1,new TGeoTranslation((Zdc_Tyvec_X_size-Zdc_VRYM_X_size), Zdc_Y_slot_position, 0));
    TyvecHV->AddNode(slot_Tyvec_AV, 1, new TGeoTranslation((Zdc_Tyvec_X_size-Zdc_VRAM_X_size-2.*Zdc_VRYM_X_size), Zdc_Y_slot_position, -(Zdc_Tyvec_Z_size-Zdc_VRAM_Z_size)));
    TyvecHV->AddNode(slot_Tyvec_AV, 2, new TGeoTranslation((Zdc_Tyvec_X_size-Zdc_VRAM_X_size-2.*Zdc_VRYM_X_size), Zdc_Y_slot_position,(Zdc_Tyvec_Z_size-Zdc_VRAM_Z_size)));
    TyvecHV->AddNode(TyvecHTubeV,1,new TGeoTranslation(0., 0., -(Zdc_Tyvec_Z_size-Zdc_VRAM_Z_size)));
    TyvecHV->AddNode(TyvecHTubeV,2,new TGeoTranslation(0., 0., (Zdc_Tyvec_Z_size-Zdc_VRAM_Z_size)));

    //insert Sc into Tyvec
    TyvecHV->AddNode(zdcSensitiveVolumeHV, 1, new TGeoTranslation(0., 0., 0.));
    zdcSensitiveVolumeHV->AddNode(slot_ScV, 1, new TGeoTranslation((Zdc_Tyvec_X_size-Zdc_VRAM_X_size-2.*Zdc_VRYM_X_size), Zdc_Y_slot_position, 0));
    zdcSensitiveVolumeHV->AddNode(ScHTubeV,1,new TGeoTranslation(0., 0., 0.));

    //insert Tyvec slices
    itInit = -(Zdc_module_Z_size-2.*Zdc_Fe_Z_size-Zdc_Tyvec_Z_size); 
    itCur=0, itDelta=2.*(Zdc_Pb_Z_size+Zdc_Tyvec_Z_size);
    for(Int_t it=1+1; it<Zdc_NTyvec_in_module+1; it++) {//60 slices (1-60), the first slice is already inserted
      itCur =itInit + (it-1)*itDelta;
      mother_volume_hole->AddNode(TyvecHV,it,new TGeoTranslation(0,0,itCur));
      cout <<"it,itCur Ty " <<it <<" " <<itCur <<endl;
    }

    //insert Pb slices    
    itInit=-(Zdc_module_Z_size-2.*Zdc_Fe_Z_size-2.*Zdc_Tyvec_Z_size-Zdc_Pb_Z_size); 
    itCur=0; itDelta=2.*(Zdc_Pb_Z_size+Zdc_Tyvec_Z_size);
    for(Int_t it=1+2; it<Zdc_NPb_in_module+2; it++) {//59 slices (2-60), the first slice is already inserted
      itCur =itInit + (it-2)*itDelta;
      mother_volume_hole->AddNode(PbHV,it,new TGeoTranslation(0,0,itCur));
      cout <<"it,itCur Pb " <<it <<" " <<itCur <<endl;
    }

}
//FillModule

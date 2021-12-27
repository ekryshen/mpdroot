// ********************************************************************
// *                                                                  *
// *  Projet:   mpd.jinr.ru                                           *
// *  Created:   30-October-2020                                     *
// *  Author: M. Martemianov / MPD collaboration                      *
// *                                                                  *
// ********************************************************************

#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoMedium.h"
#include "TFile.h"

#include "../mpdloadlibs.C"


// Media

TGeoMedium *pMedAir = 0;
TGeoMedium *pMedAluminium = 0;
TGeoMedium *pMedIron = 0;
TGeoMedium *pMedFscScint = 0;
TGeoMedium *pMedLead = 0;
TGeoMedium *pMedKapton = 0;
TGeoMedium *pMedGlue = 0;
TGeoMedium *pMedCarbComp = 0;
TGeoMedium *pMedFbGlass = 0;
TGeoMedium *pMedTiPaint = 0;
TGeoMedium *pMedDural = 0;

// Settings to check geometry quality

Bool_t fDebug = kTRUE;
Double_t fPrecision = 0.00001;

// Power frame of EMC

Double_t EMC_Length_Inner = 0.5*624.4; // half-length of lower power frame
Double_t EMC_Length_Outer = 0.5*826.0; // half-length of upper power frame
Double_t EMC_Frame_InnerRadius = 168.0, EMC_Frame_OuterRadius = 229.3; // Inner and outer radii
const Double_t EMC_Frame_Wall1 = 2.5, EMC_Frame_Wall2 = 1.0, EMC_Frame_Wall3 = 2.0; // Wall widths
const Double_t EMC_Frame_gap1 = 0.2; // additional gaps
const Double_t EMC_Frame_Ledge = 40; // z - size of frame ledge
const Double_t EMC_Frame_Angle = 2*4.734; // angle of the frame ledge
const Double_t EMC_Wall_Add = 6.0;

// Sector parameters

const Int_t EMC_NSectors = 25; // number of sectors in calorimeter
const Int_t EMC_Index_Box = 64; // number of towers along z-axis
const Int_t EMC_Box_Layers = 210; // number of layers in one calorimeter cell
const Double_t EMC_Box_Length = 41.55; // tower length (cm)
const Int_t EMC_NSector_div1 = 6, EMC_NSector_div2 = 8; // number of towers in module
const Double_t EMC_Sector_Length = 312.0; // sector length
const Double_t EMC_Sector_Wall = 0.6; // basket walls 
const Double_t EMC_Sector_InnerFrame = 170.7; // sector inner radius
const Double_t EMC_Sector_OuterRadius = 227.1; // sector outer radius
const Double_t EMC_Sector_InnerRadius = EMC_Sector_InnerFrame+EMC_Sector_Wall; // module inner radius
//// const Double_t EMC_Tower_InnerRadius = 171.56; // tower inner radius
const Double_t EMC_Tower_InnerRadius = 171.3; // tower inner radius
const Double_t EMC_Sector_UpperWall = 0.5; // upper wall of basket 
const Double_t EMC_Sector_RotWallPos = 262.5; 
const Double_t EMC_Sector_RotWall = 0.56; // inclined wall width of basket 


// Tower parameters

Double_t EMC_Box_Sc = 0.15; // scintillator width (cm)
Double_t EMC_Box_Pb = 0.03; // lead width (cm)
Double_t EMC_Box_Col = 0.005; // paint width on lead (cm)
Double_t EMC_Box_Cell = EMC_Box_Sc + EMC_Box_Pb + 2.*EMC_Box_Col;// total size of layer
Double_t EMC_Box_Plastic1 = 0.7, EMC_Box_Plastic2 = 0.8; // fixing plates width
Double_t EMC_Box_Angle1 = 0.5*0.9*TMath::DegToRad(); // cell rotation angle1 (milling angle)
Double_t EMC_Box_Angle2 = 0.5*1.2*TMath::DegToRad(); // cell rotation angle2 (milling angle)
Double_t EMC_Box_Bottom_X = 3.32, EMC_Box_Top_X = 3.97; // cell bottom/top size (cm)
Double_t EMC_Box_Bottom_Y = 3.31, EMC_Box_Top_Y = 4.00; // cell bottom/top size (cm)
Double_t* EMC_Box_A = new Double_t[EMC_Index_Box]; // cell top sizes
Double_t* EMC_Box_B = new Double_t[EMC_Index_Box]; // cell top sizes
Double_t* EMC_Box_C = new Double_t[EMC_Index_Box]; // cell bottom sizes
Double_t EMC_Box_gap = 0.02; // gap between cells -> 0.02
Double_t EMC_Module_side = EMC_Box_Length/cos(EMC_Box_Angle1);
Double_t EMC_Length; // half-length of ECal sector

Double_t EMC_Sector_angle0 = 0.;
Double_t EMC_Sector_angle = TMath::TwoPi()/EMC_NSectors;
Double_t EMC_Module_angle = EMC_Sector_angle/EMC_NSector_div1;

// Dural plates parameters

Double_t EMC_Box_angle, EMC_Plate_angle;
const Int_t EMC_NumberPlates = EMC_NSector_div2;
Double_t EMC_TowerYshift[EMC_NumberPlates];
Double_t EMC_TowerZshift[EMC_NumberPlates];

// Functions

class FairGeoMedia;
class FairGeoBuilder;
void DefineMediaParameters(FairGeoMedia* , FairGeoBuilder*);
void CreateECALFrame(); // Power frame settings
void CreateECALSectors(); // Sector geometry settings
void ComputeECALTowerParams(); // additional tower sizes
void CreateECALModShape(); // Shape for module
void CreateECALModule(); // Module (2x8 towers)
void CreateECALModuleSupport(); // Supporting elements of modules / towers
void CreateECALTower(Int_t); // tower as one trapezoid
void CreateECALCompTower(Int_t, Double_t); // tower as two trapezoids
void CreateECALCompTower(Int_t, Double_t, Double_t); // tower as three trepezoids
void CreateECALStructure(Int_t, TGeoVolume*, TList*);
void CreateECALTowerEdges(TList*, Double_t, Double_t, Double_t*, Double_t*);

const Int_t numBoxSize = 100;
Double_t* CalculateZPos(); // central position of each tower in module
Double_t* EMC_zPosition = new Double_t[numBoxSize];
	
TGeoVolume *emcChH;
TGeoVolume *emcSector; 
TGeoVolume *emcCrate; 

// List of different objects

TList* listECALModules = new TList();
TList* listECALTowers = new TList();
TList* listECALSupport = new TList();
TList* listECALTowerCuts = new TList();
TList* listECALPlateVert = new TList();


//
// Root function to create EMC/MPD geometry for root - file
//

void create_rootgeom_emc_v4() {

// Load MPD libraries


// Set working directory

  TString gPath = gSystem->Getenv("VMCWORKDIR");
  const TString geoDetectorName = "emcChamber";
  const TString geoDetectorVersion = "v4";
  TString geoFileName = gPath + "/geometry/emc_"+ geoDetectorVersion + ".root";

// Global geometry parameters

  FairGeoLoader*    geoLoad = new FairGeoLoader("TGeo","FairGeoLoader");
  FairGeoInterface* geoFace = geoLoad->getGeoInterface();
  
// Load media from media file

  TString medFile = gPath + "/geometry/media.geo";
  geoFace->setMediaFile(medFile);
  geoFace->readMedia();
  FairGeoMedia*   geoMedia = geoFace->getMedia();
  FairGeoBuilder* geoBuild = geoLoad->getGeoBuilder();
  DefineMediaParameters(geoMedia, geoBuild);


// Create geometry and global top volume

  gGeoManager = (TGeoManager*)gROOT->FindObject("FAIRGeom");
  gGeoManager->SetName(geoDetectorName + "_geom");
  TGeoVolume* top = new TGeoVolumeAssembly("TOP");
  top->SetMedium(pMedAir);
  gGeoManager->SetTopVolume(top);

// Additional parameters and updates

  EMC_Box_Top_X =   EMC_Box_Bottom_X + 2.*EMC_Box_Length*tan(EMC_Box_Angle1);
  EMC_zPosition = CalculateZPos(); // tower Z - center position
  EMC_Length = EMC_Sector_Length; 
  printf("MpdRoot : ECal half length  %f \n", EMC_Length);  
  
  Double_t rotAngleMax = (2*EMC_Index_Box-1)*EMC_Box_Angle1; 
  Double_t xPos = EMC_zPosition[EMC_Index_Box-1] - 0.5*(EMC_Box_Length + EMC_Box_Bottom_X/sin(rotAngleMax))*
		  sin(rotAngleMax) + 0.5*EMC_Box_Bottom_X/cos(rotAngleMax);

///  printf("MpdRoot : calculated last point of tower line along z : %f \n", xPos); 
///  printf("MpdRoot : tabulated %f \n", EMC_Sector_RotWallPos + EMC_Sector_Wall); 


// ECal barrel

  EMC_Length_Outer = EMC_Length_Outer - EMC_Frame_Ledge; 
  TGeoTube *EMC_ChamberShape = new TGeoTube(EMC_Frame_InnerRadius, EMC_Frame_OuterRadius, EMC_Length_Outer);
  TGeoVolume *EMC_ChamberVolume = new TGeoVolume(geoDetectorName,EMC_ChamberShape,pMedAir);
  EMC_ChamberVolume->SetVisibility(kFALSE);

// ECal right / left side

  Double_t rotation_angle = 90.;
  while  ( rotation_angle > 0 ) rotation_angle -= EMC_Sector_angle*TMath::RadToDeg();
  rotation_angle += EMC_Sector_angle*TMath::RadToDeg();
  EMC_Sector_angle0 = rotation_angle;

  TGeoTube *EMC_ChHShape = new TGeoTube(EMC_Frame_InnerRadius, EMC_Frame_OuterRadius, EMC_Length_Outer/2.);
  emcChH = new TGeoVolume("emcChH",EMC_ChHShape,pMedAir);
  emcChH->SetVisibility(kFALSE); emcChH->SetLineColor(42);
  
  TGeoCombiTrans* rotBarrelPlus = new TGeoCombiTrans();
  rotBarrelPlus->SetTranslation(0., 0., EMC_Length_Outer/2.);
  TGeoCombiTrans* rotBarrelMinus = new TGeoCombiTrans();
  rotBarrelMinus->ReflectZ(true);
  rotBarrelMinus->SetTranslation(0., 0., -EMC_Length_Outer/2.);
  EMC_ChamberVolume->AddNode(emcChH,0,rotBarrelPlus);
  EMC_ChamberVolume->AddNode(emcChH,1,rotBarrelMinus);

// Create ECAL power frame 

  CreateECALFrame();

// ECAL sectors (emcSector)

  CreateECALSectors();

// Additional ECAL tower parameters

  ComputeECALTowerParams();

// Create ECAL module support structure

  CreateECALModuleSupport();
  
//  printf("MpdRoot : ECal sector weight : %.2f kg\n", emcSector->Weight(0.01,"a")); 

// Create ECAL module shape (with glue content)

  CreateECALModShape();

// Create ECAL module

  CreateECALModule();

  top->AddNode(EMC_ChamberVolume,0,new TGeoTranslation( 0.0, 0.0, 0.0));
  top->SetVisContainers(kTRUE);

  printf("MpdRoot : ECal total weight : %.2f kg\n", EMC_ChamberVolume->Weight(0.01,"a")); 

  gGeoManager->CloseGeometry();
  gGeoManager->CheckOverlaps(fPrecision);
  gGeoManager->PrintOverlaps();
  gGeoManager->Test();
  gGeoManager->SetTopVisible(kTRUE);

// Total view draw

  TCanvas *c1 = new TCanvas("c1","ECal power frame",10,10,585,585);
  c1->cd();
   
  gGeoManager->SetVisLevel(3);
  gGeoManager->SetNsegments(20);
  emcCrate->SetVisibility(kFALSE);
  top->Draw();
  TView3D *view1 = (TView3D*)gPad->GetView();
  view1->ShowAxis();
  c1->Update();

// Draw modules one line (crate)

  TGeoCompositeShape* shapeSector = (TGeoCompositeShape*)emcSector->GetShape();
  Double_t rot_angle0 = 2.*EMC_Box_Angle2*TMath::RadToDeg();
  Double_t crate_radius = ((TGeoPcon*)emcCrate->GetShape())->GetRmin()[0];
  Double_t dZPos = ((TGeoPcon*)emcCrate->GetShape())->GetDZ();
  Double_t rotAngle = EMC_Plate_angle;
  Double_t fModRadius = crate_radius + (EMC_Tower_InnerRadius - EMC_Sector_InnerRadius);

  TCanvas *c2 = new TCanvas("c2","ECal modules one line",100,10,1000,300);
  gGeoManager->SetVisLevel(3);

  TGeoTube *chamberShapeSel = new TGeoTube(0.0, 350., 200.);
  TGeoVolume *chamberVolumeSel = new TGeoVolume("chamberVolumeSel",chamberShapeSel,pMedAir);
  chamberVolumeSel->SetVisibility(kFALSE); chamberVolumeSel->SetVisContainers(kTRUE);
  emcCrate->SetVisibility(kFALSE); 

  TGeoCombiTrans *tr_matSel = new TGeoCombiTrans();
  tr_matSel->RotateZ(-rot_angle0); tr_matSel->RotateY(-90.); tr_matSel->RotateZ(180.);
  tr_matSel->SetTranslation(dZPos, 0.0, -(crate_radius));
  chamberVolumeSel->AddNode(emcCrate, 0, tr_matSel);

  chamberVolumeSel->Draw();
  gPad->SetLeftMargin(0.01); gPad->SetBottomMargin(0.01);
  TView3D *view2 = (TView3D*)gPad->GetView();
  view2->SetParallel();
  view2->Front();
  view2->ShowAxis();

  TAxis3D *axis2 = TAxis3D::GetPadAxis();
   if (axis2) {
      axis2->SetLabelSize(0.06);
      axis2->GetXaxis()->SetLabelOffset(-0.1);
      axis2->GetXaxis()->SetTitleSize(0.07);
      axis2->GetXaxis()->SetTitleOffset(1.1);
      axis2->GetXaxis()->SetTitle("Z, cm");
    }
  c2->Update();

// Draw one module shape 

/*
 TCanvas *c3 = new TCanvas("c3","ECal module",200,10,700,580);
 
  Int_t iModule = 7; // iModule = 0, 7

  TGeoTube *chamberShapeAdd = new TGeoTube(0.0, 200., 200.);
  TGeoVolume *chamberVolumeAdd = new TGeoVolume("chamberVolumeAdd",chamberShapeAdd,pMedAir);
  chamberVolumeAdd->SetVisibility(kFALSE); chamberVolumeAdd->SetVisContainers(kTRUE);

  TGeoVolume* volModSel = (TGeoVolume*)listECALModules->At(iModule);
  TGeoVolume* volPlateSel = (TGeoVolume*)listECALSupport->At(iModule);
  
  for (Int_t iSel = 0; iSel < volModSel->GetNdaughters(); iSel++) {
    TGeoNode* selNode = volModSel->GetNode(iSel);
    TGeoVolume* selVol = (TGeoVolume*)selNode->GetVolume();
    for (Int_t iPart = 0; iPart < selVol->GetNdaughters(); iPart++) {
       TGeoVolume* selPart = (TGeoVolume*)(selVol->GetNode(iPart))->GetVolume();
       selPart->SetVisibility(kFALSE); 	
    }
  }

   TGeoCombiTrans *tr_matSelect = new TGeoCombiTrans();   
   TGeoCombiTrans *tr_plateSelect = new TGeoCombiTrans();
   tr_matSelect->RotateZ(-rot_angle0); tr_matSelect->RotateY(-90.); tr_matSelect->RotateZ(180.);
   tr_matSelect->SetTranslation(0.0, 0.0, -(fModRadius)); 
   tr_plateSelect->RotateX(90.); tr_plateSelect->RotateY(rotAngle);
   tr_plateSelect->SetTranslation(EMC_TowerZshift[iModule], 0.0, EMC_TowerYshift[iModule]);

   chamberVolumeAdd->AddNode(volModSel, 0, tr_matSelect);
   chamberVolumeAdd->AddNode(volPlateSel, 0, tr_plateSelect);
   chamberVolumeAdd->Draw();

  gPad->SetLeftMargin(0.01); gPad->SetBottomMargin(0.01);
  TView3D *view3 = (TView3D*)gPad->GetView();
  view3->SetParallel();
  view3->Front();
  view3->ShowAxis();

  TAxis3D *axis3 = TAxis3D::GetPadAxis();
   if (axis3) {
      axis3->SetLabelSize(0.06);
      axis3->GetXaxis()->SetLabelOffset(-0.1);
      axis3->GetXaxis()->SetTitleSize(0.07);
      axis3->GetXaxis()->SetTitleOffset(1.1);
      axis3->GetXaxis()->SetTitle("Z, cm");
    }

   char addTxt[256];
   sprintf(addTxt,"Module N %02d \n",iModule+1);
   TPaveText* modTextSel = new TPaveText(0.40,0.80,0.60,0.90,"blNDC");
    modTextSel->SetBorderSize(0);
    modTextSel->SetTextFont(42); modTextSel->SetLineColor(28);
    modTextSel->SetTextSize(0.05);
    modTextSel->SetFillColor(0);
    TText *modCaptionSel = modTextSel->AddText(addTxt);
    modTextSel->Draw();

  c3->Update();

// Draw one tower

  TCanvas *c4 = new TCanvas("c4","ECal tower",100,150,600,600);

  Int_t iCal = 1; // iCal = 0, 63
  c4->cd();

   TGeoVolume* towerBoxSel = (TGeoVolume*)listECALTowers->At(iCal);
   towerBoxSel->SetVisibility(kTRUE);
   for (int iNode = 0; iNode < towerBoxSel->GetNdaughters(); iNode++) {
    TGeoNode* selNode = towerBoxSel->GetNode(iNode);
    TGeoVolume* selVol = (TGeoVolume*)selNode->GetVolume();
    selVol->SetVisibility(kFALSE);
   }

   towerBoxSel->Draw(); towerBoxSel->SetVisibility(kTRUE);
   towerBoxSel->SetVisContainers(kFALSE);
   TView3D *view4 = (TView3D*)gPad->GetView();
   view4->SetParallel();

   view4->Front();
   view4->ShowAxis();

   TAxis3D *axis4 = TAxis3D::GetPadAxis();
   if (axis4) {
      axis4->SetLabelSize(0.04);
      axis4->GetXaxis()->SetLabelSize(0.04);
      axis4->GetXaxis()->SetLabelOffset(-0.1);
      axis4->GetXaxis()->SetNdivisions(4);
      axis4->GetXaxis()->SetTitleSize(0.07);
    }

   char addText[256];
   sprintf(addText,"Tower N %02d \n",iCal+1);
   TPaveText* towerTextSel = new TPaveText(0.40,0.80,0.60,0.90,"blNDC");
    towerTextSel->SetBorderSize(0);
    towerTextSel->SetTextFont(42); towerTextSel->SetLineColor(28);
    towerTextSel->SetTextSize(0.05);
    towerTextSel->SetFillColor(0);
    TText *towerCaptionSel = towerTextSel->AddText(addText);
    towerTextSel->Draw();

  c4->Update();

// Draw sector (basket) structure

  dZPos = shapeSector->GetDZ();
  rot_angle0 = 0.5*TMath::TwoPi()/EMC_NSectors*TMath::RadToDeg(); 
  crate_radius = shapeSector->GetOrigin()[0] - shapeSector->GetDX();
  fModRadius = crate_radius + (EMC_Tower_InnerRadius - EMC_Sector_InnerRadius);
  
  TCanvas *c5 = new TCanvas("c5","ECal sector",450,100,1000,500);
  gGeoManager->SetVisLevel(2);
  gGeoManager->SetNsegments(20);

  TGeoTube *chamberSectorShapeSel = new TGeoTube(0.0, 350., 200.);
  TGeoVolume *chamberSectorSel = new TGeoVolume("chamberSectorSel",chamberSectorShapeSel,pMedAir);
  chamberSectorSel->SetVisibility(kFALSE); chamberSectorSel->SetVisContainers(kTRUE);
  emcCrate->SetVisibility(kFALSE); 

  TGeoCombiTrans *tr_matrixSel = new TGeoCombiTrans();
  tr_matrixSel->RotateZ(-rot_angle0); tr_matrixSel->RotateY(-90.); tr_matrixSel->RotateZ(180.);
  tr_matrixSel->SetTranslation(dZPos, 0.0, -(crate_radius));
  chamberSectorSel->AddNode(emcSector, 0, tr_matrixSel);

  chamberSectorSel->Draw();
  gPad->SetLeftMargin(0.01); gPad->SetBottomMargin(0.01);
  TView3D *view5 = (TView3D*)gPad->GetView();
  view5->SetParallel();
  view5->Draw();
  view5->ShowAxis();

  TAxis3D *axis5 = TAxis3D::GetPadAxis();
   if (axis5) {
      axis5->SetLabelSize(0.06);
      axis5->GetXaxis()->SetLabelOffset(-0.1);
      axis5->GetXaxis()->SetTitleSize(0.07);
      axis5->GetXaxis()->SetTitleOffset(1.1);
      axis5->GetXaxis()->SetTitle("Z, cm");
    }

  c5->Update();
*/

// Write to root - file 

  TFile* geomFile = new TFile(geoFileName, "RECREATE");
  top->Write();
  geomFile->Close();

}

// Define media parameters

void DefineMediaParameters(FairGeoMedia* geoMedia, FairGeoBuilder* geoBuild) {

// Air

    FairGeoMedium* mAir = geoMedia->getMedium("air");
    if ( ! mAir ) Fatal("Main", "FairGeoMedium air not found");
    geoBuild->createMedium(mAir);
    pMedAir = gGeoManager->GetMedium("air");
    pMedAir->GetMaterial()->Print();
    if ( ! pMedAir ) Fatal("Main", "TGeoMedium air not found");

// Aluminium

    FairGeoMedium* mAluminium = geoMedia->getMedium("aluminium");
    if ( ! mAluminium  ) Fatal("Main", "FairGeoMedium aluminium not found");
    geoBuild->createMedium(mAluminium);
    pMedAluminium  = gGeoManager->GetMedium("aluminium");
    pMedAluminium->GetMaterial()->Print();
    if ( ! pMedAluminium  ) Fatal("Main", "TGeoMedium aluminium not found");

// Scintillator
   
    FairGeoMedium* mFscScint = geoMedia->getMedium("FscScint");
    if ( ! mFscScint ) Fatal("Main", "FairGeoMedium FscScint not found");
    geoBuild->createMedium(mFscScint);
    pMedFscScint = gGeoManager->GetMedium("FscScint");
    pMedFscScint->GetMaterial()->Print();
    if ( ! pMedFscScint ) Fatal("Main", "TGeoMedium FscScint not found");

// Lead    

    FairGeoMedium* mLead = geoMedia->getMedium("lead");
    if ( ! mLead ) Fatal("Main", "FairGeoMedium Lead not found");
    geoBuild->createMedium(mLead);
    pMedLead = gGeoManager->GetMedium("lead");
    pMedLead->GetMaterial()->Print();
    if ( ! pMedLead ) Fatal("Main", "TGeoMedium Lead not found");

// Plastic 

    FairGeoMedium* mKapton = geoMedia->getMedium("kapton");
    if ( ! mKapton ) Fatal("Main", "FairGeoMedium Kapton not found");
    geoBuild->createMedium(mKapton);
    pMedKapton = gGeoManager->GetMedium("kapton");
    pMedKapton->GetMaterial()->Print();
    if ( ! pMedKapton ) Fatal("Main", "TGeoMedium Kapton not found");  

// Fiberglass material

    FairGeoMedium* mFbGlass = geoMedia->getMedium("fiberglass");
    if ( ! mFbGlass ) Fatal("Main", "FairGeoMedium mFbGlass not found");
    geoBuild->createMedium(mFbGlass);
    pMedFbGlass = gGeoManager->GetMedium("fiberglass");
    pMedFbGlass->GetMaterial()->Print();
    if ( ! pMedFbGlass ) Fatal("Main", "TGeoMedium fiberglass not found");  

// ECal glue with Ti support

    FairGeoMedium* mGlue = geoMedia->getMedium("glueTi");
    if ( ! mGlue ) Fatal("Main", "FairGeoMedium mGlue not found");
    geoBuild->createMedium(mGlue);
    pMedGlue = gGeoManager->GetMedium("glueTi");
    pMedGlue->GetMaterial()->Print();
    if ( ! pMedGlue ) Fatal("Main", "TGeoMedium glueTi not found");  

// Ti02 reflecting color (EJ510)

    FairGeoMedium* mTiPaint = geoMedia->getMedium("TiO2Col");
    if ( ! mTiPaint ) Fatal("Main", "FairGeoMedium mTiPaint not found");
    geoBuild->createMedium(mTiPaint);
    pMedTiPaint = gGeoManager->GetMedium("TiO2Col");
    pMedTiPaint->GetMaterial()->Print();
    if ( ! pMedTiPaint ) Fatal("Main", "TGeoMedium TiO2Col not found");  

// Carbon composite (Graphite epoxy suprt)

    FairGeoMedium* mCarbComp = geoMedia->getMedium("carboncomp");
    if ( ! mCarbComp ) Fatal("Main", "FairGeoMedium mCarbComp not found");
    geoBuild->createMedium(mCarbComp);
    pMedCarbComp = gGeoManager->GetMedium("carboncomp");
    pMedCarbComp->GetMaterial()->Print();
    if ( ! pMedCarbComp ) Fatal("Main", "TGeoMedium carboncomp not found");  

// Dural

    FairGeoMedium* mDural = geoMedia->getMedium("dural");
    if ( ! mDural ) Fatal("Main", "FairGeoMedium mDural not found");
    geoBuild->createMedium(mDural);
    pMedDural = gGeoManager->GetMedium("dural");
    pMedDural->GetMaterial()->Print();
    if ( ! pMedDural ) Fatal("Main", "TGeoMedium dural not found");


}

// Central tower Z - position

Double_t* CalculateZPos() {

  Double_t gap = EMC_Box_gap, rotation_angle = 0.;
  Double_t* zPos = new Double_t[numBoxSize];
  Double_t* zPositions = new Double_t[numBoxSize];
  zPositions[0] = 2.*EMC_Sector_Wall + 2.*EMC_Box_gap;
  for (int ip = 0; ip < EMC_Index_Box; ip++) {
   gap = EMC_Box_gap; rotation_angle = (2*ip + 1)*EMC_Box_Angle1;
   if ( (ip == 1*EMC_Index_Box/EMC_NSector_div2) || (ip == 2*EMC_Index_Box/EMC_NSector_div2) || 
	(ip == 3*EMC_Index_Box/EMC_NSector_div2) || (ip == 4*EMC_Index_Box/EMC_NSector_div2) ||
        (ip == 5*EMC_Index_Box/EMC_NSector_div2) || (ip == 6*EMC_Index_Box/EMC_NSector_div2) || 
	(ip == 7*EMC_Index_Box/EMC_NSector_div2) ) gap = 3.*EMC_Box_gap;
    zPositions[ip+1] = zPositions[ip] + gap/cos(rotation_angle-EMC_Box_Angle1) +
            EMC_Box_Bottom_X*(cos(rotation_angle)+sin(rotation_angle)*tan(rotation_angle-EMC_Box_Angle1));
    zPos[ip] = zPositions[ip+1] - 0.5*(EMC_Box_Bottom_X*cos(rotation_angle) -
                        EMC_Box_Length*sin(rotation_angle));
   }
  return zPos;
}

// Create ECAL power frame 

void CreateECALFrame() {
 
 Double_t sizeZPos = 0.0; 
 Double_t fFrameInner = EMC_Frame_InnerRadius + EMC_Frame_Wall1;
 Double_t fFrameOuter = EMC_Frame_OuterRadius - EMC_Frame_Wall3; 
 Double_t fEdgeAngle = TMath::ASin(0.5*EMC_Frame_Wall2/fFrameOuter); 
 Double_t fEdgeLength = 0.5*EMC_Frame_Wall2/tan(fEdgeAngle) - fFrameInner;
 Double_t fFrameEdgeDiff = EMC_Length_Outer-EMC_Length_Inner;
 Double_t fFrameAngleCut = TMath::ATan(fEdgeLength/(fFrameEdgeDiff+EMC_Frame_Ledge));
 Double_t fFrameLenCut = fFrameEdgeDiff*tan(fFrameAngleCut);
 Double_t xPosTru[5] = {fFrameInner, fFrameInner+fEdgeLength, fFrameInner+fEdgeLength, 
			fFrameInner+fFrameLenCut, fFrameInner};
 Double_t yPosTru[5] = {0.0, 0.0, EMC_Length_Outer, EMC_Length_Outer,  EMC_Length_Inner};

// Parallel edges 

  TGeoXtru *frameEdgeShape1 = new TGeoXtru(2);
  frameEdgeShape1->DefinePolygon(5, xPosTru, yPosTru);
  frameEdgeShape1->DefineSection(0,-0.5*EMC_Frame_Wall2,0.,0.,1.0);
  frameEdgeShape1->DefineSection(1, 0.5*EMC_Frame_Wall2,0.,0.,1.0);

  TGeoVolume *frameEdge1 = new TGeoVolume("frameEdge1",frameEdgeShape1,pMedCarbComp);
  frameEdge1->SetLineColor(12);

  TGeoTubeSeg *frameLedgeShape = new TGeoTubeSeg(fFrameOuter, EMC_Frame_OuterRadius, 0.5*EMC_Frame_Ledge, 0., EMC_Frame_Angle);
  TGeoVolume *frameLedge = new TGeoVolume("frameLedge",frameLedgeShape, pMedCarbComp);
  frameLedge->SetLineColor(12);

// Upper and lower power parts

  TGeoPcon* frameLatticeShape1 = new TGeoPcon(0.,360.,2);
  frameLatticeShape1->DefineSection(0, 0.5*EMC_Length_Inner, EMC_Frame_InnerRadius, fFrameInner);
  frameLatticeShape1->DefineSection(1,-0.5*EMC_Length_Inner, EMC_Frame_InnerRadius, fFrameInner);
  TGeoVolume *frameLattice1 = new TGeoVolume("frameLattice1",frameLatticeShape1,pMedCarbComp);
  emcChH->AddNode(frameLattice1,0,new TGeoTranslation(0.,0.,-0.5*fFrameEdgeDiff));

  TGeoPcon* frameLatticeShape2 = new TGeoPcon(0.,360.,2);
  frameLatticeShape2->DefineSection(0, 0.5*EMC_Length_Outer, fFrameOuter, EMC_Frame_OuterRadius);
  frameLatticeShape2->DefineSection(1,-0.5*EMC_Length_Outer, fFrameOuter, EMC_Frame_OuterRadius);
  TGeoVolume *frameLattice2 = new TGeoVolume("frameLattice2",frameLatticeShape2,pMedCarbComp);
  emcChH->AddNode(frameLattice2,0,new TGeoTranslation(0.,0.,0.));
  frameLattice1->SetLineColor(12); frameLattice2->SetLineColor(12);

// Add ledge and edges

  Double_t rotation_angle = EMC_Sector_angle0;
  sizeZPos = 0.5*EMC_Length_Outer - 0.5*EMC_Frame_Ledge;
  for (int iEdge = 0; iEdge < EMC_NSectors; iEdge++) {
   TGeoCombiTrans *edge_trans1 = new TGeoCombiTrans();
   TGeoCombiTrans *ledge_trans1 = new TGeoCombiTrans();
   edge_trans1->RotateX(90.); edge_trans1->RotateZ(rotation_angle);
   edge_trans1->SetTranslation(0.,0.,-0.5*EMC_Length_Outer); 
   emcChH->AddNode(frameEdge1,iEdge,edge_trans1); 
   rotation_angle += EMC_Sector_angle*TMath::RadToDeg();
  }

}

// Create ECAL sectors and walls 

void CreateECALSectors() {

  Double_t xPos, yPos, zPos, rangle, fEMC_Pos, fEMC_Z0, fEMC_Z1, fEMC_Z2, dZPos;
  Double_t fHoleSector, fSectorAngle, fSectorRad, fSectorRadius, fGapWallAngle;
  Double_t fInnerRadius, fOuterRadius, fFrameWidth, rotation_angle, xLoc, yLoc;

  fEMC_Z0 = -EMC_Sector_Length/2., fEMC_Z1 = EMC_Sector_Length/2.;
  fFrameWidth = EMC_Sector_InnerRadius - EMC_Sector_InnerFrame;
  fHoleSector  = 0.5*EMC_Frame_Wall2 + EMC_Frame_gap1;
  fGapWallAngle = TMath::ASin(fHoleSector/EMC_Sector_InnerFrame);
  fSectorAngle  = 2.*(0.5*EMC_Sector_angle - fGapWallAngle);
  xPos = fHoleSector/tan(0.5*EMC_Sector_angle); yPos = fHoleSector;
  fSectorRad = EMC_Sector_InnerFrame - xPos;

  TGeoRotation *coneRot = new TGeoRotation("coneRot");
  coneRot->RotateZ(fGapWallAngle*TMath::RadToDeg());
  TGeoTranslation *coneTrans = new TGeoTranslation("coneTrans", xPos, yPos, 0.);
  coneRot->RegisterYourself(); coneTrans->RegisterYourself();

// Create sector (composite shape)

  TGeoPcon *emcSectorAdditional = new TGeoPcon(0.,fSectorAngle*TMath::RadToDeg(),2);
  emcSectorAdditional->SetName("emcSectorAdditional");
  emcSectorAdditional->DefineSection(0, fEMC_Z0, EMC_Sector_InnerFrame, EMC_Sector_OuterRadius);
  emcSectorAdditional->DefineSection(1, fEMC_Z1, EMC_Sector_InnerFrame, EMC_Sector_OuterRadius);
  fOuterRadius = fSectorRad + EMC_Sector_OuterRadius - EMC_Sector_InnerFrame;

  TGeoPcon *emcSectorPcone = new TGeoPcon(0.,EMC_Sector_angle*TMath::RadToDeg(),2);
  emcSectorPcone->SetName("emcSectorPcone");
  emcSectorPcone->DefineSection(0, fEMC_Z0, fSectorRad, fOuterRadius);
  emcSectorPcone->DefineSection(1, fEMC_Z1, fSectorRad, fOuterRadius);

  TGeoCompositeShape *emcSectorShape = new TGeoCompositeShape("emcSectorShape",
        "emcSectorAdditional:coneRot+emcSectorPcone:coneTrans");
  emcSector = new TGeoVolume("emcSector", emcSectorShape, pMedAir);
  emcSector->SetVisibility(kFALSE); emcSector->SetVisContainers(kTRUE);

  rotation_angle = EMC_Sector_angle0;
  dZPos = fEMC_Z0 + ((TGeoTube*)emcChH->GetShape())->GetDZ() ;
  for (Int_t iSec = 0; iSec < EMC_NSectors; iSec++) {
   TGeoCombiTrans *combi_sector = new TGeoCombiTrans();
   combi_sector->RotateZ(rotation_angle); 
   combi_sector->SetTranslation(0.,0.,-dZPos); 
   emcChH->AddNode(emcSector,iSec,combi_sector);
   rotation_angle += EMC_Sector_angle*TMath::RadToDeg();
  }

// Lower wall

  TGeoTubeSeg* emcSectorWallShape1 = new TGeoTubeSeg(EMC_Sector_InnerFrame, EMC_Sector_InnerRadius,
				     EMC_Sector_Length/2., 0.,fSectorAngle*TMath::RadToDeg());
  TGeoVolume* emcSectorWall1 = new TGeoVolume("emcSectorWall1", emcSectorWallShape1, pMedFbGlass);
  emcSectorWall1->SetLineColor(8);
  emcSector->AddNode(emcSectorWall1,0,coneRot);

// Upper wall

  TGeoTubeSeg* emcSectorWallShape2 = new TGeoTubeSeg(EMC_Sector_OuterRadius - EMC_Sector_UpperWall, 
			EMC_Sector_OuterRadius, EMC_Sector_Length/2., 0.,fSectorAngle*TMath::RadToDeg());
  TGeoVolume* emcSectorWall2 = new TGeoVolume("emcSectorWall2", emcSectorWallShape2, pMedFbGlass);
  emcSectorWall2->SetLineColor(8);
  emcSector->AddNode(emcSectorWall2,0,coneRot);

// Parrallel walls in sector

  Double_t fEdgeWidth = 2.*EMC_Sector_Wall;
  Double_t fWallInnerRadius, fWallOuterRadius, fEdgeRadius, fEdgeLength, fLastRadius;
  xPos = xPos + 0.5*EMC_Sector_Wall/tan(0.5*EMC_Sector_angle);
  yPos = yPos + 0.5*EMC_Sector_Wall;

  fWallInnerRadius = EMC_Sector_InnerRadius - xPos;
  fWallOuterRadius = fWallInnerRadius + (EMC_Sector_OuterRadius - EMC_Sector_InnerRadius - EMC_Sector_UpperWall);
  fEdgeRadius = fWallOuterRadius - (fWallInnerRadius + sqrt(xPos*xPos+yPos*yPos) - EMC_Sector_InnerRadius);
  fEdgeLength = fEdgeRadius - fWallInnerRadius;
  fLastRadius = fWallInnerRadius;  

  TGeoXtru *emcSectorWallShape3 = new TGeoXtru(2);
  Double_t xTru[5] = {fWallInnerRadius, fEdgeRadius, fEdgeRadius, fWallInnerRadius};
  Double_t yTru[5] = {fEMC_Z0+fEdgeWidth, fEMC_Z0+fEdgeWidth, fEMC_Z1, fEMC_Z1};
  emcSectorWallShape3->DefinePolygon(4, xTru, yTru);
  emcSectorWallShape3->DefineSection(0,-0.5*EMC_Sector_Wall,0.,0.,1.0);
  emcSectorWallShape3->DefineSection(1, 0.5*EMC_Sector_Wall,0.,0.,1.0);
  TGeoVolume *emcSectorWall3 = new TGeoVolume("emcSectorWall3", emcSectorWallShape3, pMedFbGlass);
  emcSectorWall3->SetLineColor(33);

  rotation_angle = 0.0;
  for (Int_t iWall = 0; iWall < 2; iWall++) {
   TGeoCombiTrans *combi_wall1 = new TGeoCombiTrans();
   combi_wall1->RotateX(90.); combi_wall1->RotateZ(rotation_angle*TMath::RadToDeg());
   combi_wall1->SetTranslation(xPos, yPos, 0.);
   emcSector->AddNode(emcSectorWall3,iWall,combi_wall1);
   rotation_angle += EMC_Sector_angle;
  }

// Horizontal and inclined walls in sector

  fWallOuterRadius = EMC_Sector_OuterRadius-EMC_Sector_UpperWall;
  
  TGeoPcon *emcSectorWallShape4 = new TGeoPcon(0.,fSectorAngle*TMath::RadToDeg(),2);
  emcSectorWallShape4->DefineSection(0, fEMC_Z0, EMC_Sector_InnerRadius, fWallOuterRadius);
  emcSectorWallShape4->DefineSection(1, fEMC_Z0+fEdgeWidth, EMC_Sector_InnerRadius, fWallOuterRadius);
  TGeoVolume *emcSectorWall4 = new TGeoVolume("emcSectorWall4", emcSectorWallShape4, pMedFbGlass);
  emcSectorWall4->SetLineColor(8);
  emcSector->AddNode(emcSectorWall4,0,coneRot);

  Double_t startAngle = 2.0*asin(0.5*EMC_Sector_Wall/EMC_Sector_InnerRadius);
  Double_t zWallPos[4], rWallMin[4], rWallMax[4], rWallPos; 
  rotation_angle = 2*EMC_Index_Box*EMC_Box_Angle1;  
  rWallPos = (EMC_Sector_Length - EMC_Sector_RotWallPos - fEdgeWidth - 
	     EMC_Sector_RotWall/cos(rotation_angle))/tan(rotation_angle) + EMC_Sector_InnerRadius;
  zWallPos[0] = fEdgeWidth + EMC_Sector_RotWallPos - EMC_Sector_Length/2.;  
  zWallPos[1] = zWallPos[0] + EMC_Sector_RotWall/cos(rotation_angle);
  zWallPos[2] = EMC_Sector_Length/2. - EMC_Sector_RotWall/cos(rotation_angle);
  zWallPos[3] = EMC_Sector_Length/2.; 
  rWallMin[0] = EMC_Sector_InnerRadius; rWallMax[0] = EMC_Sector_InnerRadius;
  rWallMin[1] = EMC_Sector_InnerRadius; 
  rWallMax[1] = EMC_Sector_InnerRadius + EMC_Sector_RotWall/sin(rotation_angle);
  rWallMin[2] = rWallPos -  EMC_Sector_RotWall/sin(rotation_angle); rWallMax[2] = rWallPos; 
  rWallMin[3] = rWallPos; rWallMax[3] = rWallPos; 
 
  TGeoPcon *emcSectorWallShape5 = new 
	TGeoPcon(startAngle*TMath::RadToDeg(),(fSectorAngle-2.*startAngle)*TMath::RadToDeg(),4);
  for (Int_t iWall = 0; iWall < 4; iWall++) 
    emcSectorWallShape5->DefineSection(iWall, zWallPos[iWall], rWallMin[iWall], rWallMax[iWall]);
  TGeoVolume *emcSectorWall5 = new TGeoVolume("emcSectorWall5", emcSectorWallShape5, pMedFbGlass);
  emcSectorWall5->SetLineColor(8);
  emcSector->AddNode(emcSectorWall5,0,coneRot);

// Four small parts to support inclined wall

  Double_t supAngle = 2.0*asin(0.5*EMC_Sector_RotWall/EMC_Sector_InnerRadius);

  TGeoPcon *emcSectorWallShape6 = new TGeoPcon(0.,supAngle*TMath::RadToDeg(),2);
  emcSectorWallShape6->DefineSection(0, zWallPos[1], EMC_Sector_InnerRadius, EMC_Sector_InnerRadius); 
  emcSectorWallShape6->DefineSection(1, EMC_Sector_Length/2., EMC_Sector_InnerRadius, rWallPos);
  TGeoVolume *emcSectorWall6 = new TGeoVolume("emcSectorWall6", emcSectorWallShape6, pMedFbGlass);
  emcSectorWall6->SetLineColor(33);

  TGeoPcon *emcSectorWallShape7 = new TGeoPcon(0.,2.*supAngle*TMath::RadToDeg(),2);
  emcSectorWallShape7->DefineSection(0, zWallPos[1], EMC_Sector_InnerRadius, EMC_Sector_InnerRadius); 
  emcSectorWallShape7->DefineSection(1, EMC_Sector_Length/2., EMC_Sector_InnerRadius, rWallPos);
  TGeoVolume *emcSectorWall7 = new TGeoVolume("emcSectorWall7", emcSectorWallShape7, pMedFbGlass);
  emcSectorWall7->SetLineColor(33);

  Double_t dStepAngle;
  Int_t iCount1 = 0, iCount2 = 0; 
  rotation_angle = startAngle + fGapWallAngle; 

  for (Int_t iWall = 0; iWall < 4; iWall++) {
   TGeoCombiTrans *combi_wall2 = new TGeoCombiTrans();
   combi_wall2->RotateZ(rotation_angle*TMath::RadToDeg());   
   dStepAngle = supAngle; 
   if ( ( iWall == 1 ) || (iWall == 2) ) { 
    dStepAngle = 2.*supAngle; 
    emcSector->AddNode(emcSectorWall7,iCount1,combi_wall2); iCount1++; 
   }
   else {
    emcSector->AddNode(emcSectorWall6,iCount2,combi_wall2); iCount2++;
   }
   rotation_angle += dStepAngle + (fSectorAngle - 2.*startAngle - 6.*supAngle)/3.;
  }

// Subsector (crate)

  Double_t dWidthGlue = EMC_Box_gap;
  Double_t fCrateRad, fCrateRad0, fCrateOuterRadius, fCrateAngle;
  xPos = fHoleSector/tan(0.5*EMC_Sector_angle); yPos = fHoleSector;
  xPos = xPos + EMC_Sector_Wall/tan(0.5*EMC_Sector_angle); 
  yPos = yPos + EMC_Sector_Wall; 

  fCrateRad0 = dWidthGlue/tan(0.5*EMC_Module_angle);
  fWallInnerRadius = EMC_Sector_InnerRadius - xPos;
  fCrateRad = fWallInnerRadius - dWidthGlue/tan(0.5*EMC_Module_angle); 
  fCrateOuterRadius = fCrateRad + fEdgeLength*fCrateRad/fLastRadius;
  rotation_angle = 2*EMC_Index_Box*EMC_Box_Angle1;
  fInnerRadius = fCrateRad + (EMC_Length - EMC_Sector_RotWallPos - fEdgeWidth)/tan(rotation_angle);

  fEMC_Z0 = -EMC_Length/2. + 2.*EMC_Sector_Wall; 
  fEMC_Z1 = -EMC_Length/2. + 2.*EMC_Sector_Wall + EMC_Sector_RotWallPos; 
  fEMC_Z2 =  EMC_Length/2.;

  TGeoPcon *emcCratePcone = new TGeoPcon(0.,EMC_Module_angle*TMath::RadToDeg(),3);
  emcCratePcone->DefineSection(0, fEMC_Z0, fCrateRad, fCrateOuterRadius);
  emcCratePcone->DefineSection(1, fEMC_Z1, fCrateRad, fCrateOuterRadius);
  emcCratePcone->DefineSection(2, fEMC_Z2, fInnerRadius, fCrateOuterRadius);
  emcCrate = new TGeoVolume("emcCrate", emcCratePcone, pMedAir);
  emcCrate->SetVisibility(kTRUE); emcCrate->SetVisContainers(kTRUE);
  emcCrate->SetLineColor(1);

  rotation_angle = 0.0;
  for (Int_t iCrate = 0; iCrate < EMC_NSector_div1; iCrate++) {
    TGeoCombiTrans *combi_crate = new TGeoCombiTrans();
    combi_crate->RotateZ(rotation_angle*TMath::RadToDeg());
    xLoc = xPos + fCrateRad0*cos(0.5*EMC_Module_angle + iCrate*EMC_Module_angle);
    yLoc = yPos + fCrateRad0*sin(0.5*EMC_Module_angle + iCrate*EMC_Module_angle); 
    combi_crate->SetTranslation(xLoc, yLoc, 0.0);
    emcSector->AddNode(emcCrate,iCrate,combi_crate);
    rotation_angle += EMC_Module_angle;
  }

 }

// Compute tower parameters

void ComputeECALTowerParams() {

  Double_t rotation_angle, arbLen;
  Double_t SideSize = (EMC_Box_Top_Y - EMC_Box_Bottom_Y)/(2.*tan(EMC_Box_Angle2));

  for (Int_t index = 0; index < EMC_Index_Box; index++) {
   rotation_angle = (2*index+1)*EMC_Box_Angle1;
   arbLen = (SideSize - EMC_Box_Bottom_X*sin(rotation_angle))*cos(EMC_Box_Angle1)/
             cos(rotation_angle-EMC_Box_Angle1);
   EMC_Box_A[index] = EMC_Box_Top_Y;
   EMC_Box_C[index] = EMC_Box_Bottom_Y + 2.*EMC_Box_Bottom_X*sin(rotation_angle)*tan(EMC_Box_Angle2);
   if (arbLen > EMC_Box_Length) {arbLen = EMC_Box_Length;
     EMC_Box_A[index] = EMC_Box_C[index] +
             2.*EMC_Box_Length*cos(rotation_angle-EMC_Box_Angle1)*tan(EMC_Box_Angle2);}
   EMC_Box_B[index] = EMC_Box_Bottom_Y + 2.*arbLen*cos(rotation_angle+EMC_Box_Angle1)*tan(EMC_Box_Angle2);
  }

}

//  Create ECAL shape of the different modules

void CreateECALModShape() {

// Set specific parameters of module shape

  Int_t shapeDim[EMC_NSector_div2], index; 
  Double_t Z12P[EMC_NSector_div2][EMC_Index_Box], R1P[EMC_NSector_div2][EMC_Index_Box], 
           R2P[EMC_NSector_div2][EMC_Index_Box];
  Double_t ZtowerInCrate[EMC_NSector_div2*EMC_NSector_div2], ZminiP[EMC_NSector_div2], xPos, yPos, fModRadius;
  Double_t fSizeMin, fSizeMax, fEMC_Z0 = -EMC_Length/2., crate_radius, xBorderMax, xLast, yLast, xPos0, yPos0; 
  Double_t xVert[4], yVert[4], rotTower1, rotTower2, dAngleRotate, dPlateLen, dCutHeight, dLenMax;
  TString nameECALModule, nameECALShape, rotCutMatrix, transName;
  Double_t hZF = EMC_Box_Bottom_X, hZB = hZF + 2.0*EMC_Box_Length*tan(EMC_Box_Angle1);
  Double_t gap, AM[EMC_Index_Box+1], Zright[EMC_Index_Box+1], Zleft[EMC_Index_Box+1], 
           Yleft[EMC_Index_Box+1], Zmin[EMC_Index_Box+1], Ymin[EMC_Index_Box+1],
           Zmax[EMC_Index_Box+1], Ymax[EMC_Index_Box+1], Z2gap[EMC_Index_Box+1];

  Z2gap[0] = 0.; AM[0]= EMC_Box_Angle1; Zright[0]=0.;
  fEMC_Z0 = fEMC_Z0 + EMC_Box_gap; // add space for glue

  for(Int_t i = 0; i < EMC_Index_Box; i++) {
  AM[i]=  2.0*EMC_Box_Angle1*i + EMC_Box_Angle1; gap = EMC_Box_gap;
  if ( (i == 1*EMC_Index_Box/EMC_NSector_div2-1) || (i == 2*EMC_Index_Box/EMC_NSector_div2-1) ||
       (i == 3*EMC_Index_Box/EMC_NSector_div2-1) || (i == 4*EMC_Index_Box/EMC_NSector_div2-1) ||
       (i == 5*EMC_Index_Box/EMC_NSector_div2-1) || (i == 6*EMC_Index_Box/EMC_NSector_div2-1) ||
       (i == 7*EMC_Index_Box/EMC_NSector_div2-1) ) gap = 3.*EMC_Box_gap;

  if(i == 0) Zright[i] = 2.*EMC_Sector_Wall + 2.*EMC_Box_gap + hZF*cos(AM[i]);
  else Zright[i] = Z2gap[i-1] + hZF*(cos(AM[i]) + sin(AM[i])*tan(AM[i] - EMC_Box_Angle1));
   ZtowerInCrate[i] = Zright[i] - 0.5*hZF*cos(AM[i]) + 0.5*EMC_Box_Length*sin(AM[i]);
   Zleft[i] =Zright[i] - hZF*cos(AM[i]); Yleft[i] = EMC_Sector_InnerRadius+sin(AM[i])*hZF;
   Zmax[i] = Zright[i] + EMC_Module_side*sin(AM[i] + EMC_Box_Angle1);
   Ymax[i]= EMC_Sector_InnerRadius + EMC_Module_side*cos(AM[i] + EMC_Box_Angle1);
   Zmin[i] = Zleft[i] + EMC_Module_side*sin(AM[i] - EMC_Box_Angle1);
   Ymin[i] = Yleft[i] + EMC_Module_side*cos(AM[i] - EMC_Box_Angle1);
   Z2gap[i] = Zright[i]+gap/cos(AM[i] + EMC_Box_Angle1);
  }

  Int_t j, k, iMod, Nmin, Nmax;
  Double_t Z0, ZPgonRmin[EMC_Index_Box], ZPgonRmax[EMC_Index_Box], YPgonRmin[EMC_Index_Box],
           YPgonRmax[EMC_Index_Box], ZPgon[EMC_Index_Box], R1Pgon[EMC_Index_Box], R2Pgon[EMC_Index_Box];

  for(Int_t jk = 0; jk < EMC_NSector_div2; jk++) {
   iMod =jk; Nmin = iMod*EMC_NSector_div2; Nmax = Nmin+EMC_NSector_div2;
   Z0 = Zleft[Nmin] - 2.*EMC_Box_gap/cos(AM[Nmin]-EMC_Box_Angle1);
   ZminiP[jk] = Z0; j=0;

  for(Int_t i = Nmin; i < Nmax;i++) {
   if(i == Nmin) {
    ZPgonRmin[j] = 0; YPgonRmin[j] = Yleft[i] - EMC_Sector_InnerRadius;
    ZPgonRmax[j] = Zmin[i]-Zleft[i]; YPgonRmax[j] = Ymin[i] - EMC_Sector_InnerRadius; j++;
   }

  if(i < Nmax - 1) {
   ZPgonRmin[j] = Zleft[i]-Z0; YPgonRmin[j] = Yleft[i] - EMC_Sector_InnerRadius;
   ZPgonRmax[j] = Zmin[i]-Z0; YPgonRmax[j] = Ymin[i] - EMC_Sector_InnerRadius; j++;
   ZPgonRmin[j] = Zright[i]-Z0; YPgonRmin[j] = 0;
   ZPgonRmax[j] = Zmax[i]-Z0; YPgonRmax[j] = Ymax[i] - EMC_Sector_InnerRadius; j++;
   ZPgonRmin[j] = Zright[i]-Z0; YPgonRmin[j] = 0;
   ZPgonRmax[j] = Zmin[i+1]-Z0; YPgonRmax[j] = Ymin[i+1] - EMC_Sector_InnerRadius; j++;
   ZPgonRmin[j] = Zright[i]-Z0; YPgonRmin[j] = 0;
   ZPgonRmax[j] = Zmin[i+1]-Z0; YPgonRmax[j] = Ymin[i+1] - EMC_Sector_InnerRadius; j++;
  }
  else {
   ZPgonRmin[j] = Zleft[i]-Z0; YPgonRmin[j] = Yleft[i]-EMC_Sector_InnerRadius;
   ZPgonRmax[j] = Zmin[i]-Z0; YPgonRmax[j] = Ymin[i] - EMC_Sector_InnerRadius; j++;
   ZPgonRmin[j] = Zright[i]-Z0; YPgonRmin[j] = 0;
   ZPgonRmax[j] = Zmax[i]-Z0; YPgonRmax[j] = Ymax[i] - EMC_Sector_InnerRadius; j++;
   ZPgonRmin[j] = Zright[i]-Z0; YPgonRmin[j] = 0;
   ZPgonRmax[j] = Zmax[i]-Z0; YPgonRmax[j] = Ymax[i] - EMC_Sector_InnerRadius; j++;
   }
  }

  Double_t ZPoly[EMC_Index_Box], R1[EMC_Index_Box], R2[EMC_Index_Box], z1, z2, y1, y2;
  for( Int_t k1 = 0; k1 < EMC_NSector_div2; k1++)
   ZminiP[k1] = Zleft[k1*EMC_NSector_div2] - 2*EMC_Box_gap/cos(AM[k1*EMC_NSector_div2] - EMC_Box_Angle1);

  for(Int_t i = 0; i < EMC_Index_Box/2; i++) {
   ZPoly[i] = ZPgonRmin[i]; R1[i] = YPgonRmin[i]; k = 0;
  while(ZPgonRmax[k] < ZPoly[i]) k++;
   if(i == 0) R2[i] = R1[i];
   else {
    if(k == 0) {z1 = ZPgonRmin[0]; z2 = ZPgonRmax[0]; y1 = YPgonRmin[0]; y2 = YPgonRmax[0];}
    else{z1 = ZPgonRmax[k-1]; z2 = ZPgonRmax[k]; y1 = YPgonRmax[k-1]; y2 = YPgonRmax[k];}
    R2[i] = y1 + (y2-y1)/(z2-z1)*(ZPoly[i]-z1);
   }
  }

  for(Int_t i=0; i < EMC_Index_Box/2; i++) {
   ZPoly[EMC_Index_Box-i-1] = ZPgonRmax[EMC_Index_Box/2-1-i];
   R2[EMC_Index_Box-i-1] = YPgonRmax[EMC_Index_Box/2-1-i]; k=0;
   while (ZPgonRmin[EMC_Index_Box/2-1-k] > ZPoly[EMC_Index_Box-1-i]) k++;

   if(i == 0) R1[EMC_Index_Box-1] = R2[EMC_Index_Box-1];
   else {
    if(k == 0) {z1 = ZPgonRmin[EMC_Index_Box/2-1]; z2 = ZPgonRmax[EMC_Index_Box/2-1];
                y1 = YPgonRmin[EMC_Index_Box/2-1]; y2 = YPgonRmax[EMC_Index_Box/2-1];}
   else {z1 = ZPgonRmin[EMC_Index_Box/2-k-1]; z2 = ZPgonRmin[EMC_Index_Box/2-k];
         y1= YPgonRmin[EMC_Index_Box/2-k-1]; y2 = YPgonRmin[EMC_Index_Box/2-k];}
         R1[EMC_Index_Box-i-1] = y1 +(y2-y1)/(z2-z1)*(ZPoly[EMC_Index_Box-i-1]-z1);}
   }

  for(Int_t i = 0; i < EMC_Index_Box - 1; i++){
   z1 = 100;
  for(Int_t i2=i+1; i2 < EMC_Index_Box-2; i2++){
   z2 = ZPoly[i2]-ZPoly[i];
   if(z2 <z1) {z1 = z2; k = i2; }}
  if(k-i >1) {y1 = ZPoly[i+1]; ZPoly[i+1] = ZPoly[k]; ZPoly[k] = y1;
   y1 = R1[i+1]; R1[i+1] = R1[k]; R1[k] = y1;
   y1 = R2[i+1]; R2[i+1] = R2[k]; R2[k] = y1;}
  }

   Int_t Ndel = 0;
   Z12P[iMod][0] = ZPoly[0]; R1P[iMod][0] = R1[0]; R2P[iMod][0] = R2[0];
   for(Int_t i = 1; i < EMC_Index_Box; i++){
   if((ZPoly[i] - ZPoly[i-1]) < 0.0001) Ndel++;
   Z12P[iMod][i-Ndel] = ZPoly[i]; R1P[iMod][i-Ndel] = R1[i]; R2P[iMod][i-Ndel] = R2[i];}
   shapeDim[jk] = EMC_Index_Box - Ndel;
  }

  dAngleRotate = EMC_Plate_angle*TMath::DegToRad();
  crate_radius = ((TGeoPcon*)emcCrate->GetShape())->GetRmin()[0];
  fModRadius = crate_radius + (EMC_Tower_InnerRadius - EMC_Sector_InnerRadius);

// Create ECAL module shape

  TGeoRotation *transRot = new TGeoRotation();
  TGeoRotation *plateRot = new TGeoRotation();

  Int_t nWallSect; 
  Double_t zGlue[2*EMC_NSector_div2], rGlueMax[2*EMC_NSector_div2]; 
  Double_t dRotAngle, dzWidth, dzPos, dLenH, dSideGlue, gPos[4], gRmin[4], gRmax[4];

  Int_t nSect = ((TGeoPcon*)emcCrate->GetShape())->GetNz();
  Double_t *crateRadMin = ((TGeoPcon*)emcCrate->GetShape())->GetRmin(); 
  Double_t *crateRadMax = ((TGeoPcon*)emcCrate->GetShape())->GetRmax(); 
  Double_t *crateZ = ((TGeoPcon*)emcCrate->GetShape())->GetZ(); 
  Double_t phiAngle = ((TGeoPcon*)emcCrate->GetShape())->GetDphi();

  Double_t dWidthGlue = EMC_Box_gap; 
  Double_t fHoleSector  = 0.5*EMC_Frame_Wall2 + EMC_Frame_gap1;
  xPos0 = fHoleSector/tan(0.5*EMC_Sector_angle); yPos0 = fHoleSector;
  xPos0 = xPos0 + (EMC_Sector_Wall + 0.5*dWidthGlue)/tan(0.5*EMC_Sector_angle);
  yPos0 = yPos0 + (EMC_Sector_Wall + 0.5*dWidthGlue);
  Double_t fGlueRad = EMC_Sector_InnerRadius - xPos0;

  for(Int_t iShape = 0; iShape < EMC_NSector_div2; iShape++) {

   nameECALShape.Form("shapeModule%02d",iShape);
   TGeoPgon *emcModuleShape = new TGeoPgon(0.,EMC_Module_angle*TMath::RadToDeg(),2,shapeDim[iShape]);
   emcModuleShape->SetName(nameECALShape.Data());

   for (Int_t iSize = 0; iSize < shapeDim[iShape]; iSize++) {
    fSizeMin = fModRadius + R1P[iShape][iSize];
    fSizeMax = fModRadius + R2P[iShape][iSize];
    emcModuleShape->DefineSection(iSize,Z12P[iShape][iSize],fSizeMin,fSizeMax);
   }

// Create vertical cut in module shape to move dural plate

   TGeoXtru *shapePlateVert = (TGeoXtru*)listECALPlateVert->At(iShape);
   for (Int_t iPoint = 0; iPoint < shapePlateVert->GetNvert(); iPoint++) {
       xVert[iPoint] = shapePlateVert->GetX(iPoint); yVert[iPoint] = shapePlateVert->GetY(iPoint);
   }

   Double_t* pRmax = emcModuleShape->GetRmax();
   xLast = emcModuleShape->GetZ(emcModuleShape->GetNz()-1);
   yLast = emcModuleShape->GetRmax(emcModuleShape->GetNz()-1);
   if (iShape == 0) dLenMax = 2.*emcModuleShape->GetDX();

   rotTower1 = 2.*EMC_NSector_div2*iShape*EMC_Box_Angle1;
   rotTower2 = 2.*EMC_NSector_div2*(iShape+1)*EMC_Box_Angle1;
   xBorderMax = emcModuleShape->GetZ((Int_t)TMath::LocMax(emcModuleShape->GetNz(),pRmax));
   if (iShape == 0) xBorderMax = - fPrecision;

   dPlateLen = xLast - EMC_TowerZshift[iShape] - (yLast - EMC_TowerYshift[iShape] - fModRadius)*tan(rotTower2);
   dPlateLen = dPlateLen*(1. - tan(rotTower2)*tan(dAngleRotate)/
               (1. + tan(rotTower2)*tan(dAngleRotate)))/cos(dAngleRotate);
   xVert[1] = dPlateLen*cos(dAngleRotate); yVert[1] = - dPlateLen*sin(dAngleRotate);
   dCutHeight = (EMC_Box_Bottom_X*sin(rotTower2+EMC_Box_Angle1) + dLenMax*cos(rotTower2)/
                                  cos(EMC_Box_Angle1) - EMC_TowerYshift[iShape] - yVert[1])/cos(rotTower2);
   xVert[2] = xVert[1] + dCutHeight*sin(rotTower2);
   yVert[2] = yVert[1] + dCutHeight*cos(rotTower2);
   yVert[3] = 2.*emcModuleShape->GetDX() - EMC_TowerYshift[iShape];
   xVert[3] = xBorderMax - EMC_TowerZshift[iShape];
   xVert[0] = xVert[3] - yVert[3]*tan(rotTower1); yVert[0] = 0.0;
   shapePlateVert->DefinePolygon(4, xVert, yVert);

    rotCutMatrix.Form("rotCutMatrix%02d",iShape);
    xPos = (fModRadius + EMC_TowerYshift[iShape])*cos(2.*EMC_Box_Angle2);
    yPos = (fModRadius + EMC_TowerYshift[iShape])*sin(2.*EMC_Box_Angle2);

    transRot->SetAngles(180., 90., 90.);
    transRot->RotateZ(2.*EMC_Box_Angle2*TMath::RadToDeg());

    TGeoCombiTrans* rotMatrix = new TGeoCombiTrans(xPos, yPos, EMC_TowerZshift[iShape], transRot);
    rotMatrix->SetName(rotCutMatrix.Data()); rotMatrix->RegisterYourself();
    transName.Form("%s-%s:%s",emcModuleShape->GetName(), shapePlateVert->GetName(), rotMatrix->GetName());

// Add modules

    TGeoCompositeShape *shapeMod = new TGeoCompositeShape("shapeMod",transName.Data());
    nameECALModule.Form("emcModule%d",iShape);
    TGeoVolume* emcModule = new TGeoVolume(nameECALModule, shapeMod, pMedGlue);
    listECALModules->Add(emcModule); emcModule->SetVisibility(kFALSE);
    emcModule->SetLineColor(1); 
    emcCrate->AddNode(emcModule, 0, new TGeoTranslation(0., 0., fEMC_Z0 + ZminiP[iShape]));

// Add plates 

    TGeoVolume* framePlate = (TGeoVolume*)listECALSupport->At(iShape);
    plateRot->SetAngles(180., 90., 90.- EMC_Plate_angle);
    plateRot->RotateZ(2.*EMC_Box_Angle2*TMath::RadToDeg());
    TGeoCombiTrans *tr_plate = new TGeoCombiTrans(xPos, yPos,
              fEMC_Z0 + ZminiP[iShape] + EMC_TowerZshift[iShape], plateRot);
    emcCrate->AddNode(framePlate, 0, tr_plate);

// Add glue between modules

   nWallSect = 4;   
   if (iShape == 0) nWallSect = 2; 
   dRotAngle = 2*iShape*EMC_NSector_div2*EMC_Box_Angle1;
   dzWidth = EMC_Box_gap; dLenH =  EMC_Box_Bottom_X*sin(dRotAngle + EMC_Box_Angle1);
   dSideGlue = EMC_Module_side - dLenH/cos(dRotAngle);
   dzPos = fEMC_Z0 + ZminiP[iShape] - dLenH*tan(dRotAngle) - dzWidth/cos(dRotAngle);      

   gPos[0] = dLenH*tan(dRotAngle); gRmin[0] = fModRadius + dLenH; gRmax[0] = gRmin[0]; 
   gPos[1] = gPos[0] + dzWidth/cos(dRotAngle);
   gRmin[1] = gRmin[0]; gRmax[1] = gRmin[0] + dzWidth/sin(dRotAngle);
   gPos[2] = gPos[0] + dSideGlue*sin(dRotAngle); 
   gRmax[2] = gRmin[0] + dSideGlue*cos(dRotAngle);
   gRmin[2] = gRmax[2] - dzWidth/sin(dRotAngle);
   gPos[3] = gPos[2] + dzWidth/cos(dRotAngle);
   gRmin[3] = gRmin[0] + dSideGlue*cos(dRotAngle); gRmax[3] = gRmin[3]; 

   TGeoPgon* emcGlueShape1 = new TGeoPgon(0.,EMC_Module_angle*TMath::RadToDeg(),2,nWallSect);

   if (iShape == 0) { 
     gRmin[0] = fModRadius+ dLenH; gRmax[0] = gRmin[0] + EMC_Module_side;
     emcGlueShape1->DefineSection(0, 0.0, gRmin[0], gRmax[0]); 
     emcGlueShape1->DefineSection(1, dzWidth, gRmin[0], gRmax[0]); 
   }
   if (iShape != 0) {
    for (Int_t iPart = 0; iPart < 4; iPart++)
     emcGlueShape1->DefineSection(iPart, gPos[iPart], gRmin[iPart], gRmax[iPart]);
   }

   TGeoVolume* emcGlueVol1 = new TGeoVolume("emcGlueVol1",emcGlueShape1, pMedGlue);
   emcGlueVol1->SetLineColor(2); emcGlueVol1->SetVisibility(kTRUE);
   emcCrate->AddNode(emcGlueVol1, iShape, new TGeoTranslation(0.,0.,dzPos));
    
   if (iShape == EMC_NSector_div2 - 1) {

    dRotAngle = 2*EMC_NSector_div2*EMC_NSector_div2*EMC_Box_Angle1;
    dzPos = fEMC_Z0 + ZminiP[iShape] + 2.*emcModuleShape->GetDZ() - EMC_Module_side*sin(dRotAngle); 

    gPos[0] = 0.0; gRmin[0] = fModRadius; gRmax[0] = gRmin[0]; 
    gPos[1] = gPos[0] + dzWidth/cos(dRotAngle);
    gRmin[1] = gRmin[0]; gRmax[1] = gRmin[0] + dzWidth/sin(dRotAngle);
    gPos[2] = gPos[0] + EMC_Module_side*sin(dRotAngle); 
    gRmax[2] = gRmin[0] + EMC_Module_side*cos(dRotAngle);
    gRmin[2] = gRmax[2] - dzWidth/sin(dRotAngle);
    gPos[3] = gPos[2] + dzWidth/cos(dRotAngle);
    gRmin[3] = gRmin[0] + EMC_Module_side*cos(dRotAngle); gRmax[3] = gRmin[3]; 

    TGeoPgon* emcGlueShape1 = new TGeoPgon(0.,EMC_Module_angle*TMath::RadToDeg(),2,nWallSect);
    for (Int_t iPart = 0; iPart < 4; iPart++)
     emcGlueShape1->DefineSection(iPart, gPos[iPart], gRmin[iPart], gRmax[iPart]);
    TGeoVolume* emcGlueVol1 = new TGeoVolume("emcGlueVol1",emcGlueShape1, pMedGlue);
    emcGlueVol1->SetLineColor(2); emcGlueVol1->SetVisibility(kTRUE);
    emcCrate->AddNode(emcGlueVol1, iShape+1, new TGeoTranslation(0.,0.,dzPos));

   }
 
   Double_t* pModRmax = emcModuleShape->GetRmax();
   Double_t* pModZ = emcModuleShape->GetZ();   
   Int_t nModIndex1 = (Int_t)TMath::LocMax(emcModuleShape->GetNz(),pModRmax);
   Int_t nModIndex2 = (Int_t)TMath::LocMax(emcModuleShape->GetNz(),pModZ);

   if (iShape == 0) {
    zGlue[0] = crateZ[0]; rGlueMax[0] = fGlueRad; 
    zGlue[iShape+1] = crateZ[0]; 
    rGlueMax[iShape+1] = fGlueRad + pModRmax[nModIndex1] -  crateRadMin[0];
   }

   else {
    zGlue[iShape+1] = fEMC_Z0 + ZminiP[iShape] + pModZ[nModIndex1]; 	
    rGlueMax[iShape+1] = fGlueRad + pModRmax[nModIndex1] - crateRadMin[0];   
   }
   
   if (iShape == EMC_NSector_div2-1) {
    dLenH = EMC_Box_Top_X*sin(dRotAngle - EMC_Box_Angle1);
    zGlue[iShape + 2] = fEMC_Z0 + ZminiP[iShape] + pModZ[nModIndex2] + dLenH*tan(dRotAngle); 
    rGlueMax[iShape + 2] = fGlueRad + pModRmax[nModIndex2] - crateRadMin[0] + dLenH;
    rGlueMax[iShape + 3] = fGlueRad;
    zGlue[iShape + 3] = dzPos; 
   } 

   
 }

// Add glue along modules 

  TGeoXtru *emcGlueShape2 = new TGeoXtru(2);
  emcGlueShape2->DefinePolygon(EMC_NSector_div2 + 3, rGlueMax, zGlue);
  emcGlueShape2->DefineSection(0,-0.5*dWidthGlue,0.,0.,1.0);
  emcGlueShape2->DefineSection(1, 0.5*dWidthGlue,0.,0.,1.0);

  dRotAngle = 0; 
  TGeoVolume *emcGlueVol2 = new TGeoVolume("emcGlueVol2", emcGlueShape2, pMedGlue);
  emcGlueVol2->SetVisibility(kTRUE); emcGlueVol2->SetLineColor(2);

  for (Int_t iGlue = 0; iGlue < EMC_NSector_div1+1; iGlue++) {
   TGeoCombiTrans *glue_trans1 = new TGeoCombiTrans();
   glue_trans1->RotateX(90.); glue_trans1->RotateZ(dRotAngle*TMath::RadToDeg());
   glue_trans1->SetTranslation(xPos0,yPos0,0.);
   emcSector->AddNode(emcGlueVol2, iGlue, glue_trans1);
   dRotAngle += EMC_Sector_angle/EMC_NSector_div1;
  }


}


// Create ECAL module

void CreateECALModule() {

  Int_t arbType, icount = 0, icont = 0, iStep = 0, iMod = 0;
  Double_t fEMC_Z0 = -EMC_Length/2., zTrans[EMC_NSector_div2];
  Double_t xPos, yPos, xTower, yTower, rPos, fModRadius, rotation_angle, tower_angle, arbLen1, arbLen2;
  Double_t SideSize = (EMC_Box_Top_Y - EMC_Box_Bottom_Y)/(2.*tan(EMC_Box_Angle2));

  Double_t crate_radius = ((TGeoPcon*)emcCrate->GetShape())->GetRmin()[0];
  fModRadius = crate_radius + (EMC_Tower_InnerRadius - EMC_Sector_InnerRadius);
  TGeoVolume* volModule;

  TGeoNode* crateNode;
  TObjArray* nodeCrate = emcCrate->GetNodes();
  TIterator* crateIter = nodeCrate->MakeIterator();
   while( (crateNode = (TGeoNode*)crateIter->Next()) ) {
   TString nameECALModule(crateNode->GetName());
    if (nameECALModule.Contains("emcModule")) {
       zTrans[icont] = crateNode->GetMatrix()->GetTranslation()[2] - fEMC_Z0; icont++;
    }
   }

  for (Int_t index = 0; index < EMC_Index_Box; index++) {
   rotation_angle = (2*index+1)*EMC_Box_Angle1; arbType = 0;
   arbLen1 = (SideSize - EMC_Box_Bottom_X*sin(rotation_angle))*cos(EMC_Box_Angle1)/
             cos(rotation_angle-EMC_Box_Angle1);
   arbLen2 = SideSize*cos(EMC_Box_Angle1)/cos(rotation_angle+EMC_Box_Angle1);
   if (arbLen1 > EMC_Box_Length) {arbLen1 = EMC_Box_Length; arbType = 1;}
   if (arbLen1 < EMC_Box_Length) arbType = 2;
   if ( (arbLen1 < EMC_Box_Length) && (arbLen2 < EMC_Box_Length) ) arbType = 3;
   if (arbType == 1) CreateECALTower(index);
   if (arbType == 2) CreateECALCompTower(index, arbLen1);
   if (arbType == 3) CreateECALCompTower(index, arbLen1, arbLen2);
   if (iStep == EMC_NSector_div2) iStep = 0;
   if (iStep == 0) {
    volModule = (TGeoVolume*)listECALModules->At(iMod); iMod++;
   }

   Double_t sizeB = EMC_Box_B[index];
   if (arbType == 2) sizeB = EMC_Box_B[index] + 2.*(EMC_Box_Length - arbLen1)*
                        cos(rotation_angle+EMC_Box_Angle1)*tan(EMC_Box_Angle2);
   if (arbType == 3) sizeB = EMC_Box_Top_Y;

   if (fDebug) printf("MpdRoot : Tower type, EMC_Box_A, B, C %d %d %.03f %.03f %.03f \n", index, arbType,
                                 EMC_Box_A[index], sizeB, EMC_Box_C[index]);

   TGeoVolume* volTower = (TGeoVolume*)listECALTowers->At(index);
   volModule->SetLineColor(42); volTower->SetLineColor(42);
   tower_angle = 2.*EMC_Box_Angle2;
   rPos = 0.5*(EMC_Box_Length*cos(rotation_angle)+EMC_Box_Bottom_X*sin(rotation_angle));
   for (Int_t iXY = 0; iXY < 2; iXY++) {
    xTower = rPos*cos(iXY*tower_angle + EMC_Box_Angle2);
    yTower = rPos*sin(iXY*tower_angle + EMC_Box_Angle2);
    xPos = fModRadius*cos(tower_angle) - (2*iXY-1)*0.5*(EMC_Box_Bottom_Y+EMC_Box_gap)*sin(tower_angle);
    yPos = fModRadius*sin(tower_angle) + (2*iXY-1)*0.5*(EMC_Box_Bottom_Y+EMC_Box_gap)*cos(tower_angle);
    TGeoCombiTrans* posMatrix = new TGeoCombiTrans();
    posMatrix->RotateY(90. - rotation_angle*TMath::RadToDeg());
    posMatrix->RotateZ((iXY*tower_angle + EMC_Box_Angle2)*TMath::RadToDeg());
    posMatrix->SetTranslation(xPos + xTower, yPos + yTower, EMC_zPosition[index] - zTrans[iMod-1]);
    volModule->AddNode(volTower, icount, posMatrix);
    icount++;

   } iStep++;
  }
}
	
// Tower consisting from one trapezoid

void CreateECALTower(Int_t index) {

   TString nameECALTower, nameTowerShape, tr1Name, tr2Name, rot1Name, rot2Name, transName;
   Double_t xTruCut1[4], xTruCut2[4], yTruCut1[4], yTruCut2[4];
   Double_t rotAngle1, rotAngle2, fWidthElemCut, dOffsetY, dOffsetZ;

   TList *arbShapes = new TList();
   TGeoArb8 *arb0 = new TGeoArb8(EMC_Box_Length/2.);
   nameTowerShape.Form("emc_boxShape%d",index+1);
   tr1Name.Form("tr1Name%d",index+1); tr2Name.Form("tr2Name%d",index+1);
   rot1Name.Form("rot1Name%d",index+1); rot2Name.Form("rot2Name%d",index+1);

   arb0->SetName(nameTowerShape.Data());
   arb0->SetVertex(0,-EMC_Box_Bottom_X/2.,-EMC_Box_Bottom_Y/2.);
   arb0->SetVertex(1,-EMC_Box_Bottom_X/2., EMC_Box_Bottom_Y/2.);
   arb0->SetVertex(2, EMC_Box_Bottom_X/2., EMC_Box_C[index]/2.);
   arb0->SetVertex(3, EMC_Box_Bottom_X/2.,-EMC_Box_C[index]/2.);
   arb0->SetVertex(4,-EMC_Box_Top_X/2.,-EMC_Box_B[index]/2.);
   arb0->SetVertex(5,-EMC_Box_Top_X/2., EMC_Box_B[index]/2.);
   arb0->SetVertex(6, EMC_Box_Top_X/2., EMC_Box_A[index]/2.);
   arb0->SetVertex(7, EMC_Box_Top_X/2.,-EMC_Box_A[index]/2.);

// Subtract upper cuts

   TGeoXtru *xElemCut = (TGeoXtru*)listECALTowerCuts->At(index);
   rotAngle1 = atan(0.5*(EMC_Box_A[index] - EMC_Box_C[index])/EMC_Box_Length);
   rotAngle2 = atan(0.5*(EMC_Box_A[index] - EMC_Box_B[index])/EMC_Box_Top_X);

   fWidthElemCut = xElemCut->GetZ(1);

   dOffsetY = 0.5*(EMC_Box_A[index] - fWidthElemCut*cos(rotAngle1)) - 0.5*EMC_Box_Length*sin(rotAngle1);
   dOffsetZ = 0.5*EMC_Box_Length*(1-cos(rotAngle1)) + 0.5*fWidthElemCut*sin(rotAngle1);

   TGeoXtru *xElemCut1 = new TGeoXtru(2); xElemCut1->SetName(tr1Name.Data());
   TGeoXtru *xElemCut2 = new TGeoXtru(2); xElemCut2->SetName(tr2Name.Data());
   for (Int_t iPoint = 0; iPoint < xElemCut->GetNvert(); iPoint++) {
     xTruCut1[iPoint] = xElemCut->GetX(iPoint); yTruCut1[iPoint] = xElemCut->GetY(iPoint);
     xTruCut2[iPoint] = xTruCut1[iPoint]; yTruCut2[iPoint] = yTruCut1[iPoint];
   }

   xTruCut1[2] = 0.5*EMC_Box_Top_Y; xTruCut1[3] = -0.5*EMC_Box_Top_Y;
   xTruCut2[2] = xTruCut1[2]; xTruCut2[3] = xTruCut1[3];
   yTruCut1[2] = yTruCut1[2]  + fPrecision;  yTruCut1[3] =  yTruCut1[3]  + fPrecision;
   yTruCut2[2] = yTruCut1[2];  yTruCut2[3] = yTruCut1[3];

   xElemCut1->DefinePolygon(4, xTruCut1, yTruCut1);
   xElemCut2->DefinePolygon(4, xTruCut2, yTruCut2);
   xElemCut1->DefineSection(0, -0.5*fWidthElemCut*cos(rotAngle2),  0.5*fWidthElemCut*sin(rotAngle2), 0., 1.);
   xElemCut1->DefineSection(1,  0.5*fWidthElemCut*cos(rotAngle2), -0.5*fWidthElemCut*sin(rotAngle2), 0., 1.);
   xElemCut2->DefineSection(0, -0.5*fWidthElemCut*cos(rotAngle2), -0.5*fWidthElemCut*sin(rotAngle2), 0., 1.);
   xElemCut2->DefineSection(1,  0.5*fWidthElemCut*cos(rotAngle2),  0.5*fWidthElemCut*sin(rotAngle2), 0., 1.);

   TGeoCombiTrans *trPlateCut1 = new TGeoCombiTrans(rot1Name.Data());
   trPlateCut1->RotateX(90.+rotAngle1*TMath::RadToDeg());
   trPlateCut1->RotateZ(180.+rotAngle2*TMath::RadToDeg());
   trPlateCut1->SetTranslation(0.0, dOffsetY, dOffsetZ);
   trPlateCut1->RegisterYourself();

   TGeoCombiTrans *trPlateCut2 = new TGeoCombiTrans(rot2Name.Data());
   trPlateCut2->RotateX(90.-rotAngle1*TMath::RadToDeg());
   trPlateCut2->RotateZ(180.-rotAngle2*TMath::RadToDeg());
   trPlateCut2->SetTranslation(0.,-dOffsetY, dOffsetZ);
   trPlateCut2->RegisterYourself();

   transName.Form("%s-%s:%s-%s:%s",nameTowerShape.Data(),
                                xElemCut1->GetName(),trPlateCut1->GetName(),
                                xElemCut2->GetName(),trPlateCut2->GetName());

   TGeoCompositeShape *arb = new TGeoCompositeShape("arb",transName.Data());
   nameECALTower.Form("emc_box%d",index+1);
   arbShapes->AddAt(arb0,0);
   TGeoVolume *emc_tower = new TGeoVolume(nameECALTower, arb, pMedTiPaint);
   CreateECALStructure(index, emc_tower, arbShapes);
   listECALTowers->Add(emc_tower);

}

// Tower consisting from two trapezoids

void CreateECALCompTower(Int_t index, Double_t Length1) {

   Double_t rotation_angle = (2*index+1)*EMC_Box_Angle1;

   TString nameECALTower, arb1Name, arb2Name, rotLower, rotUpper, rot1Name, rot2Name;
   TString transName, tr1Name, tr2Name;
   Double_t rotAngle1, rotAngle2, fWidthElemCut, fRadiusElemCut, dOffsetY, dOffsetZ, dShiftYSize;
   Double_t xTruCut1[4], yTruCut1[4], xTruCut2[4], yTruCut2[4];

   nameECALTower.Form("emc_box%d",index+1);
   Double_t fXYSize = 0.5*(EMC_Box_Bottom_X + 2.*Length1*tan(EMC_Box_Angle1));
   Double_t fBSize = EMC_Box_B[index]/2. + (EMC_Box_Length - Length1)*
                        cos(rotation_angle+EMC_Box_Angle1)*tan(EMC_Box_Angle2);

   TList *arbShapes = new TList();
   TGeoArb8 *arb1 = new TGeoArb8(Length1/2.);
   arb1Name.Form("arbLower%d",index+1); arb1->SetName(arb1Name);
   arb1->SetVertex(0,-EMC_Box_Bottom_X/2.,-EMC_Box_Bottom_Y/2.);
   arb1->SetVertex(1,-EMC_Box_Bottom_X/2., EMC_Box_Bottom_Y/2.);
   arb1->SetVertex(2, EMC_Box_Bottom_X/2., EMC_Box_C[index]/2.);
   arb1->SetVertex(3, EMC_Box_Bottom_X/2.,-EMC_Box_C[index]/2.);
   arb1->SetVertex(4,-fXYSize, -EMC_Box_B[index]/2.);
   arb1->SetVertex(5,-fXYSize, EMC_Box_B[index]/2.);
   arb1->SetVertex(6, fXYSize, EMC_Box_Top_Y/2.);
   arb1->SetVertex(7, fXYSize,-EMC_Box_Top_Y/2.);

   TGeoArb8 *arb2 = new TGeoArb8((EMC_Box_Length - Length1)/2.);
   arb2Name.Form("arbUpper%d",index+1); arb2->SetName(arb2Name);
   arb2->SetVertex(0,-fXYSize, -EMC_Box_B[index]/2.);
   arb2->SetVertex(1,-fXYSize, EMC_Box_B[index]/2.);
   arb2->SetVertex(2, fXYSize, EMC_Box_Top_Y/2.);
   arb2->SetVertex(3, fXYSize,-EMC_Box_Top_Y/2.);
   arb2->SetVertex(4,-EMC_Box_Top_X/2.,-fBSize);
  arb2->SetVertex(5,-EMC_Box_Top_X/2., fBSize);
   arb2->SetVertex(6, EMC_Box_Top_X/2., EMC_Box_Top_Y/2.);
   arb2->SetVertex(7, EMC_Box_Top_X/2.,-EMC_Box_Top_Y/2.);

   rotLower.Form("trLower%d",index+1); rotUpper.Form("trUpper%d",index+1);
   tr1Name.Form("tr1Name%d",index+1); tr2Name.Form("tr2Name%d",index+1);
   rot1Name.Form("rot1Name%d",index+1); rot2Name.Form("rot2Name%d",index+1);

   TGeoTranslation *trArb1 = new TGeoTranslation(rotLower,0.,0.,-(EMC_Box_Length - Length1)/2.);
   trArb1->RegisterYourself();
   TGeoTranslation *trArb2 = new TGeoTranslation(rotUpper,0.,0., Length1/2.);
   trArb2->RegisterYourself();

   TGeoXtru *xElemCut = (TGeoXtru*)listECALTowerCuts->At(index);
   fRadiusElemCut = xElemCut->GetY(0); fWidthElemCut = xElemCut->GetZ(1);

   rotAngle1 = atan(0.5*(EMC_Box_A[index] - EMC_Box_C[index])/Length1);
   rotAngle2 = atan(0.5*(EMC_Box_A[index] - EMC_Box_B[index])/EMC_Box_Top_X);

   TGeoXtru *xElemCut1 = new TGeoXtru(2); xElemCut1->SetName(tr1Name.Data());
   TGeoXtru *xElemCut2 = new TGeoXtru(2); xElemCut2->SetName(tr2Name.Data());
   for (Int_t iPoint = 0; iPoint < xElemCut->GetNvert(); iPoint++) {
     xTruCut1[iPoint] = xElemCut->GetX(iPoint); yTruCut1[iPoint] = xElemCut->GetY(iPoint);
     xTruCut2[iPoint] = xTruCut1[iPoint]; yTruCut2[iPoint] = yTruCut1[iPoint];
   }

   xTruCut1[2] = 0.5*EMC_Box_Top_Y; xTruCut1[3] = -0.5*EMC_Box_Top_Y;
   xTruCut2[2] = xTruCut1[2]; xTruCut2[3] = xTruCut1[3];
   yTruCut1[2] = yTruCut1[2]  + fPrecision;  yTruCut1[3] =  yTruCut1[3]  + fPrecision;
   yTruCut2[2] = yTruCut1[2];  yTruCut2[3] = yTruCut1[3];

   dShiftYSize = (EMC_Box_Length - Length1)*tan(rotAngle1);
   dOffsetY = 0.5*(EMC_Box_Top_Y - fWidthElemCut*cos(rotAngle1)) - 
		   0.5*EMC_Box_Length*sin(rotAngle1) + dShiftYSize + fPrecision;
   dOffsetZ = 0.5*EMC_Box_Length*(1-cos(rotAngle1)) + 0.5*fWidthElemCut*sin(rotAngle1) + fPrecision;

   xElemCut1->DefinePolygon(4, xTruCut1, yTruCut1);
   xElemCut2->DefinePolygon(4, xTruCut2, yTruCut2);
   xElemCut1->DefineSection(0, -0.5*fWidthElemCut*cos(rotAngle2),  0.5*fWidthElemCut*sin(rotAngle2), 0., 1.);
   xElemCut1->DefineSection(1,  0.5*fWidthElemCut*cos(rotAngle2), -0.5*fWidthElemCut*sin(rotAngle2), 0., 1.);
   xElemCut2->DefineSection(0, -0.5*fWidthElemCut*cos(rotAngle2), -0.5*fWidthElemCut*sin(rotAngle2), 0., 1.);
   xElemCut2->DefineSection(1,  0.5*fWidthElemCut*cos(rotAngle2),  0.5*fWidthElemCut*sin(rotAngle2), 0., 1.);

   TGeoCombiTrans *trPlateCut1 = new TGeoCombiTrans(rot1Name.Data());
   trPlateCut1->RotateX(90.+rotAngle1*TMath::RadToDeg());
   trPlateCut1->RotateZ(180.+rotAngle2*TMath::RadToDeg());
   trPlateCut1->SetTranslation(0.0, dOffsetY, dOffsetZ);
   trPlateCut1->RegisterYourself();

   TGeoCombiTrans *trPlateCut2 = new TGeoCombiTrans(rot2Name.Data());
   trPlateCut2->RotateX(90.-rotAngle1*TMath::RadToDeg());
   trPlateCut2->RotateZ(180.-rotAngle2*TMath::RadToDeg());
   trPlateCut2->SetTranslation(0.,-dOffsetY, dOffsetZ);
   trPlateCut2->RegisterYourself();

   transName.Form("%s:%s+%s:%s-%s:%s-%s:%s",
                arb1Name.Data(), trArb1->GetName(),
                arb2Name.Data(), trArb2->GetName(),
                xElemCut1->GetName(), trPlateCut1->GetName(),
                xElemCut1->GetName(), trPlateCut2->GetName());

   TGeoCompositeShape *arb = new TGeoCompositeShape("arb",transName);
   arbShapes->AddAt(arb1, 0); arbShapes->AddAt(arb2, 1);
   TGeoVolume *emc_tower = new TGeoVolume(nameECALTower, arb, pMedTiPaint);

   CreateECALStructure(index, emc_tower, arbShapes);
   listECALTowers->Add(emc_tower);

}

// Tower consisting from three trapezoids

void CreateECALCompTower(Int_t index, Double_t Length1, Double_t Length2) {

   TString nameECALTower, arb1Name, arb2Name, arb3Name, rotLower, rotMiddle, rotUpper;
   TString transName, rot1Name, rot2Name, tr1Name, tr2Name;
   Double_t rotAngle1, rotAngle2, fWidthElemCut, fRadiusElemCut, dOffsetY, dOffsetZ, dShiftYSize;
   Double_t xTruCut1[4], yTruCut1[4], xTruCut2[4], yTruCut2[4];

   Double_t rotation_angle = (2*index+1)*EMC_Box_Angle1;
   nameECALTower.Form("emc_box%d",index+1);
   Double_t fXYSize1 = 0.5*(EMC_Box_Bottom_X + 2.*Length1*tan(EMC_Box_Angle1));
   Double_t fXYSize2 = 0.5*(EMC_Box_Bottom_X + 2.*Length2*tan(EMC_Box_Angle1));

   TList *arbShapes = new TList();
   TGeoArb8 *arb1 = new TGeoArb8(Length1/2.);
   arb1Name.Form("arbLower%d",index+1); arb1->SetName(arb1Name);
   arb1->SetVertex(0,-EMC_Box_Bottom_X/2.,-EMC_Box_Bottom_Y/2.);
   arb1->SetVertex(1,-EMC_Box_Bottom_X/2., EMC_Box_Bottom_Y/2.);
   arb1->SetVertex(2, EMC_Box_Bottom_X/2., EMC_Box_C[index]/2.);
   arb1->SetVertex(3, EMC_Box_Bottom_X/2.,-EMC_Box_C[index]/2.);
   arb1->SetVertex(4,-fXYSize1,-EMC_Box_B[index]/2);
   arb1->SetVertex(5,-fXYSize1, EMC_Box_B[index]/2.);
   arb1->SetVertex(6, fXYSize1, EMC_Box_Top_Y/2.);
   arb1->SetVertex(7, fXYSize1,-EMC_Box_Top_Y/2.);

   TGeoArb8 *arb2 = new TGeoArb8((Length2 - Length1)/2.);
   arb2Name.Form("arbMiddle%d",index+1); arb2->SetName(arb2Name);
   arb2->SetVertex(0,-fXYSize1,-EMC_Box_B[index]/2);
   arb2->SetVertex(1,-fXYSize1, EMC_Box_B[index]/2.);
   arb2->SetVertex(2, fXYSize1, EMC_Box_Top_Y/2.);
   arb2->SetVertex(3, fXYSize1,-EMC_Box_Top_Y/2.);
   arb2->SetVertex(4,-fXYSize2,-EMC_Box_Top_Y/2.);
   arb2->SetVertex(5,-fXYSize2, EMC_Box_Top_Y/2.);
   arb2->SetVertex(6, fXYSize2, EMC_Box_Top_Y/2.);
   arb2->SetVertex(7, fXYSize2,-EMC_Box_Top_Y/2.);

   TGeoArb8 *arb3 = new TGeoArb8((EMC_Box_Length - Length2)/2.);
   arb3Name.Form("arbUpper%d",index+1); arb3->SetName(arb3Name);
   arb3->SetVertex(0,-fXYSize2,-EMC_Box_Top_Y/2.);
   arb3->SetVertex(1,-fXYSize2, EMC_Box_Top_Y/2.);
   arb3->SetVertex(2, fXYSize2, EMC_Box_Top_Y/2.);
   arb3->SetVertex(3, fXYSize2,-EMC_Box_Top_Y/2.);
   arb3->SetVertex(4,-EMC_Box_Top_X/2.,-EMC_Box_Top_Y/2.);
   arb3->SetVertex(5,-EMC_Box_Top_X/2., EMC_Box_Top_Y/2.);
   arb3->SetVertex(6, EMC_Box_Top_X/2., EMC_Box_Top_Y/2.);
   arb3->SetVertex(7, EMC_Box_Top_X/2.,-EMC_Box_Top_Y/2.);

   rotLower.Form("trLower%d",index+1); rotMiddle.Form("trMiddle%d",index+1); rotUpper.Form("trUpper%d",index+1);
   tr1Name.Form("tr1Name%d",index+1); tr2Name.Form("tr2Name%d",index+1);
   rot1Name.Form("rot1Name%d",index+1); rot2Name.Form("rot2Name%d",index+1);

   TGeoTranslation *trArb1 = new TGeoTranslation(rotLower,0.,0.,-(EMC_Box_Length - Length1)/2.);
   trArb1->RegisterYourself();
   TGeoTranslation *trArb2 = new TGeoTranslation(rotMiddle,0.,0.,-(EMC_Box_Length - Length2  - Length1)/2.);
   trArb2->RegisterYourself();
   TGeoTranslation *trArb3 = new TGeoTranslation(rotUpper,0.,0., Length2/2.);
   trArb3->RegisterYourself();

   TGeoXtru *xElemCut = (TGeoXtru*)listECALTowerCuts->At(index);
   fRadiusElemCut = xElemCut->GetY(0); fWidthElemCut = xElemCut->GetZ(1);

   rotAngle1 = atan(0.5*(EMC_Box_A[index] - EMC_Box_C[index])/Length1);
   rotAngle2 = atan(0.5*(EMC_Box_A[index] - EMC_Box_B[index])/EMC_Box_Top_X);

   TGeoXtru *xElemCut1 = new TGeoXtru(2); xElemCut1->SetName(tr1Name.Data());
   TGeoXtru *xElemCut2 = new TGeoXtru(2); xElemCut2->SetName(tr2Name.Data());
   for (Int_t iPoint = 0; iPoint < xElemCut->GetNvert(); iPoint++) {
     xTruCut1[iPoint] = xElemCut->GetX(iPoint); yTruCut1[iPoint] = xElemCut->GetY(iPoint);
     xTruCut2[iPoint] = xTruCut1[iPoint]; yTruCut2[iPoint] = yTruCut1[iPoint];
   }

   xTruCut1[2] = 0.5*EMC_Box_Top_Y; xTruCut1[3] = -0.5*EMC_Box_Top_Y;
   xTruCut2[2] = xTruCut1[2]; xTruCut2[3] = xTruCut1[3];
   yTruCut1[2] = yTruCut1[2]  + fPrecision;  yTruCut1[3] =  yTruCut1[3]  + fPrecision;
   yTruCut2[2] = yTruCut1[2];  yTruCut2[3] = yTruCut1[3];

   dShiftYSize = (EMC_Box_Length - Length1)*tan(rotAngle1);
   dOffsetY = 0.5*(EMC_Box_Top_Y - fWidthElemCut*cos(rotAngle1)) - 
	      0.5*EMC_Box_Length*sin(rotAngle1) + dShiftYSize + fPrecision;
   dOffsetZ = 0.5*EMC_Box_Length*(1-cos(rotAngle1)) + 0.5*fWidthElemCut*sin(rotAngle1) + fPrecision;

   xElemCut1->DefinePolygon(4, xTruCut1, yTruCut1);
   xElemCut2->DefinePolygon(4, xTruCut2, yTruCut2);
   xElemCut1->DefineSection(0, -0.5*fWidthElemCut*cos(rotAngle2),  0.5*fWidthElemCut*sin(rotAngle2), 0., 1.);
   xElemCut1->DefineSection(1,  0.5*fWidthElemCut*cos(rotAngle2), -0.5*fWidthElemCut*sin(rotAngle2), 0., 1.);
   xElemCut2->DefineSection(0, -0.5*fWidthElemCut*cos(rotAngle2), -0.5*fWidthElemCut*sin(rotAngle2), 0., 1.);
   xElemCut2->DefineSection(1,  0.5*fWidthElemCut*cos(rotAngle2),  0.5*fWidthElemCut*sin(rotAngle2), 0., 1.);

   TGeoCombiTrans *trPlateCut1 = new TGeoCombiTrans(rot1Name.Data());
   trPlateCut1->RotateX(90.+rotAngle1*TMath::RadToDeg());
   trPlateCut1->RotateZ(180.+rotAngle2*TMath::RadToDeg());
   trPlateCut1->SetTranslation(0.0, dOffsetY, dOffsetZ);
   trPlateCut1->RegisterYourself();

   TGeoCombiTrans *trPlateCut2 = new TGeoCombiTrans(rot2Name.Data());
   trPlateCut2->RotateX(90.-rotAngle1*TMath::RadToDeg());
   trPlateCut2->RotateZ(180.-rotAngle2*TMath::RadToDeg());
   trPlateCut2->SetTranslation(0.,-dOffsetY, dOffsetZ);
   trPlateCut2->RegisterYourself();

   transName.Form("%s:%s+%s:%s+%s:%s-%s:%s-%s:%s",
                arb1Name.Data(), trArb1->GetName(),
                arb2Name.Data(), trArb2->GetName(),
                arb3Name.Data(), trArb3->GetName(),
                xElemCut1->GetName(), trPlateCut1->GetName(),
                xElemCut2->GetName(), trPlateCut2->GetName());

   TGeoCompositeShape *arbComp = new TGeoCompositeShape("arb",transName);
   arbShapes->AddAt(arb1, 0); arbShapes->AddAt(arb2, 1); arbShapes->AddAt(arb3, 2);
   TGeoVolume *emc_tower = new TGeoVolume(nameECALTower, arbComp, pMedTiPaint);
   CreateECALStructure(index, emc_tower, arbShapes);
   listECALTowers->Add(emc_tower);

}

// Create ECAL module structure

void CreateECALStructure(Int_t index, TGeoVolume* topECALModule, TList* arbShapes) {

   TGeoArb8* topECALArb1;
   TString nameBoxSc, nameBoxPb;
   Double_t zPos0 = -EMC_Box_Length/2., xLow, xHigh, fRadiusElemCut;
   Double_t zPosSc, zPosPb, zPosPl, zPosAdd, arbLen1, arbLen2, arbLen3;
   Double_t points_low[8], points_high[8], EMC_Box_Z = 0;
   Int_t iCounter = 0, inext1 = 0, inext2 = 0;

   TGeoXtru *xElemCut = (TGeoXtru*)listECALTowerCuts->At(index);
   fRadiusElemCut = xElemCut->GetY(0);
   if (fRadiusElemCut > xElemCut->GetY(1)) fRadiusElemCut = xElemCut->GetY(1);

   if (arbShapes->GetSize() > 0) {
     topECALArb1 = (TGeoArb8*)arbShapes->At(0);
     arbLen1 = 2.*topECALArb1->GetDz();
   }

// Add first fixing plate

    TGeoArb8 *PlPlate1 = new TGeoArb8(EMC_Box_Plastic1/2.);
    zPosPl = zPos0 + EMC_Box_Plastic1/2.;

    topECALArb1->SetPlaneVertices(-arbLen1/2., points_low);
    topECALArb1->SetPlaneVertices(-arbLen1/2. + EMC_Box_Plastic1, points_high);
    for (int ivert = 0; ivert < 4; ivert++)
           PlPlate1->SetVertex(ivert, points_low[2*ivert], points_low[2*ivert+1]);
    for (int ivert = 0; ivert < 4; ivert++)
           PlPlate1->SetVertex(ivert+4, points_high[2*ivert], points_high[2*ivert+1]);
    TGeoVolume* emc_cl_pl1 = new TGeoVolume("emc_cl_pl1",PlPlate1,pMedKapton);
    emc_cl_pl1->SetLineColor(1); emc_cl_pl1->SetTransparency(0);
    topECALModule->AddNode(emc_cl_pl1, 0, new TGeoTranslation(0., 0., zPosPl));

    iCounter++; EMC_Box_Z += EMC_Box_Plastic1;
    zPos0 = zPos0 + EMC_Box_Plastic1;

   for (int ik = 0; ik <  EMC_Box_Layers; ik++) {

//  Scintillator box

      zPosSc = zPos0 + EMC_Box_Sc/2. + ik*EMC_Box_Cell;
      TGeoArb8 *ScPlate = new TGeoArb8(EMC_Box_Sc/2.);
      xLow = zPosSc - EMC_Box_Sc/2.; xHigh = zPosSc + EMC_Box_Sc/2.;
      CreateECALTowerEdges(arbShapes, xLow, xHigh, points_low, points_high);
      for (int ivert = 0; ivert < 4; ivert++) {
           ScPlate->SetVertex(ivert, points_low[2*ivert], points_low[2*ivert+1]);
      }

      for (int ivert = 0; ivert < 4; ivert++) {
           ScPlate->SetVertex(ivert+4, points_high[2*ivert], points_high[2*ivert+1]);

      }

      nameBoxSc.Form("emc_cl_sc%d",ik+1);
      TGeoVolume* emc_cl_sc = new TGeoVolume(nameBoxSc,ScPlate,pMedFscScint);
      emc_cl_sc->SetLineColor(kGreen);

      if (xHigh > fRadiusElemCut)
         topECALModule->AddNodeOverlap(emc_cl_sc, iCounter, new TGeoTranslation(0., 0., zPosSc));
      else
         topECALModule->AddNode(emc_cl_sc, iCounter, new TGeoTranslation(0., 0., zPosSc));
      iCounter++; EMC_Box_Z += EMC_Box_Sc;

// Pb box

      zPosPb = zPosSc + (EMC_Box_Sc + EMC_Box_Pb)/2. + EMC_Box_Col;
      TGeoArb8 *PbPlate = new TGeoArb8(EMC_Box_Pb/2.);
      xLow = zPosPb - EMC_Box_Pb/2.; xHigh = zPosPb + EMC_Box_Pb/2.;
      CreateECALTowerEdges(arbShapes, xLow, xHigh, points_low, points_high);

      for (int ivert = 0; ivert < 4; ivert++)  
           PbPlate->SetVertex(ivert, points_low[2*ivert], points_low[2*ivert+1]);

      for (int ivert = 0; ivert < 4; ivert++)
           PbPlate->SetVertex(ivert+4, points_high[2*ivert], points_high[2*ivert+1]);
      nameBoxPb.Form("emc_cl_pb%d",ik+1);
      TGeoVolume* emc_cl_pb = new TGeoVolume(nameBoxPb,PbPlate,pMedLead);
      emc_cl_pb->SetLineColor(kGray);

      if (xHigh > fRadiusElemCut)
        topECALModule->AddNodeOverlap(emc_cl_pb, iCounter, new TGeoTranslation(0., 0., zPosPb));
      else
        topECALModule->AddNode(emc_cl_pb, iCounter, new TGeoTranslation(0., 0., zPosPb));
      iCounter++; EMC_Box_Z += EMC_Box_Pb;

   }

//  Add last scintillator box

    zPosSc = zPos0 + (EMC_Box_Sc)/2. + EMC_Box_Layers*EMC_Box_Cell;
    TGeoArb8 *ScPlate = new TGeoArb8(EMC_Box_Sc/2.);
        xLow = zPosSc - EMC_Box_Sc/2.; xHigh = zPosSc + EMC_Box_Sc/2.;
    CreateECALTowerEdges(arbShapes, xLow, xHigh, points_low, points_high);
    for (int ivert = 0; ivert < 4; ivert++)
        ScPlate->SetVertex(ivert, points_low[2*ivert], points_low[2*ivert+1]);
    for (int ivert = 0; ivert < 4; ivert++)
        ScPlate->SetVertex(ivert+4, points_high[2*ivert], points_high[2*ivert+1]);
    nameBoxSc.Form("emc_cl_sc%d",EMC_Box_Layers+1);
    TGeoVolume* emc_cl_sc = new TGeoVolume(nameBoxSc,ScPlate,pMedFscScint);
    emc_cl_sc->SetLineColor(kGreen);
    topECALModule->AddNodeOverlap(emc_cl_sc, iCounter, new TGeoTranslation(0., 0., zPosSc));
    iCounter++; EMC_Box_Z += EMC_Box_Sc;

// Add second plastic

    TGeoArb8 *PlPlate2 = new TGeoArb8(EMC_Box_Plastic2/2.);
    zPosPl = zPosSc + (EMC_Box_Sc + EMC_Box_Plastic2)/2.;
    xLow = zPosPl - EMC_Box_Plastic2/2.; xHigh = zPosPl + EMC_Box_Plastic2/2.;
    CreateECALTowerEdges(arbShapes, xLow, xHigh, points_low, points_high);
    for (int ivert = 0; ivert < 4; ivert++)
           PlPlate2->SetVertex(ivert, points_low[2*ivert], points_low[2*ivert+1]);
     for (int ivert = 0; ivert < 4; ivert++)
           PlPlate2->SetVertex(ivert+4, points_high[2*ivert], points_high[2*ivert+1]);

     TGeoVolume* emc_cl_pl2 = new TGeoVolume("emc_cl_pl2",PlPlate2,pMedKapton);
     emc_cl_pl2->SetLineColor(1); emc_cl_pl2->SetTransparency(0);
     topECALModule->AddNodeOverlap(emc_cl_pl2, iCounter, new TGeoTranslation(0., 0., zPosPl));
     iCounter++; EMC_Box_Z += EMC_Box_Plastic2;

// Addition air space in tower

/*
    Double_t addLength = EMC_Box_Length - EMC_Box_Z;
    TGeoArb8 *PlAdd = new TGeoArb8(addLength/2.);
    zPosAdd = zPosPl + (EMC_Box_Plastic + addLength)/2.;
        xLow = zPosAdd - addLength/2.; xHigh = zPosAdd + addLength/2.;
    CreateECALTowerEdges(arbShapes, xLow, xHigh, points_low, points_high);
    for (int ivert = 0; ivert < 4; ivert++)
           PlAdd->SetVertex(ivert, points_low[2*ivert], points_low[2*ivert+1]);
     for (int ivert = 0; ivert < 4; ivert++)
           PlAdd->SetVertex(ivert+4, points_high[2*ivert], points_high[2*ivert+1]);
     TGeoVolume* emc_cl_air = new TGeoVolume("emc_cl_air",PlAdd,pMedAir);
     emc_cl_air->SetLineColor(1); emc_cl_air->SetTransparency(0);
     topECALModule->AddNode(emc_cl_air, iCounter, new TGeoTranslation(0., 0., zPosAdd));
     iCounter++; EMC_Box_Z += addLength;
*/

}

// Create low and high positions of the trapezoid intersection

void CreateECALTowerEdges(TList* arbShapes, Double_t xLow, Double_t xHigh,
                          Double_t* point_low, Double_t* point_high) {

   Double_t zPosArb0, zPosArb1, zPosArb2, dSize;
   Double_t arbLen, arbLen1, arbLen2, arbLen3;
   TGeoArb8 *topECALArb1, *topECALArb2, *topECALArb3;

   dSize = xHigh - xLow; arbLen = xLow + EMC_Box_Length/2.;
   if (arbLen > 0) arbLen = 0.5*EMC_Box_Length + xLow;

// One arb trapezoid

   if (arbShapes->GetSize() == 1) {
     topECALArb1 = (TGeoArb8*)arbShapes->At(0);
     topECALArb1->SetPlaneVertices(xLow, point_low);
     topECALArb1->SetPlaneVertices(xHigh, point_high);
   }

// Two arbs trapezoid

   if (arbShapes->GetSize() == 2) {
     topECALArb1 = (TGeoArb8*)arbShapes->At(0);
     topECALArb2 = (TGeoArb8*)arbShapes->At(1);
     arbLen1 = 2.*topECALArb1->GetDz();
     arbLen2 = 2.*topECALArb2->GetDz();

     if (arbLen < arbLen1 - dSize) {
       topECALArb1->SetPlaneVertices(xLow + (EMC_Box_Length - arbLen1)/2., point_low);
       topECALArb1->SetPlaneVertices(xHigh + (EMC_Box_Length - arbLen1)/2., point_high);
     }
     else if ( (arbLen > arbLen1 - dSize) && (arbLen < arbLen1) ) {
       topECALArb1->SetPlaneVertices(xLow + (EMC_Box_Length - arbLen1)/2., point_low);
       topECALArb2->SetPlaneVertices(xHigh + EMC_Box_Length/2. - arbLen1 - arbLen2/2., point_high);
     }
     else {
       topECALArb2->SetPlaneVertices(xLow + EMC_Box_Length/2. - arbLen1 - arbLen2/2., point_low);
       topECALArb2->SetPlaneVertices(xHigh + EMC_Box_Length/2. - arbLen1 - arbLen2/2., point_high);
     }

   }

// Three arbs trapezoid

   if (arbShapes->GetSize() == 3) {
     topECALArb1 = (TGeoArb8*)arbShapes->At(0);
     topECALArb2 = (TGeoArb8*)arbShapes->At(1);
     topECALArb3 = (TGeoArb8*)arbShapes->At(2);
     arbLen1 = 2.*topECALArb1->GetDz(); zPosArb0 = -arbLen1/2.;
     arbLen2 = 2.*topECALArb2->GetDz(); zPosArb1 = -arbLen2/2.;
     arbLen3 = 2.*topECALArb3->GetDz(); zPosArb2 = -arbLen3/2.;

     if (arbLen < arbLen1 - dSize) {
       topECALArb1->SetPlaneVertices(xLow + (EMC_Box_Length - arbLen1)/2., point_low);
       topECALArb1->SetPlaneVertices(xHigh + (EMC_Box_Length - arbLen1)/2., point_high);
     }
     else if ( (arbLen > arbLen1 - dSize) && (arbLen < arbLen1) ) {
       topECALArb1->SetPlaneVertices(xLow + (EMC_Box_Length - arbLen1)/2., point_low);
     if ( (arbLen + dSize) <  arbLen1 + arbLen2)
       topECALArb2->SetPlaneVertices(xHigh + EMC_Box_Length/2. - arbLen1 - arbLen2/2., point_high);
     else
	topECALArb3->SetPlaneVertices(xHigh + EMC_Box_Length/2. - arbLen1 - arbLen2 - arbLen3/2., point_high);
     }
     else if ( (arbLen > arbLen1) && (arbLen < arbLen1 + arbLen2 - dSize) ) {
       topECALArb2->SetPlaneVertices(xLow + EMC_Box_Length/2. - arbLen1 - arbLen2/2., point_low);
       topECALArb2->SetPlaneVertices(xHigh + EMC_Box_Length/2. - arbLen1 - arbLen2/2., point_high);
     }
     else if ( (arbLen > arbLen1 + arbLen2 - dSize) && (arbLen < arbLen1 + arbLen2) ) {
       topECALArb2->SetPlaneVertices(xLow + EMC_Box_Length/2. - arbLen1 - arbLen2/2., point_low);
       topECALArb3->SetPlaneVertices(xHigh + EMC_Box_Length/2. - arbLen1 - arbLen2 - arbLen3/2., point_high);
    }
     else {
       topECALArb3->SetPlaneVertices(xLow + EMC_Box_Length/2. - arbLen1 - arbLen2 - arbLen3/2., point_low);
       topECALArb3->SetPlaneVertices(xHigh + EMC_Box_Length/2. - arbLen1 - arbLen2-arbLen3/2., point_high);
     }   
   }

}


void CreateECALModuleSupport() {

// Dural plate parameters

  const Double_t EMC_TowerCut[EMC_NumberPlates] =
        {0.13, 0.13, 0.14, 0.16, 0.16, 0.16, 0.16, 0.16};
  const Double_t EMC_TowerUpperLen[EMC_NumberPlates] =
        {35.00, 35.20, 34.00, 32.20, 30.21, 27.70, 24.20, 20.60};
  const Double_t EMC_TowerLowerLen[EMC_NumberPlates] =
        {34.13, 34.32, 33.09, 31.26, 29.20, 26.60, 22.98, 19.19};

  const Int_t nPoints = 6;
  const Double_t EMC_WidthElem = 0.3;
  const Double_t EMC_ModuleSupprt[EMC_NumberPlates][nPoints] =
   {-0.32, 32.20, 0.8, 4.00, 11.98, 0.0,  1.1, 33.27, 1.7,  3.30, 10.90, 0.0,
      2.6, 34.75, 2.4, 4.80, 11.22, 3.4,  3.2, 37.40, 3.6,  5.20, 12.11, 3.32,
      4.6, 40.77, 4.7, 6.40, 13.17, 4.3,  6.6, 45.15, 5.8,  8.15, 14.68, 5.1,
     11.2, 52.60, 8.8, 9.40, 17.12, 4.5, 16.8, 59.84, 9.4, 13.20, 19.55, 4.85};

  const Double_t EMC_ModuleSupRot[EMC_NumberPlates][2] =
    {47.02, 46.12, 46.11, 45.21, 45.18, 44.28, 44.25, 43.30, 43.27, 42.26,
     42.23, 41.15, 41.04, 39.88, 39.71, 38.51};

// Plate holes

  const Double_t EMC_ModuleSupprtRad1 = 0.5*1.02;
  const Double_t EMC_ModuleSupprtDist1[EMC_NumberPlates] = {4.0, 4.0, 4.0, 4.4, 4.6, 4.9, 5.5, 6.4};
  const Double_t EMC_ModuleSupprtX1[EMC_NumberPlates] = {2.0, 3.0, 3.8, 3.8, 5.0, 6.3, 8.8, 8.6};
  const Double_t EMC_ModuleSupprtY1[EMC_NumberPlates] = {8.5, 7.9, 8.6, 9.0, 10.6, 11.8, 14.4, 16.6};
  const Double_t EMC_ModuleSupprtRad2 = 0.5*2.4;
  const Double_t EMC_ModuleSupprtX2[EMC_NumberPlates*EMC_NumberPlates] =
        {30.5, 25.9, 22.5, 17.9, 14.5,  9.9, 6.3, 1.7,  30.5, 25.9, 22.5, 17.9, 14.5,  9.9, 6.3, 1.7,
         30.5, 25.9, 22.5, 17.9, 14.5,  9.9, 6.3, 1.7,  32.4, 27.8, 24.0, 19.2, 15.2, 10.4, 6.4, 1.8,
         34.3, 29.7, 25.5, 20.9, 16.3, 11.7, 6.7, 2.1,  36.6, 31.8, 27.0, 22.2, 17.2, 12.6, 7.0, 2.2,
         38.8, 34.2, 29.2, 24.2, 18.6, 13.6, 7.8, 2.8,  40.2, 35.4, 30.4, 25.2, 20.2, 14.6, 9.0, 3.0};
  const Double_t EMC_ModuleSupprtY2[EMC_NumberPlates] = {3.2, 3.3, 3.4, 3.3, 4.3, 4.6, 5.0, 4.6};

// Upper plate to fix electronic part

  const Double_t EMC_ModuleUpperX[EMC_NumberPlates] = {32.2, 32.2, 32.2, 34.2, 36.2, 38.6, 41.4, 43.0};
  const Double_t EMC_ModuleUpperZ = EMC_NSector_div2;

  Int_t nPlate = 0;
  Double_t xTruShift[EMC_NumberPlates], yTruShift[EMC_NumberPlates];
  Double_t xTru[6], yTru[6], xTruCut[6], yTruCut[6], xTruVCut[4], yTruVCut[4], fShapeOrigin[3];
  Double_t dPlateLen, dAngleRotate = 0, dLenAdd, dStartNext, rotTower1, rotTower2, dCutHeight;
  Double_t dStartPos, dStartX, dLenSize, holePosX, holePosY, rotation_angle = 0;
  TString nameShape1, nameShape2, nameTrans, nameElemComp, nameElemVert;

// Create dural plate

  xTru[0] = 0.0; yTru[0] = 0.0; yTru[1] = 0.0;
  xTruCut[0] = 0.0; yTruCut[0] = 0.0;

// Estimate mean plate rotation angle

  for (Int_t iPlate = 0; iPlate < EMC_NumberPlates; iPlate++) {
   dPlateLen = EMC_ModuleSupprt[iPlate][1] - EMC_ModuleSupprt[iPlate][0];
   if (iPlate == 0) dPlateLen = EMC_ModuleSupprt[iPlate][1];
   dAngleRotate += TMath::ASin( (EMC_ModuleSupRot[iPlate][0]-EMC_ModuleSupRot[iPlate][1])/
                                dPlateLen )*TMath::RadToDeg();
  }

  dAngleRotate = ((Int_t)(dAngleRotate/EMC_NumberPlates*100))/100.*TMath::DegToRad();
  EMC_Plate_angle =  dAngleRotate*TMath::RadToDeg();

  printf("MpdRoot : Milling angle of the dural plates, grad : %f \n", dAngleRotate*TMath::RadToDeg());

  for (Int_t iPlate = 0; iPlate < EMC_NumberPlates; iPlate++) {

   nameShape1.Form("xElemShape%02d",iPlate);
   nameShape2.Form("xUpperElemShape%02d",iPlate);
   nameElemComp.Form("Plate%01d",iPlate);
   nameElemVert.Form("xElemVert%02d", iPlate);
   nameTrans.Form("%s+%s",nameShape1.Data(),nameShape2.Data());

   TGeoXtru *xElemShape = new TGeoXtru(2); xElemShape->SetName(nameShape1.Data());
   xTru[1] = EMC_ModuleSupprt[iPlate][1] - EMC_ModuleSupprt[iPlate][2];
   xTru[2] = EMC_ModuleSupprt[iPlate][1];
   yTru[2] = EMC_ModuleSupprt[iPlate][4] - EMC_ModuleSupprt[iPlate][3];
   xTru[3] = EMC_ModuleSupprt[iPlate][1];
   yTru[3] = EMC_ModuleSupprt[iPlate][4];
   xTru[4] = EMC_ModuleSupprt[iPlate][0];
   yTru[4] = EMC_ModuleSupprt[iPlate][4];
   xTru[5] = EMC_ModuleSupprt[iPlate][0];
   yTru[5] = EMC_ModuleSupprt[iPlate][4] - EMC_ModuleSupprt[iPlate][5];

   if (iPlate == 0) {
     xTru[1] = xTru[1] + EMC_ModuleSupprt[0][0];
     xTru[2] = xTru[2] + EMC_ModuleSupprt[0][0];
     xTru[3] = xTru[2];
   }

   if (iPlate < 3)
        xElemShape->DefinePolygon(5, xTru, yTru);
   else
        xElemShape->DefinePolygon(6, xTru, yTru);
    xElemShape->DefineSection(0, -0.5*EMC_WidthElem, 0., 0.);
    xElemShape->DefineSection(1,  0.5*EMC_WidthElem, 0., 0.);

   Double_t rotTower1 = 2*EMC_NumberPlates*iPlate*EMC_Box_Angle1;
   Double_t rotTower2 = 2*EMC_NumberPlates*(iPlate+1)*EMC_Box_Angle1;

   if (iPlate < 3) {xTruVCut[3] = xTru[4]; yTruVCut[3] = yTru[4];}

   TGeoXtru *xElemVertCut = new TGeoXtru(2);
   xElemVertCut->SetName(nameElemVert.Data());

   dPlateLen = xTru[1] - xTru[0];
   dCutHeight = sqrt(pow(yTru[2] - yTru[1],2) + pow(xTru[2] - xTru[1],2));

   xTruVCut[0] = 0.0; yTruVCut[0] = 0.0;
   xTruVCut[1] = xTruVCut[0] + dPlateLen*cos(dAngleRotate);
   yTruVCut[1] = yTruVCut[0] - dPlateLen*sin(dAngleRotate);
   xTruVCut[2] = xTruVCut[1] + dCutHeight*sin(rotTower2);
   yTruVCut[2] = yTruVCut[1] + dCutHeight*cos(rotTower2);
   xTruVCut[3] = xTruVCut[0] + dCutHeight*sin(rotTower1);
   yTruVCut[3] = yTruVCut[0] + dCutHeight*cos(rotTower1);
   xElemVertCut->DefinePolygon(4, xTruVCut, yTruVCut);
   xElemVertCut->DefineSection(0, xElemShape->GetZ(0), 0., 0.);
   xElemVertCut->DefineSection(1, xElemShape->GetZ(1), 0., 0.);
   listECALPlateVert->Add(xElemVertCut);

// Add upper plate

   fShapeOrigin[0] = 0.5*EMC_ModuleUpperX[iPlate] + EMC_ModuleSupprt[iPlate][0];
   fShapeOrigin[1] = EMC_ModuleSupprt[iPlate][4]+0.5*EMC_WidthElem;
   fShapeOrigin[2] = 0.0;
   TGeoBBox *xUpperElemShape = new TGeoBBox(nameShape2.Data(), 0.5*EMC_ModuleUpperX[iPlate],
             0.5*EMC_WidthElem, 0.5*EMC_ModuleUpperZ, fShapeOrigin);
   TGeoCompositeShape *emcSupportShape = new TGeoCompositeShape(nameElemComp.Data(), nameTrans.Data());
   TGeoVolume *emcElemSupportComp = new TGeoVolume(nameElemComp.Data(), emcSupportShape, pMedDural);
   emcElemSupportComp->SetVisibility(kTRUE); 
   emcElemSupportComp->SetVisContainers(kTRUE);
   emcElemSupportComp->SetVisLeaves(kTRUE);
   emcElemSupportComp->SetLineColor(4);

// Add holes to the dural plate

   TGeoTube *emcElemHoleShape1 = new TGeoTube(0., EMC_ModuleSupprtRad1, 0.5*EMC_WidthElem);
   TGeoVolume *emcElemSupportHole1 = new TGeoVolume("emcElemSupportHole1", emcElemHoleShape1, pMedGlue);
   emcElemSupportHole1->SetLineColor(12);

   TGeoTube *emcElemHoleShape2 = new TGeoTube(0., EMC_ModuleSupprtRad2, 0.5*EMC_WidthElem);
   TGeoVolume *emcElemSupportHole2 = new TGeoVolume("emcElemSupportHole2", emcElemHoleShape2, pMedAir);
   emcElemSupportHole2->SetLineColor(12);

   for (Int_t iHole = 0; iHole < EMC_NumberPlates; iHole++) {
    holePosX = EMC_ModuleSupprt[iPlate][1] - EMC_ModuleSupprtX1[iPlate] +
                (iHole + 1 - EMC_NumberPlates)*EMC_ModuleSupprtDist1[iPlate];
    holePosY = EMC_ModuleSupprt[iPlate][4] - EMC_ModuleSupprtY1[iPlate];
    if (iPlate == 0) holePosX = holePosX + EMC_ModuleSupprt[iPlate][0];
    emcElemSupportComp->AddNode(emcElemSupportHole1, iHole, new TGeoTranslation(holePosX, holePosY, 0.0));

    holePosX = EMC_ModuleSupprt[iPlate][1] - EMC_ModuleSupprtX2[iPlate*EMC_NumberPlates+iHole];
    holePosY = EMC_ModuleSupprt[iPlate][4] - EMC_ModuleSupprtY2[iPlate];
    if (iPlate == 0) holePosX = holePosX + EMC_ModuleSupprt[iPlate][0];
////    emcElemSupportComp->AddNode(emcElemSupportHole2, iHole, new TGeoTranslation(holePosX, holePosY, 0.0));

   }
   listECALSupport->Add(emcElemSupportComp);

  }

  TString nameTowerCut;
  Double_t dLengthAdd, dSideAngle;
  for (Int_t iTower = 0; iTower <  EMC_Index_Box; iTower++) {

    TGeoXtru *xElemCut = new TGeoXtru(2);
    nameTowerCut.Form("emc_boxCut%d",iTower+1);
    xElemCut->SetName(nameTowerCut.Data());
    rotation_angle = (2*iTower + 1)*EMC_Box_Angle1;

    dStartPos = dStartNext;
    if ( iTower%EMC_NumberPlates == 0)
     dStartPos = (EMC_TowerUpperLen[iTower/EMC_NumberPlates] -
                                  EMC_Box_C[iTower]*sin(EMC_Box_Angle2))*cos(EMC_Box_Angle2);
    dLenAdd = (dStartPos - EMC_Box_Bottom_X*sin(rotation_angle))/cos(rotation_angle-EMC_Box_Angle1);
    dLenSize = EMC_Box_gap*sin(dAngleRotate)/cos(rotation_angle+EMC_Box_Angle1-dAngleRotate);
    dSideAngle = atan(0.5*(EMC_Box_A[iTower]-EMC_Box_B[iTower])/EMC_Box_Top_X);

    xTruCut[0] = -(0.5*EMC_Box_Bottom_X+dLenAdd*sin(EMC_Box_Angle1));
    yTruCut[0] = dLenAdd*cos(EMC_Box_Angle1) - 0.5*EMC_Box_Length;
    xTruCut[2] =  0.5*EMC_Box_Top_X; yTruCut[2] = 0.5*EMC_Box_Length;
    xTruCut[3] = -0.5*EMC_Box_Top_X; yTruCut[3] = 0.5*EMC_Box_Length;

    if (rotation_angle < dAngleRotate) {
      dLenAdd = 2.*fabs(xTruCut[0])*tan(dAngleRotate-rotation_angle)/
                   (1.+tan(dAngleRotate-rotation_angle)*tan(EMC_Box_Angle1));
      xTruCut[1] = - xTruCut[0] - dLenAdd*tan(EMC_Box_Angle1);
      yTruCut[1] = yTruCut[0] - dLenAdd;
     }
     else {
      dLenAdd =  2*fabs(xTruCut[0])*tan(rotation_angle - dAngleRotate)/
                (1. - tan(rotation_angle - dAngleRotate)*tan(EMC_Box_Angle1));
      xTruCut[1] = - xTruCut[0] + dLenAdd*tan(EMC_Box_Angle1);
      yTruCut[1] = yTruCut[0] + dLenAdd;
     }

    dStartNext = 0.5*(EMC_Box_Length*cos(rotation_angle)+EMC_Box_Bottom_X*sin(rotation_angle))+
        yTruCut[1]*cos(rotation_angle) - xTruCut[1]*sin(rotation_angle);

/////     xTruCut[0] = -0.5*EMC_Box_Top_X; xTruCut[1] = 0.5*EMC_Box_Top_X;

     if ( ( iTower%EMC_NumberPlates == 0 ) || (iTower == 0) ){
                EMC_TowerYshift[nPlate] = EMC_TowerUpperLen[iTower/EMC_NumberPlates];
                EMC_TowerZshift[nPlate] = (EMC_TowerYshift[nPlate] - EMC_Box_Bottom_X*sin(rotation_angle))*
                           tan(rotation_angle-EMC_Box_Angle1);
     }

     if ( (iTower+1)%EMC_NumberPlates == 0) dLenSize = 0;
     dStartNext = dStartNext - dLenSize;

// Check calulated and tabulated transverse sizes

    if ( ( (iTower+1)%EMC_NumberPlates == 0) && (iTower > 0) ) {
/*
     printf("MpdRoot : calulated and tabulated transverse sizes %d %f %f \n",nPlate+1,
                dStartNext/cos(EMC_Box_Angle2)+EMC_Box_C[iTower]*sin(EMC_Box_Angle2),
                EMC_TowerLowerLen[iTower/EMC_NumberPlates]);
*/ 
     nPlate++;
    }
    xElemCut->DefinePolygon(4, xTruCut, yTruCut);
    xElemCut->DefineSection(0, 0.0, 0., 0., 1.);
    xElemCut->DefineSection(1, EMC_TowerCut[iTower/EMC_NumberPlates], 0., 0., 1.);
    listECALTowerCuts->Add(xElemCut);
  }

}

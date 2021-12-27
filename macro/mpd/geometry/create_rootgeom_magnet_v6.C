/** @file create_rootgeom_Magnet_v5_eds.C
    @author  Alexander Krylov avkrylov@jinr.ru
    @date    September 2021
    @brief   Create root geometry for Magnet at NICA
*/

#include <vector>

#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoMedium.h"
#include "TFile.h"

R__ADD_INCLUDE_PATH($VMCWORKDIR)
#include "macro/mpd/mpdloadlibs.C"

//Yoke pars
const Int_t dYokeSections = 28;
const Double_t dYokeSectionStep = 360./dYokeSections; // degree
const Double_t dYokeOuterRad = 658.3/2;
const Double_t dYokeLength = 897;
	//Plate
	const Double_t dYokePlateLength = 847;
	const Double_t dYokePlateInLength = 66.1;
	const Double_t dYokePlateOutLength = 74;
	const Double_t dYokePlateHight = 35;
	//WierdHole
	const Double_t dWierdHoleEdgePosShift = 10;
	const Double_t dWierdHoleLength = 52;
	const Double_t dWierdHoleShift = 11;
	//YokeHole									//		 ___________________________
	const Double_t dYokeHoleCenterShift = 135;	//		|	 10			 	 10		|
	const Double_t dYokeHoleRad = 2;			//		| 	o10o	270		o10o	|
	const Double_t dYokeHoleSpace = 18;			//		|	 18			 	 18		|
	const Double_t dYokeHoleShift = 10;			//		|	o10o	270		o10o	|
	const Double_t dYokeHoleHight = 0.2;
	//SupportRing
	const Double_t dSupportRingLength = 25;
	const Double_t dSupportRingInRad = 229.8;

//Cryostat pars
const Double_t dCryostatLength = 791;
const Double_t dCryostatInRad = 465.6/2;
const Double_t dCryostatOutRad = 544.3/2;
	//Suspension
	const Double_t dSuspensionLength = 247.2;
	const Double_t dSuspensionInRad = 10.25/2;
	const Double_t dSuspensionOutRad = 14.35/2;
	const Double_t dSuspensionLedge = 5.55;
	const Double_t dSuspensionShiftRad = 254;
	const Double_t dSuspensionShiftZ = 369.5;
	const Double_t dSuspensionCapLength = 13.85;
  //const Double_t dCircumscribedCircleRad = 585.6/2;
	//Chimney
	const Double_t dChimneyLength = 75;		//temporary wrong
	const Double_t dChimneyShiftZ = 373.95;
  //const Double_t dChimneyShiftX = 288.1;
	const Double_t dChimneyOutRad = 21.9/2;
	const Double_t dChimneyInRad = 18/2;

//Pole pars
const std::vector<Double_t> dPoleThickness {38, 6, 25}; // In, Mid, Out thickness
const Double_t dPoleHight = dPoleThickness[0] + dPoleThickness[1] + dPoleThickness[2];
const std::vector<Double_t> vPoleRads{191.5/2, 206.6/2, 310.6/2, 315.5/2, 344.1/2,
									  361.872/2, 373.6/2, 384.044/2, 411.4/2, 458.6/2};
	//Trim Coil
	const Double_t dTrimCoilThickness = 10.5;
	const Double_t dTrimCoilRad = 343.1/2;
	const Double_t dTrimCoilCapThickness = 0.5;
	const Double_t dTrimCoilCapInRad = 191.5/2;
	//Assembly Flange
	const Double_t dAssemblyFlangeBarrelLength = 40;
	const std::vector<Double_t> vAssemblyFlangeBarrelHight {2, 15};
	const std::vector<Double_t> vAssemblyFlangeBarrelRad {314.5/2, 321.5/2, 369.5/2, 394.5/2}; // In, Mid, Screw and Out radius
	const Double_t dAssemblyFlangeBarrelShift = 2;
	//Axial Stop
	const Double_t dAxialStopHight = 7;
	const Double_t dAxialStopBoxHigh = 11;
	const Double_t dAxialStopBoxWidth = 15;
	const Double_t dAxialStopBoxLength = 62;
	const Double_t dAxialStopScrewHight = 1.2;
	const Double_t dAxialStopScrewRad = 4;
	const Double_t dAxialStopScrewSpaceX = 29;
	const Double_t dAxialStopScrewSpaceY = 16;

//media
TGeoMedium *pMedAir = 0;	
TGeoMedium *pMedAluminium = 0;
TGeoMedium *pMedTedlar = 0;
TGeoMedium *pMedKevlar = 0;
TGeoMedium *pMedMylar = 0;
TGeoMedium *pMedN2 = 0;
TGeoMedium *pMedRohacellhf71 = 0;
TGeoMedium *pMedPolypropylene= 0;
TGeoMedium *pMedTPCmixture = 0;
TGeoMedium *pMedG10 = 0;
TGeoMedium *pMedFiberGlass = 0;
TGeoMedium *pMedCopper = 0;
TGeoMedium *pMedPlastic = 0;
TGeoMedium *pMedGold = 0;

class FairGeoMedia;
class FairGeoBuilder;

void DefineRequiredMedia(FairGeoMedia* geoMedia, FairGeoBuilder* geoBuild);
void CreateYoke(TGeoVolume* mother_volume);
void CreateYokeHoles(TGeoVolume* mother_volume);
void CreateSupportRing(TGeoVolume* mother_volume);
void CreateCryostat(TGeoVolume* mother_volume);
void CreateChimney(TGeoVolume* mother_volume);
void CreatePole(TGeoVolume* mother_volume);
void CreateAssemblyFlanges(TGeoVolume* mother_volume);
void CreateAxialStops(TGeoVolume* mother_volume);

void create_rootgeom_magnet_v6(Bool_t wrGeoWithMaterials = true)
{
    // ----  set working directory  --------------------------------------------
    TString gPath = gSystem->Getenv("VMCWORKDIR");
    // -------   Geometry file name (output)   ---------------------------------
    const TString geoDetectorName = "magnet";
    const TString geoDetectorVersion = "v6";
    TString geoFileName = gPath + "/geometry/" + geoDetectorName + "_"+ geoDetectorVersion + ".root";
    // -------------------------------------------------------------------------

    // ----  global geometry parameters  ---------------------------------------
    FairGeoLoader*    geoLoad = new FairGeoLoader("TGeo","FairGeoLoader");
    FairGeoInterface* geoFace = geoLoad->getGeoInterface();
    gGeoManager = (TGeoManager*)gROOT->FindObject("FAIRGeom");
    
    // -------   Load media from media file   ----------------------------------
    TString medFile = gPath + "/geometry/media.geo";
    geoFace->setMediaFile(medFile);
    geoFace->readMedia();
    FairGeoMedia*   geoMedia = geoFace->getMedia();
    FairGeoBuilder* geoBuild = geoLoad->getGeoBuilder();

    DefineRequiredMedia(geoMedia, geoBuild);
    // -------------------------------------------------------------------------

    // --------------   Create geometry and global top volume  ------------------------
    gGeoManager->SetName(geoDetectorName + "_geom");
	
	TGeoVolume* top = new TGeoVolumeAssembly("TOP");
    top->SetMedium(pMedAir);
    gGeoManager->SetTopVolume(top);
    //gGeoMan->SetTopVisible(1);
	
	CreateYoke(top);
	CreateCryostat(top);
	CreatePole(top);
	
    top->SetVisContainers(kTRUE);

    // ---------------   Finish   ----------------------------------------------
    gGeoManager->CloseGeometry();
    gGeoManager->CheckOverlaps(0.0001);
    gGeoManager->PrintOverlaps();

    gGeoManager->Test();

    gGeoManager->SetMaxVisNodes(150000);

    TFile* geoFile = new TFile(geoFileName, "RECREATE");
    top->Write();
    geoFile->Close();

    if (wrGeoWithMaterials)
    {
        TString geoFile_wMat = gPath + "/geometry/" + geoDetectorName + "_"+ geoDetectorVersion + "_with_materials.root";
        gGeoManager->Export(geoFile_wMat);
    }

    top->Draw("ogl");
    TGLViewer *v = (TGLViewer*)gPad->GetViewer3D();
    v->SetStyle(TGLRnrCtx::kOutline);
}

void CreateYoke(TGeoVolume* mother_volume) {
	
	TString yoke_name = "yoke";
	TGeoVolumeAssembly * yokeVA = new TGeoVolumeAssembly(yoke_name);
	
	CreateSupportRing(yokeVA);
	
	//Creating yoke beam consisting of plate-wierdhole
	TString yoke_beam_name = "beam";
	TGeoTrd1 *PlateS = new TGeoTrd1(yoke_beam_name + "plate", dYokePlateInLength/2, dYokePlateOutLength/2, dYokePlateLength/2, dYokePlateHight/2);
	
	TGeoXtru *wierdHoleS = new TGeoXtru(2);
	wierdHoleS->SetName(yoke_beam_name + "wierdhole"); 
	Double_t WierdHoleX[13] = {6, 11, 34.5, 39.1, 41, 39.1, 34.5, 11, 6,   -6,  -11, -11, -6};
	Double_t WierdHoleY[13] = {11, 6, 6,    4.6,  0,  -4.6, -6,   -6, -11, -11, -6,  6,   11};
	wierdHoleS->DefinePolygon(13, WierdHoleX, WierdHoleY);
	wierdHoleS->DefineSection(0, -(dYokePlateHight+0.1)); //+0.1 for right cut 
    wierdHoleS->DefineSection(1, dYokePlateHight+0.1);
	
	TGeoCombiTrans *wierdhole_pos = new TGeoCombiTrans(yoke_beam_name + "wierdhole_pos");
	wierdhole_pos->RotateZ(90);
	wierdhole_pos->SetTranslation(0, -dYokePlateLength/2 + dWierdHoleEdgePosShift + dWierdHoleShift, 0);
	wierdhole_pos->RegisterYourself();	
	
	TGeoCombiTrans *wierdhole_pos_refl = new TGeoCombiTrans(*wierdhole_pos);
	wierdhole_pos_refl->SetName(yoke_beam_name + "wierdhole_pos_refl");
	wierdhole_pos_refl->RotateX(180);
	wierdhole_pos_refl->RegisterYourself();
	
	for( Int_t isec = 0; isec < 28; isec++ ) {
		
		TGeoCompositeShape * yokeBeamS = new TGeoCompositeShape(yoke_beam_name, yoke_beam_name + "plate-(" + yoke_beam_name + "wierdhole:" + 
																				yoke_beam_name + "wierdhole_pos+" + yoke_beam_name + "wierdhole:" + 
																				yoke_beam_name + "wierdhole_pos_refl)");
		TGeoCombiTrans *plate_pos = new TGeoCombiTrans(yoke_beam_name + "plate_pos");
		plate_pos->RotateX(90);
		plate_pos->RotateZ(180);
		plate_pos->SetTranslation(0, dYokeOuterRad-dYokePlateHight/2, 0);
		plate_pos->RotateZ(isec*dYokeSectionStep);
		
		TGeoVolume *yokeBeamV = new TGeoVolume(yoke_beam_name, yokeBeamS);
		yokeBeamV->SetMedium(pMedAluminium);
		
		//CreateYokeHoles(yokeBeamV);  //getting overlap (underlap -> +0.2 under volume Box). Not nessessary heavyweight volume
		
		yokeVA->AddNode(yokeBeamV, isec, plate_pos);
	}
	
	mother_volume->AddNode(yokeVA, 0);
}

void CreateYokeHoles(TGeoVolume* mother_volume) {
	
	TString hole_name = "hole";
	
	TGeoPgon *holeS = new TGeoPgon(hole_name, 0., 360., 6, 2);
	holeS->DefineSection(0, dYokeHoleHight/2, 0, dYokeHoleRad);
	holeS->DefineSection(1, -dYokeHoleHight/2, 0, dYokeHoleRad);
	TGeoVolume *holeV = new TGeoVolume(hole_name, holeS);
	holeV->SetMedium(pMedAluminium);
	
	for( Int_t jsec = 0; jsec < 2; jsec++) {
		for( Int_t ksec = 0; ksec < 4; ksec++) {
				
			TGeoCombiTrans *hole_pos = new TGeoCombiTrans(hole_name + "_pos");
			hole_pos->SetTranslation(dYokePlateOutLength/2-dYokeHoleShift-ksec*dYokeHoleSpace, dYokeHoleCenterShift+jsec*dYokeHoleShift, dYokePlateHight/2+dYokeHoleHight/2);
			
			TGeoCombiTrans *hole_pos_refl = new TGeoCombiTrans(*hole_pos);
			hole_pos_refl->RotateZ(180);
			
			mother_volume->AddNode(holeV, 4*jsec+2*ksec, hole_pos);
			mother_volume->AddNode(holeV, 4*jsec+2*ksec+1, hole_pos_refl);
		}
	}
}

void CreateSupportRing(TGeoVolume* mother_volume) {
	
	TString support_ring_name = "support_ring";
	
	TGeoPgon *supportRingFullS = new TGeoPgon(support_ring_name + "full", 0., 360., dYokeSections, 2);
	supportRingFullS->DefineSection(0, dSupportRingLength/2, 0, dYokeOuterRad);
	supportRingFullS->DefineSection(1, -dSupportRingLength/2, 0, dYokeOuterRad);
	
	TGeoTube *supportRingHoleS = new TGeoTube(support_ring_name + "hole", 0., dSupportRingInRad, 900);
	TGeoCompositeShape * supportRingS = new TGeoCompositeShape(support_ring_name, support_ring_name + "full-" + support_ring_name + "hole");
	TGeoVolume *supportRingV = new TGeoVolume(support_ring_name, supportRingS);
	supportRingV->SetMedium(pMedAluminium);
	
	TGeoCombiTrans *support_ring_pos = new TGeoCombiTrans(support_ring_name + "_pos");
	support_ring_pos->SetTranslation(0., 0., dYokeLength/2 - dSupportRingLength/2);
	support_ring_pos->RotateZ(dYokeSectionStep/2);
	support_ring_pos->RegisterYourself();
	
	TGeoCombiTrans *support_ring_pos_refl = new TGeoCombiTrans(*support_ring_pos);
	support_ring_pos_refl->SetName(support_ring_name + "_pos_refl");
	support_ring_pos_refl->RotateY(180);
	support_ring_pos_refl->RegisterYourself();
	
	mother_volume->AddNode(supportRingV, 0, support_ring_pos);
	mother_volume->AddNode(supportRingV, 1, support_ring_pos_refl);
	
}

void CreateCryostat(TGeoVolume* mother_volume) {
	
	//Creating cryostat consisting of tube, suspension and chimney
	TString cryostat_name = "cryostat";
	TGeoVolumeAssembly * cryostatVA = new TGeoVolumeAssembly(cryostat_name);
	
	TGeoTube *cryostatS = new TGeoTube("tube", dCryostatInRad, dCryostatOutRad, dCryostatLength/2);
	TGeoVolume *cryostatV = new TGeoVolume("tube", cryostatS);
	cryostatV->SetMedium(pMedAluminium);
	
	cryostatVA->AddNode(cryostatV, 0);
	
	TGeoPcon *suspensionS = new TGeoPcon("suspension", 0., 360., 6);
	suspensionS->DefineSection(0, -dSuspensionLength/2, 0, dSuspensionInRad);
	suspensionS->DefineSection(1, -dSuspensionLength/2+dSuspensionLedge, 0, dSuspensionInRad);
	suspensionS->DefineSection(2, -dSuspensionLength/2+dSuspensionLedge, 0, dSuspensionOutRad);
	suspensionS->DefineSection(3, dSuspensionLength/2-dSuspensionLedge, 0, dSuspensionOutRad);
	suspensionS->DefineSection(4, dSuspensionLength/2-dSuspensionLedge, 0, dSuspensionInRad);
	suspensionS->DefineSection(5, dSuspensionLength/2, 0, dSuspensionInRad);
	TGeoVolume *suspensionV = new TGeoVolume("suspension", suspensionS);
	suspensionV->SetMedium(pMedAluminium);
	
	for( Int_t isec = 0; isec < 6; isec++) {
		
		
		TGeoCombiTrans *suspension_pos = new TGeoCombiTrans("suspension_pos");
		suspension_pos->RotateX(90);
		suspension_pos->SetTranslation(dSuspensionShiftRad, 0., dSuspensionShiftZ);
		suspension_pos->RotateZ(isec*60);
		
		TGeoCombiTrans *suspension_pos_refl = new TGeoCombiTrans(*suspension_pos);
		suspension_pos_refl->RotateX(180);
		
		cryostatVA->AddNode(suspensionV, 2*isec, suspension_pos);
		cryostatVA->AddNode(suspensionV, 2*isec+1, suspension_pos_refl);
	}
	
	CreateChimney(cryostatVA);
	
	mother_volume->AddNode(cryostatVA, 0);
}

void CreateChimney(TGeoVolume* mother_volume) {
	
	TString chimney_name = "chimney";
	
	TGeoPcon *chimneyS = new TGeoPcon(chimney_name, 0., 360., 4);
	chimneyS->DefineSection(0, 0., dChimneyInRad, dChimneyOutRad);
	chimneyS->DefineSection(1, dChimneyLength, dChimneyInRad, dChimneyOutRad);
	chimneyS->DefineSection(2, dChimneyLength, dChimneyInRad, dChimneyOutRad+5); // 5 - random
	chimneyS->DefineSection(3, dChimneyLength+2, dChimneyInRad, dChimneyOutRad+5); // 2 5 - random
	TGeoVolume *chimneyV = new TGeoVolume(chimney_name, chimneyS);
	chimneyV->SetMedium(pMedAluminium);
	
	TGeoCombiTrans *chimney_pos = new TGeoCombiTrans(chimney_name + "_pos");
	chimney_pos->RotateX(-90);
	chimney_pos->SetTranslation(0., dCryostatOutRad, dChimneyShiftZ);
	
	mother_volume->AddNode(chimneyV, 0, chimney_pos);
}

void CreatePole(TGeoVolume* mother_volume) {

	//Creating magnet pole consisting of support, trim coil and trim coil back and inside caps
	TString pole_name = "pole";
	TGeoVolumeAssembly * poleVA = new TGeoVolumeAssembly(pole_name);
	
	TGeoPcon *poleS = new TGeoPcon("support", 0., 360., 7);
	poleS->DefineSection(0, 0., vPoleRads[5], vPoleRads[6]);
	poleS->DefineSection(1, dTrimCoilThickness, vPoleRads[5], vPoleRads[7]);
	poleS->DefineSection(2, dTrimCoilThickness, vPoleRads[0], vPoleRads[7]);
	poleS->DefineSection(3, dPoleThickness[0], vPoleRads[1], vPoleRads[8]);
	poleS->DefineSection(4, dPoleThickness[0], vPoleRads[2], vPoleRads[8]);
	poleS->DefineSection(5, dPoleThickness[0]+dPoleThickness[1], vPoleRads[3], vPoleRads[9]);
	poleS->DefineSection(6, dPoleHight, vPoleRads[4], vPoleRads[9]);
	TGeoVolume *poleV = new TGeoVolume("support", poleS);
	poleV->SetMedium(pMedAluminium);
	
	poleVA->AddNode(poleV, 0);
	
	TGeoTube *trimColiS = new TGeoTube("trim_coil", vPoleRads[0]+dTrimCoilCapThickness, dTrimCoilRad, dTrimCoilThickness/2);
	TGeoVolume *trimColiV = new TGeoVolume("trim_coil", trimColiS);
	trimColiV->SetMedium(pMedAluminium);
	
	poleVA->AddNode(trimColiV, 0, new TGeoTranslation(0., 0., dTrimCoilThickness/2));
	
	TGeoTube *trimColiBackCapS = new TGeoTube("trim_coil_back_cap", vPoleRads[0], vPoleRads[5], dTrimCoilCapThickness/2);
	TGeoVolume *trimColiBackCapV = new TGeoVolume("trim_coil_back_cap", trimColiBackCapS);
	trimColiBackCapV->SetMedium(pMedAluminium);
	
	poleVA->AddNode(trimColiBackCapV, 0, new TGeoTranslation(0., 0., -dTrimCoilCapThickness/2));
	
	TGeoTube *trimColiInCapS = new TGeoTube("trim_coil_in_cap", vPoleRads[0], vPoleRads[0]+dTrimCoilCapThickness, dTrimCoilThickness/2);
	TGeoVolume *trimColiInCapV = new TGeoVolume("trim_coil_in_cap", trimColiInCapS);
	trimColiInCapV->SetMedium(pMedAluminium);
	
	poleVA->AddNode(trimColiInCapV, 0, new TGeoTranslation(0., 0., dTrimCoilThickness/2));
	
	CreateAssemblyFlanges(poleVA);
	CreateAxialStops(poleVA);
	
	TGeoCombiTrans *pole_pos = new TGeoCombiTrans(pole_name + "_pos");
	pole_pos->SetTranslation(0., 0., dYokeLength/2-dPoleHight);
	TGeoCombiTrans *pole_pos_refl = new TGeoCombiTrans(*pole_pos);
	pole_pos_refl->RotateY(180);
	
	mother_volume->AddNode(poleVA, 0, pole_pos);
	mother_volume->AddNode(poleVA, 1, pole_pos_refl);
}

void CreateAssemblyFlanges(TGeoVolume* mother_volume) {
	
	//Creating magnet pole assembly flange consisting of barrel and screws
	TString assembly_flange_name = "assembly_flange";
	TGeoVolumeAssembly * assemblyFlangeVA = new TGeoVolumeAssembly(assembly_flange_name);
	
	TGeoPcon *barrelS = new TGeoPcon("barrel", 0., 360., 4);
	barrelS->DefineSection(0, 0., vAssemblyFlangeBarrelRad[0], vAssemblyFlangeBarrelRad[3]);
	barrelS->DefineSection(1, vAssemblyFlangeBarrelHight[0],vAssemblyFlangeBarrelRad[0], vAssemblyFlangeBarrelRad[3]);
	barrelS->DefineSection(2, vAssemblyFlangeBarrelHight[0], vAssemblyFlangeBarrelRad[0], vAssemblyFlangeBarrelRad[1]);
	barrelS->DefineSection(3, vAssemblyFlangeBarrelHight[1], vAssemblyFlangeBarrelRad[0], vAssemblyFlangeBarrelRad[1]);
	TGeoVolume *barrelV = new TGeoVolume("barrel", barrelS);
	barrelV->SetMedium(pMedAluminium);
	
	TGeoCombiTrans *barrel_pos = new TGeoCombiTrans("barrel_pos");
	barrel_pos->RotateY(180);
	barrel_pos->SetTranslation(0., 0., dPoleHight+dAssemblyFlangeBarrelShift);
	
	assemblyFlangeVA->AddNode(barrelV, 0, barrel_pos);
	
	TGeoPgon *screwS = new TGeoPgon("screw", 0., 360., 6, 2);
	screwS->DefineSection(0, 0., 0., dAxialStopScrewRad);
	screwS->DefineSection(1, dAxialStopScrewHight, 0., dAxialStopScrewRad);
	TGeoVolume *screwV = new TGeoVolume("screw", screwS);
	screwV->SetMedium(pMedAluminium);
		
	for( Int_t isec = 0; isec < 40; isec++) {
		TGeoCombiTrans *screw_pos = new TGeoCombiTrans("screw_pos");
		screw_pos->SetTranslation(0., vAssemblyFlangeBarrelRad[2], dPoleHight+dAssemblyFlangeBarrelShift);
		screw_pos->RotateZ(9*isec);
		assemblyFlangeVA->AddNode(screwV, isec, screw_pos);
	}
	
	mother_volume->AddNode(assemblyFlangeVA, 0);
}

void CreateAxialStops(TGeoVolume* mother_volume) {
	
	//Creating magnet pole axial stop consisting of plate+box and screws
	TString axial_stop_name = "axial_stop";
	TGeoVolumeAssembly * axialStopVA = new TGeoVolumeAssembly(axial_stop_name);
	
	TGeoXtru *plateS = new TGeoXtru(2);
	plateS->SetName("plate");
	Double_t plateX[8] = {-18, -22, -22, -18, 18, 22, 22, 18};
	Double_t plateY[8] = {-31, -24, 24, 31, 31, 24, -24, -31};
	plateS->DefinePolygon(8, plateX, plateY);
	plateS->DefineSection(0, -dAxialStopHight/2);
    plateS->DefineSection(1, dAxialStopHight/2);
	
	TGeoBBox *boxS = new TGeoBBox("box", dAxialStopBoxWidth/2, dAxialStopBoxLength/2, dAxialStopBoxHigh/2);
	TGeoCombiTrans *box_pos = new TGeoCombiTrans("box_pos");
	box_pos->SetTranslation(0., 0., dAxialStopBoxHigh/2);
	box_pos->RegisterYourself();
	
	TGeoCompositeShape * fullVolumeS = new TGeoCompositeShape("full_volume", "plate+box:box_pos");
	TGeoVolume *fullVolumeV = new TGeoVolume("full_volume", fullVolumeS);
	fullVolumeV->SetMedium(pMedAluminium);
	
	axialStopVA->AddNode(fullVolumeV, 0);
	
	TGeoPgon *screwS = new TGeoPgon("screw", 0., 360., 6, 2);
	screwS->DefineSection(0, 0., 0., dAxialStopScrewRad);
	screwS->DefineSection(1, dAxialStopScrewHight, 0., dAxialStopScrewRad);
	TGeoVolume *screwV = new TGeoVolume("screw", screwS);
	screwV->SetMedium(pMedAluminium);
	
	Int_t count = 0;
	for( Int_t i = -1; i < 2; i+=2) {
		for ( Int_t j = -3; j < 4; j+=2) {
			TGeoCombiTrans *screw_pos = new TGeoCombiTrans("screw_pos");
			screw_pos->SetTranslation(dAxialStopScrewSpaceX/2*i, dAxialStopScrewSpaceY/2*j, dAxialStopHight/2);
			axialStopVA->AddNode(screwV, count, screw_pos);
			
			TGeoCombiTrans *axial_stop_pos = new TGeoCombiTrans(axial_stop_name + "_pos");
			axial_stop_pos->SetTranslation(0., vPoleRads[9], dPoleHight+dAxialStopHight/2);
			axial_stop_pos->RotateZ(45/2+ 45*count);
			mother_volume->AddNode(axialStopVA, count, axial_stop_pos);
			count++;
		}	
	}
}

void DefineRequiredMedia(FairGeoMedia* geoMedia, FairGeoBuilder* geoBuild) {
    //air medium
    FairGeoMedium* mAir = geoMedia->getMedium("air");
    if ( ! mAir ) Fatal("Main", "FairMedium air not found");
    geoBuild->createMedium(mAir);
    pMedAir = gGeoManager->GetMedium("air");
    if ( ! pMedAir ) Fatal("Main", "Medium air not found");

    //aluminium medium
    FairGeoMedium* mAluminium = geoMedia->getMedium("aluminium");
    if ( ! mAluminium  ) Fatal("Main", "FairMedium aluminium not found");
    geoBuild->createMedium(mAluminium);
    pMedAluminium  = gGeoManager->GetMedium("aluminium");
    if ( ! pMedAluminium  ) Fatal("Main", "Medium aluminium not found");

    //tedlar medium
    FairGeoMedium* mTedlar = geoMedia->getMedium("tedlar");
    if ( ! mTedlar  ) Fatal("Main", "FairMedium tedlar not found");
    geoBuild->createMedium(mTedlar);
    pMedTedlar = gGeoManager->GetMedium("tedlar");
    if ( ! pMedTedlar  ) Fatal("Main", "Medium tedlar not found");

    //kevlar medium
    FairGeoMedium* mKevlar = geoMedia->getMedium("kevlar");
    if ( ! mKevlar  ) Fatal("Main", "FairMedium kevlar not found");
    geoBuild->createMedium(mKevlar);
    pMedKevlar = gGeoManager->GetMedium("kevlar");
    if ( ! pMedKevlar  ) Fatal("Main", "Medium kevlar not found");

    //mylar medium
    FairGeoMedium* mMylar = geoMedia->getMedium("mylar");
    if ( ! mMylar  ) Fatal("Main", "FairMedium mylar not found");
    geoBuild->createMedium(mMylar);
    pMedMylar = gGeoManager->GetMedium("mylar");
    if ( ! pMedMylar  ) Fatal("Main", "Medium mylar not found");

    //N2 medium
    FairGeoMedium* mN2 = geoMedia->getMedium("N2");
    if ( ! mN2  ) Fatal("Main", "FairMedium N2 not found");
    geoBuild->createMedium(mN2);
    pMedN2 = gGeoManager->GetMedium("N2");
    if ( ! pMedN2  ) Fatal("Main", "Medium N2 not found");

    //rohacellhf71 medium
    FairGeoMedium* mRohacellhf71 = geoMedia->getMedium("rohacellhf71");
    if ( ! mRohacellhf71 ) Fatal("Main", "FairMedium rohacellhf71 not found");
    geoBuild->createMedium(mRohacellhf71);
    pMedRohacellhf71 = gGeoManager->GetMedium("rohacellhf71");
    if ( ! pMedRohacellhf71  ) Fatal("Main", "Medium rohacellhf71 not found");

    //polypropylene medium
    FairGeoMedium* mPolypropylene = geoMedia->getMedium("polypropylene");
    if ( ! mPolypropylene ) Fatal("Main", "FairMedium polypropylene not found");
    geoBuild->createMedium(mPolypropylene);
    pMedPolypropylene = gGeoManager->GetMedium("polypropylene");
    if ( ! pMedPolypropylene ) Fatal("Main", "Medium polypropylene not found");

    //TPCmixture medium
    FairGeoMedium* mTPCmixture = geoMedia->getMedium("TPCmixture");
    if ( ! mTPCmixture ) Fatal("Main", "FairMedium TPCmixture not found");
    geoBuild->createMedium(mTPCmixture);
    pMedTPCmixture = gGeoManager->GetMedium("TPCmixture");
    if ( ! pMedTPCmixture ) Fatal("Main", "Medium TPCmixture not found");

    //G10 medium
    FairGeoMedium* mG10 = geoMedia->getMedium("G10");
    if ( ! mG10 ) Fatal("Main", "FairMedium G10 not found");
    geoBuild->createMedium(mG10);
    pMedG10 = gGeoManager->GetMedium("G10");
    if ( ! pMedG10 ) Fatal("Main", "Medium G10 not found");

    //fiberglass medium
    FairGeoMedium* mFiberGlass = geoMedia->getMedium("fiberglass");
    if ( ! mFiberGlass ) Fatal("Main", "FairMedium fiberglass not found");
    geoBuild->createMedium(mFiberGlass);
    pMedFiberGlass = gGeoManager->GetMedium("fiberglass");
    if ( ! pMedFiberGlass ) Fatal("Main", "Medium fiberglass not found");

    //copper medium
    FairGeoMedium* mCopper = geoMedia->getMedium("copper");
    if ( ! mCopper ) Fatal("Main", "FairMedium copper not found");
    geoBuild->createMedium(mCopper);
    pMedCopper = gGeoManager->GetMedium("copper");
    if ( ! pMedCopper ) Fatal("Main", "Medium copper not found");

    //gold medium
    FairGeoMedium* mGold = geoMedia->getMedium("gold");
    if ( ! mGold) Fatal("Main", "FairMedium gold not found");
    geoBuild->createMedium(mGold);
    pMedGold = gGeoManager->GetMedium("gold");
    if ( ! pMedGold ) Fatal("Main", "Medium gold not found");

    //plastic medium
    FairGeoMedium* mPlastic = geoMedia->getMedium("plastic");
    if ( ! mPlastic) Fatal("Main", "FairMedium plastic not found");
    geoBuild->createMedium(mPlastic);
    pMedPlastic = gGeoManager->GetMedium("plastic");
    if ( ! pMedPlastic ) Fatal("Main", "Medium plastic not found");
}
/** @file create_rootgeom_TPC_v9.C
    @author  Alexander Bychkov abychkov@jinr.ru
    @date    September 2020
    @brief   Create root geometry for TPC
*/

#include <vector>
//---------------------------
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoMedium.h"
#include "TFile.h"

#include "commonFunctions.C"
#include "geometry_stage1.C"

//Detector's position
const Double_t TPC_Xpos = 0.0;
const Double_t TPC_Ypos = 0.0;
const Double_t TPC_Zpos = 0.0;

//TPC pars
Int_t Nsections = 12;
Double_t Section_step = 360./Nsections; // degree
Double_t Section_rad_step = TMath::DegToRad()*Section_step; //radian

//Detector's construct parameters

//Global dimensions for TPC in cm
const Double_t TpcInnerRadius = 27.0; // tpc inner radius
const Double_t TpcOuterRadius = 140.7; // tpc outer radius
const Double_t Container_tpc_z = 400.0; // common tpc length
const Double_t Chamber_tpc_z = 340.0;
const Double_t Chamber_tpc_z_ext_parts = 0.5;

//Barrel walls in cm
const Double_t Barrel_Al_thick = 0.01;
const Double_t Barrel_tedlar_thick = 0.01;
const std::vector<Double_t> Thick_barell_CO2_layer{6.5, 6.7}; // [inner, outer] barrel
const std::vector<Double_t> Thick_barell_inner_layer{0.3, 0.4}; // [inner, outer] barrel
const std::vector<Double_t> Thick_barell_outer_layer{0.3, 0.6}; // [inner, outer] barrel

//inner wall extension parts in cm
const Double_t Flange_in_R = 27.32;
const Double_t Flange_thickness = 5.0;
const Double_t Flange_width = 2.85;
const Double_t ExtPart_length = 20.0;
const std::vector<Double_t> Al_Ring_R{30.17, 30.5, 33.5, 33.79};
const std::vector<Double_t> Al_Ring_width{4.4, 4.7};
const Double_t AlFoil_ring_thikness = 2.5;
const Double_t AlFoil_ring_gap = 10.0;

//outer wall extension parts in cm
const Double_t Flange_out_thickness = 5.0;
const Double_t Flange_out_width = 6.47;

//Membrane
const Double_t Membrane_thickness = 0.02;
const Double_t Membrane_AlFoil_thickness = 0.005;
const Double_t Membrane_inner_holder_R_edge = 40.0;
const Double_t Membrane_inner_holder_R = Membrane_inner_holder_R_edge * TMath::Cos(Section_rad_step / 2.);
const Double_t Membrane_outer_holder_R_edge = 130.0;
const Double_t Membrane_outer_holder_R = Membrane_outer_holder_R_edge * TMath::Cos(Section_rad_step / 2.);
const Double_t Membrane_holder_thikness = 1.5;

//field cage
const Double_t Fieldcage_shift_deg = 8.;
const Double_t Fieldcage_shift_rad = Fieldcage_shift_deg * TMath::DegToRad();
const std::vector<Double_t> Fieldcage_Pin_R{2., 3.}; // [inner, outer]
const std::vector<Double_t> Fieldcage_Pin_wallthik{0.3, 0.9}; // [inner, outer]
const Double_t Fieldcage_ribbon_thikness = 0.01;
const Double_t Fieldcade_ribbon_wide = 1.3;
const Double_t Fieldcage_ribbon_gap = 0.2;

//ends
const Double_t Flange_over_inner_tube_thickness = 5.0;
const Double_t Flange_over_inner_tube_width = 0.845;
const Double_t End_Frame_thickness = 5.0;
const Double_t End_Frame_in_R = 35.4;
const Double_t End_Frame_in_R_seg = 40.15;
const Double_t End_Frame_out_R = 132.4;
const Double_t End_Frame_out_R_seg = 120.45;
const Double_t End_Frame_Spoke_width = 3.8;
const Double_t End_Frame_Spoke_thick = 4.0;

//pad plane simulation in cm
const Double_t ROC_Frame_Gap_to_End_Frame = 0.15;
const Double_t ROC_Frame_Gap_to_End_Frame_Side = 0.25;
const Double_t Sensitive_in_R_seg = End_Frame_in_R_seg + ROC_Frame_Gap_to_End_Frame;
const Double_t Sensitive_out_R_seg = End_Frame_out_R_seg - ROC_Frame_Gap_to_End_Frame;
const Double_t Pad_Plane_Pp = 0.3109;
const Double_t Pad_Thick = 0.0109;
const Double_t Pad_Width = 0.48;
const std::vector<Double_t> Pad_Height = {1.18, 1.78};
const Double_t Pad_Space_w = 0.02;
const Double_t Pad_Space_h = 0.02;
const Double_t Pad_Space_edge = 0.4;
const std::vector<Int_t> NPads_in_Row_s = {20, 21, 21, 22, 23, 23, 24, 24, 25, 26, 26, 27, 28, 28, 29, 30, 30, 31, 32, 32, 33, 33, 34, 35, 35, 36, 37}; //half
const std::vector<Int_t> NPads_in_Row_l = {38, 39, 40, 41, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62}; //half
const Double_t Pad_Plane_G10 = 0.3;
const Double_t ROC_Frame = 6.0;
const Double_t ROC_Frame_wo_plane = 5.5;
const Double_t ROC_Frame_plane = ROC_Frame - ROC_Frame_wo_plane;
const Double_t ROC_Frame_width = 1.9;
const Double_t ROC_Frame_center_width = 1.0;
const Double_t ROC_Frame_lower_chamber_height = 41.7;
const Double_t ROC_Frame_ext_width = 2.0;
const Double_t ROC_Frame_ext_thick = 1.5 - ROC_Frame_plane;
const Double_t ROC_Frame_zspace_iocards = 5.0;

//electronics moodule
const Double_t Connector_hole_length = 8.0;
const Double_t Connector_hole_width = 1.2;
const Double_t Connector_length = 3.52;
const Double_t Connector_width = 0.36;

const Double_t ElModule_CoolPlate_Thick = 0.47;
const Double_t ElModule_CoolPlateAdditional_Thick = 0.45;
const Double_t ElModule_CoolPlateAdditional_Gap = 0.58;
const Double_t ElModule_CoolPlate_ZOffsetFrameOutEdge = 0.7;
const Double_t ElModule_CoolPlate_EdgeOffsetFrame = 0.5;
const Double_t ElModule_CoolPlate_Hole_H = 7.26;
const Double_t ElModule_CoolPlate_Hole_W = 0.8;
//const Double_t ElModule_CoolPlate_InWidth = 40.7;
//const Double_t ElModule_CoolPlate_OutWidth = 32.5;
const Double_t ElModule_Cool_Chip0_WH = 1.7;
const Double_t ElModule_Cool_Chip0_Th = 0.47;
const Double_t ELModule_Chip0_Th = 0.16;
const Double_t ELModule_Chip1_Th = 0.14; //real 0.16, cool thick 0.33-0.35
const Double_t ElModule_Cool_Chip1_H = 2.6;
const Double_t ElModule_Cool_Chip1_W = 1.7;
const Double_t ElModule_Cool_Chip1_Th = 0.35;
const Double_t ELModule_Width = 2.8;
const Double_t ELModule_Height = 9.1;
const Double_t ELModule_PCBBox_Th = 0.7265;
const Double_t ELModule_PCB_Th = 0.15;
const Double_t ELModule_PCBCu_Th = 0.05;
const Double_t ElModule_Cool_Chip0_HOffsetUp = ELModule_Height / 2. - ElModule_Cool_Chip0_WH / 2. - 1.6;
const Double_t ElModule_Cool_Chip0_HOffsetDown = -(ELModule_Height / 2. - ElModule_Cool_Chip0_WH / 2. - 1.7);
const Double_t ElModule_Cool_Chip0_WOffset = -ELModule_Width / 2. + 1. + ElModule_Cool_Chip0_WH / 2.;
const Double_t ElModule_Cool_Chip1_HOffset = ELModule_Height / 2. - ElModule_Cool_Chip1_H / 2. - 2.835;
const Double_t ElModule_Cool_Chip1_WOffset = -ELModule_Width / 2. + 0.1498 + ElModule_Cool_Chip0_WH / 2.;
const Double_t ElModule_Cool_Th = 0.45;
const Double_t ElModule_Cool_Hole_H = 6.348;
const Double_t ElModule_Cool_Hole_W = 0.62;
const Double_t ElModule_Cool_Hole_HOffset = ELModule_Height / 2. - ElModule_Cool_Hole_H / 2. - 0.6;
const Double_t ElModule_Cool_Hole_WOffset = ELModule_Width / 2. - ElModule_Cool_Hole_W / 2. - 0.19; // 0.39 if ElModule_Cool_W = 3.1
const Double_t ELModule_PCB0_Connector_W = 0.48;
const Double_t ELModule_PCB0_Connector_H = 7.62;
const Double_t ELModule_PCB0_Connector_Th = 0.578;
const Double_t ElModule_PCB0_Connector_WOffset = -ELModule_Width / 2. + ELModule_PCB0_Connector_W / 2.;

//ROC pressing
const Double_t ROCPress_Backplate_thick = 0.02;
const Double_t ROCPress_Backplate_inR = 36.75;
const Double_t ROCPress_Backplate_outR = 123.85;
const Double_t ROCPress_Backplate_LowWidth = 5.45;
const Double_t ROCPress_Backplate_SideWidth = 3.95;
const Double_t ROCPress_Backplate_InSpace = 76.2;
const Double_t ROCPress_Backplate_HighWidth = ROCPress_Backplate_outR - (ROCPress_Backplate_inR + ROCPress_Backplate_LowWidth + ROCPress_Backplate_InSpace);
const Double_t ROCPress_Backplate_GapToNextSect = 0.2;
const std::vector<Double_t> ROCPress_Backplate_HighExtL = {18.761, 3.839, 5.4104, 5.2446};
const Double_t ROCPress_Backplate_HighExtX = 1.5;
const Double_t ROCPress_OutFrame_thick = 2.0;
const Double_t ROCPress_OutFrame_width = 1.8;
const Double_t ROCPress_OutFrame_LowCenterAngle = 160;
const Double_t ROCPress_OutFrame_LowCenterAngle_rad = ROCPress_OutFrame_LowCenterAngle * TMath::DegToRad();
const Double_t ROCPress_InFrame_thick = 1.7;
const Double_t ROCPress_InFrame_width = 1.9;

//Ribs
const Double_t RibProfileThickness = 0.8;
const Double_t RibProfileWidth = 7.6;
const Double_t RibProfileHeight = 4.0;
const Double_t RibsSectorRMin = 34.55;
const Double_t RibsSectorRMax = 132;
const Double_t RibsZfromTPCends = 21.4 - Chamber_tpc_z_ext_parts - ROC_Frame_zspace_iocards;

//cooling
const Double_t Thermoscreen_Ends_inR = 42.253;
const Double_t Thermoscreen_Ends_outR = Thermoscreen_Ends_inR + 84.7;
const Double_t Thermoscreen_Ends_GapToNecxtSect = 4.961;
const Double_t Thermoscreen_Thick = 2 * 0.1;
const Double_t Thermoscreen_Rad_Length = 438.2;
const Double_t Thermoscreen_Rad_R = 143.2;
const Double_t Thermoscreen_Rad_Seg_Angle = 30;
const Double_t Thermoscreen_Rad_Seg_Angle_Rad = Thermoscreen_Rad_Seg_Angle * TMath::DegToRad();
const Double_t Thermoscreen_Rad_Plate_Lenght = 101.;
const Double_t Thermoscreen_Rad_Plate_ZGap = 0.4;
const Double_t Thermoscreen_Rad_Plate_Width_Seg = 24.6;
const Double_t Thermoscreen_Rad_Hold_Width = 5.;
const Double_t Thermoscreen_Rad_Hold_Height = 5.;
const Double_t Thermoscreen_Rad_Hold_HeightIn = 3.;
const Double_t Thermoscreen_Rad_Hold_Thick = 1.;
const Double_t Thermoscreen_Rad_Rail_CentralPlate_Width = 20.;
const Double_t Thermoscreen_Rad_Rail_Bend_R = 0.4;
const Double_t Thermoscreen_Rad_Rail_BendCenters_Dist = 3.;
const Double_t Thermoscreen_Rad_Rail_CentralPlate_SideAngle = 30.;
const Double_t Thermoscreen_Rad_Rail_CentralPlate_SideAngle_Rad = 30. * TMath::DegToRad();
const Double_t Thermoscreen_Rad_Rail_CentralPlate_DistToZaxis = 140.8;
const Double_t Thermoscreen_Rad_Rail_SidePlate_Width = 23.8;
const Double_t Thermoscreen_Rad_Rail_Lenght_Cut = 101. - 68.;

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

enum ElectronicsGeometryDetailLevel
{
    Low,
    Med,
    High
};
ElectronicsGeometryDetailLevel ElectronicsGeometryDetails = High;

class FairGeoMedia;
class FairGeoBuilder;

void DefineRequiredMedia(FairGeoMedia* geoMedia, FairGeoBuilder* geoBuild);
void CreateInnerWall(TGeoVolume* mother_volume);
void CreateOuterWall(TGeoVolume* mother_volume);
void CreateTubeAlFoilRings(TGeoVolume* mother_volume, Double_t rmin, Double_t rmax, Bool_t is_outer_wall);
void CreateMembrane(TGeoVolume* mother_volume);
void CreateInnerStructure(TGeoVolume* mother_volume);
void CreateMembraneFrame(TGeoVolume* mother_volume, Double_t width, Double_t rmin, Double_t rmax);
void CreateFieldCage(TGeoVolume *mother_volume, Double_t vol_width);
void CreateEndsFrame(TGeoVolume *mother_volume, Double_t width, Double_t rmin, Double_t rmax);
void CreateSensitiveStructure(TGeoVolume *mother_volume, Double_t vol_width);
void CreateSensitiveVolumes(TGeoVolume *mother_volume, Double_t vol_width);
void CreatePadWithFrameStructure(TGeoVolume *mother_volume, Double_t vol_width);
void CreatePadPlaneWithPads(TGeoVolume *mother_volume, Double_t R_seg_in, Double_t R_seg_out);
void PlacePadsSub(TGeoVolume *mother_volume, TGeoVolume * pad_vol, Double_t R_in, const std::vector<Int_t> & nPadsInRow, Double_t padHeight);
void FillFramesPoints(unique_ptr<Double_t[]> &x, unique_ptr<Double_t[]> &y);
void PlaceElectronicsComponents(TGeoVolume * mother_vol, TGeoVolume * vol_to_place, Double_t volWidth, Double_t z_offset = 0., Bool_t ElModule = false);
void PlaceElectronicsComponentsSub(TGeoVolume * mother_vol, TGeoVolume * vol_to_place, Int_t nxlay, Int_t nylay[], 
    Double_t x0, Double_t y0, Double_t low_x, Double_t half_low_y, Double_t volWidth, Double_t z_offset = 0., Int_t index_offset = 0, Bool_t ElModule = false);
void CreateElectronicsModules(TGeoVolume *mother_volume);
void CreateEndCaps(TGeoVolume *mother_volume);
void MakeRibProfilePoints(unique_ptr<Double_t[]> &x, unique_ptr<Double_t[]> &y);
void MakeRibProfileTriConnectorsPoints(unique_ptr<Double_t[]> &x, unique_ptr<Double_t[]> &y);
void CreateRibs(TGeoVolume *mother_volume, Double_t endcap_zwidth);
void MakeROCFramesHoles(TGeoVolume * bpG10Vol, TGeoVolume * bpAlVol);
void CreateROCPressingBackPlane(TGeoVolume *mother_volume, Double_t width, Double_t rmin, Double_t rmax);
void CreateROCPressingInFrame(TGeoVolume *mother_volume, Double_t width, Double_t rmin, Double_t rmax);
void CreateROCPressingOutFrame(TGeoVolume *mother_volume, Double_t width, Double_t rmin, Double_t rmax);
void CreateThermoscreenEnds(TGeoVolume *mother_volume, Double_t endcap_zwidth);
void CreateThermoscreenRadial(TGeoVolume *mother_volume);
inline Double_t ThermoscreenRadialPlateROffset();
inline Double_t ThermoscreenRadialPlateRailBendHeight(Double_t R1, Double_t R2);
TGeoVolume * CreateMotherVolTPCandThermoscreen();

void create_rootgeom_TPC_v10(Bool_t wrGeoWithMaterials = true)
{
    // Load necessary libraries
    //gROOT->LoadMacro("$VMCWORKDIR/macro/mpd/mpdloadlibs.C");
    //mpdloadlibs(); // load libraries

    // ----  set working directory  --------------------------------------------
    TString gPath = gSystem->Getenv("VMCWORKDIR");

    // -------   Geometry file name (output)   ---------------------------------
    const TString geoDetectorName = "tpc";
    const TString geoDetectorVersion = "v10";
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
    // --------------------------------------------------------------------------
    TGeoVolume * tpc_mother_vol = CreateMotherVolTPCandThermoscreen();
    tpc_mother_vol->SetMedium(pMedAir);
    tpc_mother_vol->SetTransparency(95);

    TGeoTranslation *tpc_mother_vol_pos = new TGeoTranslation("tpc_mother_vol_pos", TPC_Xpos, TPC_Ypos, TPC_Zpos);

    top->AddNode(tpc_mother_vol, 0, tpc_mother_vol_pos);
    ///////////////////////////////////////////////////////////////////////////////////

    // Define TPC Geometry
    TString tpc_chamber_name = "tpcChamber1";
    TGeoTube *tpcChamberS = new TGeoTube(tpc_chamber_name, TpcInnerRadius, TpcOuterRadius, Container_tpc_z/2);
    TGeoVolume *tpcChamberV = new TGeoVolume(tpc_chamber_name, tpcChamberS);
    tpcChamberV->SetMedium(pMedAir);
    //tpcChamberV->SetVisibility(kTRUE);
    //tpcChamberV->SetTransparency(95);

    //TPC position
    TGeoTranslation *tpc_position= new TGeoTranslation("tpc_position", 0., 0., 0.);
    //tpc_position->RegisterYourself();

    CreateInnerWall(tpcChamberV);
    CreateOuterWall(tpcChamberV);
    
    CreateMembrane(tpcChamberV);
    
    CreateInnerStructure(tpcChamberV);
    
    CreateEndCaps(tpcChamberV);

    //Adding the tpc mother valume to the global TOP Volume
    tpc_mother_vol->AddNode(tpcChamberV, 0, tpc_position);
    
    CreateThermoscreenRadial(tpc_mother_vol);

    top->SetVisContainers(kTRUE);


    // ---------------   Finish   ----------------------------------------------
    gGeoManager->CloseGeometry();
    gGeoManager->CheckOverlaps(0.001);
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
        //TString geoFile_gdml = gPath + "/geometry/" + geoDetectorName + "_"+ geoDetectorVersion + ".gdml";
        //gGeoManager->Export(geoFile_gdml);
    }

    /*top->Draw("ogl");
    TGLViewer *v = (TGLViewer*)gPad->GetViewer3D();
    v->SetStyle(TGLRnrCtx::kOutline);*/
}

void CreateInnerWall(TGeoVolume* mother_volume)
{
    TString wall_name = "tpc01InWall";
    Double_t wall_thickness = Thick_barell_outer_layer[0] +
                              Thick_barell_CO2_layer[0] +
                              Thick_barell_inner_layer[0] +
                              Barrel_Al_thick;
    Double_t wall_rmin = TpcInnerRadius;
    Double_t wall_rmax = TpcInnerRadius + wall_thickness;

    TGeoTube *tpcWallS = new TGeoTube(wall_name, wall_rmin, wall_rmax, Chamber_tpc_z/2);
    TGeoVolume *tpcWallV = new TGeoVolume(wall_name, tpcWallS);
    tpcWallV->SetMedium(pMedN2);
    tpcWallV->SetLineColor(TColor::GetColor("#b6e1fc"));
    tpcWallV->SetTransparency(90);
    
    mother_volume->AddNode(tpcWallV, 0);

    //Create outer tube in Wall
    Double_t tubeout_rmax = wall_rmin + Thick_barell_outer_layer[0];
    Double_t tubeout_rmin = wall_rmin;
    TString tubeout_name = "tpc01OuterTube";
    
    TGeoTube *tubeOutS = new TGeoTube(tubeout_name, tubeout_rmin, tubeout_rmax, Chamber_tpc_z/2);
    TGeoVolume *tubeOutV = new TGeoVolume(tubeout_name, tubeOutS);
    tubeOutV->SetMedium(pMedKevlar);
    tubeOutV->SetLineColor(TColor::GetColor("#0000ff"));
    //tubeOutV->SetTransparency(0);
    
    tpcWallV->AddNode(tubeOutV, 0);
    
    TString tubeout_tdl_name = "tpc01OuterTubeTdl";
    
    TGeoTube *tubeOutTdlS = new TGeoTube(tubeout_tdl_name, tubeout_rmax - Barrel_tedlar_thick, tubeout_rmax, Chamber_tpc_z/2);
    TGeoVolume *tubeOutTdlV = new TGeoVolume(tubeout_tdl_name, tubeOutTdlS);
    tubeOutTdlV->SetMedium(pMedTedlar);
    tubeOutTdlV->SetLineColor(TColor::GetColor("#aa00ff"));
    
    tubeOutV->AddNode(tubeOutTdlV, 0);
    
    //Create inner tube in wall
    Double_t tubein_rmin = wall_rmax - Thick_barell_inner_layer[0] - 2 * Barrel_Al_thick;
    Double_t tubein_rmax = wall_rmax;
    TString tubein_name = "tpc01InnerTube";
    
    TGeoTube *tubeInS = new TGeoTube(tubein_name, tubein_rmin, tubein_rmax, Chamber_tpc_z/2);
    TGeoVolume *tubeInV = new TGeoVolume(tubein_name, tubeInS);
    tubeInV->SetMedium(pMedKevlar);
    tubeInV->SetLineColor(TColor::GetColor("#00ff00"));
    //tubeInV->SetTransparency(0);
    
    tpcWallV->AddNode(tubeInV, 0);
    
    CreateTubeAlFoilRings(tubeInV, tubein_rmin, tubein_rmax, false);
    
    //create flanges of inner wall
    TString tube_flange_name = "tpc01TubeFlange";
    Double_t tube_flange_rmin = Flange_in_R;
    Double_t tube_flange_rmax = tube_flange_rmin + Flange_width;

    TGeoTube * TubeFlangeDoubleS = new TGeoTube (tube_flange_name + "_double", tube_flange_rmin, tube_flange_rmax, 2 * Flange_thickness/2);
    TGeoVolume *TubeFlangeDoubleV = new TGeoVolume(tube_flange_name + "_double", TubeFlangeDoubleS);
    TubeFlangeDoubleV->SetMedium(pMedFiberGlass);
    TubeFlangeDoubleV->SetLineColor(TColor::GetColor("#c7fcec"));

    TGeoTube * TubeFlangeS = new TGeoTube (tube_flange_name, tube_flange_rmin, tube_flange_rmax, Flange_thickness/2);
    TGeoVolume *TubeFlangeV = new TGeoVolume(tube_flange_name, TubeFlangeS);
    TubeFlangeV->SetMedium(pMedFiberGlass);
    TubeFlangeV->SetLineColor(TColor::GetColor("#c7fcec"));

    TGeoCombiTrans *tube_flange_pos1 = new TGeoCombiTrans(0, 0, Chamber_tpc_z / 2 - Flange_thickness / 2, NULL);
    TGeoCombiTrans *tube_flange_pos2 = new TGeoCombiTrans(0, 0, Chamber_tpc_z / 2 - ExtPart_length, NULL);
    TGeoCombiTrans *tube_flange_pos1refl = new TGeoCombiTrans(*tube_flange_pos1);
    tube_flange_pos1refl->ReflectZ(true);
    TGeoCombiTrans *tube_flange_pos2refl = new TGeoCombiTrans(*tube_flange_pos2);
    tube_flange_pos2refl->ReflectZ(true);
    
    tpcWallV->AddNode(TubeFlangeV, 0, tube_flange_pos1);
    tpcWallV->AddNode(TubeFlangeDoubleV, 0, tube_flange_pos2);
    tpcWallV->AddNode(TubeFlangeV, 1, tube_flange_pos1refl);
    tpcWallV->AddNode(TubeFlangeDoubleV, 1, tube_flange_pos2refl);
    
    //create Al Rings of inner wall
    TString al_ring_name = "tpc01InWallAlRing";
    
    TGeoPcon * AlRingS = new TGeoPcon(al_ring_name, 0.0, 360.0, 4);
    AlRingS->DefineSection(0, 0.0, Al_Ring_R[1], Al_Ring_R[2]);
    AlRingS->DefineSection(1, -Al_Ring_width[0], Al_Ring_R[1], Al_Ring_R[2]);
    AlRingS->DefineSection(2, -Al_Ring_width[0], Al_Ring_R[0], Al_Ring_R[3]);
    AlRingS->DefineSection(3, -Al_Ring_width[1], Al_Ring_R[0], Al_Ring_R[3]);
    TGeoVolume *AlRingV = new TGeoVolume(al_ring_name, AlRingS);
    AlRingV->SetMedium(pMedAluminium);
    AlRingV->SetLineColor(TColor::GetColor("#c0c0c0"));
    
    TGeoCombiTrans *al_ring_pos = new TGeoCombiTrans(0, 0, Chamber_tpc_z / 2, NULL);
    TGeoCombiTrans *al_ring_pos_refl = new TGeoCombiTrans(*al_ring_pos);
    al_ring_pos_refl->ReflectZ(true);
    
    tpcWallV->AddNode(AlRingV, 0, al_ring_pos);
    tpcWallV->AddNode(AlRingV, 1, al_ring_pos_refl);
}

void CreateOuterWall(TGeoVolume* mother_volume)
{
    TString wall_name = "tpc01OutWall";
    Double_t wall_thickness = Thick_barell_outer_layer[1] +
                              Thick_barell_CO2_layer[1] +
                              Thick_barell_inner_layer[1] +
                              Barrel_Al_thick;
    Double_t wall_rmin = TpcOuterRadius - wall_thickness;
    Double_t wall_rmax = TpcOuterRadius;

    TGeoTube *tpcWallS = new TGeoTube(wall_name, wall_rmin, wall_rmax, Chamber_tpc_z/2);
    TGeoVolume *tpcWallV = new TGeoVolume(wall_name, tpcWallS);
    tpcWallV->SetMedium(pMedN2);
    tpcWallV->SetLineColor(TColor::GetColor("#b6e1fc"));
    
    mother_volume->AddNode(tpcWallV, 0);

    //Create outer tube in Wall
    Double_t tubeout_rmax = wall_rmax;
    Double_t tubeout_rmin = wall_rmax - Thick_barell_outer_layer[1];
    TString tubeout_name = "tpc01OuterTube";
    
    TGeoTube *tubeOutS = new TGeoTube(tubeout_name, tubeout_rmin, tubeout_rmax, Chamber_tpc_z/2);
    TGeoVolume *tubeOutV = new TGeoVolume(tubeout_name, tubeOutS);
    tubeOutV->SetMedium(pMedKevlar);
    tubeOutV->SetLineColor(TColor::GetColor("#0000ff"));
    tubeOutV->SetTransparency(0);
    
    tpcWallV->AddNode(tubeOutV, 0);
    
    TString tubeout_tdl_name = "tpc01OuterTubeTdl";
    
    TGeoTube *tubeOutTdlS = new TGeoTube(tubeout_tdl_name, tubeout_rmin, tubeout_rmin + Barrel_tedlar_thick, Chamber_tpc_z/2);
    TGeoVolume *tubeOutTdlV = new TGeoVolume(tubeout_tdl_name, tubeOutTdlS);
    tubeOutTdlV->SetMedium(pMedTedlar);
    tubeOutTdlV->SetLineColor(TColor::GetColor("#aa00ff"));
    
    tubeOutV->AddNode(tubeOutTdlV, 0);
    
    //Create inner tube in wall
    Double_t tubein_rmin = wall_rmin;
    Double_t tubein_rmax = wall_rmin + Thick_barell_inner_layer[1] + 2 * Barrel_Al_thick;
    TString tubein_name = "tpc01InnerTube";
    
    TGeoTube *tubeInS = new TGeoTube(tubein_name, tubein_rmin, tubein_rmax, Chamber_tpc_z/2);
    TGeoVolume *tubeInV = new TGeoVolume(tubein_name, tubeInS);
    tubeInV->SetMedium(pMedKevlar);
    tubeInV->SetLineColor(TColor::GetColor("#00ff00"));
    tubeInV->SetTransparency(0);
    
    tpcWallV->AddNode(tubeInV, 0);
    
    //create flanges
    TString tube_flange_name = "tpc01TubeFlange";
    Double_t tube_flange_rmin = tubein_rmax;
    Double_t tube_flange_rmax = tube_flange_rmin + Flange_out_width;

    TGeoTube * TubeFlangeS = new TGeoTube (tube_flange_name, tube_flange_rmin, tube_flange_rmax, Flange_out_thickness/2);
    TGeoVolume *TubeFlangeV = new TGeoVolume(tube_flange_name, TubeFlangeS);
    TubeFlangeV->SetMedium(pMedAluminium);
    TubeFlangeV->SetLineColor(TColor::GetColor("#c7fcec"));

    TGeoCombiTrans *tube_flange_pos1 = new TGeoCombiTrans(0, 0, Chamber_tpc_z / 2 - Flange_out_thickness / 2, NULL);
    TGeoCombiTrans *tube_flange_pos1refl = new TGeoCombiTrans(*tube_flange_pos1);
    tube_flange_pos1refl->ReflectZ(true);
    
    tpcWallV->AddNode(TubeFlangeV, 0, tube_flange_pos1);
    tpcWallV->AddNode(TubeFlangeV, 1, tube_flange_pos1refl);
    
    CreateTubeAlFoilRings(tubeInV, tubein_rmin, tubein_rmax, true);
}

void CreateTubeAlFoilRings(TGeoVolume* mother_volume, Double_t rmin, Double_t rmax, Bool_t is_outer_wall)
{
    TString foil_tube_in_name = "tpc01WallFoilTubeIn";

    TGeoTube * FoilTubeInS = new TGeoTube (foil_tube_in_name, rmin, rmin + Barrel_Al_thick, Chamber_tpc_z/2);
    TGeoVolume *FoilTubeInV = new TGeoVolume(foil_tube_in_name, FoilTubeInS);
    FoilTubeInV->SetMedium(is_outer_wall ? pMedTPCmixture : pMedN2);
    FoilTubeInV->SetLineColor(TColor::GetColor("#a0a0a0"));
    
    mother_volume->AddNode(FoilTubeInV, 0);
    
    TString foil_tube_out_name = "tpc01WallFoilTubeOut";

    TGeoTube * FoilTubeOutS = new TGeoTube (foil_tube_out_name, rmax - Barrel_Al_thick, rmax, Chamber_tpc_z/2);
    TGeoVolume *FoilTubeOutV = new TGeoVolume(foil_tube_out_name, FoilTubeOutS);
    FoilTubeOutV->SetMedium(is_outer_wall ? pMedN2 : pMedTPCmixture);
    FoilTubeOutV->SetLineColor(TColor::GetColor("#a0a0a0"));
    
    mother_volume->AddNode(FoilTubeOutV, 0);
    
    TString foil_ring_in_name = "tpc01WallFoilRingIn";

    TGeoTube * FoilRingInS = new TGeoTube (foil_ring_in_name, rmin, rmin + Barrel_Al_thick, AlFoil_ring_thikness/2);
    TGeoVolume *FoilRingInV = new TGeoVolume(foil_ring_in_name, FoilRingInS);
    FoilRingInV->SetMedium(pMedAluminium);
    FoilRingInV->SetLineColor(TColor::GetColor("#c7fcec"));
    
    FoilTubeInV->AddNode(FoilRingInV, 0);
    
    TString foil_ring_out_name = "tpc01WallFoilRingOut";

    TGeoTube * FoilRingOutS = new TGeoTube (foil_ring_out_name, rmax - Barrel_Al_thick, rmax, AlFoil_ring_thikness/2);
    TGeoVolume *FoilRingOutV = new TGeoVolume(foil_ring_out_name, FoilRingOutS);
    FoilRingOutV->SetMedium(pMedAluminium);
    FoilRingOutV->SetLineColor(TColor::GetColor("#c7fcec"));
    
    FoilTubeOutV->AddNode(FoilRingOutV, 0);
    
    Double_t z_pos = AlFoil_ring_thikness + AlFoil_ring_gap;
    Int_t ind = 1;
    while (z_pos < Chamber_tpc_z/2)
    {
        TGeoCombiTrans *ring_pos1 = new TGeoCombiTrans(0, 0, z_pos, NULL);
        TGeoCombiTrans *ring_pos2 = new TGeoCombiTrans(0, 0, z_pos, NULL);
        TGeoCombiTrans *ring_pos1refl = new TGeoCombiTrans(*ring_pos1);
        ring_pos1refl->ReflectZ(true);
        TGeoCombiTrans *ring_pos2refl = new TGeoCombiTrans(*ring_pos2);
        ring_pos2refl->ReflectZ(true);
        
        FoilTubeInV->AddNode(FoilRingInV, ind, ring_pos1);
        FoilTubeOutV->AddNode(FoilRingOutV, ind, ring_pos2);
        ind++;
        FoilTubeInV->AddNode(FoilRingInV, ind, ring_pos1refl);
        FoilTubeOutV->AddNode(FoilRingOutV, ind, ring_pos2refl);
        ind++;
        
        z_pos += AlFoil_ring_thikness + AlFoil_ring_gap;
    }
}

void CreateMembrane(TGeoVolume* mother_volume)
{
    TString membrane_name = "tpc01mb";
    
    Double_t outer_wall_thickness = Thick_barell_outer_layer[1] +
                                    Thick_barell_CO2_layer[1] +
                                    Thick_barell_inner_layer[1] +
                                    Barrel_Al_thick;

    Double_t inner_wall_thickness = Thick_barell_outer_layer[0] +
                                    Thick_barell_CO2_layer[0] +
                                    Thick_barell_inner_layer[0] +
                                    Barrel_Al_thick;
    
    Double_t membrane_rmin = TpcInnerRadius + inner_wall_thickness; // consider inner wall thickness
    Double_t membrane_rmax = TpcOuterRadius - outer_wall_thickness; // consider outer wall thickness
    
    TGeoTube *tpcMembraneS = new TGeoTube(membrane_name, membrane_rmin, membrane_rmax, Membrane_thickness/2);
    TGeoVolume *tpcMembraneV = new TGeoVolume(membrane_name, tpcMembraneS);
    tpcMembraneV->SetMedium(pMedMylar);
    tpcMembraneV->SetLineColor(TColor::GetColor("#ffc7c7"));

    mother_volume->AddNode(tpcMembraneV, 0);
    
    TString membrane_foil_name = "tpc01mbFoil";
    
    TGeoTube *tpcMembraneFoilS = new TGeoTube(membrane_foil_name, membrane_rmin, membrane_rmax, Membrane_AlFoil_thickness/2);
    TGeoVolume *tpcMembraneFoilV = new TGeoVolume(membrane_foil_name, tpcMembraneFoilS);
    tpcMembraneFoilV->SetMedium(pMedAluminium);
    tpcMembraneFoilV->SetLineColor(TColor::GetColor("#c7fcec"));
    
    tpcMembraneV->AddNode(tpcMembraneFoilV, 0, new TGeoTranslation(0, 0, Membrane_thickness/2 - Membrane_AlFoil_thickness/2));
    tpcMembraneV->AddNode(tpcMembraneFoilV, 1, new TGeoTranslation(0, 0, -Membrane_thickness/2 + Membrane_AlFoil_thickness/2));
}

void CreateInnerStructure(TGeoVolume* mother_volume)
{
    TString module_name = "tpc01mod";

    //make half of module
    Double_t half_module_width = (Chamber_tpc_z - Membrane_thickness)/2;
    Double_t half_module_width_ext = half_module_width + Chamber_tpc_z_ext_parts + ROC_Frame_zspace_iocards;
    
    Double_t outer_wall_thickness = Thick_barell_outer_layer[1] +
                                    Thick_barell_CO2_layer[1] +
                                    Thick_barell_inner_layer[1] +
                                    Barrel_Al_thick;

    Double_t inner_wall_thickness = Thick_barell_outer_layer[0] +
                                    Thick_barell_CO2_layer[0] +
                                    Thick_barell_inner_layer[0] +
                                    Barrel_Al_thick;
    
    Double_t module_rmin = TpcInnerRadius + inner_wall_thickness; // consider inner wall thickness
    Double_t module_rmax = TpcOuterRadius - outer_wall_thickness; // consider outer wall thickness

    TGeoTube *tpcHalfModuleS = new TGeoTube(module_name, module_rmin, module_rmax, half_module_width_ext/2);
    TGeoVolume *tpcHalfModuleV = new TGeoVolume(module_name, tpcHalfModuleS);
    tpcHalfModuleV->SetMedium(pMedTPCmixture);
    tpcHalfModuleV->SetLineColor(TColor::GetColor("#adff2f"));
    tpcHalfModuleV->SetTransparency(90);
    //tpcHalfModuleV->SetVisibility(kFALSE);
    
    TGeoCombiTrans *module_position = new TGeoCombiTrans(0., 0., (Membrane_thickness + half_module_width_ext) / 2, NULL);
    module_position->RotateZ(-Section_step / 2);
    TGeoCombiTrans *module_positionrefl = new TGeoCombiTrans(*module_position);
    module_positionrefl->ReflectZ(true);

    //Attach module to mother volume
    mother_volume->AddNode(tpcHalfModuleV, 0, module_position);
    mother_volume->AddNode(tpcHalfModuleV, 1, module_positionrefl);
    
    CreateMembraneFrame(tpcHalfModuleV, half_module_width, module_rmin, module_rmax);
    
    CreateFieldCage(tpcHalfModuleV, half_module_width);
    
    CreateEndsFrame(tpcHalfModuleV, half_module_width, module_rmin, module_rmax);
    
    CreateSensitiveStructure(tpcHalfModuleV, half_module_width_ext);
}

void CreateMembraneFrame(TGeoVolume* mother_volume, Double_t width, Double_t rmin, Double_t rmax)
{
    TString inner_holder_name = "tpc01mbinHold";
    
    TGeoPgon *MbInHolderPgS = new TGeoPgon(inner_holder_name + "Pg", -Fieldcage_shift_deg, 360., 12, 2);
    MbInHolderPgS->DefineSection(0, -Membrane_holder_thikness/2, rmin, Membrane_inner_holder_R);
    MbInHolderPgS->DefineSection(1, Membrane_holder_thikness/2, rmin, Membrane_inner_holder_R);
    
    TGeoTube *MbInHolderTbS = new TGeoTube(inner_holder_name + "Tb", rmin, Membrane_inner_holder_R, Membrane_holder_thikness/2);
    
    TGeoCompositeShape * MbInHolderS = new TGeoCompositeShape(inner_holder_name, inner_holder_name + "Pg" + "+" + inner_holder_name + "Tb");
    TGeoVolume *MbInHolderV = new TGeoVolume(inner_holder_name, MbInHolderS);
    MbInHolderV->SetMedium(pMedAluminium);
    MbInHolderV->SetLineColor(TColor::GetColor("#c7fcec"));
    
    mother_volume->AddNode(MbInHolderV, 0, new TGeoTranslation(0., 0., (-width - Chamber_tpc_z_ext_parts - ROC_Frame_zspace_iocards + Membrane_holder_thikness) / 2));
    
    TString outer_holder_name = "tpc01mboutHold";
    
    TGeoPgon *MbOutHolderPgS = new TGeoPgon(outer_holder_name + "Pg", -Fieldcage_shift_deg, 360., 12, 2);
    MbOutHolderPgS->DefineSection(0, -Membrane_holder_thikness/2, Membrane_outer_holder_R, rmax * TMath::Cos(Section_rad_step / 2.));
    MbOutHolderPgS->DefineSection(1, Membrane_holder_thikness/2, Membrane_outer_holder_R, rmax * TMath::Cos(Section_rad_step / 2.));
    
    TGeoTube *MbOutHolderTbS = new TGeoTube(outer_holder_name + "Tb", Membrane_outer_holder_R_edge, rmax, Membrane_holder_thikness/2);
    
    TGeoCompositeShape * MbOutHolderS = new TGeoCompositeShape(outer_holder_name, outer_holder_name + "Tb" + "+" + outer_holder_name + "Pg");
    TGeoVolume *MbOutHolderV = new TGeoVolume(outer_holder_name, MbOutHolderS);
    MbOutHolderV->SetMedium(pMedAluminium);
    MbOutHolderV->SetLineColor(TColor::GetColor("#c7fcec"));
    
    mother_volume->AddNode(MbOutHolderV, 0, new TGeoTranslation(0., 0., (-width - Chamber_tpc_z_ext_parts - ROC_Frame_zspace_iocards + Membrane_holder_thikness) / 2));
    
    TString out_pin_name = "tpc01mboutpin";

    TGeoTube *tpcOutPinS = new TGeoTube(out_pin_name, Fieldcage_Pin_R[1] - Fieldcage_Pin_wallthik[1], Fieldcage_Pin_R[1], Membrane_holder_thikness/2);
    TGeoVolume *tpcOutPinV = new TGeoVolume(out_pin_name, tpcOutPinS);
    tpcOutPinV->SetLineColor(TColor::GetColor("#c7fcec"));
    tpcOutPinV->SetMedium(pMedAluminium);
    
    for(Int_t isec = 0; isec < Nsections; isec++)
    {
        TGeoCombiTrans *out_pin_position= new TGeoCombiTrans();
        out_pin_position->SetTranslation(Membrane_outer_holder_R_edge - Fieldcage_Pin_R[1] / TMath::Cos(Section_rad_step / 2.), 0, (-width - Chamber_tpc_z_ext_parts - ROC_Frame_zspace_iocards + Membrane_holder_thikness) / 2);
        out_pin_position->RotateZ(Section_step * isec - Fieldcage_shift_deg);

        mother_volume->AddNode(tpcOutPinV, isec+1, out_pin_position);
    }
}

void CreateFieldCage(TGeoVolume *mother_volume, Double_t vol_width)
{
    Double_t cage_width = vol_width - Membrane_holder_thikness - End_Frame_thickness;
    
    TString out_pin_name = "tpc01outpin";
    TString in_pin_name = "tpc01inpin";

    TGeoTube *tpcOutPinS = new TGeoTube(out_pin_name, Fieldcage_Pin_R[1] - Fieldcage_Pin_wallthik[1], Fieldcage_Pin_R[1], cage_width/2);
    TGeoVolume *tpcOutPinV = new TGeoVolume(out_pin_name, tpcOutPinS);
    tpcOutPinV->SetLineColor(TColor::GetColor("#edaa21"));
    tpcOutPinV->SetMedium(pMedKevlar);
    
    TGeoTube *tpcOutPinPPS = new TGeoTube(out_pin_name + "pp", Fieldcage_Pin_R[1] - Fieldcage_Pin_wallthik[0], Fieldcage_Pin_R[1], cage_width/2);
    TGeoVolume *tpcOutPinPPV = new TGeoVolume(out_pin_name + "pp", tpcOutPinPPS);
    tpcOutPinPPV->SetLineColor(TColor::GetColor("#edaa21"));
    tpcOutPinPPV->SetMedium(pMedPolypropylene);
    
    tpcOutPinV->AddNode(tpcOutPinPPV, 0);
    
    TGeoTube *tpcInPinS = new TGeoTube(in_pin_name, Fieldcage_Pin_R[0] - Fieldcage_Pin_wallthik[0], Fieldcage_Pin_R[0], cage_width/2);
    TGeoVolume *tpcInPinV = new TGeoVolume(in_pin_name, tpcInPinS);
    tpcInPinV->SetLineColor(TColor::GetColor("#edaa21"));
    tpcInPinV->SetMedium(pMedPolypropylene);
    
    for(Int_t isec = 0; isec < Nsections; isec++)
    {
        TGeoCombiTrans *out_pin_position= new TGeoCombiTrans();
        out_pin_position->SetTranslation(Membrane_outer_holder_R_edge - Fieldcage_Pin_R[1] / TMath::Cos(Section_rad_step / 2.), 0, (Membrane_holder_thikness - Chamber_tpc_z_ext_parts - End_Frame_thickness - ROC_Frame_zspace_iocards) / 2);
        out_pin_position->RotateZ(Section_step * isec - Fieldcage_shift_deg);
        TGeoCombiTrans *in_pin_position= new TGeoCombiTrans();
        in_pin_position->SetTranslation(Membrane_inner_holder_R_edge - Fieldcage_Pin_R[0] / TMath::Cos(Section_rad_step / 2.), 0, (Membrane_holder_thikness - Chamber_tpc_z_ext_parts - End_Frame_thickness - ROC_Frame_zspace_iocards) / 2);
        in_pin_position->RotateZ(Section_step * isec - Fieldcage_shift_deg);

        mother_volume->AddNode(tpcOutPinV, isec+1, out_pin_position);
        mother_volume->AddNode(tpcInPinV, isec+1, in_pin_position);
    }
    
    TString FieldCageOut_name = "tpc01outfc";
    TString FieldCageIn_name = "tpc01infc";
    
    TGeoPgon *tpcFieldCageOutS = new TGeoPgon(FieldCageOut_name, -Fieldcage_shift_deg, 360., 12, 2);
    tpcFieldCageOutS->DefineSection(0, 0., Membrane_outer_holder_R, Membrane_outer_holder_R + Fieldcage_ribbon_thikness);
    tpcFieldCageOutS->DefineSection(1, Fieldcade_ribbon_wide, Membrane_outer_holder_R, Membrane_outer_holder_R + Fieldcage_ribbon_thikness);

    TGeoPgon *tpcFieldCageInS = new TGeoPgon(FieldCageIn_name, -Fieldcage_shift_deg, 360., 12, 2);
    tpcFieldCageInS->DefineSection(0, 0., Membrane_inner_holder_R, Membrane_inner_holder_R + Fieldcage_ribbon_thikness);
    tpcFieldCageInS->DefineSection(1, Fieldcade_ribbon_wide, Membrane_inner_holder_R, Membrane_inner_holder_R + Fieldcage_ribbon_thikness);

    TGeoVolume *tpcFieldCageOutV = new TGeoVolume(FieldCageOut_name, tpcFieldCageOutS);
    tpcFieldCageOutV->SetMedium(pMedMylar);
    TGeoVolume *tpcFieldCageInV = new TGeoVolume(FieldCageIn_name, tpcFieldCageInS);
    tpcFieldCageInV->SetMedium(pMedMylar);

    tpcFieldCageOutV->SetLineColor(TColor::GetColor("#edff21"));
    tpcFieldCageInV->SetLineColor(TColor::GetColor("#edff21"));

    Double_t rib_posZ = -vol_width / 2 + Membrane_holder_thikness + Fieldcage_ribbon_gap;
    Int_t rib_ind = 0;
    while (rib_posZ < vol_width / 2 - End_Frame_thickness - Fieldcade_ribbon_wide)
    {
        mother_volume->AddNode(tpcFieldCageOutV, rib_ind, new TGeoTranslation(0., 0., rib_posZ - (Chamber_tpc_z_ext_parts + ROC_Frame_zspace_iocards) / 2));
        mother_volume->AddNode(tpcFieldCageInV, rib_ind, new TGeoTranslation(0., 0., rib_posZ - (Chamber_tpc_z_ext_parts + ROC_Frame_zspace_iocards) / 2));
        
        rib_posZ += Fieldcade_ribbon_wide + Fieldcage_ribbon_gap;
        rib_ind++;
    }
}

void CreateEndsFrame(TGeoVolume *mother_volume, Double_t width, Double_t rmin, Double_t rmax)
{
    //create flanges over inner wall
    TString over_flange_name = "tpc01FlangeOverWall";
    Double_t over_flange_rmin = rmin;
    Double_t over_flange_rmax = over_flange_rmin + Flange_over_inner_tube_width;

    TGeoTube * OverFlangeS = new TGeoTube (over_flange_name, over_flange_rmin, over_flange_rmax, Flange_over_inner_tube_thickness/2);
    TGeoVolume *OverFlangeV = new TGeoVolume(over_flange_name, OverFlangeS);
    OverFlangeV->SetMedium(pMedFiberGlass);
    OverFlangeV->SetLineColor(TColor::GetColor("#c7fcec"));

    mother_volume->AddNode(OverFlangeV, 0, new TGeoTranslation(0, 0, (width - Chamber_tpc_z_ext_parts - Flange_over_inner_tube_thickness - ROC_Frame_zspace_iocards) / 2));
    
    TString inner_frame_name = "tpc01EndFrameIn";
    
    TGeoPgon *InnerFramePgS = new TGeoPgon(inner_frame_name + "Pg", 0., 360., 12, 2);
    InnerFramePgS->DefineSection(0, -End_Frame_thickness/2, End_Frame_in_R, End_Frame_in_R_seg);
    InnerFramePgS->DefineSection(1, End_Frame_thickness/2, End_Frame_in_R, End_Frame_in_R_seg);
    
    TGeoTube *InnerFrameTbS = new TGeoTube(inner_frame_name + "Tb", End_Frame_in_R, End_Frame_in_R_seg, End_Frame_thickness/2);
    
    TGeoCompositeShape * InnerFrameS = new TGeoCompositeShape(inner_frame_name, inner_frame_name + "Pg" + "+" + inner_frame_name + "Tb");
    TGeoVolume *InnerFrameV = new TGeoVolume(inner_frame_name, InnerFrameS);
    InnerFrameV->SetMedium(pMedAluminium);
    InnerFrameV->SetLineColor(TColor::GetColor("#c7fcec"));
    
    mother_volume->AddNode(InnerFrameV, 0, new TGeoTranslation(0., 0., (width - Chamber_tpc_z_ext_parts - ROC_Frame_zspace_iocards - End_Frame_thickness) / 2));
    
    TString inner_air_name = "tpc01EndAirIn";
    
    TGeoPgon *InnerAirPgS = new TGeoPgon(inner_air_name + "Pg", 0., 360., 12, 2);
    InnerAirPgS->DefineSection(0, -(Chamber_tpc_z_ext_parts + ROC_Frame_zspace_iocards) / 2., End_Frame_in_R, End_Frame_in_R_seg);
    InnerAirPgS->DefineSection(1, (Chamber_tpc_z_ext_parts + ROC_Frame_zspace_iocards) / 2., End_Frame_in_R, End_Frame_in_R_seg);
    
    TGeoTube *InnerAirTbS = new TGeoTube(inner_air_name + "Tb", rmin, End_Frame_in_R_seg, (Chamber_tpc_z_ext_parts + ROC_Frame_zspace_iocards) / 2.);
    
    TGeoCompositeShape * InnerAirS = new TGeoCompositeShape(inner_air_name, inner_air_name + "Pg" + "+" + inner_air_name + "Tb");
    TGeoVolume *InnerAirV = new TGeoVolume(inner_air_name, InnerAirS);
    InnerAirV->SetMedium(pMedAir);
    
    mother_volume->AddNode(InnerAirV, 0, new TGeoTranslation(0., 0., width / 2.));
    
    {
        TString frameIn_name = "tpc01ROCPressOutFrAirIn";
        
        TGeoVolume * tpcROCPressOutFrAirInV = new TGeoVolumeAssembly(frameIn_name);
        
        Int_t npoints_structure = 6;
        unique_ptr<Double_t[]> x(new Double_t[npoints_structure]);
        unique_ptr<Double_t[]> y(new Double_t[npoints_structure]);
        
        Double_t inR = ROCPress_Backplate_inR - (ROCPress_Backplate_GapToNextSect / 2) / TMath::Sin(Section_rad_step / 2.);
        Double_t x_shift = ROCPress_OutFrame_width / TMath::Sin(Section_rad_step / 2);
        Double_t angle = TMath::ACos(ROCPress_Backplate_HighExtX / ROCPress_Backplate_HighExtL[1]);
        Double_t angle_out = angle + TMath::PiOver2();
        Double_t angle_in = TMath::Pi() - angle - Section_rad_step / 2;
        
        x[0] = x[5] = End_Frame_in_R_seg - (ROCPress_Backplate_GapToNextSect / 2) / TMath::Sin(Section_rad_step / 2.);
        x[1] = inR * (1 + TMath::Tan(Section_rad_step / 2) / (TMath::Tan(ROCPress_OutFrame_LowCenterAngle_rad / 2) - TMath::Tan(Section_rad_step / 2)));
        x[2] = inR;
        x[3] = inR + ROCPress_OutFrame_width;
        x[4] = (inR - x_shift + ROCPress_OutFrame_width / TMath::Sin(ROCPress_OutFrame_LowCenterAngle_rad / 2)) * (1 + TMath::Tan(Section_rad_step / 2) / (TMath::Tan(ROCPress_OutFrame_LowCenterAngle_rad / 2) - TMath::Tan(Section_rad_step / 2))) + x_shift;
        y[2] = y[3] = 0.;
        y[1] = x[1] * TMath::Tan(Section_rad_step / 2);
        y[4] = (x[4] - x_shift) * TMath::Tan(Section_rad_step / 2);
        y[0] = x[0] * TMath::Tan(Section_rad_step / 2);
        y[5] = (x[5] - x_shift) * TMath::Tan(Section_rad_step / 2);
        
        TGeoXtru *ROCPressOutFrameS = new TGeoXtru(2);
        ROCPressOutFrameS->DefinePolygon(npoints_structure, x.get(), y.get());
        ROCPressOutFrameS->DefineSection(0, 0.);
        ROCPressOutFrameS->DefineSection(1, ROCPress_OutFrame_thick);
        
        TGeoVolume * ROCPressOutFrameV = new TGeoVolume(frameIn_name, ROCPressOutFrameS);
        ROCPressOutFrameV->SetMedium(pMedAluminium);
        ROCPressOutFrameV->SetLineColor(kMagenta);
        
        TGeoCombiTrans * fr_pos = new TGeoCombiTrans();
        fr_pos->SetTranslation((ROCPress_Backplate_GapToNextSect / 2) / TMath::Sin(Section_rad_step / 2.), 0., ROCPress_Backplate_thick - (Chamber_tpc_z_ext_parts + ROC_Frame_zspace_iocards) / 2.);
        TGeoCombiTrans * fr_pos_refl = new TGeoCombiTrans(*fr_pos);
        fr_pos_refl->ReflectY(true);
        
        tpcROCPressOutFrAirInV->AddNode(ROCPressOutFrameV, 0, fr_pos);
        tpcROCPressOutFrAirInV->AddNode(ROCPressOutFrameV, 1, fr_pos_refl);
        
        for(Int_t isec = 0; isec < Nsections; isec++) //Nsections
        {
            TGeoCombiTrans *pos= new TGeoCombiTrans();
            pos->SetTranslation(0., 0., 0.);
            pos->RotateZ(Section_step * isec + Section_step / 2);

            InnerAirV->AddNode(tpcROCPressOutFrAirInV, isec, pos);
        }
    }
    
    TString outer_frame_name = "tpc01EndFrameOut";
    
    TGeoPgon *OuterFramePgS = new TGeoPgon(outer_frame_name + "Pg", 0., 360., 12, 2);
    OuterFramePgS->DefineSection(0, -End_Frame_thickness/2, End_Frame_out_R_seg, End_Frame_out_R * TMath::Cos(Section_rad_step / 2.));
    OuterFramePgS->DefineSection(1, End_Frame_thickness/2, End_Frame_out_R_seg, End_Frame_out_R * TMath::Cos(Section_rad_step / 2.));
    
    TGeoTube *OuterFrameTbS = new TGeoTube(outer_frame_name + "Tb", End_Frame_out_R_seg / TMath::Cos(Section_rad_step / 2.), End_Frame_out_R, End_Frame_thickness/2);
    
    TGeoCompositeShape * OuterFrameS = new TGeoCompositeShape(outer_frame_name, outer_frame_name + "Tb" + "+" + outer_frame_name + "Pg");
    TGeoVolume *OuterFrameV = new TGeoVolume(outer_frame_name, OuterFrameS);
    OuterFrameV->SetMedium(pMedAluminium);
    OuterFrameV->SetLineColor(TColor::GetColor("#c7fcec"));
    
    mother_volume->AddNode(OuterFrameV, 0, new TGeoTranslation(0., 0., (width - Chamber_tpc_z_ext_parts - ROC_Frame_zspace_iocards - End_Frame_thickness) / 2));
    
    TString outer_air_name = "tpc01EndFrameOutAir";
    
    TGeoPgon *OuterAirPgS = new TGeoPgon(outer_air_name + "Pg", 0, 360, 12, 2);
    OuterAirPgS->DefineSection(0, -(Chamber_tpc_z_ext_parts + ROC_Frame_zspace_iocards) / 2., End_Frame_out_R_seg, End_Frame_out_R * TMath::Cos(Section_rad_step / 2.));
    OuterAirPgS->DefineSection(1, (Chamber_tpc_z_ext_parts + ROC_Frame_zspace_iocards) / 2., End_Frame_out_R_seg, End_Frame_out_R * TMath::Cos(Section_rad_step / 2.));
    
    TGeoTube *OuterAirTbS = new TGeoTube(outer_air_name + "Tb", End_Frame_out_R_seg / TMath::Cos(Section_rad_step / 2.), rmax, (Chamber_tpc_z_ext_parts + ROC_Frame_zspace_iocards) / 2.);
    
    TGeoCompositeShape * OuterAirS = new TGeoCompositeShape(outer_air_name, outer_air_name + "Tb" + "+" + outer_air_name + "Pg");
    TGeoVolume *OuterAirV = new TGeoVolume(outer_air_name, OuterAirS);
    OuterAirV->SetMedium(pMedAir);
    
    mother_volume->AddNode(OuterAirV, 0, new TGeoTranslation(0., 0., width / 2.));
    
    {
        TString frame_name = "tpc01ROCPressOutFrAirOut";
        
        TGeoVolume * tpcROCPressOutFrAirOutV = new TGeoVolumeAssembly(frame_name);
        
        Int_t npoints_structure = 12;
        unique_ptr<Double_t[]> x(new Double_t[npoints_structure]);
        unique_ptr<Double_t[]> y(new Double_t[npoints_structure]);
        
        Double_t outR = ROCPress_Backplate_outR - (ROCPress_Backplate_GapToNextSect / 2) / TMath::Sin(Section_rad_step / 2.);
        Double_t x_shift = ROCPress_OutFrame_width / TMath::Sin(Section_rad_step / 2);
        Double_t angle = TMath::ACos(ROCPress_Backplate_HighExtX / ROCPress_Backplate_HighExtL[1]);
        Double_t angle_out = angle + TMath::PiOver2();
        Double_t angle_in = TMath::Pi() - angle - Section_rad_step / 2;
        x[0] = x[1] = outR - ROCPress_Backplate_HighExtX;
        x[2] = x[3] = outR;
        x[4] = outR - ROCPress_Backplate_HighExtL[3] * TMath::Cos(angle);
        x[5] = x[6] = End_Frame_out_R_seg - (ROCPress_Backplate_GapToNextSect / 2) / TMath::Sin(Section_rad_step / 2.);
        x[7] = x[4] - ROCPress_OutFrame_width / TMath::Sin(angle_in / 2) * TMath::Cos(angle_in / 2 + Section_rad_step / 2);
        x[8] = x[9] = outR - ROCPress_OutFrame_width;
        x[10] = x[11] = x[0] - ROCPress_OutFrame_width;    
        y[0] = y[11] = 0.;
        y[4] = x[4] * TMath::Tan(Section_rad_step / 2);
        y[5] = x[5] * TMath::Tan(Section_rad_step / 2);
        y[6] = (x[6] - x_shift) * TMath::Tan(Section_rad_step / 2);
        y[7] = (x[7] - x_shift) * TMath::Tan(Section_rad_step / 2);
        y[1] = ROCPress_Backplate_HighExtL[0];
        y[2] = y[1] + ROCPress_Backplate_HighExtL[1] * TMath::Sin(angle);
        y[3] = y[2] + ROCPress_Backplate_HighExtL[2];
        y[10] = y[1] + ROCPress_OutFrame_width / TMath::Tan(angle_out / 2);
        y[9] = y[10] + ROCPress_Backplate_HighExtL[1] * TMath::Sin(angle);
        y[8] = y[9] + ROCPress_Backplate_HighExtL[2] - 2 * ROCPress_OutFrame_width / TMath::Tan(angle_out / 2);

        TGeoXtru *ROCPressOutFrameS = new TGeoXtru(2);
        ROCPressOutFrameS->DefinePolygon(npoints_structure, x.get(), y.get());
        ROCPressOutFrameS->DefineSection(0, 0.);
        ROCPressOutFrameS->DefineSection(1, ROCPress_OutFrame_thick);
        
        TGeoVolume * ROCPressOutFrameV = new TGeoVolume(frame_name, ROCPressOutFrameS);
        ROCPressOutFrameV->SetMedium(pMedAluminium);
        ROCPressOutFrameV->SetLineColor(kMagenta);
        
        TGeoCombiTrans * fr_pos = new TGeoCombiTrans();
        fr_pos->SetTranslation((ROCPress_Backplate_GapToNextSect / 2) / TMath::Sin(Section_rad_step / 2.), 0., ROCPress_Backplate_thick - (Chamber_tpc_z_ext_parts + ROC_Frame_zspace_iocards) / 2.);
        TGeoCombiTrans * fr_pos_refl = new TGeoCombiTrans(*fr_pos);
        fr_pos_refl->ReflectY(true);
        
        tpcROCPressOutFrAirOutV->AddNode(ROCPressOutFrameV, 0, fr_pos);
        tpcROCPressOutFrAirOutV->AddNode(ROCPressOutFrameV, 1, fr_pos_refl);
        
        for(Int_t isec = 0; isec < Nsections; isec++) //Nsections
        {
            TGeoCombiTrans *pos= new TGeoCombiTrans();
            pos->SetTranslation(0., 0., 0.);
            pos->RotateZ(Section_step * isec + Section_step / 2);

            OuterAirV->AddNode(tpcROCPressOutFrAirOutV, isec, pos);
        }
    }
}

void CreateSensitiveStructure(TGeoVolume *mother_volume, Double_t vol_width)
{
    TString space_name = "tpc01bp";
    Double_t space_width = Pad_Plane_Pp + Pad_Plane_G10 + ROC_Frame;

    TGeoPgon *ROCSpaceS = new TGeoPgon(space_name, -Section_step / 2, Section_step, 1, 2);
    ROCSpaceS->DefineSection(0, -space_width/2, End_Frame_in_R_seg, End_Frame_out_R_seg);
    ROCSpaceS->DefineSection(1, space_width/2 + ROC_Frame_zspace_iocards, End_Frame_in_R_seg, End_Frame_out_R_seg);
    
    TGeoVolume *ROCSpaceV = new TGeoVolume(space_name, ROCSpaceS);
    ROCSpaceV->SetMedium(pMedAir);
    
    for(Int_t isec = 0; isec < Nsections; isec++) //Nsections
    {
        TGeoCombiTrans *roc_pos= new TGeoCombiTrans();
        roc_pos->SetTranslation(0., 0., (vol_width - space_width) / 2 - ROC_Frame_zspace_iocards);
        roc_pos->RotateZ(Section_step * isec + Section_step / 2);

        mother_volume->AddNode(ROCSpaceV, isec, roc_pos);
    }
    
    CreatePadWithFrameStructure(ROCSpaceV, space_width);
    
    CreateSensitiveVolumes(mother_volume, vol_width);
}

void CreateSensitiveVolumes(TGeoVolume *mother_volume, Double_t vol_width)
{
    Double_t space_width = Pad_Plane_Pp + Pad_Plane_G10 + ROC_Frame;
    
    TString SensitiveVolume_name = "tpc01sv";
    
    TGeoPgon *SensitiveVolumeS = new TGeoPgon(SensitiveVolume_name, 0., 360., 12, 2);
    SensitiveVolumeS->DefineSection(0, -(vol_width - space_width - ROC_Frame_zspace_iocards)/2, Sensitive_in_R_seg, Sensitive_out_R_seg);
    SensitiveVolumeS->DefineSection(1, (vol_width - space_width - ROC_Frame_zspace_iocards)/2, Sensitive_in_R_seg, Sensitive_out_R_seg);
    
    std::cout << "Drift length " << vol_width - space_width << std::endl;
    
    TGeoVolume *SensitiveVolumeV = new TGeoVolume(SensitiveVolume_name, SensitiveVolumeS);
    SensitiveVolumeV->SetMedium(pMedTPCmixture);
    SensitiveVolumeV->SetLineColor(TColor::GetColor("#ff9494"));
    
    mother_volume->AddNode(SensitiveVolumeV, 1, new TGeoTranslation(0., 0., -(space_width + ROC_Frame_zspace_iocards) / 2));
}

void CreatePadWithFrameStructure(TGeoVolume *mother_volume, Double_t vol_width)
{
    //end frame spokes
    TString EndSpoke_name = "tpc01EndFrameSpoke";
    
    Int_t EndSpokePts = 4;
    unique_ptr<Double_t[]> xEndSpoke(new Double_t[EndSpokePts]);
    unique_ptr<Double_t[]> yEndSpoke(new Double_t[EndSpokePts]);
    yEndSpoke[0] = yEndSpoke[3] = 0.0;
    yEndSpoke[1] = yEndSpoke[2] = End_Frame_Spoke_width / 2;
    xEndSpoke[0] = End_Frame_in_R_seg / TMath::Cos(Section_rad_step / 2);
    xEndSpoke[1] = xEndSpoke[0] - End_Frame_Spoke_width / 2 * TMath::Tan(Section_rad_step / 2);
    xEndSpoke[3] = End_Frame_out_R_seg / TMath::Cos(Section_rad_step / 2);
    xEndSpoke[2] = xEndSpoke[3] - End_Frame_Spoke_width / 2 * TMath::Tan(Section_rad_step / 2);
    
    TGeoXtru *EndSpokeS = new TGeoXtru(2);
    EndSpokeS->DefinePolygon(EndSpokePts, xEndSpoke.get(), yEndSpoke.get());
    EndSpokeS->DefineSection(0, -End_Frame_Spoke_thick / 2);
    EndSpokeS->DefineSection(1, End_Frame_Spoke_thick / 2);
    
    TGeoVolume *EndSpokeV = new TGeoVolume(EndSpoke_name, EndSpokeS);
    EndSpokeV->SetMedium(pMedAluminium);
    EndSpokeV->SetLineColor(TColor::GetColor("#c7fcec"));
    
    TGeoCombiTrans * sp_pos = new TGeoCombiTrans();
    sp_pos->SetTranslation(0., 0., (vol_width - End_Frame_Spoke_thick) / 2 - Chamber_tpc_z_ext_parts);
    sp_pos->RotateZ(-Section_step / 2);
    TGeoCombiTrans * sp_pos_refl = new TGeoCombiTrans(*sp_pos);
    sp_pos_refl->ReflectY(true);
    
    mother_volume->AddNode(EndSpokeV, 0, sp_pos);
    mother_volume->AddNode(EndSpokeV, 1, sp_pos_refl);
    
    //pad + insulation + Al planes
    TString bpPp_name = "tpc01bpPp";
    TString bpG10_name = "tpc01bpG10";
    TString bpAl_name = "tpc01bpAl";
    
    Double_t x_offset = (End_Frame_Spoke_width / 2 + ROC_Frame_Gap_to_End_Frame_Side - ROC_Frame_ext_width) / TMath::Sin(Section_rad_step / 2);
    Double_t R_seg_in = Sensitive_in_R_seg - x_offset;
    Double_t R_seg_out = Sensitive_out_R_seg - x_offset;
    
    TGeoPgon *PpS = new TGeoPgon(bpPp_name, -Section_step / 2, Section_step, 1, 2);
    PpS->DefineSection(0, -Pad_Plane_Pp/2, R_seg_in, R_seg_out);
    PpS->DefineSection(1, Pad_Plane_Pp/2, R_seg_in, R_seg_out);
    
    TGeoVolume *PpV = new TGeoVolume(bpPp_name, PpS);
    PpV->SetMedium(pMedTPCmixture);
    PpV->SetLineColor(TColor::GetColor("#ff9494"));
    
    mother_volume->AddNode(PpV, 0, new TGeoTranslation(x_offset, 0., (-vol_width + Pad_Plane_Pp) / 2));
    
    CreatePadPlaneWithPads(PpV, R_seg_in, R_seg_out);
    
    TGeoPgon *InsulS = new TGeoPgon(bpG10_name, -Section_step / 2, Section_step, 1, 2);
    InsulS->DefineSection(0, -Pad_Plane_G10/2, R_seg_in, R_seg_out);
    InsulS->DefineSection(1, Pad_Plane_G10/2, R_seg_in, R_seg_out);
    
    TGeoVolume *InsulV = new TGeoVolume(bpG10_name, InsulS);
    InsulV->SetMedium(pMedG10);
    InsulV->SetLineColor(TColor::GetColor("#34c924"));
    
    mother_volume->AddNode(InsulV, 0, new TGeoTranslation(x_offset, 0., (-vol_width + Pad_Plane_G10) / 2 + Pad_Plane_Pp));
    
    TGeoPgon *AlPlaneS = new TGeoPgon(bpAl_name, -Section_step / 2, Section_step, 1, 2);
    AlPlaneS->DefineSection(0, -ROC_Frame_plane/2, R_seg_in, R_seg_out);
    AlPlaneS->DefineSection(1, ROC_Frame_plane/2, R_seg_in, R_seg_out);
    
    TGeoVolume *AlPlaneV = new TGeoVolume(bpAl_name, AlPlaneS);
    AlPlaneV->SetMedium(pMedAluminium);
    AlPlaneV->SetLineColor(TColor::GetColor("#84c3be")); //c4c3be
    
    mother_volume->AddNode(AlPlaneV, 0, new TGeoTranslation(x_offset, 0., (-vol_width + ROC_Frame_plane) / 2 + Pad_Plane_Pp + Pad_Plane_G10));
    
    MakeROCFramesHoles(InsulV, AlPlaneV);
    
    //Al ROC frame
    TString frame_structure_name = "tpcFrameSkeleton";
    
    const Int_t npoints_structure = 12;
    unique_ptr<Double_t[]> xFrStruct(new Double_t[npoints_structure]);
    unique_ptr<Double_t[]> yFrStruct(new Double_t[npoints_structure]);
    FillFramesPoints(xFrStruct, yFrStruct);
    
    TGeoXtru *tpcFrameStructureS = new TGeoXtru(2);
    tpcFrameStructureS->DefinePolygon(npoints_structure, xFrStruct.get(), yFrStruct.get());
    tpcFrameStructureS->DefineSection(0, -ROC_Frame_wo_plane);
    tpcFrameStructureS->DefineSection(1, 0.);
    
    TGeoVolume *tpcFrameStructureV = new TGeoVolume(frame_structure_name, tpcFrameStructureS);
    tpcFrameStructureV->SetMedium(pMedAluminium);
    tpcFrameStructureV->SetLineColor(TColor::GetColor("#84c3be"));
    
    TGeoCombiTrans * fr_pos = new TGeoCombiTrans(0., 0., vol_width / 2, NULL);
    TGeoCombiTrans * fr_pos_refl = new TGeoCombiTrans(*fr_pos);
    fr_pos_refl->ReflectY(true);
    
    mother_volume->AddNode(tpcFrameStructureV, 0, fr_pos);
    mother_volume->AddNode(tpcFrameStructureV, 1, fr_pos_refl);
    
    TString structure_ext_name = "tpcFrameSkeletonExt";
    
    const Int_t npoints_structure_ext = 4;
    unique_ptr<Double_t[]> xStructExt(new Double_t[npoints_structure]);
    unique_ptr<Double_t[]> yStructExt(new Double_t[npoints_structure]);
    xStructExt[0] = xStructExt[1] = xFrStruct[1];
    xStructExt[2] = xStructExt[3] = xFrStruct[2];
    yStructExt[0] = yFrStruct[1];
    yStructExt[1] = yFrStruct[1] + ROC_Frame_ext_width / TMath::Cos(Section_rad_step / 2);
    yStructExt[2] = yFrStruct[2] + ROC_Frame_ext_width / TMath::Cos(Section_rad_step / 2);
    yStructExt[3] = yFrStruct[2];
    
    TGeoXtru *tpcFrameStructureExtS = new TGeoXtru(2);
    tpcFrameStructureExtS->DefinePolygon(npoints_structure_ext, xStructExt.get(), yStructExt.get());
    tpcFrameStructureExtS->DefineSection(0, 0.);
    tpcFrameStructureExtS->DefineSection(1, ROC_Frame_ext_thick);
    
    TGeoVolume *tpcFrameStructureExtV = new TGeoVolume(structure_ext_name, tpcFrameStructureExtS);
    tpcFrameStructureExtV->SetMedium(pMedAluminium);
    tpcFrameStructureExtV->SetLineColor(TColor::GetColor("#84c3be"));
    
    TGeoCombiTrans * fr_ext_pos = new TGeoCombiTrans(0., 0., vol_width / 2 - ROC_Frame_wo_plane, NULL);
    TGeoCombiTrans * fr_ext_pos_refl = new TGeoCombiTrans(*fr_ext_pos);
    fr_ext_pos_refl->ReflectY(true);
    
    mother_volume->AddNode(tpcFrameStructureExtV, 0, fr_ext_pos);
    mother_volume->AddNode(tpcFrameStructureExtV, 1, fr_ext_pos_refl);
    
    CreateROCPressingInFrame(mother_volume, vol_width / 2., ROCPress_Backplate_inR + ROCPress_Backplate_LowWidth - ROCPress_InFrame_width, ROCPress_Backplate_outR - ROCPress_Backplate_HighWidth + ROCPress_InFrame_width);
    
    {
        TString frame_name = "tpc01ROCPressOutFr";
        
        Int_t npoints_structure = 4;
        unique_ptr<Double_t[]> x(new Double_t[npoints_structure]);
        unique_ptr<Double_t[]> y(new Double_t[npoints_structure]);
        
        Double_t inR = ROCPress_Backplate_inR - (ROCPress_Backplate_GapToNextSect / 2) / TMath::Sin(Section_rad_step / 2.);
        Double_t x_shift = ROCPress_OutFrame_width / TMath::Sin(Section_rad_step / 2);
        Double_t angle = TMath::ACos(ROCPress_Backplate_HighExtX / ROCPress_Backplate_HighExtL[1]);
        Double_t angle_out = angle + TMath::PiOver2();
        Double_t angle_in = TMath::Pi() - angle - Section_rad_step / 2;
        
        x[0] = x[3] = End_Frame_out_R_seg - (ROCPress_Backplate_GapToNextSect / 2) / TMath::Sin(Section_rad_step / 2.);
        x[1] = x[2] = End_Frame_in_R_seg - (ROCPress_Backplate_GapToNextSect / 2) / TMath::Sin(Section_rad_step / 2.);
        y[0] = x[0] * TMath::Tan(Section_rad_step / 2);
        y[1] = x[1] * TMath::Tan(Section_rad_step / 2);
        y[2] = (x[2] - x_shift) * TMath::Tan(Section_rad_step / 2);
        y[3] = (x[3] - x_shift) * TMath::Tan(Section_rad_step / 2);
        
        TGeoXtru *ROCPressOutFrameS = new TGeoXtru(2);
        ROCPressOutFrameS->DefinePolygon(npoints_structure, x.get(), y.get());
        ROCPressOutFrameS->DefineSection(0, 0.);
        ROCPressOutFrameS->DefineSection(1, ROCPress_OutFrame_thick);
        
        TGeoVolume * ROCPressOutFrameV = new TGeoVolume(frame_name, ROCPressOutFrameS);
        ROCPressOutFrameV->SetMedium(pMedAluminium);
        ROCPressOutFrameV->SetLineColor(kMagenta);
        
        TGeoCombiTrans * fr_pos = new TGeoCombiTrans();
        fr_pos->SetTranslation((ROCPress_Backplate_GapToNextSect / 2) / TMath::Sin(Section_rad_step / 2.), 0., vol_width / 2. - Chamber_tpc_z_ext_parts + ROCPress_Backplate_thick);
        TGeoCombiTrans * fr_pos_refl = new TGeoCombiTrans(*fr_pos);
        fr_pos_refl->ReflectY(true);
        
        mother_volume->AddNode(ROCPressOutFrameV, 0, fr_pos);
        mother_volume->AddNode(ROCPressOutFrameV, 1, fr_pos_refl);
    }
    
    
    CreateElectronicsModules(mother_volume);
}

void CreatePadPlaneWithPads(TGeoVolume *mother_volume, Double_t R_seg_in, Double_t R_seg_out)
{
    TString PpBase_name = "PpBase";
    
    TGeoPgon *PpS = new TGeoPgon(PpBase_name, -Section_step / 2., Section_step, 1, 2);
    PpS->DefineSection(0, -(Pad_Plane_Pp - Pad_Thick)/2, R_seg_in, R_seg_out);
    PpS->DefineSection(1, (Pad_Plane_Pp - Pad_Thick)/2, R_seg_in, R_seg_out);
    
    TGeoVolume *PpV = new TGeoVolume(PpBase_name, PpS);
    PpV->SetMedium(pMedG10);
    PpV->SetLineColor(TColor::GetColor("#ff2400"));
    
    mother_volume->AddNode(PpV, 0, new TGeoTranslation(0. , 0., Pad_Thick / 2));
    
    TString PadIn_name = "PadIn";
    TString PadOut_name = "PadOut";
    
    TGeoBBox * PadIn = new TGeoBBox(PadIn_name, Pad_Height[0] / 2., Pad_Width / 2., Pad_Thick / 2.);
    TGeoVolume * PadInV = new TGeoVolume(PadIn_name, PadIn);
    PadInV->SetMedium(pMedCopper);
    PadInV->SetLineColor(TColor::GetColor("#ffb841"));
    
    TGeoBBox * PadOut = new TGeoBBox(PadOut_name, Pad_Height[1] / 2., Pad_Width / 2., Pad_Thick / 2.);
    TGeoVolume * PadOutV = new TGeoVolume(PadOut_name, PadOut);
    PadOutV->SetMedium(pMedCopper);
    PadOutV->SetLineColor(TColor::GetColor("#c5e384"));
    
    PlacePadsSub(mother_volume, PadInV, R_seg_in, NPads_in_Row_s, Pad_Height[0]);
    PlacePadsSub(mother_volume, PadOutV, R_seg_in + NPads_in_Row_s.size() * (Pad_Height[0] + Pad_Space_h), NPads_in_Row_l, Pad_Height[1]);
}

void PlacePadsSub(TGeoVolume *mother_volume, TGeoVolume * pad_vol, Double_t R_in, const std::vector<Int_t> & nPadsInRow, Double_t padHeight)
{
    Double_t z_offset = (-Pad_Plane_Pp + Pad_Thick) / 2;
    
    for (Int_t row = 0; row < nPadsInRow.size(); ++row)
    {
        Int_t pad_idx = 0;
        Double_t pad_pos = (Pad_Width + Pad_Space_w) / 2.;
        
        for (pad_idx = 0; pad_idx < nPadsInRow[row]; ++pad_idx)
        {
// 2 ways to indexing pads, 0 - default
#if 0
            mother_volume->AddNode(pad_vol, row * 100 + pad_idx, new TGeoTranslation(R_in + (padHeight + Pad_Space_h) * row + padHeight / 2. + Pad_Space_edge, pad_pos, z_offset));
            mother_volume->AddNode(pad_vol, 1e4 + row * 100 + pad_idx, new TGeoTranslation(R_in + (padHeight + Pad_Space_h) * row + padHeight / 2. + Pad_Space_edge, -pad_pos, z_offset));
#else
            mother_volume->AddNode(pad_vol, row * 1e3 + nPadsInRow[row] + pad_idx, new TGeoTranslation(R_in + (padHeight + Pad_Space_h) * row + padHeight / 2. + Pad_Space_edge, pad_pos, z_offset));
            mother_volume->AddNode(pad_vol, row * 1e3 + nPadsInRow[row] - 1 - pad_idx, new TGeoTranslation(R_in + (padHeight + Pad_Space_h) * row + padHeight / 2. + Pad_Space_edge, -pad_pos, z_offset));
#endif
            
            pad_pos += Pad_Width + Pad_Space_w;
        }
        //cout << pad_idx * 2 << " pads in " << row + 1 << " row"<< endl;
    }
}

void FillFramesPoints(unique_ptr<Double_t[]> &x, unique_ptr<Double_t[]> &y)
{
    Double_t x_offset = (End_Frame_Spoke_width / 2 + ROC_Frame_Gap_to_End_Frame_Side) / TMath::Sin(Section_rad_step / 2);
    Double_t R_seg_in = Sensitive_in_R_seg - x_offset;
    Double_t R_seg_out = Sensitive_out_R_seg - x_offset;
    Double_t x_shift = ROC_Frame_width / TMath::Sin(Section_rad_step / 2);
    
    x[0] = R_seg_out;
    x[1] = x[0];
    x[2] = R_seg_in;
    x[3] = x[2];
    x[4] = x[3] + ROC_Frame_width;
    x[5] = x[4];
    x[6] = x[5] + ROC_Frame_lower_chamber_height;
    x[7] = x[6];
    x[8] = x[7] + ROC_Frame_center_width;
    x[9] = x[8];
    x[10] = R_seg_out - ROC_Frame_width;
    x[11] = x[10];
    
    y[0] = 0.;
    y[1] = x[1] * TMath::Tan(Section_rad_step / 2);
    y[2] = x[2] * TMath::Tan(Section_rad_step / 2);
    y[3] = 0.;
    y[4] = 0.;
    y[5] = (x[5] - x_shift) * TMath::Tan(Section_rad_step / 2);
    y[6] = (x[6] - x_shift) * TMath::Tan(Section_rad_step / 2);
    y[7] = 0.;
    y[8] = 0.;
    y[9] = (x[9] - x_shift) * TMath::Tan(Section_rad_step / 2);
    y[10] = (x[10] - x_shift) * TMath::Tan(Section_rad_step / 2);
    y[11] = 0.;
    
    for (Int_t i = 0; i < 12; i++)
        x[i] += x_offset;
}

void PlaceElectronicsComponents(TGeoVolume * mother_vol, TGeoVolume * vol_to_place, Double_t volWidth, Double_t z_offset, Bool_t ElModule)
{
    const Int_t npoints_structure = 12;

    unique_ptr<Double_t[]> x(new Double_t[npoints_structure]);
    unique_ptr<Double_t[]> y(new Double_t[npoints_structure]);

    FillFramesPoints(x, y);
    
    Int_t nxlay = 4;
    Int_t nylay[] = {2, 3, 4, 6};
    Double_t half_low_y = TMath::Abs(y[5] - y[4] - ElModule_CoolPlate_EdgeOffsetFrame);
    Double_t low_x = TMath::Abs(x[7] - x[4] - 2 * ElModule_CoolPlate_EdgeOffsetFrame);
    
    PlaceElectronicsComponentsSub(mother_vol, vol_to_place, nxlay, nylay, x[4] + ElModule_CoolPlate_EdgeOffsetFrame, y[4], low_x, half_low_y, volWidth, z_offset, 0, ElModule);
    
    nxlay = 3;
    Int_t nylay1[] = {5, 7, 7};
    low_x = TMath::Abs(x[11] - x[8] - 2 * ElModule_CoolPlate_EdgeOffsetFrame);
    half_low_y = TMath::Abs(y[9] - y[8] - ElModule_CoolPlate_EdgeOffsetFrame);
    
    PlaceElectronicsComponentsSub(mother_vol, vol_to_place, nxlay, nylay1, x[8] + ElModule_CoolPlate_EdgeOffsetFrame, y[8], low_x, half_low_y, volWidth, z_offset, 400, false);
}

void PlaceElectronicsComponentsSub(TGeoVolume * mother_vol, TGeoVolume * vol_to_place, Int_t nxlay, Int_t nylay[], 
    Double_t x0, Double_t y0, Double_t low_x, Double_t half_low_y, Double_t volWidth, Double_t z_offset, Int_t index_offset, Bool_t ElModule)
{
    Double_t dx = low_x / nxlay;
    Double_t dy = (half_low_y - (ELModule_Width / TMath::Cos(Section_rad_step / 2.0)) * 2) / (nylay[0] / 2.);
    
    Double_t ioffset = nylay[0] / 2.0;
    for (Int_t i = 0; i < 2; i++)
    {
        Double_t xpos = x0 + dx / 2.0;
        Double_t ypos = y0 + (ELModule_Width / TMath::Cos(Section_rad_step / 2.0)) / 2.0 + i * ELModule_Width / TMath::Cos(Section_rad_step / 2.0) + ioffset * dy;
        TGeoCombiTrans * move = new TGeoCombiTrans();
        move->RotateZ(Section_step / 2.0);
        move->SetTranslation(xpos, ypos + dx / 2.0 * TMath::Tan(Section_rad_step / 2.0), z_offset);
        mother_vol->AddNode(vol_to_place, i + 1 + nylay[0] + index_offset, move);
        TGeoCombiTrans * moveRef = new TGeoCombiTrans(*move);
        moveRef->ReflectY(true);
        mother_vol->AddNode(vol_to_place, i + 1 + nylay[0] + 2 + index_offset, moveRef);
        for (Int_t l = 1; l <= nxlay - 1; l++)
        {
            TGeoCombiTrans * movec = new TGeoCombiTrans(*move);
            movec->SetTranslation(xpos + dx * l, ypos + (dx / 2.0 + dx * l) * TMath::Tan(Section_step / 2.0 * TMath::DegToRad()), z_offset);
            mother_vol->AddNode(vol_to_place, l * 100 + i + 1 + nylay[l] + index_offset, movec);
            TGeoCombiTrans * movecr = new TGeoCombiTrans(*movec);
            movecr->ReflectY(true);
            mother_vol->AddNode(vol_to_place, l * 100 + i + 1 + nylay[l] + 2 + index_offset, movecr);
        }
    }
    
    if (ElModule) 
    {
        nylay[1] -= 1;
    }
    for (Int_t j = 0; j < nxlay; j++)
    {
        Double_t xpos = x0 + dx / 2.0 + dx * j;
        for (Int_t i = 0; i < nylay[j] / 2; i++)
        {
            Double_t ymax = half_low_y + (dx * j) * TMath::Tan(Section_step / 2.0 * TMath::DegToRad()) - (ELModule_Width / TMath::Cos(Section_rad_step / 2.0)) * 2;
            Double_t half_vol_to_place = nylay[j] / 2.0;
            Double_t curdy = ymax / half_vol_to_place;
            Double_t ypos;
            if(nylay[j] % 2 == 0)
                ypos = y0 + curdy / 2.0 + i * curdy;
            else
                ypos = y0 + (i + 1) * curdy;
            if (ypos < volWidth / 2.) ypos = volWidth / 2.;
            
            Double_t z_offset_tmp = z_offset;
            if (ElModule && (j == 0 || j == 1))
            {
                z_offset_tmp += ElModule_Cool_Chip0_Th + ELModule_Chip0_Th + ELModule_PCBBox_Th + ELModule_Chip1_Th + ElModule_Cool_Chip1_Th + ElModule_Cool_Th + ElModule_CoolPlateAdditional_Gap;
                
                Double_t coolPlateWidth = 0.;
                if (j == 0)
                    {coolPlateWidth = ELModule_Width;} // 2 modules
                if (j == 1)
                    {coolPlateWidth = ELModule_Width / 2.;} // 1 module
                
                TString tpcElModuleCool_name = "tpc01eModuleCoolPlateAdditional";
                TGeoBBox * tpcElModuleCoolS = new TGeoBBox(tpcElModuleCool_name, ELModule_Height / 2., coolPlateWidth, ElModule_CoolPlateAdditional_Thick / 2.);
                TGeoVolume * tpcElModuleCoolV = new TGeoVolume(tpcElModuleCool_name, tpcElModuleCoolS);
                tpcElModuleCoolV->SetMedium(pMedAluminium);
                tpcElModuleCoolV->SetLineColor(TColor::GetColor("#c4c3be"));
                mother_vol->AddNode(tpcElModuleCoolV, j, new TGeoTranslation(xpos, y0, (Pad_Plane_Pp + Pad_Plane_G10 + ROC_Frame) / 2.  - ElModule_CoolPlate_ZOffsetFrameOutEdge + z_offset_tmp + ElModule_CoolPlateAdditional_Thick / 2.));
                
                z_offset_tmp += ElModule_CoolPlateAdditional_Thick;
                
                if (j == 1)
                {
                    mother_vol->AddNode(vol_to_place, j * 100 + index_offset, new TGeoTranslation(xpos, y0, z_offset_tmp));
                    z_offset_tmp = z_offset;
                    ypos += Connector_hole_width / 2.; //fixme
                }
            }
            
            TGeoTranslation * move = new TGeoTranslation(xpos, ypos, z_offset_tmp);
            mother_vol->AddNode(vol_to_place, i + 1 + j * 100 + index_offset, move);
            TGeoCombiTrans * moveRef = new TGeoCombiTrans(*move);
            moveRef->ReflectY(true);
            mother_vol->AddNode(vol_to_place, i + 1 + nylay[j] / 2 + j * 100 + index_offset, moveRef);
        }
        if(nylay[j] % 2 != 0)
            mother_vol->AddNode(vol_to_place, j * 100 + index_offset, new TGeoTranslation(xpos, y0, z_offset));
    }
}

void CreateElectronicsModules(TGeoVolume *mother_volume)
{
    const Int_t nPoints_fr = 12;

    unique_ptr<Double_t[]> x_fr(new Double_t[nPoints_fr]);
    unique_ptr<Double_t[]> y_fr(new Double_t[nPoints_fr]);

    FillFramesPoints(x_fr, y_fr);
    
    TString tpcElModuleCoolPlateIn_name = "tpc01eModuleCoolPlateIn";
    TString tpcElModuleCoolPlateOut_name = "tpc01eModuleCoolPlateOut";
    
    const Int_t nPoints_coolPlane = 4;
    Double_t z0 = (Pad_Plane_Pp + Pad_Plane_G10 + ROC_Frame) / 2.;
    
    unique_ptr<Double_t[]> x_in(new Double_t[nPoints_coolPlane]);
    unique_ptr<Double_t[]> y_in(new Double_t[nPoints_coolPlane]);
    x_in[0] = x_in[1] = x_fr[6] - ElModule_CoolPlate_EdgeOffsetFrame;
    x_in[2] = x_in[3] = x_fr[5] + ElModule_CoolPlate_EdgeOffsetFrame;
    y_in[1] = y_fr[6] - ElModule_CoolPlate_EdgeOffsetFrame / TMath::Tan((TMath::PiOver2() - Section_rad_step / 2.) / 2.);
    y_in[0] = -y_in[1];
    y_in[2] = y_fr[5] - ElModule_CoolPlate_EdgeOffsetFrame / TMath::Tan((TMath::PiOver2() + Section_rad_step / 2.) / 2.);
    y_in[3] = -y_in[2];
    TGeoXtru * CoolPlateInS = new TGeoXtru(2);
    CoolPlateInS->DefinePolygon(nPoints_coolPlane, x_in.get(), y_in.get());
    CoolPlateInS->DefineSection(0, z0 - ElModule_CoolPlate_ZOffsetFrameOutEdge - ElModule_CoolPlate_Thick);
    CoolPlateInS->DefineSection(1, z0 - ElModule_CoolPlate_ZOffsetFrameOutEdge);
    
    TGeoVolume *CoolPlateInV = new TGeoVolume(tpcElModuleCoolPlateIn_name, CoolPlateInS);
    CoolPlateInV->SetMedium(pMedAluminium);
    CoolPlateInV->SetLineColor(TColor::GetColor("#c4c3be"));
    mother_volume->AddNode(CoolPlateInV, 0);
    
    unique_ptr<Double_t[]> x_out(new Double_t[nPoints_coolPlane]);
    unique_ptr<Double_t[]> y_out(new Double_t[nPoints_coolPlane]);
    x_out[0] = x_out[1] = x_fr[10] - ElModule_CoolPlate_EdgeOffsetFrame;
    x_out[2] = x_out[3] = x_fr[9] + ElModule_CoolPlate_EdgeOffsetFrame;
    y_out[1] = y_fr[10] - ElModule_CoolPlate_EdgeOffsetFrame / TMath::Tan((TMath::PiOver2() - Section_rad_step / 2.) / 2.);
    y_out[0] = -y_out[1];
    y_out[2] = y_fr[9] - ElModule_CoolPlate_EdgeOffsetFrame / TMath::Tan((TMath::PiOver2() + Section_rad_step / 2.) / 2.);
    y_out[3] = -y_out[2];
    TGeoXtru * CoolPlateOutS = new TGeoXtru(2);
    CoolPlateOutS->DefinePolygon(nPoints_coolPlane, x_out.get(), y_out.get());
    CoolPlateOutS->DefineSection(0, z0 - ElModule_CoolPlate_ZOffsetFrameOutEdge - ElModule_CoolPlate_Thick);
    CoolPlateOutS->DefineSection(1, z0 - ElModule_CoolPlate_ZOffsetFrameOutEdge);
    
    TGeoVolume *CoolPlateOutV = new TGeoVolume(tpcElModuleCoolPlateOut_name, CoolPlateOutS);
    CoolPlateOutV->SetMedium(pMedAluminium);
    CoolPlateOutV->SetLineColor(TColor::GetColor("#c4c3be"));
    mother_volume->AddNode(CoolPlateOutV, 0);
    
    TString tpcElModuleCoolPlateHole_name = "tpc01eModuleCoolPlateHole";
    
    TGeoVolume * tpcElModuleCoolPlateHole = new TGeoVolumeAssembly(tpcElModuleCoolPlateHole_name);
    
    TString tpcElModuleCoolPlateHoleB_name = "tpc01eModuleCoolPlateHoleB";
    TGeoBBox * tpcElModuleCoolPlateHoleB = new TGeoBBox(tpcElModuleCoolPlateHoleB_name, ElModule_CoolPlate_Hole_H / 2., ElModule_CoolPlate_Hole_W / 2., ElModule_CoolPlate_Thick / 2.);
    TGeoVolume *tpcElModuleCoolPlateHoleBV = new TGeoVolume(tpcElModuleCoolPlateHoleB_name, tpcElModuleCoolPlateHoleB);
    tpcElModuleCoolPlateHoleBV->SetMedium(pMedAir);
    tpcElModuleCoolPlateHoleBV->SetLineColor(kRed);
    tpcElModuleCoolPlateHole->AddNode(tpcElModuleCoolPlateHoleBV, 0, new TGeoTranslation(0., 0., z0 - ElModule_CoolPlate_ZOffsetFrameOutEdge - ElModule_CoolPlate_Thick / 2.));
    
    TString tpcElModuleCoolPlateHoleE_name = "tpc01eModuleCoolPlateHoleE";
    TGeoTubeSeg * tpcElModuleCoolPlateHoleE = new TGeoTubeSeg(tpcElModuleCoolPlateHoleE_name, 0.0, ElModule_CoolPlate_Hole_W / 2, ElModule_CoolPlate_Thick / 2, -90.0, 90.0);
    TGeoVolume *tpcElModuleCoolPlateHoleEV = new TGeoVolume(tpcElModuleCoolPlateHoleE_name, tpcElModuleCoolPlateHoleE);
    tpcElModuleCoolPlateHoleEV->SetMedium(pMedAir);
    tpcElModuleCoolPlateHoleEV->SetLineColor(kRed);
    TGeoCombiTrans * coolPlateHolePos0 = new TGeoCombiTrans(ElModule_CoolPlate_Hole_H / 2, 0., z0 - ElModule_CoolPlate_ZOffsetFrameOutEdge - ElModule_CoolPlate_Thick / 2., NULL);
    tpcElModuleCoolPlateHole->AddNode(tpcElModuleCoolPlateHoleEV, 0, coolPlateHolePos0);
    TGeoCombiTrans * coolPlateHolePos1 = new TGeoCombiTrans(ElModule_CoolPlate_Hole_H / 2, 0., z0 - ElModule_CoolPlate_ZOffsetFrameOutEdge - ElModule_CoolPlate_Thick / 2., NULL);
    coolPlateHolePos1->ReflectX(true);
    tpcElModuleCoolPlateHole->AddNode(tpcElModuleCoolPlateHoleEV, 1, coolPlateHolePos1);
    
    tpcElModuleCoolPlateHole->SetMedium(pMedAluminium);
    
    //CoolPlateOutV->AddNode(tpcElModuleCoolPlateHole, 0, new TGeoTranslation(x_out[2] + 10., 0.0, 0.0));
    Int_t nxlay = 4;
    Int_t nylay[] = {2, 3, 4, 6};
    Double_t half_low_y = TMath::Abs(y_in[2] - y_in[3]) / 2.;
    Double_t low_x = TMath::Abs(x_in[0] - x_in[2]);
    PlaceElectronicsComponentsSub(CoolPlateInV, tpcElModuleCoolPlateHole, nxlay, nylay, x_in[2], 0., low_x, half_low_y, ElModule_CoolPlate_Hole_W, 0., 0, false);
    
    nxlay = 3;
    Int_t nylay1[] = {5, 7, 7};
    low_x = TMath::Abs(x_out[0] - x_out[2]);
    half_low_y = TMath::Abs(y_out[2] - y_out[3]) / 2.;
    PlaceElectronicsComponentsSub(CoolPlateOutV, tpcElModuleCoolPlateHole, nxlay, nylay1, x_out[2], 0., low_x, half_low_y, ElModule_CoolPlate_Hole_W, 0., 400, false);
    
    
    TString tpcElModule_name = "tpc01eModule";
    
    TGeoVolume * tpcElModule = new TGeoVolumeAssembly(tpcElModule_name);
        
    const Double_t z0_Module = z0 - ElModule_CoolPlate_ZOffsetFrameOutEdge;
    
    TString tpcElModuleCoolChip0_name = "tpc01eModuleCoolChip0";
    TGeoBBox * tpcElModCoolChip0S = new TGeoBBox(tpcElModuleCoolChip0_name, ElModule_Cool_Chip0_WH / 2., ElModule_Cool_Chip0_WH / 2., ElModule_Cool_Chip0_Th / 2);
    TGeoVolume * tpcElModCoolChip0V = new TGeoVolume(tpcElModuleCoolChip0_name, tpcElModCoolChip0S);
    tpcElModCoolChip0V->SetMedium(pMedAluminium);
    tpcElModCoolChip0V->SetLineColor(TColor::GetColor("#c4c3be"));
    tpcElModule->AddNode(tpcElModCoolChip0V, 0, new TGeoTranslation(ElModule_Cool_Chip0_HOffsetUp, ElModule_Cool_Chip0_WOffset, z0_Module + ElModule_Cool_Chip0_Th / 2.));
    tpcElModule->AddNode(tpcElModCoolChip0V, 1, new TGeoTranslation(ElModule_Cool_Chip0_HOffsetDown, ElModule_Cool_Chip0_WOffset, z0_Module + ElModule_Cool_Chip0_Th / 2.));
    
    TString tpcElModulePCB0_name = "tpc01eModulePCB0";
    TGeoBBox * tpcElModPCB0S = new TGeoBBox(tpcElModulePCB0_name, ELModule_Height / 2., ELModule_Width / 2., ELModule_PCB_Th / 2.);
    TGeoVolume * tpcElModPCB0V = new TGeoVolume(tpcElModulePCB0_name, tpcElModPCB0S);
    tpcElModPCB0V->SetMedium(pMedFiberGlass);
    tpcElModPCB0V->SetLineColor(kGreen + 2);
    tpcElModule->AddNode(tpcElModPCB0V, 0, new TGeoTranslation(0., 0., z0_Module + ElModule_Cool_Chip0_Th + ELModule_Chip0_Th + ELModule_PCB_Th / 2));
    
    TString tpcElModulePCBCu_name = "tpc01eModulePCBCu";
    TGeoBBox * tpcElModPCBCuS = new TGeoBBox(tpcElModulePCBCu_name, ELModule_Height / 2., ELModule_Width / 2., ELModule_PCBCu_Th / 2.);
    TGeoVolume * tpcElModPCBCuV = new TGeoVolume(tpcElModulePCBCu_name, tpcElModPCBCuS);
    tpcElModPCBCuV->SetMedium(pMedCopper);
    tpcElModPCBCuV->SetLineColor(TColor::GetColor("#ffb841"));
    tpcElModule->AddNode(tpcElModPCBCuV, 0, new TGeoTranslation(0., 0., z0_Module + ElModule_Cool_Chip0_Th + ELModule_Chip0_Th + ELModule_PCB_Th + ELModule_PCBCu_Th / 2));
    
    TString tpcElModulePCB1_name = "tpc01eModulePCB1";
    TGeoBBox * tpcElModPCB1S = new TGeoBBox(tpcElModulePCB1_name, ELModule_Height / 2., ELModule_Width / 2., ELModule_PCB_Th / 2.);
    TGeoVolume * tpcElModPCB1V = new TGeoVolume(tpcElModulePCB1_name, tpcElModPCB1S);
    tpcElModPCB1V->SetMedium(pMedFiberGlass);
    tpcElModPCB1V->SetLineColor(kGreen + 2);
    tpcElModule->AddNode(tpcElModPCB1V, 0, new TGeoTranslation(0., 0., z0_Module + ElModule_Cool_Chip0_Th + ELModule_Chip0_Th + ELModule_PCBBox_Th - ELModule_PCB_Th / 2));
    
    TString tpcElModulePCBConnector_name = "tpc01eModulePCBConnector";
    TGeoBBox * tpcElModPCBConS = new TGeoBBox(tpcElModulePCBConnector_name, ELModule_PCB0_Connector_H / 2., ELModule_PCB0_Connector_W / 2., ELModule_PCB0_Connector_Th / 2.);
    TGeoVolume * tpcElModPCBConV = new TGeoVolume(tpcElModulePCBConnector_name, tpcElModPCBConS);
    tpcElModPCBConV->SetMedium(pMedPlastic);
    tpcElModPCBConV->SetLineColor(kWhite);
    tpcElModule->AddNode(tpcElModPCBConV, 0, new TGeoTranslation(0., ElModule_PCB0_Connector_WOffset, z0_Module + ElModule_Cool_Chip0_Th + ELModule_Chip0_Th - ELModule_PCB0_Connector_Th / 2));
    
    TString tpcElModuleCoolChip1_name = "tpc01eModuleCoolChip1";
    TGeoBBox * tpcElModCoolChip1S = new TGeoBBox(tpcElModuleCoolChip1_name, ElModule_Cool_Chip1_H / 2., ElModule_Cool_Chip1_W / 2., ElModule_Cool_Chip1_Th / 2);
    TGeoVolume * tpcElModCoolChip1V = new TGeoVolume(tpcElModuleCoolChip1_name, tpcElModCoolChip1S);
    tpcElModCoolChip1V->SetMedium(pMedAluminium);
    tpcElModCoolChip1V->SetLineColor(TColor::GetColor("#c4c3be"));
    tpcElModule->AddNode(tpcElModCoolChip1V, 0, new TGeoTranslation(ElModule_Cool_Chip1_HOffset, ElModule_Cool_Chip1_WOffset, z0_Module + ElModule_Cool_Chip0_Th + ELModule_Chip0_Th + ELModule_PCBBox_Th + ELModule_Chip1_Th + ElModule_Cool_Chip1_Th / 2.));
    
    TString tpcElModuleCool_name = "tpc01eModuleCool";
    TGeoBBox * tpcElModuleCoolS = new TGeoBBox(tpcElModuleCool_name, ELModule_Height / 2., ELModule_Width / 2., ElModule_Cool_Th / 2.);
    TGeoVolume * tpcElModuleCoolV = new TGeoVolume(tpcElModuleCool_name, tpcElModuleCoolS);
    tpcElModuleCoolV->SetMedium(pMedAluminium);
    tpcElModuleCoolV->SetLineColor(TColor::GetColor("#c4c3be"));
    tpcElModule->AddNode(tpcElModuleCoolV, 0, new TGeoTranslation(0., 0., z0_Module + ElModule_Cool_Chip0_Th + ELModule_Chip0_Th + ELModule_PCBBox_Th + ELModule_Chip1_Th + ElModule_Cool_Chip1_Th + ElModule_Cool_Th / 2.));
    
    TString tpcElModuleCoolHole_name = "tpc01eModuleCoolHole";
    TGeoBBox * tpcElModuleCoolHoleS = new TGeoBBox(tpcElModuleCoolHole_name, ElModule_Cool_Hole_H / 2., ElModule_Cool_Hole_W / 2., ElModule_Cool_Th / 2.);
    TGeoVolume * tpcElModuleCoolHoleV = new TGeoVolume(tpcElModuleCoolHole_name, tpcElModuleCoolHoleS);
    tpcElModuleCoolHoleV->SetMedium(pMedAir);
    tpcElModuleCoolHoleV->SetLineColor(kRed);
    tpcElModuleCoolV->AddNode(tpcElModuleCoolHoleV, 0, new TGeoTranslation(ElModule_Cool_Hole_HOffset, ElModule_Cool_Hole_WOffset, 0.));
    
    PlaceElectronicsComponents(mother_volume, tpcElModule, ELModule_Width, 0., true);
}

void CreateEndCaps(TGeoVolume *mother_volume) {
    TString endcap_name = "tpc01ec";

    Double_t endcap_thickness = Container_tpc_z/2 - Chamber_tpc_z/2 - Chamber_tpc_z_ext_parts - ROC_Frame_zspace_iocards;

    TGeoTube *tpcEndCapS = new TGeoTube(endcap_name, TpcInnerRadius, TpcOuterRadius, endcap_thickness/2);
    TGeoVolume *tpcEndCapV = new TGeoVolume(endcap_name, tpcEndCapS);
    tpcEndCapV->SetMedium(pMedAir);
    tpcEndCapV->SetLineColor(TColor::GetColor("#d0e3f7"));
    
    TGeoCombiTrans *endcap_position1 = new TGeoCombiTrans(0., 0., Container_tpc_z/2 - endcap_thickness/2, NULL);
    TGeoCombiTrans *endcap_position2 = new TGeoCombiTrans(*endcap_position1);
    endcap_position2->ReflectZ(true);
    
    mother_volume->AddNode(tpcEndCapV, 0, endcap_position1);
    mother_volume->AddNode(tpcEndCapV, 1, endcap_position2);
    
    //create ribs
    CreateRibs(tpcEndCapV, endcap_thickness);
    
    //cooling
    CreateThermoscreenEnds(tpcEndCapV, endcap_thickness);
}

void MakeRibProfilePoints(unique_ptr<Double_t[]> &x, unique_ptr<Double_t[]> &y)
{
    x[0] = x[7] = 0.0;
    x[1] = x[2] = x[5] = x[6] = RibProfileHeight;
    x[3] = x[4] = RibProfileThickness;
    
    y[0] = y[1] = -RibProfileWidth / 2.0;
    y[2] = y[3] = -RibProfileWidth / 2.0 + RibProfileThickness;
    y[4] = y[5] = RibProfileWidth / 2.0 - RibProfileThickness;
    y[6] = y[7] = RibProfileWidth / 2.0;
}

void MakeRibProfileTriConnectorsPoints(unique_ptr<Double_t[]> &x, unique_ptr<Double_t[]> &y)
{
    x[0] = x[2] = 0.0;
    x[1] = 2 * RibProfileHeight * TMath::Tan(Section_step / 2.0 * TMath::DegToRad());
    
    y[0] = y[1] = 0.0;
    y[2] = 2 * RibProfileHeight;
}

void CreateRibs(TGeoVolume *mother_volume, Double_t endcap_zwidth) {
    TString rib_name = "tpc01Rib";
    
    Double_t RibsSectorXlow = RibsSectorRMin * TMath::Tan(Section_step / 2.0 * TMath::DegToRad());
    Double_t RibsSectorXhigh = RibsSectorRMax * TMath::Tan(Section_step / 2.0 * TMath::DegToRad());
    TGeoTrd1 *ribPlaneS = new TGeoTrd1(rib_name, RibsSectorXlow, RibsSectorXhigh, RibProfileWidth / 2, (RibsSectorRMax - RibsSectorRMin) / 2);
    
    TGeoVolume *ribPlaneV = new TGeoVolume(rib_name, ribPlaneS);
    ribPlaneV->SetMedium(pMedAir);

    ribPlaneV->SetLineColor(kWhite);

    for(Int_t isec = 0; isec < Nsections; isec++) {
        TGeoCombiTrans *combi_trans = new TGeoCombiTrans();
        combi_trans->RotateX(-90.0);
        combi_trans->SetTranslation(0, RibsSectorRMin + (RibsSectorRMax - RibsSectorRMin) / 2, -endcap_zwidth / 2 + RibsZfromTPCends);
        combi_trans->RotateZ(isec * Section_step);

        mother_volume->AddNode(ribPlaneV, isec, combi_trans);
    }
    
    const Int_t npoints_RibStructure = 8;
    unique_ptr<Double_t[]> xRib(new Double_t[npoints_RibStructure]);
    unique_ptr<Double_t[]> yRib(new Double_t[npoints_RibStructure]);
    MakeRibProfilePoints(xRib, yRib);
    TGeoMedium *ribFrame_medium = pMedAluminium; //set medium
    
    TString rib_fr_rad = "tpc01RibRad";

    TGeoXtru *ribFrameRadialS = new TGeoXtru(2);
    ribFrameRadialS->DefinePolygon(npoints_RibStructure, xRib.get(), yRib.get());
    Double_t length = (RibsSectorRMax - RibsSectorRMin) / TMath::Cos(Section_step / 2.0 * TMath::DegToRad()) - RibProfileHeight * TMath::Tan(Section_step / 2.0 * TMath::DegToRad());
    ribFrameRadialS->DefineSection(0, -length/2 - RibProfileHeight * TMath::Tan(Section_step / 2.0 * TMath::DegToRad()) / 2);
    ribFrameRadialS->DefineSection(1, length/2 - RibProfileHeight * TMath::Tan(Section_step / 2.0 * TMath::DegToRad()) / 2);

    TGeoVolume *ribFrameRadialV = new TGeoVolume(rib_fr_rad, ribFrameRadialS);
    ribFrameRadialV->SetMedium(ribFrame_medium);
    ribFrameRadialV->SetLineColor(TColor::GetColor("#84c3be"));

    TGeoCombiTrans *combi_trans1 = new TGeoCombiTrans();
    combi_trans1->RotateY(-Section_step / 2.0);
    combi_trans1->SetTranslation(-RibsSectorXlow - (RibsSectorXhigh - RibsSectorXlow) / 2, 0, 0);
    TGeoCombiTrans *combi_trans2 = new TGeoCombiTrans(*combi_trans1);
    combi_trans2->ReflectX(true);
    
    ribPlaneV->AddNode(ribFrameRadialV, 0, combi_trans1);
    ribPlaneV->AddNode(ribFrameRadialV, 1, combi_trans2);
    
    TString rib_fr_horl = "tpc01RibHorL";
    
    TGeoXtru *ribFrameHorLS = new TGeoXtru(2);
    ribFrameHorLS->DefinePolygon(npoints_RibStructure, xRib.get(), yRib.get());
    Double_t lengthHL = RibsSectorXlow - RibProfileHeight * TMath::Cos(Section_step / 2.0 * TMath::DegToRad());
    ribFrameHorLS->DefineSection(0, -lengthHL);
    ribFrameHorLS->DefineSection(1, lengthHL);
    
    TGeoVolume *ribFrameHorV = new TGeoVolume(rib_fr_horl, ribFrameHorLS);
    ribFrameHorV->SetMedium(ribFrame_medium);
    ribFrameHorV->SetLineColor(TColor::GetColor("#84c3be"));

    TGeoCombiTrans *combi_trans3 = new TGeoCombiTrans();
    combi_trans3->RotateY(90);
    combi_trans3->SetTranslation(0, 0, -(RibsSectorRMax - RibsSectorRMin) / 2 + RibProfileHeight * (1 + TMath::Sin(Section_step / 2.0 * TMath::DegToRad())));
    TGeoCombiTrans *combi_trans4 = new TGeoCombiTrans();
    combi_trans4->RotateY(-90);
    combi_trans4->SetTranslation(0, 0, -(RibsSectorRMax - RibsSectorRMin) / 2 + RibProfileHeight * (1 + TMath::Sin(Section_step / 2.0 * TMath::DegToRad())));
    
    ribPlaneV->AddNode(ribFrameHorV, 0, combi_trans3);
    ribPlaneV->AddNode(ribFrameHorV, 1, combi_trans4);
    
    TString rib_fr_horh = "tpc01RibHorH";
    
    TGeoXtru *ribFrameHorHS = new TGeoXtru(2);
    ribFrameHorHS->DefinePolygon(npoints_RibStructure, xRib.get(), yRib.get());
    Double_t lengthHH = RibsSectorXhigh - RibProfileHeight / TMath::Cos(Section_step / 2.0 * TMath::DegToRad()) - 2 * RibProfileHeight * TMath::Tan(Section_step / 2.0 * TMath::DegToRad());
    ribFrameHorHS->DefineSection(0, -lengthHH);
    ribFrameHorHS->DefineSection(1, lengthHH);
    
    TGeoVolume *ribFrameHorHV = new TGeoVolume(rib_fr_horh, ribFrameHorHS);
    ribFrameHorHV->SetMedium(ribFrame_medium);
    ribFrameHorHV->SetLineColor(TColor::GetColor("#84c3be"));

    TGeoCombiTrans *combi_trans5 = new TGeoCombiTrans();
    combi_trans5->RotateY(90);
    combi_trans5->SetTranslation(0, 0, (RibsSectorRMax - RibsSectorRMin) / 2 - RibProfileHeight);
    TGeoCombiTrans *combi_trans6 = new TGeoCombiTrans();
    combi_trans6->RotateY(-90);
    combi_trans6->SetTranslation(0, 0, (RibsSectorRMax - RibsSectorRMin) / 2 - RibProfileHeight);
    
    ribPlaneV->AddNode(ribFrameHorHV, 0, combi_trans5);
    ribPlaneV->AddNode(ribFrameHorHV, 1, combi_trans6);
    
    TString rib_fr_triconn = "tpc01RibHorH";
    
    const Int_t npoints_RibTriConnectors = 3;
    unique_ptr<Double_t[]> xRibTriCon(new Double_t[npoints_RibTriConnectors]);
    unique_ptr<Double_t[]> yRibTriCon(new Double_t[npoints_RibTriConnectors]);
    MakeRibProfileTriConnectorsPoints(xRibTriCon, yRibTriCon);
    
    TGeoXtru *ribTriConnectorsS = new TGeoXtru(2);
    ribTriConnectorsS->DefinePolygon(npoints_RibTriConnectors, xRibTriCon.get(), yRibTriCon.get());
    ribTriConnectorsS->DefineSection(0, -RibProfileThickness / 2);
    ribTriConnectorsS->DefineSection(1, RibProfileThickness / 2);
    
    TGeoVolume *ribTriConnectorsV = new TGeoVolume(rib_fr_triconn, ribTriConnectorsS);
    ribTriConnectorsV->SetMedium(ribFrame_medium);
    ribTriConnectorsV->SetLineColor(TColor::GetColor("#84c3be"));
    
    TGeoCombiTrans *combi_trans7 = new TGeoCombiTrans();
    combi_trans7->RotateX(-90);
    combi_trans7->SetTranslation(lengthHL, -RibProfileWidth / 2 + RibProfileThickness / 2, -(RibsSectorRMax - RibsSectorRMin) / 2 + RibProfileHeight * (2 + TMath::Sin(Section_rad_step / 2.0)));
    ribPlaneV->AddNode(ribTriConnectorsV, 0, combi_trans7);
    
    TGeoCombiTrans *combi_trans8 = new TGeoCombiTrans(*combi_trans7);
    combi_trans8->ReflectX(true);
    ribPlaneV->AddNode(ribTriConnectorsV, 1, combi_trans8);
    
    TGeoCombiTrans *combi_trans9 = new TGeoCombiTrans(*combi_trans7);
    combi_trans9->ReflectY(true);
    ribPlaneV->AddNode(ribTriConnectorsV, 2, combi_trans9);
    
    TGeoCombiTrans *combi_trans10 = new TGeoCombiTrans(*combi_trans9);
    combi_trans10->ReflectX(true);
    ribPlaneV->AddNode(ribTriConnectorsV, 3, combi_trans10);
    
    TGeoCombiTrans *combi_trans11 = new TGeoCombiTrans();
    combi_trans11->RotateX(-90);
    combi_trans11->SetTranslation(lengthHH, -RibProfileWidth / 2 + RibProfileThickness / 2, (RibsSectorRMax - RibsSectorRMin) / 2);
    ribPlaneV->AddNode(ribTriConnectorsV, 4, combi_trans11);
    
    TGeoCombiTrans *combi_trans12 = new TGeoCombiTrans(*combi_trans11);
    combi_trans12->ReflectX(true);
    ribPlaneV->AddNode(ribTriConnectorsV, 5, combi_trans12);
    
    TGeoCombiTrans *combi_trans13 = new TGeoCombiTrans(*combi_trans11);
    combi_trans13->ReflectY(true);
    ribPlaneV->AddNode(ribTriConnectorsV, 6, combi_trans13);
    
    TGeoCombiTrans *combi_trans14 = new TGeoCombiTrans(*combi_trans13);
    combi_trans14->ReflectX(true);
    ribPlaneV->AddNode(ribTriConnectorsV, 7, combi_trans14);
}

void MakeROCFramesHoles(TGeoVolume * bpG10Vol, TGeoVolume * bpAlVol)
{
    Double_t x_offset = (End_Frame_Spoke_width / 2 + ROC_Frame_Gap_to_End_Frame_Side - ROC_Frame_ext_width) / TMath::Sin(Section_rad_step / 2);
    
    TString tpcElHoleG10_name = "tpc01ehg10";
    TGeoVolumeAssembly *tpcElHoleG10V = new TGeoVolumeAssembly(tpcElHoleG10_name);
    
    TString tpcElHoleBaseG10_name = "tpc01ehb";
    TGeoBBox * tpcElHoleG10B = new TGeoBBox(tpcElHoleBaseG10_name, Connector_hole_length / 2, Connector_hole_width / 2, Pad_Plane_G10 / 2);
    TGeoVolume *tpcElHoleG10BV = new TGeoVolume(tpcElHoleBaseG10_name, tpcElHoleG10B);
    tpcElHoleG10BV->SetMedium(pMedAir);
    tpcElHoleG10BV->SetLineColor(kWhite);
    tpcElHoleG10V->AddNode(tpcElHoleG10BV, 0, new TGeoTranslation(-x_offset, 0.0, 0.0));
    
    TString tpcElConnectorG10_name = "tpc01ec";
    TGeoBBox * tpcElConnectorG10 = new TGeoBBox(tpcElConnectorG10_name, Connector_length / 2, Connector_width / 2, Pad_Plane_G10 / 2);
    TGeoVolume *tpcElConnectorG10V = new TGeoVolume(tpcElConnectorG10_name, tpcElConnectorG10);
    tpcElConnectorG10V->SetMedium(pMedPlastic);
    tpcElConnectorG10V->SetLineColor(kWhite);
    Double_t move = Connector_length / 2 + (Connector_hole_length / 2.0 - Connector_length);
    tpcElHoleG10BV->AddNode(tpcElConnectorG10V, 0, new TGeoTranslation(move, 0.0, 0.0));
    tpcElHoleG10BV->AddNode(tpcElConnectorG10V, 1, new TGeoTranslation(-move, 0.0, 0.0));
    
    TString tpcElHoleEndG10_name = "tpc01ehe";
    TGeoTubeSeg * tpcElHoleG10E = new TGeoTubeSeg(tpcElHoleEndG10_name, 0.0, Connector_hole_width / 2, Pad_Plane_G10 / 2, -90.0, 90.0);
    TGeoVolume *tpcElHoleG10EV = new TGeoVolume(tpcElHoleEndG10_name, tpcElHoleG10E);
    tpcElHoleG10EV->SetMedium(pMedAir);
    tpcElHoleG10EV->SetLineColor(kWhite);
    TGeoCombiTrans * eheG10pos0 = new TGeoCombiTrans(Connector_hole_length / 2 - x_offset, 0.0, 0.0, NULL);
    tpcElHoleG10V->AddNode(tpcElHoleG10EV, 0, eheG10pos0);
    TGeoCombiTrans * eheG10pos1 = new TGeoCombiTrans(Connector_hole_length / 2 + x_offset, 0.0, 0.0, NULL);
    eheG10pos1->ReflectX(true);
    tpcElHoleG10V->AddNode(tpcElHoleG10EV, 1, eheG10pos1);
    tpcElHoleG10V->SetMedium(pMedG10);
    
    PlaceElectronicsComponents(bpG10Vol, tpcElHoleG10V, Connector_hole_width);
    
    TString tpcElHoleAl_name = "tpc01ehal";
    TGeoVolumeAssembly *tpcElHoleV = new TGeoVolumeAssembly(tpcElHoleAl_name);
    
    TString tpcElHoleBaseAl_name = "tpc01ehb";
    TGeoBBox * tpcElHoleB = new TGeoBBox(tpcElHoleBaseAl_name, Connector_hole_length / 2, Connector_hole_width / 2, ROC_Frame_plane / 2);
    TGeoVolume *tpcElHoleBV = new TGeoVolume(tpcElHoleBaseAl_name, tpcElHoleB);
    tpcElHoleBV->SetMedium(pMedAir);
    tpcElHoleBV->SetLineColor(kWhite);
    tpcElHoleV->AddNode(tpcElHoleBV, 0, new TGeoTranslation(-x_offset, 0.0, 0.0));
    
    TString tpcElConnectorAl_name = "tpc01ec";
    TGeoBBox * tpcElConnectorAl = new TGeoBBox(tpcElConnectorAl_name, Connector_length / 2, Connector_width / 2, ROC_Frame_plane / 2);
    TGeoVolume *tpcElConnectorAlV = new TGeoVolume(tpcElConnectorAl_name, tpcElConnectorAl);
    tpcElConnectorAlV->SetMedium(pMedPlastic);
    tpcElConnectorAlV->SetLineColor(kWhite);
    tpcElHoleBV->AddNode(tpcElConnectorAlV, 0, new TGeoTranslation(move, 0.0, 0.0));
    tpcElHoleBV->AddNode(tpcElConnectorAlV, 1, new TGeoTranslation(-move, 0.0, 0.0));
    
    TString tpcElHoleEndAl_name = "tpc01ehe";
    TGeoTubeSeg * tpcElHoleE = new TGeoTubeSeg(tpcElHoleEndAl_name, 0.0, Connector_hole_width / 2, ROC_Frame_plane / 2, -90.0, 90.0);
    TGeoVolume *tpcElHoleEV = new TGeoVolume(tpcElHoleEndAl_name, tpcElHoleE);
    tpcElHoleEV->SetMedium(pMedAir);
    tpcElHoleEV->SetLineColor(kWhite);
    TGeoCombiTrans * eheAlpos0 = new TGeoCombiTrans(Connector_hole_length / 2 - x_offset, 0.0, 0.0, NULL);
    tpcElHoleV->AddNode(tpcElHoleEV, 0, eheAlpos0);
    TGeoCombiTrans * eheAlpos1 = new TGeoCombiTrans(Connector_hole_length / 2 + x_offset, 0.0, 0.0, NULL);
    eheAlpos1->ReflectX(true);
    tpcElHoleV->AddNode(tpcElHoleEV, 1, eheAlpos1);
    tpcElHoleV->SetMedium(pMedAluminium);
    
    PlaceElectronicsComponents(bpAlVol, tpcElHoleV, Connector_hole_width);
}

void CreateROCPressingBackPlane(TGeoVolume *mother_volume, Double_t width, Double_t rmin, Double_t rmax)
{
    TString plane_name = "tpc01ROCPressBPlane";
    
    Int_t npoints_structure = 11;
    unique_ptr<Double_t[]> x(new Double_t[npoints_structure]);
    unique_ptr<Double_t[]> y(new Double_t[npoints_structure]);
    
    Double_t x_shift = ROCPress_Backplate_SideWidth / TMath::Sin(Section_rad_step / 2);
    Double_t angle = TMath::ACos(ROCPress_Backplate_HighExtX / ROCPress_Backplate_HighExtL[1]);
    x[0] = x[1] = rmax - ROCPress_Backplate_HighExtX;
    x[2] = x[3] = rmax;
    x[4] = rmax - ROCPress_Backplate_HighExtL[3] * TMath::Cos(angle);
    x[5] = x[6] = rmin;
    x[7] = x[8] = rmin + ROCPress_Backplate_LowWidth;
    x[9] = x[10] = rmax - ROCPress_Backplate_HighWidth;
    y[0] = y[6] = y[7] = y[10] = 0.;
    y[4] = x[4] * TMath::Tan(Section_rad_step / 2);
    y[5] = x[5] * TMath::Tan(Section_rad_step / 2);
    y[8] = (x[8] - x_shift) * TMath::Tan(Section_rad_step / 2);
    y[9] = (x[9] - x_shift) * TMath::Tan(Section_rad_step / 2);
    y[1] = ROCPress_Backplate_HighExtL[0];
    y[2] = y[1] + ROCPress_Backplate_HighExtL[1] * TMath::Sin(angle);
    y[3] = y[2] + ROCPress_Backplate_HighExtL[2];
    
    TGeoXtru *ROCPressBPlaneS = new TGeoXtru(2);
    ROCPressBPlaneS->DefinePolygon(npoints_structure, x.get(), y.get());
    ROCPressBPlaneS->DefineSection(0, 0.);
    ROCPressBPlaneS->DefineSection(1, ROCPress_Backplate_thick);
    
    TGeoVolume * ROCPressBPlaneV = new TGeoVolume(plane_name, ROCPressBPlaneS);
    ROCPressBPlaneV->SetMedium(pMedMylar);
    ROCPressBPlaneV->SetLineColor(kYellow);
    
    TGeoCombiTrans * pos = new TGeoCombiTrans(0., 0., 0., NULL);
    TGeoCombiTrans * pos_refl = new TGeoCombiTrans(*pos);
    pos_refl->ReflectY(true);
    
    mother_volume->AddNode(ROCPressBPlaneV, 0, pos);
    mother_volume->AddNode(ROCPressBPlaneV, 1, pos_refl);
}

void CreateROCPressingInFrame(TGeoVolume *mother_volume, Double_t width, Double_t rmin, Double_t rmax)
{
    TString frame_name = "tpc01ROCPressInFr";
    
    Int_t npoints_structure = 8;
    unique_ptr<Double_t[]> x(new Double_t[npoints_structure]);
    unique_ptr<Double_t[]> y(new Double_t[npoints_structure]);
    
    Double_t x_offset = (ROCPress_Backplate_SideWidth - ROCPress_InFrame_width) / TMath::Sin(Section_rad_step / 2);
    Double_t R_seg_in = rmin - x_offset;
    Double_t R_seg_out = rmax - x_offset;
    Double_t x_shift = ROCPress_InFrame_width / TMath::Sin(Section_rad_step / 2);
    x[0] = x[1] = R_seg_out;
    x[2] = x[3] = R_seg_in;
    x[4] = x[5] = R_seg_in + ROCPress_InFrame_width;
    x[6] = x[7] = R_seg_out - ROCPress_InFrame_width;
    y[0] = y[3] = y[4] = y[7] = 0.;
    y[1] = x[1] * TMath::Tan(Section_rad_step / 2);
    y[2] = x[2] * TMath::Tan(Section_rad_step / 2);
    y[5] = (x[5] - x_shift) * TMath::Tan(Section_rad_step / 2);
    y[6] = (x[6] - x_shift) * TMath::Tan(Section_rad_step / 2);
    for (Int_t i = 0; i < npoints_structure; i++)
        x[i] += x_offset;
    
    TGeoXtru *ROCPressInFrameS = new TGeoXtru(2);
    ROCPressInFrameS->DefinePolygon(npoints_structure, x.get(), y.get());
    ROCPressInFrameS->DefineSection(0, 0.);
    ROCPressInFrameS->DefineSection(1, ROCPress_InFrame_thick);
    
    TGeoVolume * ROCPressInFrameV = new TGeoVolume(frame_name, ROCPressInFrameS);
    ROCPressInFrameV->SetMedium(pMedAluminium);
    ROCPressInFrameV->SetLineColor(kBlue + 1);
    
    TGeoCombiTrans * pos = new TGeoCombiTrans(0., 0., ROCPress_Backplate_thick + width, NULL);
    TGeoCombiTrans * pos_refl = new TGeoCombiTrans(*pos);
    pos_refl->ReflectY(true);
    
    mother_volume->AddNode(ROCPressInFrameV, 0, pos);
    mother_volume->AddNode(ROCPressInFrameV, 1, pos_refl);
}

void CreateROCPressingOutFrame(TGeoVolume *mother_volume, Double_t width, Double_t rmin, Double_t rmax)
{
    TString frame_name = "tpc01ROCPressOutFr";
    
    Int_t npoints_structure = 14;
    unique_ptr<Double_t[]> x(new Double_t[npoints_structure]);
    unique_ptr<Double_t[]> y(new Double_t[npoints_structure]);
    
    Double_t x_shift = ROCPress_OutFrame_width / TMath::Sin(Section_rad_step / 2);
    Double_t angle = TMath::ACos(ROCPress_Backplate_HighExtX / ROCPress_Backplate_HighExtL[1]);
    Double_t angle_out = angle + TMath::PiOver2();
    Double_t angle_in = TMath::Pi() - angle - Section_rad_step / 2;
    x[0] = x[1] = rmax - ROCPress_Backplate_HighExtX;
    x[2] = x[3] = rmax;
    x[4] = rmax - ROCPress_Backplate_HighExtL[3] * TMath::Cos(angle);
    x[5] = rmin * (1 + TMath::Tan(Section_rad_step / 2) / (TMath::Tan(ROCPress_OutFrame_LowCenterAngle_rad / 2) - TMath::Tan(Section_rad_step / 2)));
    x[6] = rmin;
    x[7] = rmin + ROCPress_OutFrame_width;
    x[8] = (rmin - x_shift + ROCPress_OutFrame_width / TMath::Sin(ROCPress_OutFrame_LowCenterAngle_rad / 2)) * (1 + TMath::Tan(Section_rad_step / 2) / (TMath::Tan(ROCPress_OutFrame_LowCenterAngle_rad / 2) - TMath::Tan(Section_rad_step / 2))) + x_shift;
    x[9] = x[4] - ROCPress_OutFrame_width / TMath::Sin(angle_in / 2) * TMath::Cos(angle_in / 2 + Section_rad_step / 2);
    x[10] = x[11] = rmax - ROCPress_OutFrame_width;
    x[12] = x[13] = x[0] - ROCPress_OutFrame_width;
    y[0] = y[6] = y[7] = y[13] = 0.;
    y[4] = x[4] * TMath::Tan(Section_rad_step / 2);
    y[5] = x[5] * TMath::Tan(Section_rad_step / 2);
    y[8] = (x[8] - x_shift) * TMath::Tan(Section_rad_step / 2);
    y[9] = (x[9] - x_shift) * TMath::Tan(Section_rad_step / 2);
    y[1] = ROCPress_Backplate_HighExtL[0];
    y[2] = y[1] + ROCPress_Backplate_HighExtL[1] * TMath::Sin(angle);
    y[3] = y[2] + ROCPress_Backplate_HighExtL[2];
    y[12] = y[1] + ROCPress_OutFrame_width / TMath::Tan(angle_out / 2);
    y[11] = y[12] + ROCPress_Backplate_HighExtL[1] * TMath::Sin(angle);
    y[10] = y[11] + ROCPress_Backplate_HighExtL[2] - 2 * ROCPress_OutFrame_width / TMath::Tan(angle_out / 2);
    
    TGeoXtru *ROCPressOutFrameS = new TGeoXtru(2);
    ROCPressOutFrameS->DefinePolygon(npoints_structure, x.get(), y.get());
    ROCPressOutFrameS->DefineSection(0, 0.);
    ROCPressOutFrameS->DefineSection(1, ROCPress_OutFrame_thick);
    
    TGeoVolume * ROCPressOutFrameV = new TGeoVolume(frame_name, ROCPressOutFrameS);
    ROCPressOutFrameV->SetMedium(pMedAluminium);
    ROCPressOutFrameV->SetLineColor(kMagenta);
    
    TGeoCombiTrans * pos = new TGeoCombiTrans(0., 0., ROCPress_Backplate_thick + width, NULL);
    TGeoCombiTrans * pos_refl = new TGeoCombiTrans(*pos);
    pos_refl->ReflectY(true);
    
    mother_volume->AddNode(ROCPressOutFrameV, 0, pos);
    mother_volume->AddNode(ROCPressOutFrameV, 1, pos_refl);
}

void CreateThermoscreenEnds(TGeoVolume *mother_volume, Double_t endcap_zwidth)
{
    TString tscr_name = "tpc01EndThermoScr";
    Double_t x_shift = (Thermoscreen_Ends_GapToNecxtSect / 2) / TMath::Sin(Section_rad_step / 2.);
    Double_t inR = Thermoscreen_Ends_inR - x_shift;
    Double_t outR = Thermoscreen_Ends_outR - x_shift;
    
    TGeoPgon *EndTScrS = new TGeoPgon(tscr_name, -Section_step / 2, Section_step, 1, 2);
    EndTScrS->DefineSection(0, 0., inR, outR);
    EndTScrS->DefineSection(1, Thermoscreen_Thick, inR, outR);
    
    TGeoVolume *EndTScrV = new TGeoVolume(tscr_name, EndTScrS);
    EndTScrV->SetMedium(pMedAluminium);
    EndTScrV->SetLineColor(kMagenta - 1);
    
    for(Int_t isec = 0; isec < Nsections; isec++)
    {
        TGeoCombiTrans *tscr_pos= new TGeoCombiTrans();
        tscr_pos->SetTranslation(x_shift, 0., -endcap_zwidth / 2 + RibsZfromTPCends + RibProfileWidth / 2);
        tscr_pos->RotateZ(Section_step * isec);

        mother_volume->AddNode(EndTScrV, isec, tscr_pos);
    }
}

void CreateThermoscreenRadial(TGeoVolume *mother_volume)
{
    TString plate_name = "tpc01ThermoScrRad";
    
    Int_t npoints_screen = 8;
    unique_ptr<Double_t[]> xs(new Double_t[npoints_screen]);
    unique_ptr<Double_t[]> ys(new Double_t[npoints_screen]);
    xs[0] = Thermoscreen_Rad_Plate_Width_Seg / 2;
    xs[1] = xs[0] + Thermoscreen_Rad_Plate_Width_Seg * TMath::Cos(Thermoscreen_Rad_Seg_Angle_Rad / 2);
    xs[2] = xs[1] + Thermoscreen_Thick * TMath::Sin(Thermoscreen_Rad_Seg_Angle_Rad / 2);
    xs[3] = xs[0] + Thermoscreen_Thick * TMath::Tan(Thermoscreen_Rad_Seg_Angle_Rad / 4);
    xs[7] = -xs[0];
    xs[6] = -xs[1];
    xs[5] = -xs[2];
    xs[4] = -xs[3];
    ys[0] = ys[7] = 0.;
    ys[1] = ys[6] = -Thermoscreen_Rad_Plate_Width_Seg * TMath::Sin(Thermoscreen_Rad_Seg_Angle_Rad / 2);
    ys[2] = ys[5] = ys[1] + Thermoscreen_Thick * TMath::Cos(Thermoscreen_Rad_Seg_Angle_Rad / 2);
    ys[3] = ys[4] = Thermoscreen_Thick;
    
    TGeoXtru * CoolingPlateS = new TGeoXtru(2);
    CoolingPlateS->DefinePolygon(npoints_screen, xs.get(), ys.get());
    CoolingPlateS->DefineSection(0, Thermoscreen_Rad_Plate_ZGap / 2);
    CoolingPlateS->DefineSection(1, Thermoscreen_Rad_Plate_ZGap / 2 + Thermoscreen_Rad_Plate_Lenght);
    
    TGeoVolume *CoolingPlateV = new TGeoVolume(plate_name, CoolingPlateS);
    CoolingPlateV->SetMedium(pMedAluminium);
    CoolingPlateV->SetLineColor(kMagenta - 1);
    
    Int_t nsec = Nsections - 2; // 12 sectors - 2 with holders
    Int_t nzp = 4; // 4 panels per sector
    for (Int_t isec = 0; isec < nsec; ++isec)
        for (Int_t zp = 0; zp < nzp; ++zp)
        {
            TGeoCombiTrans * pos = new TGeoCombiTrans();
            pos->SetTranslation(0., ThermoscreenRadialPlateROffset(), (Thermoscreen_Rad_Plate_Lenght + Thermoscreen_Rad_Plate_ZGap) * (zp - 2));
            pos->RotateZ(-90. + Thermoscreen_Rad_Seg_Angle * (isec < nsec / 2 ? isec + 1 : isec + 2));

            mother_volume->AddNode(CoolingPlateV, 10 * isec + zp, pos);
        }
    
    TString platerail_name = "tpc01ThermoScrRadRail";
    
    Int_t npoints_screen_rail = 12;
    unique_ptr<Double_t[]> xr(new Double_t[npoints_screen_rail]);
    unique_ptr<Double_t[]> yr(new Double_t[npoints_screen_rail]);
    xr[0] = Thermoscreen_Rad_Rail_CentralPlate_Width / 2 + (Thermoscreen_Rad_Rail_Bend_R + Thermoscreen_Thick) * TMath::Tan((TMath::PiOver2() - Thermoscreen_Rad_Rail_CentralPlate_SideAngle_Rad) / 2);
    yr[0] = 0.;
    xr[5] = Thermoscreen_Rad_Rail_CentralPlate_Width / 2 + Thermoscreen_Rad_Rail_Bend_R * TMath::Tan((TMath::PiOver2() - Thermoscreen_Rad_Rail_CentralPlate_SideAngle_Rad) / 2);
    yr[5] = Thermoscreen_Thick;
    TVector2 vec(0., ThermoscreenRadialPlateRailBendHeight(Thermoscreen_Rad_Rail_Bend_R + Thermoscreen_Thick + Thermoscreen_Rad_Rail_BendCenters_Dist, Thermoscreen_Rad_Rail_Bend_R));
    vec = vec.Rotate(-Thermoscreen_Rad_Rail_CentralPlate_SideAngle_Rad);
    xr[1] = xr[0] + vec.X();
    yr[1] = yr[0] + vec.Y();
    vec.Set(0., ThermoscreenRadialPlateRailBendHeight(Thermoscreen_Rad_Rail_Bend_R + Thermoscreen_Rad_Rail_BendCenters_Dist, Thermoscreen_Rad_Rail_Bend_R + Thermoscreen_Thick));
    vec = vec.Rotate(-Thermoscreen_Rad_Rail_CentralPlate_SideAngle_Rad);
    xr[4] = xr[5] + vec.X();
    yr[4] = yr[5] + vec.Y();
    vec.Set(Thermoscreen_Rad_Rail_SidePlate_Width, 0.);
    vec = vec.Rotate(-Thermoscreen_Rad_Seg_Angle_Rad / 2);
    xr[2] = xr[1] + vec.X();
    yr[2] = yr[1] + vec.Y();
    xr[3] = xr[2] + Thermoscreen_Thick * TMath::Sin(Thermoscreen_Rad_Seg_Angle_Rad / 2);
    yr[3] = yr[2] + Thermoscreen_Thick * TMath::Cos(Thermoscreen_Rad_Seg_Angle_Rad / 2);
    
    for (Int_t pt = 0; pt < npoints_screen_rail / 2; ++pt)
    {
        xr[npoints_screen_rail - 1 - pt] = -xr[pt];
        yr[npoints_screen_rail - 1 - pt] = yr[pt];
    }
    
    TGeoXtru * CoolingPlateRailS = new TGeoXtru(2);
    CoolingPlateRailS->DefinePolygon(npoints_screen_rail, xr.get(), yr.get());
    CoolingPlateRailS->DefineSection(0, Thermoscreen_Rad_Plate_ZGap / 2);
    CoolingPlateRailS->DefineSection(1, Thermoscreen_Rad_Plate_ZGap / 2 + Thermoscreen_Rad_Plate_Lenght);
    
    TGeoVolume * CoolingPlateRailV = new TGeoVolume(platerail_name, CoolingPlateRailS);
    CoolingPlateRailV->SetMedium(pMedAluminium);
    CoolingPlateRailV->SetLineColor(kMagenta - 1);
    
    TGeoXtru * CoolingPlateRailCBaseS = new TGeoXtru(2);
    CoolingPlateRailCBaseS->DefinePolygon(npoints_screen_rail, xr.get(), yr.get());
    CoolingPlateRailCBaseS->DefineSection(0, Thermoscreen_Rad_Plate_ZGap / 2);
    CoolingPlateRailCBaseS->DefineSection(1, Thermoscreen_Rad_Plate_ZGap / 2 + Thermoscreen_Rad_Plate_Lenght - Thermoscreen_Rad_Rail_Lenght_Cut);
    
    TGeoVolume * CoolingPlateRailCBaseV = new TGeoVolume(platerail_name + "CutBase", CoolingPlateRailCBaseS);
    CoolingPlateRailCBaseV->SetMedium(pMedAluminium);
    CoolingPlateRailCBaseV->SetLineColor(kMagenta - 1);
    
    Int_t npoints_screen_rail_cut = 4;
    unique_ptr<Double_t[]> xrc(new Double_t[npoints_screen_rail_cut]);
    unique_ptr<Double_t[]> yrc(new Double_t[npoints_screen_rail_cut]);
    for (Int_t pt = 0; pt < npoints_screen_rail_cut; ++pt)
    {
        xrc[pt] = xr[pt + 1];
        yrc[pt] = yr[pt + 1];
    }
    
    TGeoXtru * CoolingPlateRailCCutS = new TGeoXtru(2);
    CoolingPlateRailCCutS->DefinePolygon(npoints_screen_rail_cut, xrc.get(), yrc.get());
    CoolingPlateRailCCutS->DefineSection(0, Thermoscreen_Rad_Plate_ZGap / 2 + Thermoscreen_Rad_Plate_Lenght - Thermoscreen_Rad_Rail_Lenght_Cut);
    CoolingPlateRailCCutS->DefineSection(1, Thermoscreen_Rad_Plate_ZGap / 2 + Thermoscreen_Rad_Plate_Lenght);
    
    TGeoVolume * CoolingPlateRailCCutV = new TGeoVolume(platerail_name + "CutSide", CoolingPlateRailCCutS);
    CoolingPlateRailCCutV->SetMedium(pMedAluminium);
    CoolingPlateRailCCutV->SetLineColor(kMagenta - 1);
    
    TGeoVolumeAssembly * CoolingPlateRailSCutV = new TGeoVolumeAssembly(platerail_name + "Cut");
    CoolingPlateRailSCutV->AddNode(CoolingPlateRailCBaseV, 0);
    CoolingPlateRailSCutV->AddNode(CoolingPlateRailCCutV, 0);
    
    for (Int_t zp = 0; zp < nzp; ++zp)
    {
        if(zp == 0 || zp == nzp - 1)
        {
            TGeoCombiTrans * pos = new TGeoCombiTrans();
            if(zp == 0)
            {
                pos->RotateY(180.);
                pos->SetTranslation(0., Thermoscreen_Rad_Rail_CentralPlate_DistToZaxis, (Thermoscreen_Rad_Plate_Lenght + Thermoscreen_Rad_Plate_ZGap) * (zp - 1));
            }
            else
            {
                pos->SetTranslation(0., Thermoscreen_Rad_Rail_CentralPlate_DistToZaxis, (Thermoscreen_Rad_Plate_Lenght + Thermoscreen_Rad_Plate_ZGap) * (zp - 2));
                pos->ReflectX(true);
            }
            pos->RotateZ(90.);
            TGeoCombiTrans * refl = new TGeoCombiTrans(*pos);
            refl->ReflectX(true);
            
            mother_volume->AddNode(CoolingPlateRailSCutV, zp, pos);
            mother_volume->AddNode(CoolingPlateRailSCutV, 10 + zp, refl);
        }
        else
        {
            TGeoCombiTrans * pos = new TGeoCombiTrans();
            pos->SetTranslation(0., Thermoscreen_Rad_Rail_CentralPlate_DistToZaxis, (Thermoscreen_Rad_Plate_Lenght + Thermoscreen_Rad_Plate_ZGap) * (zp - 2));
            pos->RotateZ(90.);
            TGeoCombiTrans * refl = new TGeoCombiTrans(*pos);
            refl->ReflectX(true);

            mother_volume->AddNode(CoolingPlateRailV, zp, pos);
            mother_volume->AddNode(CoolingPlateRailV, 10 + zp, refl);
        }
    }
    
    TString hold_name = "tpc01ThermoScrHold";
    
    Int_t npoints_structure = 12;
    unique_ptr<Double_t[]> x(new Double_t[npoints_structure]);
    unique_ptr<Double_t[]> y(new Double_t[npoints_structure]);
    y[0] = y[3] = y[4] = y[11] = Thermoscreen_Rad_Hold_Thick / 2;
    y[5] = y[6] = y[9] = y[10] = -y[0];
    y[1] = y[2] = Thermoscreen_Rad_Hold_Width / 2;
    y[7] = y[8] = -y[1];
    x[0] = x[1] = x[8] = x[9] = 0.;
    x[10] = x[11] = Thermoscreen_Rad_Hold_HeightIn;
    x[2] = x[3] = x[6] = x[7] = -Thermoscreen_Rad_Hold_Thick;
    x[4] = x[5] = -(Thermoscreen_Rad_Hold_Height - Thermoscreen_Rad_Hold_HeightIn);
    
    TGeoXtru *TScrHoldS = new TGeoXtru(2);
    TScrHoldS->DefinePolygon(npoints_structure, x.get(), y.get());
    TScrHoldS->DefineSection(0, -Thermoscreen_Rad_Length / 2);
    TScrHoldS->DefineSection(1, Thermoscreen_Rad_Length / 2);
    
    TGeoVolume * TScrHoldV = new TGeoVolume(hold_name, TScrHoldS);
    TScrHoldV->SetMedium(pMedAluminium);
    TScrHoldV->SetLineColor(kWhite);
    
    for(Int_t isec = 0; isec < Nsections; isec++)
    {
        TGeoCombiTrans * pos = new TGeoCombiTrans();
        pos->SetTranslation(Thermoscreen_Rad_R, 0., 0.);
        pos->RotateZ(90. + Thermoscreen_Rad_Seg_Angle / 2 + Thermoscreen_Rad_Seg_Angle * isec);
        
        mother_volume->AddNode(TScrHoldV, isec, pos);
    }
}

Double_t ThermoscreenRadialPlateROffset()
{
    Double_t a = Thermoscreen_Rad_Plate_Width_Seg + Thermoscreen_Rad_Plate_Width_Seg / 2 / TMath::Cos(Thermoscreen_Rad_Seg_Angle_Rad / 2);
    Double_t b = Thermoscreen_Rad_Plate_Width_Seg / TMath::Tan(Thermoscreen_Rad_Seg_Angle_Rad / 2);
    Double_t c = Thermoscreen_Rad_Plate_Width_Seg / TMath::Sin(Thermoscreen_Rad_Seg_Angle_Rad / 2);
    return c - Thermoscreen_Rad_Plate_Width_Seg / 2 * TMath::Tan(Thermoscreen_Rad_Seg_Angle_Rad / 2) + (Thermoscreen_Rad_R - b) / TMath::Cos(Thermoscreen_Rad_Seg_Angle_Rad / 2);
}

inline Double_t ThermoscreenRadialPlateRailBendHeight(Double_t R1, Double_t R2)
{
    return R1 / TMath::Cos(Thermoscreen_Rad_Rail_CentralPlate_SideAngle_Rad) + (R2 * TMath::Tan(Thermoscreen_Rad_Rail_CentralPlate_SideAngle_Rad) + R2 * TMath::Tan((TMath::PiOver2() + Thermoscreen_Rad_Seg_Angle_Rad / 2 - Thermoscreen_Rad_Rail_CentralPlate_SideAngle_Rad) / 2));
}

TGeoVolume * CreateMotherVolTPCandThermoscreen()
{
    TString mother_vol_base = "base_vol";
    TString mother_vol_sub = "subtrahend_vol";
    
    Int_t pt_per_sect = 4;
    Int_t npoints_structure = Nsections * pt_per_sect + 4;
    unique_ptr<Double_t[]> x(new Double_t[npoints_structure]);
    unique_ptr<Double_t[]> y(new Double_t[npoints_structure]);
    for (Int_t isec = 0; isec < Nsections / 2; ++isec)
    {
        std::vector<TVector2> pt;
        pt.reserve(pt_per_sect);
        Double_t xpos = Thermoscreen_Rad_Plate_Width_Seg / 2 + Thermoscreen_Thick * TMath::Tan(Thermoscreen_Rad_Seg_Angle_Rad / 4);
        pt[1].Set(Thermoscreen_Rad_Hold_Thick / 2, Thermoscreen_Rad_R + Thermoscreen_Rad_Hold_HeightIn);
        pt[2].Set(-Thermoscreen_Rad_Hold_Thick / 2, Thermoscreen_Rad_R + Thermoscreen_Rad_Hold_HeightIn);
        pt[3].Set(xpos, ThermoscreenRadialPlateROffset() + Thermoscreen_Thick);
        pt[0].Set(-xpos, ThermoscreenRadialPlateROffset() + Thermoscreen_Thick);
        pt[1] = pt[1].Rotate((isec - 2.5) * Thermoscreen_Rad_Seg_Angle_Rad);
        pt[2] = pt[2].Rotate((isec - 2.5) * Thermoscreen_Rad_Seg_Angle_Rad);
        pt[3] = pt[3].Rotate((isec - 2) * Thermoscreen_Rad_Seg_Angle_Rad);
        pt[0] = pt[0].Rotate((isec - 3) * Thermoscreen_Rad_Seg_Angle_Rad);
        for (Int_t i = 0; i < pt_per_sect; ++i)
        {
            x[pt_per_sect * isec + i + 1] = pt[i].X();
            y[pt_per_sect * isec + i + 1] = pt[i].Y();
        }
    }
    x[0] = Thermoscreen_Rad_Rail_CentralPlate_DistToZaxis + Thermoscreen_Thick;
    x[(Nsections * pt_per_sect) / 2 + 1] = - x[0];
    y[0] = y[(Nsections * pt_per_sect) / 2 + 1] = Thermoscreen_Rad_Rail_CentralPlate_Width / 2 + Thermoscreen_Rad_Rail_Bend_R * TMath::Tan((TMath::PiOver2() - Thermoscreen_Rad_Rail_CentralPlate_SideAngle_Rad) / 2);
    TVector2 vec(ThermoscreenRadialPlateRailBendHeight(Thermoscreen_Rad_Rail_Bend_R + Thermoscreen_Rad_Rail_BendCenters_Dist, Thermoscreen_Rad_Rail_Bend_R + Thermoscreen_Thick), 0);
    vec = vec.Rotate(Thermoscreen_Rad_Rail_CentralPlate_SideAngle_Rad);
    x[1] = x[0] + vec.X();
    y[1] = y[0] + vec.Y();
    x[(Nsections * pt_per_sect) / 2] = -x[1];
    y[(Nsections * pt_per_sect) / 2] = y[1];
    for (Int_t ipt = 0; ipt < npoints_structure / 2; ++ipt)
    {
        x[npoints_structure - 1 - ipt] = x[ipt];
        y[npoints_structure - 1 - ipt] = -y[ipt];
    }
    
    TGeoXtru * mother_volume_main = new TGeoXtru(2);
    mother_volume_main->DefinePolygon(npoints_structure, x.get(), y.get());
    mother_volume_main->DefineSection(0, -Thermoscreen_Rad_Length / 2);
    mother_volume_main->DefineSection(1, Thermoscreen_Rad_Length / 2);
    mother_volume_main->SetName(mother_vol_base);
    
    TGeoPcon *mother_volume_substr = new TGeoPcon (mother_vol_sub, 0., 360., 6);
    mother_volume_substr -> DefineSection(0, -Thermoscreen_Rad_Length / 2., 0., Thermoscreen_Rad_Rail_CentralPlate_DistToZaxis);
    mother_volume_substr -> DefineSection(1, -Container_tpc_z / 2., 0., Thermoscreen_Rad_Rail_CentralPlate_DistToZaxis);
    mother_volume_substr -> DefineSection(2, -Container_tpc_z / 2., 0., TpcInnerRadius);
    mother_volume_substr -> DefineSection(3, Container_tpc_z / 2., 0., TpcInnerRadius);
    mother_volume_substr -> DefineSection(4, Container_tpc_z / 2., 0., Thermoscreen_Rad_Rail_CentralPlate_DistToZaxis);
    mother_volume_substr -> DefineSection(5, Thermoscreen_Rad_Length / 2., 0.,  Thermoscreen_Rad_Rail_CentralPlate_DistToZaxis);
    
    TString top_mother_volume = "TPC_mother_volume";
    TGeoCompositeShape *TOP_mother_volumeS = new TGeoCompositeShape(top_mother_volume, mother_vol_base + "-" + mother_vol_sub);
    return new TGeoVolume(top_mother_volume, TOP_mother_volumeS);
    //return new TGeoVolume(top_mother_volume, mother_volume_main);
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

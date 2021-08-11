/*
 * Cradle for MPD setup
 *
 * Use this macro to create ROOT geometry of the cradle for the MPD setup as
 * the assembly of passive volumes
 *
 * Author: Baranov D.
 * Created: 10.08.2021
 * Updated: ***
 *
 * WARNING: all units is in cm!!!
 */

using namespace TMath;

//Cradle XYZ position
const Double_t XCradlePosition = 0.0; //cm
const Double_t YCradlePosition = 0.0; //cm
const Double_t ZCradlePosition = 0.0; //cm

//GeoManager
TGeoManager* gGeoMan = 0;

//media
TGeoMedium *pMedAir = 0;
TGeoMedium *pMedIron = 0;

class FairGeoMedia;
class FairGeoBuilder;

TGeoVolume *CreateSupportStand(TString support_stand_name);
TGeoVolume *CreateHorizontalBeam(TString support_beam_name);


void DefineRequiredMedia(FairGeoMedia* geoMedia, FairGeoBuilder* geoBuild) {

    //air medium
    FairGeoMedium* mAir = geoMedia->getMedium("air");
    if ( ! mAir ) Fatal("Main", "FairMedium air not found");
    geoBuild->createMedium(mAir);
    pMedAir = gGeoManager->GetMedium("air");
    if ( ! pMedAir ) Fatal("Main", "Medium air not found");

    //carbon medium
    FairGeoMedium* mIron = geoMedia->getMedium("iron");
    if ( ! mIron  ) Fatal("Main", "FairMedium carbon not found");
    geoBuild->createMedium(mIron);
    pMedIron = gGeoManager->GetMedium("iron");
    if ( ! pMedIron ) Fatal("Main", "Medium iron not found");
}

//create geometry of silicon detector with passive frames
void create_rootgeom_cradle() {

    // ----  set working directory  --------------------------------------------
    TString gPath = gSystem->Getenv("VMCWORKDIR");

    // -------   Geometry file name (output)   ----------------------------------
    const TString geoDetectorName = "cradle";
    const TString geoDetectorVersion = "v1";
    const TString geoFileName = gPath + "/geometry/" + geoDetectorName + "_"+ geoDetectorVersion + ".root";

    // ----  global geometry parameters  ---------------------------------------
    FairGeoLoader*    geoLoad = new FairGeoLoader("TGeo","FairGeoLoader");
    FairGeoInterface* geoFace = geoLoad->getGeoInterface();

    // -------   Load media from media file   ----------------------------------
    TString medFile = gPath + "/geometry/media.geo";
    geoFace->setMediaFile(medFile);
    geoFace->readMedia();
    FairGeoMedia*   geoMedia = geoFace->getMedia();
    FairGeoBuilder* geoBuild = geoLoad->getGeoBuilder();

    DefineRequiredMedia(geoMedia, geoBuild);
    // -------------------------------------------------------------------------

    // --------------   Create geometry and top volume  -------------------------
    gGeoMan = (TGeoManager*)gROOT->FindObject("FAIRGeom");
    gGeoMan->SetName(geoDetectorName + "_geom");
    TGeoVolume* top = new TGeoVolumeAssembly("TOP");
    top->SetMedium(pMedAir);
    gGeoMan->SetTopVolume(top);
    //gGeoMan->SetTopVisible(1);
    // --------------------------------------------------------------------------

    // Define TOP Geometry
    TGeoVolume* CRADLE = new TGeoVolumeAssembly(geoDetectorName);
    CRADLE->SetMedium(pMedAir);

    TGeoVolume *supportStand1 = CreateSupportStand("supportStand1");
    TGeoVolume *supportStand2 = CreateSupportStand("supportStand2");

    TGeoVolume *horizontalBeam1 = CreateHorizontalBeam("horizontalBeam1");
    TGeoVolume *horizontalBeam2 = CreateHorizontalBeam("horizontalBeam2");

    TGeoCombiTrans *supportStand1_transform = new TGeoCombiTrans();
        supportStand1_transform->SetDx(XCradlePosition);
        supportStand1_transform->SetDy(YCradlePosition);
        supportStand1_transform->SetDz(-(505.0*0.5 + 95.5*0.5) + ZCradlePosition);

    TGeoCombiTrans *supportStand2_transform = new TGeoCombiTrans();
        supportStand2_transform->SetDx(XCradlePosition);
        supportStand2_transform->SetDy(YCradlePosition);
        supportStand2_transform->SetDz(+(505.5*0.5 + 95.5*0.5) + ZCradlePosition);

    TGeoCombiTrans *horizontalBeam1_transform = new TGeoCombiTrans();
        horizontalBeam1_transform->SetDx(+(670.0*0.5 - 60.0*0.5 - 5.0) + XCradlePosition);
        horizontalBeam1_transform->SetDy(-425.5 + 60.0*0.5 + 5.0 + YCradlePosition);
        horizontalBeam1_transform->SetDz(0.0 + ZCradlePosition);

    TGeoCombiTrans *horizontalBeam2_transform = new TGeoCombiTrans();
        horizontalBeam2_transform->SetDx(-(670.0*0.5 - 60.0*0.5 - 5.0) + XCradlePosition);
        horizontalBeam2_transform->SetDy(-425.5 + 60.0*0.5 + 5.0 + YCradlePosition);
        horizontalBeam2_transform->SetDz(ZCradlePosition);

    CRADLE->AddNode(supportStand1, 0, supportStand1_transform);
    CRADLE->AddNode(supportStand2, 0, supportStand2_transform);

    CRADLE->AddNode(horizontalBeam1, 0, horizontalBeam1_transform);
    CRADLE->AddNode(horizontalBeam2, 0, horizontalBeam2_transform);

    top->AddNode(CRADLE, 0);
    top->SetVisContainers(kTRUE);

    // ---------------   Finish   -----------------------------------------------
    gGeoMan->CloseGeometry();
    gGeoMan->CheckOverlaps(0.0001);
    gGeoMan->PrintOverlaps();
    gGeoMan->Test();

    TFile* geoFile = new TFile(geoFileName, "RECREATE");
    top->Write();
    geoFile->Close();
    top->Draw("ogl");
}
//------------------------------------------------------------------------------

TGeoVolume *CreateSupportStand(TString support_stand_name) {

    //support stand volume
    TGeoVolume *supportStand = new TGeoVolumeAssembly(support_stand_name);
    supportStand->SetMedium(pMedAir);

    //front wall ---------------------------------------------------------------
    Double_t frontWall_ZSize = 8.0; //cm

    Double_t frontWall_ZDistBetweenParts = 95.5 - 2*frontWall_ZSize; //cm (79.5 cm)

    TGeoShape *frontWallHalfS = new TGeoXtru(2);
    frontWallHalfS->SetName(TString("frontWallHalfS")+=TString("_") + supportStand->GetName());
    {
        Double_t x_pos[12] = {
            +335.0, //0
            +335.0, //1
            +373.14, //2
            +329.1290, //3
            +314.6443, //4
            +282.8608, //5
            +236.5852, //6
            +178.2448, //7
            +110.7539, //8
            +37.5528, //9
            0.0, //10
            0.0 //11
        };
        Double_t y_pos[12] = {
            -425.5, //0
            -223.5, //1
            -41.5, //2
            -41.5, //3
            -105.2784, //4
            -173.3372, //5
            -232.4915, //6
            -279.7879, //7
            -312.7593, //8
            -329.5977, //9
            -329.5977, //10
            -425.5 //11
        };
        ((TGeoXtru*)frontWallHalfS)->DefinePolygon(12, x_pos, y_pos);
        ((TGeoXtru*)frontWallHalfS)->DefineSection(0, -frontWall_ZSize*0.5);
        ((TGeoXtru*)frontWallHalfS)->DefineSection(1, +frontWall_ZSize*0.5);
    }

    TGeoVolume *frontWallHalfV = new TGeoVolume(TString("frontWallHalfV")+=TString("_") + supportStand->GetName(), frontWallHalfS);

    //volume medium
    TGeoMedium *frontWallHalfV_medium = pMedIron;
    if(frontWallHalfV_medium) {
        frontWallHalfV->SetMedium(frontWallHalfV_medium);
    }
    else Fatal("Main", "Invalid medium for frontWallHalfV!");

    //volume visual property (transparency)
    frontWallHalfV->SetLineColor(TColor::GetColor("#ffcccc"));
    frontWallHalfV->SetTransparency(0);

    TGeoCombiTrans *frontWallHalfV_transf[4];

    frontWallHalfV_transf[0] = new TGeoCombiTrans();
    frontWallHalfV_transf[0]->SetDx(0.0);
    frontWallHalfV_transf[0]->SetDy(0.0);
    frontWallHalfV_transf[0]->SetDz(-(frontWall_ZSize*0.5 + frontWall_ZDistBetweenParts*0.5));

    frontWallHalfV_transf[1] = new TGeoCombiTrans();
    frontWallHalfV_transf[1]->ReflectX(true);
    frontWallHalfV_transf[1]->SetDx(0.0);
    frontWallHalfV_transf[1]->SetDy(0.0);
    frontWallHalfV_transf[1]->SetDz(-(frontWall_ZSize*0.5 + frontWall_ZDistBetweenParts*0.5));

    frontWallHalfV_transf[2] = new TGeoCombiTrans();
    frontWallHalfV_transf[2]->SetDx(0.0);
    frontWallHalfV_transf[2]->SetDy(0.0);
    frontWallHalfV_transf[2]->SetDz(+(frontWall_ZSize*0.5 + frontWall_ZDistBetweenParts*0.5));

    frontWallHalfV_transf[3] = new TGeoCombiTrans();
    frontWallHalfV_transf[3]->ReflectX(true);
    frontWallHalfV_transf[3]->SetDx(0.0);
    frontWallHalfV_transf[3]->SetDy(0.0);
    frontWallHalfV_transf[3]->SetDz(+(frontWall_ZSize*0.5 + frontWall_ZDistBetweenParts*0.5));

    supportStand->AddNode(frontWallHalfV, 0, frontWallHalfV_transf[0]);
    supportStand->AddNode(frontWallHalfV, 1, frontWallHalfV_transf[1]);
    supportStand->AddNode(frontWallHalfV, 2, frontWallHalfV_transf[2]);
    supportStand->AddNode(frontWallHalfV, 3, frontWallHalfV_transf[3]);
    //--------------------------------------------------------------------------

    //support base -------------------------------------------------------------
    Double_t supportBase_XSize = 680.0; //cm
    Double_t supportBase_YSize = 10.0; //cm
    Double_t supportBase_ZSize = 95.5; //cm

    Double_t supportBaseExtension_XSize = 38.5; //cm
    Double_t supportBaseExtension_YSize = 10.0; //cm
    Double_t supportBaseExtension_ZMinSize = 79.5; //cm
    Double_t supportBaseExtension_ZMaxSize = 95.5; //cm

    TGeoShape *supportBaseS = new TGeoBBox(TString("supportBaseS")+=TString("_") + supportStand->GetName(), supportBase_XSize*0.5, supportBase_YSize*0.5, supportBase_ZSize*0.5);
    TGeoShape *supportBaseExtensionS = new TGeoTrd1(TString("supportBaseExtensionS")+=TString("_") + supportStand->GetName(), supportBaseExtension_ZMinSize*0.5,supportBaseExtension_ZMaxSize*0.5, supportBaseExtension_YSize*0.5, supportBaseExtension_XSize*0.5);

    TGeoVolume *supportBaseV = new TGeoVolume(TString("supportBaseV")+=TString("_") + supportStand->GetName(), supportBaseS);
    TGeoVolume *supportBaseExtensionV = new TGeoVolume(TString("supportBaseExtensionV")+=TString("_") + supportStand->GetName(), supportBaseExtensionS);

    //volume medium
    TGeoMedium *supportBaseV_medium = pMedIron;
    if(supportBaseV_medium) {
        supportBaseV->SetMedium(supportBaseV_medium);
        supportBaseExtensionV->SetMedium(supportBaseV_medium);
    }
    else Fatal("Main", "Invalid medium for supportBaseV!");

    //volume visual property (transparency)
    supportBaseV->SetLineColor(TColor::GetColor("#ffcccc"));
    supportBaseV->SetTransparency(0);
    supportBaseExtensionV->SetLineColor(TColor::GetColor("#ffcccc"));
    supportBaseExtensionV->SetTransparency(0);

    TGeoCombiTrans *supportBaseV_transf[1];

    supportBaseV_transf[0] = new TGeoCombiTrans();
    supportBaseV_transf[0]->SetDx(0.0);
    supportBaseV_transf[0]->SetDy(-supportBase_YSize*0.5 - 425.5/*shift*/);
    supportBaseV_transf[0]->SetDz(0.0);

    TGeoCombiTrans *supportBaseExtensionV_transf[1];

    supportBaseExtensionV_transf[0] = new TGeoCombiTrans();
    supportBaseExtensionV_transf[0]->RotateY(-90.0);
    supportBaseExtensionV_transf[0]->SetDx(supportBase_XSize*0.5 + supportBaseExtension_XSize*0.5);
    supportBaseExtensionV_transf[0]->SetDy(-supportBase_YSize*0.5 - 425.5/*shift*/);
    supportBaseExtensionV_transf[0]->SetDz(0.0);

    supportStand->AddNode(supportBaseV, 0, supportBaseV_transf[0]);
    supportStand->AddNode(supportBaseExtensionV, 0, supportBaseExtensionV_transf[0]);
    //--------------------------------------------------------------------------

    //lower side wall ----------------------------------------------------------
    Double_t lowerSideWall_XSize = 5.0; //cm
    Double_t lowerSideWall_YSize = 200.0; //cm
    Double_t lowerSideWall_ZSize = 95.5; //cm

    TGeoShape *lowerSideWallS = new TGeoBBox(TString("lowerSideWallS")+=TString("_") + supportStand->GetName(), lowerSideWall_XSize*0.5, lowerSideWall_YSize*0.5, lowerSideWall_ZSize*0.5);

    TGeoVolume *lowerSideWallV = new TGeoVolume(TString("lowerSideWallV")+=TString("_") + supportStand->GetName(), lowerSideWallS);

    //volume medium
    TGeoMedium *lowerSideWallV_medium = pMedIron;
    if(lowerSideWallV_medium) {
        lowerSideWallV->SetMedium(lowerSideWallV_medium);
    }
    else Fatal("Main", "Invalid medium for lowerSideWallV!");

    //volume visual property (transparency)
    lowerSideWallV->SetLineColor(TColor::GetColor("#ffcccc"));
    lowerSideWallV->SetTransparency(0);

    TGeoCombiTrans *lowerSideWallV_transf[2];

    lowerSideWallV_transf[0] = new TGeoCombiTrans();
    lowerSideWallV_transf[0]->SetDx(-(supportBase_XSize*0.5 - lowerSideWall_XSize*0.5));
    lowerSideWallV_transf[0]->SetDy(lowerSideWall_YSize*0.5 - 425.5/*shift*/);
    lowerSideWallV_transf[0]->SetDz(0.0);

    lowerSideWallV_transf[1] = new TGeoCombiTrans();
    lowerSideWallV_transf[1]->SetDx(+(supportBase_XSize*0.5 - lowerSideWall_XSize*0.5));
    lowerSideWallV_transf[1]->SetDy(lowerSideWall_YSize*0.5 - 425.5/*shift*/);
    lowerSideWallV_transf[1]->SetDz(0.0);

    supportStand->AddNode(lowerSideWallV, 0, lowerSideWallV_transf[0]);
    supportStand->AddNode(lowerSideWallV, 1, lowerSideWallV_transf[1]);
    //--------------------------------------------------------------------------

    //upper side wall ----------------------------------------------------------
    Double_t upperSideWall_XSize = 5.0; //cm
    Double_t upperSideWall_YSize = 186.0; //cm
    Double_t upperSideWall_ZSize = 95.5; //cm

    TGeoShape *upperSideWallS = new TGeoBBox(TString("upperSideWallS")+=TString("_") + supportStand->GetName(), upperSideWall_XSize*0.5, upperSideWall_YSize*0.5, upperSideWall_ZSize*0.5);

    TGeoVolume *upperSideWallV = new TGeoVolume(TString("upperSideWallV")+=TString("_") + supportStand->GetName(), upperSideWallS);

    //volume medium
    TGeoMedium *upperSideWallV_medium = pMedIron;
    if(upperSideWallV_medium) {
        upperSideWallV->SetMedium(upperSideWallV_medium);
    }
    else Fatal("Main", "Invalid medium for upperSideWallV!");

    //volume visual property (transparency)
    upperSideWallV->SetLineColor(TColor::GetColor("#ffcccc"));
    upperSideWallV->SetTransparency(0);

    TGeoCombiTrans *upperSideWallV_transf[2];

    upperSideWallV_transf[0] = new TGeoCombiTrans();
    upperSideWallV_transf[0]->RotateZ(11.83567);
    upperSideWallV_transf[0]->SetDx(-(supportBase_XSize*0.5 - upperSideWall_XSize*0.5 + 19.07/*shift*/));
    upperSideWallV_transf[0]->SetDy(-upperSideWall_YSize*0.5 - 41.5/*shift*/ + 1.1/*shift*/);
    upperSideWallV_transf[0]->SetDz(0.0);

    upperSideWallV_transf[1] = new TGeoCombiTrans();
    upperSideWallV_transf[1]->RotateZ(-11.83567);
    upperSideWallV_transf[1]->SetDx(+(supportBase_XSize*0.5 - upperSideWall_XSize*0.5 + 19.07/*shift*/));
    upperSideWallV_transf[1]->SetDy(-upperSideWall_YSize*0.5 - 41.5/*shift*/ + 1.1/*shift*/);
    upperSideWallV_transf[1]->SetDz(0.0);

    supportStand->AddNode(upperSideWallV, 0, upperSideWallV_transf[0]);
    supportStand->AddNode(upperSideWallV, 1, upperSideWallV_transf[1]);
    //--------------------------------------------------------------------------

    //stand cover --------------------------------------------------------------
    Double_t standCover_XSize = 49.0; //cm
    Double_t standCover_YSize = 3.0; //cm
    Double_t standCover_ZSize = 95.5; //cm

    TGeoShape *standCoverS = new TGeoBBox(TString("standCoverS")+=TString("_") + supportStand->GetName(), standCover_XSize*0.5, standCover_YSize*0.5, standCover_ZSize*0.5);

    TGeoVolume *standCoverV = new TGeoVolume(TString("standCoverV")+=TString("_") + supportStand->GetName(), standCoverS);

    //volume medium
    TGeoMedium *standCoverV_medium = pMedIron;
    if(standCoverV_medium) {
        standCoverV->SetMedium(standCoverV_medium);
    }
    else Fatal("Main", "Invalid medium for standCoverV!");

    //volume visual property (transparency)
    standCoverV->SetLineColor(TColor::GetColor("#ffcccc"));
    standCoverV->SetTransparency(0);

    TGeoCombiTrans *standCoverV_transf[2];

    standCoverV_transf[0] = new TGeoCombiTrans();
    standCoverV_transf[0]->SetDx(-(-standCover_XSize*0.5 + 373.145/*shift*/ + upperSideWall_XSize));
    standCoverV_transf[0]->SetDy(-41.5/*shift*/ + standCover_YSize*0.5);
    standCoverV_transf[0]->SetDz(0.0);

    standCoverV_transf[1] = new TGeoCombiTrans();
    standCoverV_transf[1]->SetDx(+(-standCover_XSize*0.5 + 373.145/*shift*/ + upperSideWall_XSize));
    standCoverV_transf[1]->SetDy(-41.5/*shift*/ + standCover_YSize*0.5);
    standCoverV_transf[1]->SetDz(0.0);

    supportStand->AddNode(standCoverV, 0, standCoverV_transf[0]);
    supportStand->AddNode(standCoverV, 1, standCoverV_transf[1]);

    //--------------------------------------------------------------------------

    //longitudinal stiffening rib ----------------------------------------------
    Double_t longitudinalStiffeningRib_XSize = 669.0; //cm
    Double_t longitudinalStiffeningRib_YSize = 45.0; //cm
    Double_t longitudinalStiffeningRib_ZSize = 5.0; //cm

    TGeoShape *longitudinalStiffeningRibS = new TGeoBBox(TString("longitudinalStiffeningRibS")+=TString("_") + supportStand->GetName(), longitudinalStiffeningRib_XSize*0.5, longitudinalStiffeningRib_YSize*0.5, longitudinalStiffeningRib_ZSize*0.5);

    TGeoVolume *longitudinalStiffeningRibV = new TGeoVolume(TString("longitudinalStiffeningRibV")+=TString("_") + supportStand->GetName(), longitudinalStiffeningRibS);

    //volume medium
    TGeoMedium *longitudinalStiffeningRibV_medium = pMedIron;
    if(longitudinalStiffeningRibV_medium) {
        longitudinalStiffeningRibV->SetMedium(longitudinalStiffeningRibV_medium);
    }
    else Fatal("Main", "Invalid medium for longitudinalStiffeningRibV!");

    //volume visual property (transparency)
    longitudinalStiffeningRibV->SetLineColor(TColor::GetColor("#ffcccc"));
    longitudinalStiffeningRibV->SetTransparency(0);

    TGeoCombiTrans *longitudinalStiffeningRibV_transf[1];

    longitudinalStiffeningRibV_transf[0] = new TGeoCombiTrans();
    longitudinalStiffeningRibV_transf[0]->SetDx(0.0);
    longitudinalStiffeningRibV_transf[0]->SetDy(longitudinalStiffeningRib_YSize*0.5 - 425.5/*shift*/);
    longitudinalStiffeningRibV_transf[0]->SetDz(0.0);

    supportStand->AddNode(longitudinalStiffeningRibV, 0, longitudinalStiffeningRibV_transf[0]);
    //--------------------------------------------------------------------------

    //transverse stiffening rib ------------------------------------------------
    Double_t transverseStiffeningRib_XSize = 78.5; //cm
    Double_t transverseStiffeningRib_YSize = 90.0; //cm
    Double_t transverseStiffeningRib_ZSize = 5.0; //cm

    Double_t transverseStiffeningRibNotch_XSize = 6.0; //cm
    Double_t transverseStiffeningRibNotch_YSize = 45.5; //cm
    Double_t transverseStiffeningRibNotch_ZSize = 5.0; //cm

    TGeoShape *transverseStiffeningRibS = new TGeoXtru(2);
    transverseStiffeningRibS->SetName(TString("transverseStiffeningRibS")+=TString("_") + supportStand->GetName());
    {
        Double_t x_pos[8] = {
            -transverseStiffeningRibNotch_XSize*0.5,
            -transverseStiffeningRib_XSize*0.5,
            -transverseStiffeningRib_XSize*0.5,
            +transverseStiffeningRib_XSize*0.5,
            +transverseStiffeningRib_XSize*0.5,
            +transverseStiffeningRibNotch_XSize*0.5,
            +transverseStiffeningRibNotch_XSize*0.5,
            -transverseStiffeningRibNotch_XSize*0.5
        };
        Double_t y_pos[8] = {
            -transverseStiffeningRib_YSize*0.5,
            -transverseStiffeningRib_YSize*0.5,
            +transverseStiffeningRib_YSize*0.5,
            +transverseStiffeningRib_YSize*0.5,
            -transverseStiffeningRib_YSize*0.5,
            -transverseStiffeningRib_YSize*0.5,
            +(transverseStiffeningRibNotch_YSize - transverseStiffeningRib_YSize*0.5),
            +(transverseStiffeningRibNotch_YSize - transverseStiffeningRib_YSize*0.5)
        };
        ((TGeoXtru*)transverseStiffeningRibS)->DefinePolygon(8, x_pos, y_pos);
        ((TGeoXtru*)transverseStiffeningRibS)->DefineSection(0, -transverseStiffeningRib_ZSize*0.5);
        ((TGeoXtru*)transverseStiffeningRibS)->DefineSection(1, +transverseStiffeningRib_ZSize*0.5);
    }

    TGeoVolume *transverseStiffeningRibV = new TGeoVolume(TString("transverseStiffeningRibV")+=TString("_") + supportStand->GetName(), transverseStiffeningRibS);

    //volume medium
    TGeoMedium *transverseStiffeningRibV_medium = pMedIron;
    if(transverseStiffeningRibV_medium) {
        transverseStiffeningRibV->SetMedium(transverseStiffeningRibV_medium);
    }
    else Fatal("Main", "Invalid medium for transverseStiffeningRibV!");

    //volume visual property (transparency)
    transverseStiffeningRibV->SetLineColor(TColor::GetColor("#ffcccc"));
    transverseStiffeningRibV->SetTransparency(0);

    TGeoCombiTrans *transverseStiffeningRibV_transf[7];

    transverseStiffeningRibV_transf[0] = new TGeoCombiTrans();
    transverseStiffeningRibV_transf[0]->RotateY(90.0);
    transverseStiffeningRibV_transf[0]->SetDx(0.0);
    transverseStiffeningRibV_transf[0]->SetDy(transverseStiffeningRib_YSize*0.5 - 425.5/*shift*/);
    transverseStiffeningRibV_transf[0]->SetDz(0.0);

    transverseStiffeningRibV_transf[1] = new TGeoCombiTrans();
    transverseStiffeningRibV_transf[1]->RotateY(90.0);
    transverseStiffeningRibV_transf[1]->SetDx(-160.5/*shift*/);
    transverseStiffeningRibV_transf[1]->SetDy(transverseStiffeningRib_YSize*0.5 - 425.5/*shift*/);
    transverseStiffeningRibV_transf[1]->SetDz(0.0);

    transverseStiffeningRibV_transf[2] = new TGeoCombiTrans();
    transverseStiffeningRibV_transf[2]->RotateY(90.0);
    transverseStiffeningRibV_transf[2]->SetDx(+160.5/*shift*/);
    transverseStiffeningRibV_transf[2]->SetDy(transverseStiffeningRib_YSize*0.5 - 425.5/*shift*/);
    transverseStiffeningRibV_transf[2]->SetDz(0.0);

    transverseStiffeningRibV_transf[3] = new TGeoCombiTrans();
    transverseStiffeningRibV_transf[3]->RotateY(90.0);
    transverseStiffeningRibV_transf[3]->SetDx(-(160.5 + 62.5) /*shift*/);
    transverseStiffeningRibV_transf[3]->SetDy(transverseStiffeningRib_YSize*0.5 - 425.5/*shift*/);
    transverseStiffeningRibV_transf[3]->SetDz(0.0);

    transverseStiffeningRibV_transf[4] = new TGeoCombiTrans();
    transverseStiffeningRibV_transf[4]->RotateY(90.0);
    transverseStiffeningRibV_transf[4]->SetDx(+(160.5 + 62.5) /*shift*/);
    transverseStiffeningRibV_transf[4]->SetDy(transverseStiffeningRib_YSize*0.5 - 425.5/*shift*/);
    transverseStiffeningRibV_transf[4]->SetDz(0.0);

    transverseStiffeningRibV_transf[5] = new TGeoCombiTrans();
    transverseStiffeningRibV_transf[5]->RotateY(90.0);
    transverseStiffeningRibV_transf[5]->SetDx(-(160.5 + 62.5*2) /*shift*/);
    transverseStiffeningRibV_transf[5]->SetDy(transverseStiffeningRib_YSize*0.5 - 425.5/*shift*/);
    transverseStiffeningRibV_transf[5]->SetDz(0.0);

    transverseStiffeningRibV_transf[6] = new TGeoCombiTrans();
    transverseStiffeningRibV_transf[6]->RotateY(90.0);
    transverseStiffeningRibV_transf[6]->SetDx(+(160.5 + 62.5*2) /*shift*/);
    transverseStiffeningRibV_transf[6]->SetDy(transverseStiffeningRib_YSize*0.5 - 425.5/*shift*/);
    transverseStiffeningRibV_transf[6]->SetDz(0.0);

    supportStand->AddNode(transverseStiffeningRibV, 0, transverseStiffeningRibV_transf[0]);
    supportStand->AddNode(transverseStiffeningRibV, 1, transverseStiffeningRibV_transf[1]);
    supportStand->AddNode(transverseStiffeningRibV, 2, transverseStiffeningRibV_transf[2]);
    supportStand->AddNode(transverseStiffeningRibV, 3, transverseStiffeningRibV_transf[3]);
    supportStand->AddNode(transverseStiffeningRibV, 4, transverseStiffeningRibV_transf[4]);
    supportStand->AddNode(transverseStiffeningRibV, 5, transverseStiffeningRibV_transf[5]);
    supportStand->AddNode(transverseStiffeningRibV, 6, transverseStiffeningRibV_transf[6]);
    //--------------------------------------------------------------------------

    //big inclined stiffening rib ----------------------------------------------
    Double_t bigInclinedStiffeningRib_XSize = 5.0; //cm
    Double_t bigInclinedStiffeningRib_YSize = 100.0; //cm
    Double_t bigInclinedStiffeningRib_ZSize = 79.0; //cm

    TGeoShape *bigInclinedStiffeningRibS = new TGeoBBox(TString("bigInclinedStiffeningRibS")+=TString("_") + supportStand->GetName(), bigInclinedStiffeningRib_XSize*0.5, bigInclinedStiffeningRib_YSize*0.5, bigInclinedStiffeningRib_ZSize*0.5);

    TGeoVolume *bigInclinedStiffeningRibV = new TGeoVolume(TString("bigInclinedStiffeningRibV")+=TString("_") + supportStand->GetName(), bigInclinedStiffeningRibS);

    //volume medium
    TGeoMedium *bigInclinedStiffeningRibV_medium = pMedIron;
    if(bigInclinedStiffeningRibV_medium) {
        bigInclinedStiffeningRibV->SetMedium(bigInclinedStiffeningRibV_medium);
    }
    else Fatal("Main", "Invalid medium for bigInclinedStiffeningRibV!");

    //volume visual property (transparency)
    bigInclinedStiffeningRibV->SetLineColor(TColor::GetColor("#ffcccc"));
    bigInclinedStiffeningRibV->SetTransparency(0);

    TGeoCombiTrans *bigInclinedStiffeningRibV_transf[2];

    bigInclinedStiffeningRibV_transf[0] = new TGeoCombiTrans();
    bigInclinedStiffeningRibV_transf[0]->SetDx(0.0);
    bigInclinedStiffeningRibV_transf[0]->SetDy(bigInclinedStiffeningRib_YSize*0.5 - 425.5/*shift*/ - 9.3/*shift*/);
    bigInclinedStiffeningRibV_transf[0]->SetDz(0.0);
    bigInclinedStiffeningRibV_transf[0]->RotateZ(-39.0);

    bigInclinedStiffeningRibV_transf[1] = new TGeoCombiTrans();
    bigInclinedStiffeningRibV_transf[1]->SetDx(0.0);
    bigInclinedStiffeningRibV_transf[1]->SetDy(bigInclinedStiffeningRib_YSize*0.5 - 425.5/*shift*/ - 9.3/*shift*/);
    bigInclinedStiffeningRibV_transf[1]->SetDz(0.0);
    bigInclinedStiffeningRibV_transf[1]->RotateZ(+39.0);

    supportStand->AddNode(bigInclinedStiffeningRibV, 0, bigInclinedStiffeningRibV_transf[0]);
    supportStand->AddNode(bigInclinedStiffeningRibV, 1, bigInclinedStiffeningRibV_transf[1]);
    //--------------------------------------------------------------------------

    //small inclined stiffening rib --------------------------------------------
    Double_t smallInclinedStiffeningRib_XSize = 5.0; //cm
    Double_t smallInclinedStiffeningRib_YSize = 39.0; //cm
    Double_t smallInclinedStiffeningRib_ZSize = 79.0; //cm

    TGeoShape *smallInclinedStiffeningRibS = new TGeoBBox(TString("smallInclinedStiffeningRibS")+=TString("_") + supportStand->GetName(), smallInclinedStiffeningRib_XSize*0.5, smallInclinedStiffeningRib_YSize*0.5, smallInclinedStiffeningRib_ZSize*0.5);

    TGeoVolume *smallInclinedStiffeningRibV = new TGeoVolume(TString("smallInclinedStiffeningRibV")+=TString("_") + supportStand->GetName(), smallInclinedStiffeningRibS);

    //volume medium
    TGeoMedium *smallInclinedStiffeningRibV_medium = pMedIron;
    if(smallInclinedStiffeningRibV_medium) {
        smallInclinedStiffeningRibV->SetMedium(smallInclinedStiffeningRibV_medium);
    }
    else Fatal("Main", "Invalid medium for smallInclinedStiffeningRibV!");

    //volume visual property (transparency)
    smallInclinedStiffeningRibV->SetLineColor(TColor::GetColor("#ffcccc"));
    smallInclinedStiffeningRibV->SetTransparency(0);

    TGeoCombiTrans *smallInclinedStiffeningRibV_transf[2];

    smallInclinedStiffeningRibV_transf[0] = new TGeoCombiTrans();
    smallInclinedStiffeningRibV_transf[0]->SetDx(0.0);
    smallInclinedStiffeningRibV_transf[0]->SetDy(-smallInclinedStiffeningRib_YSize*0.5 - 425.5/*shift*/ + 90.7/*shift*/);
    smallInclinedStiffeningRibV_transf[0]->SetDz(0.0);
    smallInclinedStiffeningRibV_transf[0]->RotateZ(-77.58);

    smallInclinedStiffeningRibV_transf[1] = new TGeoCombiTrans();
    smallInclinedStiffeningRibV_transf[1]->SetDx(0.0);
    smallInclinedStiffeningRibV_transf[1]->SetDy(-smallInclinedStiffeningRib_YSize*0.5 - 425.5/*shift*/ + 90.7/*shift*/);
    smallInclinedStiffeningRibV_transf[1]->SetDz(0.0);
    smallInclinedStiffeningRibV_transf[1]->RotateZ(+77.58);

    supportStand->AddNode(smallInclinedStiffeningRibV, 0, smallInclinedStiffeningRibV_transf[0]);
    supportStand->AddNode(smallInclinedStiffeningRibV, 1, smallInclinedStiffeningRibV_transf[1]);
    //--------------------------------------------------------------------------

    //towing lug ---------------------------------------------------------------
    Double_t towingLug_ZSize = 5.0; //cm
    Double_t towingLug_ZDistBetweenParts = 14.0;

    Double_t towingLugHole_RMax = 6.25; //cm
    Double_t towingLugHole_ZSize = 5.0 + 0.5; //cm

    TGeoShape *towingLugS = new TGeoXtru(2);
    towingLugS->SetName(TString("towingLugS")+=TString("_") + supportStand->GetName());
    {
        Double_t x_pos[11] = {
            +5.0, //0
            +5.0, //1
            +32.0, //2
            +50.0, //3
            +62.5, //4
            +62.5, //5
            +45.9, //6
            +0.0, //7
            +0.0, //8
            +38.6, //9
            +38.6 //10
        };
        Double_t y_pos[11] = {
            -25.0, //0
            -29.0, //1
            -44.6, //2
            -44.6, //3
            -38.5, //4
            -24.6, //5
            +0.0, //6
            +0.0, //7
            -15.0, //8
            -15.0, //9
            -25.0 //10
        };
        ((TGeoXtru*)towingLugS)->DefinePolygon(11, x_pos, y_pos);
        ((TGeoXtru*)towingLugS)->DefineSection(0, -towingLug_ZSize*0.5);
        ((TGeoXtru*)towingLugS)->DefineSection(1, +towingLug_ZSize*0.5);
    }

    TGeoShape *towingLugHoleS= new TGeoTube(TString("towingLugHoleS")+=TString("_") + supportStand->GetName(), 0.0, towingLugHole_RMax, towingLugHole_ZSize*0.5);

    TGeoTranslation *towingLugS_pos = new TGeoTranslation();
        towingLugS_pos->SetName(TString("towingLugS_pos")+=TString("_") + supportStand->GetName());
        towingLugS_pos->SetDx(0.0);
        towingLugS_pos->SetDy(0.0);
        towingLugS_pos->SetDz(0.0);
        towingLugS_pos->RegisterYourself();

    TGeoTranslation *towingLugHoleS_pos = new TGeoTranslation();
        towingLugHoleS_pos->SetName(TString("towingLugHoleS_pos")+=TString("_") + supportStand->GetName());
        towingLugHoleS_pos->SetDx(+50.0);
        towingLugHoleS_pos->SetDy(-32.1);
        towingLugHoleS_pos->SetDz(0.0);
        towingLugHoleS_pos->RegisterYourself();

    TGeoCompositeShape *towingLugWithHoleS = new TGeoCompositeShape();
    towingLugWithHoleS->SetName(TString("towingLugWithHoleS")+=TString("_") + supportStand->GetName());
    {
        TString expression = "towingLugS"; expression += TString("_") + supportStand->GetName();
            expression += ":towingLugS_pos"; expression += TString("_") + supportStand->GetName();
            expression += "-towingLugHoleS"; expression += TString("_") + supportStand->GetName();
            expression += ":towingLugHoleS_pos"; expression += TString("_") + supportStand->GetName();

        towingLugWithHoleS->MakeNode(expression);
        towingLugWithHoleS->ComputeBBox(); //need to compute a bounding box
    }

    TGeoVolume *towingLugV = new TGeoVolume(TString("towingLugV")+=TString("_") + supportStand->GetName(), towingLugS);

    //volume medium
    TGeoMedium *towingLugV_medium = pMedIron;
    if(towingLugV_medium) {
        towingLugV->SetMedium(towingLugV_medium);
    }
    else Fatal("Main", "Invalid medium for towingLugV!");

    //volume visual property (transparency)
    towingLugV->SetLineColor(TColor::GetColor("#ffcccc"));
    towingLugV->SetTransparency(0);

    TGeoCombiTrans *towingLugV_transf[2];

    towingLugV_transf[0] = new TGeoCombiTrans();
    towingLugV_transf[0]->SetDx(+340.0);
    towingLugV_transf[0]->SetDy(-425.5 + 15.0 /*shift*/);
    towingLugV_transf[0]->SetDz(+(towingLug_ZDistBetweenParts*0.5 + towingLug_ZSize*0.5));

    towingLugV_transf[1] = new TGeoCombiTrans();
    towingLugV_transf[1]->SetDx(+340.0);
    towingLugV_transf[1]->SetDy(-425.5 + 15.0 /*shift*/);
    towingLugV_transf[1]->SetDz(-(towingLug_ZDistBetweenParts*0.5 + towingLug_ZSize*0.5));

    supportStand->AddNode(towingLugV, 0, towingLugV_transf[0]);
    supportStand->AddNode(towingLugV, 1, towingLugV_transf[1]);
    //--------------------------------------------------------------------------

    //brackets with spacers ----------------------------------------------------
    Double_t spacer_XSize = 24.0; //cm
    Double_t spacer_YSize = 23.0; //cm
    Double_t spacer_ZSize = 1.5; //cm

    Double_t spacer_XYChamfer = 2.0; //cm

    Double_t bracket_XSize = 20.0; //cm

    Double_t bracketVericalPart_YSize = 25.0; //cm
    Double_t bracketVericalPart_ZSize = 4.0; //cm

    Double_t bracketHorizontalPart_YSize = 4.0; //cm
    Double_t bracketHorizontalPart_ZSize = 21.0; //cm

    Double_t bracketTrianglePart_SideLength = 11.0; //cm
    Double_t bracketTrianglePart_Thickness = 4.0; //cm

    TGeoShape *spacerS = new TGeoXtru(2);
    spacerS->SetName(TString("spacerS")+=TString("_") + supportStand->GetName());
    {
        Double_t x_pos[8] = {
            +(spacer_XSize*0.5 - spacer_XYChamfer),
            +spacer_XSize*0.5,
            +spacer_XSize*0.5,
            +(spacer_XSize*0.5 - spacer_XYChamfer),
            -(spacer_XSize*0.5 - spacer_XYChamfer),
            -spacer_XSize*0.5,
            -spacer_XSize*0.5,
            -(spacer_XSize*0.5 - spacer_XYChamfer)
        };
        Double_t y_pos[8] = {
            +spacer_YSize*0.5,
            +(spacer_YSize*0.5 - spacer_XYChamfer),
            -(spacer_YSize*0.5 - spacer_XYChamfer),
            -spacer_YSize*0.5,
            -spacer_YSize*0.5,
            -(spacer_YSize*0.5 - spacer_XYChamfer),
            +(spacer_YSize*0.5 - spacer_XYChamfer),
            +spacer_YSize*0.5
        };
        ((TGeoXtru*)spacerS)->DefinePolygon(8, x_pos, y_pos);
        ((TGeoXtru*)spacerS)->DefineSection(0, -spacer_ZSize*0.5);
        ((TGeoXtru*)spacerS)->DefineSection(1, +spacer_ZSize*0.5);
    }

    TGeoShape *bracketVericalPartS = new TGeoBBox(TString("bracketVericalPartS")+=TString("_") + supportStand->GetName(), bracket_XSize*0.5, bracketVericalPart_YSize*0.5, bracketVericalPart_ZSize*0.5);
    TGeoShape *bracketHorizontalPartS = new TGeoBBox(TString("bracketHorizontalPartS")+=TString("_") + supportStand->GetName(), bracket_XSize*0.5, bracketHorizontalPart_YSize*0.5, bracketHorizontalPart_ZSize*0.5);
    TGeoShape *bracketTrianglePartS = new TGeoTrd1(TString("bracketTrianglePartS")+=TString("_") + supportStand->GetName(), 0.0, (TMath::Sqrt(bracketTrianglePart_SideLength*bracketTrianglePart_SideLength + bracketTrianglePart_SideLength*bracketTrianglePart_SideLength))*0.5, bracketTrianglePart_Thickness*0.5, (bracketTrianglePart_SideLength*TMath::Sin(45.0*TMath::DegToRad()))*0.5);

    TGeoTranslation *bracketVericalPart_pos = new TGeoTranslation();
        bracketVericalPart_pos->SetName(TString("bracketVericalPart_pos")+=TString("_") + supportStand->GetName());
        bracketVericalPart_pos->SetDx(0.0);
        bracketVericalPart_pos->SetDy(0.0);
        bracketVericalPart_pos->SetDz(0.0);
        bracketVericalPart_pos->RegisterYourself();

    TGeoTranslation *bracketHorizontalPart_pos = new TGeoTranslation();
        bracketHorizontalPart_pos->SetName(TString("bracketHorizontalPart_pos")+=TString("_") + supportStand->GetName());
        bracketHorizontalPart_pos->SetDx(0.0);
        bracketHorizontalPart_pos->SetDy(+(bracketVericalPart_YSize*0.5 - bracketHorizontalPart_YSize*0.5));
        bracketHorizontalPart_pos->SetDz(-(bracketHorizontalPart_ZSize*0.5 + bracketVericalPart_ZSize*0.5));
        bracketHorizontalPart_pos->RegisterYourself();

    TGeoCombiTrans *bracketTrianglePartS_pos = new TGeoCombiTrans();
        bracketTrianglePartS_pos->SetName(TString("bracketTrianglePartS_pos")+=TString("_") + supportStand->GetName());
        bracketTrianglePartS_pos->RotateZ(90.0);
        bracketTrianglePartS_pos->RotateY(180.0);
        bracketTrianglePartS_pos->RotateX(-45.0);
        bracketTrianglePartS_pos->SetDx(0.0);
        bracketTrianglePartS_pos->SetDy(+5.6/*shift*/);
        bracketTrianglePartS_pos->SetDz(-4.8/*shift*/);
        bracketTrianglePartS_pos->RegisterYourself();

    TGeoCompositeShape *bracketS = new TGeoCompositeShape();
    bracketS->SetName(TString("bracketS")+=TString("_") + supportStand->GetName());
    {
        TString expression = "bracketVericalPartS"; expression += TString("_") + supportStand->GetName();
            expression += ":bracketVericalPart_pos"; expression += TString("_") + supportStand->GetName();
            expression += "+bracketHorizontalPartS"; expression += TString("_") + supportStand->GetName();
            expression += ":bracketHorizontalPart_pos"; expression += TString("_") + supportStand->GetName();
            expression += "+bracketTrianglePartS"; expression += TString("_") + supportStand->GetName();
            expression += ":bracketTrianglePartS_pos"; expression += TString("_") + supportStand->GetName();

        bracketS->MakeNode(expression);
        bracketS->ComputeBBox(); //need to compute a bounding box
    }

    TGeoVolume *spacerV = new TGeoVolume(TString("spacerV")+=TString("_") + supportStand->GetName(), spacerS);
    TGeoVolume *bracketV = new TGeoVolume(TString("bracketV")+=TString("_") + supportStand->GetName(), bracketS);

    //volume medium
    TGeoMedium *spacerV_medium = pMedIron;
    if(spacerV_medium) {
        spacerV->SetMedium(spacerV_medium);
        bracketV->SetMedium(spacerV_medium);
    }
    else Fatal("Main", "Invalid medium for spacerV and bracketV!");

    //volume visual property (transparency)
    spacerV->SetLineColor(TColor::GetColor("#ffcccc"));
    spacerV->SetTransparency(0);
    bracketV->SetLineColor(TColor::GetColor("#ffcccc"));
    bracketV->SetTransparency(0);

    TGeoVolume *bracketsWithSpacersA = new TGeoVolumeAssembly(TString("bracketsWithSpacersA")+=TString("_") + supportStand->GetName());
    bracketsWithSpacersA->SetMedium(pMedAir);

    /*
    Double_t curr_angle = -90.0 + 9.31; //deg

    for(Int_t i = 0; i < 26; ++i) {

        TGeoCombiTrans *spacerV_transf[26];

        spacerV_transf[i] = new TGeoCombiTrans();
        spacerV_transf[i]->SetDx(0.0);
        spacerV_transf[i]->SetDy(-349.32 + 6.0);
        spacerV_transf[i]->SetDz(-(supportBase_ZSize*0.5 + spacer_ZSize*0.5));
        spacerV_transf[i]->RotateZ(curr_angle);

        TGeoCombiTrans *bracketV_transf[26];

        bracketV_transf[i] = new TGeoCombiTrans();
        bracketV_transf[i]->SetDx(0.0);
        bracketV_transf[i]->SetDy(-349.32 + 6.0);
        bracketV_transf[i]->SetDz(-(bracketVericalPart_ZSize*0.5 + spacer_ZSize*0.5) -(supportBase_ZSize*0.5 + spacer_ZSize*0.5));
        bracketV_transf[i]->RotateZ(curr_angle);

        curr_angle += (i % 2 == 0) ? 7.0903 : 5.7689;

        bracketsWithSpacersA->AddNode(spacerV, i, spacerV_transf[i]);
        bracketsWithSpacersA->AddNode(bracketV, i, bracketV_transf[i]);
    }

    TGeoCombiTrans *bracketsWithSpacersA_transf[2];

    bracketsWithSpacersA_transf[0] = new TGeoCombiTrans();
    bracketsWithSpacersA_transf[0]->SetDx(0.0);
    bracketsWithSpacersA_transf[0]->SetDy(0.0);
    bracketsWithSpacersA_transf[0]->SetDz(0.0);

    bracketsWithSpacersA_transf[1] = new TGeoCombiTrans();
    bracketsWithSpacersA_transf[1]->RotateY(180.0);
    bracketsWithSpacersA_transf[1]->SetDx(0.0);
    bracketsWithSpacersA_transf[1]->SetDy(0.0);
    bracketsWithSpacersA_transf[1]->SetDz(0.0);

    supportStand->AddNode(bracketsWithSpacersA, 0, bracketsWithSpacersA_transf[0]);
    supportStand->AddNode(bracketsWithSpacersA, 1, bracketsWithSpacersA_transf[1]);
    */
    Double_t curr_angle = -90.0 + 9.31; //deg

    for(Int_t i = 0; i < 26; ++i) {

        TGeoCombiTrans *spacerV_sideA_transf[26];
        spacerV_sideA_transf[i] = new TGeoCombiTrans();
        spacerV_sideA_transf[i]->SetDx(0.0);
        spacerV_sideA_transf[i]->SetDy(-349.32 + 6.0);
        spacerV_sideA_transf[i]->SetDz(-(supportBase_ZSize*0.5 + spacer_ZSize*0.5));
        spacerV_sideA_transf[i]->RotateZ(curr_angle);

        TGeoCombiTrans *spacerV_sideB_transf[26];
        spacerV_sideB_transf[i] = new TGeoCombiTrans(*spacerV_sideA_transf[i]);
        spacerV_sideB_transf[i]->ReflectZ(true);

        TGeoCombiTrans *bracketV_sideA_transf[26];
        bracketV_sideA_transf[i] = new TGeoCombiTrans();
        bracketV_sideA_transf[i]->SetDx(0.0);
        bracketV_sideA_transf[i]->SetDy(-349.32 + 6.0);
        bracketV_sideA_transf[i]->SetDz(-(bracketVericalPart_ZSize*0.5 + spacer_ZSize*0.5) -(supportBase_ZSize*0.5 + spacer_ZSize*0.5));
        bracketV_sideA_transf[i]->RotateZ(curr_angle);

        TGeoCombiTrans *bracketV_sideB_transf[26];
        bracketV_sideB_transf[i] = new TGeoCombiTrans(*bracketV_sideA_transf[i]);
        bracketV_sideB_transf[i]->ReflectZ(true);

        curr_angle += (i % 2 == 0) ? 7.0903 : 5.7689;

        supportStand->AddNode(spacerV, i, spacerV_sideA_transf[i]);
        supportStand->AddNode(bracketV, i, bracketV_sideA_transf[i]);

        supportStand->AddNode(spacerV, 26+i, spacerV_sideB_transf[i]);
        supportStand->AddNode(bracketV,26+i, bracketV_sideB_transf[i]);
    }

    //--------------------------------------------------------------------------

    //foot ---------------------------------------------------------------------
    Double_t foot_ZSize = 40.0; //cm
    Double_t footFiller_ZSize = 24.0; //cm
    Double_t footStand_ZSize = 40.0; //cm

    TGeoShape *footHalfS = new TGeoXtru(2);
    footHalfS->SetName(TString("footHalfS")+=TString("_") + supportStand->GetName());
    {
        Double_t x_pos[11] = {
            +0.0,
            +18.0,
            +18.0,
            +56.6,
            +56.6,
            +51.915,
            +37.75,
            +37.75,
            +47.75,
            +47.75,
            +0.0
        };
        Double_t y_pos[11] = {
            -6.5,
            -6.5,
            -31.5,
            -31.5,
            -29.5,
            -16.0,
            -16.0,
            -4.0,
            -4.0,
            +0.0,
            0.0
        };
        ((TGeoXtru*)footHalfS)->DefinePolygon(11, x_pos, y_pos);
        ((TGeoXtru*)footHalfS)->DefineSection(0, -foot_ZSize*0.5);
        ((TGeoXtru*)footHalfS)->DefineSection(1, +foot_ZSize*0.5);
    }

    TGeoShape *footFillerS = new TGeoXtru(2);
    footFillerS->SetName(TString("footFillerS")+=TString("_") + supportStand->GetName());
    {
        Double_t x_pos[4] = {
            +37.75,
            +51.915,
            +47.75,
            +37.75
        };
        Double_t y_pos[4] = {
            -16.0,
            -16.0,
            -4.0,
            -4.0
        };
        ((TGeoXtru*)footFillerS)->DefinePolygon(4, x_pos, y_pos);
        ((TGeoXtru*)footFillerS)->DefineSection(0, -footFiller_ZSize*0.5);
        ((TGeoXtru*)footFillerS)->DefineSection(1, +footFiller_ZSize*0.5);
    }

    TGeoShape *footStandS = new TGeoXtru(2);
    footStandS->SetName(TString("footStandS")+=TString("_") + supportStand->GetName());
    {
        Double_t x_pos[4] = {
            +33.1,
            +56.6,
            +56.6,
            +33.1,
        };
        Double_t y_pos[4] = {
            -41.1,
            -41.1,
            -31.5,
            -31.5
        };
        ((TGeoXtru*)footStandS)->DefinePolygon(4, x_pos, y_pos);
        ((TGeoXtru*)footStandS)->DefineSection(0, -footStand_ZSize*0.5);
        ((TGeoXtru*)footStandS)->DefineSection(1, +footStand_ZSize*0.5);
    }

    TGeoCombiTrans *footHalf_pos = new TGeoCombiTrans();
        footHalf_pos->SetName(TString("footHalf_pos")+=TString("_") + supportStand->GetName());
        footHalf_pos->RotateY(180.0);
        footHalf_pos->SetDx(0.0);
        footHalf_pos->SetDy(0.0);
        footHalf_pos->SetDz(0.0);
        footHalf_pos->RegisterYourself();

    TGeoCombiTrans *footFiller_pos = new TGeoCombiTrans();
        footFiller_pos->SetName(TString("footFiller_pos")+=TString("_") + supportStand->GetName());
        footFiller_pos->RotateY(180.0);
        footFiller_pos->SetDx(0.0);
        footFiller_pos->SetDy(0.0);
        footFiller_pos->SetDz(0.0);
        footFiller_pos->RegisterYourself();

    TGeoCombiTrans *footStand_pos = new TGeoCombiTrans();
        footStand_pos->SetName(TString("footStand_pos")+=TString("_") + supportStand->GetName());
        footStand_pos->RotateY(180.0);
        footStand_pos->SetDx(0.0);
        footStand_pos->SetDy(0.0);
        footStand_pos->SetDz(0.0);
        footStand_pos->RegisterYourself();

    TGeoCompositeShape *footS = new TGeoCompositeShape();
    footS->SetName(TString("footS")+=TString("_") + supportStand->GetName());
    {
        TString expression = "footHalfS"; expression += TString("_") + supportStand->GetName();
            expression += "+footFillerS"; expression += TString("_") + supportStand->GetName();
            expression += "+footHalfS"; expression += TString("_") + supportStand->GetName();
            expression += ":footHalf_pos"; expression += TString("_") + supportStand->GetName();
            expression += "+footFillerS"; expression += TString("_") + supportStand->GetName();
            expression += ":footFiller_pos"; expression += TString("_") + supportStand->GetName();
            expression += "+footStandS"; expression += TString("_") + supportStand->GetName();
            expression += "+footStandS"; expression += TString("_") + supportStand->GetName();
            expression += ":footStand_pos"; expression += TString("_") + supportStand->GetName();

        footS->MakeNode(expression);
        footS->ComputeBBox(); //need to compute a bounding box
    }

    TGeoVolume *footV = new TGeoVolume(TString("footV")+=TString("_") + supportStand->GetName(), footS);

    //volume medium
    TGeoMedium *footV_medium = pMedIron;
    if(footV_medium) {
        footV->SetMedium(footV_medium);
    }
    else Fatal("Main", "Invalid medium for footV!");

    //volume visual property (transparency)
    footV->SetLineColor(TColor::GetColor("#ffcccc"));
    footV->SetTransparency(0);

    TGeoCombiTrans *footV_transf[3];

    footV_transf[0] = new TGeoCombiTrans();
    footV_transf[0]->RotateY(90.0);
    footV_transf[0]->SetDx(0.0);
    footV_transf[0]->SetDy(-(425.5 + 10.0)/*shift*/);
    footV_transf[0]->SetDz(0.0);

    footV_transf[1] = new TGeoCombiTrans();
    footV_transf[1]->RotateY(90.0);
    footV_transf[1]->SetDx(-(supportBase_XSize*0.5 - foot_ZSize*0.5));
    footV_transf[1]->SetDy(-(425.5 + 10.0)/*shift*/);
    footV_transf[1]->SetDz(0.0);

    footV_transf[2] = new TGeoCombiTrans();
    footV_transf[2]->RotateY(90.0);
    footV_transf[2]->SetDx(+(supportBase_XSize*0.5 - foot_ZSize*0.5));
    footV_transf[2]->SetDy(-(425.5 + 10.0)/*shift*/);
    footV_transf[2]->SetDz(0.0);

    supportStand->AddNode(footV, 0, footV_transf[0]);
    supportStand->AddNode(footV, 1, footV_transf[1]);
    supportStand->AddNode(footV, 2, footV_transf[2]);
    //--------------------------------------------------------------------------

    return supportStand;
}
//------------------------------------------------------------------------------

TGeoVolume *CreateHorizontalBeam(TString horizontal_beam_name) {

    //horizontal beam volume
    TGeoVolume *horizontalBeam = new TGeoVolumeAssembly(horizontal_beam_name);
    horizontalBeam->SetMedium(pMedAir);

    //tube ---------------------------------------------------------------------
    Double_t tube_XYSize = 30.0; //cm
    Double_t tube_ZSize = 499.0; //cm
    Double_t tube_WallThickness = 1.2; //cm
    Double_t tube_RMax = 3.6;
    Double_t tube_RMin = tube_RMax - tube_WallThickness;
    Double_t tube_wallPartSize = tube_XYSize*0.5 - tube_RMax;

    TGeoShape *tubeWallSegmentCylS = new TGeoTubeSeg(TString("tubeWallSegmentCylS")+=TString("_") + horizontalBeam->GetName(), tube_RMin, tube_RMax, tube_ZSize*0.5, 0.0, 90.0);
    TGeoShape *tubeWallSegmentBoxS = new TGeoBBox(TString("tubeWallSegmentBoxS")+=TString("_") + horizontalBeam->GetName(), tube_wallPartSize*0.5, tube_WallThickness*0.5, tube_ZSize*0.5);

    TGeoCombiTrans *tubeWallSegmentCyl_pos1 = new TGeoCombiTrans();
        tubeWallSegmentCyl_pos1->SetName(TString("tubeWallSegmentCyl_pos1")+=TString("_") + horizontalBeam->GetName());
        tubeWallSegmentCyl_pos1->SetDx(+tube_wallPartSize);
        tubeWallSegmentCyl_pos1->SetDy(+tube_wallPartSize);
        tubeWallSegmentCyl_pos1->SetDz(0.0);
        tubeWallSegmentCyl_pos1->RegisterYourself();

    TGeoCombiTrans *tubeWallSegmentBoxS_pos1 = new TGeoCombiTrans();
        tubeWallSegmentBoxS_pos1->SetName(TString("tubeWallSegmentBoxS_pos1")+=TString("_") + horizontalBeam->GetName());
        tubeWallSegmentBoxS_pos1->SetDx(+tube_wallPartSize*0.5);
        tubeWallSegmentBoxS_pos1->SetDy(+tube_XYSize*0.5 - tube_WallThickness*0.5);
        tubeWallSegmentBoxS_pos1->SetDz(0.0);
        tubeWallSegmentBoxS_pos1->RegisterYourself();

    TGeoCombiTrans *tubeWallSegmentBoxS_pos2 = new TGeoCombiTrans();
        tubeWallSegmentBoxS_pos2->SetName(TString("tubeWallSegmentBoxS_pos2")+=TString("_") + horizontalBeam->GetName());
        tubeWallSegmentBoxS_pos2->RotateZ(90.0);
        tubeWallSegmentBoxS_pos2->SetDx(+tube_XYSize*0.5 - tube_WallThickness*0.5);
        tubeWallSegmentBoxS_pos2->SetDy(+tube_wallPartSize*0.5);
        tubeWallSegmentBoxS_pos2->SetDz(0.0);
        tubeWallSegmentBoxS_pos2->RegisterYourself();

    TGeoCompositeShape *tubeQuadrantS = new TGeoCompositeShape();
    tubeQuadrantS->SetName(TString("tubeQuadrantS")+=TString("_") + horizontalBeam->GetName());
    {
        TString expression = "tubeWallSegmentCylS"; expression += TString("_") + horizontalBeam->GetName();
            expression += ":tubeWallSegmentCyl_pos1"; expression += TString("_") + horizontalBeam->GetName();
            expression += "+tubeWallSegmentBoxS"; expression += TString("_") + horizontalBeam->GetName();
            expression += ":tubeWallSegmentBoxS_pos1"; expression += TString("_") + horizontalBeam->GetName();
            expression += "+tubeWallSegmentBoxS"; expression += TString("_") + horizontalBeam->GetName();
            expression += ":tubeWallSegmentBoxS_pos2"; expression += TString("_") + horizontalBeam->GetName();

        tubeQuadrantS->MakeNode(expression);
        tubeQuadrantS->ComputeBBox(); //need to compute a bounding box
    }

    TGeoVolume *tubeQuadrantV = new TGeoVolume(TString("tubeQuadrantV")+=TString("_") + horizontalBeam->GetName(), tubeQuadrantS);

    //volume medium
    TGeoMedium *tubeQuadrantV_medium = pMedIron;
    if(tubeQuadrantV_medium) {
        tubeQuadrantV->SetMedium(tubeQuadrantV_medium);
    }
    else Fatal("Main", "Invalid medium for tubeQuadrantV!");

    //volume visual property (transparency)
    tubeQuadrantV->SetLineColor(TColor::GetColor("#ffcccc"));
    tubeQuadrantV->SetTransparency(0);

    TGeoCombiTrans *tubeQuadrantV_transf[4];

    tubeQuadrantV_transf[0] = new TGeoCombiTrans();
    tubeQuadrantV_transf[0]->SetDx(0.0);
    tubeQuadrantV_transf[0]->SetDy(0.0);
    tubeQuadrantV_transf[0]->SetDz(0.0);

    tubeQuadrantV_transf[1] = new TGeoCombiTrans();
    tubeQuadrantV_transf[1]->RotateZ(90.0);
    tubeQuadrantV_transf[1]->SetDx(0.0);
    tubeQuadrantV_transf[1]->SetDy(0.0);
    tubeQuadrantV_transf[1]->SetDz(0.0);

    tubeQuadrantV_transf[2] = new TGeoCombiTrans();
    tubeQuadrantV_transf[2]->RotateZ(180.0);
    tubeQuadrantV_transf[2]->SetDx(0.0);
    tubeQuadrantV_transf[2]->SetDy(0.0);
    tubeQuadrantV_transf[2]->SetDz(0.0);

    tubeQuadrantV_transf[3] = new TGeoCombiTrans();
    tubeQuadrantV_transf[3]->RotateZ(270.0);
    tubeQuadrantV_transf[3]->SetDx(0.0);
    tubeQuadrantV_transf[3]->SetDy(0.0);
    tubeQuadrantV_transf[3]->SetDz(0.0);

    horizontalBeam->AddNode(tubeQuadrantV, 0, tubeQuadrantV_transf[0]);
    horizontalBeam->AddNode(tubeQuadrantV, 1, tubeQuadrantV_transf[1]);
    horizontalBeam->AddNode(tubeQuadrantV, 2, tubeQuadrantV_transf[2]);
    horizontalBeam->AddNode(tubeQuadrantV, 3, tubeQuadrantV_transf[3]);
    //--------------------------------------------------------------------------

    //flange -------------------------------------------------------------------
    Double_t flange_XYSize = 60.0; //cm
    Double_t flange_ZSize = 3.0; //cm
    Double_t flange_RMax = 6.0;
    Double_t flange_RMin = 0.0;
    Double_t flange_PartSize = flange_XYSize*0.5 - flange_RMax;

    TGeoShape *flangeSegmentCylS = new TGeoTubeSeg(TString("flangeSegmentCylS")+=TString("_") + horizontalBeam->GetName(), flange_RMin, flange_RMax, flange_ZSize*0.5, 0.0, 90.0);
    TGeoShape *flangeOutSegmentBoxS = new TGeoBBox(TString("flangeOutSegmentBoxS")+=TString("_") + horizontalBeam->GetName(), flange_PartSize*0.5, flange_RMax*0.5, flange_ZSize*0.5);
    TGeoShape *flangeInSegmentBoxS = new TGeoBBox(TString("flangeInSegmentBoxS")+=TString("_") + horizontalBeam->GetName(), flange_PartSize*0.5, flange_PartSize*0.5, flange_ZSize*0.5);

    TGeoCombiTrans *flangeSegmentCylS_pos1 = new TGeoCombiTrans();
        flangeSegmentCylS_pos1->SetName(TString("flangeSegmentCylS_pos1")+=TString("_") + horizontalBeam->GetName());
        flangeSegmentCylS_pos1->SetDx(+flange_PartSize);
        flangeSegmentCylS_pos1->SetDy(+flange_PartSize);
        flangeSegmentCylS_pos1->SetDz(0.0);
        flangeSegmentCylS_pos1->RegisterYourself();

    TGeoCombiTrans *flangeOutSegmentBoxS_pos1 = new TGeoCombiTrans();
        flangeOutSegmentBoxS_pos1->SetName(TString("flangeOutSegmentBoxS_pos1")+=TString("_") + horizontalBeam->GetName());
        flangeOutSegmentBoxS_pos1->SetDx(+flange_PartSize*0.5);
        flangeOutSegmentBoxS_pos1->SetDy(+flange_XYSize*0.5 - flange_RMax*0.5);
        flangeOutSegmentBoxS_pos1->SetDz(0.0);
        flangeOutSegmentBoxS_pos1->RegisterYourself();

    TGeoCombiTrans *flangeOutSegmentBoxS_pos2 = new TGeoCombiTrans();
        flangeOutSegmentBoxS_pos2->SetName(TString("flangeOutSegmentBoxS_pos2")+=TString("_") + horizontalBeam->GetName());
        flangeOutSegmentBoxS_pos2->RotateZ(90.0);
        flangeOutSegmentBoxS_pos2->SetDx(+flange_XYSize*0.5 - flange_RMax*0.5);
        flangeOutSegmentBoxS_pos2->SetDy(+flange_PartSize*0.5);
        flangeOutSegmentBoxS_pos2->SetDz(0.0);
        flangeOutSegmentBoxS_pos2->RegisterYourself();

    TGeoCombiTrans *flangeInSegmentBoxS_pos1 = new TGeoCombiTrans();
        flangeInSegmentBoxS_pos1->SetName(TString("flangeInSegmentBoxS_pos1")+=TString("_") + horizontalBeam->GetName());
        flangeInSegmentBoxS_pos1->SetDx(+flange_PartSize*0.5);
        flangeInSegmentBoxS_pos1->SetDy(+flange_PartSize*0.5);
        flangeInSegmentBoxS_pos1->SetDz(0.0);
        flangeInSegmentBoxS_pos1->RegisterYourself();

    TGeoCompositeShape *flangeQuadrantS = new TGeoCompositeShape();
    flangeQuadrantS->SetName(TString("flangeQuadrantS")+=TString("_") + horizontalBeam->GetName());
    {
        TString expression = "flangeSegmentCylS"; expression += TString("_") + horizontalBeam->GetName();
            expression += ":flangeSegmentCylS_pos1"; expression += TString("_") + horizontalBeam->GetName();
            expression += "+flangeOutSegmentBoxS"; expression += TString("_") + horizontalBeam->GetName();
            expression += ":flangeOutSegmentBoxS_pos1"; expression += TString("_") + horizontalBeam->GetName();
            expression += "+flangeOutSegmentBoxS"; expression += TString("_") + horizontalBeam->GetName();
            expression += ":flangeOutSegmentBoxS_pos2"; expression += TString("_") + horizontalBeam->GetName();
            expression += "+flangeInSegmentBoxS"; expression += TString("_") + horizontalBeam->GetName();
            expression += ":flangeInSegmentBoxS_pos1"; expression += TString("_") + horizontalBeam->GetName();

        flangeQuadrantS->MakeNode(expression);
        flangeQuadrantS->ComputeBBox(); //need to compute a bounding box
    }

    TGeoVolume *flangeQuadrantV = new TGeoVolume(TString("flangeQuadrantV")+=TString("_") + horizontalBeam->GetName(), flangeQuadrantS);

    //volume medium
    TGeoMedium *flangeQuadrantV_medium = pMedIron;
    if(flangeQuadrantV_medium) {
        flangeQuadrantV->SetMedium(flangeQuadrantV_medium);
    }
    else Fatal("Main", "Invalid medium for flangeQuadrantV!");

    //volume visual property (transparency)
    flangeQuadrantV->SetLineColor(TColor::GetColor("#ffcccc"));
    flangeQuadrantV->SetTransparency(0);

    TGeoCombiTrans *flangeQuadrantV_transf[8];

    flangeQuadrantV_transf[0] = new TGeoCombiTrans();
    flangeQuadrantV_transf[0]->SetDx(0.0);
    flangeQuadrantV_transf[0]->SetDy(0.0);
    flangeQuadrantV_transf[0]->SetDz(-(tube_ZSize*0.5+flange_ZSize*0.5));

    flangeQuadrantV_transf[1] = new TGeoCombiTrans();
    flangeQuadrantV_transf[1]->RotateZ(90.0);
    flangeQuadrantV_transf[1]->SetDx(0.0);
    flangeQuadrantV_transf[1]->SetDy(0.0);
    flangeQuadrantV_transf[1]->SetDz(-(tube_ZSize*0.5+flange_ZSize*0.5));

    flangeQuadrantV_transf[2] = new TGeoCombiTrans();
    flangeQuadrantV_transf[2]->RotateZ(180.0);
    flangeQuadrantV_transf[2]->SetDx(0.0);
    flangeQuadrantV_transf[2]->SetDy(0.0);
    flangeQuadrantV_transf[2]->SetDz(-(tube_ZSize*0.5+flange_ZSize*0.5));

    flangeQuadrantV_transf[3] = new TGeoCombiTrans();
    flangeQuadrantV_transf[3]->RotateZ(270.0);
    flangeQuadrantV_transf[3]->SetDx(0.0);
    flangeQuadrantV_transf[3]->SetDy(0.0);
    flangeQuadrantV_transf[3]->SetDz(-(tube_ZSize*0.5+flange_ZSize*0.5));

    flangeQuadrantV_transf[4] = new TGeoCombiTrans();
    flangeQuadrantV_transf[4]->SetDx(0.0);
    flangeQuadrantV_transf[4]->SetDy(0.0);
    flangeQuadrantV_transf[4]->SetDz(+(tube_ZSize*0.5+flange_ZSize*0.5));

    flangeQuadrantV_transf[5] = new TGeoCombiTrans();
    flangeQuadrantV_transf[5]->RotateZ(90.0);
    flangeQuadrantV_transf[5]->SetDx(0.0);
    flangeQuadrantV_transf[5]->SetDy(0.0);
    flangeQuadrantV_transf[5]->SetDz(+(tube_ZSize*0.5+flange_ZSize*0.5));

    flangeQuadrantV_transf[6] = new TGeoCombiTrans();
    flangeQuadrantV_transf[6]->RotateZ(180.0);
    flangeQuadrantV_transf[6]->SetDx(0.0);
    flangeQuadrantV_transf[6]->SetDy(0.0);
    flangeQuadrantV_transf[6]->SetDz(+(tube_ZSize*0.5+flange_ZSize*0.5));

    flangeQuadrantV_transf[7] = new TGeoCombiTrans();
    flangeQuadrantV_transf[7]->RotateZ(270.0);
    flangeQuadrantV_transf[7]->SetDx(0.0);
    flangeQuadrantV_transf[7]->SetDy(0.0);
    flangeQuadrantV_transf[7]->SetDz(+(tube_ZSize*0.5+flange_ZSize*0.5));

    horizontalBeam->AddNode(flangeQuadrantV, 0, flangeQuadrantV_transf[0]);
    horizontalBeam->AddNode(flangeQuadrantV, 1, flangeQuadrantV_transf[1]);
    horizontalBeam->AddNode(flangeQuadrantV, 2, flangeQuadrantV_transf[2]);
    horizontalBeam->AddNode(flangeQuadrantV, 3, flangeQuadrantV_transf[3]);
    horizontalBeam->AddNode(flangeQuadrantV, 4, flangeQuadrantV_transf[4]);
    horizontalBeam->AddNode(flangeQuadrantV, 5, flangeQuadrantV_transf[5]);
    horizontalBeam->AddNode(flangeQuadrantV, 6, flangeQuadrantV_transf[6]);
    horizontalBeam->AddNode(flangeQuadrantV, 7, flangeQuadrantV_transf[7]);
    //--------------------------------------------------------------------------

    //flange rib ---------------------------------------------------------------
    Double_t flangeRib_ZSize = 3.0; //cm

    TGeoShape *flangeRibS = new TGeoXtru(2);
    flangeRibS->SetName(TString("flangeRibS")+=TString("_") + horizontalBeam->GetName());
    {
        Double_t x_pos[6] = {
            +0.0,
            -2.5,
            -15.0,
            -15.0,
            -2.5,
            +0.0
        };
        Double_t y_pos[6] = {
            +15.0,
            +15.0,
            +2.5,
            +0.0,
            +0.0,
            +2.5
        };
        ((TGeoXtru*)flangeRibS)->DefinePolygon(6, x_pos, y_pos);
        ((TGeoXtru*)flangeRibS)->DefineSection(0, -flangeRib_ZSize*0.5);
        ((TGeoXtru*)flangeRibS)->DefineSection(1, +flangeRib_ZSize*0.5);
    }

    TGeoVolume *flangeRibV = new TGeoVolume(TString("flangeRibV")+=TString("_") + horizontalBeam->GetName(), flangeRibS);

    //volume medium
    TGeoMedium *flangeRibV_medium = pMedIron;
    if(flangeRibV_medium) {
        flangeRibV->SetMedium(flangeRibV_medium);
    }
    else Fatal("Main", "Invalid medium for flangeRibV!");

    //volume visual property (transparency)
    flangeRibV->SetLineColor(TColor::GetColor("#ffcccc"));
    flangeRibV->SetTransparency(0);

    TGeoCombiTrans *flangeRibV_transf[8];

    flangeRibV_transf[0] = new TGeoCombiTrans();
    flangeRibV_transf[0]->RotateY(90.0);
    flangeRibV_transf[0]->SetDx(0.0);
    flangeRibV_transf[0]->SetDy(+tube_XYSize*0.5);
    flangeRibV_transf[0]->SetDz(-tube_ZSize*0.5);

    flangeRibV_transf[1] = new TGeoCombiTrans();
    flangeRibV_transf[1]->RotateY(90.0);
    flangeRibV_transf[1]->SetDx(0.0);
    flangeRibV_transf[1]->SetDy(+tube_XYSize*0.5);
    flangeRibV_transf[1]->SetDz(-tube_ZSize*0.5);
    flangeRibV_transf[1]->RotateZ(90.0);

    flangeRibV_transf[2] = new TGeoCombiTrans();
    flangeRibV_transf[2]->RotateY(90.0);
    flangeRibV_transf[2]->SetDx(0.0);
    flangeRibV_transf[2]->SetDy(+tube_XYSize*0.5);
    flangeRibV_transf[2]->SetDz(-tube_ZSize*0.5);
    flangeRibV_transf[2]->RotateZ(180.0);

    flangeRibV_transf[3] = new TGeoCombiTrans();
    flangeRibV_transf[3]->RotateY(90.0);
    flangeRibV_transf[3]->SetDx(0.0);
    flangeRibV_transf[3]->SetDy(+tube_XYSize*0.5);
    flangeRibV_transf[3]->SetDz(-tube_ZSize*0.5);
    flangeRibV_transf[3]->RotateZ(270.0);

    flangeRibV_transf[4] = new TGeoCombiTrans();
    flangeRibV_transf[4]->RotateY(-90.0);
    flangeRibV_transf[4]->SetDx(0.0);
    flangeRibV_transf[4]->SetDy(+tube_XYSize*0.5);
    flangeRibV_transf[4]->SetDz(+tube_ZSize*0.5);

    flangeRibV_transf[5] = new TGeoCombiTrans();
    flangeRibV_transf[5]->RotateY(-90.0);
    flangeRibV_transf[5]->SetDx(0.0);
    flangeRibV_transf[5]->SetDy(+tube_XYSize*0.5);
    flangeRibV_transf[5]->SetDz(+tube_ZSize*0.5);
    flangeRibV_transf[5]->RotateZ(90.0);

    flangeRibV_transf[6] = new TGeoCombiTrans();
    flangeRibV_transf[6]->RotateY(-90.0);
    flangeRibV_transf[6]->SetDx(0.0);
    flangeRibV_transf[6]->SetDy(+tube_XYSize*0.5);
    flangeRibV_transf[6]->SetDz(+tube_ZSize*0.5);
    flangeRibV_transf[6]->RotateZ(180.0);

    flangeRibV_transf[7] = new TGeoCombiTrans();
    flangeRibV_transf[7]->RotateY(-90.0);
    flangeRibV_transf[7]->SetDx(0.0);
    flangeRibV_transf[7]->SetDy(+tube_XYSize*0.5);
    flangeRibV_transf[7]->SetDz(+tube_ZSize*0.5);
    flangeRibV_transf[7]->RotateZ(270.0);

    horizontalBeam->AddNode(flangeRibV, 0, flangeRibV_transf[0]);
    horizontalBeam->AddNode(flangeRibV, 1, flangeRibV_transf[1]);
    horizontalBeam->AddNode(flangeRibV, 2, flangeRibV_transf[2]);
    horizontalBeam->AddNode(flangeRibV, 3, flangeRibV_transf[3]);
    horizontalBeam->AddNode(flangeRibV, 4, flangeRibV_transf[4]);
    horizontalBeam->AddNode(flangeRibV, 5, flangeRibV_transf[5]);
    horizontalBeam->AddNode(flangeRibV, 6, flangeRibV_transf[6]);
    horizontalBeam->AddNode(flangeRibV, 7, flangeRibV_transf[7]);
    //--------------------------------------------------------------------------

    return horizontalBeam;
}

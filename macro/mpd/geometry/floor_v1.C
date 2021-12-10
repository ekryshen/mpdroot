// Geometry for pipe: aluminium + beryllium
// Demezhan Myktybekov  myktybekov@jinr.ru
// July 2020

#include <TGeoMatrix.h>
#include <TGeoMedium.h>
#include <TGeoTube.h>


//---------------------------
#ifndef __CLING__
#include <FairGeoBuilder.h>
#include <FairGeoInterface.h>
#include <FairGeoLoader.h>
#include <FairGeoMedia.h>
#include <Rtypes.h>
#include <RtypesCore.h>
#include <TError.h>
#include <TFile.h>
#include <TGLRnrCtx.h>
#include <TGLViewer.h>
#include <TGLViewerBase.h>
#include <TGeoBBox.h>
#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <TString.h>
#include <TSystem.h>
#endif

// Pipe's position
const Double_t ZMin      = -455.0;
const Double_t ZMax      = 455.0;
const Double_t YMin      = -480.0;
const Double_t YMax      = 405;
const Double_t XMin      = -405;
const Double_t XMax      = 380;
const Double_t Width     = 0.1;
const Double_t HalfWidth = Width * 0.5;
const Double_t Dz        = (ZMax - ZMin) * 0.5;
const Double_t Dy        = (YMax - YMin) * 0.5;
const Double_t Dx        = (XMax - XMin) * 0.5;
const Double_t Yshift    = (YMax + YMin) * 0.5;


// media
TGeoMedium* pMedAir       = 0;
TGeoMedium* pMedAluminium = 0;
TGeoMedium* pMedVacuum    = 0;

class FairGeoMedia;
class FairGeoBuilder;


void DefineRequiredMedia(FairGeoMedia* geoMedia, FairGeoBuilder* geoBuild);

//__________________________________________________________________________
void floor_v1() {

  // ----  set working directory  --------------------------------------------
  TString gPath = gSystem->Getenv("VMCWORKDIR");

  // -------   Module file name (output)   ---------------------------------
  const TString geoModuleName    = "flor";
  const TString geoModuleVersion = "v1";
  TString geoFileName            = gPath + "/geometry/" + geoModuleName + "_" + geoModuleVersion + ".root";
  // -------------------------------------------------------------------------

  // ----  global geometry parameters  ---------------------------------------
  FairGeoLoader* geoLoad    = new FairGeoLoader("TGeo", "FairGeoLoader");
  FairGeoInterface* geoFace = geoLoad->getGeoInterface();

  // -------   Load media from media file   ----------------------------------
  TString medFile = gPath + "/geometry/media.geo";
  geoFace->setMediaFile(medFile);
  geoFace->readMedia();
  FairGeoMedia* geoMedia   = geoFace->getMedia();
  FairGeoBuilder* geoBuild = geoLoad->getGeoBuilder();

  DefineRequiredMedia(geoMedia, geoBuild);
  // -------------------------------------------------------------------------

  // --------------   Create geometry and global top volume  ------------------------
  gGeoManager = (TGeoManager*) gROOT->FindObject("FAIRGeom");
  gGeoManager->SetName(geoModuleName + "_geom");
  TGeoVolume* top = new TGeoVolumeAssembly("TOP");
  top->SetMedium(pMedVacuum);
  gGeoManager->SetTopVolume(top);


  TGeoVolume* CRADLE = new TGeoVolumeAssembly(geoModuleName);
  CRADLE->SetMedium(pMedVacuum);
  // --------------- Add geometry stuff -------------------------------------

  // podloga
  {
    TGeoShape* wall_1 = new TGeoBBox(Dx + 2 * HalfWidth, HalfWidth, Dz + 2 * HalfWidth);
    // wall_1->SetName("floor");
    TGeoVolume* Wall1 = new TGeoVolume("fl01", wall_1, pMedVacuum);
    //  Wall1->SetTransparency(0);
    TGeoCombiTrans* emtpyTrans = new TGeoCombiTrans();
    emtpyTrans->SetDy(YMin - HalfWidth);
    CRADLE->AddNode(Wall1, 0, emtpyTrans);
  }
  // sufit
  {
    TGeoShape* wall_6 = new TGeoBBox(Dx + 2 * HalfWidth, HalfWidth, Dz + 2 * HalfWidth);
    // wall_6->SetName("ceil");
    TGeoVolume* Wall6 = new TGeoVolume("fl06", wall_6, pMedVacuum);
    // Wall6->SetTransparency(0);
    TGeoCombiTrans* emtpyTrans6 = new TGeoCombiTrans();
    emtpyTrans6->SetDy(YMax + HalfWidth);
    CRADLE->AddNode(Wall6, 0, emtpyTrans6);
  }
  // lewa sciana
  {
    TGeoShape* wall_2      = new TGeoBBox(Dx + HalfWidth * 2, Dy, HalfWidth);
    TGeoVolume* Wall2      = new TGeoVolume("fl02", wall_2, pMedVacuum);
    TGeoCombiTrans* trans2 = new TGeoCombiTrans();
    trans2->SetDz(-Dz - HalfWidth);
    trans2->SetDy(Yshift);
    CRADLE->AddNode(Wall2, 0, trans2);
  }
  {
    TGeoShape* wall_4      = new TGeoBBox(Dx + HalfWidth * 2, Dy, HalfWidth);
    TGeoVolume* Wall4      = new TGeoVolume("fl04", wall_4, pMedVacuum);
    TGeoCombiTrans* trans4 = new TGeoCombiTrans();
    trans4->SetDz(Dz + HalfWidth);
    trans4->SetDy(Yshift);
    CRADLE->AddNode(Wall4, 0, trans4);
  }
  // boczne
  {
    TGeoShape* wall_3      = new TGeoBBox(HalfWidth, Dy, Dz);
    TGeoVolume* Wall3      = new TGeoVolume("fl03", wall_3, pMedVacuum);
    TGeoCombiTrans* trans3 = new TGeoCombiTrans();
    trans3->SetDx(Dx + HalfWidth);
    trans3->SetDy(Yshift);
    CRADLE->AddNode(Wall3, 0, trans3);
  }

  {
    TGeoShape* wall_5      = new TGeoBBox(HalfWidth, Dy, Dz);
    TGeoVolume* Wall5      = new TGeoVolume("fl05", wall_5, pMedVacuum);
    TGeoCombiTrans* trans5 = new TGeoCombiTrans();
    trans5->SetDx(-(Dx + HalfWidth));
    trans5->SetDy(Yshift);
    CRADLE->AddNode(Wall5, 0, trans5);
  }
  // now led add stuff around mcord
  {  // upper part
    TGeoShape* mcordUp     = new TGeoTubeSeg(350, 350.1, 451, 0, 180);
    TGeoVolume* mcordUpV   = new TGeoVolume("fl06", mcordUp, pMedVacuum);
    TGeoCombiTrans* trans7 = new TGeoCombiTrans();
    // trans7->RotateY(90);
    CRADLE->AddNode(mcordUpV, 0, trans7);
  }
  {  // mcord down
    TGeoShape* mcordUp     = new TGeoTubeSeg(350, 350.1, 240, 180, 360);
    TGeoVolume* mcordUpV   = new TGeoVolume("fl07", mcordUp, pMedVacuum);
    TGeoCombiTrans* trans7 = new TGeoCombiTrans();
    // trans7->RotateY(90);
    CRADLE->AddNode(mcordUpV, 0, trans7);
  }

  // ---------------   Finish   ----------------------------------------------

  top->AddNode(CRADLE, 0);
  top->SetVisContainers(kTRUE);

  // ---------------   Finish   -----------------------------------------------

  gGeoManager->CloseGeometry();
  gGeoManager->CheckOverlaps(0.001);
  gGeoManager->PrintOverlaps();

  gGeoManager->Test();

  gGeoManager->SetMaxVisNodes(100000);

  TFile* geoFile = new TFile(geoFileName, "RECREATE");
  top->Write();
  geoFile->Close();

  // if (wrGeoWithMaterials)
  //{
  //    TString geoFile_wMat = gPath + "/geometry/" + geoDetectorName + "_"+ geoDetectorVersion + "_with_materials.root";
  //    gGeoManager->Export(geoFile_wMat);
  // TString geoFile_gdml = gPath + "/geometry/" + geoDetectorName + "_"+ geoDetectorVersion + ".gdml";
  // gGeoManager->Export(geoFile_gdml);
  //}

  top->Draw("ogl");
  TGLViewer* v = (TGLViewer*) gPad->GetViewer3D();
  v->SetStyle(TGLRnrCtx::kOutline);
}

//__________________________________________________________________________
void DefineRequiredMedia(FairGeoMedia* geoMedia, FairGeoBuilder* geoBuild) {

  // aluminium medium
  FairGeoMedium* mAluminium = geoMedia->getMedium("aluminium");
  if (!mAluminium) Fatal("Main", "FairMedium aluminium not found");
  geoBuild->createMedium(mAluminium);
  pMedAluminium = gGeoManager->GetMedium("aluminium");
  if (!pMedAluminium) Fatal("Main", "Medium aluminium not found");


  FairGeoMedium* mVacuum = geoMedia->getMedium("activevacuum");
  if (!mVacuum) Fatal("Main", "FairMedium vacuum not found");
  geoBuild->createMedium(mVacuum);
  pMedVacuum = gGeoManager->GetMedium("activevacuum");
  if (!pMedVacuum) Fatal("Main", "Medium Vacuum not found");
}

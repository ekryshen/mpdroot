/*
 * mcord_v4.C
 *
 *  Created on: 9 gru 2021
 *      Author: Daniel Wielanek
 *      E-mail: daniel.wielanek at pw.edu.pl
 *      Warsaw University of Technology, Faculty of Physics
 *      JINR,  Laboratory of High Energy Physics
 */

/********************************************************************************
 *    Copyright (C) 2014 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    *
 *                                                                              *
 *              This software is distributed under the terms of the             *
 *              GNU Lesser General Public Licence (LGPL) version 3,             *
 *                  copied verbatim in the file "LICENSE"                       *
 ********************************************************************************/
// in root all sizes are given in cm


#ifndef __CLING__
#include "FairGeoBuilder.h"
#include "FairGeoInterface.h"
#include "FairGeoLoader.h"
#include "FairGeoMedia.h"
#include "TFile.h"
#include "TGeoCompositeShape.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMatrix.h"
#include "TGeoMedium.h"
#include "TGeoPgon.h"
#include "TGeoTube.h"
#include "TGeoVolume.h"
#include "TGeoXtru.h"
#include "TList.h"
#include "TROOT.h"
#include "TString.h"
#include "TSystem.h"
#include <TMath.h>
#include <iostream>
using namespace std;
#endif

// Name of geometry version and output file
const TString geoVersion = "mcord";
TString FileName         = geoVersion + ".root";
// const TString FileName1  = geoVersion + "_geomanager.root";

// Names of the different used materials which are used to build the modules
// The materials are defined in the global media.geo file
const TString KeepingVolumeMedium = "air";
const TString BoxVolumeMedium     = "silicon";

const TString PvtMedium = "gold";

// Distance of the center of the first detector layer [cm];
const Float_t First_Z_Position = 10;
const Float_t Z_Distance       = 10;

// Silicon box for both module types
const Float_t Module_Size_X = 80.;
const Float_t Module_Size_Y = 80.;
const Float_t Module_Size_Z = .04;

// some global variables
TGeoManager* gGeoMan = NULL;  // Pointer to TGeoManager instance
TGeoVolume* gModules;         // Global storage for module types

// Forward declarations
void create_materials_from_media_file();
TGeoVolume* create_detector();
void position_detector();
void add_alignable_volumes();
TGeoVolume* GetSensor();
TGeoVolume* GetModule();
TGeoVolume* GetDetector();

TGeoVolume* GetSupport();
TGeoVolume* CreateLadder();
TGeoVolume* GetSector();
void DEBUG(int i) { std::cout << "******************** " << i << "  **************************" << endl; }
void mcord_v4() {
  // Load the necessary FairRoot libraries
  //  gROOT->LoadMacro("$VMCWORKDIR/gconfig/basiclibs.C");
  //  basiclibs();
  //  gSystem->Load("libGeoBase");
  //  gSystem->Load("libParBase");
  //  gSystem->Load("libBase");
  const TString geoModuleName    = "mcord";
  const TString geoModuleVersion = "v4";
  TString gPath                  = gSystem->Getenv("VMCWORKDIR");
  FileName                       = gPath + "/geometry/" + geoModuleName + "_" + geoModuleVersion + ".root";
  create_materials_from_media_file();
  gGeoMan = (TGeoManager*) gROOT->FindObject("FAIRGeom");
  gGeoMan->SetVisLevel(7);

  // Create the top volume

  TGeoVolume* top = new TGeoVolumeAssembly("TOP");
  gGeoMan->SetTopVolume(top);

  TGeoVolume* tut4 = new TGeoVolumeAssembly(geoVersion);
  top->AddNode(tut4, 1);
  DEBUG(3);
  gModules = GetDetector();
  DEBUG(4);
  // position_detector();
  DEBUG(5);
  cout << "Voxelizing." << endl;
  top->Voxelize("");
  gGeoMan->CloseGeometry();
  DEBUG(6);
  //  add_alignable_volumes();
  DEBUG(7);
  gGeoMan->CheckOverlaps(0.001);
  gGeoMan->PrintOverlaps();
  gGeoMan->Test();
  DEBUG(8);
  TFile* outfile = TFile::Open(FileName, "RECREATE");
  top->Write();
  outfile->Close();

  top->Draw("ogl");
  TGLViewer* v = (TGLViewer*) gPad->GetViewer3D();
  v->SetStyle(TGLRnrCtx::kOutline);
  /*
    TFile* outfile1 = TFile::Open(FileName1, "RECREATE");
    gGeoMan->Write();
    outfile1->Close();
  */

  // ------------------------------------------------------------------------
}

void create_materials_from_media_file() {
  // Use the FairRoot geometry interface to load the media which are already defined
  FairGeoLoader* geoLoad    = new FairGeoLoader("TGeo", "FairGeoLoader");
  FairGeoInterface* geoFace = geoLoad->getGeoInterface();
  TString geoPath           = gSystem->Getenv("VMCWORKDIR");
  TString geoFile           = geoPath + "/geometry/media.geo";
  geoFace->setMediaFile(geoFile);
  geoFace->readMedia();

  // Read the required media and create them in the GeoManager
  FairGeoMedia* geoMedia   = geoFace->getMedia();
  FairGeoBuilder* geoBuild = geoLoad->getGeoBuilder();

  FairGeoMedium* air       = geoMedia->getMedium("air");
  FairGeoMedium* aluminium = geoMedia->getMedium("aluminium");
  FairGeoMedium* pvt       = geoMedia->getMedium(PvtMedium);


  // include check if all media are found

  geoBuild->createMedium(air);
  geoBuild->createMedium(aluminium);
  geoBuild->createMedium(pvt);
}

TGeoVolume* GetDetector() {
  TGeoMedium* air = gGeoMan->GetMedium("air");
  Double_t Rmin   = 329.2;
  // TGeoTube* tube  = new TGeoTube(Rmin, Rmin + 20.0, 239.2);

  TGeoVolumeAssembly* as = new TGeoVolumeAssembly();
  TGeoVolume* vol        = gGeoMan->GetVolume(geoVersion);

  //    TGeoTube *tube = new TGeoTube(32, 4000, 3000);
  // TGeoVolume *vol = new TGeoVolume("md01",tube,air);
  // vol->SetLineColor(kGreen);
  //  vol->SetTransparency(70);
  Double_t step      = 360. / 28.;
  TGeoVolume* subVol = GetModule();
  Rmin += 7.0;
  for (int i = 0; i < 28; i++) {  // 28
    TGeoRotation* rot = new TGeoRotation();
    Double_t phi      = step * i;
    rot->RotateZ(phi);
    phi         = -phi;
    Double_t Ry = TMath::Cos(TMath::DegToRad() * phi) * Rmin;
    Double_t Rx = TMath::Sin(TMath::DegToRad() * phi) * Rmin;
    vol->AddNode(subVol, i, new TGeoCombiTrans(Rx, Ry, 0, rot));
  }
  return vol;
}


TGeoVolume* GetModule() {
  TGeoMedium* air = gGeoMan->GetMedium("air");
  // TGeoBBox* box   = new TGeoBBox("module", 0.5 * 73.5, 7.0, 478.4 * 0.5);
  // TGeoBBox* box           = new TGeoBBox("module", 0.5 * 73.5, 7.0, 474.4 * 0.5);
  TGeoVolumeAssembly* assem = new TGeoVolumeAssembly("module");
  // TGeoVolume* vol = new TGeoVolume("md01module", box, air);
  // vol->SetLineColor(kYellow);


  // TGeoVolume *sec = GetSector();

  TGeoVolume* support = GetSupport();
  TGeoVolume* ladder  = CreateLadder();

  support->SetLineColor(kRed);

  Double_t dZ = 150.;
  assem->AddNode(ladder, 0, new TGeoTranslation(0, 1.5, -dZ));
  assem->AddNode(ladder, 1, new TGeoTranslation(0, 5.5, 0));
  assem->AddNode(ladder, 2, new TGeoTranslation(0, 1.5, dZ));

  assem->AddNode(support, 0, new TGeoTranslation(0, -1., -dZ));
  assem->AddNode(support, 1, new TGeoTranslation(0, 3, 0));
  assem->AddNode(support, 2, new TGeoTranslation(0, -1., dZ));

  return assem;
}

TGeoVolume* GetSupport() {
  TGeoMedium* air      = gGeoMan->GetMedium("air");
  TGeoMedium* al       = gGeoMan->GetMedium("aluminium");
  TGeoBBox* frame      = new TGeoBBox("frame", 73.5 * 0.5, 1, 51.0);
  TGeoVolume* volFrame = new TGeoVolume("md01support", frame, air);

  Double_t dY = 0.5;
  Double_t dX = 1.5;
  Double_t Z  = 50.1;
  Double_t X  = 73.5 * 0.5;
  Double_t hX = X * 0.5;
  Double_t DY = 0.5;

  TGeoBBox* box_1 = new TGeoBBox("md01supv", dX, dY, Z);
  TGeoBBox* box_2 = new TGeoBBox("mdo1suph", X - 2.0 * dX, dY, dX);

  TGeoVolume* box_1vol = new TGeoVolume("md01supv", box_1, al);
  TGeoVolume* box_2vol = new TGeoVolume("md02supv", box_2, al);
  box_1vol->SetLineColor(kRed);
  box_2vol->SetLineColor(kCyan);

  volFrame->AddNode(box_1vol, 0, new TGeoTranslation(-dX + X, DY, 0));
  volFrame->AddNode(box_1vol, 1, new TGeoTranslation(dX - X, DY, 0));

  volFrame->AddNode(box_2vol, 0, new TGeoTranslation(0, DY, Z - dX));
  volFrame->AddNode(box_2vol, 1, new TGeoTranslation(0, DY, -Z + dX));

  TGeoXtru* tru = new TGeoXtru(2);
  Double_t LX   = X - dX * 2;
  Double_t LZ   = Z - dX * 2;
  Double_t U    = 3.0;
  Double_t x0   = 0;
  Double_t x1   = U;
  Double_t x2   = LX - U;
  Double_t x3   = LX;
  Double_t y0   = 0;
  Double_t y1   = U;
  Double_t y2   = LZ - U;
  Double_t y3   = LZ;


  Double_t x[16] = {x0, x2, x3, x3, x1, x3, x3, x2, x0, -x2, -x3, -x3, -x1, -x3, -x3, -x2};
  Double_t y[16] = {-y1, -y3, -y3, -y2, y0, y2, y3, y3, y1, y3, y3, y2, y0, -y2, -y3, -y3};

  Double_t fx[4] = {-10, -10, 10, 10};
  Double_t fy[4] = {-10, 10, 10, -11};

  tru->DefinePolygon(16, x, y);
  tru->DefineSection(0, 0, 0, 0, 1);
  tru->DefineSection(1, dY * 2, 0, 0, 1);


  TGeoVolume* truVol = new TGeoVolume("md01tru", tru, al);
  truVol->SetLineColor(kYellow);
  TGeoRotation* rot = new TGeoRotation();
  rot->RotateX(90);
  volFrame->AddNode(truVol, 0, new TGeoCombiTrans(0, 2 * dY, 0, rot));

  return volFrame;
}


TGeoVolume* GetSector() {
  TGeoMedium* air       = gGeoMan->GetMedium("air");
  TGeoMedium* aluminium = gGeoMan->GetMedium("aluminium");
  TGeoCombiTrans* com1  = new TGeoCombiTrans(0, -1.5, 0, new TGeoRotation());
  TGeoCombiTrans* com2  = new TGeoCombiTrans(0, 1.0, 0, new TGeoRotation());

  TGeoBBox* down          = new TGeoBBox("down", 73.5 * 0.5, 1.0, 51.0);
  TGeoBBox* up            = new TGeoBBox("up", 73.5 * 0.5, 1.5, 174.4 * 0.5);
  TGeoVolume* volDown     = new TGeoVolume("md01sup", down, air);
  TGeoVolume* volUp       = CreateLadder();
  TGeoVolumeAssembly* ass = new TGeoVolumeAssembly("md01sector");
  ass->AddNode(volDown, 0, com1);
  ass->AddNode(volUp, 0, com2);

  ass->SetLineColor(kYellow);
  return ass;
}


TGeoVolume* CreateSkeleton() { return nullptr; }

TGeoVolume* CreateLadder() {
  TGeoMedium* air = gGeoMan->GetMedium("air");
  Double_t dX     = 73.5 * 0.5;
  TGeoBBox* box   = new TGeoBBox("module", dX, 1.5, 174.4 * 0.5);
  TGeoVolume* vol = new TGeoVolume("md01ladder", box, air);
  vol->SetLineColor(kCyan);
  vol->SetTransparency(70);
  TGeoVolume* sensors = GetSensor();
  for (int i = 0; i < 8; i++) {
    double j = i;
    vol->AddNode(sensors, i, new TGeoTranslation(-73.5 * 0.5 + 7.0 + 8.5 * j, 0, 0));
  }
  return vol;
}

TGeoVolume* GetSensor() {
  TGeoMedium* air       = gGeoMan->GetMedium("air");
  TGeoMedium* vac       = gGeoMan->GetMedium("vacuum");
  TGeoMedium* alumunium = gGeoMan->GetMedium("aluminium");
  TGeoMedium* pvt       = gGeoMan->GetMedium(PvtMedium);
  Double_t width        = 8.0;
  Double_t dZ           = 174.4 * 0.5;
  Double_t dY           = 1.5;
  Double_t dX           = width * 0.5;
  Double_t scintDz      = 162.0 * 0.5;
  Double_t thick        = 0.30;
  Double_t dT           = thick * 0.5;
  Double_t airDz        = (dZ - scintDz - thick) * 0.5;

  Double_t height      = 3.0;
  TGeoBBox* entire_box = new TGeoBBox("sens", dX, dY, dZ);
  TGeoBBox* wall_down  = new TGeoBBox("wallA", dX, dT, dZ);
  TGeoBBox* wall_front = new TGeoBBox("wallB", dX - thick, dY - thick, dT);
  TGeoBBox* wall_up    = new TGeoBBox("wallC", dT, dY - thick, dZ);
  TGeoBBox* air_gap    = new TGeoBBox("airGap", dX - thick, dY - thick, airDz);
  TGeoBBox* scint      = new TGeoBBox("scitn", dX - thick, dY - thick, scintDz);


  TGeoVolume* wallDownVol = new TGeoVolume("md01wallDown", wall_down, alumunium);
  wallDownVol->SetLineColor(kGreen);
  TGeoVolume* wallFrontVol = new TGeoVolume("md01wallFront", wall_front, alumunium);
  wallDownVol->SetLineColor(kGreen);
  TGeoVolume* wallWallVol = new TGeoVolume("md01wallLeft", wall_up, pvt);
  wallWallVol->SetLineColor(kGreen + 2);
  TGeoVolume* airVol = new TGeoVolume("md01airVol", air_gap, air);
  airVol->SetLineColor(kBlue);
  TGeoVolume* scintVol = new TGeoVolume("md01scintVol", scint, pvt);
  scintVol->SetLineColor(kYellow);
  TGeoVolume* sensors = new TGeoVolume("md01sensor", entire_box, vac);

  sensors->AddNode(wallDownVol, 0, new TGeoCombiTrans(0, -dY + dT, 0, new TGeoRotation()));
  sensors->AddNode(wallDownVol, 1, new TGeoCombiTrans(0, dY - dT, 0, new TGeoRotation()));
  sensors->AddNode(wallFrontVol, 0, new TGeoCombiTrans(0, 0, -dZ + dT, new TGeoRotation()));
  sensors->AddNode(wallFrontVol, 1, new TGeoCombiTrans(0, 0, dZ - dT, new TGeoRotation()));
  sensors->AddNode(wallWallVol, 0, new TGeoCombiTrans(dX - dT, 0, 0, new TGeoRotation));
  sensors->AddNode(wallWallVol, 1, new TGeoCombiTrans(-dX + dT, 0, 0, new TGeoRotation));
  sensors->AddNode(airVol, 0, new TGeoCombiTrans(0, 0, -dZ + thick + airDz, new TGeoRotation));
  sensors->AddNode(airVol, 1, new TGeoCombiTrans(0, 0, dZ - thick - airDz, new TGeoRotation));
  sensors->AddNode(scintVol, 0, new TGeoCombiTrans(0, 0, 0, new TGeoRotation()));

  return sensors;
}

void position_detector() {

  gGeoMan->GetVolume(geoVersion)->AddNode(gModules, 0, 0);
  return;
  TGeoTranslation* det_trans = NULL;

  Int_t numDets = 0;
  for (Int_t detectorPlanes = 0; detectorPlanes < 40; detectorPlanes++) {
    det_trans = new TGeoTranslation("", 0., 0., First_Z_Position + (numDets * Z_Distance));
    gGeoMan->GetVolume(geoVersion)->AddNode(gModules, numDets, det_trans);
    numDets++;
  }
}

void add_alignable_volumes() {

  TString volPath;
  TString symName;
  TString detStr = "Tutorial4/det";
  TString volStr = "/TOP_1/tutorial4_1/tut4_det_";

  for (Int_t detectorPlanes = 0; detectorPlanes < 40; detectorPlanes++) {

    volPath = volStr;
    volPath += detectorPlanes;

    symName = detStr;
    symName += Form("%02d", detectorPlanes);

    cout << "Path: " << volPath << ", " << symName << endl;
    //    gGeoMan->cd(volPath);

    //   gGeoMan->SetAlignableEntry(symName.Data(),volPath.Data());
  }
  cout << "Nr of alignable objects: " << gGeoMan->GetNAlignable() << endl;
}


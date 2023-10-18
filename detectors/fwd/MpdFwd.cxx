//------------------------------------------------------------------------------------------------------------------------
/// \class MpdFwd
/// \brief Fwd detector construction and hit processing
/// \author Evgeny Kryshen (PNPI, Gatchina)
//------------------------------------------------------------------------------------------------------------------------

#include "TClonesArray.h"
#include "TObjArray.h"
#include "TLorentzVector.h"
#include "TParticle.h"
#include "TGeoManager.h"
#include "TGeoTube.h"
#include "TGeoBBox.h"
#include "TGeoVolume.h"
#include "TVirtualMC.h"
#include "TString.h"

#include "FairGeoInterface.h"
#include "FairGeoLoader.h"
#include "FairGeoNode.h"
#include "FairGeoMedia.h"
#include "FairGeoMedium.h"
#include "FairGeoRootBuilder.h"
#include "FairRootManager.h"
#include "FairRun.h"
#include "FairRuntimeDb.h"
#include "FairVolume.h"

#include "MpdMCTrack.h"
#include "MpdStack.h"
#include "MpdFwd.h"
#include "MpdFwdPoint.h"

ClassImp(MpdFwd)
MpdFwd::MpdFwd(const char *name) : FairDetector(name, kTRUE),
fFwdPoints(0),
fTrackID(0),
fPosIn(),
fPosOut(),
fTime(0),
fLength(0),
fELoss(0)
{
  fFwdPoints = new TClonesArray("MpdFwdPoint");
}

MpdFwd::~MpdFwd(){
  fFwdPoints->Delete();
  delete fFwdPoints;
}

void MpdFwd::EndOfEvent(){
  fFwdPoints->Clear();
}

void MpdFwd::Reset(){
  fFwdPoints->Clear();
}

void MpdFwd::Register(){
  FairRootManager::Instance()->Register("FwdPoint", "Fwd", fFwdPoints, kTRUE);
}

TClonesArray* MpdFwd::GetCollection(Int_t iColl) const {
  if (iColl == 0) return fFwdPoints;
  return nullptr;
}

Bool_t MpdFwd::ProcessHits(FairVolume *vol){
  if (gMC->IsTrackEntering()) {
    fELoss  = 0.;
    gMC->TrackPosition(fPosIn);
  }

  // Sum energy loss for all steps in the active volume
  fELoss += gMC->Edep();

  // Set additional parameters at exit of active volume. Create CbmMuchPoint.
  if (gMC->IsTrackExiting() || gMC->IsTrackStop() || gMC->IsTrackDisappeared()) {
    fTrackID = gMC->GetStack()->GetCurrentTrackNumber();
    gMC->TrackPosition(fPosOut);
    double trackLength = gMC->TrackLength();
    // if (fELoss == 0.) return kFALSE;
    new ((*fFwdPoints)[fFwdPoints->GetEntriesFast()]) MpdFwdPoint (fTrackID, fPosIn, fPosOut, trackLength, fELoss);
    ((MpdStack*) gMC->GetStack())->AddPoint(kFWD, gMC->GetStack()->GetCurrentParentTrackNumber());
  }
  return kTRUE;
}

void MpdFwd::ConstructGeometry(){
  FairGeoLoader*    geoLoad  = FairGeoLoader::Instance();
  FairGeoBuilder*   geoBuild = geoLoad->getGeoBuilder();
  FairGeoInterface* geoFace  = geoLoad->getGeoInterface();
  FairGeoMedia*     geoMedia = geoFace->getMedia();
  // Materials
  TGeoMedium* mVacuum   = gGeoManager->GetMedium(geoBuild->createMedium(geoMedia->getMedium("vacuum")));
  TGeoMedium* mAir      = gGeoManager->GetMedium(geoBuild->createMedium(geoMedia->getMedium("air")));
  TGeoMedium* mKapton   = gGeoManager->GetMedium(geoBuild->createMedium(geoMedia->getMedium("kapton")));
  TGeoMedium* mTungsten = gGeoManager->GetMedium(geoBuild->createMedium(geoMedia->getMedium("tungsten")));
  TGeoMedium* mGas      = gGeoManager->GetMedium(geoBuild->createMedium(geoMedia->getMedium("arco27030")));
  // Rotation matrices
  TGeoRotation* rot0   = new TGeoRotation("rot0",     0., 0., 0.);
  TGeoRotation* rotP00 = new TGeoRotation("rotP00",   0., 0., 0.);
  TGeoRotation* rotP90 = new TGeoRotation("rotP90",  90., 0., 0.);
  TGeoRotation* rotP45 = new TGeoRotation("rotP45",  45., 0., 0.);
  TGeoRotation* rotM45 = new TGeoRotation("rotM45", -45., 0., 0.);
  TGeoRotation* rotTub = new TGeoRotation("rotTub",   0.,90., 0.);
  
  // Detector parameters
  // Fwd
  Double_t eps = 1e-10;         // cm small margin
  Double_t fwdRMin = 35;        // cm
  Double_t fwdRMax = 130 + eps; // cm
  Double_t fwdHalfZ = 50;       // cm
  Double_t fwdPosZ = 265;       // cm
  // Station
  Int_t nStations = 10;
  Double_t stationRMin = 35.7; // cm
  Double_t stationRMax = 130;  // cm
  Double_t stationHalfZ = 4;   // cm
  Double_t stationDist = 2.*(fwdHalfZ-stationHalfZ)/(nStations-1.); // distance between station centers
  // Tube parameters
  Double_t tubeRMin = 0.5;               // cm
  Double_t tubeRMax = tubeRMin + 0.003;  // cm
  Double_t wireRMin = 0;                 // cm
  Double_t wireRMax = 0.0015;            // cm
  Double_t gasRMin  = wireRMax + eps;    // cm
  Double_t gasRMax  = tubeRMin - eps;    // cm
  Double_t tubeDist = 1.5;               // cm, distance between tube centers
  // Layer
  Double_t layerRMin = stationRMin + eps;
  Double_t layerRMax = stationRMax - eps;
  Double_t layerDist = 1.5; // cm, distance between layer centers
  Double_t layerHalfZ = tubeRMax + eps;
  // Sensitive plane
  Double_t planeHalfZ = 0.01; // cm

  // create two mother volumes on both sides from the IP
  TGeoTube* sFwd = new TGeoTube("sFwd", fwdRMin, fwdRMax, fwdHalfZ);
  TGeoVolume* vFwd = new TGeoVolume("vFwd", sFwd, mAir);
  Double_t* buf = NULL;
  gGeoManager->Node("vFwd", 0, "cave", 0., 0.,  fwdPosZ, 0, kTRUE, buf, 0);
  gGeoManager->Node("vFwd", 1, "cave", 0., 0., -fwdPosZ, 0, kTRUE, buf, 0);

  // create n stations in each mother volume 
  TGeoTube* sStation = new TGeoTube("sStation",stationRMin, stationRMax, stationHalfZ);
  TGeoVolume* vStation = new TGeoVolume("vStation",sStation, mAir);
  for (Int_t iStation = 0; iStation < nStations; iStation++){
    vFwd->AddNode(vStation,iStation,new TGeoCombiTrans(0,0,(-(nStations-1)/2.+iStation)*stationDist,rot0));
  }

  // create 4 layers with different tube orientations and one sensitive plane in the middle
  TGeoTube* sPlane = new TGeoTube("sPlane", layerRMin, layerRMax, planeHalfZ);
  TGeoTube* sLayer = new TGeoTube("sLayer", layerRMin, layerRMax, layerHalfZ);
  TGeoVolume* vPlane = new TGeoVolume("vPlane",sPlane, mAir);
  TGeoVolume* vLayer = new TGeoVolume("vLayer", sLayer, mAir);
  vStation->AddNode(vLayer, 0, new TGeoCombiTrans(0, 0, -1.5*layerDist, rotP00));
  vStation->AddNode(vLayer, 1, new TGeoCombiTrans(0, 0, -0.5*layerDist, rotP90));
  vStation->AddNode(vLayer, 2, new TGeoCombiTrans(0, 0,  0.5*layerDist, rotP45));
  vStation->AddNode(vLayer, 3, new TGeoCombiTrans(0, 0,  1.5*layerDist, rotM45));
  vStation->AddNode(vPlane, 0, new TGeoCombiTrans(0, 0, 0, rotP00));
  
  // Fill layers with tubes
  Int_t nTubes = TMath::Floor((layerRMax - tubeRMax)/tubeDist);
  for (Int_t iTube = 0; iTube < nTubes; iTube++){
    Double_t xPos = iTube*tubeDist;
    Double_t xMax = xPos + tubeRMax;
    Double_t xMin = xPos - tubeRMax;
    Double_t yMax = layerRMax > xMax ? TMath::Sqrt(layerRMax*layerRMax - xMax*xMax) : 0;
    Double_t yMin = layerRMin > xMin ? TMath::Sqrt(layerRMin*layerRMin - xMin*xMin) : 0;
    Double_t yPos = (yMax+yMin)/2.;
    Double_t tubeHalfLength = (yMax - yMin - eps)/2.;
    TGeoTube* sTube = new TGeoTube(Form("sTube%d",iTube), tubeRMin, tubeRMax, tubeHalfLength);
    TGeoTube* sGas  = new TGeoTube(Form("sGas%d",iTube), gasRMin, gasRMax, tubeHalfLength);
    TGeoTube* sWire = new TGeoTube(Form("sWire%d",iTube), wireRMin, wireRMax, tubeHalfLength);
    TGeoVolume* vTube = new TGeoVolume(Form("vTube%d",iTube), sTube, mKapton);
    TGeoVolume* vGas  = new TGeoVolume(Form("vGas%d",iTube), sGas, mGas);
    TGeoVolume* vWire = new TGeoVolume(Form("vWire%d",iTube), sWire, mTungsten);
    vLayer->AddNode(vTube, 0, new TGeoCombiTrans( xPos, yPos, 0, rotTub));
    vLayer->AddNode(vTube, 1, new TGeoCombiTrans(-xPos, yPos, 0, rotTub));
    vLayer->AddNode(vTube, 2, new TGeoCombiTrans( xPos,-yPos, 0, rotTub));
    vLayer->AddNode(vTube, 3, new TGeoCombiTrans(-xPos,-yPos, 0, rotTub));
    vLayer->AddNode(vGas,  0, new TGeoCombiTrans( xPos, yPos, 0, rotTub));
    vLayer->AddNode(vGas,  1, new TGeoCombiTrans(-xPos, yPos, 0, rotTub));
    vLayer->AddNode(vGas,  2, new TGeoCombiTrans( xPos,-yPos, 0, rotTub));
    vLayer->AddNode(vGas,  3, new TGeoCombiTrans(-xPos,-yPos, 0, rotTub));
    vLayer->AddNode(vWire, 0, new TGeoCombiTrans( xPos, yPos, 0, rotTub));
    vLayer->AddNode(vWire, 1, new TGeoCombiTrans(-xPos, yPos, 0, rotTub));
    vLayer->AddNode(vWire, 2, new TGeoCombiTrans( xPos,-yPos, 0, rotTub));
    vLayer->AddNode(vWire, 3, new TGeoCombiTrans(-xPos,-yPos, 0, rotTub));
    vTube->SetLineColor(kYellow);
  }
  
  // Set color for visualization
  vStation->SetLineColor(kYellow);
  vLayer->SetLineColor(kYellow);
  vPlane->SetLineColor(kGreen);
  
  // Set sensitive volumes
  AddSensitiveVolume(vPlane); // use thin plane in the middle of the station. TODO: switch to vGas
}

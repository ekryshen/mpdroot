//------------------------------------------------------------------------------------------------------------------------
/// \class MpdFwd
///
/// \brief
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
//  printf("Process hits\n");

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
//    if (fELoss == 0.) return kFALSE;
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
  FairGeoMedium*    mat = geoMedia->getMedium("vacuum");
  FairGeoMedium*    dat = geoMedia->getMedium("vacuum");
  
  FairGeoMedium*    mStation = geoMedia->getMedium("air");

  //  FairGeoMedium*    mat = geoMedia->getMedium("air");
//  FairGeoMedium*    dat = geoMedia->getMedium("air");
//  FairGeoMedium*    carbon = geoMedia->getMedium("carbon");
//  FairGeoMedium*    arco = geoMedia->getMedium("arco27030");
  FairGeoMedium*    tube = geoMedia->getMedium("kapton");
  FairGeoMedium*    wire = geoMedia->getMedium("tungsten");
  FairGeoMedium*    gas = geoMedia->getMedium("arco27030");
//  FairGeoMedium*    mFwd = geoMedia->getMedium("vacuum");
//  FairGeoMedium*    dFwd = geoMedia->getMedium("vacuum");
  Int_t imStation = geoBuild->createMedium(mStation);
  
  Int_t iMat = geoBuild->createMedium(mat);
  Int_t iDat = geoBuild->createMedium(dat);
  Int_t iTube = geoBuild->createMedium(tube);
  Int_t iWire = geoBuild->createMedium(wire);
  Int_t iGas = geoBuild->createMedium(gas);
//  Int_t iCarbon = geoBuild->createMedium(carbon);
//  Int_t iArco = geoBuild->createMedium(arco);
  Double_t* buf = NULL;
  
  Double_t a = 11; //cm 
  Double_t l = 2*162; // cm tube length
  Double_t d = 1.0; // cm, tube diameter
  Double_t ds = 0.003; // cm tube wall
  Double_t dt = 0.5; // cm, distance Between Tubes
  Double_t dl = 0.1; // cm, distance Between Layers
  Double_t dw = 0.003; // cm wire diameter
//  Double_t dx = d/2+dl/2-0.001;
  Int_t n = 200;
  Int_t k =0;
  if (n%2 == 0){
    k = dt/2 + d/2 + ds + (d+dt+2*ds)*(n/2-1);
  }else{
    k = d + 2*ds + dt + (dt+d+2*ds)*((n-1)/2-1);
  }
  
  TGeoTube* sFwd = new TGeoTube("sFwd", 0.00000001, 500., 70);
  TGeoVolume* vFwd = new TGeoVolume("vFwd",sFwd,gGeoManager->GetMedium(iMat));
  TGeoBBox* sLayer = new TGeoBBox("sLayer", k/2+0.5, 0.526, l/2+1); // half-sizes
  TGeoVolume* vLayer = new TGeoVolume("vLayer",sLayer,gGeoManager->GetMedium(iDat));
  vLayer->SetLineColor(kYellow+1);

  
  TGeoTube* sTube= new TGeoTube("sTube", d/2, d/2+ds, l/2);


  TGeoVolume* vTube = new TGeoVolume("vTube",sTube,gGeoManager->GetMedium(iTube));
  vTube->SetLineColor(kYellow+1);

  TGeoTube* sWire= new TGeoTube("sWire", 0.00000001, dw/2, l/2);
  TGeoVolume* vWire = new TGeoVolume("vWire",sWire,gGeoManager->GetMedium(iWire));
  TGeoTube* sGas = new TGeoTube("sGas", dw/2+0.00000001, d/2-0.00000001, l/2);
//  TGeoTube* sGas = new TGeoTube("sGas", 0, d/2, l);
  TGeoVolume* vGas= new TGeoVolume("vGas",sGas,gGeoManager->GetMedium(iGas));
  TGeoTube* sPlate= new TGeoTube("sPlate", 7., l/2, 0.2);
  TGeoVolume* vPlate = new TGeoVolume("vPlate",sPlate,gGeoManager->GetMedium(iMat));
  vPlate->SetLineColor(kYellow);
  
  TGeoRotation* rot0 = new TGeoRotation("rot0", 0., 0., 0.);
  TGeoRotation* rot = new TGeoRotation("rot", 0., 90., 0.);
  TGeoRotation* rot1 = new TGeoRotation("rot1", 45., 90., 0.);
  TGeoRotation* rot2 = new TGeoRotation("rot2", -45., 90., 0.);
  TGeoRotation* rot3 = new TGeoRotation("rot3", 90., 90., 0.);
  
  
//  vLayer->AddNode(vTube,0,new TGeoCombiTrans(0,0,0,rot0));
//  vLayer->AddNode(vWire,0,new TGeoCombiTrans(0,0,0,rot0));
//  vLayer->AddNode(vGas,0,new TGeoCombiTrans(0,0,0,rot0));
// for (Int_t i = 0; i<n; i++){
//      vLayer->AddNode(vTube,i,new TGeoCombiTrans(-k+i*((d/2+ds)*2+dt),0,0,rot0));
//      vLayer->AddNode(vWire,i,new TGeoCombiTrans(-k+i*((d/2+ds)*2+dt),0,0,rot0));
//      vLayer->AddNode(vGas,i,new TGeoCombiTrans(-k+i*((d/2+ds)*2+dt),0,0,rot0));
//  }
  double stationRMin = 35.7; // cm
  double stationRMax = 130;  // cm
  double stationHalfZ = 4;   // cm
 
  gGeoManager->Node("vFwd",0,"cave",0.,0.,265.,0,kTRUE,buf,0);
  gGeoManager->Node("vFwd",1,"cave",0.,0.,-265.,0,kTRUE,buf,0);

  TGeoTube* sStation = new TGeoTube("sStation",stationRMin,stationRMax,stationHalfZ);
  TGeoVolume* vStation = new TGeoVolume("vStation",sStation,gGeoManager->GetMedium(imStation));
//  vStation->AddNode(vLayer,0,new TGeoCombiTrans(0,0,-3.,rot));
//  vStation->AddNode(vLayer,1,new TGeoCombiTrans(0,0,-1,rot1));
//  vStation->AddNode(vPlate,0,new TGeoCombiTrans(0,0, 0,rot0));
//  vStation->AddNode(vLayer,2,new TGeoCombiTrans(0,0, 1,rot2));
//  vStation->AddNode(vLayer,3,new TGeoCombiTrans(0,0, 3,rot3));
  vStation->SetLineColor(kYellow);
  
  vFwd->AddNode(vStation,0,new TGeoCombiTrans(0,0,-4.5*a,rot0));
  vFwd->AddNode(vStation,1,new TGeoCombiTrans(0,0,-3.5*a,rot0));
  vFwd->AddNode(vStation,0,new TGeoCombiTrans(0,0,-2.5*a,rot0));
  vFwd->AddNode(vStation,1,new TGeoCombiTrans(0,0,-1.5*a,rot0));
  vFwd->AddNode(vStation,2,new TGeoCombiTrans(0,0,-0.5*a,rot0));
  vFwd->AddNode(vStation,3,new TGeoCombiTrans(0,0, 0.5*a,rot0));
  vFwd->AddNode(vStation,3,new TGeoCombiTrans(0,0, 1.5*a,rot0));
  vFwd->AddNode(vStation,4,new TGeoCombiTrans(0,0, 2.5*a,rot0));
  vFwd->AddNode(vStation,3,new TGeoCombiTrans(0,0, 3.5*a,rot0));
  vFwd->AddNode(vStation,4,new TGeoCombiTrans(0,0, 4.5*a,rot0));

//  vFwd->AddNode(vLayer,0,new TGeoCombiTrans(0,0,-2*a,rot));
//  vFwd->AddNode(vLayer,1,new TGeoCombiTrans(0,0,-2*a+2,rot1));
//  vFwd->AddNode(vPlate,0,new TGeoCombiTrans(0,0,-2*a+3,rot0));
//  vFwd->AddNode(vLayer,2,new TGeoCombiTrans(0,0,-2*a+4,rot2));
//  vFwd->AddNode(vLayer,3,new TGeoCombiTrans(0,0,-2*a+6,rot3));
  
//  vFwd->AddNode(vLayer,4,new TGeoCombiTrans(0,0,-a,rot));
//  vFwd->AddNode(vLayer,5,new TGeoCombiTrans(0,0,-a+2,rot1));
//  vFwd->AddNode(vPlate,1,new TGeoCombiTrans(0,0,-a+3,rot0));
//  vFwd->AddNode(vLayer,6,new TGeoCombiTrans(0,0,-a+4,rot2));
//  vFwd->AddNode(vLayer,7,new TGeoCombiTrans(0,0,-a+6,rot3));
//  
//  vFwd->AddNode(vLayer,8,new TGeoCombiTrans(0,0,0,rot));
//  vFwd->AddNode(vLayer,9,new TGeoCombiTrans(0,0,0+2,rot1));
//  vFwd->AddNode(vPlate,2,new TGeoCombiTrans(0,0,0+3,rot0));
//  vFwd->AddNode(vLayer,10,new TGeoCombiTrans(0,0,0+4,rot2));
//  vFwd->AddNode(vLayer,11,new TGeoCombiTrans(0,0,0+6,rot3));
//  
//  vFwd->AddNode(vLayer,12,new TGeoCombiTrans(0,0,a,rot));
//  vFwd->AddNode(vLayer,13,new TGeoCombiTrans(0,0,a+2,rot1));
//  vFwd->AddNode(vPlate,3,new TGeoCombiTrans(0,0,a+3,rot0));
//  vFwd->AddNode(vLayer,14,new TGeoCombiTrans(0,0,a+4,rot2));
//  vFwd->AddNode(vLayer,15,new TGeoCombiTrans(0,0,a+6,rot3));
//  
//  vFwd->AddNode(vLayer,16,new TGeoCombiTrans(0,0,2*a,rot));
//  vFwd->AddNode(vLayer,17,new TGeoCombiTrans(0,0,2*a+2,rot1));
//  vFwd->AddNode(vPlate,4,new TGeoCombiTrans(0,0,2*a+3,rot0));
//  vFwd->AddNode(vLayer,18,new TGeoCombiTrans(0,0,2*a+4,rot2));
//  vFwd->AddNode(vLayer,19,new TGeoCombiTrans(0,0,2*a+6,rot3));
  
//  TGeoTube* sFwd = new TGeoTube("sFwd", 0.0000001, 148., 50);// cm
//  TGeoVolume* vFwd = new TGeoVolume("vFwd",sFwd,gGeoManager->GetMedium(iMat));
//  TGeoTube* dsFwd = new TGeoTube("sFwd", 0.001, 147., 0.00001); // cm
//  TGeoTube* sCarbon = new TGeoTube("sCarbon", 0.0000001, 147., 0.09826275787); // cm
//  TGeoVolume* vCarbon = new TGeoVolume("vCarbon",sCarbon,gGeoManager->GetMedium(iCarbon));
//  TGeoTube* sArco= new TGeoTube("sArco", 0.0000001, 147., 0.5);
//  TGeoVolume* vArco = new TGeoVolume("vArco",sArco,gGeoManager->GetMedium(iArco));
//  TGeoTube* sWire= new TGeoTube("sWire", 0.0000001, 147., 0.5);
//  TGeoVolume* vWire = new TGeoVolume("vWire",sWire,gGeoManager->GetMedium(iWire));
//  TGeoTube* sTube= new TGeoTube("sArco", 0.0000001, 147., 0.5);
//  TGeoVolume* vTube = new TGeoVolume("vArco",sTube,gGeoManager->GetMedium(iTube));
//  gGeoManager->Node("vFwd",0,"cave",0.,0.,269.,0,kTRUE,buf,0);
//  gGeoManager->Node("vArco",0,"vFwd",0.,0.,-40.,0,kTRUE,buf,0);
//  gGeoManager->Node("dvFwd",1,"vFwd",0.,0.,-20.,0,kTRUE,buf,0);
//  gGeoManager->Node("dvFwd",2,"vFwd",0.,0.,  0.,0,kTRUE,buf,0);
//  gGeoManager->Node("dvFwd",3,"vFwd",0.,0., 20.,0,kTRUE,buf,0);
//  gGeoManager->Node("dvFwd",4,"vFwd",0.,0., 40.,0,kTRUE,buf,0);
  
//  AddSensitiveVolume(vGas);
//  AddSensitiveVolume(vArco);
//  AddSensitiveVolume(vTube);
//  AddSensitiveVolume(vWire);
 AddSensitiveVolume(vStation);
}

//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      Implementation of class TpcDetector
//      see TpcDetector.hh for details
//
// Environment:
//      Software developed for the PANDA Detector at FAIR.
//
// Author List:
//      Sebastian Neubert    TUM            (original author)
//
//
//-----------------------------------------------------------

// Panda Headers ----------------------

// This Class' Header ------------------
#include "TpcDetector.h"
#include "MpdTpcEDepParams.h"

// C/C++ Headers ----------------------
#include <iostream>
using namespace std;

// Collaborating Class Headers --------
#include "TClonesArray.h"
#include "FairRootManager.h"
#include "TpcPoint.h"
#include "TVirtualMC.h"
#include "TLorentzVector.h"
#include "TpcGeo.h"
#include "FairGeoLoader.h"
#include "FairGeoInterface.h"
#include "FairGeoNode.h"
#include "FairGeoRootBuilder.h"
#include "TList.h"
#include "FairRun.h"
#include "FairRuntimeDb.h"
#include "TpcGeoPar.h"
#include "TObjArray.h"
#include <TSystem.h>
#include "FairGeoNode.h"
#include "FairGeoVolume.h"
#include "FairVolume.h"
#include "MpdStack.h"
#include "TFile.h"
#include "TGeoManager.h"
#include "TGDMLParse.h"
#include "FairGeoMedia.h"

// Class Member definitions -----------


TpcDetector::TpcDetector(const char * Name, Bool_t Active)
  : FairDetector(Name, Active)
{

  fTpcPointCollection= new TClonesArray("TpcPoint");

  fVerboseLevel = 1;
}

TpcDetector::TpcDetector()
{
  fTpcPointCollection= new TClonesArray("TpcPoint");

  fVerboseLevel = 1;
}

TpcDetector::~TpcDetector()
{
  if (fTpcPointCollection) {
      fTpcPointCollection->Delete();
      delete fTpcPointCollection;
   }
}




void TpcDetector::EndOfEvent()
{
   if(fVerboseLevel) Print();
   fTpcPointCollection->Delete();

   //AZ - memory monitor
   ProcInfo_t proc;
   gSystem->GetProcInfo(&proc);
   cout << " User CPU time: " << proc.fCpuUser << " Memory: resident " << proc.fMemResident << ", virtual " << proc.fMemVirtual << endl;
}



void TpcDetector::Register() {

/** This will create a branch in the output tree called  TpcDetectorPoint, setting the last parameter to kFALSE means:

    this collection will not be written to the file, it will exist only during the simulation. */

   FairRootManager::Instance()->Register("TpcPoint", "Tpc", fTpcPointCollection, kTRUE);

}

//---------------------------------------------------------------------------
Bool_t TpcDetector::ProcessHits( FairVolume *v)
{
  // create Hit for every MC step

  Double_t time   = gMC->TrackTime() * 1.0e09;
  Double_t length = gMC->TrackLength();
  TLorentzVector pos;
  gMC->TrackPosition(pos);
  //  pos.Print();
  TLorentzVector mom;
  gMC->TrackMomentum(mom);
  //  mom.Print();
  Double_t eLoss = gMC->Edep();
  Int_t trackID  = gMC->GetStack()->GetCurrentTrackNumber();
  Int_t volumeID = v->getMCid();

  Float_t pX=0., pY=0., pZ=0., eTot=0.;
  gMC->TrackMomentum( pX, pY, pZ, eTot);
  if (eLoss > 0 && eTot>0.001) { //AZ - hard-coded 1 MeV cut
  //if (eLoss > 0) { //AZ
    //AZ - dE/dx simulation from the new prescription
    //cout << " xxxxxxxxxxxxx " << gMC->TrackMass() << endl;
    Double_t betgam = mom.Vect().Mag() / gMC->TrackMass();
    //cout << betgam << " " << mom.Beta() * mom.Gamma() << endl;
    if (betgam > 0.01) {
      betgam = TMath::Log10(betgam);
      MpdTpcEDepParams *deParams = MpdTpcEDepParams::Instance();
      //IR 05-APR-2021 if (betgam >= deParams->GetMinLimit() && betgam < deParams->GetMaxLimit()) {
      betgam = TMath::Min( betgam, deParams->GetMaxLimit());
      if (betgam >= deParams->GetMinLimit()) {
	// New value
	eLoss = deParams->GetEloss (betgam, gMC->TrackCharge(), gMC->TrackStep());
	//cout << " xxxxxx " << gMC->Edep() << " " << eLoss << endl;
	//AZ - fixed loss 2 keV/cm
	//eLoss = 2.e-6;
	//eLoss *= gMC->TrackStep();
	//AZ
	TpcPoint* p = AddHit(trackID, volumeID, pos.Vect(), mom.Vect(), time, length, eLoss);
	p->SetStep(gMC->TrackStep());
	
	((MpdStack*)gMC->GetStack())->AddPoint(kTPC);
	
      } 
    }
  }
  //  p->Print("");

  return kTRUE;
}

//---------------------------------------------------------------------------

TClonesArray* TpcDetector::GetCollection(Int_t iColl) const {
  if (iColl == 0) return fTpcPointCollection;
  else return nullptr;
}


TpcPoint* TpcDetector::AddHit(Int_t trackID, Int_t detID, TVector3 pos,
    TVector3 mom, Double_t time, Double_t length,
    Double_t eLoss) {
  TClonesArray& clref = *fTpcPointCollection;
  Int_t size = clref.GetEntriesFast();
  return new(clref[size]) TpcPoint(trackID, detID, pos, mom,
      time, length, eLoss);
}

void TpcDetector::Reset() {
  fTpcPointCollection->Delete();
  //ResetParameters();
}

void TpcDetector::Print() const
{
        Int_t nHits = fTpcPointCollection->GetEntriesFast();
        cout<<"-I- MpdTpc: " << nHits << " points registered in this event." << endl;

        if(fVerboseLevel > 1)
                for(Int_t i=0; i<nHits; i++) (*fTpcPointCollection)[i]->Print();
}

void TpcDetector::ConstructGeometry() {
    TString fileName = GetGeometryFileName();

    if ( fileName.EndsWith(".root") ) {
      LOG(INFO) << "Constructing TPC geometry from ROOT file " << fileName.Data();
      ConstructRootGeometry();
    }
    else if ( fileName.EndsWith(".geo") ) {
      LOG(INFO) << "Constructing TPC geometry from ASCII file " << fileName.Data();
      ConstructAsciiGeometry();
    }
    else if ( fileName.EndsWith(".gdml") )
    {
      LOG(INFO) << "Constructing TPC geometry from GDML file " << fileName.Data();
      //ConstructGDMLGeometry();
    }
    else
    {
      LOG(FATAL) << "Geometry format of TPC file " << fileName.Data() << " not supported.";
    }
}

// -----   ConstructAsciiGeometry   -------------------------------------------
void TpcDetector::ConstructAsciiGeometry() {

  FairGeoLoader*    geoLoad = FairGeoLoader::Instance();
  FairGeoInterface* geoFace = geoLoad->getGeoInterface();
  TpcGeo*       TPCGeo  = new TpcGeo();
  TPCGeo->setGeomFile(GetGeometryFileName());
  geoFace->addGeoModule(TPCGeo);

  Bool_t rc = geoFace->readSet(TPCGeo);
  if (rc) TPCGeo->create(geoLoad->getGeoBuilder());
  TList* volList = TPCGeo->getListOfVolumes();
  // store geo parameter
  FairRun *fRun = FairRun::Instance();
  FairRuntimeDb *rtdb= FairRun::Instance()->GetRuntimeDb();
  TpcGeoPar* par=(TpcGeoPar*)(rtdb->getContainer("TpcGeoPar"));
  TObjArray *fSensNodes = par->GetGeoSensitiveNodes();
  TObjArray *fPassNodes = par->GetGeoPassiveNodes();

  TListIter iter(volList);
  FairGeoNode* node   = nullptr;
  FairGeoVolume *aVol=nullptr;

  while( (node = (FairGeoNode*)iter.Next()) ) {
      aVol = dynamic_cast<FairGeoVolume*> ( node );
       if ( node->isSensitive()  ) {
           fSensNodes->AddLast( aVol );
       }else{
           fPassNodes->AddLast( aVol );
       }
  }
  par->setChanged();
  par->setInputVersion(fRun->GetRunId(),1);
  ProcessNodes( volList );
}
// ----------------------------------------------------------------------------

// -----   ConstructGDMLGeometry   -------------------------------------------
void TpcDetector::ConstructGDMLGeometry()
{
    TFile *old = gFile;
    TGDMLParse parser;
    TGeoVolume* gdmlTop;

    // Before importing GDML
    Int_t maxInd = gGeoManager->GetListOfMedia()->GetEntries() - 1;

    gdmlTop = parser.GDMLReadFile(GetGeometryFileName());

    // Cheating - reassigning media indices after GDML import (need to fix this in TGDMLParse class!!!)
    //   for (Int_t i=0; i<gGeoManager->GetListOfMedia()->GetEntries(); i++)
    //      gGeoManager->GetListOfMedia()->At(i)->Dump();
    // After importing GDML
    Int_t j = gGeoManager->GetListOfMedia()->GetEntries() - 1;
    Int_t curId;
    TGeoMedium* m;
    do {
        m = (TGeoMedium*)gGeoManager->GetListOfMedia()->At(j);
        curId = m->GetId();
        m->SetId(curId+maxInd);
        j--;
    } while (curId > 1);
    //   LOG(DEBUG) << "====================================================================";
    //   for (Int_t i=0; i<gGeoManager->GetListOfMedia()->GetEntries(); i++)
    //      gGeoManager->GetListOfMedia()->At(i)->Dump();

    Int_t newMaxInd = gGeoManager->GetListOfMedia()->GetEntries() - 1;

    gGeoManager->GetTopVolume()->AddNode(gdmlTop, 1, 0);
    ExpandNodeForGdml(gGeoManager->GetTopVolume()->GetNode(gGeoManager->GetTopVolume()->GetNdaughters()-1));

    for (Int_t k = maxInd+1; k < newMaxInd+1; k++) {
        TGeoMedium* medToDel = (TGeoMedium*)(gGeoManager->GetListOfMedia()->At(maxInd+1));
        LOG(DEBUG) << "    removing media " << medToDel->GetName() << " with id " << medToDel->GetId() << " (k=" << k << ")";
        gGeoManager->GetListOfMedia()->Remove(medToDel);
    }
    gGeoManager->SetAllIndex();

    gFile = old;
}

void TpcDetector::ExpandNodeForGdml(TGeoNode* node)
{
    LOG(DEBUG) << "----------------------------------------- ExpandNodeForGdml for node " << node->GetName();

    TGeoVolume* curVol = node->GetVolume();

    LOG(DEBUG) << "    volume: " << curVol->GetName();

    if (curVol->IsAssembly()) {
        LOG(DEBUG) << "    skipping volume-assembly";
    }
    else
    {
        TGeoMedium* curMed = curVol->GetMedium();
        TGeoMaterial* curMat = curVol->GetMaterial();
        TGeoMedium* curMedInGeoManager = gGeoManager->GetMedium(curMed->GetName());
        TGeoMaterial* curMatOfMedInGeoManager = curMedInGeoManager->GetMaterial();
        TGeoMaterial* curMatInGeoManager = gGeoManager->GetMaterial(curMat->GetName());

        // Current medium and material assigned to the volume from GDML
        LOG(DEBUG2) << "    curMed\t\t\t\t" << curMed << "\t" << curMed->GetName() << "\t" << curMed->GetId();
        LOG(DEBUG2) << "    curMat\t\t\t\t" << curMat << "\t" << curMat->GetName() << "\t" << curMat->GetIndex();

        // Medium and material found in the gGeoManager - either the pre-loaded one or one from GDML
        LOG(DEBUG2) << "    curMedInGeoManager\t\t" << curMedInGeoManager
                 << "\t" << curMedInGeoManager->GetName() << "\t" << curMedInGeoManager->GetId();
        LOG(DEBUG2) << "    curMatOfMedInGeoManager\t\t" << curMatOfMedInGeoManager
                 << "\t" << curMatOfMedInGeoManager->GetName() << "\t" << curMatOfMedInGeoManager->GetIndex();
        LOG(DEBUG2) << "    curMatInGeoManager\t\t" << curMatInGeoManager
                 << "\t" << curMatInGeoManager->GetName() << "\t" << curMatInGeoManager->GetIndex();

        TString matName = curMat->GetName();
        TString medName = curMed->GetName();

        if (curMed->GetId() != curMedInGeoManager->GetId()) {
            if (fFixedMedia.find(medName) == fFixedMedia.end()) {
                LOG(DEBUG) << "    Medium needs to be fixed";
                fFixedMedia[medName] = curMedInGeoManager;
                Int_t ind = curMat->GetIndex();
                gGeoManager->RemoveMaterial(ind);
                LOG(DEBUG) << "    removing material " << curMat->GetName()
                    << " with index " << ind;
                for (Int_t i=ind; i<gGeoManager->GetListOfMaterials()->GetEntries(); i++) {
                    TGeoMaterial* m = (TGeoMaterial*)gGeoManager->GetListOfMaterials()->At(i);
                    m->SetIndex(m->GetIndex()-1);
                }

                LOG(DEBUG) << "    Medium fixed";
            }
            else
            {
                LOG(DEBUG) << "    Already fixed medium found in the list    ";
            }
        }
        else
        {
            if (fFixedMedia.find(medName) == fFixedMedia.end()) {
                LOG(DEBUG) << "    There is no correct medium in the memory yet";

                FairGeoLoader* geoLoad = FairGeoLoader::Instance();
                FairGeoInterface* geoFace = geoLoad->getGeoInterface();
                FairGeoMedia* geoMediaBase =  geoFace->getMedia();
                FairGeoBuilder* geobuild = geoLoad->getGeoBuilder();

                FairGeoMedium* curMedInGeo = geoMediaBase->getMedium(medName);
                if (curMedInGeo == 0)
                {
                    LOG(FATAL) << "    Media not found in Geo file: " << medName;
                    //! This should not happen.
                    //! This means that somebody uses material in GDML that is not in the media.geo file.
                    //! Most probably this is the sign to the user to check materials' names in the CATIA model.
                }
                else
                {
                    LOG(DEBUG) << "    Found media in Geo file" << medName;
                    Int_t nmed = geobuild->createMedium(curMedInGeo);
                    LOG(DEBUG) << "    nmed " << nmed;
                    fFixedMedia[medName] = (TGeoMedium*)gGeoManager->GetListOfMedia()->Last();
                    gGeoManager->RemoveMaterial(curMatOfMedInGeoManager->GetIndex());
                    LOG(DEBUG) << "    removing material " << curMatOfMedInGeoManager->GetName()
                        << " with index " << curMatOfMedInGeoManager->GetIndex();
                    for (Int_t i=curMatOfMedInGeoManager->GetIndex(); i<gGeoManager->GetListOfMaterials()->GetEntries(); i++) {
                        TGeoMaterial* m = (TGeoMaterial*)gGeoManager->GetListOfMaterials()->At(i);
                        m->SetIndex(m->GetIndex()-1);
                    }
                }

                if (curMedInGeo->getSensitivityFlag()) {
                    LOG(DEBUG) << "    Adding sensitive  " << curVol->GetName();
                    AddSensitiveVolume(curVol);
                }
            }
            else
            {
                LOG(DEBUG) << "    Already fixed medium found in the list";
                LOG(DEBUG) << "!!! Sensitivity: " << fFixedMedia[medName]->GetParam(0);
                if (fFixedMedia[medName]->GetParam(0) == 1) {
                    LOG(DEBUG) << "    Adding sensitive  " << curVol->GetName();
                    AddSensitiveVolume(curVol);
                }
            }
        }

        curVol->SetMedium(fFixedMedia[medName]);
        gGeoManager->SetAllIndex();

  //      gGeoManager->GetListOfMaterials()->Print();
  //      gGeoManager->GetListOfMedia()->Print();

    }

    //! Recursevly go down the tree of nodes
    if (curVol->GetNdaughters() != 0)
    {
        TObjArray* NodeChildList = curVol->GetNodes();
        TGeoNode* curNodeChild;
        for (Int_t j=0; j<NodeChildList->GetEntriesFast(); j++)
        {
            curNodeChild = (TGeoNode*)NodeChildList->At(j);
            ExpandNodeForGdml(curNodeChild);
        }
    }
}
//-----------------------------------------------------------------------------

//Check if Sensitive-----------------------------------------------------------
Bool_t TpcDetector::CheckIfSensitive(std::string name) {
    TString tsname = name;
    if (tsname.Contains("Active") || tsname.Contains("tpc01sv")) {
      return kTRUE;
    }
    return kFALSE;

  //  if(0 == TString(name).CompareTo("DCH1DetV")) {
  //    return kTRUE;
  //  }
  //  return kFALSE;
}
//---------------------------------------------------------

ClassImp(TpcDetector)

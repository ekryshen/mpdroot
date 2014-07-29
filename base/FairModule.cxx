// -------------------------------------------------------------------------
// -----                     FairModule source file                    -----
// -----            Created 06/01/04  by M. Al-Turany                  -----
// -------------------------------------------------------------------------
/* Generated by Together */

#include "FairModule.h"


#include "FairVolume.h"
#include "FairVolumeList.h"
#include "FairBaseParSet.h"
#include "FairRun.h"

#include "FairGeoNode.h"
#include "FairRuntimeDb.h"

#include "FairGeoInterface.h"
#include "FairGeoLoader.h"
#include "FairGeoNode.h"
#include "FairGeoRootBuilder.h"
#include "FairGeoMedia.h"


#include "TString.h"
#include "TObjArray.h"
#include "TGeoVolume.h"
#include "TFile.h"
#include "TList.h"
#include "TKey.h"
#include "TGeoManager.h"
#include "TGeoVoxelFinder.h"
#include "TGeoMatrix.h"


#include "TSystem.h"
#include <fstream>
#include <iostream>


TArrayI* FairModule::volNumber=0;
Int_t FairModule::fNbOfVolumes=0;
FairVolumeList*  FairModule::vList=0;
TRefArray*    FairModule::svList=0;



//__________________________________________________________________________
void FairModule::ConstructGeometry()
{
  fLogger->Warning(MESSAGE_ORIGIN,"FairModule::ConstructGeometry() : this method has to be implimented in detector class");

}
//__________________________________________________________________________
void FairModule::ConstructOpGeometry()
{
  fLogger->Debug2(MESSAGE_ORIGIN,"FairModule::ConstructOpGeometry : this method has to be implimented in detector class");

}
//__________________________________________________________________________
FairModule::~FairModule()
{

}
//__________________________________________________________________________
FairModule::FairModule(const char* Name, const char* title ,Bool_t Active)
  :TNamed(Name, title),
   fMotherVolumeName(""),
   fgeoVer("Not defined"),
   fgeoName("Not defined"),
   fModId(-1),
   fActive(Active),
   fNbOfSensitiveVol(0),
   fVerboseLevel(0),
   flGeoPar(0),
   kGeoSaved(kFALSE),
   fLogger(FairLogger::GetLogger())
{
  if(!svList) { svList=new TRefArray(); }
  if(!vList) { vList=new FairVolumeList(); }

}

//__________________________________________________________________________

FairModule::FairModule()
  : TNamed(),
    fMotherVolumeName(""),
    fgeoVer("Not defined"),
    fgeoName("Not defined"),
    fModId(-1),
    fActive(kFALSE),
    fNbOfSensitiveVol(0),
    fVerboseLevel(0),
    flGeoPar(0),
    kGeoSaved(kFALSE),
    fLogger(FairLogger::GetLogger())
{

}

//__________________________________________________________________________
void FairModule::Streamer(TBuffer& b)
{
  TNamed::Streamer(b);


  if (b.IsReading()) {
    fgeoVer.Streamer(b);
    fgeoName.Streamer(b);
    b >> fActive;
    b >> fModId;
  } else {
    fgeoVer.Streamer(b);
    fgeoName.Streamer(b);
    b << fActive;
    b << fModId;
  }

}
//__________________________________________________________________________
void FairModule::SetGeometryFileName(TString fname, TString geoVer)
{
  fgeoVer=geoVer;
  TString FileName = fname;
  TString work = getenv("VMCWORKDIR");
  TString userwork = getenv("GEOMPATH");
  if(userwork != "") {
    fgeoName=userwork;
    if (!fgeoName.EndsWith("/")) { fgeoName+="/"; }
    //fgeoName+=fname;
    if (TString(gSystem->FindFile(fgeoName.Data(),fname)) != TString("")) {
      fgeoName=fname;
      fLogger->Info(MESSAGE_ORIGIN, "User path for detector geometry : %s ", fgeoName.Data());
    } else {
      fLogger->Warning(MESSAGE_ORIGIN, "Detector geometry was not found in user path : %s ", FileName.Data());
      fgeoName=work+"/geometry/";
      fLogger->Info(MESSAGE_ORIGIN, "Try the standard path : %s ",fgeoName.Data());
      if (TString(gSystem->FindFile(fgeoName.Data(),FileName)) != TString("")) {
        fgeoName=FileName;
        fLogger->Info(MESSAGE_ORIGIN, "Reading detector geometry from : %s ", FileName.Data());
      } else {
        fLogger->Fatal(MESSAGE_ORIGIN,"Detector geometry not found.");
      }
    }
  } else {
    fgeoName=work+"/geometry/";
    fgeoName+=fname;
  }
}
//__________________________________________________________________________
void FairModule::ProcessNodes(TList* aList)
{

  TListIter iter(aList);
  FairGeoNode* node   = NULL;
  FairGeoNode* MotherNode =NULL;
  FairVolume*  volume = NULL;
  FairRuntimeDb* rtdb= FairRun::Instance()->GetRuntimeDb();
  FairBaseParSet* par=(FairBaseParSet*)(rtdb->getContainer("FairBaseParSet"));
  TObjArray* fNodes = par->GetGeoNodes();
  while( (node = (FairGeoNode*)iter.Next()) ) {

    node->calcLabTransform();
    MotherNode=node->getMotherNode();
    volume = new FairVolume( node->getTruncName(), fNbOfVolumes++);
    volume->setRealName(node->GetName());
    vList->addVolume(volume);
    volume->setGeoNode(node);
    volume->setCopyNo(  node->getCopyNo());

    if(MotherNode!=0) {
      volume->setMotherId(node->getMCid());
      volume->setMotherCopyNo(MotherNode->getCopyNo());
    }
    FairGeoVolume* aVol=NULL;

    if ( node->isSensitive() && fActive ) {
      volume->setModId(fModId);
      volume->SetModule(this);
      svList->Add(volume);
      aVol = dynamic_cast<FairGeoVolume*> ( node );
      fNodes->AddLast( aVol );
      fNbOfSensitiveVol++;
    }

  }
}
//__________________________________________________________________________
void  FairModule::AddSensitiveVolume(TGeoVolume* v)
{

  fLogger->Debug2(MESSAGE_ORIGIN, "FairModule::AddSensitiveVolume", v->GetName());
  FairVolume*  volume = NULL;
  volume = new FairVolume(v->GetName(), fNbOfVolumes++);
  vList->addVolume(volume);
  volume->setModId(fModId);
  volume->SetModule(this);
  svList->Add(volume);
  fNbOfSensitiveVol++;

}


//__________________________________________________________________________
FairVolume* FairModule::getFairVolume(FairGeoNode* fN)
{
  FairVolume* fv;
  FairVolume*  fvol=0;
  for(Int_t i=0; i<vList->getEntries(); i++) {
    fv=vList->At(i);
    if((fv->getGeoNode())==fN) {
      fvol=fv;
      return fvol;
      break;
    }
  }
  return fvol;
}
//__________________________________________________________________________
void FairModule::ConstructRootGeometry()
{
  /** Construct the detector geometry from ROOT files, possible inputs are:
   * 1. A TGeoVolume as a mother (master) volume containing the detector geometry
   * 2. A TGeoManager with the detector geometry
   * 3. A TGeoVolume as a mother or Master volume which is the output of the CAD2ROOT geometry, in this case
   *    the materials are not proprely defined and had to be reset
   *  In all cases we have to check that the material properties are the same or is the materials defined in
   *  the current simulation session
   */
  TGeoManager* OldGeo=gGeoManager;
  TGeoManager* NewGeo=0;
  TGeoVolume* volume=0;;
  TFile* f=new TFile(GetGeometryFileName().Data());
  TList* l= f->GetListOfKeys();
  TKey* key;
  TIter next( l);
  TGeoNode* n=0;
  TGeoVolume* v1=0;
  while ((key = (TKey*)next())) {
    /**loop inside the delivered root file and try to fine a TGeoManager object
     * the first TGeoManager found will be read
     */
    if (strcmp(key->GetClassName(),"TGeoManager") != 0) { continue; }
    gGeoManager=0;
    NewGeo = (TGeoManager*)key->ReadObj();
    break;
  }
  if (NewGeo!=0) {
    /** in case a TGeoManager was found get the top most volume and the node
     */

    NewGeo->cd();
    volume=(TGeoVolume*)NewGeo->GetNode(0)->GetDaughter(0)->GetVolume();
    v1=volume->MakeCopyVolume(volume->GetShape());
    // n=NewGeo->GetTopNode();
    n=v1->GetNode(0);
    //  NewGeo=0;
    // delete NewGeo; //move it to the end of the method

  } else {
    /** The file does not contain any TGeoManager, so we assume to have a file with a TGeoVolume
     * try to look for a TGeoVolume inside the file
     */

    key=(TKey*) l->At(0);  //Get the first key in the list
    volume=dynamic_cast<TGeoVolume*> (key->ReadObj());
    if(volume!=0) { n=volume->GetNode(0); }
    if(n!=0) { v1=n->GetVolume(); }
  }

  if(v1==0) {
    fLogger->Fatal(MESSAGE_ORIGIN, "\033[5m\033[31mFairModule::ConstructRootGeometry(): could not find any geometry in File!!  \033[0m", GetGeometryFileName().Data());
  }
  gGeoManager=OldGeo;
  gGeoManager->cd();
  // If AddToVolume is empty add the volume to the top volume Cave
  // If it is defined check i´f the volume exists and if it exists add the volume from the root file
  // to the already existing volume
  TGeoVolume* Cave=NULL;
  if ( 0 == fMotherVolumeName.Length() ) {
    Cave= gGeoManager->GetTopVolume();
  } else {
    Cave = gGeoManager->GetVolume(fMotherVolumeName);
    if ( NULL == Cave ) {
      fLogger->Fatal(MESSAGE_ORIGIN,"\033[5m\033[31mFairModule::ConstructRootGeometry(): could not find the given mother volume \033[0m   %s \033[5m\033[31m where the geomanger should be added. \033[0m", fMotherVolumeName.Data());
    }
  }
  /**Every thing is OK, we have a TGeoVolume and now we add it to the simulation TGeoManager  */
  gGeoManager->AddVolume(v1);
  /** force rebuilding of voxels */
  TGeoVoxelFinder* voxels = v1->GetVoxels();
  if (voxels) { voxels->SetNeedRebuild(); }
  /**To avoid having different names of the default matrices because we could have get the volume from another
   * TGeoManager, we reset the default matrix name
   */
  TGeoMatrix* M = n->GetMatrix();
  SetDefaultMatrixName(M);

  /** NOw we can remove the matrix so that the new geomanager will rebuild it properly*/
  gGeoManager->GetListOfMatrices()->Remove(M);
  TGeoHMatrix* global = gGeoManager->GetHMatrix();
  gGeoManager->GetListOfMatrices()->Remove(global); //Remove the Identity matrix
  /**Now we can add the node to the existing cave */
  Cave->AddNode(v1,0, M);
  /** correction from O. Merle: in case of a TGeoVolume (v1) set the material properly */
  AssignMediumAtImport(v1);
  /** now go through the herachy and set the materials properly, this is important becase the CAD converter
   *  produce TGeoVolumes with materials that have only names and no properties
   */
  ExpandNode(n);
  if(NewGeo!=0) { delete NewGeo; }
  delete f;
}
//__________________________________________________________________________
void FairModule::ConstructASCIIGeometry()
{
  fLogger->Warning(MESSAGE_ORIGIN,"FairModule::ConstructASCIIGeometry() : this method has to be implimented in detector class");

}
//__________________________________________________________________________
Bool_t FairModule::CheckIfSensitive(std::string name)
{
  fLogger->Warning(MESSAGE_ORIGIN,"FairModule::CheckIfSensitive(std::string name): this method has to be implimented in detector class");
  return kFALSE;
}
//__________________________________________________________________________
void FairModule::ExpandNode(TGeoNode* fN)
{
  //FairGeoLoader* geoLoad = FairGeoLoader::Instance();
  //FairGeoInterface* geoFace = geoLoad->getGeoInterface();
  //FairGeoMedia* Media =  geoFace->getMedia();
  //FairGeoBuilder* geobuild=geoLoad->getGeoBuilder();
  TGeoMatrix* Matrix =fN->GetMatrix();
  if(gGeoManager->GetListOfMatrices()->FindObject(Matrix)) { gGeoManager->GetListOfMatrices()->Remove(Matrix); }
  TGeoVolume* v1=fN->GetVolume();
  TObjArray* NodeList=v1->GetNodes();
  for (Int_t Nod=0; Nod<NodeList->GetEntriesFast(); Nod++) {
    TGeoNode* fNode =(TGeoNode*)NodeList->At(Nod);
    TGeoMatrix* M =fNode->GetMatrix();
    //M->SetDefaultName();
    SetDefaultMatrixName(M);
    if(fNode->GetNdaughters()>0) { ExpandNode(fNode); }
    TGeoVolume* v= fNode->GetVolume();
    AssignMediumAtImport(v);
    if (!gGeoManager->FindVolumeFast(v->GetName())) {
      fLogger->Debug2(MESSAGE_ORIGIN,"Register Volume : %s ", v->GetName());
      v->RegisterYourself();
    }
    if (CheckIfSensitive(v->GetName())) {
      fLogger->Debug2(MESSAGE_ORIGIN,"Sensitive Volume : %s ", v->GetName());
      AddSensitiveVolume(v);
    }
  }
}
//__________________________________________________________________________
void FairModule::SetDefaultMatrixName(TGeoMatrix* matrix)
{
  // Copied from root TGeoMatrix::SetDefaultName() and modified (memory leak)
  // If no name was supplied in the ctor, the type of transformation is checked.
  // A letter will be prepended to the name :
  //   t - translation
  //   r - rotation
  //   s - scale
  //   c - combi (translation + rotation)
  //   g - general (tr+rot+scale)
  // The index of the transformation in gGeoManager list of transformations will
  // be appended.
  if (!gGeoManager) { return; }
  if (strlen(matrix->GetName())) { return; }
  char type = 'n';
  if (matrix->IsTranslation()) { type = 't'; }
  if (matrix->IsRotation()) { type = 'r'; }
  if (matrix->IsScale()) { type = 's'; }
  if (matrix->IsCombi()) { type = 'c'; }
  if (matrix->IsGeneral()) { type = 'g'; }
  TObjArray* matrices = gGeoManager->GetListOfMatrices();
  Int_t index = 0;
  if (matrices) { index =matrices->GetEntriesFast() - 1; }
  matrix->SetName(Form("%c%i", type, index));
}

//__________________________________________________________________________

void FairModule::AssignMediumAtImport(TGeoVolume* v)
{

  /**
   * Assign medium to the the volume v, this has to be done in all cases:
   * case 1: For CAD converted volumes they have no mediums (only names)
   * case 2: TGeoVolumes, we need to be sure that the material is defined in this session
   */
  FairGeoMedia* Media       = FairGeoLoader::Instance()->getGeoInterface()->getMedia();
  FairGeoBuilder* geobuild  = FairGeoLoader::Instance()->getGeoBuilder();

  TGeoMedium* med1=v->GetMedium();
  if(med1) {
    TGeoMaterial* mat1=v->GetMaterial();
    TGeoMaterial* newMat = gGeoManager->GetMaterial(mat1->GetName());
    if( newMat==0) {
      /**The Material is not defined in the TGeoManager, we try to create one if we have enough information about it*/
      FairGeoMedium* FairMedium=Media->getMedium(mat1->GetName());
      if (!FairMedium) {
        fLogger->Fatal(MESSAGE_ORIGIN,"Material %s is not defined in ASCII file nor in Root file we Stop creating geometry", mat1->GetName());
        //     FairMedium=new FairGeoMedium(mat1->GetName());
        //      Media->addMedium(FairMedium);
      }

      Int_t nmed=geobuild->createMedium(FairMedium);
      v->SetMedium(gGeoManager->GetMedium(nmed));
      gGeoManager->SetAllIndex();
    } else {
      /**Material is already available in the TGeoManager and we can set it */
      TGeoMedium* med2= gGeoManager->GetMedium(mat1->GetName());
      v->SetMedium(med2);
    }
  } else {
    if (strcmp(v->ClassName(),"TGeoVolumeAssembly") != 0) {
      //[R.K.-3.3.08]  // When there is NO material defined, set it to avoid conflicts in Geant
      fLogger->Fatal(MESSAGE_ORIGIN,"The volume  %s  Has no medium information and not an Assembly so we have to quit", v->GetName());
    }
  }
}

//__________________________________________________________________________
ClassImp(FairModule)




////////////////////////////////////////////////////////////////
//                                                            //
//  MpdEmcKI                                                  //
//  Cluster production for EMC                                //
//  Author List : D.Peresunko., RRCKI, 2019                   //
//                                                            //
////////////////////////////////////////////////////////////////
#include "MpdEmcKI.h"

#include "MpdEmcGeo.h"
#include "MpdEmcGeoPar.h"
#include "MpdEmcGeoUtils.h"
#include "MpdEmcPointKI.h"

#include "FairGeoInterface.h"
#include "FairGeoLoader.h"
#include "FairGeoNode.h"
#include "FairGeoVolume.h"
#include "FairRootManager.h"
#include "FairRun.h"
#include "FairRuntimeDb.h"
#include "MpdStack.h"
#include "FairVolume.h"
#include "MpdDetectorList.h"

#include "TClonesArray.h"
#include "TGeoManager.h"
#include "TGeoNode.h"
#include "TGeoTube.h"
#include "TParticle.h"
#include "TVirtualMC.h"

#include <iostream>
using std::cout;
using std::endl;

MpdEmcKI::MpdEmcKI()
  : FairDetector("EMC", kTRUE, kECAL),
    fNhits(0),
    fCurrentTrackID(-1),
    fCurrentCellID(-1),
    fCurentSuperParent(-1),
    fCurrentHit(nullptr),
    fGeom(nullptr),
    fEcalRmin(0),
    fEcalRmax(0)
{
}

MpdEmcKI::MpdEmcKI(const char* name, Bool_t active)
  : FairDetector(name, active, kECAL),
    fNhits(0),
    fCurrentTrackID(-1),
    fCurrentCellID(-1),
    fCurentSuperParent(-1),
    fCurrentHit(nullptr),
    fGeom(nullptr)
{
  fMpdEmcPointCollection = new TClonesArray("MpdEmcPointKI");
  fVerboseLevel = 1;
}

MpdEmcKI::~MpdEmcKI()
{
  if (fMpdEmcPointCollection) {
    fMpdEmcPointCollection->Delete();
    delete fMpdEmcPointCollection;
    fMpdEmcPointCollection = nullptr;
  }
}

void MpdEmcKI::Reset()
{
  fSuperParents.clear();
  fMpdEmcPointCollection->Clear();
  fNhits = 0;
  fCurrentTrackID = -1;
  fCurrentCellID = -1;
  fCurentSuperParent = -1;
  fCurrentHit = nullptr;
}
/*
void MpdEmcKI::Initialize()
{
  //Extract inner and outer radii of ECAL from current geometry

  TGeoVolume * v =gGeoManager->GetVolume("emc1Chamber1_0");
  TGeoTube *sh = dynamic_cast<TGeoTube*>(v->GetShape()) ;
//    double towerHalfSizeZ = sh->GetDz() ;
  fEcalRmin = sh->GetRmin() ;
  fEcalRmax = sh->GetRmax() ;
}
*/
Bool_t MpdEmcKI::ProcessHits(FairVolume* vol)
{
  // This method is called from the MC stepping
  // 1. Remember all particles first entered ECAL (active medium declared below: _sc* and _pl*)
  // 2. If this is _pl, do not create hit, just remember entered particle
  // 3. Collect all energy depositions in sc* by all secondaries from particle first entered ECAL

  // Check if this is first entered ECAL particle ("SuperParent")
  MpdStack* stack = static_cast<MpdStack*>(gMC->GetStack());
  const Int_t partID = stack->GetCurrentTrackNumber();
  Int_t superParent = -1;
  Bool_t isNewPartile = false;     // Create Hit even if zero energy deposition
  if (partID != fCurrentTrackID) { // not same track as before, check: same SuperParent or new one?

    auto itTr = fSuperParents.find(partID);
    if (itTr == fSuperParents.end()) {
      // Search parent
      Int_t parentID = stack->GetCurrentTrack()->GetMother(0);
      itTr = fSuperParents.find(parentID);
      if (itTr == fSuperParents.end()) { // Neither track or its parent found: new SuperParent
        // Go back to parents untill found particle entered ECAL
        TParticle* part = stack->GetCurrentTrack();
        Double_t rvert = part->R(); // Transverse radius despite the name!
        TParticle* mctr = part;
        Int_t idmother = partID;
        while (rvert > fEcalRmin && rvert < fEcalRmax) {
          // Born inside EMC - find ancestor
          Int_t tmp = mctr->GetMother(0);
          if (tmp < 0) {
            break;
          }
          idmother = tmp;
          mctr = stack->GetParticle(idmother);
          rvert = mctr->R(); // transverse radius despite the name!!!
        }

        fSuperParents[partID] = idmother;
        if (partID != idmother) { // Add intermediate tracks if necessary
          Int_t tmp = part->GetMother(0);
          while (tmp >= 0 && tmp != idmother) {
            itTr = fSuperParents.find(tmp);
            if (itTr == fSuperParents.end()) {
              fSuperParents[tmp] = idmother;
            } else {
              break;
            }
            mctr = stack->GetParticle(tmp);
            tmp = mctr->GetMother(0);
          }
        }

        superParent = idmother;
        isNewPartile = true;
        // Mark current track to keep in stack
        stack->AddPoint(kECAL, superParent);
      } else { // parent found, this track - not
        superParent = itTr->second;
        fSuperParents[partID] = superParent;
        fCurrentTrackID = partID;
      }
    } else {
      superParent = itTr->second;
      fCurrentTrackID = partID;
    }
  } else { // Same track as before
    superParent = fCurentSuperParent;
  }

  Double_t lostenergy = gMC->Edep();
  if (strstr(vol->GetName(), "cl_pl") != 0) { // Do not count energy deposited in cover
    lostenergy = 0.;
  }

  if (lostenergy < DBL_EPSILON && !isNewPartile) {
    fCurentSuperParent = superParent; // If we switched to another track/superparent, remember it
    return false;                     // do not create hits with zero energy deposition
  }

  if (!fGeom) {
    fGeom = MpdEmcGeoUtils::GetInstance();
  }

//printf("GEANT path: <%s>\n",gMC->CurrentVolPath()) ;


  Int_t chamberH, chamber, sector, crate, module, boxA, boxB, dummyA, dummyB ;
  char dummyC1,dummyC2;
/*  gMC->CurrentVolOffID(5,
                       chamberH); // 5: number of geom. levels between emc1_cl_sc* and sector: get the sector number ;
  Int_t chamber;
  gMC->CurrentVolOffID(4, chamber); // 4: number of geom. levels between emc1_cl_sc* and sector: get the sector number ;
  Int_t sector;
  gMC->CurrentVolOffID(3, sector); // 3: number of geom. levels between emc1_cl_sc* and sector: get the sector number ;
  Int_t crate;
  gMC->CurrentVolOffID(2, crate); // 2: Crate in sector: number of geom levels between emc1_cl_sc get the crate number
  Int_t box;
  gMC->CurrentVolOffID(1, box); // 1: Tower in crate: number of geom levels between emc1_cl_sc and box
*/
  //Example path: </cave_1/emcChamber_0/emcChH_1/emcSector_20/emcCrate_1/emcModule6_0/emc_box50_99/emc_cl_sc81_161>
  // /cave_1/emcChamber_0/emcChH_1/emcSector_0/emcCrate_5/emcModule6_0/emc_box53_105/emc_cl_pl1_0
  int nRead= sscanf(gMC->CurrentVolPath(),"/cave_1/emcChamber_%d/emcChH_%d/emcSector_%d/emcCrate_%d/emcModule%d_0/emc_box%d_%d/emc_cl_%c%c%d_%d",
         &chamberH,&chamber,&sector,&crate,&module,&boxA,&boxB, &dummyC1,&dummyC2,&dummyA, &dummyB) ;
  if(nRead!=11){
    LOG(FATAL)<<"Can not parse Geant path:" << gMC->CurrentVolPath() << " nRead = " << nRead<<endl ;
    return false ;
  }

  Int_t detID = fGeom->GeantToDetId(chamberH, chamber, sector, crate, module, boxA, boxB);

  if (superParent == fCurentSuperParent && detID == fCurrentCellID && fCurrentHit) {
    // continue with current hit
    fCurrentHit->AddEnergyLoss(lostenergy); // TODO implement light attenuation vs. scinitllator number
    return true;
  }

  // try to find existing Hit and add energy to it
  if (!isNewPartile) {
    for (Int_t itr = fNhits - 1; itr >= 0; itr--) {
      MpdEmcPointKI* h = static_cast<MpdEmcPointKI*>(fMpdEmcPointCollection->At(itr));
      if (h->GetTrackID() != superParent) // switched to another SuperParent, do not search further
        break;
      if (h->GetDetectorID() == detID) { // found correct hit
        h->AddEnergyLoss(lostenergy);
        fCurentSuperParent = superParent;
        fCurrentTrackID = partID;
        fCurrentCellID = detID;
        fCurrentHit = h;
        return true;
      }
    }
  }
  // Create new Hit
  Float_t posX = 0., posY = 0., posZ = 0., momX = 0, momY = 0., momZ = 0., energy = 0.;
  gMC->TrackPosition(posX, posY, posZ);
  gMC->TrackMomentum(momX, momY, momZ, energy);
  Double_t time = gMC->TrackTime() * 1.e+9; // time in ns?? To be consistent with EMCAL
  Float_t length = gMC->TrackLength();
  fCurrentHit =
    AddHit(superParent, detID, TVector3(posX, posY, posZ), TVector3(momX, momY, momZ), time, length, lostenergy);
  fCurentSuperParent = superParent;
  fCurrentTrackID = partID;
  fCurrentCellID = detID;

  return true;
}

void MpdEmcKI::EndOfEvent()
{
  Print();

  Reset();
}

void MpdEmcKI::Register()
{
  /** This will create a branch in the output tree called
      MpdEmcPointKI, setting the last parameter to kFALSE means:
      this collection will not be written to the file, it will exist
      only during the simulation.
  */

  FairRootManager::Instance()->Register("EmcPoint", "MpdEmc", fMpdEmcPointCollection, kTRUE);
}

TClonesArray* MpdEmcKI::GetCollection(Int_t iColl) const
{
  if (iColl == 0)
    return fMpdEmcPointCollection;
  else
    return NULL;
}
void MpdEmcKI::FinishEvent()
{
  // Sort Hits
  // Add duplicates if any and remove them
  if (fNhits == 0) {
    return;
  }

  fMpdEmcPointCollection->Sort();

  for (Int_t i = 0; i < fNhits; i++) {
    MpdEmcPointKI* firstHit = static_cast<MpdEmcPointKI*>(fMpdEmcPointCollection->At(i));
    if (!firstHit) { // hit already removed
      continue;
    }
    for (Int_t j = i + 1; j < fNhits; j++) {
      MpdEmcPointKI* secondHit = static_cast<MpdEmcPointKI*>(fMpdEmcPointCollection->At(j));
      if ((*firstHit) == (*secondHit)) { // same volume and same parent: add and remove second
        (*firstHit) += (*secondHit);
        fMpdEmcPointCollection->RemoveAt(j);
      } else { // hits are sorted, no need to scan further
        break;
      }
    }
    // Clean, do not store empty points
    if (firstHit->GetEnergyLoss() == 0) {
      fMpdEmcPointCollection->RemoveAt(i);
    }
  }

  // Remove empty slots and shink array
  fMpdEmcPointCollection->Compress();
  fNhits = fMpdEmcPointCollection->GetEntriesFast();
  fMpdEmcPointCollection->Expand(fNhits);
}

void MpdEmcKI::Print() const
{
  Int_t nHits = fMpdEmcPointCollection->GetEntriesFast();
  cout << "-I- MpdEmcKI: " << nHits << " points registered in this event." << endl;

  if (fVerboseLevel > 1)
    for (Int_t i = 0; i < nHits; i++)
      (*fMpdEmcPointCollection)[i]->Print();
}

void MpdEmcKI::ConstructGeometry()
{
  TString fileName = GetGeometryFileName();

  if (fileName.EndsWith(".root")) {
    LOG(INFO) << "Constructing EMC geometry from ROOT file " << fileName.Data();
    ConstructRootGeometry();
  } else if (fileName.EndsWith(".geo")) {
    LOG(INFO) << "Constructing EMC geometry from ASCII file " << fileName.Data();
    ConstructAsciiGeometry();
  } else {
    LOG(FATAL) << "Geometry format of EMC file " << fileName.Data() << " not supported.";
  }
  // Extract inner and outer radii of ECAL from current geometry
//TODO!!!!!!!!!
/*  TGeoVolume* v = gGeoManager->GetVolume("emc1Chamber1");
  TGeoTube* sh = dynamic_cast<TGeoTube*>(v->GetShape());
  //    double towerHalfSizeZ = sh->GetDz() ;
  fEcalRmin = sh->GetRmin();
  fEcalRmax = sh->GetRmax();
*/}

void MpdEmcKI::ConstructAsciiGeometry()
{
  FairGeoLoader* geoLoad = FairGeoLoader::Instance();
  FairGeoInterface* geoFace = geoLoad->getGeoInterface();
  MpdEmcGeo* Geo = new MpdEmcGeo();
  Geo->setGeomFile(GetGeometryFileName());
  geoFace->addGeoModule(Geo);

  Bool_t rc = geoFace->readSet(Geo);
  if (rc)
    Geo->create(geoLoad->getGeoBuilder());
  TList* volList = Geo->getListOfVolumes();

  // store geo parameter
  FairRun* fRun = FairRun::Instance();
  FairRuntimeDb* rtdb = FairRun::Instance()->GetRuntimeDb();
  MpdEmcGeoPar* par = (MpdEmcGeoPar*)(rtdb->getContainer("MpdEmcGeoPar"));
  TObjArray* fSensNodes = par->GetGeoSensitiveNodes();
  TObjArray* fPassNodes = par->GetGeoPassiveNodes();

  TListIter iter(volList);
  FairGeoNode* node = NULL;
  FairGeoVolume* aVol = NULL;

  while ((node = (FairGeoNode*)iter.Next())) {
    aVol = dynamic_cast<FairGeoVolume*>(node);
    if (node->isSensitive()) {
      fSensNodes->AddLast(aVol);
    } else {
      fPassNodes->AddLast(aVol);
    }
  }
  par->setChanged();
  par->setInputVersion(fRun->GetRunId(), 1);

  ProcessNodes(volList);
}

// Check sensitivity

Bool_t MpdEmcKI::CheckIfSensitive(std::string name)
{
  return (name.find("cl_sc") != std::string::npos || // Active scinillator
          name.find("cl_pl") != std::string::npos);  // forward plate to remember primary track
}

MpdEmcPointKI* MpdEmcKI::AddHit(Int_t trackID, Int_t detID, TVector3 pos, TVector3 mom, Double_t time, Double_t length,
                            Double_t ELoss)
{
  return new ((*fMpdEmcPointCollection)[fNhits++]) MpdEmcPointKI(trackID, detID, pos, mom, time, length, ELoss);
}

ClassImp(MpdEmcKI)

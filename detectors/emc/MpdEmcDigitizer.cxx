//--------------------------------------------------------------------
//
// Description:
//      MPD EMC Digitizer - takes EmcPoints and makes digits
//
//
// Author List:
//      Alexander Zinchenko LHEP, JINR, Dubna - 8-May-2016
//      Alexander Zinchenko LHEP, JINR, Dubna - 24-June-2018 - adapted to projective geometry
//
//--------------------------------------------------------------------

#include "MpdEmcDigitizer.h"
#include "MpdEmcGeoParams.h"
#include "MpdEmcDigit.h"
#include "MpdEmcPoint.h"

#include "FairMCPoint.h"
#include "MpdMCTrack.h"
#include "FairRootManager.h"
#include "FairRunAna.h"

#include <TClonesArray.h>
#include <TGeoBBox.h>
#include <TGeoManager.h>
#include <TMath.h>

using namespace std;
using namespace TMath;

FILE *lun = fopen ("file.txt","w"); //AZ

// -----   Default constructor   -------------------------------------------

MpdEmcDigitizer::MpdEmcDigitizer() :
FairTask("EMC digitizer") {

}

// -----   Destructor   ----------------------------------------------------

MpdEmcDigitizer::~MpdEmcDigitizer() {
}
// -------------------------------------------------------------------------

// -----   Public method Init   --------------------------------------------

InitStatus MpdEmcDigitizer::Init() {

    cout << "\n******************* EMC Digitizer INIT *********************" << endl;

    // Get RootManager
    FairRootManager* ioman = FairRootManager::Instance();
    if (!ioman) {
        cout << "-E- MpdEmcDigitizer::Init: "
                << "RootManager not instantiated!" << endl;
        return kFATAL;
    }

    fGeoPar = new MpdEmcGeoParams();

    // Get input array
    fPointArray = (TClonesArray*) ioman->GetObject("EmcPoint");
    if (!fPointArray) {
        cout << "-W- MpdEmcDigitizer::Init: " << "No EmcPoint array!" << endl;
        return kERROR;
    }
    
    fMcTrackArray = (TClonesArray*) ioman->GetObject("MCTrack");
    if (!fMcTrackArray) {
        cout << "-W- MpdEmcDigitizer::Init: " << "No MCTrack array!" << endl;
        return kERROR;
    }

    // Create and register output array
    fDigiArray = new TClonesArray("MpdEmcDigit", 100);
    ioman->Register("EmcDigit", "EMC", fDigiArray, kTRUE);

    cout << "-I- MpdEmcDigitizer: Intialization successfull" << endl;

    return kSUCCESS;

}

//__________________________________________________________________________

void MpdEmcDigitizer::Finish() {
    cout << "-I- MpdEmcDigitizer: Finish" << endl;
}

//__________________________________________________________________________

void MpdEmcDigitizer::Exec(Option_t* opt) 
{
  // Main processing engine

  cout << "MpdEmcDigitizer::Exec started" << endl;

  static Int_t nSecRows = -1;
  if (nSecRows < 0) 
    nSecRows = fGeoPar->GetNrows() / (fGeoPar->GetNsec() / 2 - 1); // 6 wide, 2 narrow sectors

  // Reset output Array
  if (!fDigiArray) Fatal("MpdEmcDigitizer::Exec", "No array of digits");
  fDigiArray->Delete();
  fHitMap.clear();
  Int_t nPoints = fPointArray->GetEntriesFast();

  for (Int_t iPnt = 0; iPnt < nPoints; ++iPnt) {
    FairMCPoint* pnt = (FairMCPoint*) fPointArray->UncheckedAt(iPnt);
    Int_t trId = pnt->GetTrackID();
    if (trId < 0) continue; // strange case
    MpdMCTrack* tr = (MpdMCTrack*) fMcTrackArray->UncheckedAt(trId);
    //        if (tr->GetMotherId() != -1) continue;
    Int_t pdg = tr->GetPdgCode();
    Double_t x = pnt->GetX();
    Double_t y = pnt->GetY();
    Double_t z = pnt->GetZ();
    Float_t e = pnt->GetEnergyLoss();

    gGeoManager->FindNode(x,y,z); 
    TString path(gGeoManager->GetPath());
    //cout << path << endl; 
    //if (!path.Contains("sc")) exit(0); 
    //if (!path.Contains("sc")) continue; // !!! this condition is very bad for Geant4 (not for Geant3)
    if (!path.Contains("box")) continue; // 
    if (path.Contains("_cl_")) gGeoManager->CdUp(); // to tower (channel)
    Int_t ip = path.Index("ChH"); // half-EMC
    UInt_t ih = TString(path(ip+4,1)).Atoi();
    ip = path.Index("Sector"); // sector
    Int_t ip1 = path.Index("/",ip);
    UInt_t isec = TString(path(ip+8,ip1)).Atoi();
    if (path(ip+6) == '2') {
      // Narrow sectors
      if (isec == 0) isec = 2;
      else isec = 6;
    } else {
      if (isec >= 5) isec += 2;
      else if (isec >= 2) ++isec;
    }
    ip = path.Index("Crate"); // row
    ip1 = path.Index("/",ip);
    UInt_t irow = TString(path(ip+6,ip1)).Atoi();
    ip = path.Index("box");
    if (path(ip+4) == '_') --ip; // 1-digit number
    ip1 = path.Index("/",ip);
    UInt_t ibox = TString(path(ip+6,ip1)).Atoi();
    //cout << ih << " " << isec << " " << irow << " " << ibox << endl;
    fprintf(lun,"%d %d %d %d %f %f %f %f\n",ih,isec,irow,ibox,TMath::ATan2(y,x),x,y,z);
    // Code channel ID
    TString tower = "c"; // channel
    tower += ibox; 
    tower += "r"; // row
    tower += irow;
    tower += "s"; // sector
    tower += isec;
    tower += "h"; // EMC half (Z <> 0)
    tower += ih;
    
    MpdEmcDigit* hit = SearchHit(tower);
    if (hit == NULL) {
      hit = new((*fDigiArray)[fDigiArray->GetEntriesFast()]) MpdEmcDigit(ih, isec, irow, ibox);
      fHitMap[tower] = hit;
      Double_t phi, the;
      //FindChanPhiZ(phi,z);
      FindChanPhiThe(phi,the);
      hit->SetPhiCenter(phi);
      hit->SetZCenter(the);
      Int_t iz = ibox;
      if (ih == 1) ++iz; // negative Z
      hit->SetChanZId(TMath::Abs(iz));
      Int_t iphi = irow + isec * nSecRows;
      if (isec >= 7) iphi -= nSecRows; // take into account narrow sector
      else if (isec >= 3) iphi -= nSecRows / 2; // take into account narrow sector
      hit->SetChanPhiId(iphi);
      hit->SetTimeStamp(pnt->GetTime());
    } 
    hit->IncreaseEnergy(e, trId);
    hit->SetTimeStamp(TMath::Min(pnt->GetTime(),hit->GetTimeStamp()));
    hit->SetNumTracks(hit->GetNumTracks() + 1);
    //hit->SetTrackId(trId);
    hit->SetPdg(pdg);
    //AZ hit->SetZCenter(CalcZCenter(sec, row, mod));
    //AZ hit->SetPhiCenter(CalcPhiCenter(sec, supMod, mod));

    // Redo timing - take time from the highest energy deposit
    if (TMath::Abs(e-hit->GetDy()) < 1.e-6) {
      hit->SetTimeStamp(pnt->GetTime()); // Dy - holder of Emax
      hit->SetTrackId(iPnt); // point index
    }
  }
  
  Int_t nHits = fDigiArray->GetEntriesFast();
  for (Int_t iHit = 0; iHit < nHits; ++iHit) {
    MpdEmcDigit* hit = (MpdEmcDigit*) fDigiArray->UncheckedAt(iHit);
    if (hit->GetNumTracks() > 1) {
      hit->SetPdg(0);
      //hit->SetTrackId(-1);
    }
    //hit->Print();
  }    

  // Redo track ID numbering 
  //RedoId();
}

//__________________________________________________________________________

void MpdEmcDigitizer::RedoId(TClonesArray *digis, TClonesArray *mctrs)
{
  // Redo track ID numbering 
  // Take IDs of particles produced outside EMC

  static const Double_t rmin = fGeoPar->GetRmin(), rmax = fGeoPar->GetRmax();
  TVector3 vert;
  Int_t nDigis = digis->GetEntriesFast();

  for (Int_t iDig = 0; iDig < nDigis; ++iDig) {
    MpdEmcDigit* digi = (MpdEmcDigit*) digis->UncheckedAt(iDig);
    map<Int_t,Float_t> contrib = digi->GetContrib();
    map<Int_t,Float_t> copy(contrib);
    contrib.clear();

    for (map<Int_t,Float_t>::iterator it = copy.begin(); it != copy.end(); ++it) {
      Int_t id = it->first, idm = -1;
      MpdMCTrack *mctr = (MpdMCTrack*) mctrs->UncheckedAt(id);
      mctr->GetStartVertex(vert);
      Double_t rvert = vert.Pt();
      while (rvert > rmin && rvert < rmax) {
	// Born inside EMC - find ancestor
	idm = mctr->GetMotherId();
	if (idm < 0) break;
	mctr = (MpdMCTrack*) mctrs->UncheckedAt(idm);
	mctr->GetStartVertex(vert);
	rvert = vert.Pt();
      }
      if (idm >= 0) id = idm;
      if (contrib.find(id) == contrib.end()) contrib[id] = it->second;
      else contrib[id] += it->second;
    }
  }    
}

//__________________________________________________________________________

void MpdEmcDigitizer::FindChanPhiZ(Double_t &phi, Double_t &z) 
{
  // Compute Z-position and Phi-angle of the tower (channel) center (at inner radius)

  TString path(gGeoManager->GetPath());
  Int_t ip = path.Last('/');
  Int_t ip1 = path.Last('_');
  TString volName = path(ip+1,ip1-ip-1);
  TGeoVolume *tower = gGeoManager->GetVolume(volName);
  TGeoBBox *box = (TGeoBBox*) tower->GetShape();
  Double_t xyzL[3] = {0,0,-box->GetDZ()}, xyzM[3];
  gGeoManager->LocalToMaster(xyzL,xyzM);
  z = xyzM[2];
  phi = TMath::ATan2 (xyzM[1],xyzM[0]);
  if (phi < 0) phi += TMath::TwoPi();
}

//__________________________________________________________________________

void MpdEmcDigitizer::FindChanPhiThe(Double_t &phi, Double_t &the) 
{
  // Compute Phi and Theta angles of the tower (channel) center (at inner radius)

  TString path(gGeoManager->GetPath());
  Int_t ip = path.Last('/');
  Int_t ip1 = path.Last('_');
  TString volName = path(ip+1,ip1-ip-1);
  TGeoVolume *tower = gGeoManager->GetVolume(volName);
  TGeoBBox *box = (TGeoBBox*) tower->GetShape();
  Double_t xyzL[3] = {0,0,-box->GetDZ()}, xyzM[3];
  gGeoManager->LocalToMaster(xyzL,xyzM);
  phi = TMath::ATan2 (xyzM[1],xyzM[0]);
  if (phi < 0) phi += TMath::TwoPi();
  the = TMath::ATan2 (TMath::Sqrt(xyzM[0]*xyzM[0]+xyzM[1]*xyzM[1]), xyzM[2]);
}

//__________________________________________________________________________

MpdEmcDigit* MpdEmcDigitizer::SearchHit(TString tower)
{
  // Check if the hit (tower) exists in the map

  if (fHitMap.find(tower) == fHitMap.end()) return NULL;
  return fHitMap[tower];
}

// -------------------------------------------------------------------------

ClassImp(MpdEmcDigitizer)

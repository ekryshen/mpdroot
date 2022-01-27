////////////////////////////////////////////////////////////////
//                                                            //
//  MpdEmcHitCreation                                         //
//  Hit production for EMC, v03                        	      //
//  Author List : Martemianov M., 2021                        //
//                                                            //
////////////////////////////////////////////////////////////////

#include "TClonesArray.h"
#include "FairRootManager.h"
#include "TMath.h"
#include "TGeoNode.h"
#include "TGeoManager.h"
#include "FairRunAna.h"
#include "FairEventHeader.h"

#include "MpdEmcHitCreation.h"
#include "MpdEmcGeoParams.h"
#include "MpdEmcHit.h"
#include "MpdEmcPoint.h"
#include "FairMCPoint.h"
#include "MpdMCTrack.h"

using namespace std;
using namespace TMath;

// -----   Default constructor   -------------------------------------------

MpdEmcHitCreation::MpdEmcHitCreation() : FairTask("Ideal EMC hit Creation") {}

// -----   Destructor   ----------------------------------------------------

MpdEmcHitCreation::~MpdEmcHitCreation() {}
// -------------------------------------------------------------------------

// -----   Public method Init   --------------------------------------------

InitStatus MpdEmcHitCreation::Init()
{

   cout << "******************* EMC INIT *********************" << endl;

   // Get RootManager
   FairRootManager *ioman = FairRootManager::Instance();
   if (!ioman) {
      cout << "-E- MpdEmcHitCreation::Init: "
           << "RootManager not instantiated!" << endl;
      return kFATAL;
   }

   // Get geometry EMC parameters
   fGeoPar = new MpdEmcGeoParams();

   // Get input array

   fPointArray = (TClonesArray *)ioman->GetObject("EmcPoint");
   if (!fPointArray) {
      cout << "-W- MpdEmcHitCreation::Init: "
           << "No EmcPoint array!" << endl;

      return kERROR;
   }

   fMcTrackArray = (TClonesArray *)ioman->GetObject("MCTrack");
   if (!fMcTrackArray) {
      cout << "-W- MpdEmcHitCreation::Init: "
           << "No MCTrack array!" << endl;
      return kERROR;
   }

   // Create and register output array

   fDigiArray = new TClonesArray("MpdEmcHit", 100);
   ioman->Register("MpdEmcHit", "EMC", fDigiArray, kTRUE);
   ioman->Register("MCTrack", "EMC", fMcTrackArray, kFALSE);
   cout << "-I- MpdEmcHitCreation: Intialization successfull" << endl;

   return kSUCCESS;
}

void MpdEmcHitCreation::Finish()
{

   cout << "\n-I- MpdEmcHitCreation: Finish" << endl;
}

void MpdEmcHitCreation::Exec(Option_t *opt)
{

   cout << "\n-I- MpdEmcHitCreation: Event No. " << FairRun::Instance()->GetEventHeader()->GetMCEntryNumber() << endl;

   // Reset output Array

   if (!fDigiArray) Fatal("MpdEmcHitCreation::Exec)", "No array of digits");
   fDigiArray->Delete();

   Double_t phiRow, rhoMod, zMod, thetaMod;
   for (UInt_t iPnt = 0; iPnt < fPointArray->GetEntriesFast(); ++iPnt) {

      FairMCPoint *pnt  = (FairMCPoint *)fPointArray->At(iPnt);
      Int_t        trId = pnt->GetTrackID();
      if (trId < 0) break;
      MpdMCTrack *tr = (MpdMCTrack *)fMcTrackArray->At(trId);
      //        if (tr->GetMotherId() != -1) continue;
      Int_t   pdg   = tr->GetPdgCode();
      Int_t   detID = pnt->GetDetectorID();
      Float_t x     = pnt->GetX();
      Float_t y     = pnt->GetY();
      Float_t z     = pnt->GetZ();

      // Get sector/row/tower number of one cell (tower)

      Int_t    sec         = detID / (1000 * fGeoPar->GetCrateNumber() * fGeoPar->GetNrows());
      Int_t    row         = detID / 1000;
      Int_t    tower       = detID - row * 1000;
      Double_t xTower      = fGeoPar->GetXTower(detID);
      Double_t yTower      = fGeoPar->GetYTower(detID);
      Double_t zTower      = fGeoPar->GetZTower(detID);
      Double_t phiTower    = fGeoPar->GetPhiTower(detID);
      Double_t thetaTower  = fGeoPar->GetThetaTower(detID);
      Float_t  e           = pnt->GetEnergyLoss();
      Double_t towerLength = fGeoPar->GetLengthBox();
      Double_t phi, theta = thetaTower, psi = -90.;
      phi = phiTower - 270.;
      if (phiTower < 90) phi = phiTower + 90.;

      Double_t master[3], local[3], trPos[3];
      master[0] = x;
      master[1] = y;
      master[2] = z;
      trPos[0]  = xTower;
      trPos[1]  = yTower;
      trPos[2]  = zTower;
      if (zTower < 0) {
         master[2] = -z;
         trPos[2]  = -zTower;
      }
      TGeoRotation    rotationMatrix;
      TGeoTranslation trans1(trPos[0], trPos[1], trPos[2]);
      rotationMatrix.SetAngles(phi, theta, psi);
      TGeoCombiTrans rotCombiTrans(trans1, rotationMatrix);
      rotCombiTrans.MasterToLocal(master, local);
      Double_t   cTimeDelay = (0.5 * towerLength - local[2]) * 1.6 / TMath::C() * 1.E+07;
      Double_t   timeEnergy = (cTimeDelay + pnt->GetTime()) * e;
      MpdEmcHit *hit        = SearchHit(detID);

      if (hit == NULL) {
         hit = new ((*fDigiArray)[fDigiArray->GetEntriesFast()]) MpdEmcHit(detID, e, timeEnergy);
      } else {
         hit->IncreaseEnergy(e);
         hit->IncreaseEnergyTime(timeEnergy);
      }

      hit->SetNumTracks(hit->GetNumTracks() + 1);
      hit->SetTrackId(trId);
      hit->SetPdg(pdg);
      hit->SetX(xTower);
      hit->SetY(yTower);
      hit->SetZ(zTower);
      hit->SetPhiCenter(phiTower);
      hit->SetThetaCenter(thetaTower);
   }

   for (int iHit = 0; iHit < fDigiArray->GetEntriesFast(); iHit++) {
      MpdEmcHit *hit = (MpdEmcHit *)fDigiArray->At(iHit);
      hit->SetTime(hit->GetTime() / hit->GetE());
      if (hit->GetNumTracks() > 1) {
         hit->SetPdg(0);
         hit->SetTrackId(-1);
      }
   }
}

MpdEmcHit *MpdEmcHitCreation::SearchHit(UInt_t detID)
{
   MpdEmcHit *foundHit = NULL;
   for (Int_t i = 0; i < fDigiArray->GetEntriesFast(); ++i) {
      MpdEmcHit *hit = (MpdEmcHit *)fDigiArray->At(i);
      if (hit->GetDetectorID() != detID) continue;
      foundHit = hit;
   }
   return foundHit;
}

// -------------------------------------------------------------------------

ClassImp(MpdEmcHitCreation);

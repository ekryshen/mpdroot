////////////////////////////////////////////////////////////////
//							      //
//  MpdEmcClusterCreation				      //
//  EMC clsuters in MpdEmcCluster, v03  		      //
//  Author List : Martemianov M., 2021		              //
//  							      //
////////////////////////////////////////////////////////////////

#include "TClonesArray.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRunSim.h"
#include "FairEventHeader.h"

#include "MpdEmcHit.h"
#include "MpdEmcCluster.h"
#include "MpdEmcClusterCreation.h"
#include "MpdMCTrack.h"
#include "TGeoManager.h"

//
// To produce clusters the in reco file (geometry version -> emc_v4) :
//
//  FairTask *emcHP = new MpdEmcHitCreation();
//  fRun->AddTask(emcHP);
//
//  MpdEmcClusterCreation *EmcCluster = new MpdEmcClusterCreation();
//  EmcCluster->SetClusterFrame(4, 3); -> Cluster frame
//  EmcCluster->SetEnergyThreshold(5.0); -> MeV
//  fRun->AddTask(EmcCluster);
//

namespace EMC {
UInt_t  nRows         = 300;  // number of rows (xy - plane)
UInt_t  nMods         = 128;  // number of modules in row
Float_t Distance3DCut = 20.0; // cm
Float_t HitThreshold  = 5.0;  // MeV
UInt_t  frameRow      = 5;    // row frame
UInt_t  frameMod      = 5;    // module frame
} // namespace EMC

using namespace std;
using namespace TMath;

// -----   Default constructor   -------------------------------------------

MpdEmcClusterCreation::MpdEmcClusterCreation()
   : FairTask("EMC cluster creation"), fEnergyThreshold(EMC::HitThreshold), fMaxClusterRadius(EMC::Distance3DCut),
     rowFrame(EMC::frameRow), modFrame(EMC::frameMod)
{
}

// -----   Destructor   ----------------------------------------------------

MpdEmcClusterCreation::~MpdEmcClusterCreation() {}
// -------------------------------------------------------------------------

// -----   Public method Init   --------------------------------------------

InitStatus MpdEmcClusterCreation::Init()
{

   cout << "******************* EMC INIT *********************" << endl;

   // Get RootManager

   FairRootManager *ioman = FairRootManager::Instance();
   if (!ioman) {
      cout << "-E- MpdEmcClusterCreation::Init: "
           << "RootManager not instantiated!" << endl;
      return kFATAL;
   }

   // Get input array

   fHitArray = (TClonesArray *)ioman->GetObject("MpdEmcHit");
   if (!fHitArray) {
      cout << "-W- MpdEmcHitCreation::Init: "
           << "No MpdEmcHit array!" << endl;

      return kERROR;
   }

   fMcTrackArray = (TClonesArray *)ioman->GetObject("MCTrack");
   if (!fMcTrackArray) {
      cout << "-W- MpdEmcClusterCreation::Init: "
           << "No MCTrack array!" << endl;
      return kERROR;
   }

   // Create and register output array

   fClusterArray = new TClonesArray("MpdEmcCluster", 1000);
   ioman->Register("MpdEmcCluster", "EMC", fClusterArray, kTRUE);
   ioman->Register("MCTrack", "EMC", fMcTrackArray, kTRUE);
   ////    ioman->Register("MpdEmcHit","EMC",fHitArray, kTRUE);

   cout << "-I- MpdEmcClusterCreation: Intialization successfull" << endl;

   return kSUCCESS;
}

void MpdEmcClusterCreation::Finish()
{

   cout << "\n-I- MpdEmcClusterCreation: Finish" << endl;
}

void MpdEmcClusterCreation::Exec(Option_t *opt)
{

   cout << "\n-I- MpdEmcClusterCreation: Event No. " << FairRun::Instance()->GetEventHeader()->GetMCEntryNumber()
        << endl;

   // Reset output Array

   if (!fClusterArray) Fatal("MpdEmcClusterCreation::Exec)", "No array of clusters");
   fClusterArray->Delete();

   UInt_t nHits = fHitArray->GetEntriesFast();

   MpdEmcCluster *cluster = NULL;

   vector<Int_t>   outHits, modId;
   Int_t           towerNumber, emcRow, emcMod, numHit;
   vector<Float_t> xHit, yHit, zHit, energy, time;
   Float_t         rho, phi, theta, rhoCenter, thetaCenter, phiCenter;
   Float_t         thetaPos, phiPos, rhoPos, zPos, fEnergy, fTime, valRad;

   if (nHits > 0) {

      for (UInt_t iHit = 0; iHit < nHits; ++iHit) {
         MpdEmcHit *hit = (MpdEmcHit *)fHitArray->At(iHit);
         if (hit->GetE() * 1000. > fEnergyThreshold) {
            emcRow      = hit->GetRow();
            emcMod      = hit->GetMod();
            towerNumber = 1000 * (emcRow + 1) + (emcMod + 1);
            modId.push_back(towerNumber);
            xHit.push_back(hit->GetX());
            yHit.push_back(hit->GetY());
            zHit.push_back(hit->GetZ());
            energy.push_back(hit->GetE());
            time.push_back(hit->GetTime());
         }
      }

      TVector3                fPos;
      vector<Int_t>::iterator iterMod;
      vector<Int_t>           iterTower  = modId;
      vector<Float_t>         iterEnergy = energy;
      vector<Float_t>         eH, xH, yH, zH, tH;

      while (iterEnergy.size()) {

         outHits.clear();
         numHit = 0;
         eH.clear();
         xH.clear();
         yH.clear();
         zH.clear();
         tH.clear();

         UInt_t indexMax = distance(iterEnergy.begin(), max_element(iterEnergy.begin(), iterEnergy.end()));
         emcRow          = iterTower[indexMax] / 1000;
         emcMod          = iterTower[indexMax] - emcRow * 1000;
         fEnergy         = 0.;
         fTime           = 0.;
         phiPos          = 0.;
         thetaPos        = 0.;
         SearchFrameHits(emcRow, emcMod, outHits);

         indexMax = distance(modId.begin(), find(modId.begin(), modId.end(), iterTower[indexMax]));
         rho =
            sqrt(xHit[indexMax] * xHit[indexMax] + yHit[indexMax] * yHit[indexMax] + zHit[indexMax] * zHit[indexMax]);
         phi = TMath::ATan2(yHit[indexMax], xHit[indexMax]) * RadToDeg();
         if (phi < 0) phi += 360;
         phiCenter   = phi;
         rhoCenter   = rho;
         thetaCenter = TMath::ACos(zHit[indexMax] / rho);

         for (Int_t nInd = 0; nInd < outHits.size(); nInd++) {
            iterMod = find(iterTower.begin(), iterTower.end(), outHits[nInd]);
            if (iterMod != iterTower.end()) {

               UInt_t ind = distance(modId.begin(), find(modId.begin(), modId.end(), *iterMod));
               UInt_t pos = iterMod - iterTower.begin();
               iterTower.erase(iterMod);
               iterEnergy.erase(iterEnergy.begin() + pos);

               eH.insert(eH.end(), energy[ind]);
               xH.insert(xH.end(), xHit[ind]);
               yH.insert(yH.end(), yHit[ind]);
               zH.insert(zH.end(), zHit[ind]);
               tH.insert(tH.end(), time[ind]);

               rho   = sqrt(xHit[ind] * xHit[ind] + yHit[ind] * yHit[ind] + zHit[ind] * zHit[ind]);
               theta = TMath::ACos(zHit[ind] / rho);
               phi   = TMath::ATan2(yHit[ind], xHit[ind]) * RadToDeg();
               fEnergy += energy[ind];
               thetaPos += theta * energy[ind];
               phiPos += MpdEmcMath::GetPhiDiff(phi, phiCenter) * energy[ind];
               /////        phiPos += phi*energy[ind];
               //////        fTime += time[ind]*energy[ind]; numHit++;
            }
         }

         // Cluster parameters

         thetaPos = thetaPos / fEnergy;
         phiPos   = phiCenter + phiPos / fEnergy;
         rhoPos   = rhoCenter * sin(thetaPos);
         zPos     = rhoCenter * cos(thetaPos);
         fPos.SetXYZ(phiPos, rhoPos, zPos);
         UInt_t indMax = distance(eH.begin(), max_element(eH.begin(), eH.end()));
         fTime         = tH[indMax];

         /// Calculate cluster radius

         valRad = 0.0;
         for (Int_t iHit = 0; iHit < numHit; iHit++) {
            phi = ATan2(yH[iHit], xH[iHit]) * RadToDeg();
            rho = sqrt(xH[iHit] * xH[iHit] + yH[iHit] * yH[iHit] + zH[iHit] * zH[iHit]);
            if (phi < 0) phi += 360;
            theta = ACos(zH[iHit] / rho);
            valRad += (pow(rhoCenter * (theta - thetaPos), 2) +
                       pow(rhoPos * MpdEmcMath::GetPhiDiff(phi, phiPos) * DegToRad(), 2)) *
                      eH[iHit];
         }

         cluster = new ((*fClusterArray)[fClusterArray->GetEntriesFast()]) MpdEmcCluster(fEnergy, fTime, fPos);
         cluster->SetNHits(numHit);
         cluster->SetRad(valRad / fEnergy);
      }
   }
}

// Search relative hits in the definite frame

void MpdEmcClusterCreation::SearchFrameHits(Int_t row, Int_t mod, vector<Int_t> &outMod)
{

   Int_t curRow;
   Int_t nRow = EMC::nRows, nMod = EMC::nMods;
   Int_t dRowMin = row - rowFrame, dRowMax = row + rowFrame;
   Int_t dModMin = mod - modFrame, dModMax = mod + modFrame;

   if (dRowMin <= 0) dRowMin = nRow - fabs(dRowMin);
   if (dRowMax >= nRow) dRowMax = dRowMax - nRow;
   if (dModMin <= 0) dModMin = 1;
   if (dModMax >= nMod) dModMax = nMod;

   curRow = dRowMin;
   for (Int_t iGeo1 = 0; iGeo1 < 2 * rowFrame + 1; iGeo1++) {
      if (curRow == nRow + 1) curRow = 1;
      for (Int_t iGeo2 = dModMin; iGeo2 < dModMax + 1; iGeo2++) outMod.insert(outMod.end(), 1000 * curRow + iGeo2);
      curRow++;
   }
}

// Search relative hits

void MpdEmcClusterCreation::SearchRelativeHits(Int_t row, Int_t mod, vector<Int_t> &outMod)
{

   Int_t         ind;
   Int_t         nRow = EMC::nRows;
   Int_t         nMod = EMC::nMods;
   vector<Int_t> geoMod;

   row            = row + 1;
   mod            = mod + 1;
   UInt_t rowPrev = row - 1, rowNext = row + 1;
   if (row == 1) rowPrev = nRow;
   if (row == nRow) rowNext = 1;
   geoMod = MpdEmcMath::make_vector<int>()
            << 1000 * rowPrev + (mod - 1) << 1000 * rowPrev + (mod) << 1000 * rowPrev + (mod + 1)
            << 1000 * row + (mod - 1) << 1000 * row + (mod + 1) << 1000 * rowNext + (mod - 1) << 1000 * rowNext + (mod)
            << 1000 * rowNext + (mod + 1);

   if (mod - 1 == 0) geoMod.erase(geoMod.begin() + 1, geoMod.begin() + 3);
   if (mod + 1 > nMod) geoMod.erase(geoMod.begin() + 5, geoMod.begin() + 7);

   outMod.clear();
   outMod.insert(outMod.begin(), 1000 * (row) + mod);
   for (UInt_t iHit = 0; iHit < fHitArray->GetEntriesFast(); iHit++) {
      MpdEmcHit *hit = (MpdEmcHit *)fHitArray->At(iHit);
      ind            = 1000 * (hit->GetRow() + 1) + (hit->GetMod() + 1);
      for (UInt_t iGeo = 0; iGeo < geoMod.size(); iGeo++)
         if (ind == geoMod[iGeo]) outMod.push_back(ind);
   }
}

// -------------------------------------------------------------------------

ClassImp(MpdEmcClusterCreation)

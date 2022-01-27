//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      Implementation of class MpdTpcDriftTask
//      see MpdTpcDriftTask.h for details
//
// Environment:
//      Software developed for the MPD Detector at NICA.
//
// Author List:
//      Alexandr Zinchenko LHEP, JINR, Dubna - adapted for MPD from PANDARoot
//
//-----------------------------------------------------------

// This Class' Header ------------------
#include "MpdTpcDriftTask.h"
#include "MpdTpcClusterizerTask.h"

// Collaborating Class Headers --------
#include "TpcGas.h"
#include "TpcPrimaryCluster.h"
#include "TpcDriftedElectron.h"
#include "FairRootManager.h"
#include <TClonesArray.h>
#include <TMath.h>
#include <TRandom.h>
#include <TSystem.h>

// C/C++ Headers ----------------------
#include <math.h>
#include <iostream>

// Class Member definitions -----------

//__________________________________________________________________________
MpdTpcDriftTask::MpdTpcDriftTask()
   : FairTask("TPC Drift"), fkNsec2(MpdTpcClusterizerTask::fgkNsec2), fPersistence(kFALSE), fAttach(kTRUE),
     fDiffuse(kTRUE), fDistort(kFALSE), fSec(0)
{
   fPrimBranchName        = "TpcPrimClus";
   std::string tpcGasFile = gSystem->Getenv("VMCWORKDIR");
   tpcGasFile += "/geometry/Ar-90_CH4-10.asc";
   fGas = new TpcGas(tpcGasFile, 130);
   // std::cout<<*_gas<<std::endl;
}
//__________________________________________________________________________

MpdTpcDriftTask::~MpdTpcDriftTask()
{
   delete fGas;
}
//__________________________________________________________________________

InitStatus MpdTpcDriftTask::Init()
{
   // Get ROOT Manager
   FairRootManager *ioman = FairRootManager::Instance();

   if (ioman == 0) {
      Error("MpdTpcDriftTask::Init", "RootManager not instantiated!");
      return kERROR;
   }

   // Get input collection
   // std::cout << fkNsec2 << std::endl;
   fPrimArray = new TClonesArray *[fkNsec2];
   for (Int_t i = 0; i < fkNsec2; ++i) {
      TString primBranchName = fPrimBranchName;
      primBranchName += i;
      fPrimArray[i] = (TClonesArray *)ioman->GetObject(primBranchName);
      if (fPrimArray[i] == 0x0) {
         Error("TpcDriftTask::Init", "PrimaryElectron-array not found!");
         return kERROR;
      }
   }

   // create and register output array
   // SetPersistence();
   fDriftedArray = new TClonesArray *[fkNsec2];
   for (Int_t i = 0; i < fkNsec2; ++i) {
      TString branchName = "TpcDriftEl";
      branchName += i;
      fDriftedArray[i] = new TClonesArray("TpcDriftedElectron");
      ioman->Register(branchName, "Tpc", fDriftedArray[i], fPersistence);
   }

   // TODO: Get from geom!!!!
   fzGem = 125;

   return kSUCCESS;
}
//__________________________________________________________________________

void MpdTpcDriftTask::Exec(Option_t *opt)
{

   // Sector definition - should be taken somewhere else
   static const Int_t    nSec   = fkNsec2 / 2;
   static const Double_t phiSec = TMath::TwoPi() / nSec;

   if (TString(opt) == "") return;
   sscanf(opt, "%2d", &fSec);
   // Reset output Array
   if (fDriftedArray[fSec] == 0x0) Fatal("MpdTpcDriftTask::Exec)", "No DriftedElectronArray");
   fDriftedArray[fSec]->Delete();

   // loop over incoming electrons
   Int_t nc = fPrimArray[fSec]->GetEntriesFast();
   for (Int_t ic = 0; ic < nc; ++ic) {
      TpcPrimaryCluster *pcl = (TpcPrimaryCluster *)fPrimArray[fSec]->UncheckedAt(ic);
      // create single electrons
      Int_t q = pcl->q();
      for (Int_t ie = 0; ie < q; ++ie) {
         Double_t z_length = pcl->z();
         Double_t driftl   = 0;
         Double_t dx = 0, dy = 0, dt = 0;
         Int_t    size = 0;
         if (z_length == 0) continue;

         // Calculate drift time
         if (z_length > 0)
            driftl = fzGem - z_length;
         else
            driftl = fzGem + z_length;

         if (driftl < 0) continue;
         // attachment
         if (fAttach) {
            if (exp(-driftl * fGas->k()) < gRandom->Uniform()) continue;
         }
         // diffusion
         if (fDiffuse) {
            Double_t sqrtDrift = sqrt(driftl);
            Double_t sigmat    = fGas->Dt() * sqrtDrift;
            Double_t sigmal    = fGas->Dl() * sqrtDrift;
            dt                 = (driftl + gRandom->Gaus(0, sigmal)) / fGas->VDrift();
            dx                 = gRandom->Gaus(0, sigmat);
            dy                 = gRandom->Gaus(0, sigmat);
         }
         // drift distortions
         if (fDistort) {
            // TODO: to be implemented
         }

         // Check if electron goes to another sector
         Double_t xend = pcl->x() + dx, yend = pcl->y() + dy;
         Double_t phi  = TMath::ATan2(yend, xend);
         Int_t    isec = (Int_t)(phi / phiSec);
         if (phi < 0) {
            --isec;
            isec = nSec + isec;
         }
         if (fSec > nSec - 1) isec += nSec; // z < 0

         size = fDriftedArray[isec]->GetEntriesFast();
         new ((*fDriftedArray[isec])[size]) TpcDriftedElectron(xend, yend, pcl->t() + dt, pcl);
      } // end loop over electrons
   }    // end loop over clusters

   std::cout << fDriftedArray[fSec]->GetEntriesFast() << " electrons arriving at readout No. " << fSec << std::endl;

   return;
}
//__________________________________________________________________________

void MpdTpcDriftTask::Clear(Option_t *opt)
{

   TString option = opt;
   if (option == "") return;
   if (option.Contains("fill") && !fPersistence) return;
   Int_t isec = TString(option(4, 2)).Atoi();
   // std::cout << " sec " << isec << " " << option.Data() << std::endl;
   TString name = "TpcDriftEl";
   name += isec;
   if (option.Contains("fill"))
      FairRootManager::Instance()->GetOutTree()->GetBranch(name)->Fill();
   else {
      fDriftedArray[isec]->Delete();
      fDriftedArray[isec]->Expand(1000); // shrink TClonesArray
   }
}

ClassImp(MpdTpcDriftTask)

//-----------------------------------------------------------
// File and Version Information:
// $Id$
//
// Description:
//      Implementation of class MpdTpcClusterizerTask
//      see MpdTpcClusterizerTask.h for details
//
// Environment:
//      Software developed for the MPD Detector at NICA.
//
// Author List:
//      Alexandr Zinchenko LHEP, JINR, Dubna - adapted for MPD from PANDARoot
//
//-----------------------------------------------------------

// This Class' Header ------------------
#include "MpdTpcClusterizerTask.h"

// Collaborating Class Headers --------
#include "FairRootManager.h"
#include "TpcPoint.h"
#include "TpcGas.h"
#include "TpcPrimaryCluster.h"
#include "LinearInterpolPolicy.h"

#include "TClonesArray.h"
#include "TFile.h"
#include "TMath.h"
#include "TRandom.h"
#include "TSystem.h"
#include "TTree.h"

// C/C++ Headers ----------------------
#include <iostream>
#include <math.h>
#include <vector>

using namespace std;

// Class Member definitions -----------
const Int_t MpdTpcClusterizerTask::fgkNsec2 = 24;

//__________________________________________________________________________
MpdTpcClusterizerTask::MpdTpcClusterizerTask() : FairTask("TPC Clusterizer"), fPersistence(kFALSE)
{
   fPointBranchName       = "TpcPoint";
   std::string tpcGasFile = gSystem->Getenv("VMCWORKDIR");
   tpcGasFile += "/geometry/Ar-90_CH4-10.asc";
   fGas = new TpcGas(tpcGasFile, 130);
   std::cout << *fGas << std::endl;
}
//__________________________________________________________________________

MpdTpcClusterizerTask::~MpdTpcClusterizerTask()
{
   delete fGas;
}
//__________________________________________________________________________

void MpdTpcClusterizerTask::FinishTask()
{
   for (Int_t i = 0; i < fgkNsec2; ++i) fPrimArray[i]->Delete();
}
//__________________________________________________________________________

InitStatus MpdTpcClusterizerTask::Init()
{
   // Get ROOT Manager
   FairRootManager *ioman = FairRootManager::Instance();

   if (ioman == 0) {
      Error("MpdTpcClusterizerTask::Init", "RootManager not instantiated!");
      return kERROR;
   }

   // Get input collection
   fPointArray = (TClonesArray *)ioman->GetObject(fPointBranchName);

   if (fPointArray == 0) {
      Error("TpcClusterizerTask::Init", "Point-array not found!");
      return kERROR;
   }

   // create and register output array
   // SetPersistence();
   fPrimArray    = new TClonesArray *[fgkNsec2];
   TString name0 = "TpcPrimClus";
   for (Int_t i = 0; i < fgkNsec2; ++i) {
      fPrimArray[i] = new TClonesArray("TpcPrimaryCluster");
      TString name  = "TpcPrimClus";
      name += i;
      ioman->Register(name, "Tpc", fPrimArray[i], fPersistence);
   }

   return kSUCCESS;
}
//__________________________________________________________________________

void MpdTpcClusterizerTask::Exec(Option_t *opt)
{

   const Int_t nSec = fgkNsec2 / 2; // number of TPC readout sectors
   // Reset output Array
   if (fPrimArray == 0) Fatal("MpdTpcClusterizerTask::Exec)", "No PrimClusterArray");
   for (Int_t i = 0; i < fgkNsec2; ++i) fPrimArray[i]->Delete();

   Int_t np = fPointArray->GetEntriesFast();
   if (np < 2) {
      Warning("MpdTpcClusterizerTask::Exec", "Not enough Hits in Tpc for Digitization (<2)");
      return;
   }

   TpcPoint *theLastPoint = (TpcPoint *)fPointArray->UncheckedAt(0), *point = theLastPoint;

   vector<pair<Int_t, Int_t>> sectors[fgkNsec2];

   // Sector definition - should be taken somewhere else
   const Double_t phiSec = TMath::TwoPi() / nSec, dDist = 1.2;

   // Find segments of tracks with the same ID inside one sector and one half
   // and not broken, i.e. without reentrance (short enough distance between consecutive points)
   Int_t curpointTrackID, ip0 = 0;
   curpointTrackID = point->GetTrackID();
   Double_t phi    = TMath::ATan2(point->GetY(), point->GetX());
   Int_t    isec0  = (Int_t)(phi / phiSec);
   if (phi < 0) {
      --isec0;
      isec0 = nSec + isec0;
   }
   if (point->GetZ() < 0) isec0 += nSec;
   TVector3 pos0, pos;
   point->Position(pos0);

   for (Int_t ip = 1; ip < np; ++ip) {
      point      = (TpcPoint *)fPointArray->UncheckedAt(ip);
      phi        = TMath::ATan2(point->GetY(), point->GetX());
      Int_t isec = (Int_t)(phi / phiSec);
      if (phi < 0) {
         --isec;
         isec = nSec + isec;
      }
      if (point->GetZ() < 0) isec += nSec;
      point->Position(pos);
      TVector3 dist = pos - pos0;
      if (isec != isec0 || point->GetTrackID() != curpointTrackID || dist.Mag() > dDist) {
         sectors[isec0].push_back(pair<Int_t, Int_t>(ip0, ip));
         ip0             = ip;
         isec0           = isec;
         curpointTrackID = point->GetTrackID();
      }
      pos0 = pos;
   }
   sectors[isec0].push_back(pair<Int_t, Int_t>(ip0, np));

   // Do clusterization sector-by-sector
   for (Int_t isec = 0; isec < fgkNsec2; ++isec) {
      // for (Int_t isec = 0; isec < 12; ++isec) {
      Int_t icluster = 0;
      // Loop over sectors
      Int_t nSeg = sectors[isec].size();

      for (Int_t iseg = 0; iseg < nSeg; ++iseg) {
         // Loop over track segments
         Int_t ibeg = sectors[isec][iseg].first, iend = sectors[isec][iseg].second;
         // Create "virtual" point - back-extrapolation from the first point
         theLastPoint  = (TpcPoint *)fPointArray->UncheckedAt(ibeg);
         TpcPoint virt = *theLastPoint;
         theLastPoint->Momentum(pos0);
         theLastPoint->Position(pos);
         Double_t ppp = pos0.Mag();
         if (ppp > 1.e-6) {
            ppp = theLastPoint->GetStep() / ppp;
            for (Int_t i = 0; i < 3; ++i) pos[i] -= pos0[i] * ppp;
         } else
            ++ibeg;
         virt.SetPosition(pos);
         theLastPoint = &virt;

         for (Int_t ip = ibeg; ip < iend; ++ip) {
            // Loop over points
            point = (TpcPoint *)fPointArray->UncheckedAt(ip);
            // point->Print();

            // check if hits belong to the same track
            curpointTrackID = point->GetTrackID();
            if (curpointTrackID != theLastPoint->GetTrackID()) {
               cout << " -F- TpcClusterizerTask::Exec - Different track IDs " << curpointTrackID << " "
                    << theLastPoint->GetTrackID() << endl;
               exit(0);
            }

            Double_t dE = point->GetEnergyLoss() * 1E9; // convert from GeV to eV

            // Step 0: calculate the overall ammount of charge, produced
            if (dE < 0) {
               Error("TpcClusterizerTask::Exec", "Note: particle:: negative Energy loss!");
               continue;
            }
            Int_t qTotal   = (Int_t)TMath::Abs(dE / fGas->W());
            Int_t qCluster = 0;
            Int_t nCluster = 0;

            // Step 1: Create Clusters
            Double_t pointTime = point->GetTime();
            assert(pointTime >= 0);
            // while still charge not used-up distribute charge into next cluster
            Int_t size = fPrimArray[isec]->GetEntriesFast();
            while (qTotal > 0) {
               // roll dice for next clustersize
               qCluster = fGas->GetRandomCS(gRandom->Uniform());
               if (qCluster > qTotal) qCluster = qTotal;
               qTotal -= qCluster;
               // create cluster
               new ((*fPrimArray[isec])[size++])
                  TpcPrimaryCluster(pointTime, qCluster, TVector3(0, 0, 0), curpointTrackID, ip);

               ++nCluster;
            } // finish loop for cluster creation

            // Step 2: Distribute Clusters along track segment
            LinearInterpolPolicy().Interpolate(theLastPoint, point, fPrimArray[isec], icluster, nCluster);
            icluster += nCluster;
            theLastPoint = point;
         } // for (Int ip = ibeg; ip < iend;
      }    // finish loop over segments
      std::cout << "TpcClusterizer:: " << fPrimArray[isec]->GetEntriesFast() << " clusters created" << std::endl;

      TString option = "";
      option += isec;
      ExecuteTasks(option);
      // cout << FairRootManager::Instance()->GetObject(name) << " " <<
      // FairRootManager::Instance()->GetObject(name)->ClassName() << endl; cout <<
      // FairRootManager::Instance()->GetOutTree()->GetName() << " " <<
      // FairRootManager::Instance()->GetOutTree()->GetBranch(name) << endl;
      CleanTasks();

      // Reset data in all subtasks for selected sectors
      ClearData(isec);

   } // finish loop over sectors

   return;
}
//__________________________________________________________________________

void MpdTpcClusterizerTask::ClearData(Int_t isec)
{
   // Save and clear data structures

   const Int_t nSec = fgkNsec2 / 2; // number of TPC readout sectors
   TString     option, name;
   TTask      *task = 0x0;

   if (isec > 1 && isec < nSec || isec > nSec + 1 && isec < fgkNsec2) {
      Int_t saveSec = isec - 1;
      TIter next(fTasks);
      option = "fill";
      option += saveSec;
      while ((task = (TTask *)next())) task->Clear(option);
      name = "TpcPrimClus";
      name += saveSec;
      if (fPersistence) FairRootManager::Instance()->GetOutTree()->GetBranch(name)->Fill();
      if (isec > 2 && isec < nSec - 1 || isec > nSec + 2) {
         fPrimArray[saveSec]->Delete();
         fPrimArray[saveSec]->Expand(1000); // shrink TClonesArray
      }

      option = "clea";
      option += saveSec;
      next.Reset();
      while ((task = (TTask *)next())) task->Clear(option);
   }

   if (isec == nSec - 1) {
      // Process remaining sectors at Z > 0 slightly differently
      if (fPersistence) {
         for (Int_t i = 0; i < nSec; i += (nSec - 1)) {
            name = "TpcPrimClus";
            name += i;
            FairRootManager::Instance()->GetOutTree()->GetBranch(name)->Fill();
         }
      }
      TIter next(fTasks);
      option = "fill0";
      next.Reset();
      while ((task = (TTask *)next())) task->Clear(option);
      option = "clea0";
      next.Reset();
      while ((task = (TTask *)next())) task->Clear(option);

      option = "fill";
      option += isec;
      next.Reset();
      while ((task = (TTask *)next())) task->Clear(option);
      option = "clea";
      option += isec;
      next.Reset();
      while ((task = (TTask *)next())) task->Clear(option);

      for (Int_t i = 0; i < 2; ++i) {
         fPrimArray[i]->Delete();
         fPrimArray[i]->Expand(1000); // shrink TClonesArray
      }
      for (Int_t i = nSec - 2; i < nSec; ++i) {
         fPrimArray[i]->Delete();
         fPrimArray[i]->Expand(1000); // shrink TClonesArray
      }
   }

   if (isec == fgkNsec2 - 1) {
      // Process remaining sectors at Z < 0
      if (fPersistence) {
         for (Int_t i = nSec; i < fgkNsec2; i += (nSec - 1)) {
            name = "TpcPrimClus";
            name += i;
            FairRootManager::Instance()->GetOutTree()->GetBranch(name)->Fill();
         }
      }
      TIter next(fTasks);
      option = "fill";
      option += nSec;
      next.Reset();
      while ((task = (TTask *)next())) task->Clear(option);
      option = "clea";
      option += nSec;
      next.Reset();
      while ((task = (TTask *)next())) task->Clear(option);

      option = "fill";
      option += isec;
      next.Reset();
      while ((task = (TTask *)next())) task->Clear(option);
      option = "clea";
      option += isec;
      next.Reset();
      while ((task = (TTask *)next())) task->Clear(option);

      for (Int_t i = nSec; i < nSec + 2; ++i) {
         fPrimArray[i]->Delete();
         fPrimArray[i]->Expand(1000); // shrink TClonesArray
      }
      for (Int_t i = isec - 1; i <= isec; ++i) {
         fPrimArray[i]->Delete();
         fPrimArray[i]->Expand(1000); // shrink TClonesArray
      }
   }
}

ClassImp(MpdTpcClusterizerTask);

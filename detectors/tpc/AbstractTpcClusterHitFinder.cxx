//-----------------------------------------------------------
// Description:
//      AbstractTpcClusterHitFinder interface source file
//
// Environment:
//      Software for the MPD Detector at NICA.
//
// Author:
//      Slavomir Hnatic LIT, JINR, Dubna - 5.2022
//
//-----------------------------------------------------------

#include "AbstractTpcClusterHitFinder.h"
#include "AbstractTpc2dCluster.h"
#include "AbstractTpcHit.h"

#include "QA_TpcClusterHitFinder.h"
#include "FairLogger.h"

// ROOT Headers --------
#include "TFile.h"
#include "TH2F.h"
#include "TVector3.h"

#include <iostream>

//__________________________________________________________________________

AbstractTpcClusterHitFinder::AbstractTpcClusterHitFinder(BaseTpcSectorGeo &tpcGeo, const char *name, Bool_t val)
   : FairTask(name), persistence(val)
{
   secGeo = &tpcGeo;
}

//__________________________________________________________________________

AbstractTpcClusterHitFinder::~AbstractTpcClusterHitFinder() {}

//__________________________________________________________________________

InitStatus AbstractTpcClusterHitFinder::Init()
{
   return ReadInputRegisterOutput();
}

//__________________________________________________________________________

void AbstractTpcClusterHitFinder::Exec(Option_t *opt)
{
   LOG(info) << "!!! Size of digiArray is !!! " << digiArray->GetEntriesFast();
   // digiArray->Delete();
   ClearClustersHits();
   TransformInputData();
   FindClusters();
   FindHits();
}

//__________________________________________________________________________

void AbstractTpcClusterHitFinder::Finish()
{

   if (AbstractQA::qaEngineMode == EQAMode::TPCCLUSTERHITFINDER) {
      TString infoMessage = TString("QA Engine - TpcClusterHitFinder mode: ") + ModuleNameSuffix() + TString(" module");
      LOG(info) << infoMessage;
      TString qaFileName = TString("QA_TpcClusterHitFinder") + ModuleNameSuffix() + ".root";
      TFile   f(qaFileName, "recreate");

      LOG(info) << "!!! Size of digiArray is !!! " << digiArray->GetEntriesFast();
      LOG(info) << "!!! Size of clusArray is !!! " << clusArray->GetEntriesFast();
      LOG(info) << "!!! Size of hitArray is !!! " << hitArray->GetEntriesFast();

      int              clusterNumber = clusArray->GetEntriesFast();
      TObjArray        plotList(clusterNumber);
      int              hitNumber = hitArray->GetEntriesFast();
      TObjArray        hitList(hitNumber);
      std::vector<int> hitStartIndex;

      int startMlemHitIndex = 0;
      for (int i = 0; i < clusterNumber; ++i) {
         // adding cluster
         AbstractTpc2dCluster *cluster = (AbstractTpc2dCluster *)clusArray->At(i);
         int                   padMin  = cluster->GetPadMin() - 1;
         int                   padMax  = cluster->GetPadMax() + 2;
         if (ModuleNameSuffix() == TString("Mlem")) {
            padMin -= 2;
            padMax -= 2;
         }
         int nPad       = padMax - padMin;
         int timeBinMin = cluster->GetTimeBinMin() - 1;
         int timeBinMax = cluster->GetTimeBinMax() + 2;
         int nTimeBin   = timeBinMax - timeBinMin;
         int sector     = cluster->GetSector();
         int row        = cluster->GetRow();

         LOG(info) << "Cluster " << i << ": " << sector << " " << row << " " << padMin << " " << padMax << " "
                   << timeBinMin << " " << timeBinMax;

         TString plotName =
            "Cluster " + std::to_string(i) + ", sector " + std::to_string(sector) + ", padrow " + std::to_string(row);
         TString plotTitle = "Cluster Map: " + ModuleNameSuffix();

         TH2F *clusterPlot = new TH2F(plotName, plotTitle, nPad, padMin, padMax, nTimeBin, timeBinMin, timeBinMax);
         plotList.Add(clusterPlot);
         clusterPlot->SetMaximum(100);

         // fill cluster with digits
         std::vector<AbstractTpcDigit *> clusterDigits = cluster->GetClusterDigits();
         for (int j = 0; j < clusterDigits.size(); ++j) {
            double digitPad = clusterDigits[j]->GetPad();
            if (ModuleNameSuffix() == TString("Mlem")) digitPad -= 2;
            double digitTimeBin = clusterDigits[j]->GetTimeBin();
            double digitSignal  = clusterDigits[j]->GetSignal();

            clusterPlot->Fill(digitPad, digitTimeBin, digitSignal);
            LOG(info) << digitPad << " " << digitTimeBin << " " << digitSignal;
         }

         // adding hit(s)
         if (ModuleNameSuffix() != TString("Mlem")) {
            LOG(info) << "Adding hits to cluster";
            AbstractTpcHit *hit       = (AbstractTpcHit *)hitArray->At(i);
            float           hitPad    = hit->GetPadCoordinate();
            float           hitTime   = hit->GetTimeBinCoordinate();
            float           hitSignal = hit->GetTotalSignal();
            TVector3       *currHit   = new TVector3(hitPad, hitTime, hitSignal);
            hitList.Add(currHit);
            hitStartIndex.push_back(i);

            LOG(info) << "Hit " << i << ": " << hitPad << " " << hitTime << " " << hitSignal;
            // LOG(info) << "List size " << i << ": " << ((TList *)hitMap.GetValue(i))->GetSize();
         } else {
            LOG(info) << "Adding MLEM hits to cluster";
            int kStart = startMlemHitIndex; // start Hit index for current cluster
            hitStartIndex.push_back(startMlemHitIndex);
            for (int k = kStart; k < hitArray->GetEntriesFast(); ++k) {
               AbstractTpcHit *hit       = (AbstractTpcHit *)hitArray->At(k);
               int             clusterID = hit->GetClusterID();
               if (clusterID > i) {
                  break;
               } else {
                  startMlemHitIndex++;
                  float     hitPad    = hit->GetPadCoordinate() / secGeo->PAD_WIDTH[1] + secGeo->PAD_COUNT[row];
                  float     hitTime   = hit->GetTimeBinCoordinate() / 0.55 + secGeo->TIMEBIN_COUNT;
                  float     hitSignal = hit->GetTotalSignal();
                  TVector3 *currHit   = new TVector3(hitPad, hitTime, hitSignal);
                  hitList.Add(currHit);

                  LOG(info) << "Hit " << k << " from cluster " << clusterID << ": " << hitPad << " " << hitTime << " "
                            << hitSignal;
               }
            }
         }
      }

      plotList.Write();
      // NOTE: for some unknown stupid reason, the hits are written to file starting from last hit
      // after reading it is important to reverse this back
      hitList.Write();
      f.WriteObject(&hitStartIndex, "hitStartIndex");
      f.Close();
      LOG(info) << TString("QA Engine - Plots written to file: ") + qaFileName;
   }

   // QA_ClusterHitFinder qaPlots;
   // qaPlots.WriteOut();
   //  try to output arrays with data
}

//__________________________________________________________________________

InitStatus AbstractTpcClusterHitFinder::ReadInputRegisterOutput()
{
   // get FairRoot Manager
   ioman = FairRootManager::Instance();
   if (!ioman) {
      Error("AbstractTpcClusterHitFinder::readInputRegisterOutput", "FairRootManager not instantiated!");
      return kERROR;
   }

   // Read input digi collection
   digiArray = (TClonesArray *)ioman->GetObject("MpdTpcDigit");
   if (!digiArray) {
      Error("AbstractTpcClusterHitFinder::readInputRegisterOutput", "Array of digits not found!");
      return kERROR;
   }

   LOG(info) << "!!! INITIAL Size of digiArray is !!! " << digiArray->GetEntriesFast();
   // Get event header (for vertex Z-position)
   /////////////////////////////////////////////////////////////////////////
   // This is needed ONLY for AZ's MLEM algorithm,                        //
   // placed here to be able to get rid of Init() in his implementation   //
   // should be moved away from interface to his implementation of it     //
   ioman->GetObject("MCEventHeader");
   /////////////////////////////////////////////////////////////////////////

   // Create and register output arrays
   clusArray = new TClonesArray("MpdTpc2dCluster");
   ioman->Register("TpcCluster", "Tpc", clusArray, persistence);
   hitArray = new TClonesArray("MpdTpcHit");
   ioman->Register("TpcRecPoint", "Tpc", hitArray, persistence);

   return kSUCCESS;
}

//__________________________________________________________________________

void AbstractTpcClusterHitFinder::ClearClustersHits()
{
   clusArray->Delete();
   hitArray->Delete();
}

//__________________________________________________________________________

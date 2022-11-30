#include "TpcClusterHitFinderFast.h"

// Collaborating Class Headers --------
#include "MpdTpcDigit.h"
#include "MpdTpcHit.h"
#include "TpcClustering.h"
#include "FairRootManager.h"

#include <TClonesArray.h>
#include <TMath.h>
#include <TROOT.h>
#include <TSystem.h>

// C/C++ Headers ----------------------
#include <iostream>
#include <math.h>
#include <vector>

#define nTPC_ROWS 53
#define nTPC_SECTORS 24
#define fV_DRIFT 0.0055 // cm/ns

using namespace std;
using namespace tpcClustering;
using namespace chrono;

TpcClusterHitFinderFast::TpcClusterHitFinderFast(BaseTpcSectorGeo &secGeo)
   : FairTask("TPC Cluster finder"), fPersistence(kFALSE), _nEvent(-1), _nSector(-1), _nRow(-1)
{
   _pSecGeo = dynamic_cast<TpcSectorGeoAZ *>(&secGeo);
   if (!_pSecGeo) Fatal("TpcClusterHitFinderFast::TpcClusterHitFinderFast", " !!! Wrong geometry type !!! ");
}

//__________________________________________________________________________

TpcClusterHitFinderFast::~TpcClusterHitFinderFast() {}

//__________________________________________________________________________

void TpcClusterHitFinderFast::FinishTask()
{
   // cout << "MLEM cluster finder work time = " << ((Float_t)tAllMlem) / CLOCKS_PER_SEC << endl;
   // ExecTime->Write();
   //  cout << "Average Execution Time: " << _nTime/1000 << " milliseconds" << endl;
   //  float fError = 0;
   //  for (int i = 0; i < 1000; i++) {
   //  fError += (fTime[i] - _nTime/1000) * (fTime[i] - _nTime/1000);
   // }
   // cout << "Execution Time Error: " << sqrt(fError/999) << endl;
}

//__________________________________________________________________________

InitStatus TpcClusterHitFinderFast::Init()
{
   std::cout << "TpcClusterHitFinderFast::Init started" << std::endl;
   // ExecTime = new TH1F("ExecTime", "Execution Time", 50, 0, 50);
   //_nTime = 0;

   // Get ROOT Manager
   FairRootManager *ioman = FairRootManager::Instance();
   if (ioman == 0) {
      Error("TpcClusterHitFinderFast::Init", "RootManager not instantiated!");
      return kERROR;
   }

   // Get input collection
   _pTpcBinDataAZlt = (TClonesArray *)ioman->GetObject("MpdTpcDigit");

   if (_pTpcBinDataAZlt == 0) {
      Error("TpcClusterHitFinderFast::Init", "Array of digits not found!");
      return kERROR;
   }

   _pTpcHitFinder = new TClonesArray("MpdTpcHit");
   ioman->Register("TpcRecPoint", "Tpc", _pTpcHitFinder, fPersistence);

   std::cout << "TpcClusterHitFinderFast::Init finished" << std::endl;
   return kSUCCESS;
}

//__________________________________________________________________________

void TpcClusterHitFinderFast::Exec(Option_t *opt)
{
   _pTpcHitFinder->Delete();

   std::ostringstream oss;

   FairRootManager *ioman = FairRootManager::Instance();
   _nEvent                = ioman->GetEntryNr();

   std::cout << "TpcClusterHitFinderFast::Exec started: Event " << _nEvent << std::endl;

   auto start = high_resolution_clock::now();

   EventClusters *pEventClusters = new EventClusters(_nEvent);

   int nS = _nSector, nR = _nRow;
   int nSector_ = nS;
   int nRow_    = nR;
   if (nS >= 0) pEventClusters->Add(new SectorClusters(nS));

   std::string strRows = std::string(nTPC_ROWS, '-');

   std::list<PadCluster *>                        lstPadClusters;
   std::vector<std::list<PadCluster *>::iterator> viPadClusters;

   Int_t nPoints = _pTpcBinDataAZlt->GetEntriesFast();
   for (Int_t i = 0; i < nPoints; i++) {
      MpdTpcDigit *pdigit   = (MpdTpcDigit *)_pTpcBinDataAZlt->UncheckedAt(i);
      uint         nSector  = pdigit->GetSector();
      uint         nRow     = pdigit->GetRow();
      uint         nPad     = pdigit->GetPad();
      uint         nTimeBin = pdigit->GetTimeBin();
      float        fAdc     = pdigit->GetAdc();
      int          nTrackId = pdigit->GetOrigin();
      AdcHit       adcHit(nTimeBin, fAdc, nTrackId);

      if (nSector_ < 0) { // 1st occur only!
         nSector_ = nSector;
         pEventClusters->Add(new SectorClusters(nSector));
      }
      if (nRow_ < 0 && nSector == nSector_) // 1st occur only!
         nRow_ = nRow;

      bool bSector = nS < 0 && nSector > nSector_ || nS >= 0 && nSector > nS;
      bool bRow    = nR < 0 && nRow > nRow_ || nR >= 0 && nRow > nR;

      if (bSector || bRow && nSector == nSector_) { // the row was fulfilled and ready for clustering
         if (viPadClusters.size()) {
            viPadClusters.push_back(lstPadClusters.end()); // add the end of the PadClusters
            RowClusters *pRowClusters = new RowClusters(nR < 0 ? nRow_ : nR);
            pRowClusters->Joint(viPadClusters); // joint different PadClusters to RowCluster
            pRowClusters->Split();              // split the RowClusters to HitClusters

            // Drawing clusters in oss (change to true if necessary)
            if (false) {
               DrawPadClusters(oss, viPadClusters);
               pRowClusters->Draw(oss);
               pRowClusters->Print(oss);
            }

            assert(pEventClusters->getBack() != NULL);
            pEventClusters->getBack()->Add(pRowClusters);
            if (bSector && nS < 0) pEventClusters->Add(new SectorClusters(nSector));
         }
         viPadClusters.clear();
         lstPadClusters.clear();
         if (bSector) {
            strRows[nRow_] = szDECADE[nRow_ % 10];
            if (nS < 0) {
               strRows  = std::string(nTPC_ROWS, '-');
               nRow_    = nRow;
               nSector_ = nSector;
            } else
               break;
         } else if (nR < 0) {
            strRows[nRow_] = szDECADE[nRow_ % 10];
            nRow_          = nRow;
         } else
            break;
      }

      bSector = nS < 0 && nSector == nSector_ || nS >= 0 && nSector == nS;
      bRow    = nR < 0 && nRow == nRow_ || nR >= 0 && nRow == nR;

      if (bSector && bRow) {           // add the fired pad to a PadCluster
         if (lstPadClusters.empty()) { // no one cluster is yet
            lstPadClusters.push_back(new PadCluster(nPad, adcHit));
            viPadClusters.push_back(lstPadClusters.begin());
         } else {
            PadCluster *pPadCluster = lstPadClusters.back();
            if (nPad != pPadCluster->getPad()) { // new PadCluster
               lstPadClusters.push_back(new PadCluster(nPad, adcHit));
               viPadClusters.push_back(std::prev(lstPadClusters.end()));
            } else if (pPadCluster->isOut(adcHit)) { // new cluster in the pad
               lstPadClusters.push_back(new PadCluster(nPad, adcHit));
            } else if (pPadCluster->isCombined(adcHit)) { // combined PadCluster
               PadCluster *pCombinedPadCluster = new PadCluster(nPad, pPadCluster->splitAdcHit(adcHit));
               pCombinedPadCluster->setPadClusterBefore(pPadCluster);
               pPadCluster->setPadClusterAfter(pCombinedPadCluster);
               pCombinedPadCluster->Add(adcHit);
               lstPadClusters.push_back(pCombinedPadCluster);
            } else // simple PadCluster
               pPadCluster->Add(adcHit);
         }
      }
   }
   if (viPadClusters.size()) {
      viPadClusters.push_back(lstPadClusters.end()); // add the end of the PadClusters
      RowClusters *pRowClusters = new RowClusters(nR < 0 ? nRow_ : nR);
      pRowClusters->Joint(viPadClusters);
      pRowClusters->Split();
      if (nR < 0) strRows[nRow_] = szDECADE[nRow_ % 10];
      // Drawing clusters in oss (change to true if necessary)
      if (false) {
         DrawPadClusters(oss, viPadClusters);
         pRowClusters->Draw(oss);
         pRowClusters->Print(oss);
      }
      assert(pEventClusters->getBack() != NULL);
      pEventClusters->getBack()->Add(pRowClusters);
      TpcClusterHitFinderFast::calcSector(pEventClusters);
   }
   // cout << oss.str() << endl;
   auto stop     = high_resolution_clock::now();
   auto duration = duration_cast<milliseconds>(stop - start);

   // Counting Execution Time
   // ExecTime->Fill(duration.count());
   // fTime[_nEvent] = duration.count();
   //_nTime += duration.count();

   viPadClusters.clear();
   lstPadClusters.clear();

   std::cout << "TpcClusterHitFinderFast::Exec finished for " << duration.count() << " milliseconds" << std::endl;
}

void TpcClusterHitFinderFast::calcSector(const EventClusters *pEventClusters)
{
   Int_t nClus = 0;
   Int_t nHit  = 0;
   // float r3 = 0.1 + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(0.5-0.1)));
   TVector3                           tv3GlobErr(0.05, 0.0, 0.2);
   TVector3                           p3glob(0.0, 0.0, 0.0);
   const std::list<SectorClusters *> &rlstEventClusters = pEventClusters->getSectorClusters();
   for (std::list<SectorClusters *>::const_iterator i0 = rlstEventClusters.begin(); i0 != rlstEventClusters.end();
        i0++) {
      const std::list<RowClusters *> &rlstRowClusters = (*i0)->getRowClusters();
      for (std::list<RowClusters *>::const_iterator i1 = rlstRowClusters.begin(); i1 != rlstRowClusters.end(); i1++) {
         const std::list<Cluster *> &rlstClusters = (*i1)->getClusters();
         for (std::list<Cluster *>::const_iterator i2 = rlstClusters.begin(); i2 != rlstClusters.end(); i2++) {
            const Cluster *pCluster = *i2;

            const std::list<const PadCluster *> &rlstPadClusters = (*i2)->getPadClusters();

            int   nRow  = (*i1)->getRow();
            int   nSect = (*i0)->getSector();
            float fPad  = pCluster->getPad();
            float fTime = pCluster->getTime();

            // uint nPadNum = pCluster->getPadMax() - pCluster->getPadMin() + 1;
            // uint nTimeBinNum = pCluster->getTimeBinMax() - pCluster->getTimeBinMin() + 1;
            // std::cout << "Pad: " << pCluster->getPad() << std::endl;
            // std::cout << "Row: " << (*i1)->getRow() << std::endl;
            // std::cout << "Sect: " << (*i0)->getSector() << std::endl;
            // std::cout << "TimeBinReco: " << pCluster->getTime() << std::endl;

            // --------------------------------TpcSectorGeoAZ--------------------------------------------------
            Double_t xloc  = _pSecGeo->Pad2Xloc(fPad, nRow);
            Int_t    padID = _pSecGeo->PadID(nSect % _pSecGeo->NofSectors(), nRow);
            Double_t yloc  = _pSecGeo->LocalPadPosition(padID).Y();
            Double_t zloc  = _pSecGeo->TimeBin2Z(fTime);
            TVector3 p3loc(xloc, yloc, zloc);
            // std::cout << "LocalXYZ: " << p3loc.X() << ", " << p3loc.Y() << ", " << p3loc.Z() << std::endl;
            if (nSect >= _pSecGeo->NofSectors()) p3loc[2] = -p3loc[2];
            _pSecGeo->Local2Global(nSect, p3loc, p3glob);
            // std::cout << "GlobalXYZ: " << p3glob.X() << ", " << p3glob.Y() << ", " << p3glob.Z() << std::endl;

            nHit           = _pTpcHitFinder->GetEntriesFast();
            MpdTpcHit *hit = new ((*_pTpcHitFinder)[nHit++]) MpdTpcHit(padID, p3glob, tv3GlobErr, nClus++);
            hit->SetLayer(nRow);
            hit->SetLocalPosition(p3loc); // point position
            hit->SetEnergyLoss(pCluster->getAdcSum());

            Int_t    ireg = (nRow < _pSecGeo->NofRowsReg(0)) ? 0 : 1;
            Double_t step = _pSecGeo->PadHeight(ireg);
            hit->SetStep(step);
            hit->SetModular(1); // modular geometry flag

            hit->SetPad(int(fPad));
            hit->SetBin(int(fTime));
            hit->SetNdigits(rlstPadClusters.size());
            // hit->SetRMS(0.05, 0);
            // hit->SetRMS(0.2, 1);
            //  std::cout << "Row: " << (*i1)->getRow() << " NDig: " << rlstPadClusters.size() << std::endl;
         }
      }
   }
}
ClassImp(TpcClusterHitFinderFast);

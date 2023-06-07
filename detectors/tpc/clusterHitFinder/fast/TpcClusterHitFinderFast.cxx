#include "TpcClusterHitFinderFast.h"
#include "Tpc2dClusterFast.h"
#include "MpdTpcDigit.h"
#include "MpdTpcHit.h"

// ROOT Class Headers --------
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
   : AbstractTpcClusterHitFinder(secGeo, "TPC Cluster finder Fast", kFALSE)
{
   _pSecGeo = dynamic_cast<TpcSectorGeoAZ *>(&secGeo);
   if (!_pSecGeo) Fatal("TpcClusterHitFinderFast::TpcClusterHitFinderFast", " !!! Wrong geometry type !!! ");
}

//__________________________________________________________________________

TpcClusterHitFinderFast::~TpcClusterHitFinderFast() {}

//__________________________________________________________________________

void TpcClusterHitFinderFast::FindHits()
{
   ostringstream oss;

   _nEvent = ioman->GetEntryNr();

   cout << "TpcClusterHitFinderFast::Exec started: Event " << _nEvent << endl;

   auto start = high_resolution_clock::now();

   EventClusters *pEventClusters = new EventClusters(_nEvent);

   int nS = _nSector, nR = _nRow;
   int nSector_ = nS;
   int nRow_    = nR;
   if (nS >= 0) pEventClusters->Add(new SectorClusters(nS));

   string strRows = string(nTPC_ROWS, '-');

   list<PadCluster *>                   lstPadClusters;
   vector<list<PadCluster *>::iterator> viPadClusters;

   Int_t nPoints = digiArray->GetEntriesFast();
   for (int i = 0; i < nPoints; i++) {
      MpdTpcDigit *pdigit   = (MpdTpcDigit *)digiArray->UncheckedAt(i); // Input digit
      uint         nSector  = pdigit->GetSector();
      uint         nRow     = pdigit->GetRow();
      uint         nPad     = pdigit->GetPad();
      uint         nTimeBin = pdigit->GetTimeBin();
      float        fAdc     = pdigit->GetAdc();
      // int          nTrackId = pdigit->GetOrigin();
      // AdcHit       adcHit(nTimeBin, fAdc, nTrackId);
      AdcHit adcHit(nTimeBin, fAdc);

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
            // if (false) {
            //   DrawPadClusters(oss, viPadClusters);
            //   pRowClusters->Draw(oss);
            //   pRowClusters->Print(oss);
            //}

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
            } else if (pPadCluster->isOut(adcHit)) { // new Cluster in the pad
               lstPadClusters.push_back(new PadCluster(nPad, adcHit));
            } else if (pPadCluster->isCombined(adcHit)) { // combined PadCluster
               PadCluster *pCombinedPadCluster = new PadCluster(nPad, pPadCluster->splitAdcHit(adcHit));
               pCombinedPadCluster->setPadClusterBefore(pPadCluster);
               pPadCluster->setPadClusterAfter(pCombinedPadCluster);
               pCombinedPadCluster->Add(adcHit);
               lstPadClusters.push_back(pCombinedPadCluster);
            } else
               pPadCluster->Add(adcHit); // simple PadCluster
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
      // if (false) {
      //   DrawPadClusters(oss, viPadClusters);
      //   pRowClusters->Draw(oss);
      //   pRowClusters->Print(oss);
      //}
      assert(pEventClusters->getBack() != NULL);
      pEventClusters->getBack()->Add(pRowClusters);
      TpcClusterHitFinderFast::calcSector(pEventClusters);
   }
   // Drawing clusters in oss
   // cout << oss.str() << endl;
   auto stop     = high_resolution_clock::now();
   auto duration = duration_cast<milliseconds>(stop - start);

   // Time managment
   /* Counting Execution Time
   ExecTime->Fill(duration.count());
   fTime[_nEvent] = duration.count();
   _nTime += duration.count();
   */

   viPadClusters.clear();
   lstPadClusters.clear();

   // std::cout << "TpcClusterHitFinderFast::Exec finished for " << duration.count() << " milliseconds" << std::endl;
}

void TpcClusterHitFinderFast::calcSector(const EventClusters *pEventClusters)
{
   int                           nClus = 0;
   int                           nHit  = 0;
   TVector3                      tv3GlobErr(0.05, 0.0, 0.2);
   TVector3                      p3glob(0.0, 0.0, 0.0);
   const list<SectorClusters *> &rlstEventClusters = pEventClusters->getSectorClusters();
   for (list<SectorClusters *>::const_iterator i0 = rlstEventClusters.begin(); i0 != rlstEventClusters.end(); i0++) {
      // Sectors
      const list<RowClusters *> &rlstRowClusters = (*i0)->getRowClusters();
      for (list<RowClusters *>::const_iterator i1 = rlstRowClusters.begin(); i1 != rlstRowClusters.end();
           i1++) { // Rows
         const list<Cluster *> &rlstClusters = (*i1)->getClusters();
         for (list<Cluster *>::const_iterator i2 = rlstClusters.begin(); i2 != rlstClusters.end(); i2++) { // Clusters

            const Cluster *pCluster = *i2;

            int   nRow  = (*i1)->getRow();
            int   nSect = (*i0)->getSector();
            float fPad  = pCluster->getPad();
            float fTime = pCluster->getTime();

            // --------------------------------TpcSectorGeoAZ--------------------------------------------------
            double   xloc  = _pSecGeo->Pad2Xloc(fPad, nRow);
            int      padID = _pSecGeo->PadID(nSect % _pSecGeo->NofSectors(), nRow);
            double   yloc  = _pSecGeo->LocalPadPosition(padID).Y();
            double   zloc  = _pSecGeo->TimeBin2Z(fTime);
            TVector3 p3loc(xloc, yloc, zloc);
            if (nSect >= _pSecGeo->NofSectors()) p3loc[2] = -p3loc[2];
            _pSecGeo->Local2Global(nSect, p3loc, p3glob);

            /* write cluster(s) */
            nHit                   = hitArray->GetEntriesFast();
            Tpc2dClusterFast *clus = new ((*clusArray)[nHit]) Tpc2dClusterFast(); // Add new cluster
            clus->SetClusterID(nHit);
            clus->SetSector(nSect);
            clus->SetRow(nRow);
            clus->SetTimeBinMin(pCluster->getTimeBinMin());
            clus->SetTimeBinMax(pCluster->getTimeBinMax());
            clus->SetPadMin(pCluster->getPadMin());
            clus->SetPadMax(pCluster->getPadMax());
            // write digits into cluster
            // vector<pair<int, float>>        vnTrackId;
            const list<const PadCluster *> &rlstPadClusters = (*i2)->getPadClusters();
            for (list<const PadCluster *>::const_iterator i3 = rlstPadClusters.begin(); i3 != rlstPadClusters.end();
                 i3++) {
               const vector<AdcHit> &rvAdcHits = (*i3)->getAdcHits();
               for (vector<AdcHit>::const_iterator i4 = rvAdcHits.begin(); i4 != rvAdcHits.end(); i4++) { // Digits
                  // vnTrackId = i4->getTrackID();
                  AbstractTpcDigit *pDigit =
                     new MpdTpcDigit(0, (*i3)->getPad(), nRow, i4->getTimeBin(), nSect, i4->getAdc());
                  // new MpdTpcDigit(vnTrackId[0].first, (*i3)->getPad(), nRow, i4->getTimeBin(), nSect, i4->getAdc());
                  clus->AddDigit(pDigit);
               }
            }

            /* write hit(s) */
            MpdTpcHit *hit = new ((*hitArray)[nHit++]) MpdTpcHit(padID, p3glob, tv3GlobErr, nClus++); // Add new hit
            hit->SetLayer(nRow);
            hit->SetLocalPosition(p3loc); // point position
            hit->SetEnergyLoss(pCluster->getAdcSum());

            int    ireg = (nRow < _pSecGeo->NofRowsReg(0)) ? 0 : 1;
            double step = _pSecGeo->PadHeight(ireg);
            hit->SetStep(step);
            hit->SetModular(1); // modular geometry flag

            hit->SetPad(int(fPad));
            hit->SetBin(int(fTime));
            hit->SetNdigits(rlstPadClusters.size());

            // for (vector<pair<int, float>>::const_iterator i5 = vnTrackId.begin(); i5 != vnTrackId.end(); i5++) {
            //    hit->AddTrackID((*i5).first, (*i5).second); // TrackID
            // }
         }
      }
   }
}
ClassImp(TpcClusterHitFinderFast);

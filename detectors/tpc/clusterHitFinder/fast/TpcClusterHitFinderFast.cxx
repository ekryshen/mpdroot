#include "TpcClusterHitFinderFast.h"
#include "Tpc2dClusterFast.h"
#include "MpdTpcDigit.h"
#include "MpdTpcHit.h"

// ROOT/FairRoot Class Headers --------
#include <TClonesArray.h>
#include <TMath.h>
#include <TROOT.h>
#include <TSystem.h>
#include "FairLogger.h"

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
   if (!_pSecGeo) LOG(fatal) << ("TpcClusterHitFinderFast::TpcClusterHitFinderFast", "!!! Wrong geometry type !!!");
}

//__________________________________________________________________________

TpcClusterHitFinderFast::~TpcClusterHitFinderFast() {}

//__________________________________________________________________________

void TpcClusterHitFinderFast::FindHits()
{
   _nEvent = ioman->GetEntryNr();

   LOG(info) << "TpcClusterHitFinderFast::Exec started: Event " << _nEvent;
 
  /* 
     inline void EventClusters::Load(uint nSector, uint nRow, uint nPad, uint nTimeBin, float fAdc)
     - initial loading of digits
     - partial sorting of digits to pre-clusters 
     - charge correction at the boundary between pre-clusters
     
     EventClusters::DefineClusters(std::ostringstream& debugInfo)
     - transforming pre-clusters into clusters, where they are grouped with neighboring digits 
       until cluster closure
     - on cluster closure the local coordinate of the Hit is evaluated
  */
   std::ostringstream oss; 
   EventClusters *pEventClusters = new EventClusters(_nEvent);

   int nDigits = digiArray->GetEntriesFast();
   for (int i = 0; i < nDigits; ++i) {
      MpdTpcDigit *pdigit   = (MpdTpcDigit *)digiArray->UncheckedAt(i); // Input digit
      uint         nSector  = pdigit->GetSector();
      uint         nRow     = pdigit->GetRow();
      uint         nPad     = pdigit->GetPad();
      uint         nTimeBin = pdigit->GetTimeBin();
      float        fAdc     = pdigit->GetAdc();
      //int          nTrackId = pdigit->GetOrigin();
      //AdcHit       adcHit(nTimeBin, fAdc, nTrackId);

      pEventClusters->Load(nSector, nRow, nPad, nTimeBin, fAdc);
   }

   if( !pEventClusters->isEmpty() ) {
      pEventClusters->DefineClusters(oss);
      calcSector(pEventClusters);
   }

   delete pEventClusters;
/*
   int nS = _nSector, nR = _nRow;
   int nSector_ = nS;
   int nRow_    = nR;

   list<PadCluster *>                   lstPadClusters;
   vector<list<PadCluster *>::iterator> viPadClusters;

   LOG(info) << "Sector: " << nS << _nSector << nSector_;
   LOG(info) << "Row: " << nR << _nRow << nRow_;

   Int_t nPoints = digiArray->GetEntriesFast();
   for (int i = 0; i < nPoints; i++) {
      MpdTpcDigit *pdigit   = (MpdTpcDigit *)digiArray->UncheckedAt(i); // Input digit
      uint         nSector  = pdigit->GetSector();
      uint         nRow     = pdigit->GetRow();
      uint         nPad     = pdigit->GetPad();
      uint         nTimeBin = pdigit->GetTimeBin();
      float        fAdc     = pdigit->GetAdc();
      int          nTrackId = pdigit->GetOrigin();
      AdcHit       adcHit(nTimeBin, fAdc, nTrackId);
      
      LOG(info) << nSector_;
      if (nSector_ < 0) { // 1st occur only!
         nSector_ = nSector;
         pEventClusters->Add(new SectorClusters(nSector));
         LOG(info) << "SECTOR SHIT: " << nSector << nS << nSector_;
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

            assert(pEventClusters->getBack() != NULL);
            pEventClusters->getBack()->Add(pRowClusters);
            if (bSector && nS < 0) pEventClusters->Add(new SectorClusters(nSector));
         }
         viPadClusters.clear();
         lstPadClusters.clear();
         if (bSector) {
            if (nS < 0) {
               nRow_    = nRow;
               nSector_ = nSector;
            } else
               break;
         } else if (nR < 0) {
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
      //if (nR < 0) strRows[nRow_] = szDECADE[nRow_ % 10];
      assert(pEventClusters->getBack() != NULL);
      pEventClusters->getBack()->Add(pRowClusters);
      TpcClusterHitFinderFast::calcSector(pEventClusters);
   }

   viPadClusters.clear();
   lstPadClusters.clear(); */
}

void TpcClusterHitFinderFast::calcSector(const EventClusters *pEventClusters)
{

  /* 
     - store clusters in clusArray
     - transform ClusterHit from local to global coordinate
     - store transformed hits in hitArray
  */

   int nClusters = 0;
   int nHits  = 0;
   TVector3 globErr(0.05, 0.0, 0.2);
   TVector3 p3glob(0.0, 0.0, 0.0);
   const list<SectorClusters *> &rlstEventClusters = pEventClusters->getSectorClusters();

   // loop hierarchy: i0 - sector, i1 - row, i2 - cluster, i3 - padcluster, i4 - digit
   for (list<SectorClusters *>::const_iterator i0 = rlstEventClusters.begin(); i0 != rlstEventClusters.end(); ++i0) 
   {
      int   nSector = (*i0)->getSector();
      const list<RowClusters *> &rlstRowClusters = (*i0)->getRowClusters();
      for (list<RowClusters *>::const_iterator i1 = rlstRowClusters.begin(); i1 != rlstRowClusters.end(); ++i1)
      {
         int   nRow  = (*i1)->getRow();
         const list<Cluster *> &rlstClusters = (*i1)->getClusters();
         for (list<Cluster *>::const_iterator i2 = rlstClusters.begin(); i2 != rlstClusters.end(); ++i2) 
         {
           const Cluster *pCluster = *i2; 
           const list<const PadCluster *> &rlstPadClusters = (*i2)->getPadClusters();
           // ----- write clusters -----
           Tpc2dClusterFast *clus = new ((*clusArray)[nClusters]) Tpc2dClusterFast();
           clus->SetClusterID(nHits);
           clus->SetSector(nSector);
           clus->SetRow(nRow);
           clus->SetTimeBinMin(pCluster->getTimeBinMin());
           clus->SetTimeBinMax(pCluster->getTimeBinMax());
           clus->SetPadMin(pCluster->getPadMin());
           clus->SetPadMax(pCluster->getPadMax());
           for (list<const PadCluster *>::const_iterator i3 = rlstPadClusters.begin(); i3 != rlstPadClusters.end(); 
                ++i3) 
            {
               const vector<AdcHit> &rvAdcHits = (*i3)->getAdcHits();
               for (vector<AdcHit>::const_iterator i4 = rvAdcHits.begin(); i4 != rvAdcHits.end(); ++i4) { 
                AbstractTpcDigit *pDigit = new MpdTpcDigit(0, (*i3)->getPad(), nRow, i4->getTimeBin(), nSector, i4->getAdc());
                clus->AddDigit(pDigit);
               }
            }
           ++nClusters;

           // ----- write hits -----
           float fPad  = pCluster->getPad();
           float fTime = pCluster->getTime();   
           bool isPadAreaOuter = (nRow < _pSecGeo->ROW_COUNT[0]); 
           
           //  !!! AZ geometry formulas written in BaseTpcSectorGeo notation !!!
           //  AZ formula for zloc = TimeBin2Z(ftime) = return (fTimeBinMax - fTime - 0.5 + 0.037) / fZ2TimeBin;
           //  can be by his logic simplified to fZmax + (0.037 - 0.5 - fTime) * _pSecGeo->TIMEBIN_LENGTH * fV_DRIFT
           //  where he substitutes fZmax = 170 
           double currentPadWidth = _pSecGeo->PAD_WIDTH[isPadAreaOuter];
           double xloc = currentPadWidth * (fPad - _pSecGeo->PAD_COUNT[nRow] + 0.5); 
           double yloc;
           if (nRow < _pSecGeo->ROW_COUNT[0]) yloc = _pSecGeo->PAD_HEIGHT[0] * (nRow + 0.5);
           else yloc = _pSecGeo->YPADAREA_LENGTH[0] + _pSecGeo->PAD_HEIGHT[1] * (nRow - _pSecGeo->ROW_COUNT[0] + 0.5);       
           double zloc  = 170 + (0.037 - 0.5 -fTime) * _pSecGeo->TIMEBIN_LENGTH * fV_DRIFT; 
           if (nSector >= _pSecGeo->SECTOR_COUNT_HALF) zloc = -zloc;

           TVector3 p3loc(xloc, yloc, zloc);

           // AZ geometry: to be refactored to BaseTpcSectorGeo
           _pSecGeo->Local2Global(nSector, p3loc, p3glob);

           /*
            Double_t xyz1[3], xyz2[3];

            xyzLoc.GetXYZ(xyz1);
            Local2Global(iSec, xyz1, xyz2);
            xyzGlob.SetXYZ(xyz2[0], xyz2[1], xyz2[2]);

            xyzGlob[2]     = xyzLoc[2];
            Double_t phSec = iSec * fDphi;
            Double_t cosPh = TMath::Cos(phSec);
            Double_t sinPh = TMath::Sin(phSec);
            Double_t x     = xyzLoc[1] + fYsec[0];
            Double_t y     = xyzLoc[0];
            xyzGlob[0]     = x * cosPh - y * sinPh;
            xyzGlob[1]     = x * sinPh + y * cosPh;

           */

           //LOG(info) << "Y DIFFERENCE:" << yloc-ylocAZ << " YLOC: " << yloc << " YLOC_AZ: " << ylocAZ;

           int padID = (nSector % _pSecGeo->SECTOR_COUNT_HALF) | (nRow << 5); 
           
           MpdTpcHit *hit = new ((*hitArray)[nHits]) MpdTpcHit(padID, p3glob, globErr, nClusters-1); 
           hit->SetLayer(nRow);
           hit->SetEnergyLoss(pCluster->getAdcSum());
           hit->SetNdigits(rlstPadClusters.size());
           hit->SetPadCoordinate(fPad);
           hit->SetTimeBinCoordinate(fTime);
           hit->SetPad(int(fPad));
           hit->SetBin(int(fTime));
           hit->SetStep(_pSecGeo->PAD_HEIGHT[isPadAreaOuter]);
           hit->SetModular(1); 
           hit->SetLocalPosition(p3loc); // point position

           ++nHits;
         } 
      }
   }

   LOG(info) << "Total number of clusters is: " << nClusters;
   LOG(info) << "Total number of hits is: " << nHits;

   
   

   
   

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
            //LOG(info) << "Hit pad is: " << fPad;
            float fTime = pCluster->getTime();
            
            // --------------------------------TpcSectorGeoAZ--------------------------------------------------
            int      padID = _pSecGeo->PadID(nSect % _pSecGeo->NofSectors(), nRow);
            double   xloc  = _pSecGeo->Pad2Xloc(fPad, nRow);
            double   yloc  = _pSecGeo->LocalPadPosition(padID).Y();
            double   zloc  = _pSecGeo->TimeBin2Z(fTime);
            TVector3 p3loc(xloc, yloc, zloc);
            if (nSect >= _pSecGeo->NofSectors()) p3loc[2] = -p3loc[2];
            _pSecGeo->Local2Global(nSect, p3loc, p3glob);
            
            // write cluster(s) 
            
            //nHit                   = hitArray->GetEntriesFast();
            
            //nClus                   = clusArray->GetEntriesFast();
            //LOG(info) << "Clusnumber is: " << ++nClus;
            //Tpc2dClusterFast *clus = new ((*clusArray)[nClus++]) Tpc2dClusterFast(); // Add new cluster
            /*clus->SetClusterID(nHit);
            clus->SetSector(nSect);
            clus->SetRow(nRow);
            clus->SetTimeBinMin(pCluster->getTimeBinMin());
            clus->SetTimeBinMax(pCluster->getTimeBinMax());
            clus->SetPadMin(pCluster->getPadMin());
            clus->SetPadMax(pCluster->getPadMax());
            /*
            // write digits into cluster
            //-vector<pair<int, float>>        vnTrackId;
            const list<const PadCluster *> &rlstPadClusters = (*i2)->getPadClusters();
            for (list<const PadCluster *>::const_iterator i3 = rlstPadClusters.begin(); i3 != rlstPadClusters.end();
                 i3++) {
               const vector<AdcHit> &rvAdcHits = (*i3)->getAdcHits();
               for (vector<AdcHit>::const_iterator i4 = rvAdcHits.begin(); i4 != rvAdcHits.end(); i4++) { // Digits
               //-   vnTrackId = i4->getTrackID();
               //-AbstractTpcDigit *pDigit =
               //-new MpdTpcDigit(vnTrackId[0].first, (*i3)->getPad(), nRow, i4->getTimeBin(), nSect, i4->getAdc());
                  AbstractTpcDigit *pDigit = new MpdTpcDigit(0, (*i3)->getPad(), nRow, i4->getTimeBin(), nSect, i4->getAdc());
                  clus->AddDigit(pDigit);
               }
            }
            
            // write hit(s) 
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
            hit->SetPadCoordinate(fPad);
            hit->SetTimeBinCoordinate(fTime);
            hit->SetNdigits(rlstPadClusters.size());

            for (vector<pair<int, float>>::const_iterator i5 = vnTrackId.begin(); i5 != vnTrackId.end(); i5++) {
               hit->AddTrackID((*i5).first, (*i5).second); // TrackID
            }
            */
         }
      }
   }
}
ClassImp(TpcClusterHitFinderFast);

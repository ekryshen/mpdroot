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
   if (!_pSecGeo) Fatal("TpcClusterHitFinderFast::TpcClusterHitFinderFast", " !!! Wrong geometry type !!! ");
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
   EventClusters     *pEventClusters = new EventClusters(_nEvent);

   int nDigits = digiArray->GetEntriesFast();
   for (int i = 0; i < nDigits; ++i) {
      MpdTpcDigit *pdigit   = (MpdTpcDigit *)digiArray->UncheckedAt(i); // Input digit
      uint         nSector  = pdigit->GetSector();
      uint         nRow     = pdigit->GetRow();
      uint         nPad     = pdigit->GetPad();
      uint         nTimeBin = pdigit->GetTimeBin();
      float        fAdc     = pdigit->GetAdc();

      pEventClusters->Load(nSector, nRow, nPad, nTimeBin, fAdc);
   }

   if (!pEventClusters->isEmpty()) {
      pEventClusters->DefineClusters(oss);
      calcSector(pEventClusters);
   }

   delete pEventClusters;
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
            const list<const PadCluster *> &rlstPadClusters = (*i2)->getPadClusters();
            for (list<const PadCluster *>::const_iterator i3 = rlstPadClusters.begin(); i3 != rlstPadClusters.end();
                 i3++) {
               const vector<AdcHit> &rvAdcHits = (*i3)->getAdcHits();
               for (vector<AdcHit>::const_iterator i4 = rvAdcHits.begin(); i4 != rvAdcHits.end(); i4++) { // Digits
                  AbstractTpcDigit *pDigit =
                     new MpdTpcDigit(0, (*i3)->getPad(), nRow, i4->getTimeBin(), nSect, i4->getAdc());
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
         }
      }
   }
}
ClassImp(TpcClusterHitFinderFast);

#include "TpcClusterHitFinderFast.h"
#include "Tpc2dClusterFast.h"
#include "MpdTpcDigit.h"
#include "MpdTpcHit.h"

// ROOT/FairRoot Class Headers --------
#include <TClonesArray.h>
#include <TMath.h>
#include "FairLogger.h"

// C/C++ Headers ----------------------
#include <iostream>
#include <math.h>
#include <vector>

#define fV_DRIFT 0.0055 // cm/ns

using namespace std;
using namespace tpcClustering;

//__________________________________________________________________________

void TpcClusterHitFinderFast::FindHits()
{
   int nEvent = ioman->GetEntryNr();

   LOG(info) << "TpcClusterHitFinderFast::Exec started: Event " << nEvent;

   /*
      void EventClusters::Load(uint nSector, uint nRow, uint nPad, uint nTimeBin, float fAdc, int nTrackId)
      - initial loading of digits
      - partial sorting of digits to pre-clusters
      - charge correction at the boundary between pre-clusters

      void EventClusters::DefineClusters(std::ostringstream& debugInfo)
      - converting pre-clusters into clusters, where they are grouped with neighboring digits
        until cluster closure
      - on cluster closure the local coordinate of the Hit is evaluated
   */
   std::ostringstream oss;
   EventClusters     *pEventClusters = new EventClusters(nEvent);

   int nDigits = digiArray->GetEntriesFast();
   for (int i = 0; i < nDigits; ++i) {
      MpdTpcDigit *pdigit   = (MpdTpcDigit *)digiArray->UncheckedAt(i); // Input digit
      uint         nSector  = pdigit->GetSector();
      uint         nRow     = pdigit->GetRow();
      uint         nPad     = pdigit->GetPad();
      uint         nTimeBin = pdigit->GetTimeBin();
      float        fAdc     = pdigit->GetAdc();
      int          nTrackId = pdigit->GetOrigin();

      pEventClusters->Load(nSector, nRow, nPad, nTimeBin, fAdc, nTrackId);
   }

   if (!pEventClusters->isEmpty()) {
      pEventClusters->DefineClusters(oss);
      transform(pEventClusters);
   }

   delete pEventClusters;
}

void TpcClusterHitFinderFast::transform(const EventClusters *pEventClusters)
{

   /*
      - store clusters in clusArray
      - transform ClusterHit from local to global coordinate
      - store transformed hits in hitArray
   */

   int                           nClusters         = 0;
   int                           nHits             = 0;
   const list<SectorClusters *> &rlstEventClusters = pEventClusters->getSectorClusters();

   // loop hierarchy: i0 - sector, i1 - row, i2 - cluster, i3 - padcluster, i4 - digit, i5 - trackID
   for (list<SectorClusters *>::const_iterator i0 = rlstEventClusters.begin(); i0 != rlstEventClusters.end(); ++i0) {
      int                        nSector         = (*i0)->getSector();
      const list<RowClusters *> &rlstRowClusters = (*i0)->getRowClusters();
      for (list<RowClusters *>::const_iterator i1 = rlstRowClusters.begin(); i1 != rlstRowClusters.end(); ++i1) {
         int                    nRow         = (*i1)->getRow();
         const list<Cluster *> &rlstClusters = (*i1)->getClusters();
         for (list<Cluster *>::const_iterator i2 = rlstClusters.begin(); i2 != rlstClusters.end(); ++i2) {
            const Cluster                  *pCluster        = *i2;
            const list<const PadCluster *> &rlstPadClusters = (*i2)->getPadClusters();
            vector<pair<int, float>>        vnTrackId;
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
                 ++i3) {
               const vector<AdcHit> &rvAdcHits = (*i3)->getAdcHits();
               for (vector<AdcHit>::const_iterator i4 = rvAdcHits.begin(); i4 != rvAdcHits.end(); ++i4) {
                  vnTrackId                = i4->getTrackID();
                  AbstractTpcDigit *pDigit = new MpdTpcDigit(vnTrackId[0].first, (*i3)->getPad(), nRow,
                                                             i4->getTimeBin(), nSector, i4->getAdc());
                  clus->AddDigit(pDigit);
               }
            }
            ++nClusters;

            // ----- write hits -----
            float                     fPad           = pCluster->getPad();
            float                     fTime          = pCluster->getTime();
            BaseTpcSectorGeo::PadArea currentPadArea = (nRow < secGeo->ROW_COUNT[BaseTpcSectorGeo::inner])
                                                          ? currentPadArea = BaseTpcSectorGeo::inner
                                                          : BaseTpcSectorGeo::outer;

            //  !!! AZ geometry transformations written in BaseTpcSectorGeo notation !!!
            //  AZ expression for zloc = TimeBin2Z(ftime) is equivalent to the one below
            double zloc = 170 + (0.037 - 0.5 - fTime) * secGeo->TIMEBIN_LENGTH * fV_DRIFT;
            if (nSector >= secGeo->SECTOR_COUNT_HALF) zloc = -zloc;
            TVector2 p2loc(secGeo->PadRow2Local((double)fPad, (double)(0.5 + nRow)));
            TVector3 p3loc(p2loc.X(), p2loc.Y(), zloc);

            TVector3 p3glob(p3loc.Y() + secGeo->YPADAREA_LOWEREDGE, p3loc.X(), p3loc.Z());
            double   phSec = nSector * secGeo->SECTOR_PHI_RAD;
            p3glob.RotateZ(phSec);

            /* write hit(s) */
            int        padID = (nSector % secGeo->SECTOR_COUNT_HALF) | (nRow << 5);
            MpdTpcHit *hit = new ((*hitArray)[nHits]) MpdTpcHit(padID, p3glob, TVector3(0.05, 0.0, 0.2), nClusters - 1);
            hit->SetLayer(nRow);
            hit->SetLocalPosition(p3loc);
            hit->SetEnergyLoss(pCluster->getAdcSum());
            hit->SetModular(1);
            hit->SetPad(int(fPad));
            hit->SetBin(int(fTime));
            hit->SetPadCoordinate(fPad);
            hit->SetTimeBinCoordinate(fTime);
            hit->SetNdigits(rlstPadClusters.size());
            hit->SetStep(secGeo->PAD_HEIGHT[currentPadArea]);
            for (vector<pair<int, float>>::const_iterator i5 = vnTrackId.begin(); i5 != vnTrackId.end(); ++i5)
               hit->AddTrackID((*i5).first, (*i5).second);
            ++nHits;
         }
      }
   }
}
ClassImp(TpcClusterHitFinderFast);

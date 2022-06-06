//------------------------------------------------------------------------------------------------------------------------
#include <iostream>

#include <TFile.h>
#include <TMath.h>
#include <TClonesArray.h>

#include <FairLogger.h>
#include "MpdTofHit.h"
#include "MpdTof.h"
#include "MpdTofPoint.h"
#include "MpdTofGeoUtils.h"

#include "MpdTofHitProducerQA.h"
using namespace std;
//------------------------------------------------------------------------------------------------------------------------
MpdTofHitProducerQA::MpdTofHitProducerQA(const char *flnm, bool isEndcap) : fFlnm(flnm), fIsEndcap(isEndcap)
{
   fList.SetOwner(); //  all objects will be deleted whenever the collection itself is delete.

   Add(hOccup = new TH1D(mangling("Occupancy"), "occupancy per strips;occupancy;Events", 100, -0.5, 99.5));
   Add(hMergedTimes =
          new TH2D(mangling("MergedTimes"), "Merged hits on strip times test;faster hit time, ns;slower hit time, ns",
                   1000, 0., 150., 1000, 0., 150.));
   Add(hHitPointPerEvent = new TH2D(mangling("NHitNPointPerEvent"), ";N points per event;N hits per event", 1000, 0.5,
                                    4000.5, 1000, 0.5, 1500.5));
   Add(hPointDistanceSame = new TH2D(mangling("TwoPointDistanceSame"),
                                     "distances between a pair points(same tids);distance, cm;time between hits, ns",
                                     1000, 0., 50., 1000, 0., 50.));
   Add(hPointDistanceDiff =
          new TH2D(mangling("TwoPointDistanceDiff"),
                   "distances between a pair points(different tids);distance, cm;time between hits, ns", 1000, 0., 50.,
                   1000, 0., 50.));

   Add(effHitGap2 = new TEfficiency(mangling("efficiencyHitsMiddle"), "Efficiency of middle gap hits;R, cm;efficiency",
                                    1000, -0.1, 1.));
   Add(effHitGap13 = (TEfficiency *)effHitGap2->Clone(mangling("efficiencyHitsOuter")));
   effHitGap13->SetTitle("Efficiency of outer gaps hits;R, cm;efficiency");
   Add(effCrossHit = (TEfficiency *)effHitGap2->Clone(mangling("efficiencyCrossHits")));
   effCrossHit->SetTitle("Efficiency of cross hits;R, cm;efficiency");

   // geometrical tests
   Add(hDistance =
          new TH1D(mangling("geoDistance"),
                   "Distance between strips(minus-left side, plus-right side);Distance, cm;Side", 1000, -100., 100.));
   Add(h2Strips = new TH2D(mangling("geoStripsMap"), ";Z, cm;#phi, rads", 2000, -300., 300., 500, -3.5, 3.5));
   Add(h2Detectors = new TH2D(mangling("geoDetectorMap"), ";Z, cm;#phi, rads", 2000, -300., 300., 500, -3.5, 3.5));
   Add(hNeighborPair = new TH2D(mangling("geoNeighborPair"), "Neighbor strip pairs test;stripID1;stripID2", 75, -0.5,
                                74.5, 75, -0.5, 74.5));
   Add(hEtaPhi = new TH2D(mangling("geoEtaPhi"), ";#eta;#phi, degree", 1000, -1.5, 1.5, 1000, -181., 181.));
   Add(hRZ = new TH2D(mangling("geoRZ"), ";Z, cm;R, cm", 1000, -300., 300., 1000, 150., 160.));

   // hit smearing tests
   Add(hXYSmeared = new TH2D(mangling("XYSmeared_ZPhi"), "single hits;#Delta_{Z}, cm;#Delta_{#phi}, cm", 1000, -0.7,
                             0.7, 1000, -5., 5.));
   Add(hXYSmeared2 =
          new TH2D(mangling("XYSmeared_DevR"), "single hits;#Delta, cm;#DeltaRperp, cm", 1000, 0., 2., 1000, -2., 2.));
   Add(hXYSmearedCross = new TH2D(mangling("XYSmearedCross_ZPhi"), "cross hits;#Delta_{Z}, cm;#Delta_{#phi}, cm", 1000,
                                  -2., 2., 1000, -1., 1.));
   Add(hXYSmearedCross2 = new TH2D(mangling("XYSmearedCross_DevR"), "cross hits;#Delta, cm;#DeltaRperp, cm", 1000, 0.,
                                   3., 1000, -3., 3.));

   Add(hHitPositionInsideStrip =
          new TH2D(mangling("HitPositionInsideStrip"),
                   "Hit position deviation from strip center(single hits);#Delta_{#phi}, cm;#Delta_{Z}, cm", 1000, -50.,
                   50., 1000, -0.5, 0.5));
   Add(hMCPositionInsideStrip =
          new TH2D(mangling("MCPositionInsideStrip"),
                   "MC point position deviation from strip center(single hits);#Delta_{#phi}, cm;#Delta_{Z}, cm", 1000,
                   -50., 50., 1000, -5., 5.));
   Add(hPointXZOrigin = new TH2D(mangling("PointXZOrigin"), "MC point XZ position at origin LRS;Z, cm;X, cm", 1000, -1.,
                                 1., 1000, -50., 50.));
   Add(hPointYZOrigin = new TH2D(mangling("PointYZOrigin"), "MC point YZ position at origin LRS;Z, cm;Y, cm", 1000, -1.,
                                 1., 1000, -0.2, 0.2));
   Add(hHitXZOrigin = new TH2D(mangling("HitXZOrigin"), "Hit XZ position at origin LRS;X, cm;Z, cm", 1000, -50., 50.,
                               1000, -0.05, 0.05));
   Add(hHitYZOrigin = new TH2D(mangling("HitYZOrigin"), "Hit YZ position at origin LRS;Y, cm;Z, cm", 1000, -0.1, 0.1,
                               1000, -0.05, 0.05));

   Add(hDevHitXZOrigin = new TH2D(mangling("DevHitXZOrigin"),
                                  "Hit XZ position deviation from MC point position at origin LRS;X, cm;Z, cm", 1000,
                                  -10., 10., 1000, -1., 1.));
   Add(hDevHitYZOrigin = new TH2D(mangling("DevHitYZOrigin"),
                                  "Hit YZ position deviation from MC point position at origin LRS;Y, cm;Z, cm", 1000,
                                  -0.1, 0.1, 1000, -1., 1.));

   Add(hXZCentralDetector = new TH2D(mangling("XZCentralDetector"), "central detector hits.;Z, cm;X, cm", 1000, -10.,
                                     40., 1000, -50., 50.));
   Add(hYZCentralDetector = new TH2D(mangling("YZCentralDetector"), "central detector hits.;Z, cm;Y, cm", 1000, -10.,
                                     40., 1000, 149., 152.));
   Add(hXZDetector =
          new TH2D(mangling("XYDetector"), "central detector hits.;X, cm;Y, cm", 1000, -40., 40., 1000, -0.6, 0.6));
   Add(hYZDetector =
          new TH2D(mangling("YZDetector"), "central detector hits.;Z, cm;Y, cm", 1000, -0.05, 0.05, 1000, -0.6, 0.6));
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofHitProducerQA::Finish()
{
   LOG(debug2) << "[MpdTofHitProducerQA::Finish] Update  " << fFlnm.Data() << " file. ";
   auto  ptr = gFile;
   TFile file(fFlnm.Data(), "RECREATE");
   fList.Write();
   file.Close();
   gFile = ptr;
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofHitProducerQA::PointDistanceTest(const TClonesArray *aMcPoints)
{
   TVector3 pos1, pos2;
   for (Int_t i = 0, n1 = aMcPoints->GetEntriesFast(); i < n1; i++) // cycle by tof hits
   {
      auto point1 = (MpdTofPoint *)aMcPoints->At(i);
      point1->Position(pos1);

      for (Int_t j = i, n2 = aMcPoints->GetEntriesFast(); j < n2; j++) {
         if (i == j) continue;

         auto point2 = (MpdTofPoint *)aMcPoints->At(j);
         point2->Position(pos2);
         double dT = fabs(point1->GetTime() - point2->GetTime());

         if (point1->GetTrackID() == point2->GetTrackID())
            hPointDistanceSame->Fill((pos1 - pos2).Mag(), dT);
         else
            hPointDistanceDiff->Fill((pos1 - pos2).Mag(), dT);
      }
   }
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofHitProducerQA::FillDetectorsMap(const TVector3 &A, const TVector3 &B, const TVector3 &C,
                                           const TVector3 &D) // detector edges
{
   h2Detectors->Fill(A.Z(), A.Phi());
   h2Detectors->Fill(B.Z(), B.Phi());
   h2Detectors->Fill(C.Z(), C.Phi());
   h2Detectors->Fill(D.Z(), D.Phi());
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofHitProducerQA::FillStripsMap(const TVector3 &A, const TVector3 &B, const TVector3 &C,
                                        const TVector3 &D) // strip edges
{
   h2Strips->Fill(A.Z(), A.Phi());
   h2Strips->Fill(B.Z(), B.Phi());
   h2Strips->Fill(C.Z(), C.Phi());
   h2Strips->Fill(D.Z(), D.Phi());
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofHitProducerQA::Point2HitSmearingTest(const TVector3 &mcPosition, const TVector3 &hitPosition)
{
   double delta, deltaZ, deltaR, deltaPhi;
   MpdTof::GetDelta(mcPosition, hitPosition, delta, deltaZ, deltaR, deltaPhi);

   hXYSmeared->Fill(deltaZ, deltaPhi);
   hXYSmeared2->Fill(delta, deltaR);

   hEtaPhi->Fill(mcPosition.Eta(), mcPosition.Phi() * TMath::RadToDeg());
   hRZ->Fill(mcPosition.Z(), mcPosition.Perp());
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofHitProducerQA::HitGapEfficiencyTest(bool fired, Double_t distance, Int_t gap)
{
   if (2 == gap)
      effHitGap2->Fill(fired, distance);
   else
      effHitGap13->Fill(fired, distance);
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofHitProducerQA::PositionInsideStripTest(const TVector3 &stripCenter, const TVector3 &mcPosition,
                                                  const TVector3 &hitPosition)
{
   double hitMeanR = (stripCenter.Perp() + hitPosition.Perp()) / 2.;
   hHitPositionInsideStrip->Fill((stripCenter.Phi() - hitPosition.Phi()) * hitMeanR, stripCenter.Z() - hitPosition.Z());
   // double mcMeanR = (stripCenter.Perp() + mcPosition.Perp())/2.;
   hMCPositionInsideStrip->Fill((stripCenter.Phi() - mcPosition.Phi()) * hitMeanR, stripCenter.Z() - mcPosition.Z());
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofHitProducerQA::RotationToOriginTest(const TGeoCombiTrans &matrix, const TVector3 &mcPosition,
                                               const TVector3 &hitPosition)
{
   // rotate stip to origin LRS
   Double_t localHit[3], localPoint[3], master[3] = {mcPosition.X(), mcPosition.Y(), mcPosition.Z()};
   matrix.MasterToLocal(master, localPoint);
   hPointXZOrigin->Fill(localPoint[2], localPoint[0]);
   hPointYZOrigin->Fill(localPoint[2], localPoint[1]);

   master[0] = hitPosition.X();
   master[1] = hitPosition.Y();
   master[2] = hitPosition.Z();
   matrix.MasterToLocal(master, localHit);
   hHitXZOrigin->Fill(localHit[0], localHit[2]);
   hHitYZOrigin->Fill(localHit[1], localHit[2]);

   hDevHitXZOrigin->Fill(localHit[0] - localPoint[0],
                         localHit[2] - localPoint[2]); // (Xhit - Xpoint) vs (Zhit - Zpoint)
   hDevHitYZOrigin->Fill(localHit[1] - localPoint[1],
                         localHit[2] - localPoint[2]); // (Yhit - Ypoint) vs (Zhit - Zpoint)
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofHitProducerQA::CentralDetectorTest(Int_t suid, const TVector3 &hitPosition)
{
   if (MpdTofPoint::GetSector(suid) == 4 &&
       MpdTofPoint::GetDetector(suid) ==
          11) // select unrotated detector (strips along X axis, perp. to strip along Y axis)
   {
      hXZCentralDetector->Fill(hitPosition.Z(), hitPosition.X());
      hYZCentralDetector->Fill(hitPosition.Z(), hitPosition.Y());
   }

   // back transform to the central gap LRS
   const auto &matrix = MpdTofGeoUtils::Instance()->FindStrip(MpdTofPoint::SetCentralGap(suid))->fMatrix;
   Double_t    local[3], master[3] = {hitPosition.X(), hitPosition.Y(), hitPosition.Z()};
   matrix.MasterToLocal(master, local);

   hXZDetector->Fill(local[0], local[1]);
   hYZDetector->Fill(local[2], local[1]);
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofHitProducerQA::FillCrossHitEfficiency(bool fired, Double_t distance, size_t gap, const TVector3 &position,
                                                 const TVector3 &smearedPosition)
{
   effCrossHit->Fill(fired, distance);

   if (fired) {
      double delta, deltaZ, deltaR, deltaPhi;
      MpdTof::GetDelta(position, smearedPosition, delta, deltaZ, deltaR, deltaPhi);

      hXYSmearedCross->Fill(deltaZ, deltaPhi);
      hXYSmearedCross2->Fill(delta, deltaR);
   }
}
//------------------------------------------------------------------------------------------------------------------------

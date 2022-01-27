//------------------------------------------------------------------------------------------------------------------------
#ifndef __MPD_TOF_HITPRODUCER_QA_H
#define __MPD_TOF_HITPRODUCER_QA_H 1

#include <TList.h>
#include <TString.h>
#include <TVector3.h>
#include <TEfficiency.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TGeoMatrix.h>
//------------------------------------------------------------------------------------------------------------------------
class TClonesArray;
class MpdTofHitProducerQA {
   TList fList;

   TH1D *hOccup, *htR, *hDistance;
   TH2D *hMergedTimes, *hNeighborPair, *hXYSmeared, *hXYSmeared2, *hXYSmearedCross, *hXYSmearedCross2, *hEtaPhi,
      *h2Strips, *h2Detectors, *hRZ;
   TH2D *hHitPointPerEvent = nullptr, *hPointDistanceSame = nullptr, *hPointDistanceDiff = nullptr,
        *hDevHitYZOrigin = nullptr, *hDevHitXZOrigin = nullptr;
   TH2D *hHitPositionInsideStrip = nullptr, *hMCPositionInsideStrip = nullptr, *hPointXZOrigin = nullptr,
        *hPointYZOrigin = nullptr, *hHitXZOrigin = nullptr, *hHitYZOrigin = nullptr;
   TH2D *hXZCentralDetector = nullptr, *hYZCentralDetector = nullptr, *hXZDetector = nullptr, *hYZDetector = nullptr;

   TEfficiency *effHitGap2, *effHitGap13, *effCrossHit;

   TString fFlnm;
   bool    fIsEndcap;

   const char *mangling(const char *name)
   {
      static TString nm;
      nm = fIsEndcap ? "LsETOF_" : "LsTOF_";
      nm += name;
      return nm.Data();
   }
   void Add(TH1 *hist)
   {
      hist->SetDirectory(nullptr);
      fList.Add(hist);
   }
   void Add(TEfficiency *hist)
   {
      hist->SetDirectory(nullptr);
      fList.Add(hist);
   }

public:
   MpdTofHitProducerQA(const char *flnm, bool isEndcap);
   void Finish();

   void FillDetectorsMap(const TVector3 &A, const TVector3 &B, const TVector3 &C, const TVector3 &D);
   void FillStripsMap(const TVector3 &A, const TVector3 &B, const TVector3 &C, const TVector3 &D);

   void FillCrossHitEfficiency(bool fired, Double_t distance, size_t gap, const TVector3 &position,
                               const TVector3 &smearedPosition);

   void Point2HitSmearingTest(const TVector3 &mcPosition, const TVector3 &hitPosition);
   void HitGapEfficiencyTest(bool fired, Double_t distance, Int_t gap);
   void PointDistanceTest(const TClonesArray *);
   void PositionInsideStripTest(const TVector3 &stripCenter, const TVector3 &mcPosition, const TVector3 &hitPosition);
   void RotationToOriginTest(const TGeoCombiTrans &matrix, const TVector3 &mcPosition, const TVector3 &hitPosition);
   void CentralDetectorTest(Int_t suid, const TVector3 &hitPosition);

   inline void FillOccupancy(Double_t occupancy) { hOccup->Fill(occupancy); }
   inline void FillMergedTimes(Double_t time1, Double_t time2) { hMergedTimes->Fill(time1, time2); }
   inline void FillDistanceToStrip(Double_t distance) { hDistance->Fill(distance); }
   inline void FillNeighborPairs(Int_t suid1, Int_t suid2) { hNeighborPair->Fill(suid1, suid2); }
   inline void FillNPointsHits(Double_t points, Double_t hits) { hHitPointPerEvent->Fill(points, hits); }
};
//------------------------------------------------------------------------------------------------------------------------
#endif

//------------------------------------------------------------------------------------------------------------------------
/// \class MpdTofHitProducer
///
/// \brief
/// \author Sergei Lobastov (LHE, JINR, Dubna)
//------------------------------------------------------------------------------------------------------------------------
#include <assert.h>

#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TEfficiency.h>
#include <TRandom2.h>
#include <TClonesArray.h>

#include "MpdMCTrack.h"
#include "FairLogger.h"

#include "MpdTofUtils.h"
#include "MpdTofPoint.h"
#include "MpdTofHit.h"
#include "MpdTofGeoUtils.h"
#include "MpdTofHitProducerQA.h"

#include "MpdTofHitProducer.h"
//V
#include "FairMCEventHeader.h"

using namespace std;

ClassImp(MpdTofHitProducer);
//------------------------------------------------------------------------------------------------------------------------
MpdTofHitProducer::MpdTofHitProducer(const char *name, Bool_t useMCdata, Int_t verbose, Bool_t test, const char *flnm,
                                     double stripLength)
   : MpdTofHitProducerIdeal(name, useMCdata, verbose, test, true, flnm, false), fStripLength(stripLength)
{
   pRandom = new TRandom2;
}
//------------------------------------------------------------------------------------------------------------------------
MpdTofHitProducer::~MpdTofHitProducer()
{
   delete pRandom;
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofHitProducer::AddParameters(TString &buf) const
{
   MpdTofHitProducerIdeal::AddParameters(buf);
   buf += TString::Format(", fTimeSigma=%.5g ns", fTimeSigma);
   buf += TString::Format(", fErrX=%.4g cm", fErrX);
   buf += TString::Format(", fErrZ=%.4g cm", fErrZ);
}
//------------------------------------------------------------------------------------------------------------------------
InitStatus MpdTofHitProducer::Init()
{
   InitStatus status = MpdTofHitProducerIdeal::Init();

   MpdTofGeoUtils::Instance()->ParseTGeoManager(pQA, true); //, "geometryTree.txt"); // forced

   MpdTofGeoUtils::Instance()->FindNeighborStrips(0.8, pQA, // 0.8 [cm] <--- thresh. distance between neighbor strips
                                                  true);    // forced

//V
   if (fUseMCData) {
      fEventHeaderB = (FairMCEventHeader *)FairRootManager::Instance()->GetObject("MCEventHeader.");
assert(fEventHeaderB);
   }

   LOG(info) << "[MpdTofHitProducer::Init] Initialization finished succesfully.";

   return status;
}
//------------------------------------------------------------------------------------------------------------------------
Bool_t MpdTofHitProducer::IsHitCreated(Double_t value, Int_t gap) // value - distance to the strip edge [cm]
{
   constexpr Double_t slope      = (0.98 - 0.95) / 0.2;
   Double_t           efficiency = (value > 0.2) ? 0.98 : (0.95 + slope * value);

   //-------------------------------------
   // 99% ---------
   //              \.
   //               \.
   //                \.
   // 95%             \.
   //  <-----------|--|
   //            0.2  0. value
   //-------------------------------------
   if (gap == 1 || gap == 3) efficiency /= 2.; // reduce efficiency on outer strip gap

   if (pRandom->Rndm() < efficiency) return true;
   return false;
}
//------------------------------------------------------------------------------------------------------------------------
Bool_t MpdTofHitProducer::IsCrossHitCreated(Double_t value, Int_t gap) // value - distance to the strip edge  [cm]
{
   constexpr Double_t slope      = (0.3 - 0.0) / 0.5;
   Double_t           efficiency = (value > 0.5) ? 0. : (0.3 - slope * value);

   //-------------------------------------
   // 30%               /
   //                  /
   //                 /
   //                /
   // 0%            /
   //  <-----------|----|
   //            0.5    0. value
   //-------------------------------------

   if (efficiency == 0.) return false;

   if (gap == 1 || gap == 3) efficiency /= 2.; // reduce efficiency on outer strip gap

   if (pRandom->Rndm() < efficiency) return true;
   return false;
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofHitProducer::Exec(Option_t *option)
{
   aTofHits->Clear();

   TVector3 mcPosition, hitPosition, hitPosError;
   int      nSingleHits = 0, nCrossHits = 0;

   auto SmearingAlongStrip = [this](const TGeoCombiTrans &gap, const TGeoCombiTrans &centralGap, const TVector3 &point,
                                    TVector3 &hit, TVector3 &error) {
      // rotate stip to origin LRS (dX = fStripLength, dY = 0, dZ = strip step)
      Double_t local[3], master[3] = {point.X(), point.Y(), point.Z()};
      gap.MasterToLocal(master, local);

      // smearing the hit position along x axis
      double Xsmeared, x = local[0];
      size_t nTries = 0;
      do {
         Xsmeared = pRandom->Gaus(x, fErrX);
         if (++nTries > 100) return false;
      } while (fabs(Xsmeared) > fStripLength / 2.); // check strip boundary exit

      // shift and rotate back to MRS(shift by Y from gap to central gap)
      local[0] = Xsmeared, local[2] = 0.; // smearing along X axis, shift to the strip center along Z axis(at LRS)
      centralGap.LocalToMaster(local, master);
      hit.SetXYZ(master[0], master[1], master[2]);

      // rotate error vector
      local[0] = fErrX;
      local[1] = 0.;
      local[2] = fErrZ;
      centralGap.LocalToMaster(local, master);
      error.SetXYZ(master[0], master[1], master[2]);

      return true;
   };

   if (fUseMCData) {

//V
     b_mc = fEventHeaderB->GetB(); //[fm]
     t0res_mc = 0.;
     if (b_mc > 0 && b_mc < 17)
       {
	 t0res_mc = 1.52911e-002 - 
	            5.96912e-005 * b_mc +
	            1.78144e-004 * b_mc * b_mc -
	            3.39944e-005 * b_mc * b_mc * b_mc +
	            2.46001e-006 * b_mc * b_mc * b_mc * b_mc; //[ns]
       }
//V

      for (Int_t pointIndex = 0, nTofPoint = aMcPoints->GetEntriesFast(); pointIndex < nTofPoint;
           pointIndex++) // cycle by TOF points
      {
         auto pPoint = (MpdTofPoint *)aMcPoints->UncheckedAt(pointIndex);

         if (fVerbose > 2) pPoint->Print("");

         Int_t tid  = pPoint->GetTrackID();
         Int_t suid = pPoint->GetDetectorID();
         Int_t gap  = pPoint->GetGap();

//V
         Double_t time = pRandom->Gaus(pPoint->GetTime(), t0res_mc);
         time = pRandom->Gaus(time, fTimeSigma); // default 100 ps
//V

         pPoint->Position(mcPosition);

         auto strip        = MpdTofGeoUtils::Instance()->FindStrip(suid);
         auto centralStrip = MpdTofGeoUtils::Instance()->FindStrip(MpdTofPoint::SetCentralGap(suid));
         if (!SmearingAlongStrip(strip->fMatrix, centralStrip->fMatrix, mcPosition, hitPosition, hitPosError)) {
            strip->Dump(" [MpdTofHitProducer::Exec] -E- Invalid Rotation matrix.");
            continue;
         }

         LStrip::Side_t side;
         Double_t       distance = strip->MinDistanceToEdge(&mcPosition, side); // [cm]

         bool stripFired;
         if ((stripFired = IsHitCreated(distance, gap))) // simulate hit efficiency (add flag MpdTofUtils::IsSingle)
         {
            AddHit(MpdTofPoint::ClearGap(suid), hitPosition, hitPosError, pointIndex, tid, time, MpdTofUtils::IsSingle);
            nSingleHits++;
         }

         if (pQA) {
            if (stripFired) pQA->Point2HitSmearingTest(mcPosition, hitPosition);
            pQA->HitGapEfficiencyTest(stripFired, distance, gap);
            pQA->PositionInsideStripTest(strip->center, mcPosition, hitPosition);
            pQA->RotationToOriginTest(strip->fMatrix, mcPosition, hitPosition);
            pQA->CentralDetectorTest(suid, hitPosition);
         }

         if ((stripFired =
                 IsCrossHitCreated(distance, gap))) // simulate cross hit efficiency (add flag MpdTofUtils::IsDouble)
         {
            Int_t crossSuid =
               (side == LStrip::kRight) ? strip->neighboring[LStrip::kRight] : strip->neighboring[LStrip::kLeft];

            if (LStrip::kInvalid == crossSuid) continue; // edge strip on detector

            strip        = MpdTofGeoUtils::Instance()->FindStrip(crossSuid);
            centralStrip = MpdTofGeoUtils::Instance()->FindStrip(MpdTofPoint::SetCentralGap(crossSuid));
            if (!SmearingAlongStrip(strip->fMatrix, centralStrip->fMatrix, mcPosition, hitPosition, hitPosError)) {
               strip->Dump(" [MpdTofHitProducer::Exec] -E- Invalid Rotation matrix.");
               continue;
            }

            AddHit(MpdTofPoint::ClearGap(crossSuid), hitPosition, hitPosError, pointIndex, tid, time,
                   MpdTofUtils::IsDouble);
            nCrossHits++;
         }

         if (pQA) pQA->FillCrossHitEfficiency(stripFired, distance, gap, mcPosition, hitPosition);

      } // cycle by the TOF points

   } else //  exp. data used
   {
      // FIXME: now not realized
      // AddHit(Int_t detUID, const TVector3 &posHit, const TVector3 &posHitErr, Int_t expDigitIndex, Double_t time,
      // Int_t flag)
      assert(false);
   }

   MergeHitsOnStrip(); // save only the fastest hit in the strip

   int Nhits = CompressHits(); // remove blank slotes

   if (fUseMCData && pQA) {
      pQA->FillNPointsHits(aMcPoints->GetEntriesFast(), Nhits);
      pQA->PointDistanceTest(aMcPoints);
   }

   LOG(debug1) << "[MpdTofHitProducer::Exec] single hits= " << nSingleHits << ", cross hits= " << nCrossHits
               << ", final hits= " << Nhits;
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofHitProducer::Finish()
{
   if (pQA) pQA->Finish();
}
//--------------------------------------------------------------------------------------------------------------------------------------
void MpdTofHitProducer::SetSeed(UInt_t seed)
{
   pRandom->SetSeed(seed);
}
//------------------------------------------------------------------------------------------------------------------------

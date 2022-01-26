
//------------------------------------------------------------------------------------------------------------------------
/// \class MpdEtofHitProducer
///
/// \brief
/// \author Sergei Lobastov (LHE, JINR, Dubna)
//------------------------------------------------------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <valarray>
#include <assert.h>

#include "TClonesArray.h"
#include "TRandom.h"
#include "TVector3.h"
#include "TStopwatch.h"
#include "TMath.h"

#include "FairRunAna.h"
#include "FairBaseParSet.h"
#include "FairRuntimeDb.h"
#include "FairMCApplication.h"
#include "FairDetector.h"
#include "FairRootManager.h"
#include "MpdMCTrack.h"
#include "FairGeoVolume.h"

#include "MpdEtofGeoPar.h"
#include "MpdTofPoint.h"
#include "MpdTofHit.h"
#include "MpdTofHitProducerQA.h"
#include "MpdEtof.h"
#include "MpdEtofGeoUtils.h"

#include "MpdEtofHitProducer.h"

using namespace std;

ClassImp(MpdEtofHitProducer)
   //------------------------------------------------------------------------------------------------------------------------
   MpdEtofHitProducer::MpdEtofHitProducer(const char *name, Bool_t useMCdata, Int_t verbose, Bool_t test,
                                          const char *flnm)
   : MpdEtofHitProducerIdeal(name, useMCdata, verbose, test, true, flnm), fTimeSigma(0.100), fErrPhi(0.5),
     fErrR(1. / sqrt(12.)), pRandom(new TRandom2)
{
}
//------------------------------------------------------------------------------------------------------------------------
MpdEtofHitProducer::~MpdEtofHitProducer()
{
   delete pRandom;
}
//------------------------------------------------------------------------------------------------------------------------
InitStatus MpdEtofHitProducer::Init()
{
   //	if(fOnlyPrimary) cout<<" Only primary particles are processed!!! \n"; // FIXME NOT used now ADDD

   if (fUseMCData) {
      aMcPoints = (TClonesArray *)FairRootManager::Instance()->GetObject("ETOFPoint");
      aMcTracks = (TClonesArray *)FairRootManager::Instance()->GetObject("MCTrack");
      assert(aMcPoints);
      assert(aMcTracks);
   } else {
      aExpDigits = (TClonesArray *)FairRootManager::Instance()->GetObject("??????"); // FIXME: NOW unknown name
      assert(aExpDigits);
   }

   // Create and register output array
   aTofHits = new TClonesArray("MpdTofHit");
   FairRootManager::Instance()->Register("ETOFHit", "ETof", aTofHits, kTRUE);

   MpdEtofGeoUtils::Instance()->ParseTGeoManager(true); // forced

   MpdEtofGeoUtils::Instance()->FindNeighborStrips(0.8,   // 0.8 [cm] <--- thresh. distance between neighbor strips
                                                   true); // forced

   LOG(INFO) << "MpdTofHitProducer initialization finished succesfully.";

   return kSUCCESS;
}
//------------------------------------------------------------------------------------------------------------------------
Bool_t MpdEtofHitProducer::HitExist(Double_t val) // val - rasstojanie do kraja pad
{
   constexpr Double_t slope      = (0.98 - 0.95) / 0.2;
   Double_t           efficiency = (val > 0.2) ? 0.98 : (0.95 + slope * val);

   //-------------------------------------
   // 99% ---------
   //              \.
   //               \.
   //                \.
   // 95%             \.
   //  <-----------|--|
   //            0.2  0.
   //-------------------------------------

   if (pRandom->Rndm() < efficiency) return true;
   return false;
}
//------------------------------------------------------------------------------------------------------------------------
Bool_t MpdEtofHitProducer::DoubleHitExist(Double_t val) // val - rasstojanie do kraja pad
{
   constexpr Double_t slope      = (0.3 - 0.0) / 0.5;
   Double_t           efficiency = (val > 0.5) ? 0. : (0.3 - slope * val);

   //-------------------------------------
   // 30%               /
   //                  /
   //                 /
   //                /
   // 0%            /
   //  <-----------|----|
   //            0.5    0.
   //-------------------------------------

   if (efficiency == 0.) return false;

   if (pRandom->Rndm() < efficiency) return HitExist(val);
   return false;
}
//------------------------------------------------------------------------------------------------------------------------
void MpdEtofHitProducer::Exec(Option_t *option)
{
   /*	static const TVector3 XYZ_err(0.1, 0.1, 0.1); // FIXME: dummy now,  MUST BE ROTATED!!!!

      aTofHits->Clear();

      Int_t 		UID, trackID;
      TVector3 	pos, XYZ_smeared;
         int		nSingleHits = 0, nDoubleHits = 0;

      if(fUseMCData)
      {

         for(Int_t pointIndex = 0, nTofPoint = aMcPoints->GetEntriesFast(); pointIndex < nTofPoint; pointIndex++ )  //
      cycle by TOF points
         {
            MpdTofPoint *pPoint = (MpdTofPoint*) aMcPoints->UncheckedAt(pointIndex);

            if(fVerbose > 2) pPoint->Print("");

            trackID = pPoint->GetTrackID();
            UID	= pPoint->GetDetectorID();
            Double_t time = pRandom->Gaus(pPoint->GetTime(), fTimeSigma); // 100 ps
            pPoint->Position(pos);

            // Calc Rerr, PhiErr onto MRF
            const LStrip *pStrip = MpdEtofGeoUtils::Instance()->FindStrip(UID);
            double centerPhi = pStrip->center.Phi(); // [rad]

            TVector3 posRotated(pos);
            posRotated.RotateZ(-centerPhi);

            XYZ_smeared.SetXYZ(posRotated.X(), pRandom->Gaus(posRotated.Y(), fErrPhi), posRotated.Z());
            XYZ_smeared.RotateZ(centerPhi);

            LStrip::Side_t side;
            Double_t distance = pStrip->MinDistanceToEdge(&pos, side); // [cm]

            bool passed;
            if(passed = HitExist(distance)) // check efficiency
            {
               AddHit(UID, XYZ_smeared, XYZ_err, pointIndex, trackID, time, MpdTofUtils::IsSingle);
               nSingleHits++;

               if(pHitProducerQA) pHitProducerQA->FillSingleHitPosition(pos, XYZ_smeared);
            }

            if(pQA) pQA->GetSingleHitEfficiency()->Fill(passed, distance);

               if(passed = DoubleHitExist(distance)) // check cross hit
               {
                  Int_t CrossUID = (side == LStrip::kRight) ? pStrip->neighboring[LStrip::kRight] :
      pStrip->neighboring[LStrip::kLeft];

               if(LStrip::kInvalid  == CrossUID) continue; // last strip on module

               // Calc Rerr, PhiErr onto MRF
                  pStrip = MpdEtofGeoUtils::Instance()->FindStrip(CrossUID);
               //	XYZ_smeared.SetXYZ(pRandom->Gaus(posRotated.X(), fErrX), posRotated.Y(),  pStrip->center.Z());
               //	XYZ_smeared.RotateZ(perpPhi);

                  AddHit(CrossUID, XYZ_smeared, XYZ_err, pointIndex, trackID, time, MpdTofUtils::IsDouble);
                  nDoubleHits++;

                  if(pHitProducerQA) pHitProducerQA->FillDoubleHitPosition(pos, XYZ_smeared);
               }

               if(pQA) pQA->GetDoubleHitEfficiency()->Fill(passed, distance);

            }	// cycle by the TOF points
      }
      else	//  exp. data used
      {
         // FIXME: now not realized
         //AddHit(Int_t detUID, const TVector3 &posHit, const TVector3 &posHitErr, Int_t expDigitIndex, Double_t time,
      Int_t flag) assert(false);
      }

      MergeHitsOnStrip(); // save only the fastest hit in the strip

      int nFinally = CompressHits(); // remove blank slotes

           cout<<" -I- [MpdEtofHitProducer::Exec] single hits= "<<nSingleHits<<", double hits= "<<nDoubleHits<<", final
      hits= "<<nFinally<<endl;
   */
}
//------------------------------------------------------------------------------------------------------------------------
/*void 		MpdEtofHitProducer::Finish()
{
   if(pQA) pQA->Finish();
}*/
//------------------------------------------------------------------------------------------------------------------------
void MpdEtofHitProducer::SetSeed(UInt_t seed)
{
   pRandom->SetSeed(seed);
}
//------------------------------------------------------------------------------------------------------------------------
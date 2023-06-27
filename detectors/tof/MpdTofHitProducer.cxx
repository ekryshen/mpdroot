//------------------------------------------------------------------------------------------------------------------------
/// \class MpdTofHitProducer
///
/// \brief
/// \author Sergei Lobastov (LHE, JINR, Dubna)
/// \author Victor Baryshnikov
//------------------------------------------------------------------------------------------------------------------------
#include <assert.h>

#include <TMath.h>
#include <TFile.h>
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
#include "MpdTofDigit.h"
#include "MpdTofGeoUtils.h"
#include "MpdTofHitProducerQA.h"

#include "MpdTofHitProducer.h"

using namespace std;

ClassImp(MpdTofHitProducer);
//------------------------------------------------------------------------------------------------------------------------
MpdTofHitProducer::MpdTofHitProducer(const char *name, Bool_t useMCdata, Int_t verbose, const char *CLCflnm, Bool_t test, const char *flnm,
                                     double stripLength)
   : MpdTofHitProducerIdeal(name, useMCdata, verbose, test, true, flnm, false), fStripLength(stripLength)
{
   pRandom = new TRandom2;
   if(CLCflnm) fCLCflnm = CLCflnm;
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

	if(fUseMCData)
	{
		MpdTofGeoUtils::Instance()->FindNeighborStrips(0.8, pQA, // 0.8 [cm] <--- thresh. distance between neighbor strips
                                                  true);    // forced
	}
	else
	{
		// Load cable length time correction for strips.
		if(LoadCLcorrections(fCLCflnm) != true) return kFATAL;
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
// Comparator of MpdTofDigits by strip guid (MpdTofDigit::operator<) 
struct lessPtrByGuid
{
	bool operator()(const MpdTofDigit*  lhs, const MpdTofDigit*  rhs) const 
	{
		return (*lhs) < (*rhs); 
	}
};
//------------------------------------------------------------------------------------------------------------------------
void MpdTofHitProducer::Exec(Option_t *option)
{
   aTofHits->Clear();

	// position errors at strip local reference frame
	const TVector3 errorLRF(fErrX, 0.,fErrZ);


    static const Double_t fMaxDelta = (64 * 0.5 + 10) * 0.062; // + 100 mm on the strip edge (for LR corrections) // REMOVE magic numbers

   TVector3 mcPosition, hitPosition;
   int      nSingleHits = 0, nCrossHits = 0;

   auto SmearingAlongStrip = [this](const TGeoCombiTrans &gap, const TGeoCombiTrans &centralGap, const TVector3 &point, TVector3 &hit) {
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

      return true;
   };

   if (fUseMCData) {
      for (Int_t pointIndex = 0, nTofPoint = aMcPoints->GetEntriesFast(); pointIndex < nTofPoint;
           pointIndex++) // cycle by TOF points
      {
         auto pPoint = (MpdTofPoint *)aMcPoints->UncheckedAt(pointIndex);

         if (fVerbose > 2) pPoint->Print("");

         Int_t tid  = pPoint->GetTrackID();
         Int_t suid = pPoint->GetDetectorID();
         Int_t gap  = pPoint->GetGap();

         Double_t time = pRandom->Gaus(pPoint->GetTime(), fTimeSigma); // default 100 ps
         pPoint->Position(mcPosition);

         auto strip        = MpdTofGeoUtils::Instance()->FindStrip(suid);
         auto centralStrip = MpdTofGeoUtils::Instance()->FindStrip(MpdTofPoint::SetCentralGap(suid));

         if (!SmearingAlongStrip(strip->fMatrix, centralStrip->fMatrix, mcPosition, hitPosition)) {
            strip->Dump(" [MpdTofHitProducer::Exec] -E- Invalid Rotation matrix.");
            continue;
         }

         LStrip::Side_t side;
         Double_t       distance = strip->MinDistanceToEdge(&mcPosition, side); // [cm]

         bool stripFired;
         if ((stripFired = IsHitCreated(distance, gap))) // simulate hit efficiency (add flag MpdTofUtils::IsSingle)
         {
		const TVector3 hitPosError = strip->RotateErrors(errorLRF);
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

            if (!SmearingAlongStrip(strip->fMatrix, centralStrip->fMatrix, mcPosition, hitPosition)) {
               strip->Dump(" [MpdTofHitProducer::Exec] -E- Invalid Rotation matrix.");
               continue;
            }

		const TVector3 hitPosError = strip->RotateErrors(errorLRF);
            AddHit(MpdTofPoint::ClearGap(crossSuid), hitPosition, hitPosError, pointIndex, tid, time,
                   MpdTofUtils::IsDouble);
            nCrossHits++;
         }

         if (pQA) pQA->FillCrossHitEfficiency(stripFired, distance, gap, mcPosition, hitPosition);

      } // cycle by the TOF points

   } 
   else //  exp. data used
   {
	// Sort MpdTofDigit pointers by suid
	multiset<MpdTofDigit*, lessPtrByGuid> sDigits; 	
	for(size_t index = 0, size = aExpDigits->GetEntriesFast(); index < size; index++)   						
	{
		auto digit = (MpdTofDigit*) aExpDigits->UncheckedAt(index);
		if(fVerbose > 3) digit->print("\n MpdTofDigit: ");

		sDigits.insert(digit);
	}

	for(auto it1 = sDigits.cbegin(), itEnd = sDigits.cend(); it1 != itEnd; it1 = sDigits.upper_bound(*it1)) // cycle by unique suid 
	{
		MpdTofDigit	*digit1 = nullptr, *digit2 = nullptr; // fastest digits from both edges of the strip

		if(sDigits.count(*it1) == 2) // ideal case
		{
			digit1 = *it1;

			auto it2 = std::next(it1, 1);
		 
			if(digit1->GetSide() != (*it2)->GetSide()) // set digit to different edge
			{
				digit2 = *it2;
			}
		}
		else if(sDigits.count(*it1) > 2) // additional noise digits case, select fastest for both sides
		{
			const auto range = std::equal_range(sDigits.begin(), sDigits.end(), *it1);
  			for(auto iter = range.first; iter != range.second; ++iter) // cycle by same strip digits
			{
				auto entry = *iter;
				auto digit = (entry->GetSide() == 1) ? digit1 : digit2;

				if(digit == nullptr) 				digit = entry; // first compare	
				else if(entry->GetTime() < digit->GetTime()) 	digit = entry; // change to fastest
			}
		}
		
		if(digit1 != nullptr && digit2 != nullptr) // both edges have digit
		{
			// set digit1 to left side 
			if(digit1->GetSide() == 1) std::swap(digit1, digit2); // digit1 = left side

assert(digit1->GetSector()==digit2->GetSector());
assert(digit1->GetPlane()==digit2->GetPlane());
assert(digit1->GetStrip()==digit2->GetStrip());
assert(digit1->GetSide() + digit2->GetSide() == 1);

			if(digit1->GetAmplitude() + digit2->GetAmplitude() == 0.) continue;

			Int_t sector = digit1->GetSector(), detector = digit1->GetPlane()+1, strip = digit1->GetStrip()+1;
	  		Int_t suid = MpdTofPoint::GetSuid24(sector, detector, 2, strip); // 2 = central gap

			auto pStrip = MpdTofGeoUtils::Instance()->FindStrip(suid);
			if(pStrip)
			{
				Double_t dL = (digit1->GetTime() - digit2->GetTime() - pStrip->timeCL) * fSlope; // [cm]

				TVector3 wireEdge1 = (pStrip->A + pStrip->D) *0.5, wireEdge2 = (pStrip->B + pStrip->C) *0.5;
				TVector3 hitPosition = pStrip->center + dL * (wireEdge2 - wireEdge1).Unit(); // strip center + dL shift in wire direction
		
				TVector3 hitPosError = pStrip->RotateErrors(errorLRF); 

				Double_t time = (digit1->GetTime() + digit2->GetTime()) * 0.5;
				Double_t width = digit1->GetAmplitude() + digit2->GetAmplitude();				

 				if(fabs((digit1->GetTime() - digit2->GetTime() - pStrip->timeCL)* 0.5) > fMaxDelta)  continue; // [cm] far out of strip

				auto hit = new ((*aTofHits)[aTofHits->GetEntriesFast()]) MpdTofHit(suid, hitPosition, hitPosError, -1, time, 0, width);
				if(fVerbose > 3) hit->Print("\n ");
			}
		}
	}
   }

   MergeHitsOnStrip(); // save only the fastest hit in the strip

   int Nhits = CompressHits(); // remove blank slotes

   if (fUseMCData && pQA) {
      pQA->FillNPointsHits(aMcPoints->GetEntriesFast(), Nhits);
      pQA->PointDistanceTest(aMcPoints);
   }

   if (fUseMCData) LOG(debug1) << " single hits= " << nSingleHits << ", cross hits= " << nCrossHits<< ", final hits= " << Nhits;
   else		   LOG(debug1) << " hits= " << Nhits;
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofHitProducer::Finish()
{
   if (pQA) pQA->Finish();
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofHitProducer::SetSeed(UInt_t seed)
{
   pRandom->SetSeed(seed);
}
//------------------------------------------------------------------------------------------------------------------------
Bool_t		MpdTofHitProducer::LoadCLcorrections(const char* flnm)
{
	size_t NSector = MpdTofGeoUtils::Instance()->GetNSectors(); 
        size_t NDetector = (MpdTofGeoUtils::Instance()->GetNDetectors())/NSector; 

        if (NSector <= 0 || NDetector <= 0) 
	{
                LOG(fatal)<<" No TOF detectors in geometry file for the current run! Task will be deactivated";              
                return false;
        }   

	TFile file(flnm, "READ");
	if(! file.IsOpen())
	{
		LOG(fatal)<<" can't open CL correction file="<<flnm;
		return false;
	}

	auto tree = (TTree*) file.Get("data");
	Double_t timeCL;
	tree->SetBranchAddress("correction", &timeCL);

	size_t loadedNmb = 0;
	for(size_t sector = 1; sector <= NSector; sector++)
	{
		for(size_t detector = 1; detector <= NDetector; detector++)
		{
			for(size_t strip = 1; strip <= 24; strip++)
			{
				if(tree->GetEntry((sector-1)*240 + (detector-1)*24 + (strip-1)) > 0) // return 0, if entry don't exist; return -1, if I/O error
				{
  					Int_t suid = MpdTofPoint::GetSuid24(sector, detector, 2, strip); // 2 = central gap
					auto pStrip = MpdTofGeoUtils::Instance()->FindStrip(suid);
					if(pStrip)
					{
						const_cast<LStrip*>(pStrip)->timeCL = timeCL;
						loadedNmb++;
					}
				}
			}
		}
	}  

	file.Close();

	LOG(debug1)<<" Number of loaded CL time corrections: "<<loadedNmb;

return true;
}
//------------------------------------------------------------------------------------------------------------------------


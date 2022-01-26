
//------------------------------------------------------------------------------------------------------------------------
/// \class MpdEtofHitProducerIdeal
///
/// \brief
/// \author Sergei Lobastov (LHE, JINR, Dubna)
//------------------------------------------------------------------------------------------------------------------------

#include <assert.h>
#include <iostream>

#include "TClonesArray.h"

#include "FairRootManager.h"
#include "FairLogger.h"

#include "MpdTofUtils.h"
#include "MpdTofHit.h"
#include "MpdTofPoint.h"
#include "MpdTofHitProducerQA.h"

#include "MpdEtofHitProducerIdeal.h"

ClassImp(MpdEtofHitProducerIdeal) using namespace std;
//------------------------------------------------------------------------------------------------------------------------
MpdEtofHitProducerIdeal::MpdEtofHitProducerIdeal(const char *name, Bool_t useMCdata, Int_t verbose, Bool_t test,
                                                 Bool_t merge, const char *flnm)
   : MpdTofHitProducerIdeal(name, useMCdata, verbose, test, merge, flnm, true)
{
}
//------------------------------------------------------------------------------------------------------------------------
MpdEtofHitProducerIdeal::~MpdEtofHitProducerIdeal() {}
//------------------------------------------------------------------------------------------------------------------------
InitStatus MpdEtofHitProducerIdeal::Init()
{
   /*	FairLogger::GetLogger()->Info(MESSAGE_ORIGIN, "[MpdEtofHitProducerIdeal::Init] Begin initialization.");

      if(fUseMCData)
      {
            aMcPoints = (TClonesArray*) FairRootManager::Instance()->GetObject("ETOFPoint");
            aMcTracks = (TClonesArray*) FairRootManager::Instance()->GetObject("MCTrack");
   assert(aMcPoints);
   assert(aMcTracks);
      }
      else
      {
            aExpDigits = (TClonesArray*) FairRootManager::Instance()->GetObject("??????");// FIXME: NOW unknown name
   assert(aExpDigits);
      }

           // Create and register output array
           aTofHits = new TClonesArray("MpdEtofHit");
           FairRootManager::Instance()->Register("ETOFHit", "ETof", aTofHits, kTRUE);

           FairLogger::GetLogger()->Info(MESSAGE_ORIGIN, "[MpdEtofHitProducerIdeal::Init] Initialization finished
   succesfully.");
   */
   return kSUCCESS;
}
//------------------------------------------------------------------------------------------------------------------------
void MpdEtofHitProducerIdeal::Exec(Option_t *opt)
{
   /*
      static const TVector3 XYZ_err(0., 0., 0.); // FIXME:

      aTofHits->Clear();

      Int_t 		nSingleHits = 0;
      TVector3 	pos;

      if(fUseMCData)
      {
         for(Int_t pointIndex = 0, nTofPoint = aMcPoints->GetEntriesFast(); pointIndex < nTofPoint; pointIndex++ )  //
      cycle by TOF points
         {
            MpdTofPoint *pPoint = (MpdTofPoint*) aMcPoints->UncheckedAt(pointIndex);
            pPoint->Position(pos);

            AddHit(pPoint->GetDetectorID(), pos, XYZ_err, pointIndex, pPoint->GetTrackID(), pPoint->GetTime(),
      MpdTofUtils::IsSingle); nSingleHits++;
         }
      }
      else //  exp. data used
      {
         // FIXME: now not realized
         //AddHit(Int_t detUID, const TVector3 &posHit, const TVector3 &posHitErr, Int_t expDigitIndex, Double_t time,
      Int_t flag) assert(false);
      }

      int nFinally;
      if(fDoMergeHits)
      {
         MergeHitsOnStrip(); 		// save only the fastest hit in the strip
         nFinally = CompressHits(); 	// remove blank slotes
      }
      else
         nFinally = aTofHits->GetEntriesFast();

           cout<<" -I- [MpdEtofHitProducerIdeal::Exec] single hits= "<<nSingleHits<<", final hits= "<<nFinally<<endl;
   */
}
//------------------------------------------------------------------------------------------------------------------------
/*void 			MpdEtofHitProducerIdeal::Finish()
{
   if(pHitProducerQA) pHitProducerQA->Finish();
}*/
//------------------------------------------------------------------------------------------------------------------------

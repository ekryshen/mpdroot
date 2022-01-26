//------------------------------------------------------------------------------------------------------------------------
#ifndef __MPD_ETOF_HIT_PRODUCER_IDEAL_H
#define __MPD_ETOF_HIT_PRODUCER_IDEAL_H 1

//------------------------------------------------------------------------------------------------------------------------
/// \class MpdEtofHitProducerIdeal
///
/// \brief
/// \author Sergei Lobastov (LHE, JINR, Dubna)
//------------------------------------------------------------------------------------------------------------------------

#include <TVector3.h>
#include <TList.h>

#include "MpdTofHitProducerIdeal.h"

//------------------------------------------------------------------------------------------------------------------------
class MpdEtofHitProducerIdeal : public MpdTofHitProducerIdeal {

public:
   MpdEtofHitProducerIdeal(const char *name = "ETOF Ideal HitProducer", Bool_t useMCdata = true, Int_t verbose = 1,
                           Bool_t DoTest = false, Bool_t DoMergeHits = false,
                           const char *flnm = "QA.MpdEtofHitProducerIdeal.root");
   virtual ~MpdEtofHitProducerIdeal();

   virtual InitStatus Init();
   virtual void       Exec(Option_t *option);
   //	virtual void		Finish();

   ClassDef(MpdEtofHitProducerIdeal, 3);
};
//------------------------------------------------------------------------------------------------------------------------
#endif

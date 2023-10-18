//------------------------------------------------------------------------------------------------------------------------
/// \class MpdFwdHitProducer
/// \brief Ideal hit producer
/// \author Evgeny Kryshen (PNPI, Gatchina)
//------------------------------------------------------------------------------------------------------------------------
#ifndef __HH_MPDFWDHITPRODUCER_H
#define __HH_MPDFWDHITPRODUCER_H

#include "FairTask.h"
class TClonesArray;

class MpdFwdHitProducer : public FairTask {
public:
   MpdFwdHitProducer(Double_t sigma = 0.01, Double_t timeRes = 0.05);
   virtual ~MpdFwdHitProducer();
   virtual InitStatus Init();
   virtual void       Exec(Option_t *opt);
   virtual void       Finish();
private:
   TClonesArray *fMcPoints;  //!
   TClonesArray *fFwdHits;   //!
   Double_t fSigma;          // space resolution in cm 
   Double_t fTimeRes;        // time resolution in ns

   ClassDef(MpdFwdHitProducer, 1)
};

#endif // #ifndef __HH_MPDFWDHITPRODUCER_H

//------------------------------------------------------------------------------------------------------------------------
/// \class MpdFwdHitProducer
///
/// \brief
/// \author Evgeny Kryshen (PNPI, Gatchina)
//------------------------------------------------------------------------------------------------------------------------
#ifndef __HH_MPDFWDHITPRODUCER_H
#define __HH_MPDFWDHITPRODUCER_H

#include "FairTask.h"
class TClonesArray;

class MpdFwdHitProducer : public FairTask {
public:
   MpdFwdHitProducer();
   virtual ~MpdFwdHitProducer();
   virtual InitStatus Init();
   virtual void       Exec(Option_t *opt);
   virtual void       Finish();
private:
   TClonesArray *fMcPoints;  //!
   TClonesArray *fFwdHits;   //!

   ClassDef(MpdFwdHitProducer, 1)
};

#endif // #ifndef __HH_MPDFWDHITPRODUCER_H

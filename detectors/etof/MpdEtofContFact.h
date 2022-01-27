//------------------------------------------------------------------------------------------------------------------------
#ifndef __MPD_ETOF_CONT_FACT_H
#define __MPD_ETOF_CONT_FACT_H 1

#include "FairContFact.h"

class FairContainer;
//------------------------------------------------------------------------------------------------------------------------
class MpdEtofContFact : public FairContFact {
private:
   void setAllContainers();

public:
   MpdEtofContFact();
   ~MpdEtofContFact(){};

   FairParSet *createContainer(FairContainer *);

   ClassDef(MpdEtofContFact, 0); // Factory for all TOF parameter containers
};
//------------------------------------------------------------------------------------------------------------------------
#endif

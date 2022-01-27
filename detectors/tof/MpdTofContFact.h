//------------------------------------------------------------------------------------------------------------------------
#ifndef __MPD_TOF_CONT_FACT_H
#define __MPD_TOF_CONT_FACT_H

#include "FairContFact.h"

class FairContainer;
//------------------------------------------------------------------------------------------------------------------------
class MpdTofContFact : public FairContFact {
private:
   void setAllContainers();

public:
   MpdTofContFact();
   ~MpdTofContFact(){};

   FairParSet *createContainer(FairContainer *);

   ClassDef(MpdTofContFact, 0) // Factory for all TOF parameter containers
};
//------------------------------------------------------------------------------------------------------------------------
#endif // #ifndef __MPD_TOF_CONT_FACT_H

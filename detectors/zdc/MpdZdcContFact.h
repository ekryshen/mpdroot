// -------------------------------------------------------------------------
// -----                     MpdzdcContFact header file                -----
// -------------------------------------------------------------------------

#ifndef MPDZDCCONTFACT_H
#define MPDZDCCONTFACT_H

#include "FairContFact.h"

class FairContainer;

class MpdZdcContFact : public FairContFact {
private:
   void setAllContainers();

public:
   MpdZdcContFact();
   ~MpdZdcContFact() {}
   FairParSet *createContainer(FairContainer *);

   ClassDef(MpdZdcContFact, 0) // Factory for all HYP parameter containers
};

#endif /* !MPDZDCCONTFACT_H */

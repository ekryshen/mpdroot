#ifndef MPDSTSCONTFACT_H
#define MPDSTSCONTFACT_H

#include "FairContFact.h"

class FairContainer;

class MpdStsContFact : public FairContFact {
private:
   void setAllContainers();

public:
   MpdStsContFact();
   ~MpdStsContFact() {}
   FairParSet *createContainer(FairContainer *);
   ClassDef(MpdStsContFact, 0); // Factory for all STS parameter containers
};

#endif /* !MPDSTSCONTFACT_H */

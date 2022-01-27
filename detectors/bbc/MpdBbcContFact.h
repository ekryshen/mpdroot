#ifndef MPDBBCCONTFACT_H
#define MPDBBCCONTFACT_H

#include "FairContFact.h"

class FairContainer;

class MpdBbcContFact : public FairContFact {
private:
   void setAllContainers();

public:
   MpdBbcContFact();
   ~MpdBbcContFact() {}
   FairParSet *createContainer(FairContainer *);
   ClassDef(MpdBbcContFact, 0); // Factory for all BBC parameter containers
};

#endif /* !MPDBBCCONTFACT_H */

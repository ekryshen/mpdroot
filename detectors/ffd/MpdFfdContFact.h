#ifndef MPDFFDCONTFACT_H
#define MPDFFDCONTFACT_H

#include "FairContFact.h"

class FairContainer;

class MpdFfdContFact : public FairContFact {
private:
   void setAllContainers();

public:
   MpdFfdContFact();
   ~MpdFfdContFact() {}
   FairParSet *createContainer(FairContainer *);
   ClassDef(MpdFfdContFact, 0) // Factory for all FFD parameter containers
};

#endif /* !MPDFFDCONTFACT_H */

#ifndef MPDEMCDETCONTFACT_H
#define MPDEMCDETCONTFACT_H

#include "FairContFact.h"

class FairContainer;

class MpdEmcContFact : public FairContFact {
 private:
  void setAllContainers();
 public:
  MpdEmcContFact();
  ~MpdEmcContFact() {}
  FairParSet* createContainer(FairContainer*);
  ClassDef( MpdEmcContFact,0) // Factory for all MyDet parameter containers
};

#endif  /* !MPDEMCDETCONTFACT_H */

#ifndef TPCCONTFACT_H
#define TPCCONTFACT_H

#include "FairContFact.h"

class FairContainer;
// class FairParIo;

class TpcContFact : public FairContFact {
private:
   void setAllContainers();

public:
   TpcContFact();
   ~TpcContFact() {}
   FairParSet *createContainer(FairContainer *);
   // void  activateParIo(FairParIo* io);
   ClassDef(TpcContFact, 1); // Factory for all TPC parameter containers
};

#endif

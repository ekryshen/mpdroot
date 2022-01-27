//------------------------------------------------------------------------------------------------------------------------
#include <iostream>
#include <iomanip>

#include "MpdEtofContFact.h"
#include "FairRuntimeDb.h"
#include "MpdEtofGeoPar.h"

using namespace std;

static MpdEtofContFact gMpdEtofContFact;
//------------------------------------------------------------------------------------------------------------------------
MpdEtofContFact::MpdEtofContFact()
{
   fName  = "MpdEtofContFact";
   fTitle = "Factory for parameter containers in libEtof";
   setAllContainers();
   FairRuntimeDb::instance()->addContFactory(this);
}
//------------------------------------------------------------------------------------------------------------------------
void MpdEtofContFact::setAllContainers()
{
   FairContainer *p = new FairContainer("MpdEtofGeoPar", "Etof Geometry Parameters", "TestDefaultContext");
   p->addContext("TestNonDefaultContext");

   containers->Add(p);
}
//------------------------------------------------------------------------------------------------------------------------
FairParSet *MpdEtofContFact::createContainer(FairContainer *c)
{

   const char *name = c->GetName();
   cout << "-I- container name: " << name << endl;
   FairParSet *p = NULL;
   if (strcmp(name, "MpdEtofGeoPar") == 0)
      p = new MpdEtofGeoPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());

   return p;
}
//------------------------------------------------------------------------------------------------------------------------
ClassImp(MpdEtofContFact);

//------------------------------------------------------------------------------------------------------------------------
#include <iostream>
#include <iomanip>

#include "FairRuntimeDb.h"
#include "MpdTofGeoPar.h"

#include "MpdTofContFact.h"

using namespace std;

static MpdTofContFact gMpdTofContFact;
//------------------------------------------------------------------------------------------------------------------------
MpdTofContFact::MpdTofContFact()
{
   fName  = "MpdTofContFact";
   fTitle = "Factory for parameter containers in libTof";
   setAllContainers();
   FairRuntimeDb::instance()->addContFactory(this);
}
//------------------------------------------------------------------------------------------------------------------------
void MpdTofContFact::setAllContainers()
{
   /** Creates the Container objects with all accepted contexts and adds them to
    *  the list of containers for the Tof library.*/

   FairContainer *p = new FairContainer("MpdTofGeoPar", "Tof Geometry Parameters", "TestDefaultContext");
   p->addContext("TestNonDefaultContext");

   containers->Add(p);
}
//------------------------------------------------------------------------------------------------------------------------
FairParSet *MpdTofContFact::createContainer(FairContainer *c)
{
   /** Calls the constructor of the corresponding parameter container.
    * For an actual context, which is not an empty string and not the default context
    * of this container, the name is concatinated with the context. */

   const char *name = c->GetName();
   cout << "-I- container name: " << name << endl;
   FairParSet *p = NULL;

   if (strcmp(name, "MpdTofGeoPar") == 0)
      p = new MpdTofGeoPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());

   return p;
}
//------------------------------------------------------------------------------------------------------------------------
ClassImp(MpdTofContFact);

/////////////////////////////////////////////////////////////
//
//  MpdStsContFact
//
//  Factory for the parameter containers in libSts
//
/////////////////////////////////////////////////////////////

using namespace std;
#include "MpdStsContFact.h"
#include "FairRuntimeDb.h"
#include "MpdStsGeoPar.h"
#include <iostream>
#include <iomanip>

ClassImp(MpdStsContFact)

static MpdStsContFact gMpdStsContFact;

MpdStsContFact::MpdStsContFact()
{
   // Constructor (called when the library is loaded)
   fName  = "MpdStsContFact";
   fTitle = "Factory for parameter containers in libSts";
   setAllContainers();
   FairRuntimeDb::instance()->addContFactory(this);
}

void MpdStsContFact::setAllContainers()
{
   /** Creates the Container objects with all accepted contexts and adds them to
    *  the list of containers for the Tof library.*/

   FairContainer *p = new FairContainer("MpdStsGeoPar", "Sts Geometry Parameters", "TestDefaultContext");
   p->addContext("TestNonDefaultContext");

   containers->Add(p);
}

FairParSet *MpdStsContFact::createContainer(FairContainer *c)
{
   /** Calls the constructor of the corresponding parameter container.
    * For an actual context, which is not an empty string and not the default context
    * of this container, the name is concatinated with the context. */
   const char *name = c->GetName();
   cout << " -I container name " << name << endl;
   FairParSet *p = 0;
   if (strcmp(name, "MpdStsGeoPar") == 0) {
      p = new MpdStsGeoPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
   }
   return p;
}

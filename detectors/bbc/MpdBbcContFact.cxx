/////////////////////////////////////////////////////////////
//
//  MpdBbcContFact
//
//  Factory for the parameter containers in libBbc
//
/////////////////////////////////////////////////////////////

using namespace std;
#include "MpdBbcContFact.h"
#include "FairRuntimeDb.h"
#include "MpdBbcGeoPar.h"
#include <iostream>
#include <iomanip>

ClassImp(MpdBbcContFact);

static MpdBbcContFact gMpdBbcContFact;

MpdBbcContFact::MpdBbcContFact()
{
   // Constructor (called when the library is loaded)
   fName  = "MpdBbcContFact";
   fTitle = "Factory for parameter containers in libBbc";
   setAllContainers();
   FairRuntimeDb::instance()->addContFactory(this);
}

void MpdBbcContFact::setAllContainers()
{
   /** Creates the Container objects with all accepted contexts and adds them to
    *  the list of containers for the Tof library.*/

   FairContainer *p = new FairContainer("MpdBbcGeoPar", "Bbc Geometry Parameters", "TestDefaultContext");
   p->addContext("TestNonDefaultContext");

   containers->Add(p);
}

FairParSet *MpdBbcContFact::createContainer(FairContainer *c)
{
   /** Calls the constructor of the corresponding parameter container.
    * For an actual context, which is not an empty string and not the default context
    * of this container, the name is concatinated with the context. */
   const char *name = c->GetName();
   cout << " -I container name " << name << endl;
   FairParSet *p = 0;
   if (strcmp(name, "MpdBbcGeoPar") == 0) {
      p = new MpdBbcGeoPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
   }
   return p;
}

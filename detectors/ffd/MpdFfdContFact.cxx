/////////////////////////////////////////////////////////////
//
//  MpdFfdContFact
//
//  Factory for the parameter containers in libFfd
//
/////////////////////////////////////////////////////////////

using namespace std;
#include "MpdFfdContFact.h"
#include "FairRuntimeDb.h"
#include "MpdFfdGeoPar.h"
#include <iostream>
#include <iomanip>

ClassImp(MpdFfdContFact)

   static MpdFfdContFact gMpdFfdContFact;

MpdFfdContFact::MpdFfdContFact()
{
   // Constructor (called when the library is loaded)
   fName  = "MpdFfdContFact";
   fTitle = "Factory for parameter containers in libFfd";
   setAllContainers();
   FairRuntimeDb::instance()->addContFactory(this);
}

void MpdFfdContFact::setAllContainers()
{
   /** Creates the Container objects with all accepted contexts and adds them to
    *  the list of containers for the Tof library.*/

   FairContainer *p = new FairContainer("MpdFfdGeoPar", "Ffd Geometry Parameters", "TestDefaultContext");
   p->addContext("TestNonDefaultContext");

   containers->Add(p);
}

FairParSet *MpdFfdContFact::createContainer(FairContainer *c)
{
   /** Calls the constructor of the corresponding parameter container.
    * For an actual context, which is not an empty string and not the default context
    * of this container, the name is concatinated with the context. */
   const char *name = c->GetName();
   cout << " -I container name " << name << endl;
   FairParSet *p = 0;
   if (strcmp(name, "MpdFfdGeoPar") == 0) {
      p = new MpdFfdGeoPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
   }
   return p;
}

/////////////////////////////////////////////////////////////
//
//  TpcContFact
//
//  Factory for the parameter containers in libTst
//
/////////////////////////////////////////////////////////////

using namespace std;
#include "TpcContFact.h"
#include "FairRuntimeDb.h"
#include "TpcGeoPar.h"
#include "FairParRootFileIo.h"
#include "FairParAsciiFileIo.h"
//#include "FairParIo.h"
//#include "TpcParRootFileIo.h"
//#include "TpcParAsciiFileIo.h"
#include <iostream>
#include <iomanip>

ClassImp(TpcContFact);

static TpcContFact gTpcContFact;

TpcContFact::TpcContFact()
{
   // Constructor (called when the library is loaded)
   fName  = "TpcContFact";
   fTitle = "Factory for parameter containers in libTpc";
   setAllContainers();
   FairRuntimeDb::instance()->addContFactory(this);
}

void TpcContFact::setAllContainers()
{
   /** Creates the Container objects with all accepted contexts and adds them to
    *  the list of containers for the TST library.*/

   FairContainer *p = new FairContainer("TpcGeoPar", "Tpc Geometry Parameters", "TpcDefaultContext");
   p->addContext("TpcNonDefaultContext");

   containers->Add(p);
}

FairParSet *TpcContFact::createContainer(FairContainer *c)
{
   /** Calls the constructor of the corresponding parameter container.
    * For an actual context, which is not an empty string and not the default context
    * of this container, the name is concatinated with the context. */
   const char *name = c->GetName();
   FairParSet *p    = NULL;
   if (strcmp(name, "TpcGeoPar") == 0) {
      p = new TpcGeoPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
   }
   return p;
}

/*void  TpcContFact::activateParIo(FairParIo* io) {
  // activates the input/output class for the parameters
  // needed by the Tst
  if (strcmp(io->IsA()->GetName(),"FairParRootFileIo")==0) {
    TpcParRootFileIo* p=new TpcParRootFileIo(((FairParRootFileIo*)io)->getParRootFile());
    io->setDetParIo(p);
  }
  if (strcmp(io->IsA()->GetName(),"FairParAsciiFileIo")==0) {
    TpcParAsciiFileIo* p=new TpcParAsciiFileIo(((FairParAsciiFileIo*)io)->getFile());
    io->setDetParIo(p);
  }
}*/

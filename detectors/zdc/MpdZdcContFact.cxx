/*************************************************************************************
 *
 *         Class MpdZdcContFact
 *
 *  Adopted for MPD by:   Elena Litvinenko
 *  e-mail:   litvin@nf.jinr.ru
 *  Version:  8-Apr-2008
 *
 ************************************************************************************/

using namespace std;
#include "MpdZdcContFact.h"
#include "FairRuntimeDb.h"
#include "FairParRootFileIo.h"
#include "FairParAsciiFileIo.h"
#include "MpdZdcGeoPar.h"
#include <iostream>
#include <iomanip>

ClassImp(MpdZdcContFact);

static MpdZdcContFact gMpdZdcContFact;

MpdZdcContFact::MpdZdcContFact()
{
   // Constructor (called when the library is loaded)
   fName  = "MpdZdcContFact";
   fTitle = "Factory for parameter containers in libZdc";
   setAllContainers();
   FairRuntimeDb::instance()->addContFactory(this);
}

void MpdZdcContFact::setAllContainers()
{
   /** Creates the Container objects with all accepted contexts and adds them to
    *  the list of containers for the ZDC library.*/

   FairContainer *p = new FairContainer("MpdZdcGeoPar", "Zdc Geometry Parameters", "ZdcDefaultContext");
   p->addContext("ZdcNonDefaultContext");

   containers->Add(p);

   //    p->print();
}

FairParSet *MpdZdcContFact::createContainer(FairContainer *c)
{
   /** Calls the constructor of the corresponding parameter container.
    * For an actual context, which is not an empty string and not the default context
    * of this container, the name is concatinated with the context. */
   const char *name = c->GetName();
   FairParSet *p    = NULL;
   if (strcmp(name, "MpdZdcGeoPar") == 0) {
      p = new MpdZdcGeoPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());
   }
   return p;
}

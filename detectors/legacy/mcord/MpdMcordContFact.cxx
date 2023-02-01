/*
 * MpdMcordContFact.cxx
 *
 *  Created on: 21 maj 2019
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "MpdMcordContFact.h"
#include "FairRuntimeDb.h"
#include "MpdMcordGeoPar.h"
#include <iostream>

using namespace std;
static MpdMcordContFact gMpdMcordContFact;
MpdMcordContFact::MpdMcordContFact()
{
   fName  = "MpdMcordContFact";
   fTitle = "Factory for parameter containers in libMcord";
   setAllContainers();
   FairRuntimeDb::instance()->addContFactory(this);
}

void MpdMcordContFact::setAllContainers()
{
   FairContainer *p = new FairContainer("MpdMcordGeoPar", "Mcord	 Geometry Parameters", "TestDefaultContext");
   p->addContext("TestNonDefaultContext");
   containers->Add(p);
}

MpdMcordContFact::~MpdMcordContFact()
{
   // TODO Auto-generated destructor stub
}

FairParSet *MpdMcordContFact::createContainer(FairContainer *c)
{
   const char *name = c->GetName();
   cout << "-I- container name: " << name << endl;
   FairParSet *p = NULL;

   if (strcmp(name, "MpdMcordGeoPar") == 0)
      p = new MpdMcordGeoPar(c->getConcatName().Data(), c->GetTitle(), c->getContext());

   return p;
}

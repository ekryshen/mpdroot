// -------------------------------------------------------------------------
//                            MpdFieldContFact source file             -----
//        Created 23/07/13  by P. Batyuk (MPD) <batyuk@jinr.ru>        -----
//                          from MpdFieldContFact (PNDROOT)            -----
// -------------------------------------------------------------------------

/// Last modified: 24.07.2013, P.B.

#include "MpdFieldContFact.h"
#include "MpdFieldPar.h"
#include "MpdConstPar.h"
#include "FairParSet.h"
#include "FairRuntimeDb.h"
#include "MpdMultiFieldPar.h"
#include <iostream>

using namespace std;
static MpdFieldContFact gMpdFieldContFact;

// -----   Constructor   ---------------------------------------------------
MpdFieldContFact::MpdFieldContFact() {
  fName = "MpdFieldContFact";
  fTitle = "Factory for field parameter containers";
  SetAllContainers();
  FairRuntimeDb::instance()->addContFactory(this);
}

// -----   Destructor   ----------------------------------------------------
MpdFieldContFact::~MpdFieldContFact() { }

// -----   Create containers   ---------------------------------------------
FairParSet* MpdFieldContFact::createContainer(FairContainer* container) 
{
  const char* name = container->GetName();
  cout << "create MpdFieldPar container " << name << endl;
  FairParSet* set = NULL;
  if ( strcmp(name, "MpdFieldPar") == 0 ) {
    set = new MpdFieldPar( container->getConcatName().Data(),
			   container->GetTitle(),
			   container->getContext() );
  } 
  
  else if ( strcmp(name, "MpdConstPar") == 0 ) {
    set = new MpdConstPar( container->getConcatName().Data(),
			   container->GetTitle(),
			   container->getContext() );
  }
  
  else if ( strcmp(name, "MpdMultiFieldPar") == 0 ) {
    set = new MpdMultiFieldPar( container->getConcatName().Data(),
				container->GetTitle(),	
				container->getContext() );
  }  
  return set;
}

// -----   Set all containers (private)   ----------------------------------
void MpdFieldContFact::SetAllContainers() {
  FairContainer* con1 = new FairContainer("MpdFieldPar",
					  "Field parameter container",
					  "Default field");
  containers->Add(con1);
  
  FairContainer* con2 = new FairContainer("MpdConstPar",
					  "Const Field parameter container",
					  "Default field");
  containers->Add(con2);
  
  FairContainer* con3 = new FairContainer("MpdMultiFieldPar",
					  "Multiple Field parameter container",
					  "Default field");
  containers->Add(con3);
}

ClassImp(MpdFieldContFact)

/*
 * MpdFloorContFact.cxx
 *
 *  Created on: 21 maj 2019
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "MpdFloorContFact.h"
#include "FairRuntimeDb.h"
#include "MpdFloorGeoPar.h"
#include <iostream>

using namespace std;
static MpdFloorContFact gMpdFloorContFact;
MpdFloorContFact::MpdFloorContFact() {
	fName = "MpdFloorContFact";
	fTitle = "Factory for parameter containers in libFloor";
	setAllContainers();
        FairRuntimeDb::instance()->addContFactory(this);
}

void MpdFloorContFact::setAllContainers() {
    FairContainer* p = new FairContainer("MpdFloorGeoPar", "Floor	 Geometry Parameters", "TestDefaultContext");
    p->addContext("TestNonDefaultContext");
    containers->Add(p);
}

MpdFloorContFact::~MpdFloorContFact() {
	// TODO Auto-generated destructor stub
}

FairParSet* MpdFloorContFact::createContainer(FairContainer* c) {
	const char* name = c->GetName();
	cout << "-I- container name: " << name << endl;
        FairParSet* p = NULL;

	if(strcmp(name,"MpdFloorGeoPar")==0) p = new MpdFloorGeoPar(c->getConcatName().Data(),c->GetTitle(),c->getContext());

return p;
}

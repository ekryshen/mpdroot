/*
 * MpdMcordContFact.h
 *
 *  Created on: 21 maj 2019
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef MCORD_MCORD_MPDMCORDCONTFACT_H_
#define MCORD_MCORD_MPDMCORDCONTFACT_H_

#include "FairContFact.h"

class FairContainer;

class MpdMcordContFact : public FairContFact {
private:
   void setAllContainers();

public:
   MpdMcordContFact();
   virtual ~MpdMcordContFact();
   FairParSet *createContainer(FairContainer *);

   ClassDef(MpdMcordContFact, 0)
};

#endif /* MCORD_MCORD_MPDMCORDCONTFACT_H_ */

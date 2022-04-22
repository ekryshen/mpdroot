/*
 * MpdFreezoutGenerator.h
 *
 *  Created on: 30 gru 2016
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef SRC_MPDROOT_GENERATORS_MPDFREEZOUTGENERATOR_H_
#define SRC_MPDROOT_GENERATORS_MPDFREEZOUTGENERATOR_H_
#include "FairTask.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>

/**
 * class for registering freezout coordinates in MpdROOT, for unknown reason standard registering doesn't work with
 * Geant4
 */
class MpdFreezoutGenerator : public FairTask {
   static MpdFreezoutGenerator *fgInstance;
   TClonesArray                *fFreez;

public:
   MpdFreezoutGenerator();
   static MpdFreezoutGenerator *Instance();
   InitStatus                   Init();
   /**
    *
    * @return pointer to array
    */
   TClonesArray *GetArray() const { return fFreez; };
   virtual ~MpdFreezoutGenerator();
   ClassDef(MpdFreezoutGenerator, 1);
};

#endif /* SRC_MPDROOT_GENERATORS_MPDFREEZOUTGENERATOR_H_ */

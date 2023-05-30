/*
 * MpdFreezoutGenerator.cxx
 *
 *  Created on: 30 gru 2016
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "FairRunSim.h"
#include "FairRootManager.h"
#include "MpdFreezoutGenerator.h"
MpdFreezoutGenerator *MpdFreezoutGenerator::fgInstance = NULL;
MpdFreezoutGenerator::MpdFreezoutGenerator()
{
   if (fgInstance) {
      return;
   } else {
      FairRunSim *sim = FairRunSim::Instance();
      fFreez          = new TClonesArray("TLorentzVector", 1000);
      fgInstance      = this;
      sim->AddTask(this);
   }
}

MpdFreezoutGenerator *MpdFreezoutGenerator::Instance()
{
   if (fgInstance == NULL) {
      fgInstance = new MpdFreezoutGenerator();
   }
   return fgInstance;
}

InitStatus MpdFreezoutGenerator::Init()
{
   FairRootManager *mngr = FairRootManager::Instance();
   mngr->Register("Freezouts.", "Freezouts", fFreez, kTRUE);
   return kSUCCESS;
}

MpdFreezoutGenerator::~MpdFreezoutGenerator()
{
   delete fFreez;
}

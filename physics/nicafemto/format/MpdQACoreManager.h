/*
 * MpdQACoreManager.h
 *
 *  Created on: 15 lip 2021
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef MPDROOT_PHYSICS_NICAFEMTO_NICA_HELPERS_MPDQACOREMANAGER_H_
#define MPDROOT_PHYSICS_NICAFEMTO_NICA_HELPERS_MPDQACOREMANAGER_H_

#include "NicaQACoreManager.h"

class MpdQACoreManager : public NicaQACoreManager {
public:
   MpdQACoreManager();
   virtual void        SetRecoTrackCut(NicaTrackAna *ana, NicaQACoreManager::ePidCut cut,
                                       NicaQACoreManager::eParticleType primary, TString flag);
   virtual FairRunAna *GetRunAna(TString outFile, TString simFile, TString recoFile = "", TString parFile = "");
   virtual NicaEvent  *GetFormat(eFormatType type, eAnaType ana = eAnaType::kDefault);
   virtual ~MpdQACoreManager();
   ClassDef(MpdQACoreManager, 1);
};

#endif /* MPDROOT_PHYSICS_NICAFEMTO_NICA_HELPERS_MPDQACOREMANAGER_H_ */

/*
 * MpdQACoreManager.cxx
 *
 *  Created on: 15 lip 2021
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */

#include "MpdQACoreManager.h"

#include "NicaMiniDstSource.h"
#include "NicaMpdHbtComplexEvent.h"
#include "NicaMpdHbtEvent.h"
#include "NicaMpdMiniDstEvent.h"
#include "NicaMpdMiniDstFullEvent.h"
#include "NicaMpdMiniDstMcEvent.h"

#include <FairRunAna.h>
#include <TString.h>

MpdQACoreManager::MpdQACoreManager() {
  fEta[0] = -2;
  fEta[1] = 2;
  fPt[1]  = 4;
}

FairRunAna* MpdQACoreManager::GetRunAna(TString outFile, TString simFile, TString recoFile, TString parFile) {
  if (simFile.Length() == 0) {
    simFile  = recoFile;
    recoFile = "";
  }
  FairRunAna* run         = new FairRunAna();
  NicaMiniDstSource* file = new NicaMiniDstSource(simFile);
  run->SetSource(file);
  run->SetOutputFile(outFile);
  return run;
}

NicaEvent* MpdQACoreManager::GetFormat(eFormatType type, eAnaType ana) {
  switch (type) {
    case eFormatType::kComplex: {
      if (ana == eAnaType::kHbt) return new NicaMpdHbtComplexEvent();
      return new NicaMpdMiniDstFullEvent();
    } break;
    case eFormatType::kReco: {
      if (ana == eAnaType::kHbt) return new NicaMpdHbtEvent();
      return new NicaMpdMiniDstEvent();
    } break;
    case eFormatType::kSim: {
      return new NicaMpdMiniDstMcEvent();
    };
  }
  return nullptr;
}

MpdQACoreManager::~MpdQACoreManager() {
  // TODO Auto-generated destructor stub
}

void MpdQACoreManager::SetRecoTrackCut(NicaTrackAna* ana,
                                       NicaQACoreManager::ePidCut cut,
                                       NicaQACoreManager::eParticleType primary,
                                       TString flag) {}

/*
 * v0ana.C
 *
 *  Created on: 17 gru 2021
 *      Author: Daniel Wielanek
 *      E-mail: daniel.wielanek at pw.edu.pl
 *      Warsaw University of Technology, Faculty of Physics
 *      JINR,  Laboratory of High Energy Physics
 */

#ifndef __CLING__
#include <FairRunAna.h>
#include <vector>
#include <RtypesCore.h>
#include <TString.h>
#include "../../../../external/nicafemto/features/NicaHelix.h"
#include "../../../../mpddst/MpdMiniDstSource.h"
#include "../../../physics/common/v0/MpdV0CandidateCut.h"
#include "../../../physics/common/v0/MpdV0FinderBasic.h"
#include "../../../physics/common/minidst/MpdMiniDstSource.h"
#include "../../../physics/common/v0/MpdV0DaughterCutBasic.h"
#include "../../../physics/common/v0/MpdV0FinderHelix.h"
#include "../../../physics/common/v0/MpdV0Namespace.h"
#include "../../../physics/nicafemto/nica_helpers/MpdPIDOnTheFly.h"
#include "../../../external/nicafemto/features/NicaStd.h"

#include "MpdV0DaughterCutBasic.h"
#include "MpdV0CandidateCutKalman.h"
#endif

std::vector<TString> GetListOfFiles(TString path);
void                 v0ana()
{
   auto        list = GetListOfFiles("/media/daniel/hdd/minidst");
   FairRunAna *run  = new FairRunAna();

   MpdMiniDstSource *source = new MpdMiniDstSource(list[0]);
   for (int i = 1; i < list.size(); i++) {
      source->AddFile(list[i]);
   }
   run->SetSource(source);
   run->SetOutputFile("data.root");

   MpdPIDOnTheFly *pid = new MpdPIDOnTheFly();
   run->AddTask(pid);

   MpdV0FinderHelix *helix = new MpdV0FinderHelix();

   MpdV0DaughterCutBasic proton;
   proton.SetDcaMinCut(0.6, 0.6);
   proton.SetSigmaCut(-2, 2, MpdCommonV0::ESigmaType::kProtonSigma);
   MpdV0DaughterCutBasic pion;
   pion.SetDcaMinCut(0.6, 0.6);
   pion.SetSigmaCut(-2, 2, MpdCommonV0::ESigmaType::kPionSigma);
   pion.SetChargeCut(-1);

   helix->SetFirstDaughterCut(proton);
   helix->SetSecondDaughterCut(pion);

   MpdV0CandidateCutKalman lambda;
   lambda.SetDecayLenght(4., 1E+6);
   lambda.SetDcaBetweenDaughters(0, 1.4);
   lambda.SetMaxDcaPrimVtx(0, 2.5);
   helix->SetCandicateCut(lambda);

   helix->SaveV0s(kTRUE);

   run->AddTask(helix);
   run->Init();
   run->Run();
}

NicaHelix *hh; // this is only for loading namespace

std::vector<TString> GetListOfFiles(TString path)
{
   return NicaStd::GetListOfFiles(path, "root", kTRUE, kFALSE);
}

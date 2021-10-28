/*
 * minidst_ana.C
 *
 *  Created on: 29 wrz 2020
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */

#if !defined(__CINT__) && !defined(__CLING__)
#include "FairLogger.h"
#include "MpdBasicTrackCut.h"
#include "MpdPIDOnTheFly.h"
#include "MpdTpcMonitor.h"
#include "NicaConst.h"
#include "NicaFemto1DCF.h"
#include "NicaFemtoBasicAna.h"
#include "NicaFemtoCorrFuncKt.h"
#include "NicaFemtoSourceModelGauss.h"
#include "NicaFemtoWeightGeneratorLednicky.h"
#include "NicaMiniDstSource.h"
#include "NicaMpdDstMCEventTpcPads.h"
#include "NicaMpdMiniDstFullEvent.h"
#include "NicaTrackTpcToFCut.h"
#endif

NicaFemtoBasicAna* PrepPionAna() {
  NicaFemtoBasicAna* ana = new NicaFemtoBasicAna();
  ana->SetFormat(new NicaMpdMiniDstFullEvent());
  ana->SetPdg(211);
  NicaFemtoWeightGeneratorLednicky weight;
  weight.SetPairType(ENicaFemtoPairType::kPionPlusPionPlus);
  weight.SetStrongOff();
  weight.SetCoulOn();
  weight.SetQuantumOn();
  ana->SetWeight(weight);
  NicaFemtoSourceModelGauss1D gaus;
  gaus.SetRadius(5);
  NicaFemtoFreezoutGeneratorLCMS lcms;
  lcms.SetSourceModel(gaus);
  ana->SetFreezoutGenerator(lcms);
  // analysis in LCMS frame kt bins 0.2-0.4 and 0.4-0.6
  NicaFemtoCorrFuncKt kt(NicaFemto1DCF("cf", 100, 0, .15, ENicaFemtoKinematics::kLCMS), {0.2, 0.4, 0.6});
  ana->SetCorrFctn(kt);
  ana->SetOption(NicaTwoTrackAna::BackgroundOptionMixed());
  ana->SetMixSize(5);
  return ana;
}

/**
 * prepare cuts and cut monitors, for simplicity we use NicaCutsAndMonitors
 * that initialize needed cuts and property monitors
 * @param ana
 */
void prepTrackCuts(NicaFemtoBasicAna* ana) {
  MpdBasicTrackCut cut;
  cut.MakeCutMonitors("dca+tpc+tof+kin");  // create all available cut monitors
  cut.GetTpcMonitor()->SetAxisDeDx(200, 0, 1E+5);
  cut.GetTpcMonitor()->SetAxisP(100, 0, 2);
  cut.GetChargeCut()->SetMinAndMax(1);
  cut.GetPtCut()->SetMinMax(0.1, 3);
  cut.GetEtaCut()->SetMinMax(-1.2, 1.2);
  cut.GetTpcTofCut()->SetCharge(1);
  cut.GetTpcTofCut()->SetNHits(40, 90);
  cut.GetTpcTofCut()->SetActiveSigma("pi");
  cut.GetTpcTofCut()->SetThreshold(0.5);
  cut.GetTpcTofCut()->SetSigma(-3, 3, "pi");
  cut.GetTpcTofCut()->SetM2(-0.05, 0.05);
  /**
   * we are using "complex format" therefore we have to add detector cuts as
   * "real" part of cuts
   */
  cut.SetOptionForAllCuts("re");
  ana->AddCutsAndMonitors(cut);
}

void minidst_ana(TString inFile, TString outFile = "hbt.root") {
  // prepare fairroot framework
  FairRunAna* run = new FairRunAna();
  run->SetOutputFile(outFile);

  FairLogger* log = FairLogger::GetLogger();
  log->SetColoredLog(kTRUE);

  NicaMiniDstSource* source = new NicaMiniDstSource(inFile);

  run->SetSource(source);
  // prepare analysis
  NicaFemtoBasicAna* ana = PrepPionAna();

  prepTrackCuts(ana);
  // run->AddTask(new MpdPIDOnTheFly());
  run->AddTask(ana);
  run->Init();
  run->Run(0, 5000);
}

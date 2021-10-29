/*
 * v0_finder.C
 *
 *  Created on: 29 wrz 2020
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */

/**
 * macro for lambda reconstruction by using STAR-like cuts
 */
#ifndef __CLING__

#include "FairRunAna.h"
#include "MpdPIDOnTheFly.h"
#include "MpdTpcMonitor.h"
#include "NicaConst.h"
#include "NicaCutMonitorX.h"
#include "NicaCutMonitorXY.h"
#include "NicaJobs.h"
#include "NicaMiniDstFullV0Event.h"
#include "MpdMiniDstSource.h"
#include "NicaMpdMiniDstEvent.h"
#include "NicaMpdMiniDstFullEvent.h"
#include "NicaTrackAna.h"
#include "NicaTrackBasicMCCut.h"
#include "NicaTrackChargeCut.h"
#include "NicaTrackDCACut.h"
#include "NicaTrackPCut.h"
#include "NicaTrackPdgCut.h"
#include "NicaTrackToFMass2Cut.h"
#include "NicaTrackTpcCut.h"
#include "NicaTwoTrackAna.h"
#include "NicaV0BasicFinder.h"
#include "NicaV0CandBasicCut.h"
#include "NicaV0CandidateHelix.h"
#endif

#define N_FILES 400

/**
 * positive daughter's cuts (proton)
 * @param task
 * @param tof
 */
void SetPos(NicaV0BasicFinder* task, Bool_t tof) {
  NicaTrackTpcCut tpc;
  tpc.SetCharge(1);
  tpc.SetActiveSigma(tpc.ProtonSigma());
  tpc.SetModeNotBad();
  tpc.SetNHits(20, 90);
  tpc.SetMinMax(-3, 3, tpc.ProtonSigma());
  task->AddPosDauCut(tpc);
  NicaTrackPCut p;
  task->AddPosDauCut(p);
  NicaCutMonitorXY mon(p.CutName(), 0, tpc.CutName(), tpc.DeDx());
  mon.SetXaxis(200, 0, 4);
  mon.SetYaxis(200, 0, 50E+3);
  task->AddPosDauCutMon(mon);
  NicaTrackDCACut dca;
  dca.SetMinMax(0.6, 1E+9, dca.DCA());
  task->AddPosDauCut(dca);

  NicaCutMonitorX mons(tpc.CutName(), tpc.ProtonSigma());
  mons.SetXaxis(100, -10, 10);
  // task->AddPosDauCutMon(mons);
  NicaCutMonitorX dcaMon(dca.CutName(), dca.DCA());
  dcaMon.SetXaxis(100, 0, 10);
  task->AddPosDauCutMon(dcaMon);

  if (tof) {
    NicaTrackToFMass2Cut mass;
    Double_t mP = NicaConst::ProtonMass() * NicaConst::ProtonMass();
    mass.SetMinMax(mP * 0.8, mP * 1.2);
    task->AddPosDauCut(mass);
    NicaCutMonitorXY m2(p.CutName(), 0, mass.CutName(), 0);
    m2.SetXaxis(200, 0, 4);
    m2.SetYaxis(200, 0, 1.2);
    task->AddPosDauCutMon(m2);
  }
}
/**
 * negative daughter's cuts (pion)
 * @param task
 * @param tof
 */
void SetNeg(NicaV0BasicFinder* task, Bool_t tof) {
  NicaTrackTpcCut tpc;
  tpc.SetCharge(-1);
  tpc.SetActiveSigma(tpc.PionSigma());
  tpc.SetModeNotBad();
  tpc.SetNHits(20, 90);
  tpc.SetMinMax(-3, 3, tpc.PionSigma());
  task->AddNegDauCut(tpc);
  NicaTrackPCut p;
  task->AddNegDauCut(p);
  NicaCutMonitorXY mon(p.CutName(), 0, tpc.CutName(), tpc.DeDx());
  mon.SetXaxis(200, 0, 4);
  mon.SetYaxis(200, 0, 50E+3);
  task->AddNegDauCutMon(mon);
  NicaCutMonitorX mons(tpc.CutName(), tpc.PionSigma());
  mons.SetXaxis(100, -10, 10);
  //   task->AddNegDauCutMon(mons);
  NicaTrackDCACut dca;
  dca.SetMinMax(1.5, 1E+9, dca.DCA());
  NicaCutMonitorX dcaMon(dca.CutName(), dca.DCA());
  dcaMon.SetXaxis(100, 0, 10);
  task->AddNegDauCutMon(dcaMon);
  task->AddNegDauCut(dca);
  if (tof) {
    NicaTrackToFMass2Cut mass;
    Double_t mP = NicaConst::PionPlusMass() * NicaConst::PionPlusMass();
    mass.SetMinMax(-mP, mP * 1.2);
    task->AddPosDauCut(mass);
    task->AddNegDauCut(mass);
    NicaCutMonitorXY m2(p.CutName(), 0, mass.CutName(), 0);
    task->AddNegDauCutMon(m2);
  }
}
/**
 * candidates cuts
 * @param task
 */
void SetLambCut(NicaV0BasicFinder* task) {
  // szybkie ciecie wybiera lambdy
  NicaV0CandBasicCut basF;
  basF.SetMinMax(0, 0.8, basF.DCA1to2());
  basF.SetMinMax(0, 0.4, basF.DCAToPV());
  basF.SetMinMax(7, 1E+9, basF.DecLength());
  /**
   * fast cuts are used to remove all wrong  candidates without mass
   * restriction, thanks to this we have in "failed" only lambdas with wrong
   * mass and can easy calculate the contamination "double" option is used to
   * prevent nicafemto from removing cuts with the same names
   */
  task->AddCandicateCut(basF, "fast+double");
  NicaV0CandBasicCut bas;
  bas.SetMinMax(0, 0.8, bas.DCA1to2());
  bas.SetMinMax(0, 0.4, bas.DCAToPV());
  bas.SetMinMax(7, 1E+9, bas.DecLength());
  bas.SetMinMax(NicaConst::LambdaMass() - 0.004, NicaConst::LambdaMass() + 0.004, bas.InvMass());
  NicaCutMonitorX mon(bas.CutName(), bas.InvMass());
  mon.SetXaxis(100, 1, 1.2);
  NicaCutMonitorXY mon2(bas.CutName(), bas.AlphaArm(), bas.CutName(), bas.PtArm());
  mon2.SetXaxis(100, -1, 1);
  mon2.SetYaxis(100, 0, 0.3);
  task->AddCandicateCut(bas, "double");
  task->AddCandicateCutMon(mon);
  task->AddCandicateCutMon(mon2);
}

void v0_lambda(TString inFile, TString outFile = "v0.root") {
  FairRunAna* ana = new FairRunAna();

  MpdMiniDstSource* source = new MpdMiniDstSource(inFile);

  ana->SetSource(source);
  ana->SetOutputFile(outFile);
  // two v0 findes with and without TOF
  NicaV0BasicFinder* v0_1 = new NicaV0BasicFinder();
  NicaV0BasicFinder* v0_2 = new NicaV0BasicFinder();
  SetPos(v0_1, kFALSE);
  SetPos(v0_2, kTRUE);
  SetNeg(v0_1, kFALSE);
  SetNeg(v0_2, kTRUE);
  SetLambCut(v0_1);
  SetLambCut(v0_2);
  // all trajectiories in MPD are helices - this function must be here
  v0_1->SetV0Container(NicaV0CandidateHelix());
  v0_2->SetV0Container(NicaV0CandidateHelix());
  // set format of input data
  v0_1->SetFormat(new NicaMpdMiniDstEvent());
  v0_2->SetFormat(new NicaMpdMiniDstEvent());

  // to fix n-sigma calculation (not needed for newest productions)
  MpdPIDOnTheFly* pid = new MpdPIDOnTheFly();
  ana->AddTask(pid);

  NicaTwoTrackAna* pair = new NicaTwoTrackAna();
  pair->SetFormat(new NicaMiniDstFullV0Event());

  ana->AddTask(v0_1);
  ana->AddTask(v0_2);
  ana->AddTask(pair);
  ana->Init();
  ana->Run(20000);
}

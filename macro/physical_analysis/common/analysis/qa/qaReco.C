/*
 * qaReco.C
 *
 *  Created on: 5 lis 2021
 *      Author: dr Daniel Wielanek
 *      E-mail: daniel.wielanek at pw.edu.pl
 *      Warsaw University of Technology, Faculty of Physics
 *      JINR,  Laboratory of High Energy Physics
 */

#if !defined(__CINT__) && !defined(__CLING__)
#include "FairLogger.h"
#include "FairRunAna.h"
#include "MpdBasicTrackCut.h"
#include "MpdMiniDstSource.h"
#include "MpdPIDOnTheFly.h"
#include "MpdTpcMonitor.h"
#include "NicaConst.h"
#include "NicaDataFormat.h"
#include "NicaEventImpactParameterCut.h"
#include "NicaEventMultiplicityCut.h"
#include "NicaFemto1DCF.h"
#include "NicaFemtoBasicAna.h"
#include "NicaFemtoCorrFuncKt.h"
#include "NicaFemtoSourceModelGauss.h"
#include "NicaFemtoWeightGeneratorLednicky.h"
#include "NicaMpdDstMCEvent.h"
#include "NicaMpdDstMCEventTpcPads.h"
#include "NicaMpdEvent.h"
#include "NicaMpdMiniDstEvent.h"
#include "NicaMpdMiniDstFullEvent.h"
#include "NicaQATrackTask.h"
#include "NicaTrackNullCut.h"
#include "NicaTrackPdgCut.h"
#include "NicaTrackTpcToFCut.h"
#include "TString.h"
#endif

/**
 * macro for processing reconstructed data because we use purely reconstructed events, we don't need to switch with
 * NicaDataFieldID::ReStep. See also qa.C.
 */

#define MINIDST

NicaQAPlot GetEventQA();

NicaQAPlot GetTrackQA();


#ifdef MINIDST
void qaReco(TString inFile = "$VMCWORKDIR/macro/mpd/mpddst.MiniDst.root", TString outFile = "qaReco.root") {
  FairRunAna* ana          = new FairRunAna();
  MpdMiniDstSource* source = new MpdMiniDstSource(inFile);
#else
void qaReco(TString inFile = "$VMCWORKDIR/macro/mpd/mpddst.root", TString outFile = "qaReco.root") {
  FairRunAna* ana        = new FairRunAna();
  FairFileSource* source = new FairFileSource(inFile);
#endif

  ana->SetSource(source);
  ana->SetOutputFile(outFile);

  NicaQATrackTask* trackTask = new NicaQATrackTask();
#ifdef MINIDST
  trackTask->SetFormat(new NicaMpdMiniDstEvent());
#else
  trackTask->SetFormat(new NicaMpdEvent());
#endif
  // arbitray by number of tracks
  Double_t Centrality[4] = {2000, 600, 200, 0};

  NicaEventMultiplicityCut centCut;

  trackTask->SetEventCollectionNames({"Central", "MidCentral", "Peripheral"});
  trackTask->SetTrackCollectionNames({"Positve", "Negative"});
  for (int i = 0; i < 3; i++) {
    centCut.SetMinMax(Centrality[i + 1], Centrality[i]);
    centCut.SetCollectionID(i);
    trackTask->AddCut(centCut);
    trackTask->SetQAPlot(GetEventQA());
  }

  NicaTrackChargeCut charge;
  charge.SetMinAndMax(1);
  trackTask->AddCut(charge);
  charge.SetMinAndMax(-1);
  charge.SetCollectionID(1);
  trackTask->AddCut(charge);
  trackTask->SetQAPlot(GetTrackQA());

  ana->AddTask(trackTask);
  ana->Init();
  ana->Run(4000);
}

using namespace NicaDataFieldID;

NicaQAPlot GetEventQA() {
  NicaQAPlot plot(ENicaCutUpdate::kEventUpdate);
  plot.AddTH1("centrality", EEvent::kTracksNo, 100, 0, 3000);
  plot.AddTH2("vertex", EEvent::kVertexZ, EEvent::kVertexXY, 100, -100, 100, 100, -5, 5);
  plot.AddTH2("vertexXY", EEvent::kVertexX, EEvent::kVertexY, 50, -2.5, 2.5, 50, -2.5, 2.5);
  return plot;
}

NicaQAPlot GetTrackQA() {
  NicaQAPlot plot(ENicaCutUpdate::kTrackUpdate);
  plot.AddTH2("KinematicsReco", ETrack::kEta, ETrack::kPt, 100, -2, 2, 100, 0, 4);
  plot.AddTH2("TpcMon", ETrack::kP, EExpTrack::kTpcDedx, 200, 0, 3, 1000, 0, 1E+5);
  plot.AddTH2("TofMon", ETrack::kP, EExpTrack::kTofM2, 200, 0, 3, 200, -.1, 2);
  plot.AddTH2("TofMon2", ETrack::kP, EExpTrack::kToFBeta, 200, 0, 3, 200, -0.1, 1.2);
  plot.AddTH1("NtpcHits", EExpTrack::kTpcNHits, 56, 0, 56);
  plot.AddTH1("Chi2", EExpTrack::kChi2, 40, 0, 60);
  return plot;
}

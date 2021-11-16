/*
 * qa.C
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
#include "NicaMpdMCEvent.h"
#include "NicaMpdMcEvent.h"
#include "NicaMpdMiniDstEvent.h"
#include "NicaMpdMiniDstFullEvent.h"
#include "NicaMpdMiniDstMcEvent.h"
#include "NicaQATrackTask.h"
#include "NicaTrackNullCut.h"
#include "NicaTrackPdgCut.h"
#include "NicaTrackTpcToFCut.h"
#include "TString.h"
#endif

/**
 * macro for processing simulated data, here we use only MC-data therefore, we do not have to switch to MC by
 * NicaDataFieldID::ImStep, see also qa.C
 */
#define MINIDST

NicaQAPlot GetEventQA();

NicaQAPlot GetTrackQA();

#ifdef MINIDST
void qaSim(TString inFile = "$VMCWORKDIR/macro/mpd/mpddst.MiniDst.root", TString outFile = "qaSim.root") {
  FairRunAna* ana          = new FairRunAna();
  MpdMiniDstSource* source = new MpdMiniDstSource(inFile);
#else
void qaSim(TString inFile = "$VMCWORKDIR/macro/mpd/mpddst.root", TString outFile = "qaSim.root") {
  FairRunAna* ana        = new FairRunAna();
  FairFileSource* source = new FairFileSource(inFile);
#endif

  ana->SetSource(source);
  ana->SetOutputFile(outFile);

  NicaQATrackTask* trackTask = new NicaQATrackTask();
#ifdef MINIDST
  trackTask->SetFormat(new NicaMpdMiniDstMcEvent());
#else
  trackTask->SetFormat(new NicaMpdMcEvent());
#endif
  // arbitray by number of tracks
  Double_t Centrality[4] = {0, 5, 9.10, 20};
  Int_t Pids[6]          = {211, -211, 321, -321, 2212, -2212};

  NicaEventImpactParameterCut bCut;

  trackTask->SetEventCollectionNames({"Central", "MidCentral", "Peripheral"});
  trackTask->SetTrackCollectionNames({"Pion+", "Pion-", "Kaon+", "Kaon-", "Proton", "AntiProton", "Unmatched"});
  for (int i = 0; i < 3; i++) {
    bCut.SetMinMax(Centrality[i], Centrality[i + 1]);
    bCut.SetCollectionID(i);
    trackTask->AddCut(bCut);
    trackTask->SetQAPlot(GetEventQA());
  }

  for (int j = 0; j < 6; j++) {
    NicaTrackPdgCut pid;
    pid.SetMinAndMax(Pids[j]);
    pid.SetCollectionID(j);
    trackTask->AddCut(pid);
    trackTask->SetQAPlot(GetTrackQA());
  }
  NicaTrackNullCut nullcut;
  nullcut.AcceptOnlyNotMatched();
  nullcut.SetCollectionID(6);
  trackTask->AddCut(nullcut);
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
  plot.AddTH1("imp", EMcEvent::kB, 100, 0, 20);
  plot.AddTH2("multiVsImp", EEvent::kTracksNo, EMcEvent::kB, 1000, 0, 2000, 100, 0, 20);
  return plot;
}

NicaQAPlot GetTrackQA() {
  NicaQAPlot plot(ENicaCutUpdate::kTrackUpdate);
  plot.AddTH2("KinematicsReco", ETrack::kEta, ETrack::kPt, 100, -2, 2, 100, 0, 4);
  plot.AddTH1("Rapidity", ETrack::kRapidity, 100, -2, 2);
  plot.AddTH1("pT", ETrack::kPt, 200, 0, 4);
  return plot;
}

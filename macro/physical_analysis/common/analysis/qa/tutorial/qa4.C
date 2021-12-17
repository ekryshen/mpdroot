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
#include "NicaEventVirtualCut.h"
#include "NicaFemto1DCF.h"
#include "NicaFemtoBasicAna.h"
#include "NicaFemtoCorrFuncKt.h"
#include "NicaFemtoSourceModelGauss.h"
#include "NicaFemtoWeightGeneratorLednicky.h"
#include "NicaMpdDstMCEvent.h"
#include "NicaMpdDstMCEventTpcPads.h"
#include "NicaMpdMiniDstFullEvent.h"
#include "NicaMpdMiniDstMcEvent.h"
#include "NicaQATrackTask.h"
#include "NicaTrackNullCut.h"
#include "NicaTrackPdgCut.h"
#include "NicaTrackTpcToFCut.h"
#include "TString.h"
#endif

/**
 * macro for processing "complex events", in complex events reconstructed data is represented by "real" part whereas simulated
 * data is presented by "imaginary" part To switch between "real" or "imaginary" fields, you should add NicaDataFieldID::ReStep
 * and NicaDataFieldID::ImStep. Otherwise "default" value will be presented (corresponds to real part). Note - not all fields are
 * supported:
 * List of supported enums:
 * - Works only with MC-data
 *  -#NicaDataFieldID::EMcEvent
 *  -#NicaDataFieldID::EMcTrack
 * - works only with reconstructed data
 *  -#NicaDataFieldID::EExpEvent
 *  -#NicaDataFieldID::EExpTrack
 * - works with any data
 *  -#NicaDataFieldID::ETrack
 *  -#NicaDataFieldID::EEvent
 * - works only with complex format
 *  -#NicaDataFieldID::EComplexEvent
 *  -#NicaDataFieldID::EComplexTrack
 *
 *
 *  To create www report use nica-report [root_file] [output_dir]
 *  Note: the histograms are blockedby CORS policy - you have to upload your directory to some web server or use some "simple
 * server" like "python -m SimpleHTTPServer 8069" to browse such site enter 0.0.0.0:8069
 */

NicaQAPlot GetEventQA();

NicaQAPlot GetTrackQA();

void qa4(TString inFile = "$VMCWORKDIR/macro/mpd/mpddst.MiniDst.root", TString outFile = "qa.root") {
  FairRunAna* ana          = new FairRunAna();
  MpdMiniDstSource* source = new MpdMiniDstSource(inFile);


  ana->SetSource(source);
  ana->SetOutputFile(outFile);

  NicaQATrackTask* trackTask = new NicaQATrackTask();

  trackTask->SetFormat(new NicaMpdMiniDstFullEvent());

  trackTask->SetQAPlot(GetTrackQA());
  trackTask->SetQAPlot(GetEventQA());
  trackTask->SetEventCollectionNames({"CentralE", "SemicentralE"});
  trackTask->SetTrackCollectionNames({"Pions", "Protons"});


  NicaEventImpactParameterCut central;
  central.SetMinMax(0, 3.5);
  central.SetCollectionID(0);
  NicaEventImpactParameterCut semicentral;
  semicentral.SetMinMax(5, 10);
  semicentral.SetCollectionID(1);

  trackTask->AddCut(central, "im");  // cut on "imaginary part of data"
  trackTask->AddCut(semicentral, "im");


  NicaTrackPdgCut pion, proton;
  pion.SetMinAndMax(NicaConst::PionPlusPID());
  proton.SetMinAndMax(NicaConst::ProtonPID());

  trackTask->AddCut(pion, "{0}+im");  // you can specify collection ID by {}
  trackTask->AddCut(proton, "{1}+im");


  ana->AddTask(trackTask);
  ana->Init();
  ana->Run(40000);
}

using namespace NicaDataFieldID;

NicaQAPlot GetEventQA() {
  NicaQAPlot plot(ENicaCutUpdate::kEventUpdate);
  plot.AddTH1("imp", ImStep + EMcEvent::kB, 100, 0, 20);
  return plot;
}

NicaQAPlot GetTrackQA() {
  NicaQAPlot plot(ENicaCutUpdate::kTrackUpdate);
  plot.AddTH2("KinematicsMC", ImStep + ETrack::kEta, ImStep + ETrack::kPt, 100, -2, 2, 100, 0, 4);
  plot.AddTH2("KinematicsReco", ReStep + ETrack::kEta, ReStep + ETrack::kPt, 100, -2, 2, 100, 0, 4);
  return plot;
}

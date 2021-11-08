/*
 * qa.C
 *
 *  Created on: 5 lis 2021
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#if !defined(__CINT__) && !defined(__CLING__)
#include "FairLogger.h"
#include "MpdBasicTrackCut.h"
#include "MpdMiniDstSource.h"
#include "MpdPIDOnTheFly.h"
#include "MpdTpcMonitor.h"
#include "NicaConst.h"
#include "NicaDataFormat.h"
#include "NicaEventImpactParameterCut.h"
#include "NicaFemto1DCF.h"
#include "NicaFemtoBasicAna.h"
#include "NicaFemtoCorrFuncKt.h"
#include "NicaFemtoSourceModelGauss.h"
#include "NicaFemtoWeightGeneratorLednicky.h"
#include "NicaMpdDstMCEvent.h"
#include "NicaMpdDstMCEventTpcPads.h"
#include "NicaMpdMiniDstFullEvent.h"
#include "NicaQATrackTask.h"
#include "NicaTrackPdgCut.h"
#include "NicaTrackTpcToFCut.h"
#endif

// TString CentralityNames[3] = {"Central", "MidCentral", "Peripheral"};
// TString ParticleSpieces[6] = {"Pion+", "Pion-", "Kaon+", "Kaon-", "Proton", "AntiProton"};

#define MINIDST


NicaQAPlot GetEventQA(Int_t ev);

NicaQAPlot GetTrackQA(Int_t ev, Int_t track);

void qa(TString inFile = "", TString outFile = "qa.root") {

  FairRunAna* ana = new FairRunAna();

#ifdef MINIDST
  std::ifstream fileList;
  fileList.open("/media/daniel/Baza/postdoc/2021/mpd/vhlle_prod/lista.txt");
  TString file;
  fileList >> file;
  file = Form("/media/daniel/Baza/postdoc/2021/mpd/vhlle_prod/%s", file.Data());
  cout << file << endl;
  MpdMiniDstSource* source = new MpdMiniDstSource(file);
  for (int i = 1; i < 100; i++) {
    fileList >> file;
    file = Form("/media/daniel/Baza/postdoc/2021/mpd/vhlle_prod/%s", file.Data());
    source->AddFile(file);
  }
  // MpdMiniDstSource* source = new MpdMiniDstSource(inFile);

#else
  FairFileSource* source = new FairFileSource(inFile);
  trackTask->SetFormat(new NicaMpdDstMCEvent());
#endif

  ana->SetSource(source);
  ana->SetOutputFile(outFile);

  NicaQATrackTask* trackTask = new NicaQATrackTask();
  trackTask->SetFormat(new NicaMpdMiniDstFullEvent());

  Double_t Centrality[4] = {0, 5, 9.10};
  Int_t Pids[6]          = {211, -211, 321, -321, 2212, -2212};

  NicaEventImpactParameterCut bCut;

  trackTask->SetEventCollectionNames({"Central", "MidCentral", "Peripheral"});
  trackTask->SetTrackCollectionNames({"Pion+", "Pion-", "Kaon+", "Kaon-", "Proton", "AntiProton"});
  for (int i = 0; i < 3; i++) {
    bCut.SetMinMax(Centrality[i], Centrality[i + 1]);
    bCut.SetCollectionID(i);
    trackTask->AddCut(bCut, "im");
    trackTask->SetQAPlot(GetEventQA(i));
  }

  for (int j = 0; j < 6; j++) {
    NicaTrackPdgCut pid;
    pid.SetMinAndMax(Pids[j]);
    pid.SetCollectionID(j);
    trackTask->AddCut(pid, "im");
    trackTask->SetQAPlot(GetTrackQA(0, j));
  }

  ana->AddTask(trackTask);
  ana->Init();
  ana->Run(4000);
}

using namespace NicaDataFieldID;

NicaQAPlot GetEventQA(Int_t ev) {
  NicaQAPlot plot(ENicaCutUpdate::kEventUpdate);
  plot.AddTH2("centrality", ImStep + EEvent::kTracksNo, ImStep + EMcEvent::kB, 100, 0, 3000, 100, 0, 20);
  plot.AddTH2("vertex", ImStep + EEvent::kVertexZ, ReStep + EEvent::kVertexZ, 100, -100, 100, 100, -100, 100);
  plot.AddTH1("imp", ImStep + EMcEvent::kB, 100, 0, 20);
  return plot;
}

NicaQAPlot GetTrackQA(Int_t ev, Int_t track) {
  NicaQAPlot plot(ENicaCutUpdate::kTrackUpdate);
  // plot MC values
  plot.AddTH2("KinematicsMC", ImStep + ETrack::kEta, ImStep + ETrack::kPt, 100, -2, 2, 100, 0, 4);
  // plot Reco values
  plot.AddTH2("KinematicsReco", ReStep + ETrack::kEta, ReStep + ETrack::kPt, 100, -2, 2, 100, 0, 4);
  plot.AddTH2("TpcMon", ReStep + ETrack::kP, ReStep + EExpTrack::kTpcDedx, 200, 0, 3, 1000, 0, 1E+5);
  plot.AddTH2("TofMon", ReStep + ETrack::kP, ReStep + EExpTrack::kTofM2, 200, 0, 3, 200, -.1, 2);
  plot.AddTH2("TofMon2", ReStep + ETrack::kP, ReStep + EExpTrack::kToFBeta, 200, 0, 3, 200, -0.1, 1.2);
  plot.AddTH1("NtpcHits", ReStep + EExpTrack::kTpcNHits, 56, 0, 56);
  plot.AddTH1("Chi2", ReStep + EExpTrack::kChi2, 40, 0, 20);

  // plot "complex values"
  plot.AddTH2("P reso", ReStep + ETrack::kP, EComplexTrack::kDeltaP, 100, 0, 2, 100, -0.1, 0.1);
  plot.AddTH2("Phi reso", ReStep + ETrack::kP, EComplexTrack::kDeltaPhi, 100, 0, 2, 100, -0.2, 0.2);
  plot.AddTH2("Theta reso", ReStep + ETrack::kP, EComplexTrack::kDeltaTheta, 100, 0, 2, 100, -.2, 0.2);

  return plot;
}

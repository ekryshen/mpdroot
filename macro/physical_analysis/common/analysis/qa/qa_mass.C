/*
 * qa_mass.C
 *
 *  Created on: 16 lis 2021
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
#include "NicaFemto1DCF.h"
#include "NicaFemtoBasicAna.h"
#include "NicaFemtoCorrFuncKt.h"
#include "NicaFemtoSourceModelGauss.h"
#include "NicaFemtoWeightGeneratorLednicky.h"
#include "NicaMpdDstMCEvent.h"
#include "NicaMpdDstMCEventTpcPads.h"
#include "NicaMpdMiniDstFullEvent.h"
#include "NicaQATrackTask.h"
#include "NicaTrackNullCut.h"
#include "NicaTrackPdgCut.h"
#include "NicaTrackTpcToFCut.h"
#include "TString.h"
#endif

/**
 * example macro for processing many files stored in a text file with a list with full paths to file
 * This code uses "complex events" in complex events reconstructed data is represented by "real" part whereas simulated data is
 * presented by "imaginary" part
 * See also qa.C
 */

NicaQAPlot GetEventQA();

NicaQAPlot GetTrackQA();

void qa_mass(TString inFile = "lista.txt", TString outFile = "qaMass.root") {

  FairRunAna* ana = new FairRunAna();
  std::ifstream fileList;
  fileList.open(inFile);
  TString file;
  fileList >> file;
  cout << file << endl;
  MpdMiniDstSource* source = new MpdMiniDstSource(file);
  for (int i = 1; i < 100; i++) {
    fileList >> file;
    cout << file << endl;
    source->AddFile(file);
  }

  ana->SetSource(source);
  ana->SetOutputFile(outFile);

  NicaQATrackTask* trackTask = new NicaQATrackTask();

  trackTask->SetFormat(new NicaMpdMiniDstFullEvent());

  Double_t Centrality[4] = {0, 5, 9.10, 20};
  Int_t Pids[6]          = {211, -211, 321, -321, 2212, -2212};

  NicaEventImpactParameterCut bCut;

  trackTask->SetEventCollectionNames({"Central", "MidCentral", "Peripheral"});
  trackTask->SetTrackCollectionNames({"Pion+", "Pion-", "Kaon+", "Kaon-", "Proton", "AntiProton", "Unmatched"});
  for (int i = 0; i < 3; i++) {
    bCut.SetMinMax(Centrality[i], Centrality[i + 1]);
    bCut.SetCollectionID(i);
    trackTask->AddCut(bCut, "im");
    trackTask->SetQAPlot(GetEventQA());
  }

  for (int j = 0; j < 6; j++) {
    NicaTrackPdgCut pid;
    pid.SetMinAndMax(Pids[j]);
    pid.SetCollectionID(j);
    trackTask->AddCut(pid, "im");
    trackTask->SetQAPlot(GetTrackQA());
  }
  NicaTrackNullCut nullcut;
  nullcut.AcceptOnlyNotMatched();
  nullcut.SetCollectionID(6);
  trackTask->AddCut(nullcut);
  trackTask->SetQAPlot(GetTrackQA());

  ana->AddTask(trackTask);
  ana->Init();
  ana->Run(40000);
}

using namespace NicaDataFieldID;

NicaQAPlot GetEventQA() {
  NicaQAPlot plot(ENicaCutUpdate::kEventUpdate);
  plot.AddTH2("centrality", ImStep + EEvent::kTracksNo, ImStep + EMcEvent::kB, 100, 0, 3000, 100, 0, 20);
  plot.AddTH2("vertex", ImStep + EEvent::kVertexZ, ReStep + EEvent::kVertexZ, 100, -100, 100, 100, -100, 100);
  plot.AddTH1("imp", ImStep + EMcEvent::kB, 100, 0, 20);
  plot.AddTH2("vertexXY", ReStep + EEvent::kVertexX, ReStep + EEvent::kVertexY, 50, -5, 5, 50, -5, 5);
  plot.AddTH2("multiVsImp", ReStep + EEvent::kTracksNo, ImStep + EMcEvent::kB, 1000, 0, 2000, 100, 0, 20);
  return plot;
}

NicaQAPlot GetTrackQA() {
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

/*
 * showQA.C
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
#include "NicaAnaFile.h"
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
#include "NicaMpdMiniDstEvent.h"
#include "NicaMpdMiniDstFullEvent.h"
#include "NicaQATrackTask.h"
#include "NicaTrackNullCut.h"
#include "NicaTrackPdgCut.h"
#include "NicaTrackTpcToFCut.h"
#include "TString.h"
#include <TCanvas.h>
#endif
NicaAnaFile* file;

/**
 * example macro that show how to acces the QA histograms from root macro
 */


void showQA(TString inFile = "qaReco.root") {

  NicaAnaFile* file      = new NicaAnaFile(inFile);
  TString eventLabels[3] = {"Central", "MidCentral", "Peripheral"};
  TString trackLabels[2] = {"Positve", "Negative"};

  /** get some event plot **/
  Int_t eventCol           = 0;
  NicaQAPlotReport* report = (NicaQAPlotReport*) file->GetMainObject(eventLabels[eventCol]);

  TCanvas* c = new TCanvas();
  report->Get2D(1)->Draw("colz");

  /** get some track plot **/
  Int_t trackCol = 0;

  NicaQAPlotReport* report2 =
    (NicaQAPlotReport*) file->GetMainObject(Form("%s %s", eventLabels[eventCol].Data(), trackLabels[trackCol].Data()));

  // draw first one-dimensional plot
  TCanvas* c1  = new TCanvas();
  Int_t oneDim = report2->GetSize1D();
  c1->Divide(1, oneDim);
  for (int i = 0; i < oneDim; i++) {
    c1->cd(i + 1);
    report2->Get1D(i)->Draw();
  }

  // draw all 2d plots
  TCanvas* c2 = new TCanvas();
  // first calcualte the canvas size
  Int_t twoDim = report2->GetSize2D();
  Int_t x_size = twoDim;
  Int_t y_size = 1;
  if (twoDim > 3) {
    x_size = TMath::Sqrt(twoDim);
    y_size = TMath::Ceil(twoDim / x_size);
  }

  c2->Divide(x_size, y_size);
  // then plot stuff
  for (int i = 0; i < twoDim; i++) {
    c2->cd(i + 1);
    report2->Get2D(i)->Draw("colz");
  }
}

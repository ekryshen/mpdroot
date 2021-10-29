/*
 * v0_display.C
 *
 *  Created on: 29 wrz 2020
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#if !defined(__CINT__) && !defined(__CLING__)
#include <TCanvas.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TMath.h>
#include <TStyle.h>
#include <TVirtualPad.h>
#include "NicaAnaFile.h"
#endif

/**
 * example how to read cut monitors from analysis
 */

void v0_display() {
  gStyle->SetOptStat(0);
  NicaAnaFile *f = new NicaAnaFile("v0.root");
  f->SwitchPackage(0);  // switch to first task v0 without TOF
  TH2D *armenteros_tpc =
      (TH2D *)f->GetHistogramPassed(ENicaCutUpdate::kTwoTrackUpdate, 0, 1);
  TH1D *mass_tpc =
      (TH1D *)f->GetHistogramPassed(ENicaCutUpdate::kTwoTrackUpdate, 0, 0);
  TH1D *mass_tpc_failed =
      (TH1D *)f->GetHistogramFailed(ENicaCutUpdate::kTwoTrackUpdate, 0, 0);
  mass_tpc->Add(mass_tpc_failed);
  f->SwitchPackage(1);  // switch to first task v0 with TOF
  TH2D *armenteros_tof =
      (TH2D *)f->GetHistogramPassed(ENicaCutUpdate::kTwoTrackUpdate, 0, 1);
  TH1D *mass_tof =
      (TH1D *)f->GetHistogramPassed(ENicaCutUpdate::kTwoTrackUpdate, 0, 0);
  TH1D *mass_tof_failed =
      (TH1D *)f->GetHistogramFailed(ENicaCutUpdate::kTwoTrackUpdate, 0, 0);
  mass_tof->Add(mass_tof_failed);

  TCanvas *c = new TCanvas();
  c->Divide(2, 1);
  c->cd(1);
  TVirtualPad *pad = gPad;
  pad->Divide(1, 2);
  pad->cd(1);
  armenteros_tpc->SetTitle("Arm TPC");
  armenteros_tpc->Draw("colz");
  pad->cd(2);
  armenteros_tof->SetTitle("Arm TOF");
  armenteros_tof->Draw("colz");
  c->cd(2);
  mass_tof->Sumw2();
  mass_tpc->Sumw2();
  mass_tof->SetFillStyle(0);
  mass_tpc->SetFillStyle(0);
  mass_tpc->SetLineColor(kRed);
  mass_tof->SetLineColor(kBlue);
  mass_tpc->SetMarkerColor(kRed);
  mass_tof->SetMarkerColor(kBlue);
  mass_tpc->SetMarkerStyle(kFullSquare);
  mass_tof->SetMarkerStyle(kFullSquare);
  mass_tpc->Draw();
  mass_tof->Draw("SAME");
  TLegend *l = new TLegend(0.7, 0.7, .95, 0.95);
  l->AddEntry(mass_tof, "with TOF", "LP");
  l->AddEntry(mass_tpc, "without TOF", "LP");
  l->Draw("SAME");
}

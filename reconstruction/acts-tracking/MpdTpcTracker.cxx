// This file is a part of the NICA project.
//
// Copyright (C) 2022 JINR

#include "MpdTpcTracker.h"

#include "MpdCodeTimer.h"
#include "MpdKalmanHit.h"
#include "MpdTpcKalmanTrack.h"
#include "MpdTpcHit.h"
#include "TpcPoint.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TMultiGraph.h"

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/Units.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <iostream>
#include <string>

static constexpr auto lenScalor = Acts::UnitConstants::cm  /* MPD  */ /
                                  Acts::UnitConstants::mm  /* Acts */;
static constexpr auto momScalor = Acts::UnitConstants::GeV /* MPD  */ /
                                  Acts::UnitConstants::GeV /* Acts */;

inline Mpd::Tpc::InputHitContainer convertTpcPoints(TClonesArray *tpcPoints) {
  const auto nTpcPoints = tpcPoints->GetEntriesFast();

  Mpd::Tpc::InputHitContainer hits;
  hits.reserve(nTpcPoints);

  for (Int_t i = 0; i < nTpcPoints; i++) {
    const auto *tpcPoint = static_cast<TpcPoint*>(tpcPoints->UncheckedAt(i));

    // FIXME: This is for debugging.
    // if (tpcPoint->GetTrackID() != 0) continue;

    hits.emplace_back(Mpd::Tpc::InputHit{
        tpcPoint->GetTrackID(),
        tpcPoint->GetDetectorID(),
        Acts::Vector3{
            tpcPoint->GetX() * lenScalor,
            tpcPoint->GetY() * lenScalor,
            tpcPoint->GetZ() * lenScalor
        },
        Acts::Vector3{
            tpcPoint->GetPx() * momScalor,
            tpcPoint->GetPy() * momScalor,
            tpcPoint->GetPz() * momScalor
        }
    });
  }

  return hits;
}

inline Mpd::Tpc::InputHitContainer convertTpcHits(TClonesArray *tpcHits) {
  const auto nTpcHits = tpcHits->GetEntriesFast();

  Mpd::Tpc::InputHitContainer hits;
  hits.reserve(nTpcHits);

  for (Int_t i = 0; i < nTpcHits; i++) {
    const auto *tpcHit = static_cast<MpdTpcHit*>(tpcHits->UncheckedAt(i));

    hits.emplace_back(Mpd::Tpc::InputHit{
        tpcHit->GetTrackID(),
        tpcHit->GetDetectorID(),
        Acts::Vector3{
            tpcHit->GetX() * lenScalor,
            tpcHit->GetY() * lenScalor,
            tpcHit->GetZ() * lenScalor
        },
        Acts::Vector3{0, 0, 0} // Fake momentum
    });
  }

  return hits;
}

inline TClonesArray *getArray(const char *name) {
  auto *manager = FairRootManager::Instance();
  auto *array = static_cast<TClonesArray*>(manager->GetObject(name));

  if (!array) {
    std::cout << "[MpdTpcTracker]: " << name << " not found" << std::endl;
  } else {
    std::cout << "[MpdTpcTracker]: " << name << " contains "
              << array->GetEntriesFast() << " entries" << std::endl;
  }

  return array;
}

InitStatus MpdTpcTracker::Init() {
  std::cout << "[MpdTpcTracker::Init]: Started" << std::endl;

  // Geometry must already be loaded (gGeoManager != nullptr).
  fRunner = std::make_unique<Mpd::Tpc::Runner>(
      "../../geometry/tpc_acts_tracking.json", // FIXME:
      Acts::Logging::DEBUG
  );

  std::cout << "[MpdTpcTracker::Init]: Finished" << std::endl;
  return kSUCCESS;
}

InitStatus MpdTpcTracker::ReInit() {
  std::cout << "[MpdTpcTracker::ReInit]: Do nothing" << std::endl;
  return kSUCCESS;
}

void plotOutputTracks(const int canvasX,
    const int canvasY,
    Mpd::Tpc::InputHitContainer hits,
    Mpd::Tpc::ProtoTrackContainer trajectories,
    int eventCounter) {
  TCanvas *canvas = new TCanvas("canvasName", "canvasName", canvasX, canvasY);
  TGraph *grInputHits = new TGraph();
  grInputHits->SetMarkerStyle(kFullDotMedium);

  size_t hitIndex = 0;
  for (auto hit : hits) {
    double hitX = hit.position[0];
    double hitY = hit.position[1];
    grInputHits->SetPoint(hitIndex++, hitX, hitY);
  }
  TMultiGraph *multiDotPlot = new TMultiGraph();
  multiDotPlot->Add(grInputHits, "P");
  multiDotPlot->Draw("A");

  std::string dotPlotFileNameI = "event_" + std::to_string(eventCounter) + "_i.png";
  canvas->Print(dotPlotFileNameI.c_str());
  std::cout << "File created: " << dotPlotFileNameI << std::endl;

  std::vector<TGraph*> outPlots;
  for (ActsExamples::ProtoTrack reconstructedTrack : trajectories) {
    TGraph *dotPlotReco = new TGraph();
    dotPlotReco->SetMarkerStyle(kFullDotMedium);
    dotPlotReco->SetLineWidth(3);
    dotPlotReco->SetLineColor(kRed);

    int hitCounter = 0;
    for (uint32_t hitIndex : reconstructedTrack) {
      double hitX = hits.at(hitIndex).position[0];
      double hitY = hits.at(hitIndex).position[1];
      dotPlotReco->SetPoint(hitCounter++, hitX, hitY);
    }
    multiDotPlot->Add(dotPlotReco, "PL");
    outPlots.push_back(dotPlotReco);
  }
  canvas->Clear();
  multiDotPlot->Draw("A");
  std::string dotPlotFileNameIO = "event_" + std::to_string(eventCounter) + "_io.png";
  canvas->Print(dotPlotFileNameIO.c_str());

  delete grInputHits;
  for (auto outPlot : outPlots) {
    delete outPlot;
  }
  delete multiDotPlot;
  delete canvas;

  std::cout << "File created: " << dotPlotFileNameIO << std::endl;
}

void MpdTpcTracker::Exec(Option_t *option) {
  std::cout << "[MpdTpcTracker::Exec]: Started" << std::endl;

  if (MpdCodeTimer::Active()) {
    MpdCodeTimer::Instance()->Start(Class()->GetName(), __FUNCTION__);
  }

  // Get the input points.
  fPoints = getArray(UseMcHits ? "TpcPoint" : "TpcRecPoint");
  assert(fPoints);

  // Convert the input points to the inner representation.
  auto hits = UseMcHits ? convertTpcPoints(fPoints)
                        : convertTpcHits(fPoints);

  // Run the track finding algorithm.
  const auto &trajectories = fRunner->execute(hits);

  // Get the track recognition statistics (for debugging).
  auto statistics = fRunner->getStatistics();
  auto nTracks = fRunner->getTracksNumber();

  // Build histagrams.
  static int staticEventCounter = 0;
  staticEventCounter++;

  std::string statisticsFileName = "stat_" + std::to_string(staticEventCounter) + ".root";

  TFile rootFileWithStat = TFile (statisticsFileName.c_str(), "RECREATE");
  std::string histoName = "The number of real tracks: " + std::to_string(nTracks);
  TH1* h1 = new TH1I("statistics", histoName.c_str(), 10, 0.0, 100.0);

  for (auto const& [key, val]: statistics) {
    h1->Fill(val.quality);
  }
  h1->Write();
  rootFileWithStat.Close();

  // Plot the output tracks.
  const int canvasX = 3000;
  const int canvasY = 3000;
  plotOutputTracks(canvasX, canvasY, hits, trajectories, staticEventCounter);

  // Convert the output tracks.
  (void)trajectories; // FIXME:

//  fKalmanHits = getArray("MpdKalmanHit");
//  fKalmanTracks = getArray("MpdTpcKalmanTrack");

  if (MpdCodeTimer::Active()) {
    MpdCodeTimer::Instance()->Stop(Class()->GetName(), __FUNCTION__);
  }

  std::cout << "[MpdTpcTracker::Exec]: Finished" << std::endl;
}

void MpdTpcTracker::Finish() {
  std::cout << "[MpdTpcTracker::Finish]: Do nothing" << std::endl;
}

ClassImp(MpdTpcTracker);

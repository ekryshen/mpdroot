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

void plotOutputTracks(const int canvasX,
                      const int canvasY,
                      const Mpd::Tpc::InputHitContainer &hits,
                      const Mpd::Tpc::ProtoTrackContainer &trajectories,
                      const Mpd::Tpc::SpacePointContainer &spacePoints,
                      int eventCounter);

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

void MpdTpcTracker::Exec(Option_t *option) {
  // For naming files during debug.
  static int eventCounter = 0;

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
  const auto &config = fRunner->config();
  const auto &context = fRunner->context();

  fRunner->execute(hits);

  const auto &trajectories =
      context.eventStore.get<Mpd::Tpc::ProtoTrackContainer>(
          config.trackFinding.outputTrackCandidates);

  const auto &spacePoints =
      context.eventStore.get<Mpd::Tpc::SpacePointContainer>(
          config.spacePointMaking.outputSpacePoints);

  // Get the track recognition statistics (for debugging).
  auto statistics = fRunner->getStatistics();
  auto nTracks = fRunner->getTracksNumber();

  // Build histagrams.
  std::string statisticsFileName = "stat_" + std::to_string(eventCounter) + ".root";

  TFile rootFileWithStat(statisticsFileName.c_str(), "RECREATE");
  std::string histoName = "The number of real tracks: " + std::to_string(nTracks);
  TH1 *h1 = new TH1I("statistics", histoName.c_str(), 10, 0.0, 100.0);

  for (auto const &[key, val]: statistics) {
    h1->Fill(val.quality);
  }
  h1->Write();
  rootFileWithStat.Close();

  // Plot the output tracks.
  plotOutputTracks(6000, 6000, hits, trajectories, spacePoints, eventCounter);

  // Convert the output tracks.
//  fKalmanHits = getArray("MpdKalmanHit");
//  fKalmanTracks = getArray("MpdTpcKalmanTrack");

  if (MpdCodeTimer::Active()) {
    MpdCodeTimer::Instance()->Stop(Class()->GetName(), __FUNCTION__);
  }

  std::cout << "[MpdTpcTracker::Exec]: Finished" << std::endl;
  eventCounter++;
}

void MpdTpcTracker::Finish() {
  std::cout << "[MpdTpcTracker::Finish]: Do nothing" << std::endl;
}

void plotOutputTracks(const int canvasX,
                      const int canvasY,
                      const Mpd::Tpc::InputHitContainer &hits,
                      const Mpd::Tpc::ProtoTrackContainer &trajectories,
                      const Mpd::Tpc::SpacePointContainer &spacePoints,
                      int eventCounter) {
  TCanvas *canvas = new TCanvas("canvasName", "canvasName", canvasX, canvasY);
  TMultiGraph *multiGraph = new TMultiGraph();

// SpacePoints graph
  TGraph *spacePointsGraph = new TGraph();
  spacePointsGraph->SetMarkerStyle(kFullDotMedium);
  spacePointsGraph->SetMarkerColor(kOrange);

  size_t pIndex = 0;
  for (auto spacePoint : spacePoints) {
    double x = spacePoint.x();
    double y = spacePoint.y();
    spacePointsGraph->SetPoint(pIndex++, x, y);
  }
  multiGraph->Add(spacePointsGraph, "P");

// input hits graph
  TGraph *inputHitsGraph = new TGraph();
  inputHitsGraph->SetMarkerStyle(kFullDotMedium);

  pIndex = 0;
  for (auto hit : hits) {
    double x = hit.position[0];
    double y = hit.position[1];
    inputHitsGraph->SetPoint(pIndex++, x, y);
  }
  multiGraph->Add(inputHitsGraph, "P");

// reconstructed tracks graph
  std::vector<TGraph*> graphs;
  for (ActsExamples::ProtoTrack reconstructedTrack : trajectories) {
    TGraph *outTrajectoryGraph = new TGraph();
    outTrajectoryGraph->SetMarkerStyle(kFullDotMedium);
    outTrajectoryGraph->SetLineWidth(3);
    outTrajectoryGraph->SetLineColor(kRed);

    int hitCounter = 0;
    for (uint32_t hitIndex : reconstructedTrack) {
      double x = hits.at(hitIndex).position[0];
      double y = hits.at(hitIndex).position[1];
      outTrajectoryGraph->SetPoint(hitCounter++, x, y);
    }
    multiGraph->Add(outTrajectoryGraph, "PL");
    graphs.push_back(outTrajectoryGraph);
  }
  multiGraph->Draw("A");
  std::string dotPlotFileName = "event_" + std::to_string(eventCounter) + ".png";
  canvas->Print(dotPlotFileName.c_str());

  delete inputHitsGraph;
  for (auto graph : graphs) {
    delete graph;
  }
  delete multiGraph;
  delete canvas;

  std::cout << "File created: " << dotPlotFileName << std::endl;
}

ClassImp(MpdTpcTracker);

// This file is a part of the NICA project.
//
// Copyright (C) 2022 JINR

#include "MpdTpcTracker.h"
#include "PlotUtils.h"

#include "MpdCodeTimer.h"
#include "MpdKalmanHit.h"
#include "MpdTpcKalmanTrack.h"
#include "MpdTpcHit.h"
#include "TpcPoint.h"

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/Units.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <iostream>

//===----------------------------------------------------------------------===//
// Utilities
//===----------------------------------------------------------------------===//

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

//===----------------------------------------------------------------------===//
// Conversion
//===----------------------------------------------------------------------===//

static constexpr auto lenScalor = Acts::UnitConstants::cm  /* MPD  */ /
                                  Acts::UnitConstants::mm  /* Acts */;
static constexpr auto momScalor = Acts::UnitConstants::GeV /* MPD  */ /
                                  Acts::UnitConstants::GeV /* Acts */;

/// Converts MC points to the internal representation.
inline Mpd::Tpc::InputHitContainer convertTpcPoints(TClonesArray *tpcPoints) {
  const auto nTpcPoints = tpcPoints->GetEntriesFast();

  Mpd::Tpc::InputHitContainer hits;
  hits.reserve(nTpcPoints);

  for (Int_t i = 0; i < nTpcPoints; i++) {
    const auto *tpcPoint = static_cast<TpcPoint*>(tpcPoints->UncheckedAt(i));

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

/// Converts TPC hits to the internal representation.
inline Mpd::Tpc::InputHitContainer convertTpcHits(TClonesArray *tpcHits,
                                                  TClonesArray *tpcPoints) {
  // Connect information on momentum (for estimating the algorithm efficiency).
  std::unordered_map<int, Acts::Vector3> momentum;
  const auto nTpcPoints = tpcPoints->GetEntriesFast();

  for (Int_t i = 0; i < nTpcPoints; i++) {
    const auto *tpcPoint = static_cast<TpcPoint*>(tpcPoints->UncheckedAt(i));
    momentum[tpcPoint->GetTrackID()] = Acts::Vector3{
                                         tpcPoint->GetPx() * momScalor,
                                         tpcPoint->GetPx() * momScalor,
                                         tpcPoint->GetPz() * momScalor
                                       };
  }

  const auto nTpcHits = tpcHits->GetEntriesFast();

  Mpd::Tpc::InputHitContainer hits;
  hits.reserve(nTpcHits);

  for (Int_t i = 0; i < nTpcHits; i++) {
    const auto *tpcHit = static_cast<MpdTpcHit*>(tpcHits->UncheckedAt(i));

    const auto iter = momentum.find(tpcHit->GetTrackID());
    assert(iter != momentum.end());

    hits.emplace_back(Mpd::Tpc::InputHit{
        tpcHit->GetTrackID(),
        tpcHit->GetDetectorID(),
        Acts::Vector3{
            tpcHit->GetX() * lenScalor,
            tpcHit->GetY() * lenScalor,
            tpcHit->GetZ() * lenScalor
        },
        iter->second // For statistics only.
    });
  }

  return hits;
}

//===----------------------------------------------------------------------===//
// Init / ReInit Tasks
//===----------------------------------------------------------------------===//

InitStatus MpdTpcTracker::Init() {
  std::cout << "[MpdTpcTracker::Init]: Started" << std::endl;

  // Geometry must already be loaded (gGeoManager != nullptr).
  fRunner = std::make_unique<Mpd::Tpc::Runner>(
      "../../geometry/tpc_acts_tracking.json", // FIXME:
      Acts::Logging::DEBUG
  );

  // FIXME: Should be a field.
  auto *fTracks = new TClonesArray("MpdTpcKalmanTrack");
  FairRootManager::Instance()->Register("TpcKalmanTrack", "MpdKalmanTrack", fTracks, kTRUE);

  std::cout << "[MpdTpcTracker::Init]: Finished" << std::endl;
  return kSUCCESS;
}

InitStatus MpdTpcTracker::ReInit() {
  std::cout << "[MpdTpcTracker::ReInit]: Do nothing" << std::endl;
  return kSUCCESS;
}

//===----------------------------------------------------------------------===//
// Exec Task
//===----------------------------------------------------------------------===//

void MpdTpcTracker::Exec(Option_t *option) {
  // For naming files during debug.
  static int eventCounter = -1;
  eventCounter++;

  std::cout << "[MpdTpcTracker::Exec]: Started" << std::endl;

  if (MpdCodeTimer::Active()) {
    MpdCodeTimer::Instance()->Start(Class()->GetName(), __FUNCTION__);
  }

  // Get the MC points.
  fPoints = getArray("TpcPoint");
  assert(fPoints);
  // Get the TPC hits.
  fHits = getArray("TpcRecPoint");
  assert(fHits);

  // Convert the input points to the internal representation.
  auto hits = UseMcHits ? convertTpcPoints(fPoints)
                        : convertTpcHits(fHits, fPoints);

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
  buildHistograms(statistics, nTracks, eventCounter);
  drawQualityOnP(hits, trajectories, eventCounter);

  // Plot the output tracks.
  std::shared_ptr<const Acts::TrackingGeometry> geometry =
      config.detector->getGeometry();

  plotOutputTracks(6000, 6000, geometry, spacePoints, hits, trajectories, eventCounter, false);

  // Convert the output tracks.
//  fKalmanHits = getArray("MpdKalmanHit");
//  fKalmanTracks = getArray("MpdTpcKalmanTrack");

  if (MpdCodeTimer::Active()) {
    MpdCodeTimer::Instance()->Stop(Class()->GetName(), __FUNCTION__);
  }

  std::cout << "[MpdTpcTracker::Exec]: Finished" << std::endl;
}

//===----------------------------------------------------------------------===//
// Finish Task
//===----------------------------------------------------------------------===//

void MpdTpcTracker::Finish() {
  std::cout << "[MpdTpcTracker::Finish]: Do nothing" << std::endl;
}

ClassImp(MpdTpcTracker);

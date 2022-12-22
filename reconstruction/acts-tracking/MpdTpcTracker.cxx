// This file is a part of the NICA project.
//
// Copyright (C) 2022 JINR

#include "MpdTpcTracker.h"
#include "PlotUtils.h"

#include "MpdCodeTimer.h"
#include "MpdKalmanHit.h"
#include "MpdTpcHit.h"
#include "MpdTpcTrack.h"
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
  assert(tpcPoints);
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
  assert(tpcHits && tpcPoints);

  // Connect information on momentum (for estimating the algorithm efficiency).
  std::unordered_map<int, Acts::Vector3> momentum;
  const auto nTpcPoints = tpcPoints->GetEntriesFast();

  for (Int_t i = 0; i < nTpcPoints; i++) {
    const auto *tpcPoint = static_cast<TpcPoint*>(tpcPoints->UncheckedAt(i));
    momentum[tpcPoint->GetTrackID()] = Acts::Vector3{
                                           tpcPoint->GetPx() * momScalor,
                                           tpcPoint->GetPy() * momScalor,
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

/// Converts a hit into the MpdRoot representation.
inline MpdTpcHit *convertHit(const BaseTpcSectorGeo &secGeo,
                             TClonesArray *hits, int ihit) {
  return static_cast<MpdTpcHit*>(hits->UncheckedAt(ihit));
}

/// Converts a track into the MpdRoot representation.
inline void convertTrack(const BaseTpcSectorGeo &secGeo,
                         TClonesArray *hits,
                         MpdTpcTrack *track,
                         const Mpd::Tpc::ProtoTrack &trajectory) {
  for (const auto ihit : trajectory) {
    track->GetHits()->Add(convertHit(secGeo, hits, ihit));
  }
}

/// Converts tracks into the MpdRoot representation.
inline void convertTracks(const BaseTpcSectorGeo &secGeo,
                          TClonesArray *hits,
                          TClonesArray *tracks,
                          const Mpd::Tpc::ProtoTrackContainer &trajectories) {
  tracks->Clear();
  tracks->Expand(trajectories.size());

  Int_t nTracks = 0;
  for (const auto &trajectory : trajectories) {
    auto *track = new ((*tracks)[nTracks++]) MpdTpcTrack(trajectory.size());
    convertTrack(secGeo, hits, track, trajectory);
  }
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

  fTracks = new TClonesArray("MpdTpcTrack");
  FairRootManager::Instance()->Register(
      "TpcKalmanTrack", "MpdKalmanTrack", fTracks, kTRUE);

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

  // Get the MC points and TPC hits.
  fPoints = getArray("TpcPoint");
  fHits = getArray("TpcRecPoint");

  // Convert the input points to the internal representation.
  auto hits = UseMcHits ? convertTpcPoints(fPoints)
                        : convertTpcHits(fHits, fPoints);

  // Run the track finding algorithm.
  const auto &config = fRunner->config();
  const auto &context = fRunner->context();

  fRunner->execute(hits);

  // Convert the found track to to the MpdRoot representation.
  const auto &trajectories =
      context.eventStore.get<Mpd::Tpc::ProtoTrackContainer>(
          config.trackFinding.outputTrackCandidates);

  convertTracks(fSecGeo, fHits, fTracks, trajectories);

  // Plot graphs if required.
  if constexpr(PlotGraphs) {
    const auto &spacePoints =
        context.eventStore.get<Mpd::Tpc::SpacePointContainer>(
            config.spacePointMaking.outputSpacePoints);

    // Get the track recognition statistics.
    auto statistics = fRunner->getStatistics();
    auto nTracks = fRunner->getTracksNumber();

    // Build histograms.
    buildHistograms(statistics, nTracks, eventCounter);
    plotQualityOnP(hits, trajectories, eventCounter);

    // Plot the output tracks.
    std::shared_ptr<const Acts::TrackingGeometry> geometry =
        config.detector->getGeometry();

    const Int_t lineWidth = 3;
    const Bool_t multicoloured = false;
    plotOutputTracks(6000, 6000, geometry, spacePoints, hits,
                     trajectories, eventCounter, multicoloured,
                     lineWidth, CoordSystem::XY);
  }

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

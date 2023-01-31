// This file is a part of the NICA project.
//
// Copyright (C) 2022 JINR

#include "MpdTpcTracker.h"
#include "PlotUtils.h"

#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include "MpdCodeTimer.h"
#include "MpdKalmanHit.h"
#include "MpdMCTrack.h"
#include "MpdTpcHit.h"
#include "MpdTpcTrack.h"
#include "TpcPoint.h"

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/Units.hpp>
#include <Acts/Utilities/Logger.hpp>
#include <ActsFatras/EventData/Barcode.hpp>
#include <ActsFatras/EventData/Particle.hpp>

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

    hits.emplace_back(Mpd::Tpc::InputHit{
        tpcHit->GetTrackID(),
        tpcHit->GetDetectorID(),
        Acts::Vector3{
            tpcHit->GetX() * lenScalor,
            tpcHit->GetY() * lenScalor,
            tpcHit->GetZ() * lenScalor
        },
        (iter != momentum.end()) ? iter->second : Acts::Vector3{0., 0., 0.} // FIXME:
    });
  }

  return hits;
}

/// Converts a hit into the MpdRoot representation.
inline MpdTpcHit *convertHit(TClonesArray *hits, int ihit) {
  return static_cast<MpdTpcHit*>(hits->UncheckedAt(ihit));
}

/// Converts a track into the MpdRoot representation.
inline void convertTrack(const BaseTpcSectorGeo &secGeo,
                         TClonesArray *hits,
                         MpdTpcTrack *track,
                         const Mpd::Tpc::ProtoTrack &trajectory) {
  for (const auto ihit : trajectory) {
    track->GetHits()->Add(convertHit(hits, ihit));
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
// Get input particles
//===----------------------------------------------------------------------===//

ActsExamples::SimParticleContainer getInputParticles(TClonesArray *mcTracks) {
  ActsExamples::SimParticleContainer particles;

  Int_t nMC = mcTracks->GetEntriesFast();
  for (Int_t i = 0; i < nMC; i++) {
    auto track = static_cast<MpdMCTrack*>(mcTracks->UncheckedAt(i));
    Int_t pdgCode  = track->GetPdgCode();
    Double_t x0 = track->GetStartX() * lenScalor;
    Double_t y0 = track->GetStartY() * lenScalor;
    Double_t z0 = track->GetStartZ() * lenScalor;
    Double_t t0 = track->GetStartT();
    Double_t px = track->GetPx() * momScalor;
    Double_t py = track->GetPy() * momScalor;
    Double_t pz = track->GetPz() * momScalor;
    Double_t p  = track->GetP()  * momScalor;

    ActsFatras::Barcode barCode = ActsFatras::Barcode(). // FIXME
        setVertexPrimary(1).
        setVertexSecondary(2).
        setParticle(3);

    ActsFatras::Particle particle =
        ActsFatras::Particle(barCode, Acts::PdgParticle(pdgCode)).
        setPosition4(x0, y0, z0, t0).
        setAbsoluteMomentum(p).
        setDirection(px, py, pz);

    particles.insert(std::move(particle));
  }
  return particles;
}

//===----------------------------------------------------------------------===//
// Get multimap hits -> particles
//===----------------------------------------------------------------------===//
ActsExamples::IndexMultimap<ActsFatras::Barcode> getMultimapHitToPart(
    TClonesArray *hits) {
  ActsExamples::IndexMultimap<ActsFatras::Barcode> res;
  Int_t nHits = hits->GetEntriesFast();
  res.reserve(nHits);
  for (Int_t i = 0; i < nHits; i++) {
     ActsFatras::Barcode barCode = ActsFatras::Barcode(). // FIXME
        setVertexPrimary(1).
        setVertexSecondary(2).
        setParticle(3);
     res.emplace_hint(res.end(), i, barCode);
  }
  return res;
}

//===----------------------------------------------------------------------===//
// Init / ReInit Tasks
//===----------------------------------------------------------------------===//

InitStatus MpdTpcTracker::Init() {
  std::cout << "[MpdTpcTracker::Init]: Started" << std::endl;

  // Geometry must already be loaded (gGeoManager != nullptr).
  fRunner = std::make_unique<Mpd::Tpc::Runner>(
      fSecGeo,
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

  auto mcTracks = getArray("MCTrack");
  auto inputParticles = getInputParticles(mcTracks);
  auto mapHitsToParticles = getMultimapHitToPart(fHits);

  // Run the track finding algorithm.
  const auto &config = fRunner->config();

  ActsExamples::WhiteBoard whiteBoard;
  ActsExamples::AlgorithmContext context(0, eventCounter, whiteBoard);
  fRunner->execute(hits, inputParticles, mapHitsToParticles, context);

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
    auto statistics = fRunner->getStatistics(context);
    auto nTracks = fRunner->getTracksNumber(context);

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

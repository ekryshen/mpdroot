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

/// Dump hits to file in Acts units.
void dumpData(TClonesArray *tpcHits,
              Int_t eventNumber,
              std::string outPath) {
  auto fname = outPath + "/event_" + std::to_string(eventNumber) + "_hits.txt";
  std::ofstream fout(fname);

  const auto nTpcHits = tpcHits->GetEntriesFast();

  for (size_t i = 0; i < nTpcHits; i++) {
    auto *tpcHit = static_cast<MpdTpcHit*>(tpcHits->UncheckedAt(i));
    fout << tpcHit->GetX() * lenScalor << ", " <<
            tpcHit->GetY() * lenScalor << ", " <<
            tpcHit->GetZ() * lenScalor << ", ";

    std::vector <std::pair<int, float>> pairs = tpcHit->GetTrackIDs();
    for (auto trackIdandQ : pairs) {
      auto trackID = trackIdandQ.first;
      auto q       = trackIdandQ.second;
      fout << trackID << ", " <<
              q << std::endl;
    }
  }
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
    auto *tpcHit = static_cast<MpdTpcHit*>(tpcHits->UncheckedAt(i));

    std::vector<std::pair<int, float>> trackIDs = tpcHit->GetTrackIDs();

    // Get the first MC track, correspoinding to hit.
    Int_t trackId = trackIDs[0].first;

    const auto mFind = momentum.find(trackId);
    Acts::Vector3 mom;
    if (mFind == momentum.end()) {
      mom = Acts::Vector3{0., 0., 0.};
      std::cout << "[MpdTpcTracker]: WARNING: can't find MC track " <<
          "for MpdTpcHit with index " << i << std::endl;
    } else {
      mom = mFind->second;
    }
    hits.emplace_back(Mpd::Tpc::InputHit{
        trackId,
        tpcHit->GetDetectorID(),
        Acts::Vector3{
            tpcHit->GetX() * lenScalor,
            tpcHit->GetY() * lenScalor,
            tpcHit->GetZ() * lenScalor
        },
        mom
    });
  }
  return hits;
}

/// Converts a hit into the MpdRoot representation.
inline MpdTpcHit *convertHit(TClonesArray *hits, int ihit) {
  return static_cast<MpdTpcHit*>(hits->UncheckedAt(ihit));
}

/// Converts a track into the MpdRoot representation.
inline void convertTrack(TClonesArray *hits,
                         MpdTpcTrack *track,
                         const Mpd::Tpc::ProtoTrack &trajectory) {
  for (const auto ihit : trajectory) {
    track->GetHits()->Add(convertHit(hits, ihit));
  }
}

/// Converts tracks into the MpdRoot representation.
inline void convertTracks(TClonesArray *hits,
                          TClonesArray *tracks,
                          const Mpd::Tpc::ProtoTrackContainer &trajectories) {
  tracks->Clear();
  tracks->Expand(trajectories.size());

  Int_t nTracks = 0;
  for (const auto &trajectory : trajectories) {
    auto *track = new ((*tracks)[nTracks++]) MpdTpcTrack(trajectory.size());
    convertTrack(hits, track, trajectory);
  }
}

//===----------------------------------------------------------------------===//
// Collect informations for Acts statistics
//===----------------------------------------------------------------------===//

/// Create map MC track -> ints for creating Act's Barcode
using mcTrackToBarcodeInts =
    std::map<Int_t, std::tuple<size_t, size_t, size_t, size_t, size_t>>;

mcTrackToBarcodeInts createBarcodesMap(TClonesArray *mcTracks) {

  mcTrackToBarcodeInts mcTrackToBarcode;

  Int_t nMC = mcTracks->GetEntriesFast();
  size_t primary     = 1;
  size_t secondary   = 0;
  size_t curParticle = 0;

  // Loop over primary particles.
  for (Int_t iMC = 0; iMC < nMC; iMC++) {
    auto track = static_cast<MpdMCTrack*>(mcTracks->UncheckedAt(iMC));

    Int_t motherId = track->GetMotherId();
    if (motherId != -1) {
      continue;
    }
    size_t particle   = ++curParticle;
    size_t generation = 0;
    size_t subpart    = 0;
    mcTrackToBarcode[iMC] = {primary, secondary, particle, generation, subpart};
  }

  // Loop over secondary particles.
  std::map<std::tuple<size_t, size_t>, size_t> partGenToSubpart;
  Int_t iMC;
  std::string msgPrefix = "[MpdTpcTracker]: createBarcodesMap(): ERROR: ";
  for (iMC = 0; iMC < nMC; iMC++) {
    auto track = static_cast<MpdMCTrack*>(mcTracks->UncheckedAt(iMC));
    Int_t motherId = track->GetMotherId();
    Int_t motherIdBak = motherId;
    if (motherId == -1) {
      continue;
    }
    MpdMCTrack *curTrack;
    size_t generation = 0;
    while (motherId != -1) {
      generation++;
      motherIdBak = motherId;
      curTrack = static_cast<MpdMCTrack*>(mcTracks->UncheckedAt(motherId));
      motherId = curTrack->GetMotherId();
    }
    if ((generation == 0) || (motherIdBak == -1)) {
      std::cout << msgPrefix <<
          "can't find top track for MC track with index " << iMC << std::endl;
      continue;
    }
    auto findBarcode = mcTrackToBarcode.find(motherIdBak);
    if (findBarcode == mcTrackToBarcode.end()) {
       std::cout << msgPrefix <<
          "can't process MC track with index " << iMC << std::endl;
       continue;
    }
    auto fiveNums = findBarcode->second;
    auto particle = std::get<2>(fiveNums);

    size_t subpart = partGenToSubpart[{particle, generation}];
    mcTrackToBarcode[iMC] = {primary, secondary, particle, generation, subpart};
    partGenToSubpart[{particle, generation}] += 1;
  }
  return mcTrackToBarcode;
}

/// Get input particles
ActsExamples::SimParticleContainer getInputParticles(
    TClonesArray *mcTracks,
    const mcTrackToBarcodeInts &mcTrackToBarcode) {

  ActsExamples::SimParticleContainer particles;

  Int_t nMC = mcTracks->GetEntriesFast();
  for (Int_t i = 0; i < nMC; i++) {
    auto track = static_cast<MpdMCTrack*>(mcTracks->UncheckedAt(i));

    if (auto search = mcTrackToBarcode.find(i);
             search == mcTrackToBarcode.end()) {
      std::cout << "[MpdTpcTracker.cxx]: getInputParticles(): ERROR: " <<
          "can't find Barcode for MC track with index " << i << std::endl;
      continue;
    }

    const auto [pri, sec, part, gen, sub] =  mcTrackToBarcode.at(i);
    const auto barcode = ActsFatras::Barcode().
        setVertexPrimary(pri).
        setVertexSecondary(sec).
        setParticle(part).
        setGeneration(gen).
        setSubParticle(sub);

    auto pdgCode = track->GetPdgCode();
    auto x0      = track->GetStartX() * lenScalor;
    auto y0      = track->GetStartY() * lenScalor;
    auto z0      = track->GetStartZ() * lenScalor;
    auto t0      = track->GetStartT();
    auto px      = track->GetPx() * momScalor;
    auto py      = track->GetPy() * momScalor;
    auto pz      = track->GetPz() * momScalor;
    auto p       = track->GetP()  * momScalor;

    auto particle = ActsFatras::Particle(barcode, Acts::PdgParticle(pdgCode)).
        setPosition4(x0, y0, z0, t0).
        setAbsoluteMomentum(p).
        setDirection(px, py, pz);

    particles.insert(std::move(particle));
  }
  return particles;
}

/// Create multimap TpcPoint *hits idx -> particles
ActsExamples::IndexMultimap<ActsFatras::Barcode> createTpcPointsToParticlesMap(
    TClonesArray *tpcPoints,
    const mcTrackToBarcodeInts &mcTrackToBarcode) {

  ActsExamples::IndexMultimap<ActsFatras::Barcode> mapHitsToParticles;
  Int_t n = tpcPoints->GetEntriesFast();

  for (Int_t i = 0; i < n; i++) {

    auto *tpcPoint = static_cast<TpcPoint*>(tpcPoints->UncheckedAt(i));
    auto trackId = tpcPoint->GetTrackID();

    if (auto search = mcTrackToBarcode.find(trackId);
             search == mcTrackToBarcode.end()) {
      std::cout << "[MpdTpcTracker.cxx]: createTpcPointsToParticlesMap():" <<
          " ERROR: can't find Barcode for MC track with index " <<
          trackId << std::endl;
      continue;
    }

    const auto [pri, sec, part, gen, sub] =  mcTrackToBarcode.at(trackId);
    const auto barcode = ActsFatras::Barcode().
        setVertexPrimary(pri).
        setVertexSecondary(sec).
        setParticle(part).
        setGeneration(gen).
        setSubParticle(sub);
    mapHitsToParticles.emplace(i, barcode);
  }
  return mapHitsToParticles;
}

/// Create multimap MpdTpcHit *hits idx -> particles
ActsExamples::IndexMultimap<ActsFatras::Barcode> createTpcHitsToParticlesMap(
    TClonesArray *tpcHits,
    const mcTrackToBarcodeInts &mcTrackToBarcode) {

  ActsExamples::IndexMultimap<ActsFatras::Barcode> mapHitsToParticles;
  Int_t n = tpcHits->GetEntriesFast();

  for (Int_t i = 0; i < n; i++) {

    auto *tpcHit = static_cast<MpdTpcHit*>(tpcHits->UncheckedAt(i));
    std::vector <std::pair<int, float>> trackIDs = tpcHit->GetTrackIDs();

    for (auto particle : trackIDs) {
      auto trackId = particle.first;

      if (auto search = mcTrackToBarcode.find(trackId);
               search == mcTrackToBarcode.end()) {
        std::cout << "[MpdTpcTracker.cxx]: createTpcHitsToParticlesMap():" <<
            " ERROR: can't find Barcode for MC track with index " <<
            trackId << std::endl;
        continue;
      }

      const auto [pri, sec, part, gen, sub] =  mcTrackToBarcode.at(trackId);
      const auto barcode = ActsFatras::Barcode().
          setVertexPrimary(pri).
          setVertexSecondary(sec).
          setParticle(part).
          setGeneration(gen).
          setSubParticle(sub);
      mapHitsToParticles.emplace(i, barcode);
    }
  }
  return mapHitsToParticles;
}

//===----------------------------------------------------------------------===//
// Init / ReInit Tasks
//===----------------------------------------------------------------------===//

InitStatus MpdTpcTracker::Init() {
  std::cout << "[MpdTpcTracker::Init]: Started" << std::endl;

  auto binPath = getenv("VMCWORKDIR");
  assert(binPath && "VMCWORKDIR environment variable is empty!");

  auto level = Acts::Logging::DEBUG;
  // Geometry must already be loaded (gGeoManager != nullptr).
  fRunner = std::make_unique<Mpd::Tpc::Runner>(
      fSecGeo,
      std::string(binPath) + "/geometry/tpc_acts_tracking.json",
      level
  );

  auto perfCfg = fRunner->config().perfWriterCfg(fOutPath);
  fPerfWriter = new ActsExamples::CKFPerformanceWriter(perfCfg, level);

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
  static Int_t eventCounter = -1;
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

  // Get the information necessary to collect Acts statistics.
  auto mcTracks = getArray("MCTrack");
  auto mcTrackToBarcode = createBarcodesMap(mcTracks);
  auto inputParticles = getInputParticles(mcTracks, mcTrackToBarcode);
  auto mapHitsToParticles = UseMcHits
     ? createTpcPointsToParticlesMap(fPoints, mcTrackToBarcode)
     : createTpcHitsToParticlesMap  (fHits,   mcTrackToBarcode);

  // Run the track finding algorithm.
  const auto &config = fRunner->config();

  ActsExamples::WhiteBoard whiteBoard;
  ActsExamples::AlgorithmContext context(0, eventCounter, whiteBoard);
  fRunner->execute(hits, inputParticles, mapHitsToParticles,
                   fPerfWriter, context, fOutPath);

  // Dump hits to file.
  context.eventStore.add(config.DumpDataID, config.DumpData);
  if (config.DumpData) {
    dumpData(fHits, eventCounter, fOutPath);
  }

  // Convert the found track to to the MpdRoot representation.
  const auto &trajectories =
      context.eventStore.get<Mpd::Tpc::ProtoTrackContainer>(
          config.trackFinding.outputTrackCandidates);

  convertTracks(fHits, fTracks, trajectories);

  // Plot graphs if required.
  if constexpr(PlotGraphs) {
    const auto &spacePoints =
        context.eventStore.get<Mpd::Tpc::SpacePointContainer>(
            config.spacePointMaking.outputSpacePoints);

    // Get the track recognition statistics.
    auto statistics = fRunner->getStatistics(context);
    auto nTracks = fRunner->getTracksNumber(context);

    // Build histograms.
    buildHistograms(statistics, nTracks, eventCounter, fOutPath);
    plotQualityOnP(hits, trajectories, eventCounter, fOutPath);

    // Plot the output tracks.
    std::shared_ptr<const Acts::TrackingGeometry> geometry =
        config.detector->getGeometry();

    auto lineWidth = 3;
    auto color     = true;

    auto MChits = convertTpcPoints(fPoints);
    plotRealTracks(6000, 6000, geometry, MChits, eventCounter,
        fOutPath, color, lineWidth, Projection::XY, "_MC_input_XY");

    plotRealTracks(6000, 6000, geometry, hits, eventCounter,
        fOutPath, color, lineWidth, Projection::XY, "_hits_input_XY");

    plotOutputTracks(6000, 6000, geometry, spacePoints, hits,
        trajectories, eventCounter, fOutPath,
        color, lineWidth, Projection::XY);
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

// This file is a part of the NICA project.
//
// Copyright (C) 2022 JINR

#include "MpdTpcContext.h"
#include "MpdTpcEventStorage.h"
#include "MpdTpcTrackSeeding.h"

#include <Acts/Seeding/BinFinder.hpp>
#include <Acts/Seeding/BinnedSPGroup.hpp>
#include <Acts/Seeding/Seed.hpp>
#include <Acts/Seeding/SeedFilter.hpp>
#include <Acts/Seeding/Seedfinder.hpp>
#include <Acts/Surfaces/Surface.hpp>

#include <cassert>
#include <vector>

namespace Mpd::Tpc {

TrackSeeding::TrackSeeding(Config config, Acts::Logging::Level level):
  Algorithm("TrackSeeding", level),
  m_config(std::move(config)) {

  assert(!m_config.inputSpacePoints.empty()
    && "Missing space point input collections");
  assert(!m_config.outputProtoTracks.empty()
    && "Missing proto tracks output collection");
  assert(!m_config.outputSeeds.empty()
    && "Missing seeds output collection");
  assert((m_config.gridConfig.rMax == m_config.seedFinderConfig.rMax)
    && "Inconsistent config rMax");
  assert((m_config.seedFilterConfig.deltaRMin ==
          m_config.seedFinderConfig.deltaRMin)
    && "Inconsistent config deltaRMin");
  assert((m_config.gridConfig.deltaRMax ==
          m_config.seedFinderConfig.deltaRMax)
    && "Inconsistent config deltaRMax");
  assert((m_config.gridConfig.zMin == m_config.seedFinderConfig.zMin)
    && "Inconsistent config zMin");
  assert((m_config.gridConfig.zMax == m_config.seedFinderConfig.zMax)
    && "Inconsistent config zMax");
  assert((m_config.seedFilterConfig.maxSeedsPerSpM ==
          m_config.seedFinderConfig.maxSeedsPerSpM)
    && "Inconsistent config maxSeedsPerSpM");
  assert((m_config.gridConfig.cotThetaMax ==
          m_config.seedFinderConfig.cotThetaMax)
    && "Inconsistent config cotThetaMax");
  assert((m_config.gridConfig.minPt == m_config.seedFinderConfig.minPt)
    && "Inconsistent config minPt");
  assert((m_config.gridConfig.bFieldInZ == m_config.seedFinderConfig.bFieldInZ) 
    && "Inconsistent config bFieldInZ");

  assert(((m_config.gridConfig.zBinEdges.size() - 1 == m_config.zBinNeighborsTop.size())
        || m_config.zBinNeighborsTop.empty())
    && "Inconsistent config zBinNeighborsTop");

  assert(((m_config.gridConfig.zBinEdges.size() - 1 == m_config.zBinNeighborsBottom.size())
        || m_config.zBinNeighborsBottom.empty())
    && "Inconsistent config zBinNeighborsBottom");

  if (!m_config.seedFinderConfig.zBinsCustomLooping.empty()) {
    for (size_t i = 1; i != m_config.gridConfig.zBinEdges.size(); i++) {
      auto it = std::find(m_config.seedFinderConfig.zBinsCustomLooping.begin(),
                          m_config.seedFinderConfig.zBinsCustomLooping.end(), i); 
      assert(it != m_config.seedFinderConfig.zBinsCustomLooping.end()
        && "Inconsistent config zBinsCustomLooping");
    }
  }

  if (m_config.seedFinderConfig.useDetailedDoubleMeasurementInfo) {
    m_config.seedFinderConfig.getTopHalfStripLength.connect(
      [](const void*, const ActsExamples::SimSpacePoint &sp) -> float {
        return sp.topHalfStripLength();
      });

    m_config.seedFinderConfig.getBottomHalfStripLength.connect(
      [](const void*, const ActsExamples::SimSpacePoint &sp) -> float {
        return sp.bottomHalfStripLength();
      });

    m_config.seedFinderConfig.getTopStripDirection.connect(
      [](const void*, const ActsExamples::SimSpacePoint &sp) -> Acts::Vector3 {
        return sp.topStripDirection();
      });

    m_config.seedFinderConfig.getBottomStripDirection.connect(
      [](const void*, const ActsExamples::SimSpacePoint &sp) -> Acts::Vector3 {
        return sp.bottomStripDirection();
      });

    m_config.seedFinderConfig.getStripCenterDistance.connect(
      [](const void*, const ActsExamples::SimSpacePoint &sp) -> Acts::Vector3 {
        return sp.stripCenterDistance();
      });

    m_config.seedFinderConfig.getTopStripCenterPosition.connect(
      [](const void*, const ActsExamples::SimSpacePoint &sp) -> Acts::Vector3 {
        return sp.topStripCenterPosition();
      });
  }

  m_config.seedFinderConfig.seedFilter =
      std::make_unique<Acts::SeedFilter<ActsExamples::SimSpacePoint>>(
          m_config.seedFilterConfig);
}

ProcessCode TrackSeeding::execute(const Context &context) const {
  ACTS_DEBUG("Track seeding");

  const auto &spacePoints =
    context.eventStore.get<SpacePointContainer>(m_config.inputSpacePoints);

  std::vector<const SpacePoint*> spacePointPtrs;
  spacePointPtrs.reserve(spacePoints.size());

  Acts::Extent rRangeSPExtent;

  for (const auto &sp : spacePoints) {
    spacePointPtrs.push_back(&sp);
    rRangeSPExtent.extend({sp.x(), sp.y(), sp.z()});
  }

  // Extracts covariances for a space point.
  auto extractGlobalQuantities =
    [=](const SpacePoint &sp, float, float, float) ->
        std::pair<Acts::Vector3, Acts::Vector2> {
      Acts::Vector3 position{sp.x(), sp.y(), sp.z()};
      Acts::Vector2 covariance{sp.varianceR(), sp.varianceZ()};
      return std::make_pair(position, covariance);
    };

  auto bottomBinFinder = std::make_shared<Acts::BinFinder<SpacePoint>>(
      Acts::BinFinder<SpacePoint>(m_config.zBinNeighborsBottom,
                                  m_config.nPhiNeighbors));

  auto topBinFinder = std::make_shared<Acts::BinFinder<SpacePoint>>(
      Acts::BinFinder<SpacePoint>(m_config.zBinNeighborsTop,
                                  m_config.nPhiNeighbors));

  auto gridConfig = m_config.gridConfig;
  gridConfig.bFieldInZ = 0.; // FIXME: Otherwise sqrt(negative) in Acts
  auto grid = Acts::SpacePointGridCreator::createGrid<SpacePoint>(
      gridConfig);

  auto grouping = Acts::BinnedSPGroup<SpacePoint>(
                      spacePointPtrs.begin(),
                      spacePointPtrs.end(),
                      extractGlobalQuantities,
                      bottomBinFinder,
                      topBinFinder,
                      std::move(grid),
                      m_config.seedFinderConfig
                  );

  auto finder = Acts::Seedfinder<SpacePoint>(m_config.seedFinderConfig);

  // Run the seed finder.
  static thread_local SeedContainer seeds;
  seeds.clear();
  static thread_local decltype(finder)::State state;

  for (auto group = grouping.begin(); group != grouping.end(); ++group) {
    finder.createSeedsForGroup(
        state,
        std::back_inserter(seeds),
        group.bottom(),
        group.middle(),
        group.top(),
        rRangeSPExtent
    );
  }

  // Proto tracks are groups of measurement indices from the tracks seeds.
  static thread_local ProtoTrackContainer protoTracks;
  protoTracks.clear();
  protoTracks.reserve(seeds.size());

  for (const auto &seed : seeds) {
    auto &protoTrack = protoTracks.emplace_back();
    protoTrack.reserve(seed.sp().size());

    for (auto spPtr : seed.sp()) {
      protoTrack.push_back(spPtr->measurementIndex());
    }
  }

  ACTS_DEBUG("Created " << seeds.size() << " track seeds from "
                        << spacePointPtrs.size() << " space points");

  context.eventStore.add(
      m_config.outputSeeds, SeedContainer{seeds});
  context.eventStore.add(
      m_config.outputProtoTracks, ProtoTrackContainer{protoTracks});

  return ProcessCode::SUCCESS;
}

} // namespace Mpd::Tpc

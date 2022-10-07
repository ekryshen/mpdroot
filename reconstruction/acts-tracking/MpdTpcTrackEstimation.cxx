// This file is a part of the NICA project.
//
// Copyright (C) 2022 JINR

#include "MpdTpcContext.h"
#include "MpdTpcEventStorage.h"
#include "MpdTpcTrackEstimation.h"

#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>
#include <Acts/Seeding/EstimateTrackParamsFromSeed.hpp>

#include <array>
#include <cassert>
#include <map>

namespace Mpd::Tpc {

TrackEstimation::TrackEstimation(Config config, Acts::Logging::Level level):
  Algorithm("TrackEstimation", level),
  m_config(std::move(config)) {

  // Either seeds directly or proto tracks and space points.
  assert((!m_config.inputSeeds.empty()
      || (!m_config.inputSpacePoints.empty() && m_config.inputProtoTracks.empty()))
    && "Missing seeds or space point input collection");
  assert(!m_config.inputMeasurements.empty()
    && "Missing measurements input collection");
  assert(!m_config.outputTrackParameters.empty()
    && "Missing track parameters output collection");
  assert(!m_config.outputProtoTracks.empty()
    && "Missing proto tracks output collections");
  assert(m_config.trackingGeometry
    && "Missing tracking geometry");

  // Set up the track parameters covariance (the same for all tracks).
  m_covariance(Acts::eBoundLoc0, Acts::eBoundLoc0) =
      m_config.initialVarInflation[Acts::eBoundLoc0] * config.sigmaLoc0 *
      m_config.sigmaLoc0;
  m_covariance(Acts::eBoundLoc1, Acts::eBoundLoc1) =
      m_config.initialVarInflation[Acts::eBoundLoc1] * config.sigmaLoc1 *
      m_config.sigmaLoc1;
  m_covariance(Acts::eBoundPhi, Acts::eBoundPhi) =
      m_config.initialVarInflation[Acts::eBoundPhi] * config.sigmaPhi *
      m_config.sigmaPhi;
  m_covariance(Acts::eBoundTheta, Acts::eBoundTheta) =
      m_config.initialVarInflation[Acts::eBoundTheta] * config.sigmaTheta *
      m_config.sigmaTheta;
  m_covariance(Acts::eBoundQOverP, Acts::eBoundQOverP) =
      m_config.initialVarInflation[Acts::eBoundQOverP] * config.sigmaQOverP *
      m_config.sigmaQOverP;
  m_covariance(Acts::eBoundTime, Acts::eBoundTime) =
      m_config.initialVarInflation[Acts::eBoundTime] * m_config.sigmaT0 *
      m_config.sigmaT0;
}

ProcessCode TrackEstimation::execute(Context &context) const {
  // Measurements are necesary for retrieving the geometry identifier.
  const auto &measurements =
      context.eventStore.get<MeasurementContainer>(m_config.inputMeasurements);

  // Read seeds or create them from proto tracks and space points.
  SeedContainer seeds;

  if (!m_config.inputSeeds.empty()) {
    seeds = context.eventStore.get<SeedContainer>(m_config.inputSeeds);
    ACTS_DEBUG("Read " << seeds.size() << " seeds");
  } else {
    const auto &protoTracks = context.eventStore.get<ProtoTrackContainer>(
        m_config.inputProtoTracks);
    const auto &spacePoints = context.eventStore.get<SpacePointContainer>(
        m_config.inputSpacePoints);

    seeds = createSeeds(protoTracks, spacePoints);
    ACTS_DEBUG("Read " << protoTracks.size() << " proto tracks, and created "
                       << seeds.size() << " seeds");
  }

  TrackParametersContainer trackParameters;
  trackParameters.reserve(seeds.size());

  ProtoTrackContainer tracks;
  tracks.reserve(seeds.size());

  auto bCache = m_config.magneticField->makeCache(context.mContext);

  // Loop over all found seeds to estimate track parameters.
  for (size_t iseed = 0; iseed < seeds.size(); ++iseed) {
    const auto &seed = seeds[iseed];

    // Get the bottom space point and its reference surface.
    const auto bottomSP = seed.sp().front();
    const auto index = bottomSP->measurementIndex();

    const auto geoId = std::visit(
        [](const auto &measurement) {
          return measurement.sourceLink().geometryId();
        }, measurements[index]);

    const Acts::Surface *surface = m_config.trackingGeometry->findSurface(geoId);
    if (surface == nullptr) {
      ACTS_WARNING("surface with geoID " << geoId << " is not found");
      continue;
    }

    // Get the magnetic field at the bottom space point
    const auto x = bottomSP->x();
    const auto y = bottomSP->y();
    const auto z = bottomSP->z();

    auto field = m_config.magneticField->getField({x, y, z}, bCache);
    if (!field.ok()) {
      ACTS_ERROR("Field lookup error: " << field.error());
      return ProcessCode::ABORT;
    }

    // Estimate the track parameters from seed.
    auto optParams = Acts::estimateTrackParamsFromSeed(
        context.gContext,
        seed.sp().begin(),
        seed.sp().end(),
        *surface,
        *field,
        m_config.bFieldMin
    );

    if (!optParams.has_value()) {
      ACTS_WARNING("Estimation of track parameters for seed " << iseed << " failed");
      continue;
    } else {
      const auto &params = optParams.value();
      double charge = std::copysign(1, params[Acts::eBoundQOverP]);
      trackParameters.emplace_back(surface->getSharedPtr(), params, charge, m_covariance);

      // Create a proto track for this seed.
      ProtoTrack track;
      track.reserve(3);

      for (const auto& sp : seed.sp()) {
        track.push_back(sp->measurementIndex());
      }

      tracks.emplace_back(track);
    }
  }

  ACTS_DEBUG("Estimated " << trackParameters.size()
                          << " track parameters and " << tracks.size()
                          << " tracks");

  context.eventStore.add(m_config.outputTrackParameters, std::move(trackParameters));
  context.eventStore.add(m_config.outputProtoTracks, std::move(tracks));
 
  return ProcessCode::SUCCESS;
}

SeedContainer TrackEstimation::createSeeds(
    const ProtoTrackContainer &protoTracks,
    const SpacePointContainer &spacePoints) const {
  SeedContainer seeds;
  seeds.reserve(protoTracks.size());

  std::unordered_map<Index, const SpacePoint*> spMap;

  for (const auto &sp : spacePoints) {
    spMap.emplace(sp.measurementIndex(), &sp);
  }

  for (size_t itrack = 0; itrack < protoTracks.size(); ++itrack) {
    // The list of hits and the initial start parameters.
    const auto &protoTrack = protoTracks[itrack];

    if (protoTrack.size() < 3) {
      ACTS_WARNING("Proto track " << itrack << " size is less than 3");
      continue;
    }

    // Space points on the proto track.
    std::vector<const SpacePoint*> spacePointsOnTrack;
    spacePointsOnTrack.reserve(protoTrack.size());

    // Loop over the hit index on the proto track to find the space points.
    for (const auto &hitIndex : protoTrack) {
      auto it = spMap.find(hitIndex);
      if (it != spMap.end()) {
        spacePointsOnTrack.push_back(it->second);
      }
    }

    // At least three space points are required.
    if (spacePointsOnTrack.size() < 3) {
      continue;
    }

    // Sort the space points.
    std::sort(spacePointsOnTrack.begin(), spacePointsOnTrack.end(),
      [](const SpacePoint *lhs, const SpacePoint *rhs) {
        return std::hypot(lhs->r(), lhs->z()) <
               std::hypot(rhs->r(), rhs->z());
      });

    // Loop over the found space points to find the seed with maximum deltaR
    // between the the bottom and top space points.
    // TODO: add the check of deltaZ
    bool seedFound = false;
    std::array<size_t, 3> bestSPIndices;
    double maxDeltaR = std::numeric_limits<double>::min();

    for (size_t ib = 0; ib < spacePointsOnTrack.size() - 2; ++ib) {
      for (size_t im = ib + 1; im < spacePointsOnTrack.size() - 1; ++im) {
        for (size_t it = im + 1; it < spacePointsOnTrack.size(); ++it) {
          double bmDeltaR = std::abs(spacePointsOnTrack[im]->r() - spacePointsOnTrack[ib]->r());
          double mtDeltaR = std::abs(spacePointsOnTrack[it]->r() - spacePointsOnTrack[im]->r());

          if ((bmDeltaR >= m_config.deltaRMin) && (bmDeltaR <= m_config.deltaRMax) &&
              (mtDeltaR >= m_config.deltaRMin) && (mtDeltaR <= m_config.deltaRMax)) {
            if ((bmDeltaR + mtDeltaR) > maxDeltaR) {
              maxDeltaR = bmDeltaR + mtDeltaR;
              bestSPIndices = {ib, im, it};
              seedFound = true;
            }
          }
        }
      }
    }

    if (seedFound) {
      seeds.emplace_back(*spacePointsOnTrack[bestSPIndices[0]],
                         *spacePointsOnTrack[bestSPIndices[1]],
                         *spacePointsOnTrack[bestSPIndices[2]],
                          spacePointsOnTrack[bestSPIndices[1]]->z());
    }
  }

  return seeds;
}

} // namespace Mpd::Tpc

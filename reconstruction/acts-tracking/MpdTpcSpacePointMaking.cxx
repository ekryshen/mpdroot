// This file is a part of the NICA project.
//
// Copyright (C) 2022 JINR

#include "MpdTpcSpacePointMaking.h"

#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/TrackParametrization.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/Surfaces/Surface.hpp>

#include <algorithm>
#include <cassert>

namespace Mpd::Tpc {

SpacePointMaking::SpacePointMaking(Config config, Acts::Logging::Level level):
    Algorithm("SpacePointMaking", level),
    m_config(std::move(config)) {
  assert(!m_config.inputSourceLinks.empty()
    && "Missing source link input collection");
  assert(!m_config.inputMeasurements.empty()
    && "Missing measurement input collection");
  assert(!m_config.outputSpacePoints.empty()
    && "Missing space point output collection");
  assert(m_config.trackingGeometry
    && "Missing tracking geometry");
  assert(!m_config.geometrySelection.empty()
    && "Missing space point maker geometry selection");

  // Ensure geometry selection contains valid inputs.
  for (const auto &geoId : m_config.geometrySelection) {
    assert(!geoId.approach() && !geoId.boundary() && !geoId.sensitive()
      && "Invalid geometry selection: only volumes and layers are allowed");
  }

  // Remove geometry selection duplicates (ref < cmp).
  auto isDuplicate = [](Acts::GeometryIdentifier ref,
                        Acts::GeometryIdentifier cmp) {
    // Root node contains everything.
    if (ref.volume() == 0) {
      return true;
    }
    // Different volumes are separated.
    if (ref.volume() != cmp.volume()) {
      return false;
    }
    // Each volume consists of layers.
    return (ref.layer() == cmp.layer());
  };

  // Sort geometry selection to enable unique checking.
  auto geoSelBeg = m_config.geometrySelection.begin();
  auto geoSelEnd = m_config.geometrySelection.end();
  std::sort(geoSelBeg, geoSelEnd);

  auto geoSelLastUnique = std::unique(geoSelBeg, geoSelEnd, isDuplicate);
  if (geoSelLastUnique != geoSelEnd) {
    ACTS_WARNING("Removed " << std::distance(geoSelLastUnique, geoSelEnd)
                            << " geometry selection duplicates");
    m_config.geometrySelection.erase(geoSelLastUnique, geoSelEnd);
  }

  ACTS_INFO("Space point geometry selection:");
  for (const auto &geoId : m_config.geometrySelection) {
    ACTS_INFO(" " << geoId);
  }
}

ActsExamples::ProcessCode SpacePointMaking::execute(
    const ActsExamples::AlgorithmContext &context) const {
  ACTS_DEBUG("Space point making");

  const auto &sourceLinks =
      context.eventStore.get<SourceLinkContainer>(m_config.inputSourceLinks);
  const auto &measurements =
      context.eventStore.get<MeasurementContainer>(m_config.inputMeasurements);

  SpacePointContainer spacePoints;
  spacePoints.reserve(sourceLinks.size());

  for (auto geoId : m_config.geometrySelection) {
    // Select volume/layer depending on what is set in the geometry id.
    auto range = ActsExamples::selectLowestNonZeroGeometryObject(
        sourceLinks, geoId);
    auto groupedByModule = makeGroupBy(
        range, ActsExamples::detail::GeometryIdGetter());

    for (auto [moduleGeoId, moduleSourceLinks] : groupedByModule) {
      const auto *surface = m_config.trackingGeometry->findSurface(moduleGeoId);

      if (!surface) {
        ACTS_ERROR("Could not find surface " << moduleGeoId);
        return ActsExamples::ProcessCode::ABORT;
      }

      for (auto &sourceLink : moduleSourceLinks) {
        auto index = sourceLink.get().index();
        auto [globalPos, var] = localToGlobal(
            surface, context.geoContext, measurements, index);
        ACTS_DEBUG("Space point " << index << ": ("
                                  << globalPos[0] << ", "
                                  << globalPos[1] << ", "
                                  << globalPos[2] << ")");

        // Construct space point in global coordinates.
        spacePoints.emplace_back(globalPos, var[0], var[1], index);
      }
    }
  }

  spacePoints.shrink_to_fit();

  ACTS_DEBUG("Created " << spacePoints.size() << " space points");
  context.eventStore.add(m_config.outputSpacePoints, std::move(spacePoints));

  return ActsExamples::ProcessCode::SUCCESS;
}

SpacePointMaking::GlobalPosition SpacePointMaking::localToGlobal(
    const Acts::Surface *surface,
    const Acts::GeometryContext &geoContext,
    const MeasurementContainer &measurements,
    Index index) const {
  assert(surface && "Missing surface");

  // Extract a local position/covariance independent from the concrete
  // measurement content. Since we do not know if and where the local
  // parameters are contained in the measurement parameters vector, they
  // are transformed to the bound space where we do know their location.
  // If the local parameters are not measured, this results in a zero
  // location, which is a reasonable default fall-back.
  auto [localPos, localCov] = std::visit(
      [](const auto &measurement) {
        auto expander = measurement.expander();
        Acts::BoundVector par = expander * measurement.parameters();
        Acts::BoundSymMatrix cov =
            expander * measurement.covariance() * expander.transpose();
        Acts::Vector2 lpar(par[Acts::eBoundLoc0], par[Acts::eBoundLoc1]);
        Acts::SymMatrix2 lcov = cov.block<2, 2>(0, 0);
        return std::make_pair(lpar, lcov);
      },
      measurements[index]);

  // Transform the local position to global coordinates.
  Acts::Vector3 fakeMom(1, 1, 1);
  Acts::Vector3 globalPos = surface->localToGlobal(geoContext, localPos, fakeMom);
  Acts::RotationMatrix3 rotLocalToGlobal =
      surface->referenceFrame(geoContext, globalPos, fakeMom);

  // The space point requires only the variance of the transverse and
  // longitudinal position. Reduce computations by transforming the
  // covariance directly from local to rho/z.
  //
  // Compute Jacobian from the global coordinates to rho/z:
  //
  //         rho = sqrt(x² + y²)
  // drho/d{x,y} = (1 / sqrt(x² + y²)) * 2 * {x,y}
  //             = 2 * {x,y} / r
  //       dz/dz = 1
  auto x = globalPos[Acts::ePos0];
  auto y = globalPos[Acts::ePos1];
  auto scale = 2 / std::hypot(x, y);

  Acts::ActsMatrix<2, 3> jacXyzToRhoZ = Acts::ActsMatrix<2, 3>::Zero();
  jacXyzToRhoZ(0, Acts::ePos0) = scale * x;
  jacXyzToRhoZ(0, Acts::ePos1) = scale * y;
  jacXyzToRhoZ(1, Acts::ePos2) = 1;
  // Compute Jacobian from local coordinates to rho/z.
  Acts::ActsMatrix<2, 2> J = jacXyzToRhoZ * rotLocalToGlobal.block<3, 2>(0, 0);
  // Compute rho/z variance.
  Acts::ActsVector<2> var = (J * localCov * J.transpose()).diagonal();

  return {globalPos, var};
}

} // namespace Mpd::Tpc

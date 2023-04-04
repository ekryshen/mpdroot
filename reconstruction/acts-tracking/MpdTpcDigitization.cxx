// This file is a part of the NICA project.
//
// Copyright (C) 2022 JINR

#include "MpdTpcDetector.h"
#include "MpdTpcDigitization.h"
#include "MpdTpcInputHit.h"

#include "BaseTpcSectorGeo.h"

#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/Units.hpp>
#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/Propagator/Navigator.hpp>
#include <Acts/Surfaces/Surface.hpp>

#include <cassert>
#include <type_traits>

using namespace Acts::UnitLiterals;

namespace Mpd::Tpc {

Digitization::Digitization(Config config, Acts::Logging::Level level):
    Algorithm("Digitization", level),
    m_config(std::move(config)) {
  assert(!m_config.inputSimHits.empty()
    && "Missing simulation hit input collection");
  assert(!m_config.outputSourceLinks.empty()
    && "Missing source link output collection");
  assert(!m_config.outputMeasurements.empty()
    && "Missing measurement output collection");
  assert(m_config.detector
    && "Missing detector");
}

ActsExamples::ProcessCode Digitization::execute(
    const ActsExamples::AlgorithmContext &context) const {
  ACTS_DEBUG("Digitization");

  const auto &simHits =
      context.eventStore.get<InputHitContainer>(m_config.inputSimHits);

  // Stores source link objects referenced by measurements and source links.
  std::vector<SourceLink> sourceLinkStorage;
  sourceLinkStorage.reserve(simHits.size());

  SourceLinkContainer sourceLinks;
  sourceLinks.reserve(simHits.size());

  MeasurementContainer measurements;
  measurements.reserve(simHits.size());

  for (const auto &simHit : simHits) {
    const auto &pos = simHit.position;
    const auto &mom = simHit.momentum;

    auto surface = m_config.detector->getSurface(context.geoContext, pos);
    auto moduleGeoId = surface->geometryId();
    ACTS_VERBOSE("Module identifier is " << moduleGeoId);

    double deltaR;
    if (Detector::useBaseTpcSectorGeo) {
      BaseTpcSectorGeo sectorGeo = BaseTpcSectorGeo();
      double maxPadHeight = std::max(sectorGeo.PAD_HEIGHT[0], sectorGeo.PAD_HEIGHT[1]);
      deltaR = maxPadHeight * 1_cm;
    } else {
      deltaR = Detector::DeltaR;
    }
    auto result = surface->globalToLocal(context.geoContext, pos, mom, 3. * deltaR);
    if (!result.ok()) {
      ACTS_DEBUG("Point " << pos << " is far from surface "
                          << surface->center(context.geoContext));
      assert(false);
    }
 
    ActsExamples::Index index = measurements.size();

    SourceLink sourceLink(moduleGeoId, index);
    sourceLinkStorage.push_back(std::move(sourceLink));

    auto &sourceLinkRef = sourceLinkStorage.back();
    sourceLinks.insert(sourceLinkRef);

    // FIXME: setup realistic variances for (loc0, loc1).
    Acts::SymMatrix2 cov = Acts::ActsMatrix<2, 2>::Zero();
    cov.diagonal() = Acts::Vector2(m_config.sigmaLoc0, m_config.sigmaLoc1);

    auto measurement = makeMeasurement(
        sourceLinkRef,
        result.value(),
        cov,
        Acts::eBoundLoc0,
        Acts::eBoundLoc1
    );

    measurements.push_back(measurement);
  }

  context.eventStore.add(
      m_config.outputSourceLinks, std::move(sourceLinks));
  context.eventStore.add(
      m_config.outputMeasurements, std::move(measurements));
  context.eventStore.add(
      m_config.outputSourceLinks + "_storage", std::move(sourceLinkStorage));

  return ActsExamples::ProcessCode::SUCCESS;
}

} // namespace Mpd::Tpc

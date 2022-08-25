// This file is a part of the NICA project.
//
// Copyright (C) 2022 JINR

#include "MpdTpcContext.h"
#include "MpdTpcDetector.h"
#include "MpdTpcDigitization.h"
#include "MpdTpcEventStorage.h"
#include "MpdTpcInputHit.h"

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

ProcessCode Digitization::execute(const Context &context) const {
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

    auto surface = m_config.detector->getSurface(context.gContext, pos);
    auto moduleGeoId = surface->geometryId();
    ACTS_VERBOSE("Module identifier is " << moduleGeoId);

    auto result = surface->globalToLocal(context.gContext, pos, mom, 3. * Detector::DeltaR);
    if (!result.ok()) {
      ACTS_DEBUG("Point " << pos << " is far from surface "
                          << surface->center(context.gContext));
      assert(false);
    }
 
    ActsExamples::Index index = measurements.size();

    SourceLink sourceLink(moduleGeoId, index);
    sourceLinkStorage.push_back(std::move(sourceLink));

    auto &sourceLinkRef = sourceLinkStorage.back();
    sourceLinks.insert(sourceLinkRef);

    // FIXME: setup realistic variances for (loc0, loc1).
    Acts::SymMatrix2 cov = Acts::ActsMatrix<2, 2>::Zero();
    cov.diagonal() = Acts::Vector2(0.1_mm, 0.1_mm);

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

  return ProcessCode::SUCCESS;
}

} // namespace Mpd::Tpc

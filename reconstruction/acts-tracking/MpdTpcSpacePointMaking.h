// This file is a part of the NICA project.
//
// Copyright (C) 2022 JINR

#pragma once

#include "MpdTpcAlgorithm.h"
#include "MpdTpcEventData.h"

#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <memory>
#include <string>
#include <vector>

namespace Acts {
  class Surface;
  class TrackingGeometry;
} // namespace Acts

namespace Mpd::Tpc {

/// @brief Creates space points from measurements (based on Acts examples).
///
/// This implements simple space point construction algorithm, where each
/// surface-based measurement translates into one space point using the
/// surface local-to-global transform.
///
/// The algorithm takes both the source links and measurements container as
/// input. The source link container is geometry-sorted and each element is
/// small compared to a measurement. The geometry selection is therefore much
/// easier to perform on the source links than on the unsorted measurements.
///
/// There are no explicit requirements on the content of the input measurements.
/// If no local positions are measured, the transformed global positions will
/// always be the position of the module origin.
class SpacePointMaking final : public Algorithm {
public:
  struct Config final {
    std::string inputSourceLinks;
    std::string inputMeasurements;
    std::string outputSpacePoints;

    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;

    /// For which part of the detector geometry should space points be created.
    ///
    /// Only volumes and layers can be set. Zero values can be used as wildcards
    /// to select larger parts of the hierarchy, i.e. setting only the volume
    /// selects all measurements within that volume. Adding a single identifier
    /// with all components set to zero selects all available measurements. The
    /// selection must not have duplicates.
    std::vector<Acts::GeometryIdentifier> geometrySelection;
  };

  SpacePointMaking(Config config, Acts::Logging::Level level);

  ActsExamples::ProcessCode execute(
      const ActsExamples::AlgorithmContext &context) const override;

  const Config &config() const { return m_config; }

private:
  /// Global position and the rho/z variances.
  using GlobalPosition = std::pair<Acts::Vector3, Acts::Vector2>;

  /// Gets the global position for the local measurement.
  GlobalPosition localToGlobal(const Acts::Surface *surface,
                               const Acts::GeometryContext &gContext,
                               const MeasurementContainer &measurements,
                               ActsExamples::Index index) const;

  Config m_config;
};

} // namespace Mpd::Tpc

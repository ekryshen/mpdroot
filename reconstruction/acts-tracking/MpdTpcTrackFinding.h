// This file is a part of the NICA project.
//
// Copyright (C) 2022 JINR

#pragma once

#include "MpdTpcAlgorithm.h"
#include "MpdTpcEventData.h"

#include <Acts/EventData/VectorMultiTrajectory.hpp>
#include <Acts/Propagator/EigenStepper.hpp>
#include <Acts/Propagator/Propagator.hpp>
#include <Acts/Propagator/Navigator.hpp>
#include <Acts/TrackFinding/MeasurementSelector.hpp>
#include <Acts/TrackFinding/CombinatorialKalmanFilter.hpp>
#include <Acts/TrackFinding/MeasurementSelector.hpp>
#include <Acts/TrackFitting/GainMatrixSmoother.hpp>
#include <Acts/TrackFitting/GainMatrixUpdater.hpp>

#include <Acts/Utilities/Logger.hpp>

#include <memory>
#include <string>
#include <vector>

namespace Acts {
  class MagneticFieldProvider;
  class TrackingGeometry;
  class Surface;
} // namespace Acts

namespace Mpd::Tpc {

/// @brief Finds tracks using the Kalman filter (based on Acts examples).
class TrackFinding final : public Algorithm {
public:
  using Accessor = SourceLinkAccessor;
  using Iterator = SourceLinkIterator;
  using Container = SourceLinkContainer;

  using Updater = Acts::GainMatrixUpdater;
  using Smoother = Acts::GainMatrixSmoother;
  using Selector = Acts::MeasurementSelector;

  using Stepper = Acts::EigenStepper<>;
  using Navigator = Acts::Navigator;
  using Propagator = Acts::Propagator<Stepper, Navigator>;

  using Filter = Acts::CombinatorialKalmanFilter<Propagator, Trajectory>;
  using Extensions = Acts::CombinatorialKalmanFilterExtensions<Trajectory>;

  using Options = Acts::CombinatorialKalmanFilterOptions<Iterator, Trajectory>;
  using Result = Acts::CombinatorialKalmanFilterResult<Trajectory>;
  using Results = std::vector<Acts::Result<Result>>;

  struct Config final {
    std::string inputMeasurements;
    std::string inputSourceLinks;
    std::string inputInitialTrackParameters;
    std::string outputTrajectories;
    std::string outputTrackCandidates;

    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField;
    std::shared_ptr<const Acts::Surface> referenceSurface;

    Acts::MeasurementSelector::Config measurementSelectorConfig;
    Acts::Navigator::Config navigatorConfig;
    Acts::PropagatorPlainOptions propagatorOptions;

    // Kalman filter options.
    bool multipleScattering;
    bool energyLoss;
    bool smoothing;

    bool computeSharedHits;

    /// Minimal track length (tracks w/ shorter length are ignored).
    size_t trackMinLength;
    /// Segment length bound (minimal number of new hits in a row).
    size_t newHitsInRow;
    /// Coverage ratio (new hits in a track / track size).
    double newHitsRatio;
  };

  TrackFinding(Config config, Acts::Logging::Level level);
  ProcessCode execute(const Context &context) const override;
  const Config &config() const { return m_config; }

private:
  /// Sets shared-hit flags on trajectories' states.
  void computeSharedHits(const Container &sourceLinks, Results &results) const;

  /// Post-processing: removes duplicates and constructs track candidates.
  void constructTrackCandidates(const Context &context,
                                const Container &sourceLinks,
                                Results &results) const;

  Config m_config;
};

} // namespace Mpd::Tpc

// This file is a part of the NICA project.
//
// Copyright (C) 2022 JINR

#pragma once

#include "MpdTpcAlgorithm.h"
#include "MpdTpcEventData.h"

#include <Acts/Definitions/TrackParametrization.hpp>
#include <Acts/Definitions/Units.hpp>
#include <Acts/Seeding/Seed.hpp>

#include <Rtypes.h>

#include <memory>
#include <string>

using namespace Acts::UnitLiterals;

namespace Acts {
  class TrackingGeometry;
  class MagneticFieldProvider;
} // namespace Acts

namespace Mpd::Tpc {

/// @brief Estimates track parameters for track seeds (based on Acts examples).
class TrackEstimation final : public Algorithm {
public:
  struct Config final {
    std::string inputSeeds;
    std::string inputSpacePoints;
    std::string inputProtoTracks;
    std::string inputMeasurements;
    std::string outputTrackParameters;
    std::string outputProtoTracks;

    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField;

    /// Minimum deltaR between space points in a seed.
    float deltaRMin;
    /// Maximum deltaR between space points in a seed.
    float deltaRMax;
    /// The minimum magnetic field to trigger the track parameters estimation.
    Double_t bFieldMin;
    /// Constant term of the loc0 resolution.
    Double_t sigmaLoc0;
    /// Constant term of the loc1 resolution.
    Double_t sigmaLoc1;
    /// Phi angular resolution.
    Double_t sigmaPhi;
    /// Theta angular resolution.
    Double_t sigmaTheta;
    /// q/p resolution.
    Double_t sigmaQOverP;
    /// Time resolution.
    Double_t sigmaT0;
    /// Inflate tracks.
    std::array<Double_t, 6> initialVarInflation;
  };

  TrackEstimation(Config config, Acts::Logging::Level level);

  ActsExamples::ProcessCode execute(
      const ActsExamples::AlgorithmContext &context) const override;

  const Config &config() const { return m_config; }

private:
  /// Creates seeds from proto tracks and space points.
  SeedContainer createSeeds(const ProtoTrackContainer &protoTracks,
                            const SpacePointContainer &spacePoints) const;

  /// Track parameters covariance (assumed to be the same
  /// for all estimated track parameters for the moment).
  Acts::BoundSymMatrix m_covariance = Acts::BoundSymMatrix::Zero();

  Config m_config;
};

} // namespace Mpd::Tpc

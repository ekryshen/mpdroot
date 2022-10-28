// This file is a part of the NICA project.
//
// Copyright (C) 2022 JINR

#pragma once

#include "ActsExamples/TGeoDetector/TGeoDetector.hpp"

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/Units.hpp>
#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <memory>
#include <unordered_map>

using namespace Acts::UnitLiterals;

class TGeoManager;

namespace Mpd::Tpc {

/// @brief Constructs TPC's virtual sensors (required by Acts).
class Detector final : public ActsExamples::TGeoDetector {
public:
  /// TPC sensitive volume.
  static constexpr auto GasVolume = "tpc01sv";

  /// TPC is represented as two symmetric nodes.
  static constexpr auto HasTwoNodes = true;

  static constexpr auto Rmin       =  40.29_cm;    // FIXME
  static constexpr auto Rmax       =  Rmin + 84.7_cm; //132.15_cm;   // FIXME
  static constexpr auto Zmin       = -2.*82.5_cm; // FIXME: 81.95
  static constexpr auto Zmax       =  2.*82.5_cm; // FIXME: 81.95
  static constexpr auto NumLayers  =  200u;
  static constexpr auto NumSectors =  180u;
  static constexpr auto NumZ       =  2u;
  static constexpr auto DeltaR     =  (Rmax - Rmin) / NumLayers;
  static constexpr auto DeltaZ     =  (Zmax - Zmin) / NumZ;
  static constexpr auto DeltaPhi   =  (2.*M_PI) / NumSectors;
  static constexpr auto PhiStart   =  0.5*DeltaPhi;

  Detector(const std::string &jsonFile,
           Acts::Logging::Level level = Acts::Logging::DEBUG):
      Detector("" /* No import */, jsonFile, level) {}

  Detector(const std::string &rootFile,
           const std::string &jsonFile,
           Acts::Logging::Level level = Acts::Logging::DEBUG);

  /// Constructs and returns the tracking geometry.
  std::shared_ptr<const Acts::TrackingGeometry> getGeometry();

  /// Gets the module identifier for the given position and momentum.
  std::shared_ptr<const Acts::Surface> getSurface(
      const Acts::GeometryContext &context, const Acts::Vector3 &position);

private:
  /// Adds virtual sensors to the TPC sensitive volume.
  bool editGeometry(TGeoManager *geoManager);

private:
  const Acts::Logger &logger() const { return *m_logger; }

  ActsExamples::TGeoDetector::Config m_config;
  std::unique_ptr<const Acts::Logger> m_logger;
  std::shared_ptr<const Acts::TrackingGeometry> m_trackingGeometry;
  std::unordered_map<uint64_t, std::shared_ptr<const Acts::Surface>> m_surfaces;
};

} // namespace Mpd::Tpc

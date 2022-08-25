// This file is a part of the NICA project.
//
// Copyright (C) 2022 JINR

#pragma once

#include "MpdTpcAlgorithm.h"
#include "MpdTpcDetector.h"
#include "MpdTpcEventData.h"

#include <Acts/Utilities/Logger.hpp>

#include <string>

namespace Acts {
  class TrackingGeometry;
} // namespace Acts

namespace Mpd::Tpc {

/// @brief Converts (simulation) hits to measurements and source links.
class Digitization final : public Algorithm {
public:
  struct Config final {
    std::string inputSimHits;
    std::string outputSourceLinks;
    std::string outputMeasurements;

    std::shared_ptr<Detector> detector;
  };

  Digitization(Config config, Acts::Logging::Level level);
  ProcessCode execute(const Context &context) const override;
  const Config &config() const { return m_config; }

private:
  Config m_config;
};

} // namespace Mpd::Tpc

// This file is a part of the NICA project.
//
// Copyright (C) 2022 JINR

#pragma once

#include "MpdTpcAlgorithm.h"
#include "MpdTpcEventData.h"

#include <Acts/Seeding/Seed.hpp>
#include <Acts/Seeding/SeedFilterConfig.hpp>
#include <Acts/Seeding/SeedfinderConfig.hpp>
#include <Acts/Seeding/SpacePointGrid.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <string>
#include <vector>

namespace Mpd::Tpc {

/// @brief Constructs track seeds from space points (based on Acts examples).
class TrackSeeding final : public Algorithm {
public:
  struct Config final {
    std::string inputSpacePoints;
    std::string outputSeeds;
    std::string outputProtoTracks;

    Acts::SeedFilterConfig seedFilterConfig;
    Acts::SeedfinderConfig<SpacePoint> seedFinderConfig;
    Acts::SpacePointGridConfig gridConfig;

    /// Vectors containing the map of Z bins in the top and bottom layers.
    std::vector<std::pair<int, int>> zBinNeighborsTop;
    std::vector<std::pair<int, int>> zBinNeighborsBottom;

    /// Number of phiBin neighbors at each side of the current bin.
    int nPhiNeighbors;
  };

  TrackSeeding(Config config, Acts::Logging::Level level);
  ProcessCode execute(Context &context) const override;
  const Config &config() const { return m_config; }

private:
  Config m_config;
};

} // namespace Mpd::Tpc

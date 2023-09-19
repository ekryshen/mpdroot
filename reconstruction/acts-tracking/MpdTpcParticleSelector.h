// This file is a part of the NICA project.
//
// Copyright (C) 2023 JINR

#pragma once

#include "MpdTpcAlgorithm.h"
#include "MpdTpcEventData.h"

#include "ActsExamples/TruthTracking/TruthSeedSelector.hpp"

#include <Acts/Utilities/Logger.hpp>

#include <Rtypes.h>

#include <string>
#include <vector>

namespace Mpd::Tpc {

class ParticleSelector final : public Algorithm {
public:
  struct Config final {
    ActsExamples::TruthSeedSelector::Config truthSeedSelectorConfig;

    Bool_t selectorEnabled;
    Bool_t primaryParticlesOnly;
  };

  ParticleSelector(Config config, Acts::Logging::Level level);

  ActsExamples::ProcessCode execute(
      const ActsExamples::AlgorithmContext &context) const override;

  const Config &config() const { return m_config; }

private:
  Config m_config;
};

} // namespace Mpd::Tpc

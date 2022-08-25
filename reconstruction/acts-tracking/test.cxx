// This file is a part of the NICA project.
//
// Copyright (C) 2022 JINR

#include "MpdTpcDetector.h"
#include "MpdTpcInputHit.h"
#include "MpdTpcRunner.h"

#include <Acts/Utilities/Logger.hpp>

#include <cmath>
#include <random>

int main() {
  Mpd::Tpc::Runner runner(
      "../../geometry/tpc_v9_with_materials.root",
      Acts::Logging::VERBOSE
  );

  const auto seed = 0;
  std::mt19937 gen(seed);

  std::uniform_real_distribution<double>
      rhoDist(Mpd::Tpc::Detector::Rmin, Mpd::Tpc::Detector::Rmax);
  std::uniform_real_distribution<double>
      phiDist(0., 2. * M_PI);
  std::uniform_real_distribution<double>
      zDist(Mpd::Tpc::Detector::Zmin, Mpd::Tpc::Detector::Zmax);

  const auto nHits = 10000u;

  Mpd::Tpc::InputHitContainer hits;
  hits.reserve(nHits);

  for (size_t i = 0; i < nHits; i++) {
    const auto trackId = 0;
    const auto detectorId = 0;
    const auto rho = rhoDist(gen);
    const auto phi = phiDist(gen);
    const auto x = rho * std::cos(phi);
    const auto y = rho * std::sin(phi);
    const auto z = zDist(gen);
    const auto px = 1.;
    const auto py = 1.;
    const auto pz = 1.;
    const auto time = zDist(gen);
    const auto length = zDist(gen);
    const auto energyLoss = 1.;

    hits.emplace_back(Mpd::Tpc::InputHit{
        trackId,
        detectorId,
        {x, y, z},
        {px, py, pz},
        time,
        length,
        energyLoss
    });
  }

  runner.execute(hits);

  return 0;
}

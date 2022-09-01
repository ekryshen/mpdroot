// This file is a part of the NICA project.
//
// Copyright (C) 2022 JINR

#pragma once

#include <Acts/Definitions/Algebra.hpp>

#include <cmath>
#include <ostream>
#include <vector>

namespace Mpd::Tpc {

/// @brief Represents a (simulation) hit.
struct InputHit final {
  int trackId;            ///< Simulation track. 
  int detectorId;         ///< Detector (pad).
  Acts::Vector3 position; ///< XYZ in mm.
  Acts::Vector3 momentum; ///< Pxyz in GeV.
};

inline std::ostream &operator<<(std::ostream &out, const InputHit &hit) {
  const auto x = hit.position[0];
  const auto y = hit.position[1];
  const auto z = hit.position[2];

  const auto pX = hit.momentum[0];
  const auto pY = hit.momentum[1];

  return out << "(" << x << ", " << y << ", " << z << "), "
             << "R=" << std::hypot(x, y) << ", Pt=" << std::hypot(pX, pY);
}

using InputHitContainer = std::vector<InputHit>;

} // namespace Mpd::Tpc

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
  int trackId;                  ///< Simulation track. 
  int detectorId;               ///< Detector (pad).
  Acts::Vector3 position;       ///< XYZ in mm.
  Acts::Vector3 momentum;       ///< Pxyz in FIXME.
  Acts::ActsScalar time;        ///< Time in FIXME.
  Acts::ActsScalar length;      ///< Length in mm.
  Acts::ActsScalar energyLoss;  ///< Energy loss in FIXME.
};

inline std::ostream &operator<<(std::ostream &out, const InputHit &hit) {
  const auto x = hit.position[0];
  const auto y = hit.position[1];
  const auto z = hit.position[2];

  return out << "Hit(" << x << ", " << y << ", " << z << "), "
             << "R=" << std::hypot(x, y);
}

using InputHitContainer = std::vector<InputHit>;

} // namespace Mpd::Tpc

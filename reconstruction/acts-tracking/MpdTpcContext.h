// This file is a part of the NICA project.
//
// Copyright (C) 2022 JINR

#pragma once

#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/MagneticField/MagneticFieldContext.hpp>
#include <Acts/Utilities/CalibrationContext.hpp>

#include <memory>

namespace Mpd::Tpc {

/// @brief Common context of the tracking subtasks.
struct Context {
  ActsExamples::WhiteBoard eventStore; ///< Event data storage.
  Acts::GeometryContext gContext;      ///< Geometry context.
  Acts::MagneticFieldContext mContext; ///< Magnetic field context.
  Acts::CalibrationContext cContext;   ///< Calibration context.
};

} // namespace Mpd::Tpc

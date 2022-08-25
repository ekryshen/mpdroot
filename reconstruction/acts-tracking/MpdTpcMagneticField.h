// This file is a part of the NICA project.
//
// Copyright (C) 2022 JINR

#pragma once

#include <Acts/Definitions/Units.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>

#include <memory>

using namespace Acts::UnitLiterals;

namespace Mpd::Tpc::MagneticField {

static constexpr auto Bz = 0.5_T;

/// @brief Provides information on the MPD magnetic field.
std::shared_ptr<const Acts::MagneticFieldProvider> build();

} // namespace Mpd::Tpc::MagnetricField

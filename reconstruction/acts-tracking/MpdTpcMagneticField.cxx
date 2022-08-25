// This file is part of the NICA project.
//
// Copyright (C) 2022 JINR

#include "MpdTpcMagneticField.h"

#include <Acts/MagneticField/ConstantBField.hpp>

namespace Mpd::Tpc::MagneticField {

std::shared_ptr<const Acts::MagneticFieldProvider> build() {
  static auto field = std::shared_ptr<const Acts::MagneticFieldProvider>(
      new Acts::ConstantBField({0., 0., Bz}));
  return field;
};

} // namespace Mpd::Tpc::MagneticField

// This file is a part of the NICA project.
//
// Copyright (C) 2022 JINR

#pragma once

#include "ActsExamples/Framework/ProcessCode.hpp"

#include <Acts/Utilities/Logger.hpp>

#include <memory>
#include <string>

namespace Mpd::Tpc {

class Context;

using ProcessCode = ActsExamples::ProcessCode;

/// @brief Base class for algorithms (subtasks) used in the Acts-based tracker.
class Algorithm {
public:
  Algorithm(std::string name,
            Acts::Logging::Level level = Acts::Logging::INFO):
      m_name(std::move(name)),
      m_logger(Acts::getDefaultLogger(m_name, level)) {}

  virtual ~Algorithm() = default;

  std::string name() const { return m_name; }

  /// Executes the algorithm for one event.
  virtual ProcessCode execute(Context &context) const = 0;

protected:
  const Acts::Logger &logger() const { return *m_logger; }

private:
  std::string m_name;
  std::unique_ptr<const Acts::Logger> m_logger;
};

} // namespace Mpd::Tpc

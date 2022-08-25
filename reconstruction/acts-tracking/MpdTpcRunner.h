// This file is a part of the NICA project.
//
// Copyright (C) 2022 JINR

#pragma once

#include "MpdTpcConfig.h"
#include "MpdTpcContext.h"
#include "MpdTpcEventData.h"
#include "MpdTpcEventStorage.h"
#include "MpdTpcInputHit.h"

#include <Acts/Utilities/Logger.hpp>

#include <ostream>
#include <string>

namespace Mpd::Tpc {

/// @brief Provides an API for the tracker.
class Runner final {
public:
  Runner(Acts::Logging::Level level = Acts::Logging::DEBUG):
      Runner("" /* No import */, level) {}

  Runner(const std::string &rootFile,
         Acts::Logging::Level level = Acts::Logging::DEBUG):
      m_logger(Acts::getDefaultLogger("Runner", level)),
      m_level(level),
      m_config(rootFile, m_level),
      m_store(m_level),
      m_context(m_store) {}

  const TrajectoriesContainer &execute(const InputHitContainer &hits);
  size_t countSimTracks(const InputHitContainer &hits) const;

  void logInput() const;
  void logOutput() const;

private:
  const Acts::Logger &logger() const { return *m_logger; }

  std::unique_ptr<const Acts::Logger> m_logger;
  Acts::Logging::Level m_level;

  Config m_config;
  EventStorage m_store;
  Context m_context;
};

} // namespace Mpd::Tpc

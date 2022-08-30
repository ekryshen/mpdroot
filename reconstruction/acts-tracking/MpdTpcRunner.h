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
  Runner(const std::string &jsonFile,
         Acts::Logging::Level level = Acts::Logging::DEBUG):
      Runner("" /* No import */, jsonFile, level) {}

  Runner(const std::string &rootFile,
         const std::string &jsonFile,
         Acts::Logging::Level level = Acts::Logging::DEBUG):
      m_logger(Acts::getDefaultLogger("Runner", level)),
      m_level(level),
      m_config(rootFile, jsonFile, m_level),
      m_store(m_level),
      m_context(m_store) {}

  const TrajectoriesContainer &execute(const InputHitContainer &hits);

private:
  void logInput() const;

  void logOutput() const;

  void logHit(size_t hitId, const InputHit &hit) const;

  void logHits(const InputHitContainer &hits) const;

  void logTrack(const std::string &prefix,
                size_t trackId,
                const ProtoTrack &track) const;

  void logTracks(const std::string &prefix,
                 const ProtoTrackContainer &tracks) const;

  const Acts::Logger &logger() const { return *m_logger; }

  std::unique_ptr<const Acts::Logger> m_logger;
  Acts::Logging::Level m_level;

  Config m_config;
  EventStorage m_store;
  Context m_context;
};

} // namespace Mpd::Tpc

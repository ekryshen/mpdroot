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

#include <map>
#include <ostream>
#include <string>

namespace Mpd::Tpc {

/// Debug information on a single simulation track.
struct Quality {
  /// Maximum relative length among all found tracks.
  double quality;
  /// Percent of points belonging to the original track.
  double accuracy;
  /// Length of the best found track.
  size_t length;
  /// Length of the original track.
  size_t realLength;
  /// Number of found tracks.
  size_t nTracks;
};

/// Maps simulation tracks to the recognition statistics.
using Statistics = std::map<int, Quality>;

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

  const ProtoTrackContainer &execute(const InputHitContainer &hits);

  size_t getTracksNumber() const;
  Statistics getStatistics() const;

private:

  // Logging.
  void logInput() const;
  void logOutput() const;
  void logHit(size_t hitId, const InputHit &hit) const;
  void logHits(const InputHitContainer &hits) const;
  void logTrack(const std::string &prefix,
                size_t trackId,
                const ProtoTrack &track) const;
  void logTracks(const std::string &prefix,
                 const ProtoTrackContainer &tracks) const;

  // Collecting statistics.
  void checkTrack(const InputHitContainer &hits,
                  size_t trackId,
                  const ProtoTrack &track,
                  Statistics &statistics) const;
  void checkTracks(const InputHitContainer &hits,
                   const ProtoTrackContainer &tracks,
                   Statistics &statistics) const;

  const Acts::Logger &logger() const { return *m_logger; }

  std::unique_ptr<const Acts::Logger> m_logger;
  Acts::Logging::Level m_level;

  Config m_config;
  EventStorage m_store;
  Context m_context;
};

} // namespace Mpd::Tpc

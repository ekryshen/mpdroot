// This file is a part of the NICA project.
//
// Copyright (C) 2022 JINR

#include "MpdTpcContext.h"
#include "MpdTpcConfig.h"
#include "MpdTpcDigitization.h"
#include "MpdTpcEventData.h"
#include "MpdTpcEventStorage.h"
#include "MpdTpcInputHit.h"
#include "MpdTpcRunner.h"
#include "MpdTpcSpacePointMaking.h"
#include "MpdTpcTrackEstimation.h"
#include "MpdTpcTrackFinding.h"
#include "MpdTpcTrackSeeding.h"

#include <sstream>
#include <unordered_map>
#include <vector>

namespace Mpd::Tpc {

const TrajectoriesContainer &Runner::execute(const InputHitContainer &hits) {
  // Clear the storage.
  m_store.clear();

  // Store the input hits.
  m_store.add(m_config.digitization.inputSimHits, InputHitContainer{hits});

  // Log the input hits.
  logInput();

  // Run the track finding pipeline.
  Digitization digitization(m_config.digitization, m_level);
  SpacePointMaking spacePointMaking(m_config.spacePointMaking, m_level);
  TrackSeeding trackSeeding(m_config.trackSeeding, m_level);
  TrackEstimation trackEstimation(m_config.trackEstimation, m_level);
  TrackFinding trackFinding(m_config.trackFinding, m_level);

  digitization.execute(m_context);
  spacePointMaking.execute(m_context);
  trackSeeding.execute(m_context);
  trackEstimation.execute(m_context);
  trackFinding.execute(m_context);

  // Log the output tracks.
  logOutput();

  return m_context.eventStore.get<TrajectoriesContainer>(
      m_config.trackFinding.outputTrajectories);
}

void Runner::logInput() const {
  const auto &hits = m_context.eventStore.get<InputHitContainer>(
      m_config.digitization.inputSimHits);

  std::unordered_map<int, ProtoTrack> tracks;
  tracks.reserve(hits.size());

  size_t ihit = 0;
  for (const auto &hit : hits) {
    logHit(ihit, hit);
    tracks[hit.trackId].push_back(ihit++);
  }

  for (const auto &[trackId, track] : tracks) {
    logTrack("Real track", trackId, track);
  }
}

void Runner::logOutput() const {
  const auto &hits = m_context.eventStore.get<InputHitContainer>(
      m_config.digitization.inputSimHits);
  const auto &protos = m_context.eventStore.get<ProtoTrackContainer>(
      m_config.trackSeeding.outputProtoTracks);
  const auto &tracks = m_context.eventStore.get<ProtoTrackContainer>(
      m_config.trackFinding.outputTrackCandidates);

  logTracks("Proto track", protos);
  logTracks("Found track", tracks);
  checkTracks(hits, tracks);
}

void Runner::logHit(size_t hitId, const InputHit &hit) const {
  ACTS_DEBUG("Hit " << hitId << ": " << hit);
}

void Runner::logHits(const InputHitContainer &hits) const {
  size_t ihit = 0;
  for (const auto &hit : hits) {
    logHit(ihit++, hit);
  }
}

void Runner::logTrack(const std::string &prefix,
                      size_t trackId,
                      const ProtoTrack &track) const {
  std::stringstream out;
  for (auto ihit : track) {
    out << ihit << " ";
  }
  ACTS_DEBUG(prefix << " " << trackId << ": " << out.str());
}

void Runner::logTracks(const std::string &prefix,
                       const ProtoTrackContainer &tracks) const {
  size_t itrack = 0;
  for (const auto &track : tracks) {
    logTrack(prefix, itrack++, track);
  }
}

void Runner::checkTrack(const InputHitContainer &hits,
                        size_t trackId,
                        const ProtoTrack &track) const {
  std::unordered_map<int, size_t> counts;
  counts.reserve(track.size());

  for (auto ihit : track) {
    counts[hits.at(ihit).trackId]++;
  }

  size_t maxCount = 0;
  size_t realTrack = 0;
  for (const auto &[id, count] : counts) {
    if (count > maxCount) {
      maxCount = count;
      realTrack = id;
    }
  }

  double percent = (100. * maxCount) / track.size();

  ACTS_DEBUG("Most frequent ID for " << trackId << ": "
             << realTrack << " (" << percent << "%)");
}

void Runner::checkTracks(const InputHitContainer &hits,
                         const ProtoTrackContainer &tracks) const {
  size_t itrack = 0;
  for (const auto &track : tracks) {
    checkTrack(hits, itrack++, track);
  }
}

} // namespace Mpd::Tpc

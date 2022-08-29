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

  std::unordered_map<int, std::vector<size_t>> tracks;
  tracks.reserve(hits.size());

  size_t ihit = 0;
  for (const auto &hit : hits) {
    ACTS_VERBOSE(hit);
    tracks[hit.trackId].push_back(ihit++);
  }

  for (const auto &[trackId, track] : tracks) {
    std::stringstream out;
    for (auto ihit : track) {
      out << ihit << " ";
    }

    ACTS_DEBUG("Track " << trackId << ": " << out.str());
  }

  ACTS_DEBUG("Taken " << hits.size() << " hits");
}

void Runner::logOutput() const {
  const auto &hits = m_context.eventStore.get<InputHitContainer>(
      m_config.digitization.inputSimHits);
  const auto &tracks = m_context.eventStore.get<ProtoTrackContainer>(
      m_config.trackFinding.outputTrackCandidates);

  for (const auto &track : tracks) {
    std::stringstream out;
    for (auto ihit : track) {
      out << ihit << " ";
    }

    ACTS_DEBUG("Track " << out.str());
  }

  ACTS_DEBUG("Found " << tracks.size() << " tracks");
}

} // namespace Mpd::Tpc

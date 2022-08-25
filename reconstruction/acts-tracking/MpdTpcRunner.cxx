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
#include <unordered_set>

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

size_t Runner::countSimTracks(const InputHitContainer &hits) const {
  std::unordered_set<int> trackIds;
  trackIds.reserve(hits.size());

  for (const auto &hit : hits) {
    trackIds.insert(hit.trackId);
  }

  return trackIds.size();
}

void Runner::logInput() const {
  const auto &hits = m_context.eventStore.get<InputHitContainer>(
      m_config.digitization.inputSimHits);

  for (const auto &hit : hits) {
    ACTS_VERBOSE(hit);
  }

  ACTS_DEBUG("Input: " << hits.size() << " hits");
}

void Runner::logOutput() const {
  const auto &hits = m_context.eventStore.get<InputHitContainer>(
      m_config.digitization.inputSimHits);
  const auto &trajectories = m_context.eventStore.get<TrajectoriesContainer>(
      m_config.trackFinding.outputTrajectories);

  size_t nMultiTrajectories = 0;
  size_t nTracks = 0;

  for (const auto &trajectory : trajectories) {
    const auto &tips = trajectory.tips();

    if (tips.empty()) {
      continue;
    }

    ACTS_VERBOSE("Multi-trajectory: " << nMultiTrajectories);
    nMultiTrajectories++;

    for (const auto itrack : tips) {
      ACTS_VERBOSE("Track parameters:\n" << trajectory.trackParameters(itrack));
      nTracks++;
    }

    std::stringstream out;

    const auto &fittedStates = trajectory.multiTrajectory();
    for (size_t istate = 0; istate < fittedStates.size(); istate++) {
      const auto &state = fittedStates.getTrackState(istate);

      if (!state.hasCalibrated() && !state.hasUncalibrated()) {
        continue;
      }

      const auto &sourceLink = state.hasCalibrated()
          ? static_cast<const SourceLink&>(state.calibratedSourceLink())
          : static_cast<const SourceLink&>(state.uncalibrated());

      out << sourceLink.index() << " ";
    }

    ACTS_VERBOSE("States: " << out.str());
  }

  ACTS_DEBUG("Input: " << hits.size() << " points, "
                      << countSimTracks(hits) << " tracks");
  ACTS_DEBUG("Found: " << nMultiTrajectories << " multi-trajectories, "
                      << nTracks << " tracks");
}

} // namespace Mpd::Tpc

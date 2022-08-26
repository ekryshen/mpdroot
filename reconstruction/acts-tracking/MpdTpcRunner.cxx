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
  const auto &results = m_context.eventStore.get<TrajectoriesContainer>(
      m_config.trackFinding.outputTrajectories);

  size_t nMultiTrajectories = 0;
  size_t nTrajectories = 0;

  for (const auto &trajectories : results) {
    // Last indices that identify valid trajectories.
    auto &lastIndices = trajectories.tips();

    if (lastIndices.empty()) {
      continue;
    }

    // Collection of track states of the multi-trajectory.
    auto &fittedStates = trajectories.multiTrajectory();

    ACTS_VERBOSE("Multi-trajectory: " << nMultiTrajectories);
    nMultiTrajectories++;

    for (auto lastIndex : lastIndices) {
      std::stringstream out;

      fittedStates.visitBackwards(lastIndex, [&](const auto &state) {
        if (state.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag)) {
          auto &sourceLink = static_cast<const SourceLink&>(state.uncalibrated());
          auto hitIndex = sourceLink.index();
          out << hitIndex << " ";
        } else {
          out << "* ";
        }
      });

      ACTS_VERBOSE("Trajectory (reversed): " << out.str());
      ACTS_VERBOSE("Trajectory parameters:\n"
          << trajectories.trackParameters(lastIndex));

      nTrajectories++;
    }
  }

  ACTS_DEBUG("Input: " << hits.size() << " simulation points, "
                       << countSimTracks(hits) << " simulation tracks");
  ACTS_DEBUG("Found: " << nMultiTrajectories << " multi-trajectories, "
                       << nTrajectories << " trajectories");
}

} // namespace Mpd::Tpc

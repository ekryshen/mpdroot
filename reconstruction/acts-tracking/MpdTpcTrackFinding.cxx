// This file is a part of the NICA project.
//
// Copyright (C) 2022 JINR

#include "MpdTpcContext.h"
#include "MpdTpcEventStorage.h"
#include "MpdTpcMagneticField.h"
#include "MpdTpcTrackFinding.h"

#include <Acts/Definitions/Units.hpp>
#include <Acts/EventData/TrackParameters.hpp>
#include <Acts/Geometry/TrackingGeometry.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <cassert>
#include <limits>
#include <vector>

namespace Mpd::Tpc {

TrackFinding::TrackFinding(Config config, Acts::Logging::Level level):
  Algorithm("TrackFinding", level),
  m_config(std::move(config)) {

  assert(!m_config.inputMeasurements.empty()
    && "Missing measurements input collection");
  assert(!m_config.inputSourceLinks.empty()
    && "Missing source links input collection");
  assert(!m_config.inputInitialTrackParameters.empty()
    && "Missing initial track parameters input collection");
  assert(!m_config.outputTrajectories.empty()
    && "Missing trajectories output collection");
}

ProcessCode TrackFinding::execute(const Context &context) const {
  using TrackParametersContainer = ActsExamples::TrackParametersContainer;
  using AccessorDelegate = Acts::SourceLinkAccessorDelegate<Iterator>;

  ACTS_DEBUG("Track finding");

  // Read input data.
  const auto &measurements =
      context.eventStore.get<MeasurementContainer>(m_config.inputMeasurements);
  const auto &sourceLinks =
      context.eventStore.get<Container>(m_config.inputSourceLinks);
  const auto &initialParameters =
      context.eventStore.get<TrackParametersContainer>(m_config.inputInitialTrackParameters);

  MeasurementCalibrator calibrator(measurements);
  Updater updater;
  Smoother smoother;
  Selector selector(m_config.measurementSelectorConfig);

  Extensions extensions;
  extensions.calibrator.connect
      <&MeasurementCalibrator::calibrate>(&calibrator);
  extensions.updater.connect
      <&Updater::operator()<Acts::VectorMultiTrajectory>>(&updater);
  extensions.smoother.connect
      <&Smoother::operator()<Acts::VectorMultiTrajectory>>(&smoother);
  extensions.measurementSelector.connect
      <&Selector::select<Acts::VectorMultiTrajectory>>(&selector);

  Accessor accessor;
  accessor.container = &sourceLinks;
  AccessorDelegate accessorDelegate;
  accessorDelegate.connect<&Accessor::range>(&accessor);

  Stepper stepper(m_config.magneticField);
  Navigator navigator(m_config.navigatorConfig);
  Propagator propagator(std::move(stepper), std::move(navigator));
  Filter trackFinder(std::move(propagator));

  Options options(
    context.gContext,
    context.mContext,
    context.cContext,
    accessorDelegate,
    extensions,
    Acts::LoggerWrapper{logger()},
    m_config.propagatorOptions,
    &(*m_config.referenceSurface),
    m_config.multipleScattering,
    m_config.energyLoss,
    m_config.smoothing
  );

  Results results = trackFinder.findTracks(initialParameters, options);

  // Compute shared hits from all the reconstructed tracks.
  if (m_config.computeSharedHits) {
    computeSharedHits(sourceLinks, results);
  }

  TrajectoriesContainer trajectories;
  trajectories.reserve(initialParameters.size());

  // Loop over the track finding results for all initial parameters.
  size_t trajectoryCount = 0;
  for (size_t iseed = 0; iseed < initialParameters.size(); ++iseed) {
    auto &result = results[iseed];

    if (result.ok()) {
      const auto &trajectory = result.value();

      const auto &fittedStates = trajectory.fittedStates;
      const auto &lastIndices = trajectory.lastMeasurementIndices;
      const auto &fittedParams = trajectory.fittedParameters;
 
      ACTS_DEBUG("Multi-trajectory of " << lastIndices.size() << "/"
                                        << fittedStates.size()
                                        << " trajectories/states "
                                        << " has been found");

      trajectories.emplace_back(fittedStates, lastIndices, fittedParams);
      trajectoryCount++;
    } else {
      ACTS_WARNING("Track finding failed for seed " << iseed << " with error "
                                                    << result.error());

      // Add an empty result so the output is of the same length as the input.
      trajectories.push_back(Trajectories());
    }
  }

  ACTS_DEBUG("Finalized track finding with " << trajectoryCount
                                             << " multi-trajectories");

  context.eventStore.add(m_config.outputTrajectories, std::move(trajectories));

  return ProcessCode::SUCCESS;
}

void TrackFinding::computeSharedHits(const Container &sourceLinks,
                                     Results &results) const {
  const auto nHits = sourceLinks.size();
  const auto invalid = std::numeric_limits<size_t>::max();

  std::vector<size_t> firstTrackOnHit(nHits, invalid);
  std::vector<size_t> firstStateOnHit(nHits, invalid);

  // Iterate over the seeds.
  for (size_t i = 0; i < results.size(); i++) {
    if (!results.at(i).ok()) {
      // No trajectory found for the given seed.
      continue;
    }

    // Trajectories associated w/ the given seed.
    auto &trajectories = results.at(i).value();
    // Last indices that identify valid trajectories.
    auto &lastIndices = trajectories.lastMeasurementIndices;
    // Collection of track states associated w/ the given seed.
    auto &fittedStates = trajectories.fittedStates;

    // Iterate over the valid trajectories.
    for (auto lastIndex : lastIndices) {
      fittedStates.visitBackwards(lastIndex, [&](const auto &state) {
        // Do nothing unless the state is a real measurement.
        if (!state.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag)) {
          return;
        }

        // Get the index in the measurement array.
        auto &sourceLink = static_cast<const SourceLink&>(state.uncalibrated());
        auto hitIndex = sourceLink.index();

        auto &firstTrackIndex = firstTrackOnHit.at(hitIndex);
        auto &firstStateIndex = firstStateOnHit.at(hitIndex);

        // Check if the hit has not been already used.
        if (firstTrackIndex == invalid) {
          firstTrackIndex = i;
          firstStateIndex = state.index();
          return;
        }

        // Mark the first and current states as shared.
        auto &firstStateFlags = results.at(firstTrackIndex)
            .value()
            .fittedStates.getTrackState(firstStateIndex)
            .typeFlags();

        auto &currentStateFlags = results.at(i)
            .value()
            .fittedStates.getTrackState(state.index())
            .typeFlags();

        firstStateFlags.set(Acts::TrackStateFlag::SharedHitFlag);
        currentStateFlags.set(Acts::TrackStateFlag::SharedHitFlag);
      });
    }
  }
}

} // namespace Mpd::Tpc

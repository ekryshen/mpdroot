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
      const auto &measurementIndices = trajectory.lastMeasurementIndices;
      const auto &fittedParams = trajectory.fittedParameters;
 
      ACTS_DEBUG("Trajectory of length " << fittedStates.size()
                                         << " has been found");

      trajectories.emplace_back(fittedStates, measurementIndices, fittedParams);
      trajectoryCount++;
    } else {
      ACTS_WARNING("Track finding failed for seed " << iseed << " with error "
                                                    << result.error());
      // Add an empty result so the output is of the same length as the input.
      trajectories.push_back(Trajectories());
    }
  }

  ACTS_DEBUG("Finalized track finding with " << trajectoryCount
                                             << " track candidates");

  context.eventStore.add(m_config.outputTrajectories, std::move(trajectories));

  return ProcessCode::SUCCESS;
}

void TrackFinding::computeSharedHits(const Container &sourceLinks,
                                     Results &results) const {
  std::vector<size_t> firstTrackOnTheHit(sourceLinks.size(),
                                         std::numeric_limits<size_t>::max());
  std::vector<size_t> firstStateOnTheHit(sourceLinks.size(),
                                         std::numeric_limits<size_t>::max());

  for (size_t i = 0; i < results.size(); i++) {
    if (!results.at(i).ok()) {
      continue;
    }

    auto &result = results.at(i).value();
    auto &indices = result.lastMeasurementIndices;

    for (auto index : indices) {
      result.fittedStates.visitBackwards(index, [&](const auto &state) {
        if (!state.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag)) {
          return;
        }

        auto &sourceLink = static_cast<const SourceLink&>(state.uncalibrated());
        size_t hitIndex = sourceLink.index();

        // Check if the hit has not been already used.
        if (firstTrackOnTheHit.at(hitIndex) == std::numeric_limits<size_t>::max()) {
          firstTrackOnTheHit.at(hitIndex) = i;
          firstStateOnTheHit.at(hitIndex) = state.index();
          return;
        }

        // Check if the first track state has been marked as shared.
        int indexFirstTrack = firstTrackOnTheHit.at(hitIndex);
        int indexFirstState = firstStateOnTheHit.at(hitIndex);

        if (!results.at(indexFirstTrack)
                .value()
                .fittedStates.getTrackState(indexFirstState)
                .typeFlags()
                .test(Acts::TrackStateFlag::SharedHitFlag)) {
          results.at(indexFirstTrack)
              .value()
              .fittedStates.getTrackState(indexFirstState)
              .typeFlags()
              .set(Acts::TrackStateFlag::SharedHitFlag);
        }

        // Decorate this track.
        results.at(i)
            .value()
            .fittedStates.getTrackState(state.index())
            .typeFlags()
            .set(Acts::TrackStateFlag::SharedHitFlag);
      });
    }
  }
}

} // namespace Mpd::Tpc

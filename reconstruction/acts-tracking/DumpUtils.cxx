// This file is a part of the NICA project.
//
// Copyright (C) 2022 JINR

#include "DumpUtils.h"

#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <map>

namespace Mpd::Tpc {

void dumpTracks(
    const ActsExamples::AlgorithmContext &context,
    Dump::Results &results,
    std::string spacePointsID,
    std::string outPath) {

  size_t eventNumber = context.eventNumber;
  auto fname = outPath + "/" +
      "event_" + std::to_string(eventNumber) + "_prototracks.txt";
  std::ofstream fout;

  fout.open(fname);
  const auto &spacePoints = context.eventStore.get<SpacePointContainer>(
      spacePointsID);

  std::map<Int_t, Int_t> hitIndexToSpIndexMap;

  Int_t spacePointIndex = 0;
  for (const auto &spacePoint : spacePoints) {
    Int_t hitIndex = spacePoint.measurementIndex();
    hitIndexToSpIndexMap[hitIndex] = spacePointIndex++;
  }

  fout << "# format: "
      "(x, y, z, hit-index, phi, theta, q/p, t, chi2, seed-index)+" <<
      std::endl;
  std::cout << fname << " has been created" << std::endl;

  // Iterate over the seeds.
  for (size_t itrack = 0; itrack < results.size(); itrack++) {
    if (!results.at(itrack).ok()) {
      // No trajectory found for the given seed.
      continue;
    }

    // Trajectories associated w/ the given seed.
    auto &trajectories = results.at(itrack).value();
    // Last indices that identify valid trajectories.
    auto &lastIndices = trajectories.lastMeasurementIndices;
    // Collection of track states associated w/ the given seed.
    auto &fittedStates = trajectories.fittedStates;

    // Iterate over the valid trajectories.
    for (auto lastIndex : lastIndices) {
      auto firstPass = true;
      fittedStates.visitBackwards(lastIndex, [&](const auto &state) {
        // Do nothing unless the state is a real measurement.
        if (!state.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag)) {
          return;
        }

        auto &sourceLink = static_cast<const SourceLink&>(state.uncalibrated());
        Int_t ihit = sourceLink.index();
        Int_t iSpacePoint = hitIndexToSpIndexMap.at(ihit);
        auto spacePoint = spacePoints.at(iSpacePoint);

        Double_t x = spacePoint.x();
        Double_t y = spacePoint.y();
        Double_t z = spacePoint.z();

        auto params = state.parameters();

        Double_t phi    = params[2];
        Double_t theta  = params[3];
        Double_t qOverP = params[4];
        Double_t t      = params[5];
        Double_t chi2   = state.chi2();
        size_t seedIdx  = itrack;

        if (!firstPass) {
          fout << ", ";
        }
        firstPass = false;
        fout << x        << ", " <<
                y        << ", " <<
                z        << ", " <<
                ihit     << ", " <<
                phi      << ", " <<
                theta    << ", " <<
                qOverP   << ", " <<
                t        << ", " <<
                chi2     << ", " <<
                seedIdx;
      });

      fout << std::endl;
    }
  }
}

void dumpHits(
    const Mpd::Tpc::InputHitContainer &hits,
    Int_t eventNumber,
    std::string outPath) {
  auto fname = outPath + "/event_" + std::to_string(eventNumber) + "_hits.txt";
  std::ofstream fout(fname);
  fout << "# format: x, y, z, (trackId)+" << std::endl;
  std::cout << fname << " has been created" << std::endl;

  for (const auto &hit : hits) {
    fout << hit.position[0] << ", " <<
            hit.position[1] << ", " <<
            hit.position[2] << ", " <<
            hit.trackId <<
            std::endl;
  }
}

} // namespace Mpd::Tpc

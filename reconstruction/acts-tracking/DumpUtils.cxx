// This file is a part of the NICA project.
//
// Copyright (C) 2023 JINR

#include "DumpUtils.h"
#include "MpdMCTrack.h"

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

  fout << "# format: (space-point-index, phi, theta, q/p, t, chi2)+" <<
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

        auto params = state.parameters();

        Double_t phi    = params[2];
        Double_t theta  = params[3];
        Double_t qOverP = params[4];
        Double_t t      = params[5];
        Double_t chi2   = state.chi2();

        if (!firstPass) {
          fout << ", ";
        }
        firstPass = false;
        fout << iSpacePoint << ", " <<
                phi         << ", " <<
                theta       << ", " <<
                qOverP      << ", " <<
                t           << ", " <<
                chi2;
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

void dumpSpacePoints(
    const ActsExamples::AlgorithmContext &context,
    std::string spacePointsID,
    std::string hitsID,
    std::string outPath) {

  size_t eventNumber = context.eventNumber;
  auto fname = outPath + "/" +
      "event_" + std::to_string(eventNumber) + "_space_points.txt";
  std::ofstream fout(fname);
  fout << "# format: x, y, z, (trackId)+" << std::endl;
  std::cout << fname << " has been created" << std::endl;

  const auto &spacePoints = context.eventStore.get<SpacePointContainer>(
      spacePointsID);
  const auto &hits = context.eventStore.get<InputHitContainer>(hitsID);

  for (const auto &spacePoint : spacePoints) {
    Double_t x = spacePoint.x();
    Double_t y = spacePoint.y();
    Double_t z = spacePoint.z();
    Int_t ihit = spacePoint.measurementIndex();
    auto &hit = hits.at(ihit);
    Int_t trackId = hit.trackId;

    fout << x       << ", " <<
            y       << ", " <<
            z       << ", " <<
            trackId <<
            std::endl;
  }
}

// Return map trackId -> the number of corresponding hits.
std::map<Int_t, Int_t> calcTrackIdToNHits(
    const InputHitContainer &hits,
    TClonesArray *mcTracks) {
  std::map<Int_t, Int_t> trackIdToNHitsMap;
  Int_t nMCTracks = mcTracks->GetEntriesFast();
  for (Int_t trackId = 0; trackId < nMCTracks; trackId++) {
    trackIdToNHitsMap[trackId] = 0;
  }
  for (const auto &hit : hits) {
    Int_t trackId = hit.trackId;
    trackIdToNHitsMap.at(trackId)++;
  }
  return trackIdToNHitsMap;
}

void dumpTrackIds(
    TClonesArray *mcTracks,
    Int_t eventNumber,
    std::string outPath,
    const InputHitContainer &hits,
    Bool_t dumpNHits) {

  auto fname = outPath + "/" +
      "event_" + std::to_string(eventNumber) + "_trackIds.txt";
  std::ofstream fout(fname);

  fout << "# format: trackId, primary-1-secondary-0";
  if (dumpNHits) {
    fout << ", nHits";
  }
  fout << std::endl;

  std::cout << fname << " has been created" << std::endl;

  std::map<Int_t, Int_t> trackIdToNHitsMap;
  if (dumpNHits) {
    trackIdToNHitsMap = calcTrackIdToNHits(hits, mcTracks);
  }

  Int_t nMC = mcTracks->GetEntriesFast();
  for (size_t trackId = 0; trackId < nMC; trackId++) {
    auto track = static_cast<MpdMCTrack*>(mcTracks->UncheckedAt(trackId));
    Int_t motherId = track->GetMotherId();
    Bool_t isPrimary = (motherId == -1) ? 1 : 0;
    fout << trackId << ", " << isPrimary;

    if (dumpNHits) {
      Int_t nhits = trackIdToNHitsMap.at(trackId);
      fout << ", " << nhits;
    }
    fout << std::endl;
  }
}

} // namespace Mpd::Tpc

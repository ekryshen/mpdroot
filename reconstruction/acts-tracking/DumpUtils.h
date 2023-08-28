#pragma once

#include "MpdTpcEventData.h"
#include "MpdTpcInputHit.h"

#include "ActsExamples/Framework/AlgorithmContext.hpp"

#include <Acts/TrackFinding/CombinatorialKalmanFilter.hpp>

#include <Rtypes.h>
#include <TClonesArray.h>

namespace Mpd::Tpc {

namespace Dump {

using Container = SourceLinkContainer;
using Result = Acts::CombinatorialKalmanFilterResult<Trajectory>;
using Results = std::vector<Acts::Result<Result>>;

} // namespace Dump

void dumpTracks(
    const ActsExamples::AlgorithmContext &context,
    Dump::Results &results,
    std::string spacePointsID,
    std::string outPath);

void dumpHits(const Mpd::Tpc::InputHitContainer &hits,
              Int_t eventNumber,
              std::string outPath);

void dumpSpacePoints(
    const ActsExamples::AlgorithmContext &context,
    std::string spacePointsID,
    std::string hitsID,
    std::string outPath);

void dumpTrackIds(
    TClonesArray *mcTracks,
    Int_t eventNumber,
    std::string outPath);
} // namespace Mpd::Tpc

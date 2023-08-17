#pragma once

#include "MpdTpcEventData.h"
#include "MpdTpcInputHit.h"

#include "ActsExamples/Framework/AlgorithmContext.hpp"

#include <Acts/TrackFinding/CombinatorialKalmanFilter.hpp>

#include <Rtypes.h>

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

} // namespace Mpd::Tpc

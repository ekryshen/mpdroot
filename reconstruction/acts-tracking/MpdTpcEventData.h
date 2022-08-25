// This file is a part of the NICA project.
//
// Copyright (C) 2022 JINR

#pragma once

#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"

#include <Acts/Seeding/Seed.hpp>

#include <vector>

namespace Mpd::Tpc {

using SpacePoint = ActsExamples::SimSpacePoint;
using SpacePointContainer = ActsExamples::SimSpacePointContainer;

using MeasurementCalibrator = ActsExamples::MeasurementCalibrator;
using MeasurementContainer = ActsExamples::MeasurementContainer;

using Index = ActsExamples::Index;

using SourceLink = ActsExamples::IndexSourceLink;
using SourceLinkAccessor = ActsExamples::IndexSourceLinkAccessor;
using SourceLinkIterator = ActsExamples::IndexSourceLinkAccessor::Iterator;
using SourceLinkContainer = ActsExamples::IndexSourceLinkContainer;

using SeedContainer = std::vector<Acts::Seed<SpacePoint>>;

using ProtoTrack = ActsExamples::ProtoTrack;
using ProtoTrackContainer = ActsExamples::ProtoTrackContainer;

using Trajectory = Acts::VectorMultiTrajectory;
using Trajectories = ActsExamples::Trajectories;
using TrajectoriesContainer = ActsExamples::TrajectoriesContainer;

} // namespace Mpd::Tpc

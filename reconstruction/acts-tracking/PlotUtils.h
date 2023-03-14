// This file is a part of the NICA project.
//
// Copyright (C) 2022 JINR

#pragma once

enum CoordSystem {
  XY,
  ZY,
  ZR
};

void buildHistograms(
    const Mpd::Tpc::Statistics &statistics,
    Int_t nTracks,
    Int_t eventCounter,
    std::string outPath);

void plotQualityOnP(
    const Mpd::Tpc::InputHitContainer &hits,
    const Mpd::Tpc::ProtoTrackContainer &trajectories,
    Int_t eventCounter,
    std::string outPath);

void plotOutputTracks(
    Int_t canvasX,
    Int_t canvasY,
    const std::shared_ptr<const Acts::TrackingGeometry> &geometry,
    const Mpd::Tpc::SpacePointContainer &spacePoints,
    const Mpd::Tpc::InputHitContainer &hits,
    const Mpd::Tpc::ProtoTrackContainer &trajectories,
    Int_t eventCounter,
    std::string outPath,
    bool multicoloured,
    Int_t lineWidth,
    CoordSystem coordSystem);

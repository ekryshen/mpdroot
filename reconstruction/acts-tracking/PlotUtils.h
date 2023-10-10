// This file is a part of the NICA project.
//
// Copyright (C) 2022 JINR

#pragma once

enum Projection {
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

void plotRealTracks(
    Int_t canvasX,
    Int_t canvasY,
    const std::shared_ptr<const Acts::TrackingGeometry> &geometry,
    const Mpd::Tpc::InputHitContainer &hits,
    Int_t eventCounter,
    std::string outPath,
    Bool_t multicoloured,
    Int_t lineWidth,
    Projection projection,
    Double_t zmin,
    Double_t zmax,
    Double_t rmax,
    Bool_t plotLines = false,
    std::string namePostfix = "",
    Bool_t grid = false,
    Bool_t realAspectRatio = false,
    // Size in pixels
    Int_t txtSize = 15,
    // Print trackId label for every txtStep point
    Int_t txtStep = 10);

void plotOutputTracks(
    Int_t canvasX,
    Int_t canvasY,
    const std::shared_ptr<const Acts::TrackingGeometry> &geometry,
    const Mpd::Tpc::SpacePointContainer &spacePoints,
    const Mpd::Tpc::InputHitContainer &hits,
    const Mpd::Tpc::ProtoTrackContainer &trajectories,
    Int_t eventCounter,
    std::string outPath,
    Bool_t multicoloured,
    Int_t lineWidth,
    Projection projection,
    std::string namePostfix = "",
    Bool_t plotLabels = true,
    Int_t txtSize = 15,
    Int_t txtStep = 10);


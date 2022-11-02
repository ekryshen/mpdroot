// This file is a part of the NICA project.
//
// Copyright (C) 2022 JINR

#pragma once

void buildHistograms(const Mpd::Tpc::Statistics &statistics,
                     const int nTracks,
                     const int eventCounter);

void drawQualityOnP(const Mpd::Tpc::InputHitContainer &hits,
                    const Mpd::Tpc::ProtoTrackContainer &trajectories,
                    const int eventCounter);

void plotOutputTracks(const int canvasX,
                      const int canvasY,
                      const std::shared_ptr<const Acts::TrackingGeometry> &geometry,
                      const Mpd::Tpc::SpacePointContainer &spacePoints,
                      const Mpd::Tpc::InputHitContainer &hits,
                      const Mpd::Tpc::ProtoTrackContainer &trajectories,
                      const int eventCounter,
                      const bool multicoloured);

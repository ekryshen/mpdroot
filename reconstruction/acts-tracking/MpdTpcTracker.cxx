// This file is a part of the NICA project.
//
// Copyright (C) 2022 JINR

#include "MpdTpcTracker.h"

#include "MpdCodeTimer.h"
#include "MpdKalmanHit.h"
#include "MpdTpcKalmanTrack.h"
#include "TpcPoint.h"

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/Units.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <iostream>

MpdTpcTracker::MpdTpcTracker(const char *title):
    FairTask(title) {}

MpdTpcTracker::~MpdTpcTracker() {}

inline TClonesArray *getArray(const char *name) {
  auto *manager = FairRootManager::Instance();
  auto *array = static_cast<TClonesArray*>(manager->GetObject(name));

  if (!array) {
    std::cout << "[MpdTpcTracker::Init]: No " << name << " array!" << std::endl;
  }

  return array;
}

InitStatus MpdTpcTracker::Init() {
  std::cout << "[MpdTpcTracker::Init]: Started" << std::endl;

  // Geometry must already be loaded (gGeoManager != nullptr).
  fRunner = std::make_unique<Mpd::Tpc::Runner>(Acts::Logging::DEBUG);

  fPoints = getArray("TpcPoint");
//  fKalmanHits = getArray("MpdKalmanHit");
//  fKalmanTracks = getArray("MpdTpcKalmanTrack");

  if (!fPoints /*|| !fKalmanHits || !fKalmanTracks*/) {
    return kERROR;
  }

  std::cout << "[MpdTpcTracker::Init]: Finished" << std::endl;
  return kSUCCESS;
}

InitStatus MpdTpcTracker::ReInit() {
  std::cout << "[MpdTpcTracker::ReInit]: Do nothing" << std::endl;
  return kSUCCESS;
}

void MpdTpcTracker::Exec(Option_t *option) {
  std::cout << "[MpdTpcTracker::Exec]: Started" << std::endl;

  if (MpdCodeTimer::Active()) {
    MpdCodeTimer::Instance()->Start(Class()->GetName(), __FUNCTION__);
  }

  // Convert the input TPC points.
  const auto nPoints = fPoints->GetEntriesFast();

  Mpd::Tpc::InputHitContainer hits;
  hits.reserve(nPoints);

  const auto lengthScalor = Acts::UnitConstants::cm / Acts::UnitConstants::mm;

  for (Int_t i = 0; i < nPoints; i++) {
    const auto *point = static_cast<TpcPoint*>(fPoints->UncheckedAt(i));

    hits.emplace_back(Mpd::Tpc::InputHit{
        point->GetTrackID(),
        point->GetDetectorID(),
        Acts::Vector3{
            point->GetX() * lengthScalor,
            point->GetY() * lengthScalor,
            point->GetZ() * lengthScalor
        },
        Acts::Vector3{
            point->GetPx(),    // TODO: Scaling
            point->GetPy(),    // TODO: Scaling
            point->GetPz()     // TODO: Scaling
        },
        point->GetTime(),      // TODO: Scaling
        point->GetLength() * lengthScalor,
        point->GetEnergyLoss() // TODO: Scaling
    });
  }

  // Run the track finding algorithm.
  const auto &trajectories = fRunner->execute(hits);

  // Convert the output tracks.
  /* base
  track.fID;
  track.fNhits;
  track.fTrackDir;
  track.fTrackType;
  track.fLastLay;
  track.fNofWrong;
  track.fNode;
  track.fNodeNew;
  track.fPartID;
  track.fPos;
  track.fPosNew;
  track.fPosAtHit;
  track.fChi2;
  track.fChi2Vertex;
  track.fLength;
  track.fLengAtHit;
  track.fParam = nullptr;
  track.fParamNew    = nullptr;
  track.fParamAtHit  = nullptr;
  track.fCovar       = nullptr;
  track.fWeightAtHit = nullptr;
  track.fVertex      = nullptr;
  track.fHits        = new TObjArray(70);
  track.fStepMap;
  track.fFlag;
  track.fWeight = nullptr;
  fHits->Add(track.fHits->UncheckedAt(i));

  Int_t              GetNofTrHits() const { return fTrHits->GetEntriesFast(); } ///< get number of track hits
  TClonesArray      *GetTrHits() const { return fTrHits; }  // MpdKalmanHit
  */

  /*
  auto nTracks = 0;
  for (const auto &trajectory : trajectories) {
    if (trajectory.tips().empty()) {
      continue;
    }

    auto *track = new ((*fKalmanTracks)[nTracks++]) MpdTpcKalmanTrack();
    const auto &fittedStates = trajectory.multiTrajectory();

    for (size_t i = 0; i < fittedStates.size(); i++) {
      auto *hit = new ((*track->GetTrHits())[i]) MpdKalmanHit();
      const auto &state = fittedStates.getTrackState(i);
      const auto index = static_cast<const Mpd::Tpc::SourceLink&>(
          state.uncalibrated()).index();
      // TODO: Fill hit
    }
  }
  */ 
  if (MpdCodeTimer::Active()) {
    MpdCodeTimer::Instance()->Stop(Class()->GetName(), __FUNCTION__);
  }

  std::cout << "[MpdTpcTracker::Exec]: Finished" << std::endl;
}

void MpdTpcTracker::Finish() {
  std::cout << "[MpdTpcTracker::Finish]: Do nothing" << std::endl;
}

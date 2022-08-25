// This file is a part of the NICA project.
//
// Copyright (C) 2022 JINR

#ifndef MPDTPCTRACKER_H
#define MPDTPCTRACKER_H

#include "MpdTpcRunner.h"

#include <FairTask.h>

#include <TClonesArray.h>

/// @brief Acts-based track finder for TPC.
class MpdTpcTracker final : public FairTask {
public:
  static constexpr auto Title = "TPC Acts-based tracker";
   
  explicit MpdTpcTracker(const char *title = Title);
  virtual ~MpdTpcTracker();

  InitStatus Init() override;
  InitStatus ReInit() override;
  void Exec(Option_t *option) override;
  void Finish() override;

private:
  Mpd::Tpc::Runner fRunner;

  TClonesArray *fPoints;
  TClonesArray *fKalmanHits;
  TClonesArray *fKalmanTracks;

  ClassDef(MpdTpcTracker, 0);
};

#endif // MPDTPCTRACKER_H

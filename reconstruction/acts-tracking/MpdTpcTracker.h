// This file is a part of the NICA project.
//
// Copyright (C) 2022 JINR

#ifndef MPDTPCTRACKER_H
#define MPDTPCTRACKER_H

#include "MpdTpcRunner.h"

#include "BaseTpcSectorGeo.h"

#include <FairTask.h>
#include <TClonesArray.h>

#include <memory>

class BaseTpcSectorGeo;

/// @brief Acts-based track finder for TPC.
class MpdTpcTracker final : public FairTask {
public:
  static constexpr auto TaskTitle  = "TPC Acts-based tracker";
  static constexpr auto UseMcHits  = false;
  static constexpr auto PlotGraphs = true;
   
  explicit MpdTpcTracker(const BaseTpcSectorGeo &secGeo,
                         const char *title = TaskTitle):
      FairTask(title), fSecGeo(secGeo) {}

  virtual ~MpdTpcTracker() {}

  InitStatus Init() override;
  InitStatus ReInit() override;
  void Exec(Option_t *option) override;
  void Finish() override;

private:
  std::unique_ptr<Mpd::Tpc::Runner> fRunner;

  TClonesArray *fPoints;
  TClonesArray *fHits;
  TClonesArray *fTracks;

  const BaseTpcSectorGeo &fSecGeo;

  ClassDef(MpdTpcTracker, 0);
};

#endif // MPDTPCTRACKER_H

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
                         std::string outPath,
                         const char *title = TaskTitle):
      FairTask(title), fSecGeo(secGeo), fOutPath(std::move(outPath)) {}

  virtual ~MpdTpcTracker() {
    delete fTracks; // TODO: investigate Clear()
}

  InitStatus Init() override;
  InitStatus ReInit() override;
  void Exec(Option_t *option) override;
  void Finish() override;

private:
  std::unique_ptr<Mpd::Tpc::Runner> fRunner;

  TClonesArray *fPoints;
  TClonesArray *fHits;
  /// Needed to collect performance statistics.
  TClonesArray *fMCTracks;

  TClonesArray *fTracks;

  TEfficiency *fEffPt;    // Add delete memory
  TEfficiency *fEffEta;   // Add delete memory

  Int_t fNTruth;
  Int_t fNFake;
  Int_t fNRealTracks;

  const BaseTpcSectorGeo &fSecGeo;
  std::string fOutPath;
  ActsExamples::CKFPerformanceWriter *fPerfWriter;

  ClassDef(MpdTpcTracker, 0);
};

#endif // MPDTPCTRACKER_H

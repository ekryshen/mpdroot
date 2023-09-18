// This file is a part of the NICA project.
//
// Copyright (C) 2022 JINR

#pragma once

#include "MpdTpcDetector.h"
#include "MpdTpcDigitization.h"
#include "MpdTpcMagneticField.h"
#include "MpdTpcParticleSelector.h"
#include "MpdTpcSpacePointMaking.h"
#include "MpdTpcTrackEstimation.h"
#include "MpdTpcTrackFinding.h"
#include "MpdTpcTrackSeeding.h"

#include "ActsExamples/Io/Performance/CKFPerformanceWriter.hpp"
#include "ActsExamples/TruthTracking/TruthSeedSelector.hpp"

#include <Acts/Definitions/Units.hpp>
#include <Acts/Utilities/Logger.hpp>

#include <algorithm>
#include <memory>
#include <string>

using namespace Acts::UnitLiterals;

namespace Mpd::Tpc {

/// @brief Configuration of the Acts-based tracker.
struct Config final {

  //===--------------------------------------------------------------------===//
  // Objects identifiers
  //===--------------------------------------------------------------------===//

  static constexpr auto SimHitsID           = "simhits";
  static constexpr auto SelectedSimHitsID   = "selectedsimhits";
  static constexpr auto SourceLinksID       = "sourcelinks";
  static constexpr auto MeasurementsID      = "measurements";
  static constexpr auto SpacePointsID       = "spacepoints";
  static constexpr auto SeedsID             = "seeds";
  static constexpr auto ProtoTracksID       = "prototracks";
  static constexpr auto EstTrackParamsID    = "estimatedparameters";
  static constexpr auto EstProtoTracksID    = "estimatedprototracks";
  static constexpr auto TrajectoriesID      = "trajectories";
  static constexpr auto TrackCandidatesID   = "trackcandidates";
  static constexpr auto ParticlesID         = "inputparticles";
  static constexpr auto SelectedParticlesID = "selectedparticles";
  static constexpr auto HitParticlesMapID   = "hitparticlesmap";
  static constexpr auto SelectedHitParticlesMapID = "selectedhitparticlesmap";

  //===--------------------------------------------------------------------===//
  // Common options (for selector and for seeding)
  //===--------------------------------------------------------------------===//

  static constexpr auto Rmin                =  Detector::Rmin;      // ~ 0.4 m
  static constexpr auto Rmax                =  Detector::Rmax;      // ~ 1.4 m
  static constexpr auto Zmin                =  Detector::Zmin;      // ~-1.7 m
  static constexpr auto Zmax                =  Detector::Zmax;      // ~ 1.7 m
  static constexpr auto EtaMax              =  2.5;
  static constexpr auto PtMin               =  0._GeV;              // 0.02 < Pt < 10 GeV

  //===--------------------------------------------------------------------===//
  // Particle selector
  //===--------------------------------------------------------------------===//

  static constexpr auto SelectorEnabled     =  false;
  static constexpr auto PrimaryParticlesOnly = false;
  // Truth particle kinematic cuts.
  static constexpr auto SelectorRmin        =  Rmin;
  static constexpr auto SelectorRmax        =  Rmax;
  static constexpr auto SelectorZmin        =  Zmin;
  static constexpr auto SelectorZmax        =  Zmax;
  static constexpr auto SelectorPhiMin      = -3.15;
  static constexpr auto SelectorPhiMax      =  3.15;
  static constexpr auto SelectorEtaMin      =  std::numeric_limits<Double_t>::lowest();
  static constexpr auto SelectorEtaMax      =  std::numeric_limits<Double_t>::max();
  static constexpr auto AbsEtaMin           =  0;
  static constexpr auto AbsEtaMax           =  EtaMax;
  static constexpr auto SelectorPtMin       =  PtMin;
  static constexpr auto SelectorPtMax       =  10._GeV;
  static constexpr auto KeepNeutral         =  false;

  /// Requirement on number of recorded hits.
  static constexpr auto NHitsMin            =  9;
  static constexpr auto NHitsMax            =  std::numeric_limits<size_t>::max();

  //===--------------------------------------------------------------------===//
  // Track seeding
  //===--------------------------------------------------------------------===//

  static constexpr auto SeedRmin            =  Rmin;
  static constexpr auto SeedRmax            =  Rmax;
  static constexpr auto SeedZmin            =  Zmin;
  static constexpr auto SeedZmax            =  Zmax;
  static constexpr auto CollisionZmin       = -30._cm;              // Close to 0
  static constexpr auto CollisionZmax       =  30._cm;              // Close to 0
  Double_t CotThetaMax;                                             // = 1/2 * (e^eta - e^(-eta)); eta = EtaMax
                                                                    // CotThetaMax = max(dZ/dR)=1.7 ~ 1.3 eta (eta < 1.2)
  static constexpr auto SeedBinSizeR        =  10._mm;              // 10._mm for MC (pads are ~12-18 mm)
  static constexpr auto SeedDeltaRmin       =  10._mm;              // 02._mm for MC
  static constexpr auto SeedDeltaRmax       =  60._mm;              // 20._mm for MC
  static constexpr auto SeedDeltaZmax       =  20._cm;              // FIXME 10._cm for MC
  static constexpr auto MaxSeedsPerSpM      =  3;                   // FIXME
  static constexpr auto SigmaScattering     =  5;                   // FIXME
  static constexpr auto MaxPtScattering     =  5._GeV;              // Max Pt for scattering
  static constexpr auto RadLengthPerSeed    =  0.05;                // OK
  static constexpr auto SeedPtMin           =  PtMin;
  static constexpr auto Bz                  =  MagneticField::Bz;   // 0.5 T
  static constexpr auto BeamX               =  0._mm;               // Center
  static constexpr auto BeamY               =  0._mm;               // Center
  static constexpr auto ImpactMax           =  30._mm;              // FIXME

  //===--------------------------------------------------------------------===//
  // Track parameter estimation
  //===--------------------------------------------------------------------===//

  static constexpr auto Bmin                =  MagneticField::Bz;   // 0.5 T
  static constexpr auto EstDeltaRmin        =  SeedDeltaRmin;       // OK
  static constexpr auto EstDeltaRmax        =  SeedDeltaRmax;       // OK
  static constexpr auto SigmaLoc0           =  0.5_mm;              // OK
  static constexpr auto SigmaLoc1           =  0.5_mm;              // OK
  static constexpr auto SigmaPhi            =  0.5_degree;          // OK
  static constexpr auto SigmaTheta          =  0.5_degree;          // OK
  static constexpr auto SigmaQOverP         =  0.1 / 1._GeV;        // OK
  static constexpr auto SigmaT0             =  2000._s;             // FIXME
  static constexpr auto InitVarInflatLoc0   =  1.;                  // No inflation
  static constexpr auto InitVarInflatLoc1   =  1.;                  // No inflation
  static constexpr auto InitVarInflatPhi    =  1.;                  // No inflation
  static constexpr auto InitVarInflatTheta  =  1.;                  // No inflation
  static constexpr auto InitVarInflatQOverP =  1.;                  // No inflation
  static constexpr auto InitVarInflatT0     =  1.;                  // No inflation

  //===--------------------------------------------------------------------===//
  // Track finding
  //===--------------------------------------------------------------------===//

  static constexpr auto ResolvePassive      = false;                // FIXME
  static constexpr auto ResolveMaterial     = true;                 // FIXME
  static constexpr auto ResolveSensitive    = true;                 // FIXME
  static constexpr auto PropagationMaxSteps = 1000u;                // FIXME
  static constexpr auto MultipleScattering  = true;                 // FIXME
  static constexpr auto EnergyLoss          = true;                 // FIXME
  static constexpr auto Smoothing           = true;                 // FIXME
  /// Maximum local Chi2 contribution.
  static constexpr auto Chi2max             = 30.0;                 // FIXME
  /// Maximum number of associated measurements on a single surface.
  static constexpr auto NmaxPerSurface      = 5u;                   // FIXME
  static constexpr auto ComputeSharedHits   = false;
  static constexpr auto TrackMinLength      = 4u;
  static constexpr auto NewHitsInRow        = 3u;
  static constexpr auto NewHitsRatio        = 0.25;
  static constexpr auto PostProcess         = false;

  //===--------------------------------------------------------------------===//
  // Performance writer Acts implementation
  //===--------------------------------------------------------------------===//

  static constexpr auto ActsPerfFilePath    = "performance_ckf.root";
  /// Min reco-truth matching probability.
  static constexpr auto ActsTruthMatchProbMin = 0.5;
  static constexpr auto ActsMeasurementsMin = NHitsMin;
  static constexpr auto PerfPtMin           = PtMin;

  // Parameters for EffPlotTool, FakeRatePlotTool,
  // DuplicationPlotTool, TrackSummaryPlotTool.

  static constexpr auto EtaName             = "#eta";
  static constexpr auto EtaNBins            =  40;
  static constexpr auto PerfPlotToolEtaMin  = -EtaMax;
  static constexpr auto PerfPlotToolEtaMax  =  EtaMax;

  static constexpr auto PhiName             = "#phi";
  static constexpr auto PhiNBins            =  100;
  static constexpr auto PerfPlotToolPhiMin  =  SelectorPhiMin;
  static constexpr auto PerfPlotToolPhiMax  =  SelectorPhiMax;

  static constexpr auto PtName              = "pT [GeV/c]";
  static constexpr auto PtNBins             =  40;
  static constexpr auto PerfPlotToolPtMin   =  0._GeV;
  static constexpr auto PerfPlotToolPtMax   =  2.5_GeV;

  static constexpr auto NumName             = "N";
  static constexpr auto NumNBins            =  30;
  static constexpr auto NumMin              = -0.5;
  static constexpr auto NumMax              =  29.5;

  //===--------------------------------------------------------------------===//
  // Performance writer own implementation
  //===--------------------------------------------------------------------===//

  static constexpr auto OPerfFilePath       = "eff.root";
  ///Whether calculate efficiency only for certain trackIds.
  static constexpr auto OnlyCertainTracks   = false;
  /// Path to files with saved trackIds.
  /// Used for calculation efficiency for tracks with certain trackIds only.
  static constexpr auto PathWithTrackIds    = "";
  /// Min reco-truth matching probability.
  static constexpr auto OTruthMatchProbMin  = ActsTruthMatchProbMin;
  static constexpr auto OMeasurementsMin    = ActsMeasurementsMin;

  //===--------------------------------------------------------------------===//
  // System parameters
  //===--------------------------------------------------------------------===//

  /// Whether to dump data for tracks post-processing.
  static constexpr auto DumpData            = false;

  //===--------------------------------------------------------------------===//
  // Constructor
  //===--------------------------------------------------------------------===//

  Config(const BaseTpcSectorGeo &secGeo,
         const std::string &jsonFile,
         Acts::Logging::Level level = Acts::Logging::DEBUG):
      Config(secGeo, "" /* No import */, jsonFile, level) {}

  Config(const BaseTpcSectorGeo &secGeo,
         const std::string &rootFile,
         const std::string &jsonFile,
         Acts::Logging::Level level = Acts::Logging::DEBUG);

  ActsExamples::CKFPerformanceWriter::Config perfWriterCfg(
    std::string outPath) const;

  std::shared_ptr<Detector> detector;

  ParticleSelector::Config particleSelector;
  Digitization::Config digitization;
  SpacePointMaking::Config spacePointMaking;
  TrackSeeding::Config trackSeeding;
  TrackEstimation::Config trackEstimation;
  TrackFinding::Config trackFinding;
};

} // namespace Mpd::Tpc

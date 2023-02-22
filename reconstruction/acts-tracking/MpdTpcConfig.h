// This file is a part of the NICA project.
//
// Copyright (C) 2022 JINR

#pragma once

#include "MpdTpcDetector.h"
#include "MpdTpcDigitization.h"
#include "MpdTpcMagneticField.h"
#include "MpdTpcSpacePointMaking.h"
#include "MpdTpcTrackEstimation.h"
#include "MpdTpcTrackFinding.h"
#include "MpdTpcTrackSeeding.h"

#include "ActsExamples/Io/Performance/CKFPerformanceWriter.hpp"
#include "ActsExamples/TruthTracking/TruthSeedSelector.hpp"

#include <Acts/Definitions/Units.hpp>
#include <Acts/Utilities/Logger.hpp>

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
  static constexpr auto SelectedID          = "selectedparticles";
  static constexpr auto MeasParticlesMapID  = "measurementparticles";

  //===--------------------------------------------------------------------===//
  // Track seeding
  //===--------------------------------------------------------------------===//

  static constexpr auto Rmin                =  Detector::Rmin;      // ~ 0.4 m
  static constexpr auto Rmax                =  Detector::Rmax;      // ~ 1.4 m
  static constexpr auto Zmin                =  Detector::Zmin;      // ~-1.7 m
  static constexpr auto Zmax                =  Detector::Zmax;      // ~ 1.7 m
  static constexpr auto CollisionZmin       = -30._cm;              // Close to 0
  static constexpr auto CollisionZmax       =  30._cm;              // Close to 0
  static constexpr auto CotThetaMax         =  2.;                  // max(dZ/dR)=1.7 ~ 1.3 eta (eta < 1.2)
  static constexpr auto SeedBinSizeR        =  20._mm;              // 10._mm for MC (pads are ~12-18 mm)
  static constexpr auto SeedDeltaRmin       =  10._mm;              // 02._mm for MC
  static constexpr auto SeedDeltaRmax       =  60._mm;              // 20._mm for MC
  static constexpr auto SeedDeltaZmax       =  10._cm;              // FIXME 10._cm for MC
  static constexpr auto MaxSeedsPerSpM      =  3;                   // FIXME
  static constexpr auto SigmaScattering     =  5;                   // FIXME
  static constexpr auto MaxPtScattering     =  5._GeV;              // Max Pt for scattering
  static constexpr auto RadLengthPerSeed    =  0.05;                // OK
  static constexpr auto MinPt               =  0.02_GeV;            // 0.02 < Pt < 10 GeV
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

  //===--------------------------------------------------------------------===//
  // Particle selector
  //===--------------------------------------------------------------------===//

  /// Minimum distance from the origin in the transverse plane.
  static constexpr auto RhoMin              = 0.;
  /// Maximum distance from the origin in the transverse plane.
  static constexpr auto RhoMax              = std::numeric_limits<double>::max();
  /// Minimum absolute distance from the origin along z.
  static constexpr auto ZminSelector        = Zmin;
  /// Maximum absolute distance from the origin along z.
  static constexpr auto ZmaxSelector        = Zmax;
  // Truth particle kinematic cuts.
  static constexpr auto PhiMin              = std::numeric_limits<double>::lowest();
  static constexpr auto PhiMax              = std::numeric_limits<double>::max();
  static constexpr auto EtaMin              = -1.2;
  static constexpr auto EtaMax              =  1.2;
  static constexpr auto AbsEtaMin           = std::numeric_limits<double>::lowest();
  static constexpr auto AbsEtaMax           = std::numeric_limits<double>::max();
  static constexpr auto PtMinSelector       = 0.0;
  static constexpr auto PtMaxSelector       = std::numeric_limits<double>::max();
  /// Keep neutral particles.
  static constexpr auto KeepNeutral         = false;
  /// Requirement on number of recorded hits.
  static constexpr auto NHitsMin            = 9;
  static constexpr auto NHitsMax            = std::numeric_limits<size_t>::max();

  //===--------------------------------------------------------------------===//
  // Performance writer
  //===--------------------------------------------------------------------===//

  // File name will be PerfFilePrefixPath + "_event_<n>.root"
  static constexpr auto PerfFilePrefixPath  = "performance_ckf";
  /// Min reco-truth matching probability.
  static constexpr auto TruthMatchProbMin   = 0.5;
  /// Min number of measurements.
  static constexpr auto MeasurementsMin     = 9u;
  /// Min transverse momentum.
  static constexpr auto PtMinPerf           = MinPt;

  // Parameters for EffPlotTool, FakeRatePlotTool,
  // DuplicationPlotTool, TrackSummaryPlotTool.

  // Pseudorapidity.
  static constexpr auto PlotToolEtaName     = "#eta";
  static constexpr auto PlotToolEtaNBins    = 40;
  static constexpr auto PlotToolEtaMin      = -4;
  static constexpr auto PlotToolEtaMax      = 4;

  static constexpr auto PlotToolPhiName     = "#phi";
  static constexpr auto PlotToolPhiNBins    = 100;
  static constexpr auto PlotToolPhiMin      = -3.15;
  static constexpr auto PlotToolPhiMax      = 3.15;

  static constexpr auto PlotToolPtName      = "pT [GeV/c]";
  static constexpr auto PlotToolPtNBins     = 40;
  static constexpr auto PlotToolPtMin       = 0._GeV;
  static constexpr auto PlotToolPtMax       = 2.5_GeV;

  static constexpr auto PlotToolNumName     = "N";
  static constexpr auto PlotToolNumNBins    = 30;
  static constexpr auto PlotToolNumMin      = -0.5;
  static constexpr auto PlotToolNumMax      = 29.5;

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

  std::shared_ptr<Detector> detector;

  Digitization::Config digitization;
  ActsExamples::TruthSeedSelector::Config truthSeedSelector;
  SpacePointMaking::Config spacePointMaking;
  TrackSeeding::Config trackSeeding;
  TrackEstimation::Config trackEstimation;
  TrackFinding::Config trackFinding;
  ActsExamples::CKFPerformanceWriter::Config perfWriting;
};

} // namespace Mpd::Tpc

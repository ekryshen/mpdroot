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

  //===--------------------------------------------------------------------===//
  // Track seeding
  //===--------------------------------------------------------------------===//

  static constexpr auto Rmin                =  Detector::Rmin;      // ~ 0.4 m
  static constexpr auto Rmax                =  Detector::Rmax;      // ~ 1.4 m 
  static constexpr auto Zmin                =  Detector::Zmin;      // ~-1.7 m
  static constexpr auto Zmax                =  Detector::Zmax;      // ~ 1.7 m
  static constexpr auto CollisionZmin       = -30._cm;              // Close to 0
  static constexpr auto CollisionZmax       =  30._cm;              // Close to 0
  static constexpr auto CotThetaMax         =  1.7;                 // max(dZ/dR)=1.7 ~ 1.3 eta (eta < 1.2)
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
  static constexpr auto ImpactMax           =  10._mm;              // FIXME

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
  static constexpr auto Chi2max             = 50.0;                 // FIXME
  /// Maximum number of associated measurements on a single surface.
  static constexpr auto NmaxPerSurface      = 5u;                   // FIXME
  static constexpr auto ComputeSharedHits   = false;
  static constexpr auto TrackMinLength      = 4u;
  static constexpr auto NewHitsInRow        = 3u;
  static constexpr auto NewHitsRatio        = 0.25;

  Config(const std::string &jsonFile,
         Acts::Logging::Level level = Acts::Logging::DEBUG):
      Config("" /* No import */, jsonFile, level) {}

  Config(const std::string &rootFile,
         const std::string &jsonFile,
         Acts::Logging::Level level = Acts::Logging::DEBUG);

  std::shared_ptr<Detector> detector;

  Digitization::Config digitization;
  SpacePointMaking::Config spacePointMaking;
  TrackSeeding::Config trackSeeding;
  TrackEstimation::Config trackEstimation;
  TrackFinding::Config trackFinding;
};

} // namespace Mpd::Tpc

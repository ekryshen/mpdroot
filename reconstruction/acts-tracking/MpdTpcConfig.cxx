// This file is part of the NICA project.
//
// Copyright (C) 2022 JINR

#include "MpdTpcConfig.h"

#include <Acts/Surfaces/PerigeeSurface.hpp>

namespace Mpd::Tpc {

Config::Config(const BaseTpcSectorGeo &secGeo,
               const std::string &rootFile,
               const std::string &jsonFile,
               Acts::Logging::Level level):
    detector(std::make_shared<Detector>(secGeo, rootFile, jsonFile, level)) {
  auto trackingGeometry = detector->getGeometry();
  assert(trackingGeometry && "Missing tracking geometry");

  auto magneticField = MagneticField::build();
  assert(magneticField && "Missing magnetic field");

  auto referenceSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(
      Acts::Vector3{0., 0., 0.});
  assert(referenceSurface && "Missing reference surface");

  // Digitization.
  digitization.inputSimHits = SimHitsID;
  digitization.outputSourceLinks = SourceLinksID;
  digitization.outputMeasurements = MeasurementsID;
  digitization.sigmaLoc0 = SigmaLoc0;
  digitization.sigmaLoc1 = SigmaLoc1;
  digitization.detector = detector;

  // Space point making.
  spacePointMaking.inputSourceLinks = SourceLinksID;
  spacePointMaking.inputMeasurements = MeasurementsID;
  spacePointMaking.outputSpacePoints = SpacePointsID;
  spacePointMaking.trackingGeometry = trackingGeometry;
  spacePointMaking.geometrySelection = {Acts::GeometryIdentifier().setVolume(0)};

  // Track seeding.
  trackSeeding.inputSpacePoints = SpacePointsID;
  trackSeeding.outputSeeds = SeedsID;
  trackSeeding.outputProtoTracks = ProtoTracksID;
  trackSeeding.seedFilterConfig.deltaRMin = SeedDeltaRmin;
  trackSeeding.seedFilterConfig.maxSeedsPerSpM = MaxSeedsPerSpM;
  trackSeeding.seedFinderConfig.rMin = Rmin;
  trackSeeding.seedFinderConfig.rMax = Rmax;
  trackSeeding.seedFinderConfig.binSizeR = SeedBinSizeR;
  trackSeeding.seedFinderConfig.deltaRMin = SeedDeltaRmin;
  trackSeeding.seedFinderConfig.deltaRMinTopSP = SeedDeltaRmin;
  trackSeeding.seedFinderConfig.deltaRMinBottomSP = SeedDeltaRmin;
  trackSeeding.seedFinderConfig.deltaRMax = SeedDeltaRmax;
  trackSeeding.seedFinderConfig.deltaRMaxTopSP = SeedDeltaRmax;
  trackSeeding.seedFinderConfig.deltaRMaxBottomSP = SeedDeltaRmax;
  trackSeeding.seedFinderConfig.collisionRegionMin = CollisionZmin;
  trackSeeding.seedFinderConfig.collisionRegionMax = CollisionZmax;
  trackSeeding.seedFinderConfig.zMin = Zmin;
  trackSeeding.seedFinderConfig.zMax = Zmax;
  trackSeeding.seedFinderConfig.deltaZMax = SeedDeltaZmax;
  trackSeeding.seedFinderConfig.maxSeedsPerSpM = MaxSeedsPerSpM;
  trackSeeding.seedFinderConfig.cotThetaMax = CotThetaMax;
  trackSeeding.seedFinderConfig.maxPtScattering = MaxPtScattering;
  trackSeeding.seedFinderConfig.sigmaScattering = SigmaScattering;
  trackSeeding.seedFinderConfig.radLengthPerSeed = RadLengthPerSeed;
  trackSeeding.seedFinderConfig.minPt = MinPt;
  trackSeeding.seedFinderConfig.bFieldInZ = Bz;
  trackSeeding.seedFinderConfig.beamPos = {BeamX, BeamY};
  trackSeeding.seedFinderConfig.impactMax = ImpactMax;
  trackSeeding.gridConfig.rMax = Rmax;
  trackSeeding.gridConfig.deltaRMax = SeedDeltaRmax;
  trackSeeding.gridConfig.zMin = Zmin;
  trackSeeding.gridConfig.zMax = Zmax;
  trackSeeding.gridConfig.cotThetaMax = CotThetaMax;
  trackSeeding.gridConfig.minPt = MinPt;
  trackSeeding.gridConfig.bFieldInZ = Bz;

  // Track parameter estimation.
  trackEstimation.inputSeeds = SeedsID;
  trackEstimation.inputProtoTracks = ProtoTracksID;
  trackEstimation.inputSpacePoints = SpacePointsID;
  trackEstimation.inputMeasurements = MeasurementsID;
  trackEstimation.outputTrackParameters = EstTrackParamsID;
  trackEstimation.outputProtoTracks = EstProtoTracksID;
  trackEstimation.trackingGeometry = trackingGeometry;
  trackEstimation.magneticField = magneticField;
  trackEstimation.bFieldMin = Bmin;
  trackEstimation.deltaRMin = EstDeltaRmin;
  trackEstimation.deltaRMax = EstDeltaRmax;
  trackEstimation.sigmaLoc0 = SigmaLoc0;
  trackEstimation.sigmaLoc1 = SigmaLoc1;
  trackEstimation.sigmaPhi = SigmaPhi;
  trackEstimation.sigmaTheta = SigmaTheta;
  trackEstimation.sigmaQOverP = SigmaQOverP;
  trackEstimation.sigmaT0 = SigmaT0;
  trackEstimation.initialVarInflation = {
      InitVarInflatLoc0,
      InitVarInflatLoc1,
      InitVarInflatPhi,
      InitVarInflatTheta,
      InitVarInflatQOverP,
      InitVarInflatT0
  };

  // Track finding.
  trackFinding.inputMeasurements = MeasurementsID;
  trackFinding.inputSourceLinks = SourceLinksID;
  trackFinding.inputInitialTrackParameters = EstTrackParamsID;
  trackFinding.outputTrajectories = TrajectoriesID;
  trackFinding.outputTrackCandidates = TrackCandidatesID;
  trackFinding.trackingGeometry = trackingGeometry;
  trackFinding.magneticField = magneticField;
  trackFinding.referenceSurface = referenceSurface;
  trackFinding.navigatorConfig.trackingGeometry = trackingGeometry;
  trackFinding.navigatorConfig.resolvePassive = ResolvePassive;
  trackFinding.navigatorConfig.resolveMaterial = ResolveMaterial;
  trackFinding.navigatorConfig.resolveSensitive = ResolveSensitive;
  trackFinding.propagatorOptions.maxSteps = PropagationMaxSteps;
  trackFinding.multipleScattering = MultipleScattering;
  trackFinding.energyLoss = EnergyLoss;
  trackFinding.smoothing = Smoothing;
  trackFinding.measurementSelectorConfig = {{
      Acts::GeometryIdentifier(), {
          {/* Single |eta| bin */},
          {Chi2max},
          {NmaxPerSurface}
      }
  }};
  trackFinding.computeSharedHits = ComputeSharedHits;
  trackFinding.trackMinLength = TrackMinLength;
  trackFinding.newHitsInRow = NewHitsInRow;
  trackFinding.newHitsRatio = NewHitsRatio;

  // Performance writing.
  perfWriting.inputTrajectories = TrajectoriesID;
  perfWriting.inputParticles = ParticlesID;
  perfWriting.inputMeasurementParticlesMap = MeasParticlesMapID;
  perfWriting.filePath = PerfFilePath;
  perfWriting.truthMatchProbMin = TruthMatchProbMin;
  perfWriting.nMeasurementsMin = MeasurementsMin;
  perfWriting.ptMin = PtMin;
}

} // namespace Mpd::Tpc

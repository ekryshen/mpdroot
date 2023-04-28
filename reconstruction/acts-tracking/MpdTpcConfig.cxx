// This file is part of the NICA project.
//
// Copyright (C) 2022 JINR

#include "MpdTpcConfig.h"

#include "ActsExamples/Utilities/Helpers.hpp"

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

  // Selecting particles.
  particleSelector.truthSeedSelectorConfig.inputParticles = ParticlesID;
  particleSelector.truthSeedSelectorConfig.inputMeasurementParticlesMap =
      HitParticlesMapID;
  particleSelector.truthSeedSelectorConfig.outputParticles =
      SelectedParticlesID;
  particleSelector.outputHitsParticlesMap = SelectedHitParticlesMapID;
  particleSelector.inputSimHits = SimHitsID;
  particleSelector.outputSimHits = SelectedSimHitsID;
  particleSelector.primaryParticlesOnly = PrimaryParticlesOnly;
  particleSelector.truthSeedSelectorConfig.rhoMin = Rmin;
  particleSelector.truthSeedSelectorConfig.rhoMax = Rmax;
  particleSelector.truthSeedSelectorConfig.zMin = Zmin;
  particleSelector.truthSeedSelectorConfig.zMax = Zmax;
  particleSelector.truthSeedSelectorConfig.phiMin = PhiMin;
  particleSelector.truthSeedSelectorConfig.phiMax = PhiMax;
  particleSelector.truthSeedSelectorConfig.etaMin = EtaMin;
  particleSelector.truthSeedSelectorConfig.etaMax = EtaMax;
  particleSelector.truthSeedSelectorConfig.absEtaMin = AbsEtaMin;
  particleSelector.truthSeedSelectorConfig.absEtaMax = AbsEtaMax;
  particleSelector.truthSeedSelectorConfig.ptMin = PtMin;
  particleSelector.truthSeedSelectorConfig.ptMax = PtMax;
  particleSelector.truthSeedSelectorConfig.keepNeutral = KeepNeutral;
  particleSelector.truthSeedSelectorConfig.nHitsMin = NHitsMin;
  particleSelector.truthSeedSelectorConfig.nHitsMax = NHitsMax;

  // Digitization.
  digitization.inputSimHits = SelectedSimHitsID;
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
  trackSeeding.seedFinderConfig.minPt = PtMin;
  trackSeeding.seedFinderConfig.bFieldInZ = Bz;
  trackSeeding.seedFinderConfig.beamPos = {BeamX, BeamY};
  trackSeeding.seedFinderConfig.impactMax = ImpactMax;
  trackSeeding.gridConfig.rMax = Rmax;
  trackSeeding.gridConfig.deltaRMax = SeedDeltaRmax;
  trackSeeding.gridConfig.zMin = Zmin;
  trackSeeding.gridConfig.zMax = Zmax;
  trackSeeding.gridConfig.cotThetaMax = CotThetaMax;
  trackSeeding.gridConfig.minPt = PtMin;
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
  trackFinding.spacePointsID = SpacePointsID;
  trackFinding.dumpData = DumpData;
}

ActsExamples::CKFPerformanceWriter::Config Config::perfWriterCfg(
    std::string outPath) const {
  ActsExamples::CKFPerformanceWriter::Config result;
  result.inputTrajectories = TrajectoriesID;
  result.inputParticles = SelectedParticlesID;
  result.inputMeasurementParticlesMap = SelectedHitParticlesMapID;
  result.truthMatchProbMin = TruthMatchProbMin;
  result.nMeasurementsMin = MeasurementsMin;
  result.ptMin = PtMinPerf;

  ActsExamples::EffPlotTool::Config effConfig;
  effConfig.varBinning["Eta"] = ActsExamples::PlotHelpers::Binning(
      PlotToolEtaName,
      PlotToolEtaNBins,
      PlotToolEtaMin,
      PlotToolEtaMax);
  effConfig.varBinning["Phi"] = ActsExamples::PlotHelpers::Binning(
      PlotToolPhiName,
      PlotToolPhiNBins,
      PlotToolPhiMin,
      PlotToolPhiMax);
  effConfig.varBinning["Pt"]  = ActsExamples::PlotHelpers::Binning(
      PlotToolPtName,
      PlotToolPtNBins,
      PlotToolPtMin,
      PlotToolPtMax);
  result.effPlotToolConfig = effConfig;

  ActsExamples::FakeRatePlotTool::Config fakeConfig;
  fakeConfig.varBinning["Eta"] = ActsExamples::PlotHelpers::Binning(
      PlotToolEtaName,
      PlotToolEtaNBins,
      PlotToolEtaMin,
      PlotToolEtaMax);
  fakeConfig.varBinning["Phi"] = ActsExamples::PlotHelpers::Binning(
      PlotToolPhiName,
      PlotToolPhiNBins,
      PlotToolPhiMin,
      PlotToolPhiMax);
  fakeConfig.varBinning["Pt"]  = ActsExamples::PlotHelpers::Binning(
      PlotToolPtName,
      PlotToolPtNBins,
      PlotToolPtMin,
      PlotToolPtMax);
  fakeConfig.varBinning["Num"] = ActsExamples::PlotHelpers::Binning(
      PlotToolNumName,
      PlotToolNumNBins,
      PlotToolNumMin,
      PlotToolNumMax);
  result.fakeRatePlotToolConfig = fakeConfig;

  ActsExamples::DuplicationPlotTool::Config duplicationConfig;
  duplicationConfig.varBinning["Eta"] = ActsExamples::PlotHelpers::Binning(
      PlotToolEtaName,
      PlotToolEtaNBins,
      PlotToolEtaMin,
      PlotToolEtaMax);
  duplicationConfig.varBinning["Phi"] = ActsExamples::PlotHelpers::Binning(
      PlotToolPhiName,
      PlotToolPhiNBins,
      PlotToolPhiMin,
      PlotToolPhiMax);
  duplicationConfig.varBinning["Pt"] = ActsExamples::PlotHelpers::Binning(
      PlotToolPtName,
      PlotToolPtNBins,
      PlotToolPtMin,
      PlotToolPtMax);
  duplicationConfig.varBinning["Num"] = ActsExamples::PlotHelpers::Binning(
      PlotToolNumName,
      PlotToolNumNBins,
      PlotToolNumMin,
      PlotToolNumMax);
  result.duplicationPlotToolConfig = duplicationConfig;

  ActsExamples::TrackSummaryPlotTool::Config trackSummaryConfig;
  trackSummaryConfig.varBinning["Eta"] = ActsExamples::PlotHelpers::Binning(
      PlotToolEtaName,
      PlotToolEtaNBins,
      PlotToolEtaMin,
      PlotToolEtaMax);
  trackSummaryConfig.varBinning["Phi"] = ActsExamples::PlotHelpers::Binning(
      PlotToolPhiName,
      PlotToolPhiNBins,
      PlotToolPhiMin,
      PlotToolPhiMax);
  trackSummaryConfig.varBinning["Pt"] = ActsExamples::PlotHelpers::Binning(
      PlotToolPtName,
      PlotToolPtNBins,
      PlotToolPtMin,
      PlotToolPtMax);
  trackSummaryConfig.varBinning["Num"] = ActsExamples::PlotHelpers::Binning(
      PlotToolNumName,
      PlotToolNumNBins,
      PlotToolNumMin,
      PlotToolNumMax);
  result.trackSummaryPlotToolConfig = trackSummaryConfig;
  result.filePath = outPath + "/" + Config::PerfFilePath;

  return result;
}

} // namespace Mpd::Tpc

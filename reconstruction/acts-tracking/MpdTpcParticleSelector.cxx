// This file is a part of the NICA project.
//
// Copyright (C) 2023 JINR

#include "MpdTpcInputHit.h"
#include "MpdTpcParticleSelector.h"

#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Utilities/Range.hpp"

#include <Acts/Utilities/Helpers.hpp>
#include <ActsFatras/EventData/Barcode.hpp>
#include <ActsFatras/EventData/Particle.hpp>

#include <map>
#include <set>

namespace Mpd::Tpc {

  ParticleSelector::ParticleSelector(Config config,
                                     Acts::Logging::Level level)
      : Algorithm("ParticleSelector", level),
        m_config(std::move(config)) {
    assert(!m_config.truthSeedSelectorConfig.inputParticles.empty()
      && "Missing input truth particles collection");
    assert(!m_config.truthSeedSelectorConfig.inputMeasurementParticlesMap.empty()
      && "Missing input hit-particles map collection");
    assert(!m_config.truthSeedSelectorConfig.outputParticles.empty()
      && "Missing output truth particles collection");
  }

  // Get barcodes of particles.
  std::set<ActsFatras::Barcode> getParticleIds(
      const ActsExamples::SimParticleContainer &particles) {
    std::set<ActsFatras::Barcode> result;
    for (auto particle : particles) {
      ActsFatras::Barcode particleId = particle.particleId();
      result.insert(particleId);
    }
    return result;
  }

  // Make a selection of hits and hitParticlesMap
  // based on already selected particles.
  std::pair<InputHitContainer, ActsExamples::IndexMultimap<ActsFatras::Barcode>>
      selectHitsAndHitParticlesMap(
          const InputHitContainer &hits,
          const ActsExamples::IndexMultimap<ActsFatras::Barcode>
              &hitParticlesMap,
          const std::set<ActsFatras::Barcode> &particlesIds) {

    ActsExamples::IndexMultimap<ActsFatras::Barcode> selectedHitParticlesMap;
    InputHitContainer selectedHits;

    Int_t iselectedHit = 0;
    Int_t ihit = -1;
    for (auto hit : hits) {
      ihit++;

      Bool_t needAddHit = false;

      // Loop over multimap particlesIds[ihit]
      auto pairs = hitParticlesMap.equal_range(ihit);
      for (auto hitParticle = pairs.first; hitParticle != pairs.second;
          hitParticle++) {
        ActsFatras::Barcode particleId = hitParticle->second;
        if (particlesIds.count(particleId) == 0) {
          continue;
        }
        selectedHitParticlesMap.emplace(iselectedHit++, particleId);
        needAddHit = true;
      }
      if (needAddHit) {
        selectedHits.push_back(hit);
      }
    }
    return std::make_pair(selectedHits, selectedHitParticlesMap);
  }

  ActsExamples::ProcessCode ParticleSelector::execute(
      const ActsExamples::AlgorithmContext &context) const {
    using HitParticlesMap = ActsExamples::IndexMultimap<ActsFatras::Barcode>;

    auto &seedConfig = m_config.truthSeedSelectorConfig;

    // prepare input collections
    const auto& inputParticles =
        context.eventStore.get<ActsExamples::SimParticleContainer>(
            seedConfig.inputParticles);
    const auto& hitParticlesMap =
        context.eventStore.get<HitParticlesMap>(
            seedConfig.inputMeasurementParticlesMap);
    // compute particle_id -> {hit_id...} map from the
    // hit_id -> {particle_id...} map on the fly.
    const auto& particleHitsMap =
        ActsExamples::invertIndexMultimap(hitParticlesMap);

    // prepare output collection
    ActsExamples::SimParticleContainer selectedParticles;
    selectedParticles.reserve(inputParticles.size());

    auto within = [](Double_t x, Double_t min, Double_t max) {
      return (min <= x) and (x < max);
    };
    auto isValidparticle = [&](const auto& p) {
      const auto eta = Acts::VectorHelpers::eta(p.unitDirection());
      const auto phi = Acts::VectorHelpers::phi(p.unitDirection());
      const auto rho = Acts::VectorHelpers::perp(p.position());
      // find the corresponding hits for this particle
      const auto& hits = ActsExamples::makeRange(
          particleHitsMap.equal_range(p.particleId()));
      // number of recorded hits
      size_t nHits = hits.size();

      auto particleId = p.particleId();

      return
          within(rho, 0., seedConfig.rhoMax) and
          within(p.position().z(), seedConfig.zMin, seedConfig.zMax) and
          within(std::abs(eta), seedConfig.absEtaMin, seedConfig.absEtaMax) and
          within(eta, seedConfig.etaMin, seedConfig.etaMax) and
          within(phi, seedConfig.phiMin, seedConfig.phiMax) and
          within(p.transverseMomentum(), seedConfig.ptMin, seedConfig.ptMax) and
          within(nHits, seedConfig.nHitsMin, seedConfig.nHitsMax) and
          (seedConfig.keepNeutral or (p.charge() != 0)) and
          (not m_config.primaryParticlesOnly or (particleId.generation() == 0));
    };

    InputHitContainer selectedHits;
    ActsExamples::IndexMultimap<ActsFatras::Barcode> selectedMap;

    if (m_config.selectorEnabled) {
      for (const auto& particle : inputParticles) {
        if (isValidparticle(particle)) {
          selectedParticles.insert(particle);
        }
      }
    } else {
      // Create prototracks for all input particles
      selectedParticles = std::move(inputParticles);
    }
    context.eventStore.add(seedConfig.outputParticles,
        std::move(selectedParticles));

    return ActsExamples::ProcessCode::SUCCESS;
  }

} // namespace Mpd::Tpc

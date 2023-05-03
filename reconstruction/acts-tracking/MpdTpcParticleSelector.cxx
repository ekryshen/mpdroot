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
    assert(!m_config.inputSimHits.empty()
      && "Missing input hits collection");
    assert(!m_config.outputSimHits.empty()
      && "Missing output hits collection");
    assert(!m_config.outputHitsParticlesMap.empty()
      && "Missing output hit-particles map collection");
  }

  std::pair<InputHitContainer, ActsExamples::IndexMultimap<ActsFatras::Barcode>>
      selectHitsAndHitParticlesMap(
          const InputHitContainer &hits,
          const ActsExamples::IndexMultimap<ActsFatras::Barcode>
              &hitParticlesMap,
          const ActsExamples::SimParticleContainer &particles) {
    InputHitContainer selectedHits;
    ActsExamples::IndexMultimap<ActsFatras::Barcode> selectedHitParticlesMap;

    size_t nHits = hits.size();
    selectedHits.reserve(nHits);

    size_t newIndex = 0;
    std::set<size_t> hitAdded;
    std::map<size_t, size_t> oldNewIndexMap;

    for (const auto &[hitIndex, particleId] : hitParticlesMap) {
      if (particles.find(particleId) == particles.end()) {
        continue;
      }
      if (hitAdded.count(hitIndex) == 0) {
        auto hit = hits.at(hitIndex);
        selectedHits.push_back(hit);
        hitAdded.insert(hitIndex);
        oldNewIndexMap[hitIndex] = newIndex++;
      }
      selectedHitParticlesMap.emplace(oldNewIndexMap[hitIndex], particleId);
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

    auto within = [](double x, double min, double max) {
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

    // create prototracks for all input particles
    for (const auto& particle : inputParticles) {
      if (isValidparticle(particle)) {
        selectedParticles.insert(particle);
      }
    }

    const auto& hits = context.eventStore.get<InputHitContainer>(
        m_config.inputSimHits);

    auto [selectedHits, selectedHitParticlesMap] = selectHitsAndHitParticlesMap(
        hits,
        hitParticlesMap,
        selectedParticles);

    context.eventStore.add(seedConfig.outputParticles,
        std::move(selectedParticles));
    context.eventStore.add(m_config.outputSimHits,std::move(selectedHits));
    context.eventStore.add(m_config.outputHitsParticlesMap,
        std::move(selectedHitParticlesMap));

    return ActsExamples::ProcessCode::SUCCESS;
  }

} // namespace Mpd::Tpc

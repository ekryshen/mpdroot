// This file is a part of the NICA project.
//
// Copyright (C) 2022 JINR

#pragma once

#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <Acts/Utilities/Logger.hpp>

#include <memory>
#include <string>

namespace Mpd::Tpc {

/// @brief Base class for algorithms (subtasks) used in the Acts-based tracker.
class Algorithm : public ActsExamples::IAlgorithm {
public:
  Algorithm(std::string name,
            Acts::Logging::Level level = Acts::Logging::INFO):
      m_name(std::move(name)),
      m_logger(Acts::getDefaultLogger(m_name, level)) {}

  virtual ~Algorithm() = default;

  std::string name() const override { return m_name; }

  /// Executes the algorithm for one event.
  ActsExamples::ProcessCode execute(
    const ActsExamples::AlgorithmContext &context) const override = 0;

  /// Initialize the algorithm
  ActsExamples::ProcessCode initialize() const override {
    return ActsExamples::ProcessCode::SUCCESS;
  };

  /// Finalize the algorithm
  ActsExamples::ProcessCode finalize() const override {
    return ActsExamples::ProcessCode::SUCCESS;
  };

protected:
  const Acts::Logger &logger() const { return *m_logger; }

private:
  std::string m_name;
  std::unique_ptr<const Acts::Logger> m_logger;
};

} // namespace Mpd::Tpc

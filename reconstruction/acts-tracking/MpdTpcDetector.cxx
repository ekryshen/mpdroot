// This file is a part of the NICA project.
//
// Copyright (C) 2022 JINR

#include "MpdTpcDetector.h"

#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Utilities/Logger.hpp>
#include <Acts/Visualization/EventDataView3D.hpp>
#include <Acts/Visualization/ObjVisualization3D.hpp>

#include <TGeoManager.h>
#include <TGeoNode.h>
#include <TGeoVolume.h>
#include <TMath.h>
#include <TString.h>

#include <cassert>
#include <cmath>
#include <iostream>
#include <sstream>

namespace Mpd::Tpc {

/// Converts Acts length (in mm) to ROOT units.
inline Double_t toRootLength(Acts::ActsScalar actsLength) {
  // ROOT uses both cm (TGeoManager::kRootUnits) and mm (TGeoManager::kG4Units).
  const auto unitScalor = 1. / 1._cm;
  return unitScalor * actsLength; 
}

/// Produces name by concatenating prefix and suffix.
template<typename T>
inline std::string makeName(const TString &prefix, const T &suffix) {
  std::stringstream out;
  out << prefix << "_" << suffix;
  return out.str();
}

/// Creates a sensor volume.
inline TGeoVolume *makeSensorVolume(TGeoManager *geoManager,
                                    TGeoMedium *medium,
                                    Double_t dX,
                                    Double_t dY,
                                    Double_t dZ,
                                    Int_t ivolume) {
  auto sensorName = makeName("tpc_sensor", ivolume);
  auto *sensorVolume = geoManager->MakeBox(
      sensorName.c_str(), medium, dX, dY, dZ);

  sensorVolume->SetLineColor(kBlack);
  sensorVolume->SetFillColor(kYellow);
  sensorVolume->SetVisibility(kTRUE);

  sensorVolume->Print();
  sensorVolume->CheckOverlaps();

  return sensorVolume;
}

/// Adds a sensor at the given position.
inline void addSensor(TGeoVolume *gasVolume,
                      TGeoVolume *sensorVolume,
                      Double_t x,
                      Double_t y,
                      Double_t z,
                      Double_t phi,
                      Int_t &isensor) {
  auto phiDegree = TMath::RadToDeg() * phi;

  auto rotName = makeName("tpc_rotation", isensor);
  auto *rot = new TGeoRotation(rotName.c_str(), phiDegree, 0., 0.);
  auto *loc = new TGeoCombiTrans(x, y, z, rot);
 
  gasVolume->AddNode(sensorVolume, isensor, loc);
  isensor++;
}

inline void addSensorsPhi(TGeoVolume *gasVolume,
                          TGeoVolume *sensorVolume,
                          Double_t r,
                          Double_t z,
                          Int_t &isensor) {
  for (size_t isec = 0; isec < Detector::NumSectors; isec++) {
    for (size_t iphi = 0; iphi < Detector::NumPhi; iphi++) {
      auto phi = (isec * Detector::SecDeltaPhi) +
                 (iphi * Detector::IntDeltaPhi) + Detector::PhiStart;

      auto x = r * std::cos(phi);
      auto y = r * std::sin(phi);

      addSensor(gasVolume, sensorVolume, x, y, z, phi, isensor);
    }
  }
}

/// Adds sensors along the Z axis.
inline void addSensorsZPhi(TGeoVolume *gasVolume,
                           TGeoVolume *sensorVolume,
                           Double_t r,
                           Int_t &isensor) {
  auto nz = Detector::HasTwoNodes ? Detector::NumZ / 2
                                  : Detector::NumZ;
  auto mz = Detector::HasTwoNodes ? Detector::Zmin / 2.
                                  : Detector::Zmin;

  for (size_t iz = 0; iz < nz; iz++) {
    auto z = toRootLength(mz + ((iz + 0.5) * Detector::DeltaZ));
    addSensorsPhi(gasVolume, sensorVolume, r, z, isensor);
  }
}

void addSensors(TGeoManager *geoManager, TGeoVolume *gasVolume) {
  auto *medium = gasVolume->GetMedium();
  assert(medium && "Missing TPC gas medium");

  const auto tanHalfDeltaPhi = std::tan(0.5 * Detector::IntDeltaPhi);

  auto isensor = 0;
  for(size_t ilayer = 0; ilayer < Detector::NumLayers; ilayer++) {
    // Set r*phi so as sensors touch at the inner edges.
    auto minR = toRootLength(Detector::Rmin + (ilayer * Detector::DeltaR));
    auto rPhi = 2. * tanHalfDeltaPhi * minR;

    auto sizeX = toRootLength(Detector::DeltaR);
    auto sizeY = rPhi;
    auto sizeZ = toRootLength(Detector::DeltaZ);

    auto eps = 0.00001;
    auto deltaX = (0.5 - eps) * sizeX; // Rotation error
    auto deltaY = (0.5 - eps) * sizeY; // Rotation error
    auto deltaZ = 0.5 * sizeZ;
 
    auto *sensorVolume = makeSensorVolume(
        geoManager, medium, deltaX, deltaY, deltaZ, ilayer);

    auto r = minR + (0.5 * sizeX);
    addSensorsZPhi(gasVolume, sensorVolume, r, isensor);
  }
}

bool Detector::editGeometry(TGeoManager *geoManager) {
  TGeoVolume *gasVolume = geoManager->GetVolume(Detector::GasVolume);
  if (!gasVolume) {
    ACTS_ERROR("Volume '" << Detector::GasVolume << "' not found");
    return false;
  }

  ACTS_DEBUG("TPC gas volume: ");
  gasVolume->Print();

  addSensors(geoManager, gasVolume);

  if (m_logger->doPrint(Acts::Logging::DEBUG)) {
    geoManager->Export("tpc_acts_tracking.root");
  }

  geoManager->CloseGeometry();
  return true;
}

Detector::Detector(const std::string &rootFile,
                   const std::string &jsonFile,
                   Acts::Logging::Level level):
    m_logger(Acts::getDefaultLogger("Detector", level)) {
  geoEditor = [this](auto *gm){ return editGeometry(gm); };

  m_config.surfaceLogLevel = level;
  m_config.layerLogLevel = level;
  m_config.volumeLogLevel = level;

  m_config.fileName = rootFile;
  m_config.readJson(jsonFile);
}

std::shared_ptr<const Acts::TrackingGeometry> Detector::getGeometry() {
  // Use the cached value if available.
  if (m_trackingGeometry != nullptr) {
    return m_trackingGeometry;
  }

  const auto &geometry = finalize(m_config, nullptr);
  m_trackingGeometry = geometry.first;

  // Debug output.
  if (m_logger->doPrint(Acts::Logging::VERBOSE)) {
    Acts::GeometryContext context;
    Acts::ObjVisualization3D visualizer;

    m_trackingGeometry->visitSurfaces([&](const Acts::Surface *surface) {
      // Visual representation.
      Acts::GeometryView3D::drawSurface(visualizer, *surface, context);

      // Textual representation.
      std::ostream &os = std::cout;
      os << std::endl << "geometryId=" << surface->geometryId() << std::endl;
      surface->toStream(context, os);
    });

    visualizer.write("tpc_acts_tracking");
  }

  return m_trackingGeometry;
}

const Acts::TrackingVolume *findVolume(
    const Acts::TrackingVolume *master, const std::string &name) {
  if (master->volumeName() == name) {
    return master;
  }
  if(master->volumeName().empty() || master->volumeName()[0] != '{') {
    return nullptr;
  }

  for(const auto &child : master->confinedVolumes()->arrayObjects()) {
    if(child->volumeName() == name) {
      return child.get();
    }
    if(auto found = findVolume(child.get(), name)) {
      return found;
    }
  }

  return nullptr;
}

inline Acts::ActsScalar getPhi(Acts::ActsScalar x,
                               Acts::ActsScalar y,
                               Acts::ActsScalar deltaPhi) {
  auto phi = std::acos(x / std::hypot(x, y));
  if (y < 0) { phi = 2*M_PI - phi; }
  phi += deltaPhi;

  const auto eps = 0.001;
  if (phi >= 2*M_PI - eps) { phi -= 2*M_PI; }
  if (phi < 0) { phi = 0; }

  return phi;
}

inline uint64_t toBits(const Acts::Vector3 &position,
                       Acts::ActsScalar deltaPhi) {
  const auto x = position[0];
  const auto y = position[1];
  const auto z = position[2];
  assert(Detector::Zmin <= z && z <= Detector::Zmax);

  const auto r = std::hypot(x, y);
  assert(Detector::Rmin <= r && r <= Detector::Rmax);

  const auto eps = 0.001;
  const auto phi = getPhi(x, y, deltaPhi) + eps;
  assert(0 <= phi && phi < 2*M_PI);

  auto rBits = static_cast<uint64_t>((r - Detector::Rmin) / Detector::DeltaR);
  auto pBits = static_cast<uint64_t>(phi / Detector::IntDeltaPhi);
  auto zBits = static_cast<uint64_t>((z - Detector::Zmin) / Detector::DeltaZ);

  return (rBits << 32) | (pBits << 16) | zBits;
}

void fillSurfaces(
    const Acts::GeometryContext &context,
    const std::shared_ptr<const Acts::TrackingGeometry> &geometry,
    std::unordered_map<uint64_t,
                       std::shared_ptr<const Acts::Surface>> &surfaces) {
  auto topVolume = geometry->highestTrackingVolume();
  assert(topVolume);

  auto tpcVolumeName = std::string(Detector::GasVolume) + "::Barrel";
  auto tpcVolume = findVolume(topVolume, tpcVolumeName);
  assert(tpcVolume);

  auto tpcLayerVector = tpcVolume->confinedLayers()->arrayObjects();
  for(size_t i = 0; i < tpcLayerVector.size(); i++) {
    auto surfaceArray = tpcLayerVector.at(i)->surfaceArray();
    if(surfaceArray == nullptr) {
      continue;
    }

    auto surfaceVector = surfaceArray->surfaces();
    for(size_t j = 0; j < surfaceVector.size(); j++) {
      auto surface = surfaceVector.at(j)->getSharedPtr();
      auto center = surface->center(context);
      auto bits = toBits(center, 0.5 * Detector::IntDeltaPhi);

      assert(surfaces.find(bits) == surfaces.end()
          && "There are surfaces with the same bit keys");
      surfaces.insert({bits, surface});
    }
  }
}

std::shared_ptr<const Acts::Surface> Detector::getSurface(
    const Acts::GeometryContext &gcontext,
    const Acts::Vector3 &position) {
  if (m_surfaces.empty()) {
    fillSurfaces(gcontext, m_trackingGeometry, m_surfaces);
    ACTS_DEBUG("Constructed " << m_surfaces.size() << " surfaces");
  }

  ACTS_VERBOSE("Finding a surface for " << position);
  auto bits = toBits(position, 0.5 * Detector::IntDeltaPhi);
  auto iter = m_surfaces.find(bits);
  assert(iter != m_surfaces.end() && "Surface not found");

  return iter->second;
}

} // namespace Mpd::Tpc

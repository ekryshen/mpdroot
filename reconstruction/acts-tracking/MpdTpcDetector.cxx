// This file is a part of the NICA project.
//
// Copyright (C) 2022 JINR

#include "MpdTpcDetector.h"

#include "BaseTpcSectorGeo.h"

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
//sensorVolume->CheckOverlaps();

  return sensorVolume;
}

//===----------------------------------------------------------------------===//
// Virtual surfaces arranged in cylindrical layers
//===----------------------------------------------------------------------===//

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
    auto phi = (isec * Detector::DeltaPhi) + Detector::PhiStart;

    auto x = r * std::cos(phi);
    auto y = r * std::sin(phi);

    addSensor(gasVolume, sensorVolume, x, y, z, phi, isensor);
  }
}

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

void addCylinderSensors(TGeoManager *geoManager, TGeoVolume *gasVolume) {
  auto *medium = gasVolume->GetMedium();
  assert(medium && "Missing TPC gas medium");

  const auto tanHalfDeltaPhi = std::tan(0.5 * Detector::DeltaPhi);

  auto isensor = 0;
  for(size_t ilayer = 0; ilayer < Detector::NumLayers; ilayer++) {
    // Set r*phi so as sensors touch at the inner edges.
    auto minR = toRootLength(Detector::Rmin + (ilayer * Detector::DeltaR));
    auto rPhi = 2. * tanHalfDeltaPhi * minR;

    auto sizeX = toRootLength(Detector::DeltaR);
    auto sizeY = rPhi;
    auto sizeZ = toRootLength(Detector::DeltaZ);

    auto eps = 0.00001; // FIXME: if 0.02, then no overlapping
    auto deltaX = (0.5 - eps) * sizeX; // Rotation error
    auto deltaY = (0.5 - eps) * sizeY; // Rotation error
    auto deltaZ = 0.5 * sizeZ;

    auto *sensorVolume = makeSensorVolume(
        geoManager, medium, deltaX, deltaY, deltaZ, ilayer);

    auto r = minR + (0.5 * sizeX);
    addSensorsZPhi(gasVolume, sensorVolume, r, isensor);
  }
}

inline Acts::ActsScalar getPhi(Acts::ActsScalar x,
                               Acts::ActsScalar y,
                               Acts::ActsScalar deltaPhi,
                               Double_t eps) {
  auto phi = std::acos(x / std::hypot(x, y));
  if (y < 0) { phi = 2*M_PI - phi; }
  phi += deltaPhi;

  if (phi >= 2*M_PI - eps) { phi -= 2*M_PI; }
  if (phi < 0) { phi = 0; }

  return phi;
}

inline uint64_t toBitsForCylinder(const Acts::Vector3 &position,
                                  Acts::ActsScalar deltaPhi) {
  const auto x = position[0];
  const auto y = position[1];
  const auto z = position[2];
  assert(Detector::Zmin <= z && z <= Detector::Zmax);

  const auto r = std::hypot(x, y);
  assert(Detector::Rmin <= r && r <= Detector::Rmax);

  const auto eps = 0.001;
  const auto epsGetPhi = 0.001;
  const auto phi = getPhi(x, y, deltaPhi, epsGetPhi) + eps;
  assert(0 <= phi && phi < 2*M_PI);

  auto rBits = static_cast<uint64_t>((r - Detector::Rmin) / Detector::DeltaR);
  auto pBits = static_cast<uint64_t>(phi / Detector::DeltaPhi);
  auto zBits = static_cast<uint64_t>((z - Detector::Zmin) / Detector::DeltaZ);

  return (rBits << 32) | (pBits << 16) | zBits;
}

//===----------------------------------------------------------------------===//
// Virtual surfaces based real geometry
//===----------------------------------------------------------------------===//

inline Double_t halfPhiSector(const BaseTpcSectorGeo &secGeo) {
  return M_PI / secGeo.SECTOR_COUNT_HALF;
}

inline Int_t getSector(const BaseTpcSectorGeo &secGeo,
                       Double_t x, Double_t y) {

  const Double_t halfPhiSec = halfPhiSector(secGeo);
  const Double_t phi = getPhi(x, y, halfPhiSec, 0);
  return static_cast<Int_t>(phi / (2.*halfPhiSec));
}

// Rotate vector (x, y) on phi
inline std::pair<Double_t, Double_t> rotateXY(Double_t x,
                                              Double_t y,
                                              Double_t phi) {
  return std::make_pair(
      x * std::cos(phi) - y * std::sin(phi),
      x * std::sin(phi) + y * std::cos(phi));
}

// Shift and rotate of (x, y) coordinates
// phi - is the rotation angle
// (x0, y0) - shift coordinates
inline std::pair<Double_t, Double_t> coordinateTransform(
    Double_t x,  Double_t y,
    Double_t x0, Double_t y0,
    Double_t phi) {
  return rotateXY(x - x0, y - y0, phi);
}

std::tuple <Int_t, Int_t, Int_t> getSectorRowPad (
    const BaseTpcSectorGeo &secGeo,
    Double_t x,
    Double_t y,
    Double_t phi) { // The -1*angle by which Acts rotates the geometry

  std::tie(x, y) = rotateXY(x, y, phi - halfPhiSector(secGeo));

  Int_t iSector = getSector(secGeo, x, y);

  // iSector * 30grad;
  Double_t iSectorPhi = const_cast<BaseTpcSectorGeo&>
      (secGeo).SectorAxisAngleRad(iSector); // FIXME:

  Double_t r0 = secGeo.YPADAREA_LOWEREDGE;
  Double_t x0 = r0 * std::cos(iSectorPhi);
  Double_t y0 = r0 * std::sin(iSectorPhi);

  // Angle of rotation of the coordinate axes (pad, row) with respect to (x, y)
  Double_t phiCoord = iSectorPhi - M_PI / 2.;

  auto [iPad, iRow] = coordinateTransform(x, y, x0, y0, -phiCoord);

  Double_t lengthPads0 = secGeo.YPADAREA_LENGTH[0];

  if (iRow < lengthPads0) {
    iRow = iRow / secGeo.PAD_HEIGHT[0];
    iPad = iPad / secGeo.PAD_WIDTH[0];
  } else {
    iRow = secGeo.ROW_COUNT[0] + (iRow - lengthPads0) / secGeo.PAD_HEIGHT[1];
    iPad = iPad / secGeo.PAD_WIDTH[1];
  }
  assert(iRow >= 0 && "The point is under the sector");
  Int_t maxPad = secGeo.PAD_COUNT[static_cast<Int_t>(iRow)];

  // Fixes the fact that the fast cluster finder algorithm generates
  // MC points not in the center, but on the boundary of the pad.
  // As a result of calculations the current pad can go beyond the
  // maximum pad for a small eps
  Double_t eps = 0.001;
  if ((maxPad < iPad) && (iPad < maxPad + eps)) {
    iPad -= eps;
  }
  if ((-maxPad - eps < iPad) && (iPad < -maxPad)) {
    iPad += eps;
  }

  // Shift to be >= 0
  Double_t iPadShifted = iPad + maxPad;

  return std::make_tuple(iSector,
                         static_cast<Int_t>(iRow),
                         static_cast<Int_t>(iPadShifted));
}

// Checks if row and pad number are within bounds and fix if they are not
void checkAndFixSecRowPad(const BaseTpcSectorGeo &secGeo,
                          const Acts::Vector3 &position,
                          Int_t sector,
                          Int_t &row, Int_t &pad) {
  Double_t x = position[0];
  Double_t y = position[1];
  Double_t z = position[2];

  std::string msgPrefix =
      "MpdTpcDetector.cxx::checkAndFixSecRowPad() ERROR: point " +
      std::to_string(x) + ", " +
      std::to_string(y) + ", " +
      std::to_string(z);
  std::string msgMiddle =
      "sector = " + std::to_string(sector) + ", " +
      "row = "    + std::to_string(row)    + ", " +
      "pad = "    + std::to_string(pad)    + ". ";
  Int_t maxRows = secGeo.ROW_COUNT[0] + secGeo.ROW_COUNT[1];
  if (row < 0) {
    row = 0;
    std::cout << msgPrefix << " is under first row. " <<
        msgMiddle << "Setting row = " << row << std::endl;
  } else if (row >= maxRows) {
    row = maxRows - 1;
    std::cout << msgPrefix << " is far from last row. " <<
        msgMiddle << "Setting row = " << row << std::endl;
  }
  Int_t maxPads = 2*(secGeo.PAD_COUNT[row]);
  if (pad < 0) {
    pad = 0;
    std::cout << msgPrefix << " goes beyond pads. " <<
        msgMiddle << "Setting pad = " << pad << std::endl;
  } else if (pad >= maxPads) {
    pad = maxPads - 1;
    std::cout << msgPrefix << " goes beyond pads. " <<
        msgMiddle << "Setting pad = " << pad << std::endl;
  }
}

inline uint64_t toBitsForSector(const BaseTpcSectorGeo &secGeo,
                                const Acts::Vector3 &position,
                                Acts::ActsScalar deltaPhi) {
  // From mm to cm
  Double_t x = toRootLength(position[0]);
  Double_t y = toRootLength(position[1]);

  Double_t z;

  if (Detector::getZFromSectorGeo) {
    z = toRootLength(position[2]);
    assert(-secGeo.Z_MAX <= z && z <= secGeo.Z_MAX);
  } else {
    z = position[2];
    assert(Detector::Zmin <= z && z <= Detector::Zmax);
  }
  auto [iSector, iRow, iPad] = getSectorRowPad(secGeo, x, y, deltaPhi);

  checkAndFixSecRowPad(secGeo, position, iSector, iRow, iPad);

  auto sBits = static_cast<uint64_t>(iSector);
  auto rBits = static_cast<uint64_t>(iRow);
  auto pBits = static_cast<uint64_t>(iPad);

  const auto deltaZ = Detector::getZFromSectorGeo
      ? secGeo.Z_MAX : Detector::DeltaZ;

  const auto zMin = Detector::getZFromSectorGeo
      ? -secGeo.Z_MAX : Detector::Zmin;

  auto zBits = static_cast<uint64_t>((z - zMin) / deltaZ);

  return (sBits << 48) | (rBits << 32) | (pBits << 16) | zBits;
}

void addSensorsPad(const BaseTpcSectorGeo &secGeo,
                   TGeoVolume *sensorVol, TGeoVolume *gasVolume,
                   Double_t phi, Int_t iRow, TGeoRotation *rot,
                   Int_t &iSensor) {
  const Int_t nRowsInner = secGeo.ROW_COUNT[0];

  const Double_t yPadAreaLowerEdge = secGeo.YPADAREA_LOWEREDGE;
  const Double_t h0 = secGeo.PAD_HEIGHT[0]; // cm
  const Double_t h1 = secGeo.PAD_HEIGHT[1]; // cm
  const Double_t x0 = (iRow < nRowsInner) ?
      yPadAreaLowerEdge + h0 * (iRow + 0.5) :
      yPadAreaLowerEdge + h0 * nRowsInner + h1 * (iRow - nRowsInner + 0.5);

  const Double_t w = (iRow < nRowsInner) ? secGeo.PAD_WIDTH[0] :
                                           secGeo.PAD_WIDTH[1];
  const Int_t nPadsInRow = secGeo.PAD_COUNT[iRow];

  for (Int_t iPadLeftRight = -1; iPadLeftRight < 2; iPadLeftRight +=2) {
    for (size_t iPad = 0; iPad < nPadsInRow; iPad++) {

      Double_t y0 = iPadLeftRight * w * (iPad + 0.5);
      Double_t z0 = 0;

      auto [xRotated0, yRotated0] = rotateXY(x0, y0, phi);
      auto loc = new TGeoCombiTrans(xRotated0, yRotated0, z0, rot);

      gasVolume->AddNode(sensorVol, iSensor++, loc);
    }
  }
}

void addSensorsRowPad(const BaseTpcSectorGeo &secGeo,
                      TGeoVolume *sensorVol0,
                      TGeoVolume *sensorVol1,
                      TGeoVolume *gasVolume,
                      Double_t phi,
                      TGeoRotation *rot,
                      Int_t &iSensor) {
  const Int_t nRows = secGeo.ROW_COUNT[0] + secGeo.ROW_COUNT[1];

  for (size_t iRow = 0; iRow < nRows; iRow++) {
    auto curVol = (iRow < secGeo.ROW_COUNT[0]) ? sensorVol0 :
                                                 sensorVol1;
    addSensorsPad(secGeo, curVol, gasVolume,
                  phi, iRow, rot, iSensor);
  }
}

void addSensorsSecRowPad(const BaseTpcSectorGeo &secGeo,
                         TGeoVolume *sensorVol0,
                         TGeoVolume *sensorVol1,
                         TGeoVolume *gasVolume) {
  const size_t nSectors = secGeo.SECTOR_COUNT_HALF;

  TGeoRotation *rot;
  Int_t iSensor = 0;
  for (size_t iSector = 0; iSector < nSectors; iSector++) {
    // Setting phi = iSector * 30_grad - 15_grad
    Double_t phi =
        const_cast<BaseTpcSectorGeo&>(secGeo).SectorAxisAngleRad(iSector) + // FIXME
        secGeo.SECTOR_PHI0_RAD;

    Double_t phiDegree = phi * TMath::RadToDeg();

    auto rotName = makeName("tpc_rotation_sec", iSector);
    rot = new TGeoRotation(rotName.c_str(), phiDegree, 0., 0.);

    addSensorsRowPad(secGeo, sensorVol0, sensorVol1, gasVolume, phi, rot, iSensor);
  }
}

void addSectorSensors(const BaseTpcSectorGeo &secGeo,
                      TGeoManager *geoManager,
                      TGeoVolume *gasVolume) {
  const Double_t dZ = Detector::getZFromSectorGeo
      ? secGeo.Z_MAX : toRootLength(Detector::DeltaZ);

  Int_t iVolume = 0;
  TGeoMedium *medium = gasVolume->GetMedium();

  const Double_t w0 = secGeo.PAD_WIDTH[0];  // cm
  const Double_t h0 = secGeo.PAD_HEIGHT[0]; // cm
  TGeoVolume *sensorVol0 = makeSensorVolume(
      geoManager,
      medium,
      0.5 * h0,
      0.5 * w0,
      0.5 * dZ,
      iVolume++);

  const Double_t w1 = secGeo.PAD_WIDTH[1];  // cm
  const Double_t h1 = secGeo.PAD_HEIGHT[1]; // cm
  TGeoVolume *sensorVol1 = makeSensorVolume(
      geoManager,
      medium,
      0.5 * h1,
      0.5 * w1,
      0.5 * dZ,
      iVolume++);

  assert(Detector::HasTwoNodes && "addSectorSensors() Detector::HasTwoNodes must be true");

  addSensorsSecRowPad(secGeo, sensorVol0, sensorVol1, gasVolume);
}

//===----------------------------------------------------------------------===//
// Class members and additional functions
//===----------------------------------------------------------------------===//

inline uint64_t toBits(const BaseTpcSectorGeo &secGeo,
                       const Acts::Vector3 &position) {
  auto bits = (Detector::geometryType == Detector::sectorBased) ?
      toBitsForSector(secGeo, position, halfPhiSector(secGeo)) :
      toBitsForCylinder(position, 0.5 * Detector::DeltaPhi);
  return bits;
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

void fillSurfaces(
    const BaseTpcSectorGeo &secGeo,
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
      auto bits = toBits(secGeo, center);
      assert(surfaces.find(bits) == surfaces.end()
          && "There are surfaces with the same bit keys");
      surfaces.insert({bits, surface});
    }
  }
}

Detector::Detector(const BaseTpcSectorGeo &secGeo,
                   const std::string &rootFile,
                   const std::string &jsonFile,
                   Acts::Logging::Level level):
  m_secGeo(secGeo),
  m_logger(Acts::getDefaultLogger("Detector", level)) {
  geoEditor = [this](auto *gm){ return editGeometry(gm); };

  m_config.surfaceLogLevel = level;
  m_config.layerLogLevel = level;
  m_config.volumeLogLevel = level;

  m_config.fileName = rootFile;
  m_config.readJson(jsonFile);
}

Bool_t Detector::editGeometry(TGeoManager *geoManager) {
  TGeoVolume *gasVolume = geoManager->GetVolume(Detector::GasVolume);
  if (!gasVolume) {
    ACTS_ERROR("Volume '" << Detector::GasVolume << "' not found");
    return false;
  }

  ACTS_DEBUG("TPC gas volume: ");
  gasVolume->Print();

  if (Detector::geometryType == Detector::sectorBased) {
    addSectorSensors(m_secGeo, geoManager, gasVolume);
  } else {
    addCylinderSensors(geoManager, gasVolume);
  }

  geoManager->CloseGeometry();
  gasVolume->CheckOverlaps();

  return true;
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

std::shared_ptr<const Acts::Surface> Detector::getSurface(
    const Acts::GeometryContext &gcontext,
    const Acts::Vector3 &position) {

  if (m_surfaces.empty()) {
    fillSurfaces(m_secGeo, gcontext, m_trackingGeometry, m_surfaces);
    ACTS_DEBUG("Constructed " << m_surfaces.size() << " surfaces");
  }

  ACTS_VERBOSE("Finding a surface for " << position);

  auto bits = toBits(m_secGeo, position);
  auto iter = m_surfaces.find(bits);
  assert(iter != m_surfaces.end() && "Surface not found");

  return iter->second;
}

} // namespace Mpd::Tpc

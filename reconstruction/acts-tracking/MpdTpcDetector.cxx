// This file is a part of the NICA project.
//
// Copyright (C) 2022 JINR

#include "MpdTpcDetector.h"
#include "geometry/BaseTpcSectorGeo.h"

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
                               Acts::ActsScalar deltaPhi) {
  auto phi = std::acos(x / std::hypot(x, y));
  if (y < 0) { phi = 2*M_PI - phi; }
  phi += deltaPhi;

  const auto eps = 0.001;
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
  const auto phi = getPhi(x, y, deltaPhi) + eps;
  assert(0 <= phi && phi < 2*M_PI);

  auto rBits = static_cast<uint64_t>((r - Detector::Rmin) / Detector::DeltaR);
  auto pBits = static_cast<uint64_t>(phi / Detector::DeltaPhi);
  auto zBits = static_cast<uint64_t>((z - Detector::Zmin) / Detector::DeltaZ);

  return (rBits << 32) | (pBits << 16) | zBits;
}

//===----------------------------------------------------------------------===//
// Virtual surfaces based real geometry
//===----------------------------------------------------------------------===//

//inline Acts::ActsScalar getPhi(Acts::ActsScalar x, Acts::ActsScalar y) {
//  auto phi = std::acos(x / std::hypot(x, y));
//  if (y < 0) { phi = 2*M_PI - phi; }
//  return phi;
//}

int getSector(const BaseTpcSectorGeo &secGeo, double x, double y) {
  int nSectors = secGeo.SECTOR_COUNT_HALF;

  const double phiSector = 2*M_PI / nSectors;
  const double phi = getPhi(x, y, phiSector/2.);

//  double phiForSector = phi + phiSector/2.;

//  if (phiForSector < 0) {
//    phiForSector += 2*M_PI;
//  } else if (phiForSector >= 2*M_PI) {
//    phiForSector -= 2*M_PI;
//  }

  return /*phiForSector*/ phi / phiSector;
}

void getSectorRowPad(const BaseTpcSectorGeo &secGeo,
                     double x,
                     double y,
                     int &sector,
                     double &row,
                     double &pad) {
  int iSector = getSector(secGeo, x, y);

// iSector * 30grad;
  double iSectorPhi = const_cast<BaseTpcSectorGeo&>(secGeo).SectorAxisAngleRad(iSector); // FIXME:

// distance from the Z axis to sector with pads
  double yPadAreaLowerEdge = secGeo.YPADAREA_LOWEREDGE;

  double x0 = yPadAreaLowerEdge * std::cos(iSectorPhi);
  double y0 = yPadAreaLowerEdge * std::sin(iSectorPhi);

// Angle of rotation of the coordinate axes (pad, row) with respect to (x, y)
  double phiCoord = iSectorPhi - M_PI / 2.;

  double iPad = (x - x0) * std::cos(-phiCoord) - (y - y0) * std::sin(-phiCoord); // FIXME
  double iRow = (x - x0) * std::sin(-phiCoord) + (y - y0) * std::cos(-phiCoord);

  double lengthPads0 = secGeo.YPADAREA_LENGTH[0];

  if (iRow < lengthPads0) {
    iRow = iRow / secGeo.PAD_HEIGHT[0] ;
    iPad = iPad / secGeo.PAD_WIDTH[0];
  } else {
    iRow = secGeo.ROW_COUNT[0] + (iRow - lengthPads0) / secGeo.PAD_HEIGHT[1] ;
    iPad = iPad / secGeo.PAD_WIDTH[1];
  }

  int maxRows = secGeo.ROW_COUNT[0] + secGeo.ROW_COUNT[1];
  if (iRow < 0) {
    iRow = 0.5;
    std::cout << "MpdTpcDetector.cxx::getSectorRowPad() ERROR: " <<
        "point " << x << ", " << y << "z, is under first row. " <<
        "Setting row = " << iRow << std::endl;
  } else if (iRow >= maxRows) {
    iRow = maxRows - 0.5;
    std::cout << "MpdTpcDetector.cxx::getSectorRowPad() ERROR: " <<
        "point " << x << ", " << y << ", z is far from last row. " <<
        "Setting row = " << iRow << std::endl;
  }

// Additional pads in row, as LMEM algorithm generates points outside pads in row
  int extraPads = 0;

// Shift to be >= 0
  double iPadShifted = iPad + secGeo.PAD_COUNT[(int)iRow] + extraPads;

  int maxPads = secGeo.PAD_COUNT[(int)iRow] + extraPads;
  if (fabs(iPad) >= maxPads) {
    int sign = (0 <= iPad) - (iPad < 0);
    iPad = sign * (maxPads - 0.5);
    iPadShifted = iPad + maxPads;
    std::cout << "MpdTpcDetector.cxx::getSectorRowPad() ERROR: " <<
        "point " << x << ", " << y << ", z goes beyond pads. " <<
        "row = " << iRow << ". " <<
        "Setting pad = " << iPad << std::endl;
  }

  sector = iSector;
  row    = iRow;
  pad    = iPadShifted;
}

// rotate vector (xIn, yIn) on phi
inline void rotateXY(double xIn,
                     double yIn,
                     double phi,
                     double &xOut,
                     double &yOut) {
  double x = xIn * std::cos(phi) - yIn * std::sin(phi);
  double y = xIn * std::sin(phi) + yIn * std::cos(phi);
  xOut = x;
  yOut = y;
}

inline uint64_t toBitsForSector(const BaseTpcSectorGeo &secGeo,
                                const Acts::Vector3 &position,
                                const Acts::ActsScalar deltaPhi) {
  // from mm to cm
  double x = toRootLength(position[0]);
  double y = toRootLength(position[1]);
  double z;

  if (Detector::getZFromSectorGeo) {
    z = toRootLength(position[2]);
    assert(-secGeo.Z_MAX <= z && z <= secGeo.Z_MAX);
  } else {
    z = position[2];
    assert(Detector::Zmin <= z && z <= Detector::Zmax);
  }

  rotateXY(x, y, deltaPhi, x, y);

  int iSector;
  double iRow;
  double iPad;

  getSectorRowPad(secGeo, x, y, iSector, iRow, iPad);

  auto sBits = static_cast<uint64_t>(iSector);
  auto rBits = static_cast<uint64_t>(iRow);
  auto pBits = static_cast<uint64_t>(iPad);

  const auto deltaZ = Detector::getZFromSectorGeo
      ? secGeo.DRIFT_LENGTH : Detector::DeltaZ;
  const auto zMax = Detector::getZFromSectorGeo
      ? secGeo.Z_MAX : Detector::Zmax;

  auto zBits = static_cast<uint64_t>((z + zMax) / deltaZ);

  return (sBits << 48) | (rBits << 32) | (pBits << 16) | zBits;
}

void addSensorsPad(const BaseTpcSectorGeo &secGeo,
                   TGeoVolume *sensorVol0, TGeoVolume *sensorVol1, TGeoVolume *gasVolume,
                   int iSector, double phi,
                   int iRow, TGeoRotation *rot,
                   int &iSensor, int extraPads) {
  const int nRowsInner = secGeo.ROW_COUNT[0];

  const double yPadAreaLowerEdge = secGeo.YPADAREA_LOWEREDGE;
  const double h0 = secGeo.PAD_HEIGHT[0]; // cm
  const double h1 = secGeo.PAD_HEIGHT[1]; // cm
  const double x0 = (iRow < nRowsInner) ?
      yPadAreaLowerEdge + h0 * (iRow + 0.5) :
      yPadAreaLowerEdge + h0 * nRowsInner + h1 * (iRow - nRowsInner + 0.5);

  const double w = (iRow < nRowsInner) ? secGeo.PAD_WIDTH[0] :
                                         secGeo.PAD_WIDTH[1];

  const auto nPadsInRow = secGeo.PAD_COUNT[iRow];

  for (int iPadLeftRight = -1; iPadLeftRight < 2; iPadLeftRight +=2) {
    for (auto iPad = 0; iPad < nPadsInRow + extraPads ; iPad++) {
      int maxPads = nPadsInRow + extraPads;

      double y0 = iPadLeftRight * w * (iPad + 0.5);
      double z0 = 0;

      double xRotated0;
      double yRotated0;
      rotateXY(x0, y0, phi, xRotated0, yRotated0);

      auto loc = new TGeoCombiTrans(xRotated0, yRotated0, z0, rot);

      gasVolume->AddNode(sensorVol0, iSensor++, loc);
    }
  }
}

void addSensorsRowPad(const BaseTpcSectorGeo &secGeo,
                      TGeoVolume *sensorVol0,
                      TGeoVolume *sensorVol1,
                      TGeoVolume *gasVolume,
                      int iSector,
                      double phi,
                      TGeoRotation *rot,
                      int &iSensor) {
  const int nRowsInner = secGeo.ROW_COUNT[0];
  const int nRowsOuter = secGeo.ROW_COUNT[1];

  const int extraPads = 0;

  for (auto iRow = 0; iRow < nRowsInner + nRowsOuter; iRow++) {
    addSensorsPad(secGeo, sensorVol0, sensorVol1, gasVolume,
                  iSector, phi, iRow, rot, iSensor, extraPads);
  }
}

void addSensorsSecRowPad(const BaseTpcSectorGeo &secGeo,
                         TGeoVolume *sensorVol0,
                         TGeoVolume *sensorVol1,
                         TGeoVolume *gasVolume) {
  const int nSectors = secGeo.SECTOR_COUNT_HALF;

  TGeoRotation *rot;
  int iSensor = 0;

  for (int iSector = 0; iSector < nSectors; iSector++) {

    double phi = const_cast<BaseTpcSectorGeo&>(secGeo).SectorAxisAngleRad(iSector) +  // = iSector * 30_grad // FIXME
                 secGeo.SECTOR_PHI0_RAD;               //   - 15_grad
    double phiDegree = phi * TMath::RadToDeg();

    auto rotName = makeName("tpc_rotation_sec", iSector);
    rot = new TGeoRotation(rotName.c_str(), phiDegree, 0., 0.);

    addSensorsRowPad(secGeo, sensorVol0, sensorVol1, gasVolume, iSector, phi, rot, iSensor);
  }
}

void addSectorSensors(const BaseTpcSectorGeo &secGeo,
                      TGeoManager *geoManager,
                      TGeoVolume *gasVolume) {
  auto dZ = Detector::getZFromSectorGeo ? secGeo.DRIFT_LENGTH :
                                          toRootLength(Detector::DeltaZ);
  Int_t iVolume = 0;
  auto *medium = gasVolume->GetMedium();

  const double w0 = secGeo.PAD_WIDTH[0];  // cm
  const double h0 = secGeo.PAD_HEIGHT[0]; // cm
  TGeoVolume *sensorVol0 = makeSensorVolume(
      geoManager,
      medium,
      0.5 * h0,
      0.5 * w0,
      0.5 * dZ,
      iVolume++ );

  const double w1 = secGeo.PAD_WIDTH[1];  // cm
  const double h1 = secGeo.PAD_HEIGHT[1]; // cm
  TGeoVolume *sensorVol1 = makeSensorVolume(
      geoManager,
      medium,
      0.5 * h1,
      0.5 * w1,
      0.5 * dZ,
      iVolume++ );

  assert(Detector::HasTwoNodes && "addSectorSensors() Detector::HasTwoNodes must be true");

  addSensorsSecRowPad(secGeo, sensorVol0, sensorVol1, gasVolume);
}

//===----------------------------------------------------------------------===//
// Class members and additional functions
//===----------------------------------------------------------------------===//

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

      auto bits = Detector::useBaseTpcSectorGeo ? // FIXME:
          toBitsForSector(secGeo, center, 0) :
          toBitsForCylinder(center, 0.5 * Detector::DeltaPhi);
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

bool Detector::editGeometry(TGeoManager *geoManager) {
  TGeoVolume *gasVolume = geoManager->GetVolume(Detector::GasVolume);
  if (!gasVolume) {
    ACTS_ERROR("Volume '" << Detector::GasVolume << "' not found");
    return false;
  }

  ACTS_DEBUG("TPC gas volume: ");
  gasVolume->Print();

  if (Detector::useBaseTpcSectorGeo) {
    addSectorSensors(m_secGeo, geoManager, gasVolume);
  } else {
    addCylinderSensors(geoManager, gasVolume);
  }

  geoManager->CloseGeometry();
  gasVolume->CheckOverlaps();

  if (m_logger->doPrint(Acts::Logging::DEBUG)) {
    geoManager->Export("tpc_acts_tracking.root");
  }

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
  auto bits = Detector::useBaseTpcSectorGeo ?
      toBitsForSector  (m_secGeo, position, 0) : // FIXME:
      toBitsForCylinder(position, 0.5 * Detector::DeltaPhi);
  auto iter = m_surfaces.find(bits);
  assert(iter != m_surfaces.end() && "Surface not found");

  return iter->second;
}

} // namespace Mpd::Tpc

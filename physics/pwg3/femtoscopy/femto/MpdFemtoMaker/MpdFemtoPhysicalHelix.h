/**
 * \class MpdFemtoPhysicalHelix
 * \author Grigory Nigmatkulov, May 07 2018
 *
 * Parametrization of a physical helix (modification of StPhysicalHelix).
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * \email nigmatkulov@gmail.com
 */

#ifndef MpdFemtoPhysicalHelix_h
#define MpdFemtoPhysicalHelix_h

// ROOT headers
#include "TVector3.h"

// PicoDst headers
#include "MpdFemtoHelix.h"

//_________________

class MpdFemtoPhysicalHelix : public MpdFemtoHelix {
 public:
  /// Default constructor
  MpdFemtoPhysicalHelix();
  /// Constructor with momentum, origin, signed Magnetic Field
  /// and Charge of particle (+/- 1)
  MpdFemtoPhysicalHelix(const TVector3&,
			const TVector3&,
			double, double);
  /// Constructor with Curvature, dip angle, phase, origin, h
  MpdFemtoPhysicalHelix(double, double, double,
			const TVector3&, Int_t h = -1);
  /// Destructor
  ~MpdFemtoPhysicalHelix();

  /// Requires: signed Magnetic Field
  TVector3 momentum(double) const; // returns the momentum at origin
  TVector3 momentumAt(double, double) const; // returns momemtum at S
  Int_t charge(double) const; // returns charge of particle
  /// 2d DCA to x,y point signed relative to curvature
  double curvatureSignedDistance(double x, double y);
  /// 2d DCA to x,y point signed relative to rotation
  double geometricSignedDistance(double x, double y);
  /// 3d DCA to 3d point signed relative to curvature
  double curvatureSignedDistance(const TVector3&);
  /// 3d DCA to 3d point signed relative to rotation
  double geometricSignedDistance(const TVector3&);

  ClassDef(MpdFemtoPhysicalHelix, 1)
};

#endif

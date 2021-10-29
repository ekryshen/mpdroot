/**
 * \class MpdFemtoTriplet
 * \brief The class allows to perform three-particle analysis
 *
 * The MpdFemtoTriplet class allows to perform three-particle
 * analysis and calculate quantities like: qInv, mInv, etc...
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * \email nigmatkulov@gmail.com
 */

#ifndef MpdFemtoTriplet_h
#define MpdFemtoTriplet_h

// C++ headers
#include <utility>

// MpdFemtoMaker headers
#include "MpdFemtoParticle.h"
#include "MpdFemtoTypes.h"

// ROOT headers
#include "TLorentzVector.h"

//_________________
class MpdFemtoTriplet {
 public:
  /// Default constructor
  MpdFemtoTriplet();
  /// Constructor with three particles
  MpdFemtoTriplet(MpdFemtoParticle*, MpdFemtoParticle*, MpdFemtoParticle*);
  /// Copy constructor
  MpdFemtoTriplet(const MpdFemtoTriplet&);
  /// Default destructor
  virtual ~MpdFemtoTriplet();

  /// Retrieve the first track
  MpdFemtoParticle* track1() const {
    return mTrack1;
  }
  /// Retrieve the second track
  MpdFemtoParticle* track2() const {
    return mTrack2;
  }
  /// Retrieve the third track
  MpdFemtoParticle* track3() const {
    return mTrack3;
  }
  /// Set first track

  void setTrack1(const MpdFemtoParticle* trkPtr) {
    mTrack1 = (MpdFemtoParticle*) trkPtr;
  }
  /// Set second track
  void setTrack2(const MpdFemtoParticle* trkPtr) {
    mTrack2 = (MpdFemtoParticle*) trkPtr;
  }
  /// Set third track
  void setTrack3(const MpdFemtoParticle* trkPtr) {
    mTrack3 = (MpdFemtoParticle*) trkPtr;
  }

  /// Four-momentum of three particles
  TLorentzVector fourMomentum() const;
  /// Relative momentum of three particles
  double qInv() const;
  /// Relative momentum the first and second particles
  double qInv12() const;
  /// Relative momentum the third and second particles
  double qInv23() const;
  /// Relative momentum the first and third particles
  double qInv31() const;
  /// Half of transverse momentum of three particles
  double kT() const;
  /// Invariant mass of three particles
  double mInv() const;

  /// Track-splitting quantity of three particles
  double quality() const;
  /// Track-splitting quantity of three particles
  double splittinLevel() const {
    return quality();
  }

  // the following two methods calculate the "nominal" separation of the tracks
  // at the inner field cage (EntranceSeparation) and when they exit the TPC,
  // which may be at the outer field cage, or at the endcaps.
  // "nominal" means that the tracks are assumed to start at (0,0,0).  Making this
  // assumption is important for the Event Mixing-- it is not a mistake

  /// Nominal TPC exit separation
  double nominalTpcExitSeparation() const;
  /// Nominal TPC entrance separation
  double nominalTpcEntranceSeparation() const;
  /// Nominal TPC average separation
  double nominalTpcAverageSeparation() const;

 private:

  /// The first particle from triplet
  MpdFemtoParticle* mTrack1;
  /// The second particle from triplet
  MpdFemtoParticle* mTrack2;
  /// The third particle from triplet
  MpdFemtoParticle* mTrack3;

  ClassDef(MpdFemtoTriplet, 1)
};

#endif // MpdFemtoTriplet_h

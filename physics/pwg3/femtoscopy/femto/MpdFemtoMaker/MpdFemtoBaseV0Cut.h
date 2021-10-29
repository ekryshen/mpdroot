/**
 * \class MpdFemtoBaseV0Cut
 * \brief Base class for V0 cuts
 *
 * The pure virtual base class for the V0 cut. All V0 cuts
 * must inherit from this one.
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * \email nigmatkulov@gmail.com
 */

#ifndef MpdFemtoBaseV0Cut_h
#define MpdFemtoBaseV0Cut_h

// MpdFemtoMaker headers
// Infrastructure
#include "MpdFemtoTypes.h"
#include "MpdFemtoV0.h"
// Base
#include "MpdFemtoBaseParticleCut.h"

//_________________
class MpdFemtoBaseV0Cut : public MpdFemtoBaseParticleCut {
 public:
  /// Default constructor
  MpdFemtoBaseV0Cut() {
    /* empty */
  }
  /// Copy constructor
  MpdFemtoBaseV0Cut(const MpdFemtoBaseV0Cut& copy);
  /// Assignment operator
  MpdFemtoBaseV0Cut& operator=(const MpdFemtoBaseV0Cut& c);
  /// Default destructor
  virtual ~MpdFemtoBaseV0Cut() {
    /* empty */
  }

  /// Returns true if cut has been passed and false if not
  virtual bool pass(const MpdFemtoV0*) = 0;

  /// Return V0 type
  virtual MpdFemtoParticleType type() {
    return hbtV0;
  }
  /// Clone V0 cut
  virtual MpdFemtoBaseV0Cut* clone() {
    return nullptr;
  }

  ClassDef(MpdFemtoBaseV0Cut, 0)
};

//_________________
inline MpdFemtoBaseV0Cut::MpdFemtoBaseV0Cut(const MpdFemtoBaseV0Cut& c) : MpdFemtoBaseParticleCut(c) {
  /* empty */
}

//_________________
inline MpdFemtoBaseV0Cut& MpdFemtoBaseV0Cut::operator=(const MpdFemtoBaseV0Cut& c) {
  if (this != &c) {
    MpdFemtoBaseParticleCut::operator=(c);
  }
  return *this;
}

#endif // #define MpdFemtoV0Cut_h

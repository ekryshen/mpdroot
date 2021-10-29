/**
 * \class MpdFemtoBaseKinkCut
 * \brief Base class for kink cuts
 *
 * The the pure virtual base class for the kink cut. All kink cuts
 * must inherit from this one
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * \email nigmatkulov@gmail.com
 */

#ifndef MpdFemtoBaseKinkCut_h
#define MpdFemtoBaseKinkCut_h

// MpdFemtoMaker headers
// Base
#include "MpdFemtoBaseParticleCut.h"
// Infrastructure
#include "MpdFemtoTypes.h"
#include "MpdFemtoKink.h"

//_________________
class MpdFemtoBaseKinkCut : public MpdFemtoBaseParticleCut {
 public:
  /// Default constructor

  MpdFemtoBaseKinkCut() {
    /* empty */
  }
  /// Copy constructor
  MpdFemtoBaseKinkCut(const MpdFemtoBaseKinkCut& copy);
  /// Assignment operator
  MpdFemtoBaseKinkCut& operator=(const MpdFemtoBaseKinkCut& copy);
  /// Default destructor
  virtual ~MpdFemtoBaseKinkCut() {
    /* empty */
  }

  /// Returns true is cut has been passed, and false if not
  virtual bool pass(const MpdFemtoKink*) = 0;
  /// Return kink type
  virtual MpdFemtoParticleType type() {
    return hbtKink;
  }
  /// Clone kink cut
  virtual MpdFemtoBaseKinkCut* clone() {
    return nullptr;
  }

  ClassDef(MpdFemtoBaseKinkCut, 0)
};

//_________________
inline MpdFemtoBaseKinkCut::MpdFemtoBaseKinkCut(const MpdFemtoBaseKinkCut& c) : MpdFemtoBaseParticleCut(c) {
  /* empty */
}

//_________________
inline MpdFemtoBaseKinkCut& MpdFemtoBaseKinkCut::operator=(const MpdFemtoBaseKinkCut& c) {
  if (this != &c) {
    MpdFemtoBaseParticleCut::operator=(c);
  }
  return *this;
}

#endif // #define MpdFemtoKinkCut_h

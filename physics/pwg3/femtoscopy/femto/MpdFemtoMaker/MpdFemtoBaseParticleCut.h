/**
 * \class MpdFemtoBaseParticleCut
 * \brief The pure virtual base class for the particle cut.
 *
 * All particle cuts must inherit from this one
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * \email nigmatkulov@gmail.comm
 */

#ifndef MpdFemtoBaseParticleCut_h
#define MpdFemtoBaseParticleCut_h

// C++ headers
#include <iostream>

// MpdFemtoMaker headers
#include "MpdFemtoTypes.h"
#include "MpdFemtoCutMonitorHandler.h"

// ROOT headers
#include "TList.h"
#include "TObjString.h"
#include "TString.h"

// Forward declaration
class MpdFemtoBaseAnalysis;

//_________________
class MpdFemtoBaseParticleCut : public MpdFemtoCutMonitorHandler {
 public:
  /// Default constructor
  MpdFemtoBaseParticleCut();
  /// Copy constructor
  MpdFemtoBaseParticleCut(const MpdFemtoBaseParticleCut& copy);
  /// Assignment operator
  MpdFemtoBaseParticleCut& operator=(const MpdFemtoBaseParticleCut& c);
  /// Default destructor
  virtual ~MpdFemtoBaseParticleCut() {
    /* no-op */
  }

  /// User-written method to return string describing cuts
  virtual MpdFemtoString report() = 0;
  /// User-written list of settings which is stored in the result file
  virtual TList *listSettings();

  /// Return mass of the particle to be selected
  double mass() {
    return mMass;
  }
  /// Set mass of the particle to be selected
  virtual void setMass(const double& mass) {
    mMass = mass;
  }

  /// Declare event start
  virtual void eventBegin(const MpdFemtoEvent*) {
    /* no-op */
  }
  /// Declare event end
  virtual void eventEnd(const MpdFemtoEvent*) {
    /* no-op */
  }
  /// Clone base particle cut
  virtual MpdFemtoBaseParticleCut* clone() {
    return nullptr;
  }

  /// Return particle type
  virtual MpdFemtoParticleType type() = 0;

  /// The following allows "back-pointing" from the CorrFctn
  /// to the "parent" Analysis
  friend class MpdFemtoBaseAnalysis;
  /// Return a pointer to the analysis
  MpdFemtoBaseAnalysis* hbtAnalysis() {
    return mBaseAnalysis;
  }
  /// Set analysis
  void setAnalysis(MpdFemtoBaseAnalysis* ana) {
    mBaseAnalysis = ana;
  }

 protected:

  /// Particle mass
  double mMass;
  /// Pointer to the base analysis
  MpdFemtoBaseAnalysis* mBaseAnalysis; //!<!

  ClassDef(MpdFemtoBaseParticleCut, 0)
};

//_________________
inline MpdFemtoBaseParticleCut::MpdFemtoBaseParticleCut() : MpdFemtoCutMonitorHandler(), mMass(0), mBaseAnalysis(nullptr) {
  /* empty */
}

//_________________
inline MpdFemtoBaseParticleCut::MpdFemtoBaseParticleCut(const MpdFemtoBaseParticleCut& c) :
			       MpdFemtoCutMonitorHandler(c), mMass(c.mMass), mBaseAnalysis(c.mBaseAnalysis) {
  /* empty */
}

//_________________
inline MpdFemtoBaseParticleCut& MpdFemtoBaseParticleCut::operator=(const MpdFemtoBaseParticleCut& c) {
  if (this != &c) {
    MpdFemtoCutMonitorHandler::operator=(c);
    mBaseAnalysis = c.mBaseAnalysis;
    mMass = c.mMass;
  }
  return *this;
}

//_________________
inline TList *MpdFemtoBaseParticleCut::listSettings() {
  TList *listOfSettings = new TList();
  listOfSettings->Add(new TObjString(Form("MpdFemtoBaseParticleCut::mass = %5.3f", mMass)));
  return listOfSettings;
}

#endif // #define MpdFemtoBaseParticleCut_h

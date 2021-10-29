/**
 * \class MpdFemtoBaseCorrFctn
 * \brief The pure-virtual base class for correlation functions
 *
 * All correlation function classes must inherit from this one.
 * This class has a optional pointers to the "parent" analysis and
 * a pair cut
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * \email nigmatkulov@gmail.com
 */

#ifndef MpdFemtoBaseCorrFctn_h
#define MpdFemtoBaseCorrFctn_h

// C++ headers
#include <iostream>

// MpdFemtoMaker headers
// Base
#include "MpdFemtoBaseAnalysis.h"
#include "MpdFemtoBasePairCut.h"
// Real
#include "MpdFemtoEvent.h"
#include "MpdFemtoPair.h"
#include "MpdFemtoTriplet.h"

// ROOT headers
#include "TList.h"

//_________________
class MpdFemtoBaseCorrFctn {
 public:
  /// Default constuctor
  MpdFemtoBaseCorrFctn();
  /// Copy constructor
  MpdFemtoBaseCorrFctn(const MpdFemtoBaseCorrFctn& copy);
  /// Assignment operator
  MpdFemtoBaseCorrFctn& operator=(const MpdFemtoBaseCorrFctn& copy);
  /// Default destructor

  virtual ~MpdFemtoBaseCorrFctn() {
    /* no-op */
  }

  /// Report
  virtual MpdFemtoString report() = 0;

  /// Add real pair
  virtual void addRealPair(MpdFemtoPair*);
  /// Add mixed pair
  virtual void addMixedPair(MpdFemtoPair*);

  /// Add real triplet
  virtual void addRealTriplet(const MpdFemtoTriplet*);
  /// Add mixed triplet
  virtual void addMixedTriplet(const MpdFemtoTriplet*);

  /// Not Implemented - Add particle with options
  virtual void addFirstParticle(MpdFemtoParticle *particle, bool mixing);
  /// Not Implemented - Add particle
  virtual void addSecondParticle(MpdFemtoParticle *particle);

  /// Event start
  virtual void eventBegin(const MpdFemtoEvent*) {
    /* no-op */
  }
  
  /// Event end
  virtual void eventEnd(const MpdFemtoEvent*) {
    /* no-op */
  }
  /// Finish method
  virtual void finish() = 0;

  /// Return output list
  virtual TList* getOutputList() = 0;

  /// Return correlation function clone
  virtual MpdFemtoBaseCorrFctn* clone() const = 0;
  
  /// Return pointer to the pair cut
  virtual MpdFemtoBasePairCut* getPairCut() {
    return mPairCut;
  }

  // The following allows "back-pointing" from the CorrFctn to the "parent" Analysis
  friend class MpdFemtoBaseAnalysis;
  
  /// Return a pointer to the analysis
  MpdFemtoBaseAnalysis* hbtAnalysis() {
    return mBaseAnalysis;
  }
  
  /// Set analysis
  void setAnalysis(MpdFemtoBaseAnalysis* ana) {
    mBaseAnalysis = ana;
  }
  
  /// Set pair cut
  void setPairSelectionCut(MpdFemtoBasePairCut *cut) {
    mPairCut = cut;
  }

 protected:
  /// Pointer to the base analysis
  MpdFemtoBaseAnalysis* mBaseAnalysis; //!
  /// Pointer to the base pair cut
  MpdFemtoBasePairCut* mPairCut; //!

 private:

  ClassDef(MpdFemtoBaseCorrFctn, 0)
};

#endif // #define MpdFemtoBaseCorrFctn_h

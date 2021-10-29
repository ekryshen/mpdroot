/**
 * \class MpdFemtoLikeSignAnalysis
 * \brief The class is a base class for like-sign analysis with z-binning
 *
 * The class provides posibility to perform femtoscopic analysis
 * using vertex binning
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * \email nigmatkulov@gmail.com
 */

#ifndef MpdFemtoLikeSignAnalysis_h
#define MpdFemtoLikeSignAnalysis_h

// MpdFemtoMaker headers
// Base
#include "MpdFemtoBaseAnalysis.h"
#include "MpdFemtoBaseEventCut.h"
#include "MpdFemtoBaseParticleCut.h"
#include "MpdFemtoBasePairCut.h"
#include "MpdFemtoBaseLikeSignCorrFctn.h"
// Infrastructure
#include "MpdFemtoTypes.h"
#include "MpdFemtoAnalysis.h"
#include "MpdFemtoCorrFctnCollection.h"

//_________________
class MpdFemtoLikeSignAnalysis : public MpdFemtoAnalysis {
 public:
  /// Constructor
  MpdFemtoLikeSignAnalysis(const unsigned int& bins = 20, const double& min = -100., const double& max = 100.);
  /// Copy constructor
  MpdFemtoLikeSignAnalysis(const MpdFemtoLikeSignAnalysis& copy);
  /// Copy constructor
  MpdFemtoLikeSignAnalysis& operator=(const MpdFemtoLikeSignAnalysis& copy);
  /// Default destructor
  virtual ~MpdFemtoLikeSignAnalysis();

  /// Function that calls processing
  virtual void processEvent(const MpdFemtoEvent*);
  /// Make report
  virtual MpdFemtoString report();
  /// Make report
  virtual MpdFemtoString Report() {
    return report();
  }
  /// Return number of events were lower than min z
  virtual unsigned int overflow() {
    return mOverFlow;
  }
  /// Return number of events were lower than min z
  virtual unsigned int Overflow() {
    return overflow();
  }
  /// Return number of events were higher than max z
  virtual unsigned int underflow() {
    return mUnderFlow;
  }
  /// Return number of events were higher than max z
  virtual unsigned int Underflow() {
    return mUnderFlow;
  }

 protected:
  /// Min/Max z-vertex position allowed to be processed
  double mVertexZ[2];
  /// Number of mixing bins in z-vertex in EventMixing Buffer
  unsigned int mVertexBins;
  /// Number of events encountered which had too large z-vertex
  unsigned int mOverFlow;
  /// Number of events encountered which had too small z-vertex
  unsigned int mUnderFlow;

  ClassDef(MpdFemtoLikeSignAnalysis, 0)
};

#endif // #define MpdFemtoLikeSignAnalysis_h

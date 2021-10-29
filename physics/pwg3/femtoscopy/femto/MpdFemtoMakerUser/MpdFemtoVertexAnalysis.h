/**
 * \class MpdFemtoVertexAnalysis
 * \brief Femtoscopic analysis which mixes events with vertex z-binning
 *
 * Femtoscopic analysis which mixes events with respect to the z-position
 * of the primary vertex
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * \email nigmatkulov@gmail.com
 */

#ifndef MpdFemtoVertexAnalysis_h
#define MpdFemtoVertexAnalysis_h

// MpdFemtoMaker headers
#include "MpdFemtoAnalysis.h"

//_________________
class MpdFemtoVertexAnalysis : public MpdFemtoAnalysis {
 public:
  /// Standard constructor
  MpdFemtoVertexAnalysis(unsigned int bins = 10, float zMin = -100., float zMax = +100.);
  /// Copy constructor
  MpdFemtoVertexAnalysis(const MpdFemtoVertexAnalysis& copy);
  /// Assignment constructor
  MpdFemtoVertexAnalysis& operator=(const MpdFemtoVertexAnalysis& copy);
  /// Destructor
  virtual ~MpdFemtoVertexAnalysis();

  /// Main processor
  virtual void processEvent(const MpdFemtoEvent* thisEvent);
  /// Returns reports of all cuts applied and correlation functions being done
  virtual MpdFemtoString report();
  /// Returns number of events that have z greater than the zMax
  virtual unsigned int overflow() {
    return mOverFlow;
  }
  /// Returns number of events that have z smaller than the zMin
  virtual unsigned int underflow() {
    return mUnderFlow;
  }

  /// Return list of cut settings for the analysis
  virtual TList* listSettings();
  /// Return a TList of objects to be written as output
  virtual TList* getOutputList();

 protected:

  /// min/max z-vertex position allowed to be processed
  float mVertexZ[2];
  /// Number of mixing bins in z-vertex in EventMixing Buffer
  unsigned int mVertexBins;
  /// Number of events encountered which had too large z-vertex
  unsigned int mOverFlow;
  /// Number of events encountered which had too small z-vertex
  unsigned int mUnderFlow;

  ClassDef(MpdFemtoVertexAnalysis, 0)
};

#endif

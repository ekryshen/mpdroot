/**
 * \class MpdFemtoReactionPlaneAnalysis
 * \brief Class that allows the azimuthal analysis
 *
 * Analysis with respect to the z position of the primary vertex and event
 * total multiplicity and uses only events in certain reaction plane angle bin
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * \email nigmatkulov@gmail.com
 */

#ifndef MpdFemtoReactionPlaneAnalysis_h
#define MpdFemtoReactionPlaneAnalysis_h

// MpdFemtoMaker headers
#include "MpdFemtoAnalysis.h"

// ROOT headers
#include "TList.h"

// Forward declaration
class MpdFemtoPicoEventCollectionVectorHideAway;

//_________________
class MpdFemtoReactionPlaneAnalysis : public MpdFemtoAnalysis {
 public:
  /// Default constructor
  MpdFemtoReactionPlaneAnalysis(unsigned int binsVertex = 10, double minVertex = -100., double maxVertex = +100.,
				unsigned int binsMult = 10, double minMult = -1.e9, double maxMult = +1.e9,
				unsigned short binsRP = 10, double minRP = -1.e9, double maxRP = +1.e9);
  /// Copy constructor
  MpdFemtoReactionPlaneAnalysis(const MpdFemtoReactionPlaneAnalysis& copy);
  /// Assignment operator
  MpdFemtoReactionPlaneAnalysis& operator=(const MpdFemtoReactionPlaneAnalysis& copy);
  /// Destructor
  virtual ~MpdFemtoReactionPlaneAnalysis();

  /// Process event
  virtual void processEvent(const MpdFemtoEvent* thisEvent);
  /// Reports of all cuts applied and correlation functions being done
  virtual MpdFemtoString report();
  /// Number of events that overflow Z max
  virtual unsigned int overflowVertexZ() const {
    return mOverFlowVertexZ;
  }
  /// Number of events that underflow Z min
  virtual unsigned int underflowVertexZ() const {
    return mUnderFlowVertexZ;
  }
  /// Number of events that overflow max multiplicity
  virtual unsigned int overflowMult() const {
    return mOverFlowMult;
  }
  /// Number of events that underflow min multiplicity
  virtual unsigned int underflowMult() const {
    return mUnderFlowMult;
  }
  /// Reaction plane angle
  double getCurrentReactionPlane() const {
    return mCurrentRP;
  }

  /// Return a TList of analysis settings.
  ///
  /// The TList comprises TObjStrings describing the settings provided by the
  /// MpdFemtoAnalysis::ListSettings class followed by all event-mixing
  /// binning parameters.
  virtual TList *listSettings();
  /// Obtain number of objects to be written as an output
  virtual TList *getOutputList();

protected:

  /// Min/max z-vertex position allowed to be processed
  double mVertexZ[2];
  /// Number of VERTEX mixing bins in z-vertex in EventMixing Buffer
  unsigned int mVertexZBins;
  /// Number of events encountered which had too large z-vertex
  unsigned int mOverFlowVertexZ;
  /// Number of events encountered which had too small z-vertex
  unsigned int mUnderFlowVertexZ;

  /// Min/max multiplicity allowed for event to be processed
  double mMult[2];
  /// Number of MULTIPLICITY mixing bins in z-vertex in EventMixing Buffer
  unsigned int mMultBins;
  /// Number of events encountered which had too large multiplicity
  unsigned int mOverFlowMult;
  /// Number of events encountered which had too small multiplicity
  unsigned int mUnderFlowMult;

  /// Number of reaction plane angle orientation bins
  unsigned short mRPBins;
  /// Reaction plane angle range
  double mRP[2];
  /// Reaction plane angle of the current event
  double mCurrentRP;

  ClassDef(MpdFemtoReactionPlaneAnalysis, 0)
};

#endif // #define MpdFemtoReactionPlaneAnalysis_h

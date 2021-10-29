/**
 * \class MpdFemtoPicoEventCollectionVectorHideAway
 * \brief A helper class for managing many mixing buffers
 *
 * A helper class for managing many mixing buffers with up to
 * three variables used for binning.
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * \email nigmatkulov@gmail.com
 */

#ifndef MpdFemtoPicoEventCollectionVectorHideAway_h
#define MpdFemtoPicoEventCollectionVectorHideAway_h

// C++ headers
#include <cmath>
#include <limits>
#include <float.h>

// MpdFemtoMaker headers
// Infrastructure
#include "MpdFemtoPicoEvent.h"
#include "MpdFemtoPicoEventCollection.h"
#include "MpdFemtoPicoEventCollectionVector.h"

//_________________
class MpdFemtoPicoEventCollectionVectorHideAway {
 public:
  /// Constructor
  MpdFemtoPicoEventCollectionVectorHideAway(int bx = 1, double lx = -FLT_MAX, double ux = FLT_MAX,
					    int by = 1, double ly = -FLT_MAX, double uy = FLT_MAX,
					    int bz = 1, double lz = -FLT_MAX, double uz = FLT_MAX);
  /// Copy constructor
  MpdFemtoPicoEventCollectionVectorHideAway(const MpdFemtoPicoEventCollectionVectorHideAway &copy);
  /// Copy constructor
  MpdFemtoPicoEventCollectionVectorHideAway& operator=(const MpdFemtoPicoEventCollectionVectorHideAway& copy);
  /// Default destructor
  virtual ~MpdFemtoPicoEventCollectionVectorHideAway();

  /// Return pico event collection
  MpdFemtoPicoEventCollection* picoEventCollection(int, int, int);
  /// Return pico event collection
  MpdFemtoPicoEventCollection* picoEventCollection(double x, double y = 0, double z = 0);

  /// Return bin number on the x axis
  unsigned int binXNumber(double x) const {
    return (int) floor((x - mMinX) / mStepX);
  }
  /// Return bin number on the y axis
  unsigned int binYNumber(double y) const {
    return (int) floor((y - mMinY) / mStepY);
  }
  /// Return bin number on the z axis
  unsigned int binZNumber(double z) const {
    return (int) floor((z - mMinZ) / mStepZ);
  }

 private:

  /// Total number of bins
  int mBinsTot;
  /// Number of bins on x axis
  int mBinsX;
  /// Number of bins on y axis
  int mBinsY;
  /// Number of bins on z axis
  int mBinsZ;
  /// Minimal value on x axis
  double mMinX;
  /// Minimal value on y axis
  double mMinY;
  /// Minimal value on z axis
  double mMinZ;
  /// Maximal value on x axis
  double mMaxX;
  /// Maximal value on y axis
  double mMaxY;
  /// Maximal value on z axis
  double mMaxZ;
  /// Steps on x axis
  double mStepX;
  /// Steps on y axis
  double mStepY;
  /// Steps on z axis
  double mStepZ;
  /// Pico event collection
  MpdFemtoPicoEventCollection* mCollection;
  /// Collection vector
  MpdFemtoPicoEventCollectionVector mCollectionVector;

  ClassDef(MpdFemtoPicoEventCollectionVectorHideAway, 0)
};

#endif // #define MpdFemtoPicoEventCollectionVectorHideAway_h

/**
 * \class MpdFemtoBaseLikeSignCorrFctn
 * \brief A pure virtual base class for the like sign correlation function
 *
 * All like sign correlation functions  must inherit from this one
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * \email nigmatkulov@gmail.com
 */

#ifndef MpdFemtoBaseLikeSignCorrFctn_h
#define MpdFemtoBaseLikeSignCorrFctn_h

// Forward declaration
class MpdFemtoPair;

// MpdFemtoMaker headers
#include "MpdFemtoBaseCorrFctn.h"

// ROOT headers
#include "TList.h"

//_________________
class MpdFemtoBaseLikeSignCorrFctn : public MpdFemtoBaseCorrFctn {
 public:
  /// Default constructor

  MpdFemtoBaseLikeSignCorrFctn() {
    /* empty */
  }
  /// Copy constructor
  MpdFemtoBaseLikeSignCorrFctn(const MpdFemtoBaseLikeSignCorrFctn& copy);
  /// Assignment operator
  MpdFemtoBaseLikeSignCorrFctn& operator=(const MpdFemtoBaseLikeSignCorrFctn& copy);
  /// Default destructor

  virtual ~MpdFemtoBaseLikeSignCorrFctn() {
    /* empty */
  }

  /// Add positive pair
  virtual void addLikeSignPositivePair(const MpdFemtoPair*) = 0;
  /// Add negative pair
  virtual void addLikeSignNegativePair(const MpdFemtoPair*) = 0;

  /// Clone like-sign correlation function
  virtual MpdFemtoBaseLikeSignCorrFctn* clone() const = 0;

  /// Return output list
  virtual TList* getOutputList() = 0;

  /// The following allows "back-pointing" from the CorrFctn
  /// to the "parent" Analysis
  friend class MpdFemtoLikeSignAnalysis;
};

//_________________
inline MpdFemtoBaseLikeSignCorrFctn::MpdFemtoBaseLikeSignCorrFctn(const MpdFemtoBaseLikeSignCorrFctn& c) : MpdFemtoBaseCorrFctn(c) {
  mBaseAnalysis = nullptr;
}

//_________________
inline MpdFemtoBaseLikeSignCorrFctn& MpdFemtoBaseLikeSignCorrFctn::operator=(const MpdFemtoBaseLikeSignCorrFctn& c) {
  if (this != &c) {
    MpdFemtoBaseCorrFctn::operator=(c);
  }
  return *this;
}

#endif // #define MpdFemtoLikeSignCorrFctn_h

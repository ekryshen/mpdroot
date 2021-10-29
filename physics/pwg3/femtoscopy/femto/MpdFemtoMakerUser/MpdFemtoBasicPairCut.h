/**
 * \class MpdFemtoBasicPairCut
 * \brief A pair cut that keeps most of the pair observables
 *
 * MpdFemtoBasicPairCut keeps most of the pair observables and allows
 * to cut on them
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * \email nigmatkulov@gmail.com
 */

#ifndef MpdFemtoBasicPairCut_h
#define MpdFemtoBasicPairCut_h

// MpdFemtoMaker headers
#include "MpdFemtoBasePairCut.h"
#include "MpdFemtoTypes.h"
#include "SystemOfUnits.h"

// ROOT headers
#include "TVector3.h"

// C++ headers
#include <string>
#include <cstdio>
#include <iostream>

// Forward declaration
class MpdFemtoPair;

//_________________
class MpdFemtoBasicPairCut : public MpdFemtoBasePairCut {
 public:
  /// Default constructor
  MpdFemtoBasicPairCut();
  /// Copy constructor
  MpdFemtoBasicPairCut(const MpdFemtoBasicPairCut& copy);
  /// Assignment operator
  MpdFemtoBasicPairCut& operator=(const MpdFemtoBasicPairCut& copy);
  /// Destructor
  virtual ~MpdFemtoBasicPairCut();

  /// List of settings
  virtual TList *listSettings();

  /// Prepare report
  virtual MpdFemtoString report();

  /// Check if pair passes the cut
  virtual bool pass(const MpdFemtoPair*);

  /// Set anti-splitting cut (min,max) distributed in [-0.5, 1.0]
  void setQuality(const float& lo, const float& hi) {
    mQuality[0] = lo;
    mQuality[1] = hi;
  }
  /// Set half of pair transverse momentum cut (min, max)
  void setKt(const float& lo, const float& hi) {
    mKt[0] = lo;
    mKt[1] = hi;
  }
  /// Set pair transverse momentum cut (min, max)
  void setPt(const float& lo, const float& hi) {
    mPt[0] = lo;
    mPt[1] = hi;
  }
  /// Set cut on opening angle (min, max)
  void setOpeningAngle(const float& lo, const float& hi) {
    mOpeningAngle[0] = lo;
    mOpeningAngle[1] = hi;
  }
  /// Set cut on the distance between two track an the most inner TPC row (min, max)
  void setEntranceSeparation(const float& lo, const float& hi) {
    mEntranceSeparation[0] = lo;
    mEntranceSeparation[1] = hi;
  }
  /// Set rapidity cut (min, max)
  void setRapidity(const float& lo, const float& hi) {
    mRapidity[0] = lo;
    mRapidity[1] = hi;
  }
  /// Set pseudorapidity cut (min, max)
  void setEta(const float& lo, const float& hi) {
    mEta[0] = lo;
    mEta[1] = hi;
  }
  /// Set cut on \Delta \Eta (min, max)
  void setDEta(const float& lo, const float& hi) {
      mDEta[0] = lo;
      mDEta[1] = hi;  
  }
  /// Set cut on the invariant relative momentum (min, max)
  void setQinv(const float& lo, const float& hi) {
    mQinv[0] = lo;
    mQinv[1] = hi;
  }
  /// Set cut on the invariant mass (min, max)
  void setMinv(const float& lo, const float& hi) {
    mMinv[0] = lo;
    mMinv[1] = hi;
  }
  /// Set angle to primary vertex (min, max)
  void setAngleToPrimaryVertex(const float& lo, const float& hi) {
    mAngleToPrimaryVertex[0] = lo;
    mAngleToPrimaryVertex[1] = hi;
  }
  /// Set fraction of merged rows (min, max) distributed in [0., 1.]
  /// and -1. when tracks go through different sectors
  void setFracOfMergedRow(const float& lo, const float& hi) {
    mFracOfMergedRow[0] = lo;
    mFracOfMergedRow[1] = hi;
  }
  /// Set cut on closest row at DCA (min, max)
  void setClosestRowAtDCA(const float& lo, const float& hi) {
    mClosestRowAtDCA[0] = lo;
    mClosestRowAtDCA[1] = hi;
  }
  /// Set cut on weighted average separation between tracks in TPC (min, max)
  void setWeightedAvSep(const float& lo, const float& hi) {
    mWeightedAvSep[0] = lo;
    mWeightedAvSep[1] = hi;
  }
  /// Set cut on average separation between tracks in TPC (min, max)
  void setAverageSeparation(const float& lo, const float& hi) {
    mAverageSeparation[0] = lo;
    mAverageSeparation[1] = hi;
  }
  /// Set cut on R=\sqrt( dEta^2 + dPhi^2) (low)
  void setRValue(const float& lo) {
    mRValueLo = lo;
  }
  /// Set cut on delta phi* min
  void setDPhiStarMin(const float& lo, const float& hi) {
    mDPhiStarMin[0] = lo;
    mDPhiStarMin[1] = hi;
  }
  /// Set verbose mode
  /// \param true print information for each pair
  /// \param false do not print (default)
  void setVerboseMode(const bool& isVerbose) {
    mVerbose = isVerbose;
  }

 private:

  /// Primary vertex position
  TVector3 mPrimaryVertex;

  /// Quality - splitting level [-0.5, 1]
  float mQuality[2];
  /// Half of pair transverse momentum
  float mKt[2];
  /// Pair transverse momentum
  float mPt[2];
  /// Opening angle of a pair
  float mOpeningAngle[2];
  /// Rapidity of pair
  float mRapidity[2];
  /// Pseudorapidity of a pair
  float mEta[2];
  /// Invariant relative momentum
  float mQinv[2];
  /// Invariant mass
  float mMinv[2];
  /// Distance between two tracks at entrance to the TPC (inner radius)
  float mEntranceSeparation[2];
  /// Angle to primary vertex
  float mAngleToPrimaryVertex[2];
  /// Fraction of merged rows [0.,1.] and -1. corresponds to the case
  /// when tracks go to the different sectors
  float mFracOfMergedRow[2];
  /// Closest row at DCA
  float mClosestRowAtDCA[2];
  /// Weighted average spatial separation between two tracks estimated over several
  /// points in TPC (see MpdFemtoMaker/MpdFemtoPair for more details)
  float mWeightedAvSep[2];
  /// Average spatial separation between two tracks estimated over several
  /// points in TPC (see MpdFemtoMaker/MpdFemtoPair for more details)
  float mAverageSeparation[2];
  /// R = \sqrt( dEta^2 + dPhi^2)
  float mRValueLo;
  /// Minimal azimuthal angle between two tracks estimated over several radial points in TPC
  float mDPhiStarMin[2];
  /// Minimal and maximal value of \Delta \Eta
  float mDEta[2];
  /// Number of pairs that passed cut
  long mNPairsPassed;
  /// Number of pairs that failed cut
  long mNPairsFailed;
  /// Print information for each pair (default is false)
  bool mVerbose;

 protected:

  ClassDef(MpdFemtoBasicPairCut, 2)
};

#endif

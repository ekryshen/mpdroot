//
// Base class fo the femtoscopic weight generators
//

// MpdFemtoMaker headers
// Base
#include "MpdFemtoBaseModelWeightGenerator.h"
// Infrastructure
#include "MpdFemtoPair.h"
#include "MpdFemtoModelHiddenInfo.h"

ClassImp(MpdFemtoBaseModelWeightGenerator);

// Constants
const int MpdFemtoBaseModelWeightGenerator::fgkPairTypeNone = 0;
const int MpdFemtoBaseModelWeightGenerator::fgkPionPlusPionPlus = 1;
const int MpdFemtoBaseModelWeightGenerator::fgkPionPlusPionMinus = 2;
const int MpdFemtoBaseModelWeightGenerator::fgkKaonPlusKaonPlus = 3;
const int MpdFemtoBaseModelWeightGenerator::fgkKaonPlusKaonMinus = 4;
const int MpdFemtoBaseModelWeightGenerator::fgkProtonProton = 5;
const int MpdFemtoBaseModelWeightGenerator::fgkProtonAntiproton = 6;
const int MpdFemtoBaseModelWeightGenerator::fgkPionPlusKaonPlus = 7;
const int MpdFemtoBaseModelWeightGenerator::fgkPionPlusKaonMinus = 8;
const int MpdFemtoBaseModelWeightGenerator::fgkPionPlusProton = 9;
const int MpdFemtoBaseModelWeightGenerator::fgkPionPlusAntiproton = 10;
const int MpdFemtoBaseModelWeightGenerator::fgkKaonPlusProton = 11;
const int MpdFemtoBaseModelWeightGenerator::fgkKaonPlusAntiproton = 12;
const int MpdFemtoBaseModelWeightGenerator::fgkLambdaLambda = 13;
const int MpdFemtoBaseModelWeightGenerator::fgkAntilambdaAntilambda = 14;
const int MpdFemtoBaseModelWeightGenerator::fgkLambdaAntilambda = 15;

//_________________
MpdFemtoBaseModelWeightGenerator::MpdFemtoBaseModelWeightGenerator() :
mPairType(0),
mKStarOut(0), mKStarSide(0), mKStarLong(0), mKStar(0),
mRStarOut(0), mRStarSide(0), mRStarLong(0), mRStar(0) {
  /* empty */
}

//_________________
MpdFemtoBaseModelWeightGenerator::MpdFemtoBaseModelWeightGenerator(const MpdFemtoBaseModelWeightGenerator &aModel) :
mPairType(0),
mKStarOut(0), mKStarSide(0), mKStarLong(0), mKStar(0),
mRStarOut(0), mRStarSide(0), mRStarLong(0), mRStar(0) {
  // Copy constructor
  mPairType = aModel.mPairType;
}

//_________________
MpdFemtoBaseModelWeightGenerator& MpdFemtoBaseModelWeightGenerator::operator=(const MpdFemtoBaseModelWeightGenerator &aModel) {
  // Assignment operator
  if (this != &aModel) {
    mPairType = aModel.mPairType;
  }
  return *this;
}

//_____________________________________________
MpdFemtoBaseModelWeightGenerator::~MpdFemtoBaseModelWeightGenerator() {
  /* empty */
}

//_________________
int MpdFemtoBaseModelWeightGenerator::pairTypeFromPair(MpdFemtoPair *aPair) {
  // Get the type of pair from PID of particles in the pair
  MpdFemtoModelHiddenInfo *inf1 = (MpdFemtoModelHiddenInfo*) aPair->track1()->hiddenInfo();
  MpdFemtoModelHiddenInfo *inf2 = (MpdFemtoModelHiddenInfo*) aPair->track2()->hiddenInfo();

  int tPairType = fgkPairTypeNone;

  const int ktPid1 = inf1->pdgPid();
  const int ktPid2 = inf2->pdgPid();

  if (((ktPid1 == 211) && (ktPid2 == 211)) ||
      ((ktPid1 == -211) && (ktPid2 == -211))) {
    tPairType = fgkPionPlusPionPlus;
  } else if (((ktPid1 == -211) && (ktPid2 == 211)) ||
	     ((ktPid1 == 211) && (ktPid2 == -211))) {
    tPairType = fgkPionPlusPionMinus;
  } else if (((ktPid1 == 321) && (ktPid2 == 321)) ||
	     ((ktPid1 == -321) && (ktPid2 == -321))) {
    tPairType = fgkKaonPlusKaonPlus;
  } else if (((ktPid1 == -321) && (ktPid2 == 321)) ||
	     ((ktPid1 == 321) && (ktPid2 == -321))) {
    tPairType = fgkKaonPlusKaonMinus;
  } else if (((ktPid1 == 2212) && (ktPid2 == 2212)) ||
	     ((ktPid1 == -2212) && (ktPid2 == -2212))) {
    tPairType = fgkProtonProton;
  } else if (((ktPid1 == -2212) && (ktPid2 == 2212)) ||
	     ((ktPid1 == 2212) && (ktPid2 == -2212))) {
    tPairType = fgkProtonAntiproton;
  } else if (((ktPid1 == 211) && (ktPid2 == 321)) ||
	     ((ktPid1 == -211) && (ktPid2 == -321))) {
    tPairType = fgkPionPlusKaonPlus;
  } else if (((ktPid1 == -211) && (ktPid2 == 321)) ||
	     ((ktPid1 == 211) && (ktPid2 == -321))) {
    tPairType = fgkPionPlusKaonMinus;
  } else if (((ktPid1 == 211) && (ktPid2 == 2212)) ||
	     ((ktPid1 == -211) && (ktPid2 == -2212))) {
    tPairType = fgkPionPlusProton;
  } else if (((ktPid1 == -211) && (ktPid2 == 2212)) ||
	     ((ktPid1 == 211) && (ktPid2 == -2212))) {
    tPairType = fgkPionPlusAntiproton;
  } else if (((ktPid1 == 321) && (ktPid2 == 2212)) ||
	     ((ktPid1 == -321) && (ktPid2 == -2212))) {
    tPairType = fgkKaonPlusProton;
  } else if (((ktPid1 == -321) && (ktPid2 == 2212)) ||
	     ((ktPid1 == 321) && (ktPid2 == -2212))) {
    tPairType = fgkKaonPlusAntiproton;
  } else if (((ktPid1 == 3122) && (ktPid2 == 3122))) {
    tPairType = fgkLambdaLambda;
  } else if (((ktPid1 == -3122) && (ktPid2 == -3122))) {
    tPairType = fgkAntilambdaAntilambda;
  } else if (((ktPid1 == 3122) && (ktPid2 == -3122)) ||
	     ((ktPid1 == -3122) && (ktPid2 == 3122))) {
    tPairType = fgkLambdaAntilambda;
  }

  return tPairType;
}

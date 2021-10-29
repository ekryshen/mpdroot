//
// Manager for model studies
//

// MpdFemtoMaker headers
#include "MpdFemtoModelManager.h"
#include "MpdFemtoModelHiddenInfo.h"

// C++ headers
#include <iostream>

ClassImp(MpdFemtoModelManager);

//_____________________________________________
MpdFemtoModelManager::MpdFemtoModelManager() :
mFreezeOutGenerator(nullptr),
mWeightGenerator(nullptr),
mCreateCopyHiddenInfo(false) {
  // Default constructor
  /* empty */
}

//_________________
MpdFemtoModelManager::MpdFemtoModelManager(const MpdFemtoModelManager& copy) :
mFreezeOutGenerator(nullptr),
mWeightGenerator(nullptr),
mCreateCopyHiddenInfo(copy.mCreateCopyHiddenInfo) {
  // Copy constructor
  // Copy freeze-out coordinate generator if exists
  if (copy.mFreezeOutGenerator) {
    mFreezeOutGenerator = copy.mFreezeOutGenerator->clone();
  }

  // Copy femtoscopic weight generator if exists
  if (copy.mWeightGenerator) {
    mWeightGenerator = copy.mWeightGenerator->clone();
  }
}

//_________________
MpdFemtoModelManager::~MpdFemtoModelManager() {
  if (mFreezeOutGenerator) delete mFreezeOutGenerator;
  if (mWeightGenerator) delete mWeightGenerator;
}

//_________________
MpdFemtoModelManager& MpdFemtoModelManager::operator=(const MpdFemtoModelManager& copy) {
  // Check if the instance is not a copy of the object.
  // If not then copy it.
  if (this != &copy) {

    // Copy freeze-out coordinate generator if exists
    if (copy.mFreezeOutGenerator) {
      mFreezeOutGenerator = copy.mFreezeOutGenerator->clone();
    } else {
      mFreezeOutGenerator = nullptr;
    }

    // Copy femtoscopic weight generator if exists
    if (copy.mWeightGenerator) {
      mWeightGenerator = copy.mWeightGenerator->clone();
    } else {
      mWeightGenerator = nullptr;
    }

    // Copy hidden info copy status
    mCreateCopyHiddenInfo = copy.mCreateCopyHiddenInfo;
  }

  return *this;
}

//_____________________________________________
double MpdFemtoModelManager::weight(MpdFemtoPair *aPair) {
  if (!mWeightGenerator) {
    std::cout << "No weight generator set! Cannot calculate weight" << std::endl;
    return 1.0;
    // exit(0);
  }

  // Return femtoscopic weight for a given pair
  if (mCreateCopyHiddenInfo) {

    // Try to guess particle masses and pid from the weight generator
    int tPid1 = 0, tPid2 = 0;

    // pi+pi+
    if (mWeightGenerator->pairType() == MpdFemtoBaseModelWeightGenerator::pionPlusPionPlus()) {
      tPid1 = 211;
      tPid2 = 211;
    }

    // pi+pi-
    if (mWeightGenerator->pairType() == MpdFemtoBaseModelWeightGenerator::pionPlusPionMinus()) {
      tPid1 = 211;
      tPid2 = -211;
    }

    // K+K+
    if (mWeightGenerator->pairType() == MpdFemtoBaseModelWeightGenerator::kaonPlusKaonPlus()) {
      tPid1 = 321;
      tPid2 = 321;
    }

    // K+K-
    if (mWeightGenerator->pairType() == MpdFemtoBaseModelWeightGenerator::kaonPlusKaonMinus()) {
      tPid1 = 321;
      tPid2 = -321;
    }

    // pp
    if (mWeightGenerator->pairType() == MpdFemtoBaseModelWeightGenerator::protonProton()) {
      tPid1 = 2212;
      tPid2 = 2212;
    }

    // p antip
    if (mWeightGenerator->pairType() == MpdFemtoBaseModelWeightGenerator::protonAntiproton()) {
      tPid1 = 2212;
      tPid2 = -2212;
    }

    // pi+K+
    if (mWeightGenerator->pairType() == MpdFemtoBaseModelWeightGenerator::pionPlusKaonPlus()) {
      tPid1 = 211;
      tPid2 = 321;
    }

    // pi+K-
    if (mWeightGenerator->pairType() == MpdFemtoBaseModelWeightGenerator::pionPlusKaonMinus()) {
      tPid1 = 211;
      tPid2 = -321;
    }

    // pi+p
    if (mWeightGenerator->pairType() == MpdFemtoBaseModelWeightGenerator::pionPlusProton()) {
      tPid1 = 211;
      tPid2 = 2212;
    }

    // pi+antip
    if (mWeightGenerator->pairType() == MpdFemtoBaseModelWeightGenerator::pionPlusAntiproton()) {
      tPid1 = 211;
      tPid2 = -2212;
    }

    // K+p
    if (mWeightGenerator->pairType() == MpdFemtoBaseModelWeightGenerator::kaonPlusProton()) {
      tPid1 = 321;
      tPid2 = 2212;
    }

    // K+antip
    if (mWeightGenerator->pairType() == MpdFemtoBaseModelWeightGenerator::kaonPlusAntiproton()) {
      tPid1 = 321;
      tPid2 = -2212;
    }

    // LambdaLambda
    if (mWeightGenerator->pairType() == MpdFemtoBaseModelWeightGenerator::lambdaLambda()) {
      tPid1 = 3122;
      tPid2 = 3122;
    }

    // AntiLambdaAntiLambda
    if (mWeightGenerator->pairType() == MpdFemtoBaseModelWeightGenerator::antilambdaAntilambda()) {
      tPid1 = -3122;
      tPid2 = -3122;
    }

    // LambdaAntiLambda
    if (mWeightGenerator->pairType() == MpdFemtoBaseModelWeightGenerator::lambdaAntilambda()) {
      tPid1 = 3122;
      tPid2 = -3122;
    }

    // If hiddenInfo of the first particle from the pair is not set, then set it
    if (!(aPair->track1()->hiddenInfo())) {
      MpdFemtoModelHiddenInfo *inf1 = new MpdFemtoModelHiddenInfo();
      inf1->setTrueMomentum(aPair->track1()->track()->p());
      inf1->setPDGPid(tPid1);
      aPair->track1()->setHiddenInfo(inf1);
      delete inf1;
    }

    // If hiddenInfo of the second particle from the pair is not set, then set it
    if (!(aPair->track2()->hiddenInfo())) {
      MpdFemtoModelHiddenInfo *inf2 = new MpdFemtoModelHiddenInfo();
      inf2->setTrueMomentum(aPair->track2()->track()->p());
      inf2->setPDGPid(tPid2);
      aPair->track2()->setHiddenInfo(inf2);
      delete inf2;
    }
  }

  // If freeze-out coordinate generator exists, then perform generation
  if (mFreezeOutGenerator) {
    mFreezeOutGenerator->generateFreezeOut(aPair);
  }

  return mWeightGenerator->generateWeight(aPair);
}

//_________________
void MpdFemtoModelManager::createCopyHiddenInfo(bool aCopy) {
  mCreateCopyHiddenInfo = aCopy;
}

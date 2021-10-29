//
// The pure-virtual base class for correlation functions
//

#include "MpdFemtoBaseCorrFctn.h"

//_________________
MpdFemtoBaseCorrFctn::MpdFemtoBaseCorrFctn() :
mBaseAnalysis(nullptr),
mPairCut(nullptr) {
  /* no-op */
}

//_________________
MpdFemtoBaseCorrFctn::MpdFemtoBaseCorrFctn(const MpdFemtoBaseCorrFctn& c) :
mBaseAnalysis(c.mBaseAnalysis),
mPairCut(c.mPairCut) {
  /* no-op */
}

//_________________
MpdFemtoBaseCorrFctn& MpdFemtoBaseCorrFctn::operator=(const MpdFemtoBaseCorrFctn& c) {
  if (this != &c) {
    mBaseAnalysis = c.mBaseAnalysis;
    mPairCut = c.mPairCut;
  }
  return *this;
}

//________________
void MpdFemtoBaseCorrFctn::addRealPair(MpdFemtoPair* /* pair */) {
  std::cout << " MpdFemtoBaseCorrFctn::addRealPair - Not implemented" << std::endl;
}

//_________________
void MpdFemtoBaseCorrFctn::addMixedPair(MpdFemtoPair* /* pair */) {
  std::cout << "MpdFemtoBaseCorrFctn::addMixedPair - Not implemented" << std::endl;
}

//_________________
void MpdFemtoBaseCorrFctn::addRealTriplet(const MpdFemtoTriplet* /* triplet */) {
  std::cout << "MpdFemtoBaseCorrFctn::addRealTriplet - Not implemented" << std::endl;
}

//_________________
void MpdFemtoBaseCorrFctn::addMixedTriplet(const MpdFemtoTriplet* /* triplet */) {
  std::cout << "MpdFemtoBaseCorrFctn::addMixedTriplet - Not implemented" << std::endl;
}

//_________________
void MpdFemtoBaseCorrFctn::addFirstParticle(MpdFemtoParticle* /* part */, bool /* isMixing*/) {
  std::cout << "MpdFemtoBaseCorrFctn::addFirstParticle -- Not implemented" << std::endl;
}

//_________________
void MpdFemtoBaseCorrFctn::addSecondParticle(MpdFemtoParticle* /* part */) {
  std::cout << "MpdFemtoBaseCorrFctn::addSecondParticle -- Not implemented" << std::endl;
}

ClassImp(MpdFemtoBaseCorrFctn)

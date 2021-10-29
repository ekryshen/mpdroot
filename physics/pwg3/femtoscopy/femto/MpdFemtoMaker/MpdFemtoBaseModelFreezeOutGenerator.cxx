//
// A base class for freeze-out coordinate generator
//

// MpdFemtoMaker headers
#include "MpdFemtoBaseModelFreezeOutGenerator.h"

//_________________
MpdFemtoBaseModelFreezeOutGenerator::MpdFemtoBaseModelFreezeOutGenerator() : mRandom(0) {
  /* emtpy */
}

//_________________
MpdFemtoBaseModelFreezeOutGenerator::MpdFemtoBaseModelFreezeOutGenerator(const MpdFemtoBaseModelFreezeOutGenerator &/* aModel */) :
mRandom(0) {
  /* empty */
}

//_________________
MpdFemtoBaseModelFreezeOutGenerator& MpdFemtoBaseModelFreezeOutGenerator::operator=(const MpdFemtoBaseModelFreezeOutGenerator& aGen) {
  // Assignment operator
  if (this != &aGen) {
    if (aGen.mRandom) {
      mRandom = new TRandom3(*aGen.mRandom);
    } else {
      mRandom = 0;
    }
  }

  return *this;
}

//_________________
MpdFemtoBaseModelFreezeOutGenerator::~MpdFemtoBaseModelFreezeOutGenerator() {
  if (mRandom) delete mRandom;
}

//_________________
MpdFemtoBaseModelFreezeOutGenerator* MpdFemtoBaseModelFreezeOutGenerator::clone() const {
  return nullptr;
}

ClassImp(MpdFemtoBaseModelFreezeOutGenerator);

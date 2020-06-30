//
// Stores collection of particles for processing
//

// MpdFemtoMaker headers
#include "MpdFemtoPicoEvent.h"

//_________________
MpdFemtoPicoEvent::MpdFemtoPicoEvent() {
  mFirstParticleCollection = new MpdFemtoParticleCollection();
  mSecondParticleCollection = new MpdFemtoParticleCollection();
  mThirdParticleCollection = new MpdFemtoParticleCollection();
  // Default constructor
  /* empty */
}

//_________________
MpdFemtoPicoEvent::~MpdFemtoPicoEvent() {
  // Destructor
  MpdFemtoParticleIterator iter;

  if (mFirstParticleCollection) {
    for (iter = mFirstParticleCollection->begin(); iter != mFirstParticleCollection->end(); iter++) {
      delete *iter;
    }
    mFirstParticleCollection->clear();
    delete mFirstParticleCollection;
    mFirstParticleCollection = nullptr;
  }

  if (mSecondParticleCollection) {
    for (iter = mSecondParticleCollection->begin(); iter != mSecondParticleCollection->end(); iter++) {
      delete *iter;
    }
    mSecondParticleCollection->clear();
    delete mSecondParticleCollection;
    mSecondParticleCollection = nullptr;
  }

  if (mThirdParticleCollection) {
    if (mThirdParticleCollection->size() != 0) {
      for (iter = mThirdParticleCollection->begin(); iter != mThirdParticleCollection->end(); iter++) {
	delete *iter;
      }
    }
    mThirdParticleCollection->clear();
    delete mThirdParticleCollection;
    mThirdParticleCollection = nullptr;
  }
}

//_________________
MpdFemtoPicoEvent::MpdFemtoPicoEvent(const MpdFemtoPicoEvent& pico) :
mFirstParticleCollection(nullptr),
mSecondParticleCollection(nullptr),
mThirdParticleCollection(nullptr) {
  // Copy constructor
  MpdFemtoParticleIterator iter;

  mFirstParticleCollection = new MpdFemtoParticleCollection;
  if (pico.mFirstParticleCollection) {
    for (iter = pico.mFirstParticleCollection->begin();
	 iter != pico.mFirstParticleCollection->end(); iter++) {
      mFirstParticleCollection->push_back(*iter);
    }
  }

  mSecondParticleCollection = new MpdFemtoParticleCollection;
  if (pico.mSecondParticleCollection) {
    for (iter = pico.mSecondParticleCollection->begin();
	 iter != pico.mSecondParticleCollection->end(); iter++) {
      mSecondParticleCollection->push_back(*iter);
    }
  }

  mThirdParticleCollection = new MpdFemtoParticleCollection;
  if (pico.mThirdParticleCollection) {
    for (iter = pico.mThirdParticleCollection->begin();
	 iter != pico.mThirdParticleCollection->end(); iter++) {
      mThirdParticleCollection->push_back(*iter);
    }
  }
}

//_________________
MpdFemtoPicoEvent& MpdFemtoPicoEvent::operator=(const MpdFemtoPicoEvent& pico) {
  // Assignment operator
  if (this != &pico) {

    MpdFemtoParticleIterator iter;

    // Clean collections
    if (mFirstParticleCollection) {
      for (iter = mFirstParticleCollection->begin();
	   iter != mFirstParticleCollection->end(); iter++) {
	delete *iter;
      }
      mFirstParticleCollection->clear();
      delete mFirstParticleCollection;
      mFirstParticleCollection = nullptr;
    } //if (mFirstParticleCollection)

    if (mSecondParticleCollection) {
      for (iter = mSecondParticleCollection->begin();
	   iter != mSecondParticleCollection->end(); iter++) {
	delete *iter;
      }
      mSecondParticleCollection->clear();
      delete mSecondParticleCollection;
      mSecondParticleCollection = nullptr;
    } //if (mSecondParticleCollection)

    if (mThirdParticleCollection) {
      for (iter = mThirdParticleCollection->begin();
	   iter != mThirdParticleCollection->end(); iter++) {
	delete *iter;
      }
      mThirdParticleCollection->clear();
      delete mThirdParticleCollection;
      mThirdParticleCollection = nullptr;
    } //if (mThirdParticleCollection)

    // Copy collections
    mFirstParticleCollection = new MpdFemtoParticleCollection;
    if (pico.mFirstParticleCollection) {
      for (iter = pico.mFirstParticleCollection->begin();
	   iter != pico.mFirstParticleCollection->end(); iter++) {
	mFirstParticleCollection->push_back(*iter);
      }
    } //if (pico.mFirstParticleCollection)

    mSecondParticleCollection = new MpdFemtoParticleCollection;
    if (pico.mSecondParticleCollection) {
      for (iter = pico.mSecondParticleCollection->begin();
	   iter != pico.mSecondParticleCollection->end(); iter++) {
	mSecondParticleCollection->push_back(*iter);
      }
    } //if (pico.mSecondParticleCollection)

    mThirdParticleCollection = new MpdFemtoParticleCollection;
    if (pico.mThirdParticleCollection) {
      for (iter = pico.mThirdParticleCollection->begin();
	   iter != pico.mThirdParticleCollection->end(); iter++) {
	mThirdParticleCollection->push_back(*iter);
      }
    } //if (pico.mThirdParticleCollection)

  } //if( this != &pico)

  return *this;
}

ClassImp(MpdFemtoPicoEvent)

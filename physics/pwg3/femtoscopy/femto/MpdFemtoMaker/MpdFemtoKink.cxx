//
// Main class holding kink information
//

// MpdFemtoMaker headers
#include "MpdFemtoKink.h"
#include "MpdFemtoTrack.h"

ClassImp(MpdFemtoKink)

//_________________
MpdFemtoKink::MpdFemtoKink() :
mDcaParentDaughter(999),
mHitDistanceParentDaughter(999),
mHitDistanceParentVertex(999),
mDeltaEnergy{},
mDecayAngle(0),
mDecayAngleCM(0),
mDaughter(),
mParent(),
mPositionX(999),
mPositionY(999),
mPositionZ(999) {
  // Default constructor
  /* empty */
}

//_________________
MpdFemtoKink::MpdFemtoKink(const MpdFemtoKink& k) :
mDcaParentDaughter(k.mDcaParentDaughter),
mHitDistanceParentDaughter(k.mHitDistanceParentDaughter),
mHitDistanceParentVertex(k.mHitDistanceParentVertex),
mDecayAngle(k.mDecayAngle),
mDecayAngleCM(k.mDecayAngleCM),
mDaughter(k.mDaughter),
mParent(k.mParent),
mPositionX(k.mPositionX),
mPositionY(k.mPositionY),
mPositionZ(k.mPositionZ) {
  // Copy constructor
  memcpy(mDeltaEnergy, k.mDeltaEnergy, sizeof (mDeltaEnergy));
}

//_________________
MpdFemtoKink& MpdFemtoKink::operator=(const MpdFemtoKink &k) {
  if (this != &k) {
    mDcaParentDaughter = k.mDcaParentDaughter;
    mHitDistanceParentDaughter = k.mHitDistanceParentDaughter;
    mHitDistanceParentVertex = k.mHitDistanceParentVertex;
    mDeltaEnergy[0] = k.mDeltaEnergy[0];
    mDeltaEnergy[1] = k.mDeltaEnergy[1];
    mDeltaEnergy[2] = k.mDeltaEnergy[2];
    mDecayAngle = k.mDecayAngle;
    mDecayAngleCM = k.mDecayAngleCM;
    mDaughter = k.mDaughter;
    mParent = k.mParent;
    mPositionX = k.mPositionX;
    mPositionY = k.mPositionY;
    mPositionZ = k.mPositionZ;
  }

  return *this;
}

//_________________
MpdFemtoKink::~MpdFemtoKink() {
  // Destructor
  /* empty */
}

//_________________
void MpdFemtoKink::setDaugherHbtTrack(const MpdFemtoTrack& trk) {
  mDaughter = trk;
}

//_________________
void MpdFemtoKink::setParetntHbtTrack(const MpdFemtoTrack& trk) {
  mParent = trk;
}

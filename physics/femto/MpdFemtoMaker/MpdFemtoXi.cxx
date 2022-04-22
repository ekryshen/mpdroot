//
// Main class holding xi information
//

// C++ headers
#include <limits>

// MpdFemtoMaker headers
#include "MpdFemtoXi.h"

ClassImp(MpdFemtoXi);

//_________________
MpdFemtoXi::MpdFemtoXi() :
MpdFemtoV0(),
mCharge(false),
mDecayVertexXiX(0), mDecayVertexXiY(0), mDecayVertexXiZ(0),
mDcaXiDaughters(999), mDcaXiToPrimVertex(999),
mMomXiX(0), mMomXiY(0), mMomXiZ(0),
mChi2Xi(std::numeric_limits<unsigned char>::max()),
mClXi(0),
mBachelor(nullptr),
mMomBacAtDca2DecayPointX(0), mMomBacAtDca2DecayPointY(0), mMomBacAtDca2DecayPointZ(0) {
  // Default constructor
  if (!mBachelor) {
    mBachelor = new MpdFemtoTrack();
  }
}

//_________________
MpdFemtoXi::MpdFemtoXi(const MpdFemtoXi& xi) :
MpdFemtoV0(xi),
mCharge(xi.mCharge),
mDecayVertexXiX(xi.mDecayVertexXiX),
mDecayVertexXiY(xi.mDecayVertexXiY),
mDecayVertexXiZ(xi.mDecayVertexXiZ),
mDcaXiDaughters(xi.mDcaXiDaughters),
mDcaXiToPrimVertex(xi.mDcaXiToPrimVertex),
mMomXiX(xi.mMomXiX), mMomXiY(xi.mMomXiY), mMomXiZ(xi.mMomXiZ),
mChi2Xi(xi.mChi2Xi),
mBachelor(nullptr),
mMomBacAtDca2DecayPointX(xi.mMomBacAtDca2DecayPointX),
mMomBacAtDca2DecayPointY(xi.mMomBacAtDca2DecayPointY),
mMomBacAtDca2DecayPointZ(xi.mMomBacAtDca2DecayPointZ) {
  // Copy bachelor info
  if (!mBachelor) {
    mBachelor = new MpdFemtoTrack(*xi.mBachelor);
  } else {
    mBachelor = xi.mBachelor;
  }
}

//_________________
MpdFemtoXi& MpdFemtoXi::operator=(const MpdFemtoXi& xi) {
  // Assignment operator
  if (this != &xi) {
    mCharge = xi.mCharge;
    mDecayVertexXiX = xi.mDecayVertexXiX;
    mDecayVertexXiY = xi.mDecayVertexXiY;
    mDecayVertexXiZ = xi.mDecayVertexXiZ;
    mDcaXiDaughters = xi.mDcaXiDaughters;
    mDcaXiToPrimVertex = xi.mDcaXiToPrimVertex;
    mMomXiX = xi.mMomXiX;
    mMomXiY = xi.mMomXiY;
    mMomXiZ = xi.mMomXiZ;
    mChi2Xi = xi.mChi2Xi;
    if (!mBachelor) {
      mBachelor = new MpdFemtoTrack(*xi.mBachelor);
    } else {
      mBachelor = xi.mBachelor;
    }
    mMomBacAtDca2DecayPointX = xi.mMomBacAtDca2DecayPointX;
    mMomBacAtDca2DecayPointY = xi.mMomBacAtDca2DecayPointY;
    mMomBacAtDca2DecayPointZ = xi.mMomBacAtDca2DecayPointZ;
  }

  return *this;
}

//_________________
MpdFemtoXi::~MpdFemtoXi() {
  if (mBachelor) delete mBachelor;
}

//_________________
void MpdFemtoXi::updateXi() {
  // TODO:
  // Momentum and the logic of Xi decay point MUST be written here!!!!
  // Not yet done due to the lack of time.
  // For example, mMomBacAtDca2DecayPointX(Y,Z) should be calculated here
}

//_________________
float MpdFemtoXi::alphaArmXi() const {
  float MomV0AlongXi = momV0() * momXi() / ptotXi();
  float MomBacAlongXi = momBac() * momXi() / ptotXi();
  return ( (MomBacAlongXi - MomV0AlongXi) / (MomBacAlongXi + MomV0AlongXi));
}

//_________________
float MpdFemtoXi::ptArmXi() const {
  float MomBacAlongXi = momBac() * momXi() / ptotXi();
  return TMath::Sqrt(ptot2Bac() - MomBacAlongXi * MomBacAlongXi);
}

//_________________
void MpdFemtoXi::setChi2Xi(const float& chi2) {
  if (chi2 < 0) {
    mChi2Xi = 0;
  } else {
    mChi2Xi = ((chi2 * 10) > std::numeric_limits<unsigned char>::max() ?
	       std::numeric_limits<unsigned char>::max() :
	       (unsigned char) (chi2 * 10));
  }
}

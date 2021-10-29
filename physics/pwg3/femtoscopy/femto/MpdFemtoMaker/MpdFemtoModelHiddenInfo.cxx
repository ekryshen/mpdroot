//
// The hidden information for model calculations
//

// C++ headers
#include <iostream>

// MpdFemtoMaker headers
#include "MpdFemtoModelHiddenInfo.h"

//_________________
MpdFemtoModelHiddenInfo::MpdFemtoModelHiddenInfo() :
mTrueMomentumX(nullptr), mTrueMomentumY(nullptr), mTrueMomentumZ(nullptr),
mMotherMomentumX(nullptr), mMotherMomentumY(nullptr), mMotherMomentumZ(nullptr),
mEmissionPointX(nullptr), mEmissionPointY(nullptr), mEmissionPointZ(nullptr), mEmissionPointT(nullptr),
mPDGPid(0), mMotherPdg(0),
mTrueMomentumPosX(nullptr), mTrueMomentumPosY(nullptr), mTrueMomentumPosZ(nullptr),
mEmissionPointPosX(nullptr), mEmissionPointPosY(nullptr), mEmissionPointPosZ(nullptr),
mEmissionPointPosT(nullptr),
mPDGPidPos(0),
mTrueMomentumNegX(nullptr), mTrueMomentumNegY(nullptr), mTrueMomentumNegZ(nullptr),
mEmissionPointNegX(nullptr), mEmissionPointNegY(nullptr), mEmissionPointNegZ(nullptr),
mEmissionPointNegT(nullptr),
mPDGPidNeg(0),
mOrigin(0) {
  // Default constructor
}

//_________________
MpdFemtoModelHiddenInfo::MpdFemtoModelHiddenInfo(const MpdFemtoModelHiddenInfo &aInfo) :
MpdFemtoHiddenInfo(aInfo),
mTrueMomentumX(nullptr), mTrueMomentumY(nullptr), mTrueMomentumZ(nullptr),
mMotherMomentumX(nullptr), mMotherMomentumY(nullptr), mMotherMomentumZ(nullptr),
mEmissionPointX(nullptr), mEmissionPointY(nullptr), mEmissionPointZ(nullptr), mEmissionPointT(nullptr),
mPDGPid(aInfo.mPDGPid), mMotherPdg(aInfo.mMotherPdg),
mTrueMomentumPosX(nullptr), mTrueMomentumPosY(nullptr), mTrueMomentumPosZ(nullptr),
mEmissionPointPosX(nullptr), mEmissionPointPosY(nullptr), mEmissionPointPosZ(nullptr),
mEmissionPointPosT(nullptr),
mPDGPidPos(aInfo.mPDGPidPos),
mTrueMomentumNegX(nullptr), mTrueMomentumNegY(nullptr), mTrueMomentumNegZ(nullptr),
mEmissionPointNegX(nullptr), mEmissionPointNegY(nullptr), mEmissionPointNegZ(nullptr),
mEmissionPointNegT(nullptr),
mPDGPidNeg(aInfo.mPDGPidNeg),
mOrigin(aInfo.mOrigin) {
  // Copy constructor
  setTrueMomentum(aInfo.trueMomentum());
  setMotherMomentum(aInfo.motherMomentum());
  setEmissionPoint(aInfo.emissionPoint());
  setTrueMomentumPos(aInfo.trueMomentumPos());
  setEmissionPointPos(aInfo.emissionPointPos());
  setTrueMomentumNeg(aInfo.trueMomentumNeg());
  setEmissionPointNeg(aInfo.emissionPointNeg());
  setOrigin(aInfo.origin());
}

//_________________
MpdFemtoModelHiddenInfo& MpdFemtoModelHiddenInfo::operator=(const MpdFemtoModelHiddenInfo& aInfo) {
  // Assignment operator
  if (this != &aInfo) {

    setTrueMomentum(aInfo.trueMomentum());
    setMotherMomentum(aInfo.motherMomentum());
    setEmissionPoint(aInfo.emissionPoint());
    mPDGPid = aInfo.pdgPid();
    mMotherPdg = aInfo.motherPdgCode();

    setTrueMomentumPos(aInfo.trueMomentumPos());
    setEmissionPointPos(aInfo.emissionPointPos());
    mPDGPidPos = aInfo.pdgPidPos();

    setTrueMomentumNeg(aInfo.trueMomentumNeg());
    setEmissionPointNeg(aInfo.emissionPointNeg());
    mPDGPidNeg = aInfo.pdgPidNeg();

    mOrigin = aInfo.origin();
  }

  return *this;
}

//_________________
MpdFemtoModelHiddenInfo::~MpdFemtoModelHiddenInfo() {
  // Destructor
  if (mTrueMomentumX) delete mTrueMomentumX;
  if (mTrueMomentumY) delete mTrueMomentumY;
  if (mTrueMomentumZ) delete mTrueMomentumZ;
  if (mMotherMomentumX) delete mMotherMomentumX;
  if (mMotherMomentumY) delete mMotherMomentumY;
  if (mMotherMomentumZ) delete mMotherMomentumZ;
  if (mEmissionPointX) delete mEmissionPointX;
  if (mEmissionPointY) delete mEmissionPointY;
  if (mEmissionPointZ) delete mEmissionPointZ;
  if (mEmissionPointT) delete mEmissionPointT;

  if (mTrueMomentumPosX) delete mTrueMomentumPosX;
  if (mTrueMomentumPosY) delete mTrueMomentumPosY;
  if (mTrueMomentumPosZ) delete mTrueMomentumPosZ;
  if (mEmissionPointPosX) delete mEmissionPointPosX;
  if (mEmissionPointPosY) delete mEmissionPointPosY;
  if (mEmissionPointPosZ) delete mEmissionPointPosZ;
  if (mEmissionPointPosT) delete mEmissionPointPosT;

  if (mTrueMomentumNegX) delete mTrueMomentumNegX;
  if (mTrueMomentumNegY) delete mTrueMomentumNegY;
  if (mTrueMomentumNegZ) delete mTrueMomentumNegZ;
  if (mEmissionPointNegX) delete mEmissionPointNegX;
  if (mEmissionPointNegY) delete mEmissionPointNegY;
  if (mEmissionPointNegZ) delete mEmissionPointNegZ;
  if (mEmissionPointNegT) delete mEmissionPointNegT;
}

//_________________
MpdFemtoHiddenInfo* MpdFemtoModelHiddenInfo::getParticleHiddenInfo() const {
  MpdFemtoModelHiddenInfo* mBuf = new MpdFemtoModelHiddenInfo(*this);
  return mBuf;
}

//_____________________________________________
TVector3 MpdFemtoModelHiddenInfo::trueMomentum() const {
  TVector3 mom;
  if (!mTrueMomentumX || !mTrueMomentumY || !mTrueMomentumZ) {
    if (!mTrueMomentumX) {
      // std::cout << "mTrueMomentumX does not exist" << std::endl;
    }
    if (!mTrueMomentumY) {
      // std::cout << "mTrueMomentumY does not exist" << std::endl;
    }
    if (!mTrueMomentumZ) {
      // std::cout << "mTrueMomentumZ does not exist" << std::endl;
    }
    // std::cout << "Will return dummy (0,0,0) momentum" << std::endl;

    mom.SetXYZ(0, 0, 0);
  } else {
    mom.SetXYZ(*mTrueMomentumX, *mTrueMomentumY, *mTrueMomentumZ);
  }
  return mom;
}

//_________________
TVector3 MpdFemtoModelHiddenInfo::motherMomentum() const {
  TVector3 mom;
  if (!mMotherMomentumX || !mMotherMomentumY || !mMotherMomentumZ) {
    if (!mMotherMomentumX) {
      // std::cout << "mMotherMomentumX does not exist" << std::endl;
    }
    if (!mMotherMomentumY) {
      // std::cout << "mMotherMomentumY does not exist" << std::endl;
    }
    if (!mMotherMomentumZ) {
      // std::cout << "mMotherMomentumZ does not exist" << std::endl;
    }
    mom.SetXYZ(0., 0., 0.);
  } else {
    mom.SetXYZ(*mMotherMomentumX, *mMotherMomentumY, *mMotherMomentumZ);
  }
  return mom;
}

//_________________
TLorentzVector MpdFemtoModelHiddenInfo::emissionPoint() const {
  TLorentzVector point;
  if (!mEmissionPointX || !mEmissionPointY || !mEmissionPointZ || !mEmissionPointT) {
    if (!mEmissionPointX) {
      // std::cout << "mEmissionPointX does not exist" << std::endl;
    }
    if (!mEmissionPointY) {
      // std::cout << "mEmissionPointY does not exist" << std::endl;
    }
    if (!mEmissionPointZ) {
      // std::cout << "mEmissionPointZ does not exist" << std::endl;
    }
    if (!mEmissionPointT) {
      // std::cout << "mEmissionPointT does nto exist" << std::endl;
    }
    point.SetXYZT(0, 0, 0, 0);
  } else {
    point.SetXYZT(*mEmissionPointX, *mEmissionPointY, *mEmissionPointZ, *mEmissionPointT);
  }
  return point;
}

//_________________
void MpdFemtoModelHiddenInfo::setTrueMomentum(TVector3 *aMom) {
  if (!mTrueMomentumX) mTrueMomentumX = new float;
  if (!mTrueMomentumY) mTrueMomentumY = new float;
  if (!mTrueMomentumZ) mTrueMomentumZ = new float;

  if (aMom) {
    *mTrueMomentumX = aMom->X();
    *mTrueMomentumY = aMom->Y();
    *mTrueMomentumZ = aMom->Z();
  } else {
    *mTrueMomentumX = 0;
    *mTrueMomentumY = 0;
    *mTrueMomentumZ = 0;
  }
}

//_________________
void MpdFemtoModelHiddenInfo::setTrueMomentum(const TVector3& aMom) {
  // Set momentum from vector
  if (!mTrueMomentumX) mTrueMomentumX = new float;
  if (!mTrueMomentumY) mTrueMomentumY = new float;
  if (!mTrueMomentumZ) mTrueMomentumZ = new float;

  *mTrueMomentumX = aMom.X();
  *mTrueMomentumY = aMom.Y();
  *mTrueMomentumZ = aMom.Z();
}

//_________________
void MpdFemtoModelHiddenInfo::setTrueMomentum(const double& aPx, const double& aPy,
					      const double& aPz) {
  // Set momentum from components
  if (mTrueMomentumX == nullptr) mTrueMomentumX = new float;
  if (mTrueMomentumY == nullptr) mTrueMomentumY = new float;
  if (mTrueMomentumZ == nullptr) mTrueMomentumZ = new float;
  *mTrueMomentumX = aPx;
  *mTrueMomentumY = aPy;
  *mTrueMomentumZ = aPz;
}

//_________________
void MpdFemtoModelHiddenInfo::setMotherMomentum(TVector3 *aMom) {
  if (mMotherMomentumX == nullptr) mMotherMomentumX = new float;
  if (mMotherMomentumY == nullptr) mMotherMomentumY = new float;
  if (mMotherMomentumZ == nullptr) mMotherMomentumZ = new float;

  if (aMom) {
    *mMotherMomentumX = aMom->X();
    *mMotherMomentumY = aMom->Y();
    *mMotherMomentumZ = aMom->Z();
  } else {
    *mMotherMomentumX = 0;
    *mMotherMomentumY = 0;
    *mMotherMomentumZ = 0;
  }
}

//_________________
void MpdFemtoModelHiddenInfo::setMotherMomentum(const TVector3& aMom) {
  // Set momentum from vector
  if (mMotherMomentumX == nullptr) mMotherMomentumX = new float;
  if (mMotherMomentumY == nullptr) mMotherMomentumY = new float;
  if (mMotherMomentumZ == nullptr) mMotherMomentumZ = new float;
  *mMotherMomentumX = aMom.X();
  *mMotherMomentumY = aMom.Y();
  *mMotherMomentumZ = aMom.Z();
}

//_________________
void MpdFemtoModelHiddenInfo::setMotherMomentum(const double& aPx, const double& aPy,
						const double& aPz) {
  // Set momentum from components
  if (mMotherMomentumX == nullptr) mMotherMomentumX = new float;
  if (mMotherMomentumY == nullptr) mMotherMomentumY = new float;
  if (mMotherMomentumZ == nullptr) mMotherMomentumZ = new float;
  *mMotherMomentumX = aPx;
  *mMotherMomentumY = aPy;
  *mMotherMomentumZ = aPz;
}

//_________________
void MpdFemtoModelHiddenInfo::setEmissionPoint(TLorentzVector *aPos) {
  // Set emission point pararmeters
  if (mEmissionPointX == nullptr) mEmissionPointX = new float;
  if (mEmissionPointY == nullptr) mEmissionPointY = new float;
  if (mEmissionPointZ == nullptr) mEmissionPointZ = new float;
  if (mEmissionPointT == nullptr) mEmissionPointT = new float;
  *mEmissionPointX = aPos->X();
  *mEmissionPointY = aPos->Y();
  *mEmissionPointZ = aPos->Z();
  *mEmissionPointT = aPos->T();
}

//_________________
void MpdFemtoModelHiddenInfo::setEmissionPoint(const TLorentzVector& aPos) {
  // Set position from vector
  if (mEmissionPointX == nullptr) mEmissionPointX = new float;
  if (mEmissionPointY == nullptr) mEmissionPointY = new float;
  if (mEmissionPointZ == nullptr) mEmissionPointZ = new float;
  if (mEmissionPointT == nullptr) mEmissionPointT = new float;
  *mEmissionPointX = aPos.X();
  *mEmissionPointY = aPos.Y();
  *mEmissionPointZ = aPos.Z();
  *mEmissionPointT = aPos.T();
}

//_________________
void MpdFemtoModelHiddenInfo::setEmissionPoint(const double& aRx, const double& aRy,
					       const double& aRz, const double& aT) {
  // Set position from components
  if (mEmissionPointX == nullptr) mEmissionPointX = new float;
  if (mEmissionPointY == nullptr) mEmissionPointY = new float;
  if (mEmissionPointZ == nullptr) mEmissionPointZ = new float;
  if (mEmissionPointT == nullptr) mEmissionPointT = new float;
  *mEmissionPointX = aRx;
  *mEmissionPointY = aRy;
  *mEmissionPointZ = aRz;
  *mEmissionPointT = aT;
}

//_________________
TVector3 MpdFemtoModelHiddenInfo::trueMomentumPos() const {
  TVector3 mom;
  if (!mTrueMomentumPosX || !mTrueMomentumPosY || !mTrueMomentumPosZ) {
    if (!mTrueMomentumPosX) {
      // std::cout << "mTrueMomentumPosX does not exist" << std::endl;
    }
    if (!mTrueMomentumPosY) {
      // std::cout << "mTrueMomentumPosY does not exist" << std::endl;
    }
    if (!mTrueMomentumPosZ) {
      // std::cout << "mTrueMomentumPosZ does not exist" << std::endl;
    }
    mom.SetXYZ(0, 0, 0);
  } else {
    mom.SetXYZ(*mTrueMomentumPosX, *mTrueMomentumPosY, *mTrueMomentumPosZ);
  }
  return mom;
}

//_________________
TLorentzVector MpdFemtoModelHiddenInfo::emissionPointPos() const {
  TLorentzVector point;
  if (!mEmissionPointPosX || !mEmissionPointPosY || !mEmissionPointPosZ || !mEmissionPointPosT) {
    if (!mEmissionPointPosX) {
      // std::cout << "mEmissionPointPosX does not exist" << std::endl;
    }
    if (!mEmissionPointPosY) {
      // std::cout << "mEmissionPointPosY does not exist" << std::endl;
    }
    if (!mEmissionPointPosZ) {
      // std::cout << "mEmissionPointPosZ does not exist" << std::endl;
    }
    if (!mEmissionPointPosT) {
      // std::cout << "mEmissionPointPosT does not exist" << std::endl;
    }
    point.SetXYZT(0, 0, 0, 0);
  } else {
    point.SetXYZT(*mEmissionPointPosX, *mEmissionPointPosY, *mEmissionPointPosY, *mEmissionPointPosT);
  }
  return point;
}

//_________________
void MpdFemtoModelHiddenInfo::setTrueMomentumPos(TVector3 *aMom) {
  // Set momentum of the positive daughter
  if (!mTrueMomentumPosX) mTrueMomentumPosX = new float;
  if (!mTrueMomentumPosY) mTrueMomentumPosY = new float;
  if (!mTrueMomentumPosZ) mTrueMomentumPosZ = new float;

  if (aMom) {
    *mTrueMomentumPosX = aMom->X();
    *mTrueMomentumPosY = aMom->Y();
    *mTrueMomentumPosZ = aMom->Z();
  } else {
    delete mTrueMomentumPosX;
    mTrueMomentumPosX = nullptr;
    delete mTrueMomentumPosY;
    mTrueMomentumPosY = nullptr;
    delete mTrueMomentumPosZ;
    mTrueMomentumPosZ = nullptr;
  }
}

//_________________
void MpdFemtoModelHiddenInfo::setTrueMomentumPos(const TVector3& aMom) {
  // Set momentum from vector
  if (!mTrueMomentumPosX) mTrueMomentumPosX = new float;
  if (!mTrueMomentumPosY) mTrueMomentumPosY = new float;
  if (!mTrueMomentumPosZ) mTrueMomentumPosZ = new float;
  *mTrueMomentumPosX = aMom.X();
  *mTrueMomentumPosY = aMom.Y();
  *mTrueMomentumPosZ = aMom.Z();
}

//_________________
void MpdFemtoModelHiddenInfo::setTrueMomentumPos(const double& aPx, const double& aPy,
						 const double& aPz) {
  // Set momentum from components
  if (!mTrueMomentumPosX) mTrueMomentumPosX = new float;
  if (!mTrueMomentumPosY) mTrueMomentumPosY = new float;
  if (!mTrueMomentumPosZ) mTrueMomentumPosZ = new float;
  *mTrueMomentumPosX = aPx;
  *mTrueMomentumPosY = aPy;
  *mTrueMomentumPosZ = aPz;
}

//_________________
void MpdFemtoModelHiddenInfo::setEmissionPointPos(TLorentzVector *aPos) {
  // Set emission point
  if (!mEmissionPointPosX) mEmissionPointPosX = new float;
  if (!mEmissionPointPosY) mEmissionPointPosY = new float;
  if (!mEmissionPointPosZ) mEmissionPointPosZ = new float;
  if (!mEmissionPointPosT) mEmissionPointPosT = new float;

  if (aPos) {
    *mEmissionPointPosX = aPos->X();
    *mEmissionPointPosY = aPos->Y();
    *mEmissionPointPosZ = aPos->Z();
    *mEmissionPointPosT = aPos->T();
  } else {
    delete mEmissionPointPosX;
    mEmissionPointPosX = nullptr;
    delete mEmissionPointPosY;
    mEmissionPointPosY = nullptr;
    delete mEmissionPointPosZ;
    mEmissionPointPosZ = nullptr;
    delete mEmissionPointPosT;
    mEmissionPointPosT = nullptr;
  }
}

//_________________
void MpdFemtoModelHiddenInfo::setEmissionPointPos(const TLorentzVector& aPos) {
  // Set position from vector
  if (!mEmissionPointPosX) mEmissionPointPosX = new float;
  if (!mEmissionPointPosY) mEmissionPointPosY = new float;
  if (!mEmissionPointPosZ) mEmissionPointPosZ = new float;
  if (!mEmissionPointPosT) mEmissionPointPosT = new float;

  *mEmissionPointPosX = aPos.X();
  *mEmissionPointPosY = aPos.Y();
  *mEmissionPointPosZ = aPos.Z();
  *mEmissionPointPosT = aPos.T();
}

//_________________
void MpdFemtoModelHiddenInfo::setEmissionPointPos(const double& aRx, const double& aRy,
						  const double& aRz, const double& aT) {
  // Set position from components
  if (!mEmissionPointPosX) mEmissionPointPosX = new float;
  if (!mEmissionPointPosY) mEmissionPointPosY = new float;
  if (!mEmissionPointPosZ) mEmissionPointPosZ = new float;
  if (!mEmissionPointPosT) mEmissionPointPosT = new float;

  *mEmissionPointPosX = aRx;
  *mEmissionPointPosY = aRy;
  *mEmissionPointPosZ = aRz;
  *mEmissionPointPosT = aT;
}

//_________________
TVector3 MpdFemtoModelHiddenInfo::trueMomentumNeg() const {
  TVector3 mom;
  if (!mTrueMomentumNegX || !mTrueMomentumNegY || !mTrueMomentumNegZ) {
    if (!mTrueMomentumNegX) {
      // std::cout << "mTrueMomentumNegX does not exist" << std::endl;
    }
    if (!mTrueMomentumNegY) {
      // std::cout << "mTrueMomentumNegY does not exist" << std::endl;
    }
    if (!mTrueMomentumNegZ) {
      // std::cout << "mTrueMomentumNegZ does not exist" << std::endl;
    }
    mom.SetXYZ(0, 0, 0);
  } else {
    mom.SetXYZ(*mTrueMomentumNegX, *mTrueMomentumPosY, *mTrueMomentumNegZ);
  }
  return mom;
}

//_________________
TLorentzVector MpdFemtoModelHiddenInfo::emissionPointNeg() const {
  TLorentzVector point;
  if (!mEmissionPointNegX || !mEmissionPointNegY || !mEmissionPointNegZ || !mEmissionPointNegT) {
    if (!mEmissionPointNegX) {
      // std::cout << "mEmissionPointNegX does not exist" << std::endl;
    }
    if (!mEmissionPointNegY) {
      // std::cout << "mEmissionPointNegY does not exist" << std::endl;
    }
    if (!mEmissionPointNegZ) {
      // std::cout << "mEmissionPointNegZ does not exist" << std::endl;
    }
    if (!mEmissionPointNegT) {
      // std::cout << "mEmissionPointNegT does not exist" << std::endl;
    }
    point.SetXYZT(0, 0, 0, 0);
  } else {
    point.SetXYZT(*mEmissionPointNegX, *mEmissionPointNegY, *mEmissionPointNegY, *mEmissionPointNegT);
  }
  return point;
}

//_________________
void MpdFemtoModelHiddenInfo::setTrueMomentumNeg(TVector3 *aMom) {
  // Set momentum of the positive daughter
  if (!mTrueMomentumNegX) mTrueMomentumNegX = new float;
  if (!mTrueMomentumNegY) mTrueMomentumNegY = new float;
  if (!mTrueMomentumNegZ) mTrueMomentumNegZ = new float;

  if (aMom) {
    *mTrueMomentumNegX = aMom->X();
    *mTrueMomentumNegY = aMom->Y();
    *mTrueMomentumNegZ = aMom->Z();
  } else {
    delete mTrueMomentumNegX;
    mTrueMomentumNegX = nullptr;
    delete mTrueMomentumNegY;
    mTrueMomentumNegY = nullptr;
    delete mTrueMomentumNegZ;
    mTrueMomentumNegZ = nullptr;
  }
}

//_________________
void MpdFemtoModelHiddenInfo::setTrueMomentumNeg(const TVector3& aMom) {
  // Set momentum from vector
  if (!mTrueMomentumNegX) mTrueMomentumNegX = new float;
  if (!mTrueMomentumNegY) mTrueMomentumNegY = new float;
  if (!mTrueMomentumNegZ) mTrueMomentumNegZ = new float;
  *mTrueMomentumNegX = aMom.X();
  *mTrueMomentumNegY = aMom.Y();
  *mTrueMomentumNegZ = aMom.Z();
}

//_________________
void MpdFemtoModelHiddenInfo::setTrueMomentumNeg(const double& aPx, const double& aPy,
						 const double& aPz) {
  // Set momentum from components
  if (!mTrueMomentumNegX) mTrueMomentumNegX = new float;
  if (!mTrueMomentumNegY) mTrueMomentumNegY = new float;
  if (!mTrueMomentumNegZ) mTrueMomentumNegZ = new float;
  *mTrueMomentumNegX = aPx;
  *mTrueMomentumNegY = aPy;
  *mTrueMomentumNegZ = aPz;
}

//_________________
void MpdFemtoModelHiddenInfo::setEmissionPointNeg(TLorentzVector *aNeg) {
  // Set emission point
  if (!mEmissionPointNegX) mEmissionPointNegX = new float;
  if (!mEmissionPointNegY) mEmissionPointNegY = new float;
  if (!mEmissionPointNegZ) mEmissionPointNegZ = new float;
  if (!mEmissionPointNegT) mEmissionPointNegT = new float;

  if (aNeg) {
    *mEmissionPointNegX = aNeg->X();
    *mEmissionPointNegY = aNeg->Y();
    *mEmissionPointNegZ = aNeg->Z();
    *mEmissionPointNegT = aNeg->T();
  } else {
    delete mEmissionPointNegX;
    mEmissionPointNegX = nullptr;
    delete mEmissionPointNegY;
    mEmissionPointNegY = nullptr;
    delete mEmissionPointNegZ;
    mEmissionPointNegZ = nullptr;
    delete mEmissionPointNegT;
    mEmissionPointNegT = nullptr;
  }
}

//_________________
void MpdFemtoModelHiddenInfo::setEmissionPointNeg(const TLorentzVector& aNeg) {
  // Set position from vector
  if (!mEmissionPointNegX) mEmissionPointNegX = new float;
  if (!mEmissionPointNegY) mEmissionPointNegY = new float;
  if (!mEmissionPointNegZ) mEmissionPointNegZ = new float;
  if (!mEmissionPointNegT) mEmissionPointNegT = new float;

  *mEmissionPointNegX = aNeg.X();
  *mEmissionPointNegY = aNeg.Y();
  *mEmissionPointNegZ = aNeg.Z();
  *mEmissionPointNegT = aNeg.T();
}

//_________________
void MpdFemtoModelHiddenInfo::setEmissionPointNeg(const double& aRx, const double& aRy, const double& aRz, const double& aT) {
  // Set position from components
  if (!mEmissionPointNegX) mEmissionPointNegX = new float;
  if (!mEmissionPointNegY) mEmissionPointNegY = new float;
  if (!mEmissionPointNegZ) mEmissionPointNegZ = new float;
  if (!mEmissionPointNegT) mEmissionPointNegT = new float;

  *mEmissionPointNegX = aRx;
  *mEmissionPointNegY = aRy;
  *mEmissionPointNegZ = aRz;
  *mEmissionPointNegT = aT;
}

ClassImp(MpdFemtoModelHiddenInfo)

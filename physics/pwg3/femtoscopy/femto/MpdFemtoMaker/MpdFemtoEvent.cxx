//
// MpdFemtoEvent holds the event information specific and particle collections
//

// MpdFemtoMaker headers
// Infrastructure
#include "MpdFemtoEvent.h"
#include "MpdFemtoTrack.h"
#include "MpdFemtoV0.h"
#include "MpdFemtoXi.h"
#include "MpdFemtoKink.h"
#include "PhysicalConstants.h"
#include "SystemOfUnits.h"
// Base
#include "MpdFemtoBaseTrackCut.h"
#include "MpdFemtoBaseV0Cut.h"
#include "MpdFemtoBaseXiCut.h"
#include "MpdFemtoBaseKinkCut.h"

ClassImp(MpdFemtoEvent)

//_________________
MpdFemtoEvent::MpdFemtoEvent() :
mEventNumber(0), mRunNumber(0), mMagneticField(0),
mRefMult(0), mRefMultPos(0), mRefMultCorr(0),
mRefMultCorrWeight(0), mRefMult2(0), mRefMult2Pos(0),
mGRefMult(0), mGRefMultPos(0), mBTofTrayMultiplicity(0),
mNBTOFMatch(0), mNBEMCMatch(0), mNumberOfPrimaryTracks(0),
mNumberOfGlobalTracks(0), mZdcSumAdcEast(0), mZdcSumAdcWest(0),
mZdcCoincidenceRate(0), mBbcCoincidenceRate(0), mSphericity(-1), mSphericity2(-1),
mEventPlaneAngle(0), mEventPlaneResolution(0), mCent16(-1),
mPrimaryVertexPositionX(-999), mPrimaryVertexPositionY(-999), mPrimaryVertexPositionZ(-999),
mVpdVz(0), mRanking(-1e5), mImpactPar(0) {

  // Default constuctor
  if (!mTriggerIds.empty()) mTriggerIds.clear();

  mTrackCollection = new MpdFemtoTrackCollection;
  mV0Collection = new MpdFemtoV0Collection;
  mXiCollection = new MpdFemtoXiCollection;
  mKinkCollection = new MpdFemtoKinkCollection;
}

//___________________
MpdFemtoEvent::MpdFemtoEvent(const MpdFemtoEvent& ev, MpdFemtoBaseTrackCut* tCut,
			     MpdFemtoBaseV0Cut* vCut, MpdFemtoBaseXiCut* xCut,
			     MpdFemtoBaseKinkCut* kCut) {

  // Copy constructor with track and v0 cuts
  mEventNumber = ev.mEventNumber;
  mRunNumber = ev.mRunNumber;
  mMagneticField = ev.mMagneticField;

  mRefMult = ev.mRefMult;
  mRefMultPos = ev.mRefMultPos;
  mRefMultCorr = ev.mRefMultCorr;
  mRefMultCorrWeight = ev.mRefMultCorrWeight;
  mRefMult2 = ev.mRefMult2;
  mRefMult2Pos = ev.mRefMult2Pos;
  mGRefMult = ev.mGRefMult;
  mGRefMultPos = ev.mGRefMultPos;
  mBTofTrayMultiplicity = ev.mBTofTrayMultiplicity;
  mNBTOFMatch = ev.mNBTOFMatch;
  mNBEMCMatch = ev.mNBEMCMatch;
  mNumberOfPrimaryTracks = ev.mNumberOfPrimaryTracks;
  mNumberOfGlobalTracks = ev.mNumberOfGlobalTracks;
  mZdcSumAdcEast = ev.mZdcSumAdcEast;
  mZdcSumAdcWest = ev.mZdcSumAdcWest;
  mZdcCoincidenceRate = ev.mZdcCoincidenceRate;
  mBbcCoincidenceRate = ev.mBbcCoincidenceRate;

  mSphericity = ev.mSphericity;
  mSphericity2 = ev.mSphericity2;
  mEventPlaneAngle = ev.mEventPlaneAngle;
  mEventPlaneResolution = ev.mEventPlaneResolution;
  mCent16 = ev.mCent16;

  mPrimaryVertexPositionX = ev.mPrimaryVertexPositionX;
  mPrimaryVertexPositionY = ev.mPrimaryVertexPositionY;
  mPrimaryVertexPositionZ = ev.mPrimaryVertexPositionZ;
  mVpdVz = ev.mVpdVz;
  mRanking = ev.mRanking;
  mImpactPar = ev.mImpactPar;
  mTriggerIds = ev.mTriggerIds;

  // Create empty collections
  mTrackCollection = new MpdFemtoTrackCollection;
  mV0Collection = new MpdFemtoV0Collection;
  mXiCollection = new MpdFemtoXiCollection;
  mKinkCollection = new MpdFemtoKinkCollection;

  // Copy tracks from one collection to the another one
  for (MpdFemtoTrackIterator tIter = ev.mTrackCollection->begin();
       tIter != ev.mTrackCollection->end(); tIter++) {
    if (!tCut || tCut->pass(*tIter)) {
      MpdFemtoTrack* trackCopy = new MpdFemtoTrack(**tIter);
      mTrackCollection->push_back(trackCopy);
    }
  }

  // Copy V0s from one collection to the another one
  for (MpdFemtoV0Iterator vIter = ev.mV0Collection->begin();
       vIter != ev.mV0Collection->end(); vIter++) {
    if (!vCut || vCut->pass(*vIter)) {
      MpdFemtoV0* v0Copy = new MpdFemtoV0(**vIter);
      mV0Collection->push_back(v0Copy);
    }
  }

  // Copy Xis from one collection to the another one
  for (MpdFemtoXiIterator xIter = ev.mXiCollection->begin();
       xIter != ev.mXiCollection->end(); xIter++) {
    if (!xCut || xCut->pass(*xIter)) {
      MpdFemtoXi* xiCopy = new MpdFemtoXi(**xIter);
      mXiCollection->push_back(xiCopy);
    }
  }

  // Copy kinks from one collection to the another one
  for (MpdFemtoKinkIterator kIter = ev.mKinkCollection->begin();
       kIter != ev.mKinkCollection->end(); kIter++) {
    if (!kCut || kCut->pass(*kIter)) {
      MpdFemtoKink* kinkCopy = new MpdFemtoKink(**kIter);
      mKinkCollection->push_back(kinkCopy);
    }
  }
}

//_________________
MpdFemtoEvent& MpdFemtoEvent::operator=(const MpdFemtoEvent& ev) {

  // Assignment operator
  if (this != &ev) {
    mEventNumber = ev.mEventNumber;
    mRunNumber = ev.mRunNumber;
    mMagneticField = ev.mMagneticField;

    mRefMult = ev.mRefMult;
    mRefMultPos = ev.mRefMultPos;
    mRefMultCorr = ev.mRefMultCorr;
    mRefMultCorrWeight = ev.mRefMultCorrWeight;
    mRefMult2 = ev.mRefMult2;
    mRefMult2Pos = ev.mRefMult2Pos;
    mGRefMult = ev.mGRefMult;
    mGRefMultPos = ev.mGRefMultPos;
    mBTofTrayMultiplicity = ev.mBTofTrayMultiplicity;
    mNBTOFMatch = ev.mNBTOFMatch;
    mNBEMCMatch = ev.mNBEMCMatch;
    mNumberOfPrimaryTracks = ev.mNumberOfPrimaryTracks;
    mNumberOfGlobalTracks = ev.mNumberOfGlobalTracks;
    mZdcSumAdcEast = ev.mZdcSumAdcEast;
    mZdcSumAdcWest = ev.mZdcSumAdcWest;
    mZdcCoincidenceRate = ev.mZdcCoincidenceRate;
    mBbcCoincidenceRate = ev.mBbcCoincidenceRate;
    mSphericity = ev.mSphericity;
    mSphericity2 = ev.mSphericity2;
    mEventPlaneAngle = ev.mEventPlaneAngle;
    mEventPlaneResolution = ev.mEventPlaneResolution;
    mCent16 = ev.mCent16;
    mPrimaryVertexPositionX = ev.mPrimaryVertexPositionX;
    mPrimaryVertexPositionY = ev.mPrimaryVertexPositionY;
    mPrimaryVertexPositionZ = ev.mPrimaryVertexPositionZ;
    mRanking = ev.mRanking;
    mImpactPar = ev.mImpactPar;

    for (unsigned int iIter = 0; iIter < ev.mTriggerIds.size(); iIter++) {
      mTriggerIds.push_back(ev.mTriggerIds.at(iIter));
    }

    // Clear collections if exist
    if (mTrackCollection) {
      for (MpdFemtoTrackIterator iter = mTrackCollection->begin();
	   iter != mTrackCollection->end(); iter++) {
	delete *iter;
      }
      mTrackCollection->clear();
    } else {
      mTrackCollection = new MpdFemtoTrackCollection;
    }

    if (mV0Collection) {
      for (MpdFemtoV0Iterator iter = mV0Collection->begin();
	   iter != mV0Collection->end(); iter++) {
	delete *iter;
      }
      mV0Collection->clear();
    } else {
      mV0Collection = new MpdFemtoV0Collection;
    }

    if (mKinkCollection) {
      for (MpdFemtoKinkIterator iter = mKinkCollection->begin();
	   iter != mKinkCollection->end(); iter++) {
	delete *iter;
      }
      mKinkCollection->clear();
    } else {
      mKinkCollection = new MpdFemtoKinkCollection;
    }

    if (mXiCollection) {
      for (MpdFemtoXiIterator iter = mXiCollection->begin();
	   iter != mXiCollection->end(); iter++) {
	delete *iter;
      }
      mXiCollection->clear();
    } else {
      mXiCollection = new MpdFemtoXiCollection;
    }

    // Copy collections
    for (MpdFemtoTrackIterator iIter = ev.mTrackCollection->begin();
	 iIter != ev.mTrackCollection->end(); iIter++) {
      MpdFemtoTrack *trackCopy = new MpdFemtoTrack(**iIter);
      mTrackCollection->push_back(trackCopy);
    }

    for (MpdFemtoV0Iterator iIter = ev.mV0Collection->begin();
	 iIter != ev.mV0Collection->end(); iIter++) {
      MpdFemtoV0 *v0Copy = new MpdFemtoV0(**iIter);
      mV0Collection->push_back(v0Copy);
    }

    for (MpdFemtoKinkIterator iIter = ev.mKinkCollection->begin();
	 iIter != ev.mKinkCollection->end(); iIter++) {
      MpdFemtoKink *kinkCopy = new MpdFemtoKink(**iIter);
      mKinkCollection->push_back(kinkCopy);
    }

    for (MpdFemtoXiIterator iIter = ev.mXiCollection->begin();
	 iIter != ev.mXiCollection->end(); iIter++) {
      MpdFemtoXi *xiCopy = new MpdFemtoXi(**iIter);
      mXiCollection->push_back(xiCopy);
    }
  } //if ( this != &ev )

  return *this;
}

//___________________
MpdFemtoEvent::~MpdFemtoEvent() {

  // Remove track collection
  for (MpdFemtoTrackIterator iter = mTrackCollection->begin();
       iter != mTrackCollection->end(); iter++) {
    delete *iter;
  }
  mTrackCollection->clear();
  delete mTrackCollection;

  // Remove V0 collection
  for (MpdFemtoV0Iterator V0iter = mV0Collection->begin();
       V0iter != mV0Collection->end(); V0iter++) {
    delete *V0iter;
  }
  mV0Collection->clear();
  delete mV0Collection;

  // Remove Xi collection
  for (MpdFemtoXiIterator XiIter = mXiCollection->begin();
       XiIter != mXiCollection->end(); XiIter++) {
    delete *XiIter;
  }
  mXiCollection->clear();
  delete mXiCollection;

  // Remove Kink collection
  for (MpdFemtoKinkIterator kinkIter = mKinkCollection->begin();
       kinkIter != mKinkCollection->end(); kinkIter++) {
    delete *kinkIter;
  }
  mKinkCollection->clear();
  delete mKinkCollection;
}

//___________________
void MpdFemtoEvent::rotateZ(const double& angle) {

  MpdFemtoTrackIterator iter;
  MpdFemtoV0Iterator V0iter;

  MpdFemtoPhysicalHelix helix;
  TVector3 p;
  TVector3 o;

  mEventPlaneAngle += angle;
  std::cout << " MpdFemtoEvent::rotateZ(const double angle) - angle="
            << angle << " rad    ";
  std::cout << angle / degree << " deg " << std::endl;
  for (iter = mTrackCollection->begin(); iter != mTrackCollection->end(); iter++) {

    p = (*iter)->p();
    p.RotateZ(angle);
    (*iter)->setP(p);
  }
  for (V0iter = mV0Collection->begin(); V0iter != mV0Collection->end(); V0iter++) {

    p = (*V0iter)->decayPoint();
    p.RotateX(angle);
    (*V0iter)->setDecayPoint(p);

    p = (*V0iter)->momV0();
    p.RotateX(angle);
    (*V0iter)->setV0Mom(p);

    p = (*V0iter)->momPos();
    p.RotateX(angle);
    (*V0iter)->setMomPos(p);

    p = (*V0iter)->momNeg();
    p.RotateX(angle);
    (*V0iter)->setMomNeg(p);
  }
}

//_________________
bool MpdFemtoEvent::isTrigger(const unsigned int& id) const {
  return std::find(mTriggerIds.begin(), mTriggerIds.end(), id) != mTriggerIds.end();
}

//_________________
short MpdFemtoEvent::cent9() const {
  if (mCent16 == 15) {
    return 8; //0-5%
  } else if (mCent16 == 14) {
    return 7; //5-10%
  } else if ((mCent16 == 13) || (mCent16 == 12)) {
    return 6; //10-20%
  } else if ((mCent16 == 11) || (mCent16 == 10)) {
    return 5; //20-30%
  } else if ((mCent16 == 9) || (mCent16 == 8)) {
    return 4; //30-40%
  } else if ((mCent16 == 7) || (mCent16 == 6)) {
    return 3; //40-50%
  } else if ((mCent16 == 5) || (mCent16 == 4)) {
    return 2; //50-60%
  } else if ((mCent16 == 3) || (mCent16 == 2)) {
    return 1; //60-70%
  } else if ((mCent16 == 1) || (mCent16 == 0)) {
    return 0; //70-80%
  } else {
    return -1; //not defined
  }
}

//_________________
void MpdFemtoEvent::addTriggerId(const unsigned int& id) {
  // Check if the trigger id is not in the mTriggerIds vector
  if (std::find(mTriggerIds.begin(),
		mTriggerIds.end(),
		id) == mTriggerIds.end()) {
    mTriggerIds.push_back(id);
  }
}

//_________________
void MpdFemtoEvent::setSphericity(const float& sph) {
  // Set sphericity value for |eta|<0.5, nHits>10, pT>0.15 and DCA<3
  if (sph < 0) {
    mSphericity = -1;
  } else {
    mSphericity = ((sph * 100.) > 100. ? -1. : (char) (sph * 100.));
  }
}

//_________________
void MpdFemtoEvent::setSphericity2(const float& sph) {
  // Set sphericity value for |eta|<0.5, nHits>10, pT>0.15 and DCA<3
  if (sph < 0) {
    mSphericity2 = -1;
  } else {
    mSphericity2 = ((sph * 100.) > 100. ? -1. : (char) (sph * 100.));
  }
}

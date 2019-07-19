//
// Main class holding track information
//

// MpdFemtoMaker headers
#include "MpdFemtoTrack.h"

//________________
MpdFemtoTrack::MpdFemtoTrack() :
mId(0), mFlag(0), mNHits(0), mNHitsPoss(0), mNHitsDedx(0), mChi2(0), mDedx(0),
mNSigmaElectron(-30), mNSigmaPion(-30), mNSigmaKaon(-30), mNSigmaProton(-30),
mPidProbElectron(0), mPidProbPion(0), mPidProbKaon(0), mPidProbProton(0),
mMap{}

, mTofBeta(0),
mPrimaryPx(0), mPrimaryPy(0), mPrimaryPz(0), mGlobalPx(0), mGlobalPy(0), mGlobalPz(0),
mDcaX(-999), mDcaY(-999), mDcaZ(-999),
mPrimaryVertexX(0), mPrimaryVertexY(0), mPrimaryVertexZ(0),
/* mXfr(0), mYfr(0), mZfr(0), mTfr(0), mPdgId(0), */
mHiddenInfo(nullptr) {
    // Default constructor
}

//________________

MpdFemtoTrack::MpdFemtoTrack(const MpdFemtoTrack& t) {

    // Copy constructor
    mId = t.mId;
    mFlag = t.mFlag;
    mNHits = t.mNHits;
    mNHitsPoss = t.mNHitsPoss;
    mNHitsDedx = t.mNHitsDedx;
    mChi2 = t.mChi2;
    mDedx = t.mDedx;
    mNSigmaElectron = t.mNSigmaElectron;
    mNSigmaPion = t.mNSigmaPion;
    mNSigmaKaon = t.mNSigmaKaon;
    mNSigmaProton = t.mNSigmaProton;
    mPidProbElectron = t.mPidProbElectron;
    mPidProbPion = t.mPidProbPion;
    mPidProbKaon = t.mPidProbKaon;
    mPidProbProton = t.mPidProbProton;
    mDcaX = t.mDcaX;
    mDcaY = t.mDcaY;
    mDcaZ = t.mDcaZ;
    mMap[0] = t.mMap[0];
    mMap[1] = t.mMap[1];
    mTofBeta = t.mTofBeta;
    mPrimaryPx = t.mPrimaryPx;
    mPrimaryPy = t.mPrimaryPy;
    mPrimaryPz = t.mPrimaryPz;
    mGlobalPx = t.mGlobalPx;
    mGlobalPy = t.mGlobalPy;
    mGlobalPz = t.mGlobalPz;
    mPrimaryVertexX = t.mPrimaryVertexX;
    mPrimaryVertexY = t.mPrimaryVertexY;
    mPrimaryVertexZ = t.mPrimaryVertexZ;
    /*
    mXfr = t.mXfr;
    mYfr = t.mYfr;
    mZfr = t.mZfr;
    mTfr = t.mTfr;
    mPdgId = t.mPdgId;
     */

    if (t.validHiddenInfo()) {
        mHiddenInfo = t.getHiddenInfo()->clone();
    } else {
        mHiddenInfo = nullptr;
    }
}

//_________________

MpdFemtoTrack& MpdFemtoTrack::operator=(const MpdFemtoTrack& trk) {

    // Assignment operator
    if (this != &trk) {
        mId = trk.mId;
        mFlag = trk.mFlag;
        mNHits = trk.mNHits;
        mNHitsPoss = trk.mNHitsPoss;
        mNHitsDedx = trk.mNHitsDedx;
        mChi2 = trk.mChi2;
        mDedx = trk.mDedx;
        mNSigmaElectron = trk.mNSigmaElectron;
        mNSigmaPion = trk.mNSigmaPion;
        mNSigmaKaon = trk.mNSigmaKaon;
        mNSigmaProton = trk.mNSigmaProton;
        mPidProbElectron = trk.mPidProbElectron;
        mPidProbPion = trk.mPidProbPion;
        mPidProbKaon = trk.mPidProbKaon;
        mPidProbProton = trk.mPidProbProton;
        mMap[0] = trk.mMap[0];
        mMap[1] = trk.mMap[1];
        mTofBeta = trk.mTofBeta;
        mPrimaryPx = trk.mPrimaryPx;
        mPrimaryPy = trk.mPrimaryPy;
        mPrimaryPz = trk.mPrimaryPz;
        mGlobalPx = trk.mGlobalPx;
        mGlobalPy = trk.mGlobalPy;
        mGlobalPz = trk.mGlobalPz;
        mDcaX = trk.mDcaX;
        mDcaY = trk.mDcaY;
        mDcaZ = trk.mDcaZ;
        mPrimaryVertexX = trk.mPrimaryVertexX;
        mPrimaryVertexY = trk.mPrimaryVertexY;
        mPrimaryVertexZ = trk.mPrimaryVertexZ;
        mBField = trk.mBField;

        /*
        mXfr = trk.mXfr;
        mYfr = trk.mYfr;
        mZfr = trk.mZfr;
        mTfr = trk.mTfr;
        mPdgId = trk.mPdgId;
         */

        if (mHiddenInfo) delete mHiddenInfo;
        mHiddenInfo = trk.validHiddenInfo() ? trk.getHiddenInfo()->clone() : nullptr;
    }

    return *this;
}

//_________________

MpdFemtoTrack::~MpdFemtoTrack() {
    if (mHiddenInfo) delete mHiddenInfo;
}

//_________________

float MpdFemtoTrack::massSqr() const {
    // Set squared mass
    float massSqr = -999.;
    if (isTofTrack()) {
        massSqr = gPtot2() * (invBeta2() - 1.);
    }
    return massSqr;
}

//_________________

void MpdFemtoTrack::setNSigmaElectron(const float& ns) {
    // Set nSigma(e)
    mNSigmaElectron = (TMath::Abs(ns * 1000.) > std::numeric_limits<short>::max() ?
            ((ns > 0) ? std::numeric_limits<short>::max() : std::numeric_limits<short>::min()) :
            (short) (ns * 1000.));
}

//_________________

void MpdFemtoTrack::setNSigmaPion(const float& ns) {
    // Set nSigma(pi)
    mNSigmaPion = (TMath::Abs(ns * 1000.) > std::numeric_limits<short>::max() ?
            ((ns > 0) ? std::numeric_limits<short>::max() : std::numeric_limits<short>::min()) :
            (short) (ns * 1000.));
}

//_________________

void MpdFemtoTrack::setNSigmaKaon(const float& ns) {
    // Set nSigma(K)
    mNSigmaKaon = (TMath::Abs(ns * 1000.) > std::numeric_limits<short>::max() ?
            ((ns > 0) ? std::numeric_limits<short>::max() : std::numeric_limits<short>::min()) :
            (short) (ns * 1000.));
}

//_________________

void MpdFemtoTrack::setNSigmaProton(const float& ns) {
    // Set nSigma(p)
    mNSigmaProton = (TMath::Abs(ns * 1000.) > std::numeric_limits<short>::max() ?
            ((ns > 0) ? std::numeric_limits<short>::max() : std::numeric_limits<short>::min()) :
            (short) (ns * 1000.));
}

//_________________

void MpdFemtoTrack::setChi2(const float& x) {
    // Set chi2
    mChi2 = ((x * 1000.) > std::numeric_limits<unsigned short>::max() ?
            std::numeric_limits<unsigned short>::max() :
            (unsigned short) (x * 1000.));
}

//_________________

void MpdFemtoTrack::setPidProbElectron(const float& prob) {
    // Set probability(e)
    mPidProbElectron = ((prob * 10000.) > std::numeric_limits<unsigned short>::max() ?
            std::numeric_limits<unsigned short>::max() :
            (unsigned short) (prob * 10000.));
}

//_________________

void MpdFemtoTrack::setPidProbPion(const float& prob) {
    // Set probability(pi)
    if (prob < 0) {
        mPidProbPion = std::numeric_limits<unsigned short>::max();
    } else {
        mPidProbPion = ((prob * 10000.) > std::numeric_limits<unsigned short>::max() ?
                std::numeric_limits<unsigned short>::max() :
                (unsigned short) (prob * 10000.));
    }
}

//_________________

void MpdFemtoTrack::setPidProbKaon(const float& prob) {
    // Set probability(K)
    if (prob < 0) {
        mPidProbKaon = std::numeric_limits<unsigned short>::max();
    } else {
        mPidProbKaon = ((prob * 10000.) > std::numeric_limits<unsigned short>::max() ?
                std::numeric_limits<unsigned short>::max() :
                (unsigned short) (prob * 10000.));
    }
}

//_________________

void MpdFemtoTrack::setPidProbProton(const float& prob) {
    // Set probability(p)
    if (prob < 0) {
        mPidProbProton = std::numeric_limits<unsigned short>::max();
    } else {
        mPidProbProton = ((prob * 10000.) > std::numeric_limits<unsigned short>::max() ?
                std::numeric_limits<unsigned short>::max() :
                (unsigned short) (prob * 10000.));
    }
}

//_________________

void MpdFemtoTrack::setDedx(const double& dEdx) {
    // Set dE/dx (from GeV/cm)
    if (dEdx < 0) {
        mDedx = 0;
    } else {
        mDedx = ((dEdx * 1e9) > std::numeric_limits<unsigned short>::max() ?
                std::numeric_limits<unsigned short>::max() :
                (unsigned short) (dEdx * 1e9));
    }
}

//_________________

void MpdFemtoTrack::setDedxFromKeV(const double& dEdx) {
    // Set dE/dx (from keV/cm)
    if (dEdx < 0) {
        mDedx = 0;
    } else {
        mDedx = ((dEdx * 1e3) > std::numeric_limits<unsigned short>::max() ?
                std::numeric_limits<unsigned short>::max() :
                (unsigned short) (dEdx * 1e3));
    }
}

//_________________

void MpdFemtoTrack::setBeta(const float& beta) {
    // Set velocity (from TOF)
    if (beta <= 0) {
        mTofBeta = -666;
    } else {
        mTofBeta = ((beta * 20000.) > std::numeric_limits<unsigned short>::max() ?
                std::numeric_limits<unsigned short>::max() :
                (unsigned short) (beta * 20000.));
    }
}

//_________________

MpdFemtoPhysicalHelix MpdFemtoTrack::helix() const {
    // Return helix of the primary track
    return MpdFemtoPhysicalHelix(pMom(), primaryVertex(),
            mBField * kilogauss,
            static_cast<float> (charge()));
}

//_________________

MpdFemtoPhysicalHelix MpdFemtoTrack::gHelix() const {
    // Return helix of the global track
    return MpdFemtoPhysicalHelix(gMom(), origin(),
            mBField * kilogauss,
            static_cast<float> (charge()));
}

ClassImp(MpdFemtoTrack)
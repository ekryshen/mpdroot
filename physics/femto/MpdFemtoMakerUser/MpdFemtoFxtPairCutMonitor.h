/**
 * \class MpdFemtoFxtPairCutMonitor
 * \brief Pair cut monitor for basic analysis
 *
 * The class provides histograms for monitoring pair cuts
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * email: nigmatkulov@gmail.com
 */

#ifndef MpdFemtoFxtPairCutMonitor_h
#define MpdFemtoFxtPairCutMonitor_h

// MpdFemtoMaker headers
// Base
#include "MpdFemtoBaseCutMonitor.h"

// ROOT headers
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"

// Forward declaration
class MpdFemtoPair;

//_________________
class MpdFemtoFxtPairCutMonitor : public MpdFemtoBaseCutMonitor {
 public:
  /// Constructor
  MpdFemtoFxtPairCutMonitor(const char*);
  /// Copy constructor
  MpdFemtoFxtPairCutMonitor(const MpdFemtoFxtPairCutMonitor& cutMoni);
  /// Assignment operator
  MpdFemtoFxtPairCutMonitor operator=(const MpdFemtoFxtPairCutMonitor& c);
  /// Destructor
  virtual ~MpdFemtoFxtPairCutMonitor();

  /// Fill pair cut histograms
  virtual void fill(const MpdFemtoPair* pair);

  /// Write all histograms
  void writeOutHistos();

  /// Ouput list with histograms
  virtual TList* getOutputList();

  // These dummy Fill() functions were introduced to remove a compiler
  // warning related to overloaded base-class Fill() functions being
  // hidden by a single version of Fill() in this derived class
  void fill(const MpdFemtoEvent* /* d */) {
    ;
  }

  void fill(const MpdFemtoTrack* /* d */) {
    ;
  }

  void fill(const MpdFemtoV0* /* d */) {
    ;
  }

  void fill(const MpdFemtoXi* /* xi */) {
    ;
  }

  void fill(const MpdFemtoKink* /* d */) {
    ;
  }

  void fill(const MpdFemtoParticleCollection* /* d */) {
    ;
  }

  void fill(const MpdFemtoEvent* /* d1 */, const MpdFemtoParticleCollection* /* d2 */) {
    ;
  }

  void fill(const MpdFemtoParticleCollection* /* d1 */, const MpdFemtoParticleCollection* /* d2 */) {
    ;
  }

  /// Half transverse momentum of a pair
  TH1F* kt() {
    return (mKt) ? mKt : nullptr;
  }
  /// Half transverse momentum of a pair
  TH1F* Kt() {
    return Kt();
  }
  /// pT1 vs. pT2 vs. qInv
  TH3F* pt1pt2qInv() {
    return (mPt1Pt2Qinv) ? mPt1Pt2Qinv : nullptr;
  }
  /// pT1 vs. pT2 vs. qInv
  TH3F* Pt1Pt2Qinv() {
    return pt1pt2qInv();
  }
  /// Fraction of merged rows vs. qInv
  TH2F* fmrVsQinv() {
    return (mFMRvsQinv) ? mFMRvsQinv : nullptr;
  }
  /// Fraction of merged rows vs. qInv
  TH2F* FMRvsQinv() {
    return fmrVsQinv();
  }
  /// Splitting level (quality) vs. qInv
  TH2F* slVsQinv() {
    return (mSLvsQinv) ? mSLvsQinv : nullptr;
  }
  /// Splitting level (quality) vs. qInv
  TH2F* SLvsQinv() {
    return slVsQinv();
  }
  /// Average separation vs. qInv
  TH2F* aveSepVsQinv() {
    return (mAveSepVsQinv) ? mAveSepVsQinv : nullptr;
  }
  /// Average separation vs. qInv
  TH2F* AveSepVsQinv() {
    return aveSepVsQinv();
  }
  /// sqrt(deta^2+dphi^2) vs. qInv
  TH2F* rValueVsQinv() {
    return (mRValueVsQinv) ? mRValueVsQinv : nullptr;
  }
  /// sqrt(deta^2+dphi^2) vs. qInv
  TH2F* RValueVsQinv() {
    return rValueVsQinv();
  }
  /// delta eta vs. delta phi
  TH2F* dEtaVsdPhi() {
    return (mDEtaDPhi) ? mDEtaDPhi : nullptr;
  }
  /// delta eta vs. delta phi
  TH2F* DEtaDPhi() {
    return dEtaVsdPhi();
  }
  /// Invariant mass distribution
  TH1F* mInv() {
    return (mMinv) ? mMinv : nullptr;
  }
  /// Invariant mass distribution
  TH1F* Minv() {
    return mInv();
  }
  /// Transverse momentum vs. pseudorapidity
  TH2F* ptVsEta() {
    return (mPtVsEta) ? mPtVsEta : nullptr;
  }
  /// Transverse momentum vs. pseudorapidity
  TH2F* PtVsEta() {
    return ptVsEta();
  }
  /// delta eta vs. delta phi* main
  TH2F* dEtaVsdPhiStarMin() {
    return (mDEtaVsDPhiStarMin) ? mDEtaVsDPhiStarMin : nullptr;
  }

  /// Set particle mass
  void setParticleMass(const double& mass) {
    mPartMass = mass;
  }
  /// Set particle mass
  void SetParticleMass(const double& mass) {
    setParticleMass(mass);
  }

 private:

  /// Half transverse momentum of a pair
  TH1F* mKt;
  /// pT1 vs. pT2 vs. qInv
  TH3F* mPt1Pt2Qinv;
  /// Fraction of merged rows vs. qInv
  TH2F* mFMRvsQinv;
  /// Splitting level vs. qInv
  TH2F* mSLvsQinv;
  /// Average two-track separation vs. qInv
  TH2F* mAveSepVsQinv;
  /// R value vs. qInv
  TH2F* mRValueVsQinv;
  /// dEta vs. dPhi
  TH2F* mDEtaDPhi;
  /// Invaraint mass
  TH1F* mMinv;
  /// pT vs. eta
  TH2F* mPtVsEta;
  /// dEta vs. dPhiStarMin
  TH2F *mDEtaVsDPhiStarMin;
  /// Particle mass
  double mPartMass;

  ClassDef(MpdFemtoFxtPairCutMonitor, 3);
};

#endif // #define MpdFemtoFxtPairCutMonitor_h

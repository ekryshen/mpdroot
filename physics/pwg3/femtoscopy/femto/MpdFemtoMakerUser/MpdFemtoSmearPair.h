/**
 * \class MpdFemtoSmearPair
 * \brief Class that returns a MpdFemtoPair with smeared momentum
 *
 * Given a MpdFemtoPair, this class provides a corresponding
 * MpdFemtoPair whose constituent MpdFemtoParticles' momenta have been
 * smeared according to a parameterization of the STAR momentum
 * resolution.
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * \email nigmatkulov@gmail.com
 */

#ifndef MpdFemtoSmearPair_h
#define MpdFemtoSmearPair_h

// MpdFemtoMaker headers
#include "MpdFemtoPair.h"

// ROOT headers
#include "TLorentzVector.h"
#include "TRandom3.h"

//_________________
class MpdFemtoSmearPair {
 public:
  /// Default constructor
  MpdFemtoSmearPair();
  /// Constructor that takes unsmeared pair
  MpdFemtoSmearPair(const MpdFemtoPair* unSmearedPair);
  /// Destructor
  virtual ~MpdFemtoSmearPair() {
    /* empty */
  }

  /// Essentially same as c'tor
  void setUnsmearedPair(const MpdFemtoPair* unSmearedPair);

  /// Access to the smeared pair
  MpdFemtoPair& smearedPair() {
    return mSmearedPair;
  }

  //========= Resolution parameters ==========

  /// pT resolution parameterized by d(pT) = Frac*pT
  void setFractionalPtRes(double val) {
    mFracPtRes = val;
  }
  /// phi resolution parameterized d(phi)= by a+b*P^alpha
  void setPhiRes_a(double val) {
    mPhi_a = val;
  }
  /// phi resolution parameterized d(phi)= by a+b*P^alpha
  void setPhiRes_b(double val) {
    mPhi_b = val;
  }
  /// phi resolution parameterized d(phi)= by a+b*P^alpha
  void setPhiRes_alpha(double val) {
    mPhi_alpha = val;
  }
  /// phi resolution parameterized by d(theta) = a+b*P^alpha
  void setThetaRes_a(double val) {
    mTheta_a = val;
  }
  /// phi resolution parameterized by d(theta) = a+b*P^alpha
  void setThetaRes_b(double val) {
    mTheta_b = val;
  }
  /// phi resolution parameterized by d(theta) = a+b*P^alpha
  void setThetaRes_alpha(double val) {
    mTheta_alpha = val;
  }

  /// Return smeared four-momentum of the pair
  TLorentzVector smearedMomentum(TLorentzVector input);

private:

  /// Pair with smeared parameters
  MpdFemtoPair mSmearedPair;
  /// The first particle
  MpdFemtoParticle mParticle1;
  /// The second particle
  MpdFemtoParticle mParticle2;

  /// Pair resolution parameters
  double mFracPtRes;
  /// Parameter phi_a
  double mPhi_a;
  /// Parameter phi_b
  double mPhi_b;
  /// Parameter phi_alpha
  double mPhi_alpha;
  /// Parameter theta_a
  double mTheta_a;
  /// Parameter theta_b
  double mTheta_b;
  /// Parameter theta_alpha
  double mTheta_alpha;

  /// Random value
  TRandom3* mHbtRandom;

  /// Setup
  void setup();

  ClassDef(MpdFemtoSmearPair, 0)
};

#endif // #define MpdFemtoSmearPair_h

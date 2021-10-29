/**
 * \class MpdFemtoHiddenInfo
 * \brief A pure virtual base class for the hidden (Monte Carlo) data.
 *
 * Hidden info stores additional information, which is not in a standard track.
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * \email nigmatkulov@gmail.com
 */

#ifndef MpdFemtoHiddenInfo_h
#define MpdFemtoHiddenInfo_h

// MpdFemtoMaker headers
#include "MpdFemtoTypes.h"

// ROOT headers
#include "TLorentzVector.h"

//_________________
class MpdFemtoHiddenInfo {
 public:
  /// Default constructor
  MpdFemtoHiddenInfo() {
    /* empty */
  }
  /// Default destructor
  virtual ~MpdFemtoHiddenInfo() {
    /* empty */
  }

  // !!! MANDATORY !!!
  // --- Copy the hidden info from MpdFemtoTrack to MpdFemtoParticle

  /// Set emission point
  virtual void setEmissionPoint(const double& x, const double& y, const double& z, const double& t) = 0;
  /// Set emission point
  virtual void SetEmissionPoint(const double& x, const double& y, const double& z, const double& t) {
    setEmissionPoint(x, y, z, t);
  }
  /// Retrieve emission point
  virtual TLorentzVector emissionPoint() const = 0;
  /// Retrieve emission point
  virtual TLorentzVector EmissionPoint() const {
    return emissionPoint();
  }
  /// Return PDG code of the particle
  virtual int pdgId() const = 0;
  /// Return PDG code of the particle
  virtual int PdgId() const {
    return pdgId();
  }
  /// Retrieve hidden information
  virtual MpdFemtoHiddenInfo* getParticleHiddenInfo() const = 0;
  /// Retrieve hidden information
  virtual MpdFemtoHiddenInfo* GetParticleHiddenInfo() const {
    return getParticleHiddenInfo();
  }
  /// Clone hidden information
  virtual MpdFemtoHiddenInfo* clone() const;
  /// Clone hidden information
  virtual MpdFemtoHiddenInfo* Clone() const {
    return clone();
  }
};

//_________________
inline MpdFemtoHiddenInfo* MpdFemtoHiddenInfo::clone() const {
  return getParticleHiddenInfo();
}

#endif // #define MpdFemtoHiddenInfo_h

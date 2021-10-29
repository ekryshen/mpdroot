/**
 * \class MpdFemtoKtPairCut
 * \brief Simple cut to select pairs within the given kT range
 *
 * A cut to select pairs with the given reduced transverse momentum, kT,
 * emission angle w.r.t. reaction plane
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 11, 2020
 * \email nigmatkulov@gmail.com
 */


#ifndef MpdFemtoKtPairCut_h
#define MpdFemtoKtPairCut_h

// MpdFemto 
#include "MpdFemtoBasePairCut.h"

//_________________
class MpdFemtoKtPairCut : public MpdFemtoBasePairCut {
  
 public:
  /// Default cut
  MpdFemtoKtPairCut();
  /// Constructor that accepts kTmin, kTmax, phiMin, phiMax values
  MpdFemtoKtPairCut(double ktLo = 0., double ktHi = 1e6,
		    double phiLow = -1e6, double phiHi = 1e6);
  /// Copy constructor
  MpdFemtoKtPairCut(const MpdFemtoKtPairCut& c);
  /// Assignment operator
  MpdFemtoKtPairCut& operator=(const MpdFemtoKtPairCut& c);
  /// Destructor
  virtual ~MpdFemtoKtPairCut();


  /// Construct report
  virtual MpdFemotString report();
  /// Construct report
  virtual MpdFemotString Report()                     { return report(); }
  /// Prepare a report using list of settings
  virtual TList *listSettings();
  /// Prepare a report using list of settings
  virtual TList *ListSettings()                       { return listSettings(); }
  /// Clone pair cut
  MpdFemotPairCut* clone()                            { return new MpdFemtoKtPairCut(*this); }
  /// Clone pair cut
  MpdFemotPairCut* Clone()                            { return clone(); }
  /// Set kT range of the pair
  /// \par lo Low kT value
  /// \oar hi High kT value
  void setKtRange(double lo, double hi)               { fKtRange[0] = lo; fKtRange[1] = hi; }
  /// Set kT range of the pair
  /// \par lo Low kT value
  /// \oar hi High kT value
  void SetKtRange(double lo, double hi)               { setKtRange(lo, hi); }
  /// Set emsission angle range w.r.t. reaction plane
  /// \par lo Low value of angle
  /// \par hi High value of angle
  void setPhiRange(double lo, double hi);
  /// Set emsission angle range w.r.t. reaction plane
  /// \par lo Low value of angle
  /// \par hi High value of angle
  void SetPhiRange(double lo, double hi)              { setPhiRange(lo, hi); }
  /// Check if the pair will pass the cut
  virtual bool pass(const MpdFemotPair* pair);
  /// Check if the pair will pass the cut
  virtual bool pass(const MpdFemotPair* pair, double aRPAngle);

 protected:

  /// Reduced pair transverse momentum range
  float fKtRange[2];
  /// Emission angle range w.r.t. reaction plane
  float fPhiRange[2];

  ClassDef(MpdFemtoKtPairCut, 0);
};

#endif // MpdFemtoKtPairCut_h

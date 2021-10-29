/**
 * \class MpdFemtoCorrFctn3DLCMSSym
 * \brief 3D correlation function in Bertsch-Pratt coordinate system
 *
 * Three-dimensional correlation function of identical particles
 * in Bertsch-Pratt coordinate system (out-side-long)
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * \email nigmatkulov@gmail.com
 */

#ifndef MpdFemtoCorrFctn3DLCMSSym_h
#define MpdFemtoCorrFctn3DLCMSSym_h

// Forward declarations
class TH3F;
class MpdFemtoPairCut;

// MpdFemtoMaker headers
#include "MpdFemtoBaseCorrFctn.h"

//_________________
class MpdFemtoCorrFctn3DLCMSSym : public MpdFemtoBaseCorrFctn {
 public:

  /// Build the correlation function with parameters.
  ///
  /// \param title The title with which to give the output
  /// \param nbins The number of bins in each direction of , and q
  MpdFemtoCorrFctn3DLCMSSym(const char* title, const int nbins, const float QHi);
  /// Copy Constructor
  MpdFemtoCorrFctn3DLCMSSym(const MpdFemtoCorrFctn3DLCMSSym& aCorrFctn);
  /// Assignment operator
  MpdFemtoCorrFctn3DLCMSSym& operator=(const MpdFemtoCorrFctn3DLCMSSym& aCorrFctn);
  /// Destructor
  virtual ~MpdFemtoCorrFctn3DLCMSSym();

  /// Make report
  virtual MpdFemtoString report();
  /// Add real pair (from the same event)
  virtual void addRealPair(MpdFemtoPair* aPair);
  /// Add mixed pair (from event mixing)
  virtual void addMixedPair(MpdFemtoPair* aPair);

  /// Finish
  virtual void finish();

  /// Return numerator
  TH3F* numerator() {
    return (mNumerator) ? mNumerator : nullptr;
  }
  /// Return denominator
  TH3F* denominator() {
    return (mDenominator) ? mDenominator : nullptr;
  }
  /// Return numerator weighted by qinv
  TH3F* numeratorW() {
    return (mNumeratorW) ? mNumeratorW : nullptr;
  }
  /// Return denominator weighted by qinv
  TH3F* denominatorW() {
    return (mDenominatorW) ? mDenominatorW : nullptr;
  }

  /// Write all histograms
  void writeOutHistos();
  /// Retrieve list of histograms
  virtual TList* getOutputList();

  /// Set which system to use:
  /// \param true LCMS (default)
  /// \param false PRF
  void setUseLCMS(bool useLCMS) {
    mUseLCMS = useLCMS;
  }
  /// Retrieve which system is used for calculations
  bool isUseLCMS() {
    return mUseLCMS;
  }

  /// Clone correlation function
  virtual MpdFemtoBaseCorrFctn* clone() const {
    return new MpdFemtoCorrFctn3DLCMSSym(*this);
  }

private:

  /// Numerator
  TH3F* mNumerator;
  /// Denominator
  TH3F* mDenominator;
  /// Qinv-weighted numerator
  TH3F* mNumeratorW;
  /// Qinv-weighted denominator
  TH3F* mDenominatorW;

  /// False - use PRF, True - use LCMS
  bool mUseLCMS;

  ClassDef(MpdFemtoCorrFctn3DLCMSSym, 1)
};

#endif // MpdFemtoCorrFctn3DLCMSSym_h

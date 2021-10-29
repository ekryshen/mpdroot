/**
 * \class MpdFemtoModelBPLCMS3DCorrFctnKt
 * \brief Three-dimensional Bertsch-Pratt correlation function in LCMS for the model estimations
 *
 * Three-dimensional Bertsch-Pratt correlation function which will store
 * numerators (from same event) and denominatros (from mixed events) with
 * and without femtoscopic weights
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * \email nigmatkulov@gmail.com
 */

#ifndef MpdFemtoModelBPLCMS3DCorrFctnKt_h
#define MpdFemtoModelBPLCMS3DCorrFctnKt_h

// MpdFemtoMaker headers
// Base
#include "MpdFemtoBaseCorrFctn.h"

// C++ headers
#include <vector>

// Forward declarations
class MpdFemtoPair;
class MpdFemtoModelManager;
class TH3F;

//_________________
class MpdFemtoModelBPLCMS3DCorrFctnKt : public MpdFemtoBaseCorrFctn {
 public:
  /// Parametrized constructor
  ///
  /// \param title  Name of the histogram
  /// \param nBins  Number of bins (will be used for out, side anMpd long projections)
  /// \param qLo    Minimum value of the q
  /// \param qHi    Maximum value of the q
  /// \param ktBins Number of kT bins used in the analysis
  /// \param ktLo   Minimum value of kT
  /// \param ktHi   Maximum value of kT
  ///
  MpdFemtoModelBPLCMS3DCorrFctnKt(const char* title = "hBPLCMSCorrFctn",
				  const int& nBins = 80, const double& qLo = -0.4, const double& qHi = 0.4,
				  const int& ktBins = 10, const double& ktLo = 0.05, const double& ktHi = 1.05,
				  const bool isUseDenominator = true);
  /// Copy constructor
  MpdFemtoModelBPLCMS3DCorrFctnKt(const MpdFemtoModelBPLCMS3DCorrFctnKt& corrFctn);
  /// Assignment operator
  MpdFemtoModelBPLCMS3DCorrFctnKt& operator=(const MpdFemtoModelBPLCMS3DCorrFctnKt& corrFctn);
  /// Destructor
  virtual ~MpdFemtoModelBPLCMS3DCorrFctnKt();

  /// Method that allows front-loading model manager
  virtual void connectToManager(MpdFemtoModelManager *manager);

  /// Set kT range (nbins, ktLow, ktHi)
  void setKtRange(const int& nKtBins = 10, const float& kTLow = 0.05, const float& kTHi = 1.05);

  /// Make report
  virtual MpdFemtoString report();

  /// Add real pair
  virtual void addRealPair(MpdFemtoPair* pair);
  /// Add mixed pair
  virtual void addMixedPair(MpdFemtoPair* pair);

  /// Begin event
  virtual void eventBegin(const MpdFemtoEvent* event);
  /// Event end
  virtual void eventEnd(const MpdFemtoEvent* event);
  /// Finish
  virtual void finish();

  /// Return i-th histogram for numerator
  TH3F* numerator(int i) {
    return ( mNumerator.at(i) ) ? mNumerator.at(i) : nullptr;
  }

  /// Return i-th histogram for denominator
  TH3F* denominator(int i) {
    return ( mDenominator.at(i) ) ? mDenominator.at(i) : nullptr;
  }

  /// Return i-th weighted numerator
  TH3F* numeratorWeighted(int i) {
    return ( mNumeratorWeighted.at(i) ) ? mNumeratorWeighted.at(i) : nullptr;
  }

  /// Return i-th weighted numerator
  TH3F* denominatorWeighted(int i) {
    return ( mDenominatorWeighted.at(i) ) ? mDenominatorWeighted.at(i) : nullptr;
  }

  /// Return i-th qInv weighted numerator
  TH3F* numeratorQinvWeighted(int i) {
    return ( mNumeratorQinvWeighted.at(i) ) ? mNumeratorQinvWeighted.at(i) : nullptr;
  }

  /// Return i-th weighted numerator
  TH3F* denominatorQinvWeighted(int i) {
    return ( mDenominatorQinvWeighted.at(i) ) ? mDenominatorQinvWeighted.at(i) : nullptr;
  }

  /// Write histograms
  virtual void writeOutHistos();

  /// Return output list
  virtual TList* getOutputList();

  /// Set use/not use denominators
  void setUseDenominator(bool useDen) { mIsUseDenominator = useDen; }
  
  /// Clone correlation function
  virtual MpdFemtoModelBPLCMS3DCorrFctnKt* clone() const {
    return new MpdFemtoModelBPLCMS3DCorrFctnKt(*this);
  }

protected:

  /// Pointer to the model manager that performs femtoscopic
  /// weight calculation
  MpdFemtoModelManager *mManager;

  /// Numerator made with pairs from the same event with femtoscopic weight
  std::vector< TH3F* > mNumeratorWeighted;
  /// Numerator made with pairs from the same event without weights
  std::vector< TH3F* > mNumerator;
  /// Numerator made with pairs from the same event with Qinv weight
  std::vector< TH3F* > mNumeratorQinvWeighted;
  /// Denominator made with mixed pairs with femtoscopic weight
  std::vector< TH3F* > mDenominatorWeighted;
  /// Denominator made with mixed pairs without weights
  std::vector< TH3F* > mDenominator;
  /// Denominator made with pairs from the same event with Qinv weight
  std::vector< TH3F* > mDenominatorQinvWeighted;

  /// Number of the kT bins
  int mNKtBins;
  /// kT step
  float mKtStep;
  /// kT range
  float mKtRange[2];

  /// Use denominator (default true)
  bool mIsUseDenominator;

  ClassDef(MpdFemtoModelBPLCMS3DCorrFctnKt, 2)
};

#endif // MpdFemtoModelBPLCMS3DCorrFctnKt_h

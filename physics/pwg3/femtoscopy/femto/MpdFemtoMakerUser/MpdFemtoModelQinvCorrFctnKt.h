/**
 * \class MpdFemtoModelQinvCorrFctnKt
 * \brief One-dimensional correlation function with kT-binning for models
 *
 * One-dimensional correlation function with correlation functions binned
 * according to the user request
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * \email nigmatkulov@gmail.com
 */

#ifndef MpdFemtoModelQinvCorrFctnKt_h
#define MpdFemtoModelQinvCorrFctnKt_h

// MpdFemtoMaker headers
#include "MpdFemtoBaseCorrFctn.h"
#include "MpdFemtoBasePairCut.h"

// C++ headers
#include <vector>

// Forward declarations
class MpdFemtoPair;
class MpdFemtoModelManager;
class TH1D;

//_________________
class MpdFemtoModelQinvCorrFctnKt : public MpdFemtoBaseCorrFctn {
 public:

  /// Parametrized constructor
  ///
  /// \param title  Name of the histogram
  /// \param nBins  Number of qInv bins 
  /// \param qLo    Minimum value of the qInv
  /// \param qHi    Maximum value of the qInv
  /// \param ktBins Number of kT bins used in the analysis
  /// \param ktLo   Minimum value of kT
  /// \param ktHi   Maximum value of kT
  ///
  MpdFemtoModelQinvCorrFctnKt(const char* title = "hQinv",
			      const int& nbins = 100, const double& QinvLo = 0., const double& QinvHi = 1.,
			      const int& ktBins = 10, const double& KtLo = 0.05, const double& KtHi = 1.05,
			      const bool useDenominator = true);
  /// Copy constructor
  MpdFemtoModelQinvCorrFctnKt(const MpdFemtoModelQinvCorrFctnKt& copy);
  /// Assignment operator
  MpdFemtoModelQinvCorrFctnKt& operator=(const MpdFemtoModelQinvCorrFctnKt& copy);
  /// Destructor
  virtual ~MpdFemtoModelQinvCorrFctnKt();

  /// Method that allows front-loading model manager
  virtual void connectToManager(MpdFemtoModelManager *manager);

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

  /// Write histograms
  void writeOutHistos();

  /// Return output list
  virtual TList* getOutputList();

  /// Set kT bins
  void setKtRange(const int& nbins = 10, const double& ktLo = 0.05, const double& ktHi = 1.05);

  /// Return i-th numerator
  TH1D *numerator(unsigned int i) {
    return ( mNumerator.at(i)) ? mNumerator.at(i) : nullptr;
  }
  /// Return i-th weighted numerator
  TH1D *numeratorWeighted(unsigned int i) {
    return ( mNumeratorWeighted.at(i)) ? mNumeratorWeighted.at(i) : nullptr;
  }
  /// Return i-th denominator
  TH1D *denominator(unsigned int i) {
    return ( mDenominator.at(i)) ? mDenominator.at(i) : nullptr;
  }
  /// Return i-th denominator
  TH1D *denominatorWeighted(unsigned int i) {
    return ( mDenominatorWeighted.at(i)) ? mDenominatorWeighted.at(i) : nullptr;
  }

  /// Set write denominator
  void setUseDenominator(bool isUse) { mIsUseDenominator = isUse; }

  /// Clone correlation function
  virtual MpdFemtoModelQinvCorrFctnKt* clone() const {
    return new MpdFemtoModelQinvCorrFctnKt(*this);
  }

 private:

  /// Pointer to the model manager that performs femtoscopic
  /// weight calculation
  MpdFemtoModelManager *mManager;

  /// Numerator
  std::vector< TH1D* > mNumerator;
  /// Numerator weighted
  std::vector< TH1D* > mNumeratorWeighted;
  /// Denominator
  std::vector< TH1D* > mDenominator;
  /// Denominator weighted
  std::vector< TH1D* > mDenominatorWeighted;

  /// Number of kT bins to study
  int mNKtBins;
  /// Range of kT to plot [min,max]
  double mKtRange[2];

  /// kT bin width
  double mKtStep;

  /// Use and write denominators
  bool mIsUseDenominator;

  ClassDef(MpdFemtoModelQinvCorrFctnKt, 1)
};

#endif // #define MpdFemtoModelQinvCorrFctnKt_h

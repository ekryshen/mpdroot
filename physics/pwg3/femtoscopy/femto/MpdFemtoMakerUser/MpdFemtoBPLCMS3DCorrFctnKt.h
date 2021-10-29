/**
 * \class MpdFemtoBPLCMS3DCorrFctnKt
 * \brief A class to calculate 3D correlation for pairs of identical particles.
 *
 * It also stored the weighted qinv per bin denominator histogram for the coulomb
 * correction.
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * \email nigmatkulov@gmail.com
 */

#ifndef MpdFemtoBPLCMS3DCorrFctnKt_h
#define MpdFemtoBPLCMS3DCorrFctnKt_h

// MpdFemtoMaker headers
#include "MpdFemtoBaseCorrFctn.h"
#include "MpdFemtoBasePairCut.h"
#include "MpdFemtoString.h"

// ROOT headers
#include "TH3.h"

//_________________
class MpdFemtoBPLCMS3DCorrFctnKt : public MpdFemtoBaseCorrFctn {
 public:
  /// Parametrized constructor
  ///
  /// \param title  Name of the histogram
  /// \param nBins  Number of bins (will be used for out, side and long projections)
  /// \param qLo    Minimum value of the q
  /// \param qHi    Maximum value of the q
  /// \param ktBins Number of kT bins used in the analysis
  /// \param ktLo   Minimum value of kT
  /// \param ktHi   Maximum value of kT
  ///
  MpdFemtoBPLCMS3DCorrFctnKt(const char* title = "hBPLCMSCorrFctn",
			     const int& nBins = 80, const double& qLo = -0.4, const double& qHi = 0.4,
			     const int& ktBins = 10, const double& ktLo = 0.05, const double& ktHi = 1.05);
  /// Copy constructor
  MpdFemtoBPLCMS3DCorrFctnKt(const MpdFemtoBPLCMS3DCorrFctnKt& copy);
  /// Assignment operator
  MpdFemtoBPLCMS3DCorrFctnKt& operator=(const MpdFemtoBPLCMS3DCorrFctnKt& aCorrFctn);
  /// Destructor
  virtual ~MpdFemtoBPLCMS3DCorrFctnKt();

  /// Set histogram parameters (nbins, qLow, qHi)
  void setHistoParameters(const int& nBins = 80, const float& qLo = -0.4, const float& qHi = 0.4);
  /// Set kT range (nbins, ktLow, ktHi)
  void setKtRange(const int& nKtBins = 10, const float& kTLow = 0.05, const float& kTHi = 1.05);

  /// Report
  virtual MpdFemtoString report();
  /// Add real pair
  virtual void addRealPair(MpdFemtoPair* pair);
  /// Add mixed pair
  virtual void addMixedPair(MpdFemtoPair* pair);

  /// Event begin
  virtual void eventBegin(const MpdFemtoEvent* /* event*/) {
    /* empty */
  }
  /// Event end
  virtual void eventEnd(const MpdFemtoEvent* /* event*/) {
    /* empty */
  }
  /// Finish
  virtual void finish();

  /// Return i-th histogram for numerator
  TH3F* numerator(int i) {
    return (mNumerator[i]) ? mNumerator[i] : nullptr;
  }
  /// Return i-th histogram for denominator
  TH3F* denominator(int i) {
    return (mDenominator[i]) ? mDenominator[i] : nullptr;
  }
  /// Return i-th histogram with the q_{inv} weight
  TH3F* qInvWeight(int i) {
    return (mQinvHisto[i]) ? mQinvHisto[i] : nullptr;
  }

  /// Write histograms
  void writeOutHistos();
  /// Get output list
  virtual TList* getOutputList();
  /// Clone correalation function
  virtual MpdFemtoBaseCorrFctn* clone() const {
    return new MpdFemtoBPLCMS3DCorrFctnKt(*this);
  }

 private:

  /// Numerator
  TH3F* mNumerator[20];
  /// Denominator
  TH3F* mDenominator[20];
  /// qInv weighted denominator
  TH3F* mQinvHisto[20];

  /// Number of bins in the outward direction
  int mNQbins;
  /// qOut range
  float mQRange[2];

  /// Number of the kT bins
  int mNKtBins;
  /// kT step
  float mKtStep;
  /// kT range
  float mKtRange[2];

  ClassDef(MpdFemtoBPLCMS3DCorrFctnKt, 1)
};

#endif // #ifndef MpdFemtoBPLCMS3DCorrFctnKt_h

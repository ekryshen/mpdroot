/**
 * \class MpdFemtoQinvCorrFctnKt
 * \brief One-dimensional correlation function with kT-binning
 *
 * One-dimensional correlation function with correlation functions binned
 * according to the user request
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * \email nigmatkulov@gmail.com
 */

#ifndef MpdFemtoQinvCorrFctnKt_h
#define MpdFemtoQinvCorrFctnKt_h

// MpdFemtoMaker headers
#include "MpdFemtoBaseCorrFctn.h"
#include "MpdFemtoBasePairCut.h"

// C++ headers
#include <vector>

// Forward declaration
class TH1D;

//_________________
class MpdFemtoQinvCorrFctnKt : public MpdFemtoBaseCorrFctn {
 public:

  /// Constructor
  MpdFemtoQinvCorrFctnKt(const char* title = "hQinv",
			 const int& nbins = 100, const double& QinvLo = 0., const double& QinvHi = 1.,
			 const int& nCFs = 10, const double& KtLo = 0.05, const double& KtHi = 1.05);
  /// Copy constructor
  MpdFemtoQinvCorrFctnKt(const MpdFemtoQinvCorrFctnKt& copy);
  /// Assignment operator
  MpdFemtoQinvCorrFctnKt& operator=(const MpdFemtoQinvCorrFctnKt& copy);
  /// Destructor
  virtual ~MpdFemtoQinvCorrFctnKt();

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
  /// Return i-th denominator
  TH1D *denominator(unsigned int i) {
    return ( mDenominator.at(i)) ? mDenominator.at(i) : nullptr;
  }

  /// Clone correlation function
  virtual MpdFemtoQinvCorrFctnKt* clone() const {
    return new MpdFemtoQinvCorrFctnKt(*this);
  }

private:

  /// Numerators
  std::vector< TH1D* > mNumerator;
  /// Denominators
  std::vector< TH1D* > mDenominator;

  /// Number of kT bins to study
  int mNKtBins;
  /// Range of kT to plot [min,max]
  double mKtRange[2];

  /// kT bin width
  double mKtStep;

  ClassDef(MpdFemtoQinvCorrFctnKt, 2)
};

#endif // #define MpdFemtoQinvCorrFctnKt_h

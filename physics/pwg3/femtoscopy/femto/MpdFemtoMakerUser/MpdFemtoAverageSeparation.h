/**
 * \class MpdFemtoAverageSeparation
 * \brief One-dimensional distr. of Average Separation
 *
 *
 * \author Malinina Ludmila (SINP MSU & JINR)
 * \date  Sep 2020, 2020
 * \email lmalinin@cern.ch
 */

#ifndef MpdFemtoAverageSeparation_h
#define MpdFemtoAverageSeparation_h

// MpdFemtoMaker headers
// Base
#include "MpdFemtoBaseCorrFctn.h"

// ROOT headers
#include "TH1.h"
#include "TVector3.h"

// C++ headers
#include <vector>

// Forward declarations
class MpdFemtoPair;

//_________________
class MpdFemtoAverageSeparation : public MpdFemtoBaseCorrFctn {
 public:

  /// Parametrized constructor
  ///
  /// \param title     Name of the histogram
  /// \param nBinsEta  Number of Average Separation bins
  /// \param mEtaLo    Minimum value of Average Separation
  /// \param mEtaHi    Maximum value of Average Separation
  ///
  MpdFemtoAverageSeparation(const char* title = "hAverageSepartion",
				    const int& nBinsAvSep = 200, const double& AvSep = 0, const double& aVSep = 100);
  
  /// Destructor
  virtual ~MpdFemtoAverageSeparation();

  
  void setHistoParameters(const int& nBinsAvSep = 200, const double& AvSep = 0., const double& avSep = 100.);
			 
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
  TH1F* numerator(unsigned int i) {
    return ( mNumerator.at(i)) ? mNumerator.at(i) : nullptr;
  }
  /// Return i-th histogram for denominator
  TH1F* denominator(unsigned int i) {
    return ( mDenominator.at(i)) ? mDenominator.at(i) : nullptr;
  }

  /// Write histograms
  virtual void writeOutHistos();

  /// Return output list
  virtual TList* getOutputList();
  /// Clone correlation function
  virtual MpdFemtoAverageSeparation* clone() const {
    return new MpdFemtoAverageSeparation(*this);
  }

protected:

  /// Numerator made with pairs from the same event
  std::vector< TH1F* > mNumerator;
  /// Denominator made with mixed pairs
  std::vector< TH1F* > mDenominator;

  /// Number of eta bins
  int mAvSepBins;
  /// Eta range
  double mAvSepRange[2];
 

  ClassDef(MpdFemtoAverageSeparation, 1)
};

#endif // MpdFemtoAverageSeparation_h

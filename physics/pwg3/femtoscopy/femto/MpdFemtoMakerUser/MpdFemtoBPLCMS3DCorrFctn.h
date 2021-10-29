/*
 * \class MpdFemtoBPLCMS3DCorrFctn
 * \brief A class that calculates 3D correlation function of identical particles
 *
 * It also stored the weighted qinv per bin denominator histogram for the coulomb
 * correction.
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 11, 2020
 * \email nigmatkulov@gmail.com
*/

#ifndef MpdFemtoBPLCMS3DCorrFctn_h
#define MpdFemtoBPLCMS3DCorrFctn_h

// MpdFemtoMaker headers
#include "MpdFemtoBaseCorrFctn.h"
#include "MpdFemtoBasePairCut.h"

// ROOT headers
class TH3F;
class TList;

//_________________
class MpdFemtoBPLCMS3DCorrFctn : public MpdFemtoBaseCorrFctn {
  
 public:
  /// Parametrized constructor
  MpdFemtoBPLCMS3DCorrFctn(const char* title = "hBPLCMSCorrFctn", const int& nbins = 80,
			   const double& qLo = -0.4, const double& qHi = 0.4);
  /// Copy constructor
  MpdFemtoBPLCMS3DCorrFctn(const MpdFemtoBPLCMS3DCorrFctn& aCorrFctn);
  /// Destructor
  virtual ~MpdFemtoBPLCMS3DCorrFctn();
  /// Assignment operator
  MpdFemtoBPLCMS3DCorrFctn& operator=(const MpdFemtoBPLCMS3DCorrFctn& aCorrFctn);

  /// Construct report
  virtual MpdFemtoString report();
  /// Add pair from the real event
  virtual void addRealPair( MpdFemtoPair* aPair);
  /// Add pair from the mixed event
  virtual void addMixedPair( MpdFemtoPair* aPair);
  /// Finish
  virtual void finish();

  /// Return numerator
  TH3F* numerator()                                  { return ( mNumerator ) ? mNumerator : nullptr; }
  /// Return denominator
  TH3F* denominator()                                { return ( mDenominator ) ? mDenominator : nullptr; }
  /// Return qInv-weighted denominator
  TH3F* qInvHisto()                                  { return ( mQinvHisto ) ? mQinvHisto : nullptr; }

  /// Write histograms to the file
  void writeOutHistos();                              
  /// Return list of output list
  virtual TList* getOutputList();
  /// Clone BPLCMS histograms
  virtual MpdFemtoBaseCorrFctn* clone() const { return new MpdFemtoBPLCMS3DCorrFctn( *this ); }

 private:

  /// Numerator
  TH3F* mNumerator;
  /// Denominator
  TH3F* mDenominator;
  /// Denominator with qInv weights
  TH3F* mQinvHisto;

  ClassDef(MpdFemtoBPLCMS3DCorrFctn, 1);

};

#endif
